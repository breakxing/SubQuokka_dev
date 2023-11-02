#include <omp.h>
#include <unistd.h> //pread, pwrite

#include <complex>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>
#include <algorithm>
#include "circuit.hpp"
#include <functional>
#include "MPI_Runner.hpp"


#define outer_loop_m0(inner_loop) {\
    for (long long i = 0; i < env.thread_state; i += env.chunk_state) {\
        inner_loop\
        for (auto &x : task.fd_offset_using) {\
            x += env.chunk_size;\
        }\
    }\
}

#define outer_loop_m1(inner_loop, m0) {\
    long long half_off = env.qubit_offset[m0];\
    long long half_off_size = env.qubit_size[m0];\
    long long off = half_off * 2;\
    for (long long i = 0; i < env.thread_state; i += off) {\
        for (long long j = 0; j < half_off; j += env.chunk_state) {\
            inner_loop\
            for (auto &x : task.fd_offset_using) {\
                x += env.chunk_size;\
            }\
        }\
        for (auto &x : task.fd_offset_using) {\
            x += half_off_size;\
        }\
    }\
}

#define outer_loop_m2(inner_loop, m0, m1) {\
    long long half_off_m0 = env.qubit_offset[m0];\
    long long half_off_m1 = env.qubit_offset[m1];\
    long long half_off_size_m0 = env.qubit_size[m0];\
    long long half_off_size_m1 = env.qubit_size[m1];\
    long long off_m0 = half_off_m0 * 2;\
    long long off_m1 = half_off_m1 * 2;\
    for (long long i = 0; i < env.thread_state; i += off_m1) {\
        for (long long j = 0; j < half_off_m1; j += off_m0) {\
            for (long long k = 0; k < half_off_m0; k += env.chunk_state) {\
                inner_loop\
                for (auto &x : task.fd_offset_using) {\
                    x += env.chunk_size;\
                }\
            }\
            for (auto &x : task.fd_offset_using) {\
                x += half_off_size_m0;\
            }\
        }\
        for (auto &x : task.fd_offset_using) {\
            x += half_off_size_m1;\
        }\
    }\
}

#define innerloop {\
    int fsz = task.fd_using.size();\
    int osz = task.fd_offset_using.size();\
    for (int ii = 0; ii < fsz; ii++) {\
        for (int jj = 0; jj < osz; jj++) {\
            if(pread(env.fd_arr[task.fd_using[ii]], &(task.buffer[(ii * osz + jj) * env.chunk_state]), env.chunk_size, task.fd_offset_using[jj]))\
                ;\
        }\
    }\
    g->run(task.buffer);\
    for (int ii = 0; ii < fsz; ii++) {\
        for (int jj = 0; jj < osz; jj++) {\
            if(pwrite(env.fd_arr[task.fd_using[ii]], &(task.buffer[(ii * osz + jj) * env.chunk_state]), env.chunk_size, task.fd_offset_using[jj]))\
                ;\
        }\
    }\
}

#define innerloop_sub {\
    int fsz = task.fd_using.size();\
    int osz = task.fd_offset_using.size();\
    for (int ii = 0; ii < fsz; ii++) {\
        for (int jj = 0; jj < osz; jj++) {\
            if(pread(env.fd_arr[task.fd_using[ii]], &(task.buffer[(ii * osz + jj) * env.chunk_state]), env.chunk_size, task.fd_offset_using[jj]))\
                ;\
        }\
    }\
    for (auto &g : subcircuit){\
        g->run(task.buffer);\
    }\
    for (int ii = 0; ii < fsz; ii++) {\
        for (int jj = 0; jj < osz; jj++) {\
            if(pwrite(env.fd_arr[task.fd_using[ii]], &(task.buffer[(ii * osz + jj) * env.chunk_state]), env.chunk_size, task.fd_offset_using[jj]))\
                ;\
        }\
    }\
}

using namespace std;

template<class T>
struct aligned_buffer_allocator: public std::allocator<T>{
    T* allocate (std::size_t n);
    void deallocate (T* p, std::size_t n);
};

inline int findCorrespond(int tid, int targ) {
    int targMask = 1 << targ;
    return tid + targMask;
}

thread_MPI_task::thread_MPI_task(int tid) {
    request_size = 16;
    request.resize(request_size);
    for (int i = 0; i < seg.chunk; i++) {
        fd_table.push_back(tid);
        fd_offset_table.push_back(0);
    }

    for (int i = 0; i < seg.middle; i++) {
        fd_table.push_back(tid);
        fd_offset_table.push_back(env.qubit_size[i + seg.chunk]);
    }

    for (int i = 0; i < seg.file; i++) {
        int corr = findCorrespond(tid, i);

        if (corr > tid) {
            fd_table.push_back(corr - tid);
            fd_offset_table.push_back(0);
        }
    }
    for(int i = 0;i < seg.mpi;i++)
    {
        fd_table.push_back(tid);
        fd_offset_table.push_back(0);
    }
}

MPI_Runner::MPI_Runner() {
    int provided;
    if(MPI_Init_thread(NULL, NULL,MPI_THREAD_MULTIPLE,&provided)!=MPI_SUCCESS) exit(-1);
    if(provided < MPI_THREAD_MULTIPLE)
    {
        printf("The threading support level is lesser than that demanded.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if(MPI_Comm_rank(MPI_COMM_WORLD, &env.rank)!=MPI_SUCCESS)
    {
        exit(-1);
    }
    for (int i = 0; i < env.num_thread; i++) {
        thread_tasks.push_back(thread_MPI_task(i));
        if(env.is_MPI)
            thread_tasks[i].buffer.resize(16 * env.chunk_state);
        else if (env.is_directIO) {
            thread_tasks[i].buffer.resize(8 * env.chunk_state); // 64 for VSWAP_6_6, 8 for normal use.
            // posix_memalign(thread_tasks[i].buffer, 4096, 8 * env.chunk_state * sizeof(complex<double>)); // 64 for VSWAP_6_6, 8 for normal use.
        }
        else
            thread_tasks[i].buffer.resize(8 * env.chunk_state); // 64 for VSWAP_6_6, 8 for normal use.
    }
}
void MPI_Runner::setFD(thread_MPI_task &task, Gate *&gate) {
    task.fd_using.resize(0);
    task.fd_offset_using.resize(0);
    auto targ = gate->targs;
    int file_count = gate->file_count;
    int middle_count = gate->middle_count;
    int chunk_count = gate->chunk_count;
    int nonchunk_count = file_count + middle_count;

    int nonchunk_qbit[6];
    for (int i = 0; i < nonchunk_count; i++){
        nonchunk_qbit[i] = targ[i + chunk_count];
    }
    
    vector<int> temp_fd_table;
    vector<int> temp_fd;
    for (int i = middle_count; i < middle_count + file_count; i++) {
        temp_fd_table.push_back(task.fd_table[nonchunk_qbit[i]]);
    }

    int file_range = 1 << file_count;
    for (int i = 0; i < file_range; i++) {
        int fd = task.fd_table[0]; // tid
        for (int j = 0; j < file_count; j++) {
            if (i & (1<<j)) {
                fd += temp_fd_table[j];
            }
        }
        temp_fd.push_back(fd);
    }

    vector<long long> temp_fd_offset_table;
    vector<long long> temp_fd_offset;
    for (int i = 0; i < middle_count; i++) {
        temp_fd_offset_table.push_back(task.fd_offset_table[nonchunk_qbit[i]]);
    }

    int middle_range = 1 << middle_count;
    for (int i = 0; i < middle_range; i++) {
        long long fd_off = 0;
        for (int j = 0; j < middle_count; j++) {
            if (i & (1<<j)) {
                fd_off += temp_fd_offset_table[j];
            }
        }
        temp_fd_offset.push_back(fd_off);
    }

    task.fd_using = temp_fd;
    task.fd_offset_using = temp_fd_offset;

    return;
}

void MPI_Runner::setFD_sub(thread_MPI_task &task) {
    task.fd_using.resize(0);
    task.fd_offset_using.resize(0);

    int fd = task.fd_table[0]; // tid

    task.fd_using.push_back(fd);
    task.fd_offset_using.push_back(0);

    return;
}

static inline bool skipThread(int tid, vector<int> &targ) {
    for (int q : targ) {
        if (isFile(q) && ((tid & (1 << (q - seg.middle - seg.chunk))) > 0))
            return true;
    }
    return false;
}

void outer_loop_m3(thread_MPI_task &task, function<void()> &inner_loop, int m0, int m1, int m2) {
    long long half_off_m0 = env.qubit_offset[m0];
    long long half_off_m1 = env.qubit_offset[m1];
    long long half_off_m2 = env.qubit_offset[m2];
    
    long long half_off_size_m0 = env.qubit_size[m0];
    long long half_off_size_m1 = env.qubit_size[m1];
    long long half_off_size_m2 = env.qubit_size[m2];
    
    long long off_m0 = half_off_m0 * 2;
    long long off_m1 = half_off_m1 * 2;
    long long off_m2 = half_off_m2 * 2;
    for (long long i = 0; i < env.thread_state; i += off_m2) {
        for (long long j = 0; j < half_off_m2; j += off_m1) {
            for (long long k = 0; k < half_off_m1; k += off_m0) {
                for (long long m = 0; m < half_off_m0; m += env.chunk_state) {
                    inner_loop();
                    for (auto &x : task.fd_offset_using) {
                        x += env.chunk_size;
                    }
                }
                for (auto &x : task.fd_offset_using) {
                    x += half_off_size_m0;
                }
            }
            for (auto &x : task.fd_offset_using) {
                x += half_off_size_m1;
            }
        }
        for (auto &x : task.fd_offset_using) {
            x += half_off_size_m2;
        }
    }
}

void outer_loop_m4up_helper(thread_MPI_task &task, function<void()> &inner_loop, vector<int> &m, int level){
    if (level == -1){
        inner_loop();
        for (auto &x : task.fd_offset_using) {
            x += env.chunk_size;
        }
        return;
    }
    long long range = env.qubit_offset[m[level]];
    long long incre = level == 0 ? env.chunk_state : 2 * env.qubit_offset[m[level-1]];
    long long range_size = range * sizeof(complex<double>);

    for (long long i = 0; i < range; i += incre) {
        outer_loop_m4up_helper(task, inner_loop, m, level-1);
    }
    for (auto &x : task.fd_offset_using) {
        x += range_size;
    }
    return;
}

void outer_loop_m4up(thread_MPI_task &task, function<void()> &inner_loop, vector<int> m) {
    long long incre = 2 * env.qubit_offset[m.back()];

    for (long long i = 0; i < env.thread_state; i += incre) {
        outer_loop_m4up_helper(task, inner_loop, m, m.size()-1);
    }
}

void MPI_Runner::run(vector<Gate *> &circuit) {
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto task = thread_tasks[tid];
        task.tid = tid;
        // static int count = 0;
        for (auto &g : circuit) {
            // if(tid == 0)
            //     cout << ++count << "th gate doing...... " << endl;
            vector<int> targ = g->targs;  // increasing
            int mpi_count = g->mpi_count;
            int file_count = g->file_count;
            int middle_count = g->middle_count;
            int chunk_count = g->chunk_count;
            task.fd_using.resize(0);
            task.fd_offset_using.resize(0);
            if(mpi_count == 0)
            {
                setFD(task, g);
                if (file_count > 0 && skipThread(tid, targ)) 
                {
                    #pragma omp barrier
                    continue;
                }
                switch(middle_count){
                    case 0:
                        outer_loop_m0(innerloop)
                        break;

                    case 1:
                        outer_loop_m1(innerloop, targ[chunk_count])
                        break;

                    case 2:
                        outer_loop_m2(innerloop, targ[chunk_count], targ[chunk_count+1])
                        break;
                    default:
                        exit(-1);
                }
            }
            else if(mpi_count <= 2)
            {
                if(g->name == "SWAP_Gate")
                    MPI_Swap(task,g);
                else
                    MPI_gate_scheduler(task,g);
            }
            else
                exit(-1);
#pragma omp barrier
        }
    #pragma omp barrier
    #pragma omp single
    MPI_Barrier(MPI_COMM_WORLD);
    #pragma omp barrier
    }
    if(MPI_Finalize()!=MPI_SUCCESS)
    {
        cout<<"Error"<<endl;
        exit(-1);
    }
}

void MPI_Runner::run(vector<vector<Gate *>> &subcircuits) {
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto task = thread_tasks[tid];

        for (auto &subcircuit : subcircuits) {
            Gate *g = subcircuit[0];
            if (g->type == VSWAP) {
                vector<int> targ = subcircuit[0]->targs;  // increasing
                int file_count = subcircuit[0]->file_count;
                int middle_count = subcircuit[0]->middle_count;
                int chunk_count = subcircuit[0]->chunk_count;

                if (file_count > 0 && skipThread(tid, targ)) {
                    #pragma omp barrier
                    continue;
                }

                setFD(task, subcircuit[0]);

                switch(middle_count){
                    case 0:
                        outer_loop_m0(innerloop)
                        break;

                    case 1:
                        outer_loop_m1(innerloop, targ[chunk_count])
                        break;

                    case 2:
                        outer_loop_m2(innerloop, targ[chunk_count], targ[chunk_count+1])
                        break;
                        
                    // case 3:
                    //     ol = bind(outer_loop_m3, ref(task), ref(il), targ[chunk_count], targ[chunk_count+1], targ[chunk_count+2]);
                    //     break;
                    
                    default:
                        vector<int> m;
                        for (int i = chunk_count; i < chunk_count+middle_count; i++) {
                            m.push_back(targ[i]);
                        }
                        // ol = bind(outer_loop_m4up, ref(task), ref(il), m);
                        break;
                }
                // ol();
            }
            else {
                setFD_sub(task);
                outer_loop_m0(innerloop_sub)
            }
            #pragma omp barrier
        }
    }
}
void MPI_Runner::MPI_Swap(thread_MPI_task &task,Gate * &g)
{
    long long mpi_targ_mask_0 = 1 << (g->targs[1] - seg.N);
    long long mpi_targ_mask_1 = 1 << (g->targs[0] - seg.N);
    task.partner_using = {int(((long long)env.rank) ^ mpi_targ_mask_0)};
    task.fd_using = {env.fd_arr[task.tid],env.fd_arr[task.tid]};
    task.fd_offset_using = {0};
    long long loop_bound = env.thread_size;
    long long loop_stride = (env.thread_state == env.chunk_state)? env.chunk_size : (env.chunk_size << 1);
    bool firstround = true;
    if(isMpi(g->targs[0]))
    {
        if(((env.rank >> (g->targs[1] - seg.N)) & 1) == ((env.rank >> (g->targs[0] - seg.N)) & 1)) return;
        task.partner_using = {int(((long long)env.rank) ^ mpi_targ_mask_0 ^ mpi_targ_mask_1)};
        loop_stride = env.chunk_size;
    }
    else if(isFile(g->targs[0]))
    {
        long long file_mask = 1 << (g->targs[0] - seg.middle - seg.chunk);
        if(env.thread_state == env.chunk_state && (((env.rank >> (g->targs[1] - seg.N)) & 1) == ((task.tid >> (g->targs[0] - seg.middle - seg.chunk)) & 1))) return;
        task.fd_using = {env.fd_arr[task.tid & (~file_mask)],env.fd_arr[task.tid | file_mask]};
    }
    else if(isMiddle(g->targs[0]))
    {
        loop_stride = env.qubit_size[g->targs[0] + 1];
    }
    for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
    {
        if(isMiddle(g->targs[0]))
        {
            task.fd_offset_using = env.rank < task.partner_using[0]? vector<long long>{cur_offset + env.qubit_size[g->targs[0]]} : vector<long long>{cur_offset};
            for(long long j = 0;j < env.qubit_size[g->targs[0]]; j += env.chunk_size)
            {
                _thread_swap(task,g,firstround);
                update_offset(task,env.chunk_size);
                firstround = false;
            }
        }
        else
        {
            _thread_swap(task,g,firstround);
            update_offset(task,loop_stride);
            firstround = false;
        }
    }
    
}
void MPI_Runner::_thread_swap(thread_MPI_task &task,Gate * &g,bool firstround)
{
    int partner_rank = task.partner_using[0];
    long long file_mask = 1 << (g->targs[0] - seg.middle - seg.chunk);
    int member_th = env.rank > partner_rank? 1 : 0;
    int thread_member_th;//For all thread drive such as D F,otherwise it is 0
    int recv_tag;
    if(isFile(g->targs[0]))
    {
        thread_member_th = (env.thread_state == env.chunk_state)? 0 : (task.tid & (file_mask)? 1 : 0);
        recv_tag = (env.thread_state != env.chunk_state)? task.tid : task.tid ^ file_mask;
    }
    else
    {
        thread_member_th = 0;
        recv_tag = task.tid;
    }
        MPI_Irecv(&task.buffer[(!member_th) * env.chunk_state], env.chunk_state, MPI_DOUBLE_COMPLEX, partner_rank,recv_tag, MPI_COMM_WORLD,&task.request[!member_th]);
        if(!firstround)
            MPI_Wait(&task.request[member_th],MPI_STATUS_IGNORE);
        if(pread(task.fd_using[!member_th],&task.buffer[member_th * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + thread_member_th * env.chunk_size));
        MPI_Isend(&task.buffer[member_th * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request[member_th]);
        MPI_Wait(&task.request[!member_th],MPI_STATUS_IGNORE);
        if(isChunk(g->targs[0]))
            g->run(task.buffer);
        if(pwrite(task.fd_using[!member_th],&task.buffer[(!member_th) * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + thread_member_th * env.chunk_size));
}
void MPI_Runner::MPI_gate_scheduler(thread_MPI_task &task,Gate * &g)
{
    long long mpi_targ_mask_0 = 1 << (g->targs[1] - seg.N);
    long long mpi_targ_mask_2 = 1 << (g->targs[0] - seg.N);
    long long loop_bound = env.thread_size;
    long long loop_stride;
    long long threshold;
    int rank0;
    int rank1;
    int rank2;
    int rank3;
    if(g->targs.size() == 1)
    {
        task.partner_using = {int(((long long)env.rank) ^ mpi_targ_mask_2)};
        task.fd_using = {env.fd_arr[task.tid]};
        task.fd_offset_using = {0};
        loop_stride = env.chunk_size << 1;
        if(env.thread_state == env.chunk_state && env.rank > task.partner_using[0])
            _thread_no_exec_MPI(task,1);
        else
        {
            for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
            {
                _thread_read1_recv1(task,g);
                update_offset(task,loop_stride);
            }
        }
    }
    else
    {
        if(isMpi(g->targs[0]))
        {
            rank0 = env.rank & (~(mpi_targ_mask_0 | mpi_targ_mask_2));
            rank1 = rank0 | mpi_targ_mask_2;
            rank2 = rank0 | mpi_targ_mask_0;
            rank3 = rank0 | mpi_targ_mask_0 | mpi_targ_mask_2;
            task.partner_using = {rank0,rank1,rank2,rank3};
            task.fd_using = {env.fd_arr[task.tid],env.fd_arr[task.tid],env.fd_arr[task.tid],env.fd_arr[task.tid]};
            task.fd_offset_using = {0,env.chunk_size,env.chunk_size << 1,env.chunk_size * 3};
            threshold = min(env.thread_state / env.chunk_state,(long long)4);
            loop_stride = (1 << int(log2(threshold))) * env.chunk_size;
            if(env.rank > task.partner_using[threshold - 1])
            {
                _thread_no_exec_MPI(task,1 << int(log2(threshold)));
            }
            else
            {
                for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
                {
                    _thread_read1_recv3(task,g,1 << int(log2(threshold)));
                    update_offset(task,loop_stride);
                }
            }

        }
        else if(isFile(g->targs[0]))
        {
            long long file_mask = 1 << (g->targs[0] - seg.middle - seg.chunk);
            rank0 = env.rank ^ mpi_targ_mask_0;
            rank1 = env.rank ^ mpi_targ_mask_0;
            task.partner_using = {rank0,rank1};
            task.fd_using = {env.fd_arr[task.tid & (~file_mask)],env.fd_arr[task.tid | file_mask]};
            task.fd_offset_using = {0,0};
            threshold = min(env.thread_state / env.chunk_state,(long long)4);
            loop_stride = (threshold > 2)? (env.chunk_size << 2) : (env.chunk_size << 1);
            int num_worker = (threshold == 1)? 1 : 2;
            if((threshold <= 2) && (task.tid & ((1 << (g->targs[0] - seg.middle - seg.chunk)))))
                return;
            if(threshold == 1 && env.rank > task.partner_using[0])
            {
                _thread_no_exec_MPI(task,2);
                return;
            }
            for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
            {
                _thread_read2_recv2(task,g,num_worker,env.chunk_size);
                update_offset(task,loop_stride);
            }
        }
        else if(isMiddle(g->targs[0]))
        {
            rank0 = env.rank ^ mpi_targ_mask_0;
            rank1 = env.rank ^ mpi_targ_mask_0;
            task.partner_using = {rank0,rank1};
            task.fd_using = {env.fd_arr[task.tid],env.fd_arr[task.tid]};
            loop_bound = (env.thread_state == env.qubit_offset[g->targs[0] + 1])? env.qubit_size[g->targs[0]] : env.thread_size;
            loop_stride = (env.thread_state == env.qubit_offset[g->targs[0] + 1])? (env.chunk_size << 1) : (env.qubit_size[g->targs[0] + 1] << 1);
            if(env.thread_state == env.qubit_offset[g->targs[0] + 1])
            {
                task.fd_offset_using = {0,env.qubit_size[g->targs[0]]};
                if(env.chunk_state == env.qubit_offset[g->targs[0]])
                {
                    if(env.rank > task.partner_using[0])
                        _thread_no_exec_MPI(task,2);
                    else
                        _thread_read2_recv2(task,g,1,0);
                }
                else
                {
                    for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
                    {   
                        _thread_read2_recv2(task,g,2,loop_stride >> 1);
                        update_offset(task,loop_stride);
                    }
                }
            }
            else
            {
                for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
                {
                    task.fd_offset_using = {cur_offset,cur_offset + env.qubit_size[g->targs[0]]};
                    for(long long j = 0;j < env.qubit_size[g->targs[0]];j += env.chunk_size)
                    {
                        _thread_read2_recv2(task,g,2,loop_stride >> 1);
                        update_offset(task,env.chunk_size);
                    }
                }
            }
        }
        else
        {
            rank0 = env.rank ^ mpi_targ_mask_0;
            task.partner_using = {rank0};
            task.fd_using = {env.fd_arr[task.tid]};
            task.fd_offset_using = {0};
            loop_stride = env.chunk_size << 1;
            if(env.thread_state == env.chunk_state && env.rank > task.partner_using[0])
                _thread_no_exec_MPI(task,1);
            else
            {
                for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
                {
                    _thread_read1_recv1(task,g);
                    update_offset(task,loop_stride);
                }
            }
        }
    }
}
void MPI_Runner::_thread_read1_recv1(thread_MPI_task &task,Gate * &g)
{
    int num_worker = (env.chunk_state == env.thread_state)? 1 : 2;
    int member_th = env.rank > task.partner_using[0]? 1 : 0;
    int partner_rank = task.partner_using[0];
    MPI_Irecv(&task.buffer[(!member_th) * env.chunk_state], env.chunk_state, MPI_DOUBLE_COMPLEX, partner_rank,task.tid, MPI_COMM_WORLD,&task.request[!member_th]);
    if(num_worker == 2)
    {
        if(pread(task.fd_using[0],&task.buffer[env.chunk_state << 1],env.chunk_size,task.fd_offset_using[0] + (!member_th) * env.chunk_size));
        MPI_Isend(&task.buffer[env.chunk_state << 1],env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request[2]);
    }
    if(pread(task.fd_using[0],&task.buffer[member_th * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + member_th * env.chunk_size));
    MPI_Wait(&task.request[!member_th],MPI_STATUS_IGNORE);
    g->run(task.buffer);
    MPI_Isend(&task.buffer[(!member_th) * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request[!member_th]);
    if(pwrite(task.fd_using[0],&task.buffer[member_th * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + member_th * env.chunk_size));
    if(num_worker == 2)
    {
        MPI_Recv(&task.buffer[env.chunk_state << 1],env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        if(pwrite(task.fd_using[0],&task.buffer[env.chunk_state << 1],env.chunk_size,task.fd_offset_using[0] + (!member_th) * env.chunk_size));
    }
}
void MPI_Runner::_thread_read2_recv2(thread_MPI_task &task,Gate * &g,int num_worker,long long stride)
{
    int member_th = env.rank > task.partner_using[0]? 1 : 0;//for rank
    int thread_member_th = isMpi(g->targs[1]) && isFile(g->targs[0]) && (task.tid & (1 << (g->targs[0] - seg.middle - seg.chunk)))? 1 : 0;//For all thread drive such as D F,otherwise it is 0
    MPI_Irecv(&task.buffer[((!member_th) << 1) * env.chunk_state], env.chunk_state, MPI_DOUBLE_COMPLEX, task.partner_using[0],task.tid, MPI_COMM_WORLD,&task.request[(!member_th) << 1]);
    MPI_Irecv(&task.buffer[(((!member_th) << 1) | 1) * env.chunk_state], env.chunk_state, MPI_DOUBLE_COMPLEX, task.partner_using[1],task.tid, MPI_COMM_WORLD,&task.request[((!member_th) << 1) | 1]);
    if(num_worker == 2)
    {
        if(pread(task.fd_using[0],&task.buffer[env.chunk_state << 2],env.chunk_size,task.fd_offset_using[0] + (!member_th) * stride + thread_member_th * 2 * stride));
        MPI_Isend(&task.buffer[env.chunk_state << 2],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,&task.request[4]);
        if(pread(task.fd_using[1],&task.buffer[5 * env.chunk_state],env.chunk_size,task.fd_offset_using[1] + (!member_th) * stride + thread_member_th * 2 * stride));
        MPI_Isend(&task.buffer[5 * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[1],task.tid,MPI_COMM_WORLD,&task.request[5]);
    }
    if(pread(task.fd_using[0],&task.buffer[(member_th << 1) * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + member_th * stride + thread_member_th * 2 * stride));
    if(pread(task.fd_using[1],&task.buffer[((member_th << 1) | 1) * env.chunk_state],env.chunk_size,task.fd_offset_using[1] + member_th * stride + thread_member_th * 2 * stride));
    MPI_Waitall(2,&task.request[(!member_th) << 1],MPI_STATUS_IGNORE);
    g->run(task.buffer);
    MPI_Isend(&task.buffer[((!member_th) << 1) * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,&task.request[(!member_th) << 1]);
    MPI_Isend(&task.buffer[(((!member_th) << 1) | 1) * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[1],task.tid,MPI_COMM_WORLD,&task.request[((!member_th) << 1) | 1]);
    if(num_worker == 2)
    {
        MPI_Irecv(&task.buffer[6 * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,&task.request[6]);
        MPI_Irecv(&task.buffer[7 * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[1],task.tid,MPI_COMM_WORLD,&task.request[7]);
    }
    if(pwrite(task.fd_using[0],&task.buffer[(member_th << 1) * env.chunk_state],env.chunk_size, task.fd_offset_using[0] + member_th * stride + thread_member_th * 2 * stride));
    if(pwrite(task.fd_using[1],&task.buffer[((member_th << 1) | 1) * env.chunk_state],env.chunk_size,task.fd_offset_using[1] + member_th * stride + thread_member_th * 2 * stride));

    if(num_worker == 2)
    {
        vector<bool>used(2,false);
        int next_buffer_idx = Get_Next_Undone_Buffer_index(task.request,used,2,6);
        while(next_buffer_idx != -1)
        {
            if(next_buffer_idx == 6 && pwrite(task.fd_using[0],&task.buffer[6 * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + (!member_th) * stride + thread_member_th * 2 * stride));
            if(next_buffer_idx == 7 && pwrite(task.fd_using[1],&task.buffer[7 * env.chunk_state],env.chunk_size,task.fd_offset_using[1] + (!member_th) * stride + thread_member_th * 2 * stride));
            next_buffer_idx = Get_Next_Undone_Buffer_index(task.request,used,2,6);
        }
    }
}
void MPI_Runner::_thread_read1_recv3(thread_MPI_task &task,Gate * &g,int num_worker)
{
    int member_th;
    for(int i = 0;i < 4;i++)
    {
        if(env.rank == task.partner_using[i])
            member_th = i;
        else
            MPI_Irecv(&task.buffer[i * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[i],task.tid,MPI_COMM_WORLD,&task.request[i]);
    }
    for(int i = 0;i < num_worker;i++)
    {
        if(i != member_th)
        {
            if(pread(task.fd_using[i],&task.buffer[(i + 4) * env.chunk_state],env.chunk_size,task.fd_offset_using[i]));
            MPI_Isend(&task.buffer[(i + 4) * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[i],task.tid,MPI_COMM_WORLD,&task.request[i + 4]);
        }
    }
    if(pread(task.fd_using[member_th],&task.buffer[member_th * env.chunk_state],env.chunk_size,task.fd_offset_using[member_th]));
    for(int i = 0;i < 4;i++)
    {
        if(i != member_th)
            MPI_Wait(&task.request[i],MPI_STATUS_IGNORE);
    }
    g->run(task.buffer);
    for(int i = 0;i < 4;i++)
        if(i != member_th) MPI_Isend(&task.buffer[i * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[i],task.tid,MPI_COMM_WORLD,&task.request[i]);
    if(pwrite(task.fd_using[member_th],&task.buffer[member_th * env.chunk_state],env.chunk_size,task.fd_offset_using[member_th]));
    unordered_map<int,int>m;
    for(int i = 0,rd_off = 8;i < num_worker;i++)
    {
        if(i != member_th)
        {
            MPI_Irecv(&task.buffer[rd_off * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[i],task.tid,MPI_COMM_WORLD,&task.request[rd_off]);
            m[rd_off++] = i;
        }
    }
    vector<bool>used(num_worker - 1,false);
    int next_buffer_idx = Get_Next_Undone_Buffer_index(task.request,used,num_worker - 1,8);
    while(next_buffer_idx != -1)
    {
        if(pwrite(task.fd_using[m[next_buffer_idx]],&task.buffer[next_buffer_idx * env.chunk_state],env.chunk_size,task.fd_offset_using[m[next_buffer_idx]]));
        next_buffer_idx = Get_Next_Undone_Buffer_index(task.request,used,num_worker - 1,8);
    }
}
void MPI_Runner::_thread_no_exec_MPI(thread_MPI_task &task,int round)
{
    vector<bool>used(round,false);
    for(int i = 0;i < round;i++)
    {
        MPI_Irecv(&task.buffer[(round + i) * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[i],task.tid,MPI_COMM_WORLD,&task.request[round + i]);
    }
    for(int i = 0;i < round;i++)
    {
        if(pread(task.fd_using[i], &task.buffer[i * env.chunk_state], env.chunk_size, task.fd_offset_using[i]));
        MPI_Isend(&task.buffer[i * env.chunk_state], env.chunk_state, MPI_DOUBLE_COMPLEX, task.partner_using[i], task.tid, MPI_COMM_WORLD,&task.request[i]);
    }
    int next_buffer_idx = Get_Next_Undone_Buffer_index(task.request,used,round,round);
    while(next_buffer_idx != -1)
    {
        if(pwrite(task.fd_using[next_buffer_idx - round], &task.buffer[next_buffer_idx * env.chunk_state], env.chunk_size, task.fd_offset_using[next_buffer_idx - round]));
        next_buffer_idx = Get_Next_Undone_Buffer_index(task.request,used,round,round);
    }
}
int MPI_Runner::Get_Next_Undone_Buffer_index(vector<MPI_Request> &request,vector<bool>&used,int round,int start)
{
    int res = -1;
    int fl;
    int done_cnt = 0;
    while(done_cnt != round)
    {
        done_cnt = 0;
        for(int i = 0;i < round;i++)
        {
            MPI_Test(&request[start + i],&fl,MPI_STATUS_IGNORE);
            if(fl)
            {
                if(!used[i])
                {
                    res = start + i;
                    used[i] = true;
                    return res;
                }
                else
                {
                    done_cnt++;
                }
            }
        }
    }
    return res;
}
void MPI_Runner::update_offset(thread_MPI_task &task,long long off)
{
    for(auto &x:task.fd_offset_using)
    {
        x+=off;
    }
}