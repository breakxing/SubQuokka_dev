#include <omp.h>
#include <unistd.h> //pread, pwrite

#include <complex>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

#include "circuit.hpp"

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
    for (int i = 0; i < seg.chunk; i++) {
        fd_table.push_back(tid);
        fd_offset_table.push_back(0);
    }

    for (int i = 0; i < seg.middle; i++) {
        fd_table.push_back(tid);
        fd_offset_table.push_back(env.qubit_size[i + seg.chunk]);
    }

    long long file_size = env.file_size;
    for (int i = 0; i < seg.file; i++) {
        int corr = findCorrespond(tid, i);

        if (corr > tid) {
            fd_table.push_back(corr - tid);
            fd_offset_table.push_back(0);
        }
    }
}
void MPI_Runner::_thread_read1_recv1_MPI(thread_MPI_task &task,Gate * &g,int member_th,int num_worker,long long stride,int partner_rank,long long cur_offset)
{
    MPI_Irecv(&task.buffer[(!member_th) * env.chunk_state], env.chunk_state, MPI_DOUBLE_COMPLEX, partner_rank,task.tid, MPI_COMM_WORLD,&task.request[!member_th]);
    if(num_worker == 2)
    {
        if(pread(env.fd_arr[task.tid],&task.buffer[env.chunk_state << 1],env.chunk_size,cur_offset + (!member_th) * stride));
        MPI_Isend(&task.buffer[env.chunk_state << 1],env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request[2]);
    }
    if(pread(env.fd_arr[task.tid],&task.buffer[member_th * env.chunk_state],env.chunk_size,cur_offset + member_th * stride));
    MPI_Wait(&task.request[!member_th],MPI_STATUS_IGNORE);
    g->run(task.buffer);
    MPI_Isend(&task.buffer[(!member_th) * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request[!member_th]);
    if(pwrite(env.fd_arr[task.tid],&task.buffer[member_th * env.chunk_state],env.chunk_size,cur_offset + member_th * stride));
    if(num_worker == 2)
    {
        MPI_Recv(&task.buffer[env.chunk_state << 1],env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        if(pwrite(env.fd_arr[task.tid],&task.buffer[env.chunk_state << 1],env.chunk_size,cur_offset + (!member_th) * stride));
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
        if (env.is_directIO) {
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
            
            if (file_count > 0 && skipThread(tid, targ)) {
#pragma omp barrier
                continue;
            }
            if(mpi_count == 0)
            {
                setFD(task, g);
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
            else if(mpi_count == 1)
            {
                int mpi_targ_offset = 1 << (g->targs[0] - seg.N);
                int member_0 = env.rank & (~mpi_targ_offset);
                int partner_rank = env.rank ^ mpi_targ_offset;
                if(targ.size() == 1)
                {
                    if(env.chunk_state == env.thread_state)
                    {
                        if(tid == 0)
                            cout<<"AAA\n";
                        if(env.rank == member_0)
                        {
                            _thread_read1_recv1_MPI(task,g,0,1,env.chunk_size,partner_rank,0);
                        }
                        else
                        {
                            if(pread(env.fd_arr[tid], &task.buffer[0], env.chunk_size, 0));
                            MPI_Isend(&task.buffer[0], env.chunk_state, MPI_DOUBLE_COMPLEX, partner_rank, tid, MPI_COMM_WORLD,&task.request[0]);
                            MPI_Recv(&task.buffer[env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,tid,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                            if(pwrite(env.fd_arr[tid], &task.buffer[env.chunk_state], env.chunk_size, 0));
                        }
                    }
                    else
                    {
                        if(tid == 0)
                            cout<<"BBB\n";
                        for(long long cur_offset = 0;cur_offset < env.thread_size;cur_offset += 2 * env.chunk_size)
                        {
                            if((env.rank & mpi_targ_offset) == 0)
                            {
                                _thread_read1_recv1_MPI(task,g,0,2,env.chunk_size,partner_rank,cur_offset);
                            }
                            else
                            {
                                _thread_read1_recv1_MPI(task,g,1,2,env.chunk_size,partner_rank,cur_offset);
                            }
                        }
                    }
                }
                else if(targ.size() == 2)
                {

                }
                else
                {

                }
            }
            else if(mpi_count == 2)
            {

            }
            else
                exit(-1);
            
#pragma omp barrier
        }
    }
}

void MPI_Runner::run(vector<vector<Gate *>> &subcircuits) {
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto task = thread_tasks[tid];

        static int count = 0;
        for (auto &subcircuit : subcircuits) {
            Gate *g = subcircuit[0];
            if (g->type == VSWAP) {
                vector<int> targ = subcircuit[0]->targs;  // increasing
                int file_count = subcircuit[0]->file_count;
                int middle_count = subcircuit[0]->middle_count;
                int chunk_count = subcircuit[0]->chunk_count;
                int nonchunk_count = subcircuit[0]->nonchunk_count;

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