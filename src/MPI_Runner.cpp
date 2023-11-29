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
            if((g->name == "Z_Gate" || g->name == "Phase_Gate") && (ii * osz + jj == 0) && (!isChunk(targ[0])))\
                continue;\
            else if((g->name == "SWAP_Gate") && (ii * osz + jj == 0 || ii * osz + jj == 3) && (!isChunk(targ[0])))\
                continue;\
            else if(g->name == "CPhase_Gate" && (ii * osz + jj != 3) && (!isChunk(targ[1])))\
                continue;\
            if(pread(env.fd_arr[task.fd_using[ii]], &(task.buffer[(ii * osz + jj) * env.chunk_state]), env.chunk_size, task.fd_offset_using[jj]))\
                ;\
        }\
    }\
    g->run(task.buffer);\
    for (int ii = 0; ii < fsz; ii++) {\
        for (int jj = 0; jj < osz; jj++) {\
            if((g->name == "Z_Gate" || g->name == "Phase_Gate") && (ii * osz + jj == 0) && (!isChunk(targ[0])))\
                continue;\
            else if((g->name == "SWAP_Gate") && (ii * osz + jj == 0 || ii * osz + jj == 3) && (!isChunk(targ[0])))\
                continue;\
            else if(g->name == "CPhase_Gate" && (ii * osz + jj != 3) && (!isChunk(targ[1])))\
                continue;\
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
        thread_tasks[i].buffer.resize(16 * env.chunk_state);
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
            if(mpi_count == 0)
            {
                if(file_count != 0 && g->type != THREE_QUBIT)
                    all_thread_drive_scheduler(task,g);
                else
                {
                    if (!(file_count > 0 && skipThread(tid, targ))) 
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
            #pragma omp master
            MPI_Barrier(MPI_COMM_WORLD);
            #pragma omp barrier
        }
    #pragma omp barrier
    #pragma omp master
    MPI_Barrier(MPI_COMM_WORLD);
    #pragma omp barrier
    }
}

void MPI_Runner::run(vector<vector<Gate *>> &subcircuits) {
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto task = thread_tasks[tid];
        task.tid = tid;
        for (auto &subcircuit : subcircuits) {
            Gate *g = subcircuit[0];
            if (g->type == VSWAP) {
                vector<int> targ = subcircuit[0]->targs;  // increasing
                int file_count = subcircuit[0]->file_count;
                int middle_count = subcircuit[0]->middle_count;
                int chunk_count = subcircuit[0]->chunk_count;
                if(file_count != 0)
                {
                    if(g->name == "VSWAP_Gate_1_1")
                        all_thread_drive_scheduler(task,g);
                    else
                        all_thread_drive_vs2_2(task,g);
                }
                else
                {
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
            }
            else {
                if(subcircuit[0]->mpi_count)
                {
                    if(subcircuit[0]->name != "SWAP_Gate" || subcircuit.size() != 1)
                    {
                        if(task.tid == 0)
                            cout<<"ERROR"<<endl;
                        exit(-1);
                    }
                    else
                        MPI_Swap(task,g);
                }
                else
                {
                    setFD_sub(task);
                    outer_loop_m0(innerloop_sub)
                }
            }
            #pragma omp barrier
            #pragma omp master
            MPI_Barrier(MPI_COMM_WORLD);
            #pragma omp barrier
        }
        #pragma omp barrier
        #pragma omp master
        MPI_Barrier(MPI_COMM_WORLD);
        #pragma omp barrier
    }
}
void MPI_Runner::all_thread_drive_scheduler(thread_MPI_task &task,Gate * &g)
{
    vector<int> targ = g->targs;
    int tid = task.fd_table[0];
    long long func_loop_size;
    if(targ.size() == 1)
    {
        long long file_mask = 1 << (g->targs[0] - seg.middle - seg.chunk);
        int fd0 = tid & (~file_mask);
        int fd1 = tid | file_mask;
        if(env.thread_state == env.chunk_state && tid == fd1) return;
        task.fd_using = {env.fd_arr[fd0],env.fd_arr[fd1]};
        func_loop_size = (env.thread_state == env.chunk_state)? env.thread_size : env.thread_size >> 1;
        task.fd_offset_using = (tid == fd0)? vector<long long>{0,0} : vector<long long>{func_loop_size,func_loop_size};
        task.gate_buffer_using = {0,1};
        inner_all_thread(task,g,func_loop_size,2);
    }
    else if(targ.size() == 2)
    {
        long long file_mask_left = 1 << (g->targs[1] - seg.middle - seg.chunk);
        long long file_mask_right = 1 << (g->targs[0] - seg.middle - seg.chunk);
        if(isFile(targ[0]))
        {
            int fd0 = tid & (~file_mask_left) & (~file_mask_right);
            int fd1 = fd0 | file_mask_right;
            int fd2 = fd0 | file_mask_left;
            int fd3 = fd1 | file_mask_left;
            task.fd_using = {env.fd_arr[fd0],env.fd_arr[fd1],env.fd_arr[fd2],env.fd_arr[fd3]};
            int shift = min(int(log2(env.thread_state / env.chunk_state)),2);
            func_loop_size = env.thread_size >> shift;
            vector<long long>off0 = {0}; off0.resize(4,off0[0]);
            vector<long long>off1 = {func_loop_size}; off1.resize(4,off1[0]);
            vector<long long>off2 = {func_loop_size << 1};off2.resize(4,off2[0]);
            vector<long long>off3 = {func_loop_size * 3};off3.resize(4,off3[0]);
            unordered_map<int,vector<long long>>m = {{fd0,off0},{fd1,off1},{fd2,off2},{fd3,off3}};
            task.fd_offset_using = m[tid];
            if((shift == 0) && (tid != fd0)) return;
            else if((shift == 1) && (tid != fd0 && tid != fd1)) return;
            task.gate_buffer_using = {0,1,2,3};
            inner_all_thread(task,g,func_loop_size,4);
        }
        else if(isMiddle(targ[0]))
        {
            int fd0 = tid & (~file_mask_left);
            int fd1 = tid | file_mask_left;
            task.fd_using = {env.fd_arr[fd0],env.fd_arr[fd0],env.fd_arr[fd1],env.fd_arr[fd1]};
            long long base1 = 0;
            long long base2 = env.qubit_size[targ[0]];
            bool cond1 = env.thread_state == env.qubit_offset[targ[0] + 1];
            bool cond2 = env.chunk_state == env.qubit_offset[targ[0]];
            func_loop_size = (cond1 && !cond2)?env.qubit_size[targ[0] - 1] : env.qubit_size[targ[0]];
            if(cond1 && cond2 && tid == fd1) return;
            task.fd_offset_using = {base1,base2,base1,base2};
            if(cond1 && !cond2 && tid == fd1) task.fd_offset_using = {base1 + (env.qubit_size[targ[0] - 1]),base2 + (env.qubit_size[targ[0] - 1]),base1 + (env.qubit_size[targ[0] - 1]),base2 + (env.qubit_size[targ[0] - 1])};
            if(!cond1 && tid == fd1) task.fd_offset_using = {base1 + (env.thread_size >> 1),base2 + (env.thread_size >> 1),base1 + (env.thread_size >> 1),base2 + (env.thread_size >> 1)};
            task.gate_buffer_using = {0,1,2,3};
            if(cond1)
                inner_all_thread(task,g,func_loop_size,4);
            else
            {
                for(long long cur_offset = 0;cur_offset < (env.thread_size >> 1);cur_offset += env.qubit_size[targ[0] + 1])
                {
                    inner_all_thread(task,g,env.qubit_size[targ[0]],4);
                    task.fd_offset_using[0] += env.qubit_size[targ[0]];
                    task.fd_offset_using[1] += env.qubit_size[targ[0]];
                    task.fd_offset_using[2] += env.qubit_size[targ[0]];
                    task.fd_offset_using[3] += env.qubit_size[targ[0]];
                }
            }
        }
        else
        {
            int fd0 = tid & (~file_mask_left);
            int fd1 = tid | file_mask_left;
            if(env.thread_state == env.chunk_state && tid == fd1) return;
            task.fd_using = {env.fd_arr[fd0],env.fd_arr[fd1]};
            func_loop_size = (env.thread_state == env.chunk_state)? env.thread_size : env.thread_size >> 1;
            unordered_map<int,vector<long long>>m = {{fd0,{0,0}},{fd1,{func_loop_size,func_loop_size}}};
            task.fd_offset_using = m[tid];
            task.gate_buffer_using = {0,1};
            inner_all_thread(task,g,func_loop_size,2);
        }
    }
}
void MPI_Runner::all_thread_drive_vs2_2(thread_MPI_task &task,Gate * &g)
{
    vector<int> targ = g->targs;
    int tid = task.fd_table[0];
    long long func_loop_size;
    if(isFile(targ[2]))//L L D D
    {
        long long file_mask_left = 1 << (g->targs[3] - seg.middle - seg.chunk);
        long long file_mask_right = 1 << (g->targs[2] - seg.middle - seg.chunk);
        int fd0 = tid & (~file_mask_left) & (~file_mask_right);
        int fd1 = fd0 | file_mask_right;
        int fd2 = fd0 | file_mask_left;
        int fd3 = fd1 | file_mask_left;
        task.fd_using = {env.fd_arr[fd0],env.fd_arr[fd1],env.fd_arr[fd2],env.fd_arr[fd3]};
        int shift = min(int(log2(env.thread_state / env.chunk_state)),2);
        func_loop_size = env.thread_size >> shift;
        vector<long long>off0 = {0}; off0.resize(4,off0[0]);
        vector<long long>off1 = {func_loop_size}; off1.resize(4,off1[0]);
        vector<long long>off2 = {func_loop_size << 1};off2.resize(4,off2[0]);
        vector<long long>off3 = {func_loop_size * 3};off3.resize(4,off3[0]);
        unordered_map<int,vector<long long>>m = {{fd0,off0},{fd1,off1},{fd2,off2},{fd3,off3}};
        task.fd_offset_using = m[tid];
        if((shift == 0) && (tid != fd0)) return;
        else if((shift == 1) && (tid != fd0 && tid != fd1)) return;
        task.gate_buffer_using = {0,1,2,3};
        inner_all_thread(task,g,func_loop_size,4);
    }
    else//L L M D
    {
        long long file_mask = 1 << (g->targs[3] - seg.middle - seg.chunk);
        int fd0 = tid & (~file_mask);
        int fd1 = tid | file_mask;
        task.fd_using = {env.fd_arr[fd0],env.fd_arr[fd0],env.fd_arr[fd1],env.fd_arr[fd1]};
        long long base1 = 0;
        long long base2 = env.qubit_size[targ[2]];
        bool cond1 = env.thread_state == env.qubit_offset[targ[2] + 1];
        bool cond2 = env.chunk_state == env.qubit_offset[targ[2]];
        func_loop_size = (cond1 && !cond2)?env.qubit_size[targ[2] - 1] : env.qubit_size[targ[2]];
        if(cond1 && cond2 && tid == fd1) return;
        task.fd_offset_using = {base1,base2,base1,base2};
        if(cond1 && !cond2 && tid == fd1) task.fd_offset_using = {base1 + (env.qubit_size[targ[2] - 1]),base2 + (env.qubit_size[targ[2] - 1]),base1 + (env.qubit_size[targ[2] - 1]),base2 + (env.qubit_size[targ[2] - 1])};
        if(!cond1 && tid == fd1) task.fd_offset_using = {base1 + (env.thread_size >> 1),base2 + (env.thread_size >> 1),base1 + (env.thread_size >> 1),base2 + (env.thread_size >> 1)};
        task.gate_buffer_using = {0,1,2,3};
        if(cond1)
            inner_all_thread(task,g,func_loop_size,4);
        else
        {
            for(long long cur_offset = 0;cur_offset < (env.thread_size >> 1);cur_offset += env.qubit_size[targ[2] + 1])
            {
                inner_all_thread(task,g,env.qubit_size[targ[2]],4);
                task.fd_offset_using[0] += env.qubit_size[targ[2]];
                task.fd_offset_using[1] += env.qubit_size[targ[2]];
                task.fd_offset_using[2] += env.qubit_size[targ[2]];
                task.fd_offset_using[3] += env.qubit_size[targ[2]];
            }
        }
    }
}
void MPI_Runner::inner_all_thread(thread_MPI_task &task,Gate * &g,long long func_loop_size,int round)
{
    long long stride = env.chunk_size;
    if(g->name == "CPhase_Gate" && isMpi(g->targs[1]) && isFile(g->targs[0]))
        stride = env.chunk_size << 1;
    for(long long i = 0;i < func_loop_size;i += stride)
    {
        for(int j = 0;j < round;j++)
        {
            if((g->name == "Z_Gate" || g->name == "Phase_Gate") && (j == 0) && (!isChunk(g->targs[0]) && (!isMpi(g->targs[0])))) continue;
            else if((g->name == "SWAP_Gate") && (j == 0 || j == 3) && (!isChunk(g->targs[0]))) continue;
            else if(g->name == "CPhase_Gate" && (j != 3) && (!isChunk(g->targs[1]) && (!isMpi(g->targs[1])))) continue;
            if(pread(task.fd_using[j],&task.buffer[task.gate_buffer_using[j] * env.chunk_state],env.chunk_size,task.fd_offset_using[j]));
        }
        g->run(task.buffer);
        for(int j = 0;j < round;j++)
        {
            if((g->name == "Z_Gate" || g->name == "Phase_Gate") && (j == 0) && (!isChunk(g->targs[0])) && (!isMpi(g->targs[0]))) continue;
            else if((g->name == "SWAP_Gate") && (j == 0 || j == 3) && (!isChunk(g->targs[0]))) continue;
            else if(g->name == "CPhase_Gate" && (j != 3) && (!isChunk(g->targs[1]) && (!isMpi(g->targs[1])))) continue;
            if(pwrite(task.fd_using[j],&task.buffer[task.gate_buffer_using[j] * env.chunk_state],env.chunk_size,task.fd_offset_using[j]));
            task.fd_offset_using[j] += stride;
        }
    }
}
void MPI_Runner::update_offset(thread_MPI_task &task,long long &off)
{
    for(auto &x:task.fd_offset_using)
        x+=off;
}