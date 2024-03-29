#include <omp.h>
#include <unistd.h> //pread, pwrite

#include <complex>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>
#include <algorithm>
#include "circuit.hpp"
#include <algorithm>
#include "MEM_Runner.hpp"
using namespace std;

void process_mem_usage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}
thread_MEM_task::thread_MEM_task(int tid,int MPI_buffer_size)
{
    buffer2.resize(MPI_buffer_size * env.chunk_state);
    buffer3.resize(MPI_buffer_size * env.chunk_state);
    buffer4.resize(MPI_buffer_size * env.chunk_state);
}
MEM_Runner::MEM_Runner() {
    if(env.is_MPI)
    {
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
        for(int i = 0;i < env.num_thread;i++)
        {
            thread_tasks.push_back(thread_MEM_task(i,env.MPI_buffer_size));
        }
    }
    //double vm, rss;
    //process_mem_usage(vm, rss);
    //cout << "[State Vector Before]VM: " << vm << "; RSS: " << rss << endl;
    buffer.resize((1ULL<<seg.N));
    if(env.rank == 0)
    buffer[0] = complex<double>(1, 0);

    //process_mem_usage(vm, rss);
    //cout << "[State Vector After] VM: " << vm << "; RSS: " << rss << endl;
    
    // buffer.resize((1ULL<<seg.N));
    // buffer[0] = complex<double>(1, 0);
}

inline long long bit_string(long long &task, const vector<int> &targ){
    long long res = task;
    for (auto &x : targ)
    {
        long long mask = 1ULL << x;
        mask -= 1;
        res = ((res >> x) << (x+1)) | (res & mask);
    }
    return res;    
}

void MEM_Runner::run(vector<Gate *> &circuit) {
    for (auto &g : circuit) {
        vector<int> targ = g->targs;  // increasing
        if(g->mpi_count)
        {
            if(g->name == "Z_Gate_MEM" || g->name == "Phase_Gate_MEM" || g->name == "RZ_Gate_MEM") MPI_one_qubit_gate_diagonal(g);
            else if(g->name == "CPhase_Gate_MEM" || g->name == "RZZ_Gate_MEM") MPI_two_qubit_gate_diagonal(g);
            else
            {
                #pragma omp parallel
                {
                    int tid = omp_get_thread_num();
                    auto task = thread_tasks[tid];
                    task.tid = tid;
                    if(g->name == "SWAP_Gate_MEM") MPI_Swap_1_1(task,g);
                    else MPI_gate_scheduler(task,g);
                }
            }
        }
        else
        {
            long long  n_task = 1ULL << (seg.N - targ.size());
            int chunk_count = g->chunk_count;
            long long incre = env.chunk_state / (1 << chunk_count);
            #pragma omp parallel for schedule(static)
            for (long long task = 0; task < n_task; task += incre) {
                long long idx = bit_string(task, targ);
                g->_run(buffer, idx);
            }
        }
    }
    // ed = omp_get_wtime();
    // cout << ed-st << "s" << endl;
}

void MEM_Runner::run(vector<vector<Gate *>> &subcircuits) {
    
    // int i = 0;

    for (auto &subcircuit : subcircuits) {
        Gate g = *subcircuit[0];
        vector<int> targ = subcircuit[0]->targs;  // increasing
        
        // printf("%d\n", i++);
        
        if (g.type == VSWAP) {
            if(subcircuit[0]->mpi_count)
            {
                #pragma omp parallel
                {
                    int tid = omp_get_thread_num();
                    auto task = thread_tasks[tid];
                    task.tid = tid;
                    if(subcircuit[0]->mpi_count == 1) MPI_Swap_1_1(task,subcircuit[0]);
                    else MPI_Swap_2_2(task,subcircuit[0]);
                }
            }
            else
            {
                long long  n_task = 1ULL << (seg.N - targ.size());
                int chunk_count = subcircuit[0]->chunk_count;
                long long incre = env.chunk_state / (1 << chunk_count);
                #pragma omp parallel for schedule(static)
                for (long long task = 0; task < n_task; task += incre) {
                    long long idx = bit_string(task, targ);
                    for(auto &g : subcircuit) {
                        g->_run(buffer, idx);
                    }
                }
            }
        }
        else{
            long long  n_task = 1ULL << seg.N;
            long long incre = env.chunk_state;
            // #pragma omp parallel for schedule(guided, 1)
            // #pragma omp parallel for schedule(dynamic)
            #pragma omp parallel for schedule(static)
            for (long long idx = 0; idx < n_task; idx += incre) {
                for(auto &g : subcircuit) {
                    g->_run(buffer, idx);
                }
            }
        }
    }
}