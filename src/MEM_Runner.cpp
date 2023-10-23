#include <omp.h>
#include <unistd.h> //pread, pwrite

#include <complex>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>
#include <algorithm>
#include "circuit.hpp"

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

MEM_Runner::MEM_Runner() {
    //double vm, rss;
    //process_mem_usage(vm, rss);
    //cout << "[State Vector Before]VM: " << vm << "; RSS: " << rss << endl;
    buffer.resize((1ULL<<seg.N));
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

        long long  n_task = 1ULL << (seg.N - targ.size());
        int chunk_count = g->chunk_count;
        long long incre = env.chunk_state / (1 << chunk_count);
        #pragma omp parallel for schedule(static)
        for (long long task = 0; task < n_task; task += incre) {
            long long idx = bit_string(task, targ);
            g->_run(buffer, idx);
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