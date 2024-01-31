#ifndef MEM_RUNNER_HPP
#define MEM_RUNNER_HPP
#include <mpi.h>
#include "simulator.hpp"
void process_mem_usage(double& vm_usage, double& resident_set);
class thread_MEM_task
{
public:
    vector<complex<double>>buffer2;
    vector<int>partner_using;
    vector<MPI_Request>request_send;
    vector<MPI_Request>request_recv;
    int tid;
    thread_MEM_task()=default;
    thread_MEM_task(int,int);
};
class MEM_Runner: public Simulator::circuitRunner {
public:
    vector<thread_MEM_task> thread_tasks;
    std::vector<std::complex<double>> buffer;
    MEM_Runner();
    void run(std::vector<Gate *> &) override;
    void run(std::vector<std::vector<Gate *>> &) override;
    void MPI_gate_scheduler(thread_MEM_task &,Gate* &g);
    void _mpi_one_gate_inner(thread_MEM_task &,Gate* &g,long long);
};


#endif
