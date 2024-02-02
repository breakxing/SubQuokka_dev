#ifndef MEM_RUNNER_HPP
#define MEM_RUNNER_HPP
#include <mpi.h>
#include "simulator.hpp"
void process_mem_usage(double& vm_usage, double& resident_set);
class thread_MEM_task
{
public:
    vector<complex<double>>buffer2;
    vector<complex<double>>buffer3;
    vector<complex<double>>buffer4;
    vector<int>partner_using;
    MPI_Request request1_send;
    MPI_Request request1_recv;
    MPI_Request request2_send;
    MPI_Request request2_recv;
    MPI_Request request3_send;
    MPI_Request request3_recv;
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
    void MPI_one_qubit_gate_diagonal(Gate* &);
    void MPI_two_qubit_gate_diagonal(Gate* &);
    void MPI_Swap_1_1(thread_MEM_task &,Gate* &g);
    void MPI_Swap_2_2(thread_MEM_task &,Gate* &g);
};


#endif
