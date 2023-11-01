#ifndef MPI_Runner_HPP
#define MPI_Runner_HPP

#include "simulator.hpp"
#include <mpi.h>
class thread_MPI_task
{
public:
    std::vector<std::complex<double>> buffer;
    std::vector<int> fd_table;
    std::vector<long long> fd_offset_table;
    int tid;
    vector<MPI_Request> request;
    int request_size;
    std::vector<int> fd_using;
    std::vector<long long> fd_offset_using;
    std::vector<int>partner_using;
    // std::vector<std::vector<Gate>> subcircuits;
    thread_MPI_task()=default;
    thread_MPI_task(int);
};
class MPI_Runner: public Simulator::circuitRunner {
    // using fd_table and fd_offset_table to set fd_using and fd_offset_using of each thread;
    std::vector<thread_MPI_task> thread_tasks;
    void setFD(thread_MPI_task &, Gate * &);
    void setFD_sub(thread_MPI_task &);
    void _thread_read1_recv1(thread_MPI_task &,Gate * &);
    void MPI_gate_scheduler(thread_MPI_task &,Gate * &);
    void _thread_read2_recv2(thread_MPI_task &,Gate * &,int,long long);
    void _thread_read1_recv3(thread_MPI_task &,Gate * &,int);
    void _thread_pure_send_recv_MPI(thread_MPI_task &,int round);
    int Get_Next_Undone_Buffer_index(vector<MPI_Request> &,vector<bool>&,int ,int);
    void update_offset(thread_MPI_task &,long long);

public:
    MPI_Runner();
    void run(std::vector<Gate *> &) override;
    void run(std::vector<std::vector<Gate *>> &) override;
};

#endif