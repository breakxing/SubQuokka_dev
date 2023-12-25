#ifndef MPI_Runner_HPP
#define MPI_Runner_HPP

#include "simulator.hpp"
#include <mpi.h>
class thread_MPI_task
{
public:
    std::vector<std::complex<double>> buffer1;
    vector<complex<double>>buffer2;
    std::vector<int> fd_table;
    std::vector<long long> fd_offset_table;
    int tid;
    vector<MPI_Request> request;
    int request_size;
    std::vector<int> fd_using;
    std::vector<long long> fd_offset_using;
    std::vector<int>partner_using;
    std::vector<int>gate_buffer_using;
    // std::vector<std::vector<Gate>> subcircuits;
    thread_MPI_task()=default;
    thread_MPI_task(int,int);
};
class MPI_Runner: public Simulator::circuitRunner {
    // using fd_table and fd_offset_table to set fd_using and fd_offset_using of each thread;
    std::vector<thread_MPI_task> thread_tasks;
    void setFD(thread_MPI_task &, Gate * &);
    void setFD_sub(thread_MPI_task &);
    void MPI_Swap(thread_MPI_task &,Gate * &);
    void _thread_MPI_swap(thread_MPI_task &,Gate * &,bool &);
    void MPI_gate_scheduler(thread_MPI_task &,Gate * &);
    void _two_gate_mpi_read1_recv1(thread_MPI_task &,Gate * &);
    void _mpi_one_gate_inner(thread_MPI_task &,Gate * &);
    void _thread_read2_recv2(thread_MPI_task &,Gate * &,int,long long);
    void _thread_read1_recv3(thread_MPI_task &,Gate * &,int);
    void _thread_no_exec_MPI(thread_MPI_task &,int round);
    int Get_Next_Undone_Buffer_index(vector<MPI_Request> &,vector<bool>&,int ,int);
    void update_offset(thread_MPI_task &,long long&);
    void inner_all_thread(thread_MPI_task &,Gate * &,long long,int);
    void all_thread_drive_scheduler(thread_MPI_task &,Gate * &);
    void all_thread_drive_vs2_2(thread_MPI_task &,Gate * &);

    bool skip_read_write(Gate * &,const int &);
    void MPI_one_qubit_gate_diagonal(thread_MPI_task &,Gate * &);
    void MPI_two_qubit_gate_diagonal(thread_MPI_task &,Gate * &);
    void MPI_special_gate_inner(thread_MPI_task &,Gate * &,long long,int);

    void MPI_vs2_2(thread_MPI_task &,Gate * &);
    void aaa(vector<complex<double>> &buffer1,vector<complex<double>> &buffer2,vector<complex<double>> &buffer3,vector<complex<double>> &buffer4,vector<int>&);

public:
    MPI_Runner();
    void run(std::vector<Gate *> &) override;
    void run(std::vector<std::vector<Gate *>> &) override;
    int MPI_buffer_size;
};

#endif
