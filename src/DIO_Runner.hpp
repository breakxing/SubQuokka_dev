#ifndef DIO_RUNNER_HPP
#define DIO_RUNNER_HPP

#include "simulator.hpp"
#include <mpi.h>
class thread_DIO_task
{
public:
    complex<double> *buffer1;
    complex<double> *buffer2;
    complex<double> *buffer3;
    complex<double> *buffer4;
    std::vector<int> fd_table;
    std::vector<long long> fd_offset_table;
    int tid;
    vector<MPI_Request> request;
    int request_size;
    std::vector<int> fd_using;
    std::vector<long long> fd_offset_using;
    std::vector<int>partner_using;
    std::vector<int>gate_buffer_using;
    bool has_non_blocking;
    // std::vector<std::vector<Gate>> subcircuits;
    thread_DIO_task()=default;
    thread_DIO_task(int,int);
};

class DIO_Runner: public Simulator::circuitRunner {
    // using fd_table and fd_offset_table to set fd_using and fd_offset_using of each thread;
    std::vector<thread_DIO_task> thread_tasks;
    void setFD(thread_DIO_task &, Gate * &);
    void setFD_sub(thread_DIO_task &);
    void MPI_Swap(thread_DIO_task &,Gate * &);
    void _thread_MPI_swap(thread_DIO_task &,Gate * &,bool &);
    void MPI_gate_scheduler(thread_DIO_task &,Gate * &);
    void _two_gate_mpi_read1_recv1(thread_DIO_task &,Gate * &);
    void _mpi_one_gate_inner(thread_DIO_task &,Gate * &);
    void _thread_read2_recv2(thread_DIO_task &,Gate * &,int,long long);
    void _thread_read1_recv3(thread_DIO_task &,Gate * &,int);
    void _thread_no_exec_MPI(thread_DIO_task &,int round);
    int Get_Next_Undone_Buffer_index(vector<MPI_Request> &,vector<bool>&,int ,int);
    void update_offset(thread_DIO_task &,long long&);
    void inner_all_thread(thread_DIO_task &,Gate * &,long long,int);
    void all_thread_drive_scheduler(thread_DIO_task &,Gate * &);
    void all_thread_drive_vs2_2(thread_DIO_task &,Gate * &);

    bool skip_read_write(Gate * &,const int &);
    void MPI_one_qubit_gate_diagonal(thread_DIO_task &,Gate * &);
    void MPI_two_qubit_gate_diagonal(thread_DIO_task &,Gate * &);
    void MPI_special_gate_inner(thread_DIO_task &,Gate * &,long long,int);

    void MPI_vs2_2(thread_DIO_task &,Gate * &);
public:
    DIO_Runner();
    ~DIO_Runner();
    void run(std::vector<Gate *> &) override;
    void run(std::vector<std::vector<Gate *>> &) override;
    int MPI_buffer_size;
};

#endif
