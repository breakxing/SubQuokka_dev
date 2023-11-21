#ifndef IO_RUNNER_HPP
#define IO_RUNNER_HPP

#include "simulator.hpp"

class thread_IO_task
{
public:
    std::vector<std::complex<double>> buffer;
    std::vector<int> fd_table;
    std::vector<long long> fd_offset_table;
    
    std::vector<int> fd_using;
    std::vector<long long> fd_offset_using;

    // std::vector<std::vector<Gate>> subcircuits;
    thread_IO_task()=default;
    thread_IO_task(int);
};

class IO_Runner: public Simulator::circuitRunner {
    // using fd_table and fd_offset_table to set fd_using and fd_offset_using of each thread;
    std::vector<thread_IO_task> thread_tasks;
    void setFD(thread_IO_task &, Gate * &);
    void setFD_sub(thread_IO_task &);
    void inner_all_thread(thread_IO_task &,Gate * &,long long,int,bool);
    void all_thread_drive_scheduler(thread_IO_task &,Gate * &);
public:
    IO_Runner();
    void run(std::vector<Gate *> &) override;
    void run(std::vector<std::vector<Gate *>> &) override;
};

#endif
