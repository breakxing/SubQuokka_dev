#ifndef DIO_RUNNER_HPP
#define DIO_RUNNER_HPP

#include "simulator.hpp"

class thread_DIO_task
{
public:
    std::complex<double> *buffer;
    std::vector<int> fd_table;
    std::vector<long long> fd_offset_table;
    
    std::vector<int> fd_using;
    std::vector<long long> fd_offset_using;

    // std::vector<std::vector<Gate>> subcircuits;
    thread_DIO_task()=default;
    thread_DIO_task(int);
};

class DIO_Runner: public Simulator::circuitRunner {
    // using fd_table and fd_offset_table to set fd_using and fd_offset_using of each thread;
    std::vector<thread_DIO_task> thread_tasks;
    void setFD(thread_DIO_task &, Gate * &);
    void setFD_sub(thread_DIO_task &);
public:
    DIO_Runner();
    ~DIO_Runner();
    void run(std::vector<Gate *> &) override;
    void run(std::vector<std::vector<Gate *>> &) override;
};

#endif
