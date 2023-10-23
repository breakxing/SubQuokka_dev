#ifndef MEM_RUNNER_HPP
#define MEM_RUNNER_HPP

#include "simulator.hpp"
void process_mem_usage(double& vm_usage, double& resident_set);

class MEM_Runner: public Simulator::circuitRunner {
public:
    std::vector<std::complex<double>> buffer;
    MEM_Runner();
    void run(std::vector<Gate *> &) override;
    void run(std::vector<std::vector<Gate *>> &) override;
};


#endif
