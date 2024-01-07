#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP
#include <vector>
#include <iostream>
#include "circuit.hpp"

struct SEG{
    int N; // total_qubit
    int mpi;
    int file;
    int middle;
    int chunk;
};

struct ENV{
    int num_file;
    int num_thread;
    int half_num_thread;
    int quarter_num_thread;
    int eighth_num_thread;
    int rank;
    long long file_state;
    long long thread_state;
    long long chunk_state;

    long long file_size;
    long long thread_size;
    long long chunk_size;
    
    std::vector<long long> qubit_offset;
    std::vector<long long> qubit_size;
    std::vector<int> fd_arr;
    std::string runner_type;
    int is_subcircuit;
    int is_directIO;
    int is_MPI;
    int MPI_testing;
    int MPI_buffer_size;
    std::string dumpfile;
    // int SetOfSaveState;
};

extern struct SEG seg;
extern struct ENV env;

bool isMpi(int);
bool isFile(int);
bool isMiddle(int);
bool isChunk(int);

class Simulator
{
    std::vector<Gate *> circuit; // 全部的 gateMap
    std::vector<std::vector<Gate *>> subcircuits; // 由 gateMap 組成的 subcircuit

    void setupIni(std::string);
    void setupCir(std::string);

    Gate* setGate_IO(std::string &);
    void setupCircuit_IO(std::string);
    void setupSubCircuits_IO(std::string);
    
    Gate* setGate_DIO(std::string &);
    void setupCircuit_DIO(std::string);
    void setupSubCircuits_DIO(std::string);
    
    Gate* setGate_MEM(std::string &);
    void setupCircuit_MEM(std::string);
    void setupSubCircuits_MEM(std::string);

    void setupStateFile();
public:
    class circuitRunner;
    friend class circuitRunner;
    friend class IO_Runner;
    friend class MPI_Runner;
    circuitRunner *Runner;
    
    Simulator(std::string, std::string);
    virtual ~Simulator();
    void run();
    // run(Circuit);
};

class Simulator::circuitRunner {
public:
    // run the subcircuit with parameters in th Simulator;
    virtual void run(std::vector<Gate *> &) {
        std::cout << "[circuitRunner]: Not implemented" << std::endl;
        exit(1);
    }
    virtual void run(std::vector<std::vector<Gate *>> &) {
        std::cout << "[circuitRunner]: Not implemented" << std::endl;
        exit(1);
    }
    virtual ~circuitRunner(){};
};
#endif