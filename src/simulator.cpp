#include <dirent.h>
#include <fcntl.h>
#include <libgen.h>
#include <omp.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cassert>
#include <complex>
#include <fstream>
#include <functional>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "INIReader.h"
#include "circuit.hpp"

#include "simulator.hpp"

#include "IO_Runner.hpp"
#include "DIO_Runner.hpp"
#include "MEM_Runner.hpp"
#include "MPI_Runner.hpp"

using namespace std;

struct SEG seg;
struct ENV env;

// check qubit is in mpi seg;
bool isMpi(int target) {
    return target >= seg.N;
}

// check qubit is in file seg;
bool isFile(int target) {
    return (target >= (seg.middle + seg.chunk)) && (target < seg.N);
}

// check qubit is in middle seg;
bool isMiddle(int target) {
    return (target < (seg.middle + seg.chunk)) && (target >= seg.chunk);
}

// check qubit is in chunk seg;
bool isChunk(int target) {
    return target < seg.chunk && target >= 0;
}

inline int make_dir(string dir) {
    DIR *mydir = NULL;
    if (!(mydir = opendir(dir.c_str()))) {  // check is dir or not
        if (mkdir(dir.c_str(), (S_IRWXU | S_IRWXG | S_IRWXO))) {
            cerr << "[DIR]: " << dir << " create failed!" << endl;
            exit(1);
        }
        // cout << "[DIR]: " << dir << " created sucess!" << endl;
    } else {
        // cout << "[DIR]: " << dir << " exist!" << endl;
        closedir(mydir);
    }
    return 0;
}

inline int file_exists(string state_path) {
    struct stat buffer;
    return (stat(state_path.c_str(), &buffer) == 0);
}

inline void create_file(string state_path) {
    int fd;
    if (env.runner_type == "DirectIO") {
        fd = open(state_path.c_str(), O_RDWR | O_CREAT | O_DIRECT, 0666);
    }
    else {
        fd = open(state_path.c_str(), O_RDWR | O_CREAT, 0666);
    }
    if (fd < 0) {
        cerr << "[FILE]: " << state_path << " create failed!, fd: " << fd << endl;
        exit(1);
    } else {
        env.fd_arr.push_back(fd);
        // cout << "[FILE]: " << state_path << " create success!, fd: " << fd << endl;
    }
}

void setSEG(INIReader &reader) {
    string section("system");
    seg.N = reader.GetInteger(section, "total_qbit", 0);
    seg.mpi = reader.GetInteger(section, "mpi_qbit", 0);
    seg.N = seg.N - seg.mpi;
    seg.file = reader.GetInteger(section, "file_qbit", 0);
    seg.chunk = reader.GetInteger(section, "chunk_qbit", 0);
    seg.middle = seg.N - seg.file - seg.chunk;

    // assert for invalid case
    
    assert(seg.N >= (seg.file + seg.chunk));
}

inline void setupFiles(string state_paths) {
    stringstream ss(state_paths);
    vector<string> tokens;
    string token;
    while (getline(ss, token, ',')) {
        // cout << token << endl; // print path of state files
        tokens.push_back(token);
    }

    for (int i = 0; i < env.num_file; i++) {
        string state_path = tokens[i];
        // cout << state_path << endl;
        int found = state_path.find_last_of("/\\");
        string state_dir = state_path.substr(0, found);
        if (!make_dir(state_dir)) {
            if (file_exists(state_path)) {
                remove(state_path.c_str());
            }
            create_file(state_path);
        }
    }
}

void setENV(INIReader &reader) {
    string section("system");

    // max_depth = reader.GetInteger(section, "max_depth", 1000);
    // IsDensity = reader.GetInteger(section, "is_density", 0);
    // SkipInithread_state = reader.GetInteger(section, "skip_init_state", 0);
    // SetOfSaveState = reader.GetInteger(section, "set_of_save_state", 1);
    env.dumpfile = reader.Get(section, "dump_file", "");
    env.runner_type = reader.Get(section, "runner_type", "IO");
    env.is_subcircuit = reader.GetInteger(section, "is_subcircuit", 1);
    env.is_MPI = (seg.mpi > 0);
    env.num_file = (1ULL << seg.file);
    env.num_thread = env.num_file;
    env.half_num_thread = env.num_thread / 2;
    env.quarter_num_thread = env.num_thread / 4;
    env.eighth_num_thread = env.num_thread / 8;

    env.file_state = (1ULL << (seg.N - seg.file));
    env.thread_state = env.file_state;
    env.chunk_state = (1ULL << seg.chunk);

    env.file_size = env.file_state * sizeof(complex<double>);
    env.thread_size = env.thread_state * sizeof(complex<double>);
    env.chunk_size = env.chunk_state * sizeof(complex<double>);

    long long temp_offset = 1;
    long long temp_size = sizeof(complex<double>);
    for (int i = 0; i <= seg.N; i++) {
        env.qubit_offset.push_back(temp_offset);
        env.qubit_size.push_back(temp_size);
        temp_offset *= 2;
        temp_size *= 2;
    }

    // printf("is density: %d\n", IsDensity);
    // if(env.runner_type != "MEM") {
    //     string state_paths;
    //     state_paths = reader.GetString(section, "state_paths", "");
    //     if (state_paths.empty()) {
    //         cerr << "[Config File]: state_paths not found." << endl;
    //         exit(1);
    //     }
    //     setupFiles(state_paths);
    // }
}

// initialize the seg and env from the *.ini file
void Simulator::setupIni(string ini) {
    INIReader reader(ini);
    if (reader.ParseError() < 0) {
        cerr << "Error parsing config.ini!" << endl;
        exit(1);
    }
    
    setSEG(reader);
    setENV(reader);
    //double vm, rss;
    //process_mem_usage(vm, rss);
    //cout << "[Setup SEG and ENV]VM: " << vm << "; RSS: " << rss << endl;
    omp_set_num_threads(env.num_thread);

    if (env.runner_type == "IO") {
        std::cout<<"IO Mode\n";
        Runner = new IO_Runner();
    } else if (env.runner_type == "DirectIO") {
        std::cout<<"DIO Mode\n";
        Runner = new DIO_Runner();
    } else if (env.runner_type == "MEM") {
        std::cout<<"MEM Mode\n";
        Runner = new MEM_Runner();
    } else if (env.runner_type == "RDMA") {
        cerr << "[Config File]: RDMA runner not found." << endl;
        exit(1);
    } else if (env.runner_type == "MPI") {
        std::cout<<"MPI Mode\n";
        Runner = new MPI_Runner();
        string section("system");
        string state_paths;
        state_paths = reader.GetString(section, "state_paths", "");
        std::istringstream ss(state_paths);
        string new_paths = "";
        string token;
        while(getline(ss,token,','))
        {
            token.insert(7,to_string(env.rank));
            new_paths += token + ",";
        }
        if (state_paths.empty()) {
            cerr << "[Config File]: state_paths not found." << endl;
            exit(1);
        }
        new_paths.pop_back();
        setupFiles(new_paths);
    } else if (env.runner_type == "GPU") {
        cerr << "[Config File]: GPU runner not found." << endl;
        exit(1);
    }
    return;
}

void Simulator::setupCir(string cir){
    if (env.runner_type == "IO") {
        if(env.is_subcircuit)
            setupSubCircuits_IO(cir);
        else
            setupCircuit_IO(cir);
    } else if (env.runner_type == "DirectIO") {
        if(env.is_subcircuit)
            setupSubCircuits_DIO(cir);
        else {
            setupCircuit_DIO(cir);
        }
    } else if (env.runner_type == "MEM") {
        if(env.is_subcircuit)
            setupSubCircuits_MEM(cir);
        else {
            setupCircuit_MEM(cir);
        }
    } else if (env.runner_type == "RDMA") {
        cerr << "[Config File]: RDMA runner not found." << endl;
        exit(1);
    } else if (env.runner_type == "MPI") {
        if(env.is_subcircuit)
        {
            cerr << "Not implement yet" << endl;
            exit(1);
        }
        else
            setupCircuit_MPI(cir);
    } else if (env.runner_type == "GPU") {
        cerr << "[Config File]: GPU runner not found." << endl;
        exit(1);
    }
}

template<class T>
T *setGate_1(stringstream &ss){
    int target;
    ss >> target;
    vector<int> targ{target};
    return new T(targ);
}

template<class T>
T *setGate_2(stringstream &ss){
    int target0;
    int target1;
    ss >> target0 >> target1;
    vector<int> targ{target0, target1};
    return new T(targ);
}


template<class T>
T *setGate_3(stringstream &ss){
    int target0;
    int target1;
    int target2;
    ss >> target0 >> target1 >> target2;
    vector<int> targ{target0, target1, target2};
    return new T(targ);
}

template<class T>
T *setGate_vswap(stringstream &ss, int swap_size){
    vector<int> targ;
    for (int i = 0; i < swap_size*2; i++) {
        int temp;
        ss >> temp;
        targ.push_back(temp);
    }
    return new T(targ);
}

template <class T>
T *setGate_Phase(stringstream &ss, int n_qubits) {
    int target;
    vector<int> targ;
    for (int i = 0; i < n_qubits; i++) {
        ss >> target;
        targ.push_back(target);
    }
    double phi;
    ss >> phi;
    return new T(targ, phi);
}

template <class T>
Gate *setUnitary(stringstream &ss, int n_qubits) {
    int target;
    vector<int> targ;
    for (int i = 0; i < n_qubits; i++) {
        ss >> target;
        targ.push_back(target);
    }

    int n_coeff = 1ULL << (2*n_qubits);
    double re;
    double im;
    vector<complex<double>> coeff(n_coeff);
    for (int i = 0; i < n_coeff; i++) {
        ss >> re >> im;
        coeff[i] = complex<double>(re, im);
    }

    return new T(targ, coeff);
}

Gate *Simulator::setGate_MPI(string &line) {
    // static int count = 0;
    // cerr << "[setGate] " << count++ << "th gate." << endl;
    stringstream ss;
    ss << line;
    string gate_ops;
    ss >> gate_ops;
    if (gate_ops == "U1") {
        return setUnitary<U1_Gate>(ss, ONE_QUBIT);
    } else if (gate_ops == "U2") {
        return setUnitary<U2_Gate>(ss, TWO_QUBIT);
    } else if (gate_ops == "U3") {
        return setUnitary<U3_Gate>(ss, THREE_QUBIT);
    } else if (gate_ops == "H"){
        return setGate_1<H_Gate>(ss);
    // } else if (gate_ops == "S"){
    //     return setGate1<S_Gate>(ss);
    // } else if (gate_ops == "T"){
    //     return setGate1<T_Gate>(ss);
    } else if (gate_ops == "X"){
        return setGate_1<X_Gate>(ss);
    } else if (gate_ops == "Y"){
        return setGate_1<Y_Gate>(ss);
    } else if (gate_ops == "Z"){
        return setGate_1<Z_Gate>(ss);
    } else if (gate_ops == "P"){
        return setGate_Phase<Phase_Gate>(ss, ONE_QUBIT);
    } else if (gate_ops == "RX"){
        return setGate_Phase<RX_Gate>(ss, ONE_QUBIT);
    } else if (gate_ops == "RY"){
        return setGate_Phase<RY_Gate>(ss, ONE_QUBIT);
    } else if (gate_ops == "RZ"){
        return setGate_Phase<RZ_Gate>(ss, ONE_QUBIT);
    } else if (gate_ops == "SWAP") {
        return setGate_2<SWAP_Gate>(ss);
    } else if (gate_ops == "CP"){
        return setGate_Phase<CPhase_Gate>(ss, TWO_QUBIT);
    } else if (gate_ops == "RZZ"){
        return setGate_Phase<RZZ_Gate>(ss, TWO_QUBIT);
    } /*else if (gate_ops == "VSWAP_1_1") {
        return setGate_vswap<VSWAP_Gate_1_1>(ss, 1);
    } else if (gate_ops == "VSWAP_2_2") {
        return setGate_vswap<VSWAP_Gate_2_2>(ss, 2);
    } else if (gate_ops == "VSWAP_3_3") {
        return setGate_vswap<VSWAP_Gate_3_3>(ss, 3);
    } else if (gate_ops == "VSWAP_4_4") {
        return setGate_vswap<VSWAP_Gate_4_4>(ss, 4);
    } else if (gate_ops == "VSWAP_6_6") {
        return setGate_vswap<VSWAP_Gate_6_6>(ss, 6);
    }*/ else {
        cerr << "[setGate]: Not implemented yet." << endl;
        exit(1);
    }
}


// create a Gate by a string description
Gate *Simulator::setGate_IO(string &line) {
    // static int count = 0;
    // cerr << "[setGate] " << count++ << "th gate." << endl;
    stringstream ss;
    ss << line;
    string gate_ops;
    ss >> gate_ops;
    if (gate_ops == "U1") {
        return setUnitary<U1_Gate>(ss, ONE_QUBIT);
    } else if (gate_ops == "U2") {
        return setUnitary<U2_Gate>(ss, TWO_QUBIT);
    } else if (gate_ops == "U3") {
        return setUnitary<U3_Gate>(ss, THREE_QUBIT);
    } else if (gate_ops == "H"){
        return setGate_1<H_Gate>(ss);
    // } else if (gate_ops == "S"){
    //     return setGate1<S_Gate>(ss);
    // } else if (gate_ops == "T"){
    //     return setGate1<T_Gate>(ss);
    } else if (gate_ops == "X"){
        return setGate_1<X_Gate>(ss);
    } else if (gate_ops == "Y"){
        return setGate_1<Y_Gate>(ss);
    } else if (gate_ops == "Z"){
        return setGate_1<Z_Gate>(ss);
    } else if (gate_ops == "P"){
        return setGate_Phase<Phase_Gate>(ss, ONE_QUBIT);
    } else if (gate_ops == "RX"){
        return setGate_Phase<RX_Gate>(ss, ONE_QUBIT);
    } else if (gate_ops == "RY"){
        return setGate_Phase<RY_Gate>(ss, ONE_QUBIT);
    } else if (gate_ops == "RZ"){
        return setGate_Phase<RZ_Gate>(ss, ONE_QUBIT);
    } else if (gate_ops == "SWAP") {
        return setGate_2<SWAP_Gate>(ss);
    } else if (gate_ops == "CP"){
        return setGate_Phase<CPhase_Gate>(ss, TWO_QUBIT);
    } else if (gate_ops == "RZZ"){
        return setGate_Phase<RZZ_Gate>(ss, TWO_QUBIT);
    } else if (gate_ops == "VSWAP_1_1") {
        return setGate_vswap<VSWAP_Gate_1_1>(ss, 1);
    } else if (gate_ops == "VSWAP_2_2") {
        return setGate_vswap<VSWAP_Gate_2_2>(ss, 2);
    } else if (gate_ops == "VSWAP_3_3") {
        return setGate_vswap<VSWAP_Gate_3_3>(ss, 3);
    } else if (gate_ops == "VSWAP_4_4") {
        return setGate_vswap<VSWAP_Gate_4_4>(ss, 4);
    } else if (gate_ops == "VSWAP_6_6") {
        return setGate_vswap<VSWAP_Gate_6_6>(ss, 6);
    } else {
        cerr << "[setGate]: Not implemented yet." << endl;
        exit(1);
    }
}

void printGateName(Gate *gate) {
    cout << gate->name << endl
         << flush;
    for (auto t : gate->targs) {
        cout << t << " " << flush;
    }
    cout << endl << flush;
}
void Simulator::setupCircuit_MPI(string cir) {
    ifstream cirfile;
    cirfile.open(cir);

    string line;
    while (getline(cirfile, line)) {
        circuit.push_back(setGate_MPI(line));
        // printGateName(circuit.back());
    }
    cirfile.close();
}
// create a circuit by circuit file
void Simulator::setupCircuit_IO(string cir) {
    // read from circuit file
    ifstream cirfile;
    cirfile.open(cir);

    string line;
    while (getline(cirfile, line)) {
        circuit.push_back(setGate_IO(line));
        // printGateName(circuit.back());
    }
    cirfile.close();
}

// create a subcircuit from the circuit
void Simulator::setupSubCircuits_IO(string cir) {
    // read from circuit file
    ifstream cirfile;
    cirfile.open(cir);
    string line;
    while (getline(cirfile, line)) {
        int sc_size = stoi(line);
        vector<Gate *> subcircuit;
        for (int i = 0; i < sc_size; i++)
        {
            getline(cirfile, line);
            subcircuit.push_back(setGate_IO(line));
            // printGateName(subcircuit.back());
        }
        subcircuits.push_back(subcircuit);
    }
    cirfile.close();
}

/*DIO version*/
Gate *Simulator::setGate_DIO(string &line) {
    stringstream ss;
    ss << line;
    string gate_ops;
    ss >> gate_ops;
    if (gate_ops == "U1") {
        return setUnitary<U1_Gate_DIO>(ss, ONE_QUBIT);
    } else if (gate_ops == "U2") {
        return setUnitary<U2_Gate_DIO>(ss, TWO_QUBIT);
    } else if (gate_ops == "H"){
        return setGate_1<H_Gate_DIO>(ss);
    } else if (gate_ops == "X"){
        return setGate_1<X_Gate_DIO>(ss);
    } else if (gate_ops == "Y"){
        return setGate_1<Y_Gate_DIO>(ss);
    } else if (gate_ops == "Z"){
        return setGate_1<Z_Gate_DIO>(ss);
    } else if (gate_ops == "P"){
        return setGate_Phase<Phase_Gate_DIO>(ss, ONE_QUBIT);
    } else if (gate_ops == "RX"){
        return setGate_Phase<RX_Gate_DIO>(ss, ONE_QUBIT);
    } else if (gate_ops == "RY"){
        return setGate_Phase<RY_Gate_DIO>(ss, ONE_QUBIT);
    } else if (gate_ops == "RZ"){
        return setGate_Phase<RZ_Gate_DIO>(ss, ONE_QUBIT);
    } else if (gate_ops == "SWAP") {
        return setGate_2<SWAP_Gate_DIO>(ss);
    } else if (gate_ops == "CP"){
        return setGate_Phase<CPhase_Gate_DIO>(ss, TWO_QUBIT);
    } else if (gate_ops == "RZZ"){
        return setGate_Phase<RZZ_Gate_DIO>(ss, TWO_QUBIT);
    } else if (gate_ops == "VSWAP_1_1") {
        return setGate_vswap<VSWAP_Gate_1_1_DIO>(ss, 1);
    } else if (gate_ops == "VSWAP_2_2") {
        return setGate_vswap<VSWAP_Gate_2_2_DIO>(ss, 2);
    } else {
        cerr << "[setGate]: Not implemented yet." << endl;
        exit(1);
    }
}

// create a circuit by circuit file
void Simulator::setupCircuit_DIO(string cir) {
    // read from circuit file
    ifstream cirfile;
    cirfile.open(cir);

    string line;
    while (getline(cirfile, line)) {
        circuit.push_back(setGate_DIO(line));
        // printGateName(circuit.back());
    }
    cirfile.close();
}

// create a subcircuit from the circuit
void Simulator::setupSubCircuits_DIO(string cir) {
    // read from circuit file
    ifstream cirfile;
    cirfile.open(cir);
    string line;
    while (getline(cirfile, line)) {
        int sc_size = stoi(line);
        vector<Gate *> subcircuit;
        for (int i = 0; i < sc_size; i++)
        {
            getline(cirfile, line);
            subcircuit.push_back(setGate_DIO(line));
            // printGateName(subcircuit.back());
        }
        subcircuits.push_back(subcircuit);
    }
    cirfile.close();
}

/*MEM version*/
Gate *Simulator::setGate_MEM(string &line) {
    // static int count = 0;
    // cerr << "[setGate] " << count++ << "th gate." << endl;
    stringstream ss;
    ss << line;
    string gate_ops;
    ss >> gate_ops;
    if (gate_ops == "H") {
        return setGate_1<H_Gate_MEM>(ss);
    } else if (gate_ops == "X") {
        return setGate_1<X_Gate_MEM>(ss);
    } else if (gate_ops == "Y") {
        return setGate_1<Y_Gate_MEM>(ss);
    } else if (gate_ops == "Z") {
        return setGate_1<Z_Gate_MEM>(ss);
    } else if (gate_ops == "P") {
        return setGate_Phase<Phase_Gate_MEM>(ss, ONE_QUBIT);
    } else if (gate_ops == "RX") {
        return setGate_Phase<RX_Gate_MEM>(ss, ONE_QUBIT);
    } else if (gate_ops == "RY") {
        return setGate_Phase<RY_Gate_MEM>(ss, ONE_QUBIT);
    } else if (gate_ops == "RZ") {
        return setGate_Phase<RZ_Gate_MEM>(ss, ONE_QUBIT);
    } else if (gate_ops == "U1") {
        return setUnitary<U1_Gate_MEM>(ss, ONE_QUBIT);
    } else if (gate_ops == "RX") {
        return setGate_Phase<RX_Gate_MEM>(ss, ONE_QUBIT);
    } else if (gate_ops == "U2") {
        return setUnitary<U2_Gate_MEM>(ss, TWO_QUBIT);
    } else if (gate_ops == "SWAP") {
        return setGate_2<SWAP_Gate_MEM>(ss);
    } else if (gate_ops == "CP") {
        return setGate_Phase<CPhase_Gate_MEM>(ss, TWO_QUBIT);
    } else if (gate_ops == "RZZ") {
        return setGate_Phase<RZZ_Gate_MEM>(ss, TWO_QUBIT);
    } else if (gate_ops == "U3") {
        return setUnitary<U3_Gate_MEM>(ss, THREE_QUBIT);
    } else if (gate_ops == "VSWAP_1_1") {
        return setGate_vswap<VSWAP_Gate_1_1_MEM>(ss, 1);
    } else if (gate_ops == "VSWAP_2_2") {
        return setGate_vswap<VSWAP_Gate_2_2_MEM>(ss, 2);
    } else if (gate_ops == "VSWAP_3_3") {
        return setGate_vswap<VSWAP_Gate_3_3_MEM>(ss, 3);
    } else if (gate_ops == "VSWAP_4_4") {
        return setGate_vswap<VSWAP_Gate_4_4_MEM>(ss, 4);
    } else if (gate_ops == "VSWAP_6_6") {
        return setGate_vswap<VSWAP_Gate_6_6_MEM>(ss, 6);
    } else {
        cerr << "[setGate]: Not implemented yet." << endl;
        exit(1);
    }
}

/*MEM version*/
// create a circuit by circuit file
void Simulator::setupCircuit_MEM(string cir) {
    // read from circuit file
    ifstream cirfile;
    cirfile.open(cir);

    string line;
    while (getline(cirfile, line)) {
        circuit.push_back(setGate_MEM(line));
        //double vm, rss;
        //process_mem_usage(vm, rss);
        //cout << "[Setup MEM]VM: " << vm << "; RSS: " << rss << endl;
        // printGateName(circuit.back());
    }
    cirfile.close();
}

/*MEM version*/
// create a subcircuit from the circuit
void Simulator::setupSubCircuits_MEM(string cir) {
    // read from circuit file
    ifstream cirfile;
    cirfile.open(cir);
    string line;
    while (getline(cirfile, line)) {
        int sc_size = stoi(line);
        vector<Gate *> subcircuit;
        for (int i = 0; i < sc_size; i++)
        {
            getline(cirfile, line);
            subcircuit.push_back(setGate_MEM(line));
            // printGateName(subcircuit.back());
        }
        subcircuits.push_back(subcircuit);
    }
    cirfile.close();
}

// initialize state file to |0>
void Simulator::setupStateFile() {
    vector<complex<double>> buffer(env.chunk_state, complex<double>(0.0, 0.0));

#pragma omp parallel
    {
        int t = omp_get_thread_num();
        long long fd_off = 0;
        for (long long sz = 0; sz < env.file_state; sz += env.chunk_state) {
            if (pwrite(env.fd_arr[t], static_cast<void *>(buffer.data()), env.chunk_size, fd_off))
                ;
            fd_off += env.chunk_size;
        }
    }
    if(env.rank == 0)
    {
        buffer[0] = complex<double>(1.0, 0.0);
        if (pwrite(env.fd_arr[0], static_cast<void *>(buffer.data()), sizeof(complex<double>), 0))
            ;
    }
}

bool checkValue(vector<complex<double>> &sv){
    complex<double> q = sv[0];
    for(auto &x : sv){
        if(x != q)
            return false;
    }
    return true;
}



void Simulator::run() {

    cerr << "[Simulator]: RUN!!!!!!\n"
         << flush;
    
    std::cerr << env.is_subcircuit << std::endl;

    double start = omp_get_wtime();
    if(env.is_subcircuit)
        Runner->run(subcircuits);
    else
        Runner->run(circuit);
    
    double end = omp_get_wtime();
    std::cout << end - start << "s" << std::endl;

    // if(checkValue(((MEM_Runner *)Runner)->buffer))
    //     cout << "[Check value]: all entries the same, " << ((MEM_Runner *)Runner)->buffer[0] << endl;

    // cerr << "before dump" << endl;
    // cerr << "env.dumpfile: " << env.dumpfile << endl;
    if(env.runner_type == "MEM" && env.dumpfile != ""){
        // cerr << "inside dump" << endl;
        ofstream resfile;
        resfile.open(env.dumpfile);
        for(auto &x :static_cast<MEM_Runner *>(Runner)->buffer)
        resfile << fixed << setprecision(11) << x.real() << endl << x.imag() << endl;
        resfile.close();
    }
    return;
}

Simulator::Simulator(string ini, string cir) {
    setupIni(ini);
    cerr << "[Simulator]: Ini setup." << endl;
        
    setupCir(cir);
    cerr << "[Simulator]: Circuit setup." << endl;
    
    if(env.runner_type != "MEM")
        setupStateFile();

    cerr << "[Simulator]: StateFile setup." << endl;

    cerr << "[Simulator]: Ready to go." << endl;

    // cout << "[Simulator]: Early return." << endl << flush;
    // exit(0);
}

Simulator::~Simulator() {
    for(auto &g:circuit)
        delete g;
    for(auto &gg:subcircuits)
    {
        for(auto &g:gg)
            delete g;
    }
    delete Runner;
}

// void inner_loop_1io_sub(thread_IO_task &task, vector<Gate *> &subcircuit) {
//     for (int j = 0; j < task.fd_using.size(); j++) {
//         pread(env.fd_arr[task.fd_using[j]], &(task.buffer[j * env.chunk_state]), env.chunk_size, task.fd_offset_using[j]);
//     }

//     for (auto &g : subcircuit) {
//         g->run(task.buffer);
//     }

//     for (int j = 0; j < task.fd_using.size(); j++) {
//         pwrite(env.fd_arr[task.fd_using[j]], &(task.buffer[j * env.chunk_state]), env.chunk_size, task.fd_offset_using[j]);
//     }
// }
