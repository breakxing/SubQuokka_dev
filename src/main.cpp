#include <iostream>
#include <string>
#include "common.hpp"
#include "circuit.hpp"
#include "simulator.hpp"

// #include <stdlib.h>
#include <omp.h>
// #include "init.h"
// #include "gate.h"
// #include <mpi.h>


int main(int argc, char *argv[]) {
    std::string ini, cir;
    readargs(argc, argv, ini, cir);
    
    Simulator Quokka(ini, cir);
    
    Quokka.run();
    
    return 0;
}
