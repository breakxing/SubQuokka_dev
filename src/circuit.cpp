#include <iostream>
#include <complex>
#include <vector>
#include <string>
#include <functional>
#include "simulator.hpp"
#include "circuit.hpp"
#include <cassert>

using namespace std;
using namespace placeholders;

#define bind_gate_1(gate_name) \
    if(isChunk(targ[0]))\
        run = bind(&gate_name::run_chunk, this, _1);\
    else\
        run = bind(&gate_name::run_nonchunk, this, _1);

#define bind_gate_mpi_1(gate_name)\
    run_one_qubit_mpi_io = bind(&gate_name::run_mpi_chunk, this, _1,_2,_3,_4,_5);

#define bind_gate_mpi_1_dio(gate_name)\
    run_one_qubit_mpi_dio = bind(&gate_name::run_mpi_chunk, this, _1,_2,_3,_4,_5);

#define bind_gate_2(gate_name) \
    if(isChunk(targ[1]))\
        run = bind(&gate_name::run_chunk_chunk, this, _1);\
    else if (isChunk(targ[0]))\
        run = bind(&gate_name::run_nonchunk_chunk, this, _1);\
    else\
        run = bind(&gate_name::run_nonchunk_nonchunk, this, _1);

#define bind_gate_special_swap(gate_name)\
    run_one_qubit_mpi_io = bind(&gate_name::run_mpi_nonchunk_chunk, this, _1,_2,_3,_4,_5);

#define bind_gate_special_swap_dio(gate_name)\
    run_one_qubit_mpi_dio = bind(&gate_name::run_mpi_nonchunk_chunk, this, _1,_2,_3,_4,_5);

#define bind_gate_1_dio(gate_name) \
    if(isChunk(targ[0]))\
        run_dio = bind(&gate_name::run_chunk, this, _1);\
    else\
        run_dio = bind(&gate_name::run_nonchunk, this, _1);

#define bind_gate_2_dio(gate_name) \
    if(isChunk(targ[1]))\
        run_dio = bind(&gate_name::run_chunk_chunk, this, _1);\
    else if (isChunk(targ[0]))\
        run_dio = bind(&gate_name::run_nonchunk_chunk, this, _1);\
    else\
        run_dio = bind(&gate_name::run_nonchunk_nonchunk, this, _1);

#define bind_gate_1_mem(gate_name) \
    if(isChunk(targ[0]))\
        _run = bind(&gate_name::run_chunk, this, _1, _2);\
    else\
        _run = bind(&gate_name::run_nonchunk, this, _1, _2);

#define bind_gate_mpi_1_mem(gate_name) \
    run_one_qubit_mpi_mem = bind(&gate_name::run_mpi_mem, this, _1,_2,_3,_4,_5);

#define bind_gate_mpi_1_mem_diagonal(gate_name) \
    run_one_qubit_mpi_mem_diagonal = bind(&gate_name::run_mpi_mem, this, _1,_2,_3,_4);

#define bind_gate_swap_mpi_mpi(gate_name) \
    run_one_qubit_mpi_mem = bind(&gate_name::run_mpi_mpi, this, _1,_2,_3,_4,_5);

#define bind_gate_swap_mpi_chunk(gate_name) \
    run_one_qubit_mpi_mem = bind(&gate_name::run_mpi_chunk, this, _1,_2,_3,_4,_5);

#define bind_gate_swap_mpi_file_middle(gate_name) \
    run_one_qubit_mpi_mem = bind(&gate_name::run_mpi_file_middle, this, _1,_2,_3,_4,_5);

#define bind_gate_mpi_mpi_u2_mem(gate_name) \
    run_mem_u2 = bind(&gate_name::run_mpi_nonchunk,this,_1,_2,_3,_4,_5,_6,_7,_8,_9);

#define bind_gate_2_mem(gate_name) \
    if(isChunk(targ[1]))\
        _run = bind(&gate_name::run_chunk_chunk, this, _1, _2);\
    else if (isChunk(targ[0]))\
        _run = bind(&gate_name::run_nonchunk_chunk, this, _1, _2);\
    else\
        _run = bind(&gate_name::run_nonchunk_nonchunk, this, _1, _2);

#define chunk_gate(init_type, update_type1, update_type2, gate_func) \
    init_type\
    for (int i = 0; i < env.chunk_state; i += off) {\
        for (int j = 0; j < half_off; j++){\
            gate_func\
            update_type1\
        }\
        update_type2\
    }\

#define nonchunk_gate(init_type, update_type1, gate_func) \
    init_type\
    for (int i = 0; i < env.chunk_state; i ++) {\
        gate_func\
        update_type1\
    }

#define mpi_gate(init_type, round, update_type1, gate_func) \
    init_type\
    for (int i = 0; i < round * env.chunk_state; i ++) {\
        gate_func\
        update_type1\
    }

#define mpi_gate_mem(init_type, length, update_type1, gate_func) \
    init_type\
    for (int i = 0; i < length; i ++) {\
        gate_func\
        update_type1\
    }

#define mpi_nonchunk_chunk(init_type,round, update_type1, update_type2, gate_func) \
    init_type\
    for (int j = 0; j < round * env.chunk_state; j += off_0) {\
        for (int k = 0; k < half_off_0; k++) {\
            gate_func\
            update_type1\
        }\
        update_type2\
    }

#define mpi_nonchunk_chunk_mem(init_type,length, update_type1, update_type2, gate_func) \
    init_type\
    for (int j = 0; j < length; j += off_0) {\
        for (int k = 0; k < half_off_0; k++) {\
            gate_func\
            update_type1\
        }\
        update_type2\
    }

#define mpi_nonchunk_nonchunk_mem(init_type,length, update_type1, gate_func) \
    init_type\
    for (int k = 0; k < length; k++) {\
        gate_func\
        update_type1\
    }

#define chunk_gate_mem(init_type, update_type1, update_type2, gate_func) \
    init_type\
    for (long long i = 0; i < env.chunk_state; i += off) {\
        for (long long j = 0; j < half_off; j++){\
            gate_func\
            update_type1\
        }\
        update_type2\
    }\

#define nonchunk_gate_mem(init_type, update_type1, gate_func) \
    init_type\
    for (long long i = 0; i < env.chunk_state; i ++) {\
        gate_func\
        update_type1\
    }

#define chunk_chunk_gate(init_type, update_type1, update_type2, update_type3, gate_func) \
    init_type\
    for (int i = 0; i < env.chunk_state; i += off_1) {\
        for (int j = 0; j < half_off_1; j += off_0) {\
            for (int k = 0; k < half_off_0; k++) {\
                gate_func\
                update_type1\
            }\
            update_type2\
        }\
        update_type3\
    }

#define nonchunk_chunk_gate(init_type, update_type1, update_type2, gate_func) \
    init_type\
    for (int j = 0; j < env.chunk_state; j += off_0) {\
        for (int k = 0; k < half_off_0; k++) {\
            gate_func\
            update_type1\
        }\
        update_type2\
    }

#define nonchunk_nonchunk_gate(init_type, update_type1, gate_func) \
    init_type\
    for (int k = 0; k < env.chunk_state; k++) {\
        gate_func\
        update_type1\
    }

#define chunk_chunk_gate_mem(init_type, update_type1, update_type2, update_type3, gate_func) \
    init_type\
    for (long long i = 0; i < env.chunk_state; i += off_1) {\
        for (long long j = 0; j < half_off_1; j += off_0) {\
            for (long long k = 0; k < half_off_0; k++) {\
                gate_func\
                update_type1\
            }\
            update_type2\
        }\
        update_type3\
    }

#define nonchunk_chunk_gate_mem(init_type, update_type1, update_type2, gate_func) \
    init_type\
    for (long long j = 0; j < env.chunk_state; j += off_0) {\
        for (long long k = 0; k < half_off_0; k++) {\
            gate_func\
            update_type1\
        }\
        update_type2\
    }

#define nonchunk_nonchunk_gate_mem(init_type, update_type1, gate_func) \
    init_type\
    for (long long k = 0; k < env.chunk_state; k++) {\
        gate_func\
        update_type1\
    }

#define chunk_chunk_chunk_gate
#define nonchunk_chunk_chunk_gate
#define nonchunk_nonchunk_chunk_gate
#define nonchunk_nonchunk_nonchunk_gate

#define init_off1(type, value0) \
type off0 = value0; complex<double> q0;

#define init_off2(type, value0, value1) \
type off0 = value0; type off1 = value1;\
complex<double> q0; complex<double> q1;

#define init_off4(type, value0, value1, value2, value3) \
type off0 = value0; type off1 = value1; type off2 = value2; type off3 = value3;\
complex<double> q0; complex<double> q1; complex<double> q2; complex<double> q3;

#define init_off8(type, value0, value1, value2, value3, value4, value5, value6, value7) \
type off0 = value0; type off1 = value1; type off2 = value2; type off3 = value3; type off4 = value4; type off5 = value5; type off6 = value6; type off7 = value7;\
complex<double> q0; complex<double> q1; complex<double> q2; complex<double> q3; complex<double> q4; complex<double> q5; complex<double> q6; complex<double> q7;

#define update_off1(incre) off0 += incre;
#define update_off2(incre) off0 += incre; off1 += incre;
#define update_off4(incre) off0 += incre; off1 += incre; off2 += incre; off3 += incre;


#define HGATE \
    q0 = buffer[off0]; q1 = buffer[off1];\
    buffer[off0] = 0.70710678118 * (q0 + q1);\
    buffer[off1] = 0.70710678118 * (q0 - q1);

#define HGATEMPI \
    q0 = buffer1[off0]; q1 = buffer2[off1];\
    buffer1[off0] = 0.70710678118 * (q0 + q1);\
    buffer2[off1] = 0.70710678118 * (q0 - q1);

#define XGATE \
    q0 = buffer[off0]; q1 = buffer[off1];\
    buffer[off0] = q1; buffer[off1] = q0;

#define XGATEMPI \
    q0 = buffer1[off0]; q1 = buffer2[off1];\
    buffer1[off0] = q1; buffer2[off1] = q0;

#define YGATE \
    q0 = buffer[off0]; q1 = buffer[off1];\
    buffer[off0] = q1 * (-1i); buffer[off1] = q0 * 1i;

#define YGATEMPI \
    q0 = buffer1[off0]; q1 = buffer2[off1];\
    buffer1[off0] = q1 * (-1i); buffer2[off1] = q0 * 1i;

#define ZGATE \
    buffer[off0] *= -1;

#define U1GATE \
    q0 = buffer[off0]; q1 = buffer[off1];\
    buffer[off0] = coeff[0] * q0 + coeff[1] * q1;\
    buffer[off1] = coeff[2] * q0 + coeff[3] * q1;

#define U1GATEMPI \
    q0 = buffer1[off0]; q1 = buffer2[off1];\
    buffer1[off0] = coeff[0] * q0 + coeff[1] * q1;\
    buffer2[off0] = coeff[2] * q0 + coeff[3] * q1;

#define PGATE \
    buffer[off0] *= phi;

#define RXGATE \
    q0 = buffer[off0]; q1 = buffer[off1];\
    buffer[off0] = q0 * cos_Phi_2 + q1 * i_sin_Phi_2;\
    buffer[off1] = q0 * i_sin_Phi_2 + q1 * cos_Phi_2;

#define RXGATEMPI \
    q0 = buffer1[off0]; q1 = buffer2[off1];\
    buffer1[off0] = q0 * cos_Phi_2 + q1 * i_sin_Phi_2;\
    buffer2[off1] = q0 * i_sin_Phi_2 + q1 * cos_Phi_2;

#define RYGATE \
    q0 = buffer[off0]; q1 = buffer[off1];\
    buffer[off0] = q0 * cos_Phi_2 - q1 * sin_Phi_2;\
    buffer[off1] = q0 * sin_Phi_2 + q1 * cos_Phi_2;

#define RYGATEMPI \
    q0 = buffer1[off0]; q1 = buffer2[off1];\
    buffer1[off0] = q0 * cos_Phi_2 - q1 * sin_Phi_2;\
    buffer2[off1] = q0 * sin_Phi_2 + q1 * cos_Phi_2;

#define RZGATE \
    buffer[off0] *= exp_n_iPhi_2;\
    buffer[off1] *= exp_p_iPhi_2;

#define RZ0 \
    buffer[off0] *= exp_n_iPhi_2;

#define RZ1 \
    buffer[off0] *= exp_p_iPhi_2;

#define RZGATEMPI \
    buffer1[off0] *= exp_n_iPhi_2;\
    buffer2[off1] *= exp_p_iPhi_2;

#define U2Gate \
    q0 = buffer[off0];\
    q1 = buffer[off1];\
    q2 = buffer[off2];\
    q3 = buffer[off3];\
    buffer[off0] = coeff[ 0] * q0 + coeff[ 1] * q1 + coeff[ 2] * q2 + coeff[ 3] * q3;\
    buffer[off1] = coeff[ 4] * q0 + coeff[ 5] * q1 + coeff[ 6] * q2 + coeff[ 7] * q3;\
    buffer[off2] = coeff[ 8] * q0 + coeff[ 9] * q1 + coeff[10] * q2 + coeff[11] * q3;\
    buffer[off3] = coeff[12] * q0 + coeff[13] * q1 + coeff[14] * q2 + coeff[15] * q3;

#define U2MPICHUNK \
    q0 = buffer1[off0];\
    q1 = buffer1[off1];\
    q2 = buffer2[off2];\
    q3 = buffer2[off3];\
    buffer1[off0] = coeff[ 0] * q0 + coeff[ 1] * q1 + coeff[ 2] * q2 + coeff[ 3] * q3;\
    buffer1[off1] = coeff[ 4] * q0 + coeff[ 5] * q1 + coeff[ 6] * q2 + coeff[ 7] * q3;\
    buffer2[off2] = coeff[ 8] * q0 + coeff[ 9] * q1 + coeff[10] * q2 + coeff[11] * q3;\
    buffer2[off3] = coeff[12] * q0 + coeff[13] * q1 + coeff[14] * q2 + coeff[15] * q3;

#define U2GENERAL \
    q0 = buffer1[off0];\
    q1 = buffer2[off1];\
    q2 = buffer3[off2];\
    q3 = buffer4[off3];\
    buffer1[off0] = coeff[ 0] * q0 + coeff[ 1] * q1 + coeff[ 2] * q2 + coeff[ 3] * q3;\
    buffer2[off1] = coeff[ 4] * q0 + coeff[ 5] * q1 + coeff[ 6] * q2 + coeff[ 7] * q3;\
    buffer3[off2] = coeff[ 8] * q0 + coeff[ 9] * q1 + coeff[10] * q2 + coeff[11] * q3;\
    buffer4[off3] = coeff[12] * q0 + coeff[13] * q1 + coeff[14] * q2 + coeff[15] * q3;

#define SWAPGate \
    q0 = buffer[off0];\
    q1 = buffer[off1];\
    buffer[off0] = q1;\
    buffer[off1] = q0;

#define SWAPGateMPIChunk \
    q0 = buffer1[off0]; q1 = buffer2[off1];\
    buffer1[off0] = q1;\
    buffer2[off1] = q0;

#define CPGATE \
    buffer[off0] *= exp_iPhi;

#define RZZGATE \
    buffer[off0] *= exp_n_iPhi_2;\
    buffer[off1] *= exp_p_iPhi_2;\
    buffer[off2] *= exp_p_iPhi_2;\
    buffer[off3] *= exp_n_iPhi_2;

#define RZZN \
    buffer[off0] *= exp_n_iPhi_2;

#define RZZP \
    buffer[off0] *= exp_p_iPhi_2;

Gate::Gate(vector<int> targ): targs(targ){
    for(auto &x: targ){
        if(isMpi(x))
            mpi_count++;
        else if(isFile(x))
            file_count++;
        else if(isMiddle(x))
            middle_count++;
        else if(isChunk(x))
            chunk_count++;
    }
    nonchunk_count = mpi_count + file_count + middle_count;

    type = ERROR;
}

/*-----ONE_QUBIT_GATE-----*/

ONE_QUBIT_GATE::ONE_QUBIT_GATE(vector<int> targ): Gate(targ){
    type = ONE_QUBIT;
    if(isMpi(targ[0]))
        half_off = env.chunk_state;
    else if(isChunk(targ[0]))
        half_off = env.qubit_offset[targ[0]];
    else
        half_off = env.chunk_state;

    off = half_off * 2;
}

H_Gate::H_Gate(vector<int> targ): ONE_QUBIT_GATE(targ) {
    name = "H_Gate";
    if(mpi_count)
    {
        bind_gate_mpi_1(H_Gate)
    }   
    else
    {
        bind_gate_1(H_Gate)
    }
}

void H_Gate::run_chunk(vector<complex<double>> &buffer){
    chunk_gate(init_off2(int, 0, half_off), update_off2(1), update_off2(half_off), HGATE)
}

void H_Gate::run_nonchunk(vector<complex<double>> &buffer){
    nonchunk_gate(init_off2(int, 0, half_off), update_off2(1), HGATE)
}

void H_Gate::run_mpi_chunk(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,int offset0,int offset1,int round){
    mpi_gate(init_off2(int, offset0 * env.chunk_state, offset1 * env.chunk_state),round, update_off2(1), HGATEMPI)
}

X_Gate::X_Gate(vector<int> targ): ONE_QUBIT_GATE(targ) {
    name = "X_Gate";
    if(mpi_count)
    {
        bind_gate_mpi_1(X_Gate)
    }   
    else
    {
        bind_gate_1(X_Gate)
    }
}

void X_Gate::run_chunk(vector<complex<double>> &buffer){
    chunk_gate(init_off2(int, 0, half_off), update_off2(1), update_off2(half_off), XGATE)
}

void X_Gate::run_nonchunk(vector<complex<double>> &buffer){
    nonchunk_gate(init_off2(int, 0, half_off), update_off2(1), XGATE)
}

void X_Gate::run_mpi_chunk(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,int offset0,int offset1,int round){
    mpi_gate(init_off2(int, offset0 * env.chunk_state, offset1 * env.chunk_state),round, update_off2(1), XGATEMPI)
}

Y_Gate::Y_Gate(vector<int> targ): ONE_QUBIT_GATE(targ) {
    name = "Y_Gate";
    if(mpi_count)
    {
        bind_gate_mpi_1(Y_Gate)
    }   
    else
    {
        bind_gate_1(Y_Gate)
    }
}

void Y_Gate::run_chunk(vector<complex<double>> &buffer){
    chunk_gate(init_off2(int, 0, half_off), update_off2(1), update_off2(half_off), YGATE)
}

void Y_Gate::run_nonchunk(vector<complex<double>> &buffer){
    nonchunk_gate(init_off2(int, 0, half_off), update_off2(1), YGATE)
}

void Y_Gate::run_mpi_chunk(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,int offset0,int offset1,int round){
    mpi_gate(init_off2(int, offset0 * env.chunk_state, offset1 * env.chunk_state),round, update_off2(1), YGATEMPI)
}


Z_Gate::Z_Gate(vector<int> targ): ONE_QUBIT_GATE(targ) {
    name = "Z_Gate";
    bind_gate_1(Z_Gate)
}

void Z_Gate::run_chunk(vector<complex<double>> &buffer){
    chunk_gate(init_off1(int, half_off), update_off1(1), update_off1(half_off), ZGATE)
}

void Z_Gate::run_nonchunk(vector<complex<double>> &buffer){
    nonchunk_gate(init_off1(int, half_off), update_off1(1), ZGATE)
}

U1_Gate::U1_Gate(vector<int> targ, vector<complex<double>> coeff): ONE_QUBIT_GATE(targ), coeff(coeff) {
    name = "U1_Gate";
    if(mpi_count)
    {
        bind_gate_mpi_1(U1_Gate)
    }   
    else
    {
        bind_gate_1(U1_Gate)
    }
}

void U1_Gate::run_chunk(vector<complex<double>> &buffer){
    chunk_gate(init_off2(int, 0, half_off), update_off2(1), update_off2(half_off), U1GATE)
}

void U1_Gate::run_nonchunk(vector<complex<double>> &buffer){
    nonchunk_gate(init_off2(int, 0, half_off), update_off2(1), U1GATE)
}

void U1_Gate::run_mpi_chunk(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,int offset0,int offset1,int round){
    mpi_gate(init_off2(int, offset0 * env.chunk_state, offset1 * env.chunk_state),round, update_off2(1), U1GATEMPI)
}

Phase_Gate::Phase_Gate(vector<int> targ, double phi): ONE_QUBIT_GATE(targ), phi(cos(phi), sin(phi)) {
    name = "Phase_Gate";
    bind_gate_1(Phase_Gate)
}

void Phase_Gate::run_chunk(vector<complex<double>> &buffer){
    chunk_gate(init_off1(int, half_off), update_off1(1), update_off1(half_off), PGATE)
}

void Phase_Gate::run_nonchunk(vector<complex<double>> &buffer){
    nonchunk_gate(init_off1(int, half_off), update_off1(1), PGATE)
}

RX_Gate::RX_Gate(vector<int> targ, double phi): ONE_QUBIT_GATE(targ), cos_Phi_2(cos(phi/2), 0), i_sin_Phi_2(0, -sin(phi/2)) {
    name = "RX_Gate";
    if(mpi_count)
    {
        bind_gate_mpi_1(RX_Gate)
    }   
    else
    {
        bind_gate_1(RX_Gate)
    }
}

void RX_Gate::run_chunk(vector<complex<double>> &buffer){
    chunk_gate(init_off2(int, 0, half_off), update_off2(1), update_off2(half_off), RXGATE)
}

void RX_Gate::run_nonchunk(vector<complex<double>> &buffer){
    nonchunk_gate(init_off2(int, 0, half_off), update_off2(1), RXGATE)
}

void RX_Gate::run_mpi_chunk(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,int offset0,int offset1,int round){
    mpi_gate(init_off2(int, offset0 * env.chunk_state, offset1 * env.chunk_state),round, update_off2(1), RXGATEMPI)
}

RY_Gate::RY_Gate(vector<int> targ, double phi): ONE_QUBIT_GATE(targ), cos_Phi_2(cos(phi/2), 0), sin_Phi_2(sin(phi/2), 0) {
    name = "RY_Gate";
    if(mpi_count)
    {
        bind_gate_mpi_1(RY_Gate)
    }   
    else
    {
        bind_gate_1(RY_Gate)
    }
}

void RY_Gate::run_chunk(vector<complex<double>> &buffer){
    chunk_gate(init_off2(int, 0, half_off), update_off2(1), update_off2(half_off), RYGATE)
}

void RY_Gate::run_nonchunk(vector<complex<double>> &buffer){
    nonchunk_gate(init_off2(int, 0, half_off), update_off2(1), RYGATE)
}

void RY_Gate::run_mpi_chunk(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,int offset0,int offset1,int round){
    mpi_gate(init_off2(int, offset0 * env.chunk_state, offset1 * env.chunk_state),round, update_off2(1), RYGATEMPI)
}

RZ_Gate::RZ_Gate(vector<int> targ, double phi): ONE_QUBIT_GATE(targ), exp_p_iPhi_2(cos(phi/2), sin(phi/2)), exp_n_iPhi_2(cos(-phi/2), sin(-phi/2)) {
    name = "RZ_Gate";
    bind_gate_1(RZ_Gate)
}

void RZ_Gate::run_chunk(vector<complex<double>> &buffer){
    chunk_gate(init_off2(int, 0, half_off), update_off2(1), update_off2(half_off), RZGATE)
}

void RZ_Gate::run_nonchunk(vector<complex<double>> &buffer){
    nonchunk_gate(init_off2(int, 0, half_off), update_off2(1), RZGATE)
}
/*-----TWO_QUBIT_GATE-----*/
TWO_QUBIT_GATE::TWO_QUBIT_GATE(vector<int> targ): Gate(targ) {
    type = TWO_QUBIT;
    sort(targ.begin(), targ.end());
    // cerr << targ[0] << " " << targ[1] << endl;
    if(isChunk(targ[1])){
        half_off_0 = env.qubit_offset[targ[0]];
        half_off_1 = env.qubit_offset[targ[1]];
    }
    else if (isChunk(targ[0])){
        half_off_0 = env.qubit_offset[targ[0]];
        half_off_1 = env.chunk_state;
    }
    else {
        half_off_0 = env.chunk_state;
        half_off_1 = env.chunk_state * 2;
    }
    
    off_0 = half_off_0 * 2;
    off_1 = half_off_1 * 2;
}

U2_Gate::U2_Gate(vector<int> targ, vector<complex<double>> matrix): TWO_QUBIT_GATE(targ){
    name = "U2_Gate";
    assert(targ[0] < targ[1]);
    for (int i = 0; i < 16; i++) {
        int res = 0;
        int b1 = i & 1;
        int b2 = (i >> 1) & 1;
        int b3 = (i >> 2) & 1;
        int b4 = (i >> 3) & 1;
        res = (b3 << 3) | (b4 << 2) | (b1 << 1) | b2;
        coeff[i] = matrix[res];
    }
    if(mpi_count == 1 && chunk_count == 1)
    {
        bind_gate_mpi_1(U2_Gate)
    }
    else
    {
        bind_gate_2(U2_Gate)
    }
}

void U2_Gate::run_chunk_chunk(vector<complex<double>> &buffer){
    chunk_chunk_gate(init_off4(int, 0, half_off_0, half_off_1, (half_off_0 + half_off_1)),
                    update_off4(1), 
                    update_off4(half_off_0), 
                    update_off4(half_off_1), 
                    U2Gate)
}

void U2_Gate::run_nonchunk_chunk(vector<complex<double>> &buffer){
    nonchunk_chunk_gate(init_off4(int, 0, half_off_0, half_off_1, (half_off_0 + half_off_1)),
                    update_off4(1), 
                    update_off4(half_off_0), 
                    U2Gate)
}

void U2_Gate::run_nonchunk_nonchunk(vector<complex<double>> &buffer){
    nonchunk_nonchunk_gate(init_off4(int, 0, half_off_0, half_off_1, (half_off_0 + half_off_1)),
                    update_off4(1), 
                    U2Gate)
}

void U2_Gate::run_mpi_chunk(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,int offset0,int offset1,int round){
    mpi_nonchunk_chunk(init_off4(int,offset0 * env.chunk_state,offset0 * env.chunk_state + half_off_0,offset1 * env.chunk_state, offset1 * env.chunk_state + half_off_0),
                    round,
                    update_off4(1),
                    update_off4(half_off_0),
                    U2MPICHUNK)
}

SWAP_Gate::SWAP_Gate(vector<int> targ): TWO_QUBIT_GATE(targ) {
    name = "SWAP_Gate";
    if(mpi_count == 1 && chunk_count == 1)
    {
        bind_gate_special_swap(SWAP_Gate)
    }
    else
    {
        bind_gate_2(SWAP_Gate)
    }
}

void SWAP_Gate::run_chunk_chunk(vector<complex<double>> &buffer){
    chunk_chunk_gate(init_off2(int, half_off_0, half_off_1),
                    update_off2(1), 
                    update_off2(half_off_0), 
                    update_off2(half_off_1), 
                    SWAPGate)
}

void SWAP_Gate::run_nonchunk_chunk(vector<complex<double>> &buffer){
    nonchunk_chunk_gate(init_off2(int, half_off_0, half_off_1),
                    update_off2(1), 
                    update_off2(half_off_0), 
                    SWAPGate)
}

void SWAP_Gate::run_nonchunk_nonchunk(vector<complex<double>> &buffer){
    nonchunk_nonchunk_gate(init_off2(int, half_off_0, half_off_1),
                    update_off2(1), 
                    SWAPGate)
}

void SWAP_Gate::run_mpi_nonchunk_chunk(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,int offset0,int offset1,int round){
    mpi_nonchunk_chunk(init_off2(int, offset0 * env.chunk_state + half_off_0, offset1 * env.chunk_state),round, update_off2(1),update_off2(half_off_0),SWAPGateMPIChunk)
}

VSWAP_Gate_1_1::VSWAP_Gate_1_1(vector<int> targ): SWAP_Gate(targ) {
    type = VSWAP;
    name = "VSWAP_Gate_1_1";
}

MPI_VSWAP_Gate_1_1_IO::MPI_VSWAP_Gate_1_1_IO(vector<int> targ): SWAP_Gate(targ) {
    type = VSWAP;
    name = "MPI_VSWAP_Gate_1_1_IO";
}


CPhase_Gate::CPhase_Gate(vector<int> targ, double phi): TWO_QUBIT_GATE(targ) {
    name = "CPhase_Gate";
    exp_iPhi = complex<double>(cos(phi), sin(phi));
    bind_gate_2(CPhase_Gate)
}

void CPhase_Gate::run_chunk_chunk(vector<complex<double>> &buffer){
    chunk_chunk_gate(init_off1(int, (half_off_0 + half_off_1)),
                    update_off1(1), 
                    update_off1(half_off_0), 
                    update_off1(half_off_1), 
                    CPGATE)
}

void CPhase_Gate::run_nonchunk_chunk(vector<complex<double>> &buffer){
    nonchunk_chunk_gate(init_off1(int, (half_off_0 + half_off_1)),
                    update_off1(1), 
                    update_off1(half_off_0), 
                    CPGATE)
}

void CPhase_Gate::run_nonchunk_nonchunk(vector<complex<double>> &buffer){
    nonchunk_nonchunk_gate(init_off1(int, (half_off_0 + half_off_1)),
                    update_off1(1), 
                    CPGATE)
}

RZZ_Gate::RZZ_Gate(vector<int> targ, double phi): TWO_QUBIT_GATE(targ) {
    name += "RZZ_Gate";
    exp_p_iPhi_2 = complex<double>(cos(phi/2), sin(phi/2));
    exp_n_iPhi_2 = complex<double>(cos(-phi/2), sin(-phi/2));
    bind_gate_2(RZZ_Gate)
}

void RZZ_Gate::run_chunk_chunk(vector<complex<double>> &buffer){
    chunk_chunk_gate(init_off4(int, 0, half_off_0, half_off_1, (half_off_0 + half_off_1)),
                    update_off4(1), 
                    update_off4(half_off_0), 
                    update_off4(half_off_1), 
                    RZZGATE)
}

void RZZ_Gate::run_nonchunk_chunk(vector<complex<double>> &buffer){
    nonchunk_chunk_gate(init_off4(int, 0, half_off_0, half_off_1, (half_off_0 + half_off_1)),
                    update_off4(1), 
                    update_off4(half_off_0), 
                    RZZGATE)
}

void RZZ_Gate::run_nonchunk_nonchunk(vector<complex<double>> &buffer){
    nonchunk_nonchunk_gate(init_off4(int, 0, half_off_0, half_off_1, (half_off_0 + half_off_1)),
                    update_off4(1), 
                    RZZGATE)
}

/*-----THREE_QUBIT_GATE-----*/
U3_Gate::U3_Gate(vector<int> targ, vector<complex<double>> coeff): Gate(targ), coeff(coeff) {
    type = THREE_QUBIT;
    if(isChunk(targ[2])){
        half_target_offset0 = env.qubit_offset[targ[0]];
        half_target_offset1 = env.qubit_offset[targ[1]];
        half_target_offset2 = env.qubit_offset[targ[2]];
        
        target_offset0 = half_target_offset0 * 2;
        target_offset1 = half_target_offset1 * 2;
        target_offset2 = half_target_offset2 * 2;
        
        run = bind(&U3_Gate::run_chunk_chunk_chunk, this, _1);
        name += "U3_Gate_chunk_chunk_chunk";
    }
    if(isChunk(targ[1])){
        half_target_offset0 = env.qubit_offset[targ[0]];
        half_target_offset1 = env.qubit_offset[targ[1]];
        half_target_offset2 = env.chunk_state;

        
        target_offset0 = half_target_offset0 * 2;
        target_offset1 = half_target_offset1 * 2;
        target_offset2 = half_target_offset2 * 2;
        
        run = bind(&U3_Gate::run_nonchunk_chunk_chunk, this, _1);
        name += "U3_Gate_nonchunk_chunk_chunk";
    }
    else if (isChunk(targ[0])){
        half_target_offset0 = env.qubit_offset[targ[0]];
        half_target_offset1 = env.chunk_state;
        half_target_offset2 = env.chunk_state * 2;
        
        target_offset0 = half_target_offset0 * 2;
        target_offset1 = half_target_offset1 * 2;
        target_offset2 = half_target_offset2 * 2;

        run = bind(&U3_Gate::run_nonchunk_nonchunk_chunk, this, _1);
        name += "U3_Gate_nonchunk_nonchunk_chunk";
    }
    else {
        half_target_offset0 = env.chunk_state;
        half_target_offset1 = env.chunk_state * 2;
        half_target_offset2 = env.chunk_state * 4;
        
        target_offset0 = half_target_offset0 * 2;
        target_offset1 = half_target_offset1 * 2;
        target_offset2 = half_target_offset2 * 2;

        run = bind(&U3_Gate::run_nonchunk_nonchunk_nonchunk, this, _1);
        name += "U3_Gate_nonchunk_nonchunk_nonchunk";
    }
}

void U3_Gate::run_chunk_chunk_chunk(vector<complex<double>> &buffer){
    int off000 = 0;
    int off001 = half_target_offset0;
    int off010 = half_target_offset1;
    int off011 = half_target_offset1 + half_target_offset0;
    int off100 = half_target_offset2;
    int off101 = half_target_offset2 + half_target_offset0;
    int off110 = half_target_offset2 + half_target_offset1;
    int off111 = half_target_offset2 + half_target_offset1 + half_target_offset0;
    
    complex<double> q000;
    complex<double> q001;
    complex<double> q010;
    complex<double> q011;
    complex<double> q100;
    complex<double> q101;
    complex<double> q110;
    complex<double> q111;

    for (int i = 0; i < env.chunk_state; i += target_offset2) {
        for (int j = 0; j < half_target_offset2; j += target_offset1) {
            for (int k = 0; k < half_target_offset1; k += target_offset0) {
                for (int m = 0; m < half_target_offset0; m++) {
                    q000 = buffer[off000];
                    q001 = buffer[off001];
                    q010 = buffer[off010];
                    q011 = buffer[off011];
                    q100 = buffer[off100];
                    q101 = buffer[off101];
                    q110 = buffer[off110];
                    q111 = buffer[off111];
                    
                    buffer[off000] = coeff[ 0] * q000 + coeff[ 1] * q001 + coeff[ 2] * q010 + coeff[ 3] * q011 + coeff[ 4] * q100 + coeff[ 5] * q101 + coeff[ 6] * q110 + coeff[ 7] * q111;
                    buffer[off001] = coeff[ 8] * q000 + coeff[ 9] * q001 + coeff[10] * q010 + coeff[11] * q011 + coeff[12] * q100 + coeff[13] * q101 + coeff[14] * q110 + coeff[15] * q111;
                    buffer[off010] = coeff[16] * q000 + coeff[17] * q001 + coeff[18] * q010 + coeff[19] * q011 + coeff[20] * q100 + coeff[21] * q101 + coeff[22] * q110 + coeff[23] * q111;
                    buffer[off011] = coeff[24] * q000 + coeff[25] * q001 + coeff[26] * q010 + coeff[27] * q011 + coeff[28] * q100 + coeff[29] * q101 + coeff[30] * q110 + coeff[31] * q111;
                    buffer[off100] = coeff[32] * q000 + coeff[33] * q001 + coeff[34] * q010 + coeff[35] * q011 + coeff[36] * q100 + coeff[37] * q101 + coeff[38] * q110 + coeff[39] * q111;
                    buffer[off101] = coeff[40] * q000 + coeff[41] * q001 + coeff[42] * q010 + coeff[43] * q011 + coeff[44] * q100 + coeff[45] * q101 + coeff[46] * q110 + coeff[47] * q111;
                    buffer[off110] = coeff[48] * q000 + coeff[49] * q001 + coeff[50] * q010 + coeff[51] * q011 + coeff[52] * q100 + coeff[53] * q101 + coeff[54] * q110 + coeff[55] * q111;
                    buffer[off111] = coeff[56] * q000 + coeff[57] * q001 + coeff[58] * q010 + coeff[59] * q011 + coeff[60] * q100 + coeff[61] * q101 + coeff[62] * q110 + coeff[63] * q111;
                    
                    off000++;
                    off001++;
                    off010++;
                    off011++;
                    off100++;
                    off101++;
                    off110++;
                    off111++;
                }
                off000 += half_target_offset0;
                off001 += half_target_offset0;
                off010 += half_target_offset0;
                off011 += half_target_offset0;
                off100 += half_target_offset0;
                off101 += half_target_offset0;
                off110 += half_target_offset0;
                off111 += half_target_offset0;
            }
            off000 += half_target_offset1;
            off001 += half_target_offset1;
            off010 += half_target_offset1;
            off011 += half_target_offset1;
            off100 += half_target_offset1;
            off101 += half_target_offset1;
            off110 += half_target_offset1;
            off111 += half_target_offset1;
        }
        off000 += half_target_offset2;
        off001 += half_target_offset2;
        off010 += half_target_offset2;
        off011 += half_target_offset2;
        off100 += half_target_offset2;
        off101 += half_target_offset2;
        off110 += half_target_offset2;
        off111 += half_target_offset2;
    }
}

void U3_Gate::run_nonchunk_chunk_chunk(vector<complex<double>> &buffer){
    int off000 = 0;
    int off001 = half_target_offset0;
    int off010 = half_target_offset1;
    int off011 = half_target_offset1 + half_target_offset0;
    int off100 = half_target_offset2;
    int off101 = half_target_offset2 + half_target_offset0;
    int off110 = half_target_offset2 + half_target_offset1;
    int off111 = half_target_offset2 + half_target_offset1 + half_target_offset0;
    
    complex<double> q000;
    complex<double> q001;
    complex<double> q010;
    complex<double> q011;
    complex<double> q100;
    complex<double> q101;
    complex<double> q110;
    complex<double> q111;

    for (int j = 0; j < env.chunk_state; j += target_offset1) {
        for (int k = 0; k < half_target_offset1; k += target_offset0) {
            for (int m = 0; m < half_target_offset0; m++) {
                q000 = buffer[off000];
                q001 = buffer[off001];
                q010 = buffer[off010];
                q011 = buffer[off011];
                q100 = buffer[off100];
                q101 = buffer[off101];
                q110 = buffer[off110];
                q111 = buffer[off111];
                
                buffer[off000] = coeff[ 0] * q000 + coeff[ 1] * q001 + coeff[ 2] * q010 + coeff[ 3] * q011 + coeff[ 4] * q100 + coeff[ 5] * q101 + coeff[ 6] * q110 + coeff[ 7] * q111;
                buffer[off001] = coeff[ 8] * q000 + coeff[ 9] * q001 + coeff[10] * q010 + coeff[11] * q011 + coeff[12] * q100 + coeff[13] * q101 + coeff[14] * q110 + coeff[15] * q111;
                buffer[off010] = coeff[16] * q000 + coeff[17] * q001 + coeff[18] * q010 + coeff[19] * q011 + coeff[20] * q100 + coeff[21] * q101 + coeff[22] * q110 + coeff[23] * q111;
                buffer[off011] = coeff[24] * q000 + coeff[25] * q001 + coeff[26] * q010 + coeff[27] * q011 + coeff[28] * q100 + coeff[29] * q101 + coeff[30] * q110 + coeff[31] * q111;
                buffer[off100] = coeff[32] * q000 + coeff[33] * q001 + coeff[34] * q010 + coeff[35] * q011 + coeff[36] * q100 + coeff[37] * q101 + coeff[38] * q110 + coeff[39] * q111;
                buffer[off101] = coeff[40] * q000 + coeff[41] * q001 + coeff[42] * q010 + coeff[43] * q011 + coeff[44] * q100 + coeff[45] * q101 + coeff[46] * q110 + coeff[47] * q111;
                buffer[off110] = coeff[48] * q000 + coeff[49] * q001 + coeff[50] * q010 + coeff[51] * q011 + coeff[52] * q100 + coeff[53] * q101 + coeff[54] * q110 + coeff[55] * q111;
                buffer[off111] = coeff[56] * q000 + coeff[57] * q001 + coeff[58] * q010 + coeff[59] * q011 + coeff[60] * q100 + coeff[61] * q101 + coeff[62] * q110 + coeff[63] * q111;
                
                off000++;
                off001++;
                off010++;
                off011++;
                off100++;
                off101++;
                off110++;
                off111++;
            }
            off000 += half_target_offset0;
            off001 += half_target_offset0;
            off010 += half_target_offset0;
            off011 += half_target_offset0;
            off100 += half_target_offset0;
            off101 += half_target_offset0;
            off110 += half_target_offset0;
            off111 += half_target_offset0;
        }
        off000 += half_target_offset1;
        off001 += half_target_offset1;
        off010 += half_target_offset1;
        off011 += half_target_offset1;
        off100 += half_target_offset1;
        off101 += half_target_offset1;
        off110 += half_target_offset1;
        off111 += half_target_offset1;
    }
}

void U3_Gate::run_nonchunk_nonchunk_chunk(vector<complex<double>> &buffer){
    int off000 = 0;
    int off001 = half_target_offset0;
    int off010 = half_target_offset1;
    int off011 = half_target_offset1 + half_target_offset0;
    int off100 = half_target_offset2;
    int off101 = half_target_offset2 + half_target_offset0;
    int off110 = half_target_offset2 + half_target_offset1;
    int off111 = half_target_offset2 + half_target_offset1 + half_target_offset0;
    
    complex<double> q000;
    complex<double> q001;
    complex<double> q010;
    complex<double> q011;
    complex<double> q100;
    complex<double> q101;
    complex<double> q110;
    complex<double> q111;

    for (int k = 0; k < env.chunk_state; k += target_offset0) {
        for (int m = 0; m < half_target_offset0; m++) {
            q000 = buffer[off000];
            q001 = buffer[off001];
            q010 = buffer[off010];
            q011 = buffer[off011];
            q100 = buffer[off100];
            q101 = buffer[off101];
            q110 = buffer[off110];
            q111 = buffer[off111];
            
            buffer[off000] = coeff[ 0] * q000 + coeff[ 1] * q001 + coeff[ 2] * q010 + coeff[ 3] * q011 + coeff[ 4] * q100 + coeff[ 5] * q101 + coeff[ 6] * q110 + coeff[ 7] * q111;
            buffer[off001] = coeff[ 8] * q000 + coeff[ 9] * q001 + coeff[10] * q010 + coeff[11] * q011 + coeff[12] * q100 + coeff[13] * q101 + coeff[14] * q110 + coeff[15] * q111;
            buffer[off010] = coeff[16] * q000 + coeff[17] * q001 + coeff[18] * q010 + coeff[19] * q011 + coeff[20] * q100 + coeff[21] * q101 + coeff[22] * q110 + coeff[23] * q111;
            buffer[off011] = coeff[24] * q000 + coeff[25] * q001 + coeff[26] * q010 + coeff[27] * q011 + coeff[28] * q100 + coeff[29] * q101 + coeff[30] * q110 + coeff[31] * q111;
            buffer[off100] = coeff[32] * q000 + coeff[33] * q001 + coeff[34] * q010 + coeff[35] * q011 + coeff[36] * q100 + coeff[37] * q101 + coeff[38] * q110 + coeff[39] * q111;
            buffer[off101] = coeff[40] * q000 + coeff[41] * q001 + coeff[42] * q010 + coeff[43] * q011 + coeff[44] * q100 + coeff[45] * q101 + coeff[46] * q110 + coeff[47] * q111;
            buffer[off110] = coeff[48] * q000 + coeff[49] * q001 + coeff[50] * q010 + coeff[51] * q011 + coeff[52] * q100 + coeff[53] * q101 + coeff[54] * q110 + coeff[55] * q111;
            buffer[off111] = coeff[56] * q000 + coeff[57] * q001 + coeff[58] * q010 + coeff[59] * q011 + coeff[60] * q100 + coeff[61] * q101 + coeff[62] * q110 + coeff[63] * q111;
            
            off000++;
            off001++;
            off010++;
            off011++;
            off100++;
            off101++;
            off110++;
            off111++;
        }
        off000 += half_target_offset0;
        off001 += half_target_offset0;
        off010 += half_target_offset0;
        off011 += half_target_offset0;
        off100 += half_target_offset0;
        off101 += half_target_offset0;
        off110 += half_target_offset0;
        off111 += half_target_offset0;
    }
}

void U3_Gate::run_nonchunk_nonchunk_nonchunk(vector<complex<double>> &buffer){
    int off000 = 0;
    int off001 = half_target_offset0;
    int off010 = half_target_offset1;
    int off011 = half_target_offset1 + half_target_offset0;
    int off100 = half_target_offset2;
    int off101 = half_target_offset2 + half_target_offset0;
    int off110 = half_target_offset2 + half_target_offset1;
    int off111 = half_target_offset2 + half_target_offset1 + half_target_offset0;
    
    complex<double> q000;
    complex<double> q001;
    complex<double> q010;
    complex<double> q011;
    complex<double> q100;
    complex<double> q101;
    complex<double> q110;
    complex<double> q111;

    for (int m = 0; m < env.chunk_state; m++) {
        q000 = buffer[off000];
        q001 = buffer[off001];
        q010 = buffer[off010];
        q011 = buffer[off011];
        q100 = buffer[off100];
        q101 = buffer[off101];
        q110 = buffer[off110];
        q111 = buffer[off111];
        
        buffer[off000] = coeff[ 0] * q000 + coeff[ 1] * q001 + coeff[ 2] * q010 + coeff[ 3] * q011 + coeff[ 4] * q100 + coeff[ 5] * q101 + coeff[ 6] * q110 + coeff[ 7] * q111;
        buffer[off001] = coeff[ 8] * q000 + coeff[ 9] * q001 + coeff[10] * q010 + coeff[11] * q011 + coeff[12] * q100 + coeff[13] * q101 + coeff[14] * q110 + coeff[15] * q111;
        buffer[off010] = coeff[16] * q000 + coeff[17] * q001 + coeff[18] * q010 + coeff[19] * q011 + coeff[20] * q100 + coeff[21] * q101 + coeff[22] * q110 + coeff[23] * q111;
        buffer[off011] = coeff[24] * q000 + coeff[25] * q001 + coeff[26] * q010 + coeff[27] * q011 + coeff[28] * q100 + coeff[29] * q101 + coeff[30] * q110 + coeff[31] * q111;
        buffer[off100] = coeff[32] * q000 + coeff[33] * q001 + coeff[34] * q010 + coeff[35] * q011 + coeff[36] * q100 + coeff[37] * q101 + coeff[38] * q110 + coeff[39] * q111;
        buffer[off101] = coeff[40] * q000 + coeff[41] * q001 + coeff[42] * q010 + coeff[43] * q011 + coeff[44] * q100 + coeff[45] * q101 + coeff[46] * q110 + coeff[47] * q111;
        buffer[off110] = coeff[48] * q000 + coeff[49] * q001 + coeff[50] * q010 + coeff[51] * q011 + coeff[52] * q100 + coeff[53] * q101 + coeff[54] * q110 + coeff[55] * q111;
        buffer[off111] = coeff[56] * q000 + coeff[57] * q001 + coeff[58] * q010 + coeff[59] * q011 + coeff[60] * q100 + coeff[61] * q101 + coeff[62] * q110 + coeff[63] * q111;
        
        off000++;
        off001++;
        off010++;
        off011++;
        off100++;
        off101++;
        off110++;
        off111++;
    }
}


VSWAP_Gate_2_2::VSWAP_Gate_2_2(vector<int> targ): Gate(targ) {
    type = VSWAP;
    name = "VSWAP_Gate_2_2";

    half_off_0 = env.qubit_offset[targ[0]];
    half_off_1 = env.qubit_offset[targ[1]];
    half_off_2 = env.chunk_state;
    half_off_3 = 2 * env.chunk_state;
        
    off_0 = half_off_0 * 2;
    off_1 = half_off_1 * 2;
    off_2 = half_off_2 * 2;
    off_3 = half_off_3 * 2;
        
    run = bind(&VSWAP_Gate_2_2::run_nonchunk_chunk, this, _1);
}

MPI_VSWAP_Gate_2_2_IO::MPI_VSWAP_Gate_2_2_IO(vector<int> targ): Gate(targ) {
    type = VSWAP;
    name = "MPI_VSWAP_Gate_2_2_IO";
    half_off_0 = env.qubit_offset[targ[0]];
    half_off_1 = env.qubit_offset[targ[1]];
    half_off_2 = 0;
    half_off_3 = 0;
        
    off_0 = half_off_0 * 2;
    off_1 = half_off_1 * 2;
    off_2 = half_off_2 * 2;
    off_3 = half_off_3 * 2;
    run_mpi_vswap2_2_io = bind(&MPI_VSWAP_Gate_2_2_IO::run_vswap_2_2,this,_1,_2,_3,_4,_5);
}

void MPI_VSWAP_Gate_2_2_IO::run_vswap_2_2(vector<complex<double>> &buffer1,vector<complex<double>> &buffer2,vector<complex<double>> &buffer3,vector<complex<double>> &buffer4,int round)
{
    int off0001 = half_off_0;
    int off0010 = half_off_1;
    int off0011 = half_off_0 + half_off_1;
    int off0100 = half_off_2;
    int off0110 = half_off_1 + half_off_2;
    int off0111 = half_off_0 + half_off_1 + half_off_2;
    int off1000 = half_off_3;
    int off1001 = half_off_0 + half_off_3;
    int off1011 = half_off_0 + half_off_1 + half_off_3;
    int off1100 = half_off_2 + half_off_3;
    int off1101 = half_off_0 + half_off_2 + half_off_3;
    int off1110 = half_off_1 + half_off_2 + half_off_3;

    complex<double> q0001;
    complex<double> q0010;
    complex<double> q0011;
    complex<double> q0100;
    complex<double> q0110;
    complex<double> q0111;
    complex<double> q1000;
    complex<double> q1001;
    complex<double> q1011;
    complex<double> q1100;
    complex<double> q1101;
    complex<double> q1110;

    for (int i = 0; i < round * env.chunk_state; i += off_1) {
        for (int j = 0; j < half_off_1; j += off_0) {
            for (int k = 0; k < half_off_0; k++) {
                q0001 = buffer1[off0001];
                q0010 = buffer1[off0010];
                q0011 = buffer1[off0011];
                q0100 = buffer2[off0100];
                q0110 = buffer2[off0110];
                q0111 = buffer2[off0111];
                q1000 = buffer3[off1000];
                q1001 = buffer3[off1001];
                q1011 = buffer3[off1011];
                q1100 = buffer4[off1100];
                q1101 = buffer4[off1101];
                q1110 = buffer4[off1110];
                buffer1[off0001] = q0100;
                buffer1[off0010] = q1000;
                buffer1[off0011] = q1100;
                buffer2[off0100] = q0001;
                buffer2[off0110] = q1001;
                buffer2[off0111] = q1101;
                buffer3[off1000] = q0010;
                buffer3[off1001] = q0110;
                buffer3[off1011] = q1110;
                buffer4[off1100] = q0011;
                buffer4[off1101] = q0111;
                buffer4[off1110] = q1011;
                off0001++;
                off0010++;
                off0011++;
                off0100++;
                off0110++;
                off0111++;
                off1000++;
                off1001++;
                off1011++;
                off1100++;
                off1101++;
                off1110++;
            }
            off0001 += half_off_0;
            off0010 += half_off_0;
            off0011 += half_off_0;
            off0100 += half_off_0;
            off0110 += half_off_0;
            off0111 += half_off_0;
            off1000 += half_off_0;
            off1001 += half_off_0;
            off1011 += half_off_0;
            off1100 += half_off_0;
            off1101 += half_off_0;
            off1110 += half_off_0;
        }
        off0001 += half_off_1;
        off0010 += half_off_1;
        off0011 += half_off_1;
        off0100 += half_off_1;
        off0110 += half_off_1;
        off0111 += half_off_1;
        off1000 += half_off_1;
        off1001 += half_off_1;
        off1011 += half_off_1;
        off1100 += half_off_1;
        off1101 += half_off_1;
        off1110 += half_off_1;
    }
}

void VSWAP_Gate_2_2::run_nonchunk_chunk(vector<complex<double>> &buffer){
    int off0001 = half_off_0;
    int off0010 = half_off_1;
    int off0011 = half_off_0 + half_off_1;
    int off0100 = half_off_2;
    int off0110 = half_off_1 + half_off_2;
    int off0111 = half_off_0 + half_off_1 + half_off_2;
    int off1000 = half_off_3;
    int off1001 = half_off_0 + half_off_3;
    int off1011 = half_off_0 + half_off_1 + half_off_3;
    int off1100 = half_off_2 + half_off_3;
    int off1101 = half_off_0 + half_off_2 + half_off_3;
    int off1110 = half_off_1 + half_off_2 + half_off_3;

    complex<double> q0001;
    complex<double> q0010;
    complex<double> q0011;
    complex<double> q0100;
    complex<double> q0110;
    complex<double> q0111;
    complex<double> q1000;
    complex<double> q1001;
    complex<double> q1011;
    complex<double> q1100;
    complex<double> q1101;
    complex<double> q1110;

    for (int i = 0; i < env.chunk_state; i += off_1) {
        for (int j = 0; j < half_off_1; j += off_0) {
            for (int k = 0; k < half_off_0; k++) {
                q0001 = buffer[off0001];
                q0010 = buffer[off0010];
                q0011 = buffer[off0011];
                q0100 = buffer[off0100];
                q0110 = buffer[off0110];
                q0111 = buffer[off0111];
                q1000 = buffer[off1000];
                q1001 = buffer[off1001];
                q1011 = buffer[off1011];
                q1100 = buffer[off1100];
                q1101 = buffer[off1101];
                q1110 = buffer[off1110];
                buffer[off0001] = q0100;
                buffer[off0010] = q1000;
                buffer[off0011] = q1100;
                buffer[off0100] = q0001;
                buffer[off0110] = q1001;
                buffer[off0111] = q1101;
                buffer[off1000] = q0010;
                buffer[off1001] = q0110;
                buffer[off1011] = q1110;
                buffer[off1100] = q0011;
                buffer[off1101] = q0111;
                buffer[off1110] = q1011;
                off0001++;
                off0010++;
                off0011++;
                off0100++;
                off0110++;
                off0111++;
                off1000++;
                off1001++;
                off1011++;
                off1100++;
                off1101++;
                off1110++;
            }
            off0001 += half_off_0;
            off0010 += half_off_0;
            off0011 += half_off_0;
            off0100 += half_off_0;
            off0110 += half_off_0;
            off0111 += half_off_0;
            off1000 += half_off_0;
            off1001 += half_off_0;
            off1011 += half_off_0;
            off1100 += half_off_0;
            off1101 += half_off_0;
            off1110 += half_off_0;
        }
        off0001 += half_off_1;
        off0010 += half_off_1;
        off0011 += half_off_1;
        off0100 += half_off_1;
        off0110 += half_off_1;
        off0111 += half_off_1;
        off1000 += half_off_1;
        off1001 += half_off_1;
        off1011 += half_off_1;
        off1100 += half_off_1;
        off1101 += half_off_1;
        off1110 += half_off_1;
    }
}

VSWAP_Gate_3_3::VSWAP_Gate_3_3(vector<int> targ): Gate(targ) {
    type = VSWAP;
    swap_out = {targ[0], targ[1], targ[2]};
    swap_in  = {targ[3], targ[4], targ[5]};

    half_targ_off0 = env.qubit_offset[targ[0]];
    half_targ_off1 = env.qubit_offset[targ[1]];
    half_targ_off2 = env.qubit_offset[targ[2]];
    half_targ_off3 = env.chunk_state;
    half_targ_off4 = 2 * env.chunk_state;
    half_targ_off5 = 4 * env.chunk_state;
        
    targ_off0 = half_targ_off0 * 2;
    targ_off1 = half_targ_off1 * 2;
    targ_off2 = half_targ_off2 * 2;
    targ_off3 = half_targ_off3 * 2;
    targ_off4 = half_targ_off4 * 2;
    targ_off5 = half_targ_off5 * 2;
        
    run = bind(&VSWAP_Gate_3_3::run_nonchunk_chunk, this, _1);
    name += "VSWAP_Gate_3_3";
}

void VSWAP_Gate_3_3::run_nonchunk_chunk(vector<complex<double>> &buffer){
    int off000001 = half_targ_off0;
    int off000010 = half_targ_off1;
    int off000011 = half_targ_off0 + half_targ_off1;
    int off000100 = half_targ_off2;
    int off000101 = half_targ_off0 + half_targ_off2;
    int off000110 = half_targ_off1 + half_targ_off2;
    int off000111 = half_targ_off0 + half_targ_off1 + half_targ_off2;
    int off001000 = half_targ_off3;
    int off001010 = half_targ_off1 + half_targ_off3;
    int off001011 = half_targ_off0 + half_targ_off1 + half_targ_off3;
    int off001100 = half_targ_off2 + half_targ_off3;
    int off001101 = half_targ_off0 + half_targ_off2 + half_targ_off3;
    int off001110 = half_targ_off1 + half_targ_off2 + half_targ_off3;
    int off001111 = half_targ_off0 + half_targ_off1 + half_targ_off2 + half_targ_off3;
    int off010000 = half_targ_off4;
    int off010001 = half_targ_off0 + half_targ_off4;
    int off010011 = half_targ_off0 + half_targ_off1 + half_targ_off4;
    int off010100 = half_targ_off2 + half_targ_off4;
    int off010101 = half_targ_off0 + half_targ_off2 + half_targ_off4;
    int off010110 = half_targ_off1 + half_targ_off2 + half_targ_off4;
    int off010111 = half_targ_off0 + half_targ_off1 + half_targ_off2 + half_targ_off4;
    int off011000 = half_targ_off3 + half_targ_off4;
    int off011001 = half_targ_off0 + half_targ_off3 + half_targ_off4;
    int off011010 = half_targ_off1 + half_targ_off3 + half_targ_off4;
    int off011100 = half_targ_off2 + half_targ_off3 + half_targ_off4;
    int off011101 = half_targ_off0 + half_targ_off2 + half_targ_off3 + half_targ_off4;
    int off011110 = half_targ_off1 + half_targ_off2 + half_targ_off3 + half_targ_off4;
    int off011111 = half_targ_off0 + half_targ_off1 + half_targ_off2 + half_targ_off3 + half_targ_off4;
    int off100000 = half_targ_off5;
    int off100001 = half_targ_off0 + half_targ_off5;
    int off100010 = half_targ_off1 + half_targ_off5;
    int off100011 = half_targ_off0 + half_targ_off1 + half_targ_off5;
    int off100101 = half_targ_off0 + half_targ_off2 + half_targ_off5;
    int off100110 = half_targ_off1 + half_targ_off2 + half_targ_off5;
    int off100111 = half_targ_off0 + half_targ_off1 + half_targ_off2 + half_targ_off5;
    int off101000 = half_targ_off3 + half_targ_off5;
    int off101001 = half_targ_off0 + half_targ_off3 + half_targ_off5;
    int off101010 = half_targ_off1 + half_targ_off3 + half_targ_off5;
    int off101011 = half_targ_off0 + half_targ_off1 + half_targ_off3 + half_targ_off5;
    int off101100 = half_targ_off2 + half_targ_off3 + half_targ_off5;
    int off101110 = half_targ_off1 + half_targ_off2 + half_targ_off3 + half_targ_off5;
    int off101111 = half_targ_off0 + half_targ_off1 + half_targ_off2 + half_targ_off3 + half_targ_off5;
    int off110000 = half_targ_off4 + half_targ_off5;
    int off110001 = half_targ_off0 + half_targ_off4 + half_targ_off5;
    int off110010 = half_targ_off1 + half_targ_off4 + half_targ_off5;
    int off110011 = half_targ_off0 + half_targ_off1 + half_targ_off4 + half_targ_off5;
    int off110100 = half_targ_off2 + half_targ_off4 + half_targ_off5;
    int off110101 = half_targ_off0 + half_targ_off2 + half_targ_off4 + half_targ_off5;
    int off110111 = half_targ_off0 + half_targ_off1 + half_targ_off2 + half_targ_off4 + half_targ_off5;
    int off111000 = half_targ_off3 + half_targ_off4 + half_targ_off5;
    int off111001 = half_targ_off0 + half_targ_off3 + half_targ_off4 + half_targ_off5;
    int off111010 = half_targ_off1 + half_targ_off3 + half_targ_off4 + half_targ_off5;
    int off111011 = half_targ_off0 + half_targ_off1 + half_targ_off3 + half_targ_off4 + half_targ_off5;
    int off111100 = half_targ_off2 + half_targ_off3 + half_targ_off4 + half_targ_off5;
    int off111101 = half_targ_off0 + half_targ_off2 + half_targ_off3 + half_targ_off4 + half_targ_off5;
    int off111110 = half_targ_off1 + half_targ_off2 + half_targ_off3 + half_targ_off4 + half_targ_off5;

    complex<double> q000001;
    complex<double> q000010;
    complex<double> q000011;
    complex<double> q000100;
    complex<double> q000101;
    complex<double> q000110;
    complex<double> q000111;
    complex<double> q001000;
    complex<double> q001010;
    complex<double> q001011;
    complex<double> q001100;
    complex<double> q001101;
    complex<double> q001110;
    complex<double> q001111;
    complex<double> q010000;
    complex<double> q010001;
    complex<double> q010011;
    complex<double> q010100;
    complex<double> q010101;
    complex<double> q010110;
    complex<double> q010111;
    complex<double> q011000;
    complex<double> q011001;
    complex<double> q011010;
    complex<double> q011100;
    complex<double> q011101;
    complex<double> q011110;
    complex<double> q011111;
    complex<double> q100000;
    complex<double> q100001;
    complex<double> q100010;
    complex<double> q100011;
    complex<double> q100101;
    complex<double> q100110;
    complex<double> q100111;
    complex<double> q101000;
    complex<double> q101001;
    complex<double> q101010;
    complex<double> q101011;
    complex<double> q101100;
    complex<double> q101110;
    complex<double> q101111;
    complex<double> q110000;
    complex<double> q110001;
    complex<double> q110010;
    complex<double> q110011;
    complex<double> q110100;
    complex<double> q110101;
    complex<double> q110111;
    complex<double> q111000;
    complex<double> q111001;
    complex<double> q111010;
    complex<double> q111011;
    complex<double> q111100;
    complex<double> q111101;
    complex<double> q111110;

    for (int i = 0; i < env.chunk_state; i += targ_off2) {
        for (int j = 0; j < half_targ_off2; j += targ_off1) {
            for (int k = 0; k < half_targ_off1; k += half_targ_off0) {
                for (int m = 0; m < half_targ_off0; m++) {
                    q000001 = buffer[off000001];
                    q000010 = buffer[off000010];
                    q000011 = buffer[off000011];
                    q000100 = buffer[off000100];
                    q000101 = buffer[off000101];
                    q000110 = buffer[off000110];
                    q000111 = buffer[off000111];
                    q001000 = buffer[off001000];
                    q001010 = buffer[off001010];
                    q001011 = buffer[off001011];
                    q001100 = buffer[off001100];
                    q001101 = buffer[off001101];
                    q001110 = buffer[off001110];
                    q001111 = buffer[off001111];
                    q010000 = buffer[off010000];
                    q010001 = buffer[off010001];
                    q010011 = buffer[off010011];
                    q010100 = buffer[off010100];
                    q010101 = buffer[off010101];
                    q010110 = buffer[off010110];
                    q010111 = buffer[off010111];
                    q011000 = buffer[off011000];
                    q011001 = buffer[off011001];
                    q011010 = buffer[off011010];
                    q011100 = buffer[off011100];
                    q011101 = buffer[off011101];
                    q011110 = buffer[off011110];
                    q011111 = buffer[off011111];
                    q100000 = buffer[off100000];
                    q100001 = buffer[off100001];
                    q100010 = buffer[off100010];
                    q100011 = buffer[off100011];
                    q100101 = buffer[off100101];
                    q100110 = buffer[off100110];
                    q100111 = buffer[off100111];
                    q101000 = buffer[off101000];
                    q101001 = buffer[off101001];
                    q101010 = buffer[off101010];
                    q101011 = buffer[off101011];
                    q101100 = buffer[off101100];
                    q101110 = buffer[off101110];
                    q101111 = buffer[off101111];
                    q110000 = buffer[off110000];
                    q110001 = buffer[off110001];
                    q110010 = buffer[off110010];
                    q110011 = buffer[off110011];
                    q110100 = buffer[off110100];
                    q110101 = buffer[off110101];
                    q110111 = buffer[off110111];
                    q111000 = buffer[off111000];
                    q111001 = buffer[off111001];
                    q111010 = buffer[off111010];
                    q111011 = buffer[off111011];
                    q111100 = buffer[off111100];
                    q111101 = buffer[off111101];
                    q111110 = buffer[off111110];

                    buffer[off000001] = q001000;
                    buffer[off000010] = q010000;
                    buffer[off000011] = q011000;
                    buffer[off000100] = q100000;
                    buffer[off000101] = q101000;
                    buffer[off000110] = q110000;
                    buffer[off000111] = q111000;
                    buffer[off001000] = q000001;
                    buffer[off001010] = q010001;
                    buffer[off001011] = q011001;
                    buffer[off001100] = q100001;
                    buffer[off001101] = q101001;
                    buffer[off001110] = q110001;
                    buffer[off001111] = q111001;
                    buffer[off010000] = q000010;
                    buffer[off010001] = q001010;
                    buffer[off010011] = q011010;
                    buffer[off010100] = q100010;
                    buffer[off010101] = q101010;
                    buffer[off010110] = q110010;
                    buffer[off010111] = q111010;
                    buffer[off011000] = q000011;
                    buffer[off011001] = q001011;
                    buffer[off011010] = q010011;
                    buffer[off011100] = q100011;
                    buffer[off011101] = q101011;
                    buffer[off011110] = q110011;
                    buffer[off011111] = q111011;
                    buffer[off100000] = q000100;
                    buffer[off100001] = q001100;
                    buffer[off100010] = q010100;
                    buffer[off100011] = q011100;
                    buffer[off100101] = q101100;
                    buffer[off100110] = q110100;
                    buffer[off100111] = q111100;
                    buffer[off101000] = q000101;
                    buffer[off101001] = q001101;
                    buffer[off101010] = q010101;
                    buffer[off101011] = q011101;
                    buffer[off101100] = q100101;
                    buffer[off101110] = q110101;
                    buffer[off101111] = q111101;
                    buffer[off110000] = q000110;
                    buffer[off110001] = q001110;
                    buffer[off110010] = q010110;
                    buffer[off110011] = q011110;
                    buffer[off110100] = q100110;
                    buffer[off110101] = q101110;
                    buffer[off110111] = q111110;
                    buffer[off111000] = q000111;
                    buffer[off111001] = q001111;
                    buffer[off111010] = q010111;
                    buffer[off111011] = q011111;
                    buffer[off111100] = q100111;
                    buffer[off111101] = q101111;
                    buffer[off111110] = q110111;
                    
                    off000001++;
                    off000010++;
                    off000011++;
                    off000100++;
                    off000101++;
                    off000110++;
                    off000111++;
                    off001000++;
                    off001010++;
                    off001011++;
                    off001100++;
                    off001101++;
                    off001110++;
                    off001111++;
                    off010000++;
                    off010001++;
                    off010011++;
                    off010100++;
                    off010101++;
                    off010110++;
                    off010111++;
                    off011000++;
                    off011001++;
                    off011010++;
                    off011100++;
                    off011101++;
                    off011110++;
                    off011111++;
                    off100000++;
                    off100001++;
                    off100010++;
                    off100011++;
                    off100101++;
                    off100110++;
                    off100111++;
                    off101000++;
                    off101001++;
                    off101010++;
                    off101011++;
                    off101100++;
                    off101110++;
                    off101111++;
                    off110000++;
                    off110001++;
                    off110010++;
                    off110011++;
                    off110100++;
                    off110101++;
                    off110111++;
                    off111000++;
                    off111001++;
                    off111010++;
                    off111011++;
                    off111100++;
                    off111101++;
                    off111110++;
                }
                off000001 += half_targ_off0;
                off000010 += half_targ_off0;
                off000011 += half_targ_off0;
                off000100 += half_targ_off0;
                off000101 += half_targ_off0;
                off000110 += half_targ_off0;
                off000111 += half_targ_off0;
                off001000 += half_targ_off0;
                off001010 += half_targ_off0;
                off001011 += half_targ_off0;
                off001100 += half_targ_off0;
                off001101 += half_targ_off0;
                off001110 += half_targ_off0;
                off001111 += half_targ_off0;
                off010000 += half_targ_off0;
                off010001 += half_targ_off0;
                off010011 += half_targ_off0;
                off010100 += half_targ_off0;
                off010101 += half_targ_off0;
                off010110 += half_targ_off0;
                off010111 += half_targ_off0;
                off011000 += half_targ_off0;
                off011001 += half_targ_off0;
                off011010 += half_targ_off0;
                off011100 += half_targ_off0;
                off011101 += half_targ_off0;
                off011110 += half_targ_off0;
                off011111 += half_targ_off0;
                off100000 += half_targ_off0;
                off100001 += half_targ_off0;
                off100010 += half_targ_off0;
                off100011 += half_targ_off0;
                off100101 += half_targ_off0;
                off100110 += half_targ_off0;
                off100111 += half_targ_off0;
                off101000 += half_targ_off0;
                off101001 += half_targ_off0;
                off101010 += half_targ_off0;
                off101011 += half_targ_off0;
                off101100 += half_targ_off0;
                off101110 += half_targ_off0;
                off101111 += half_targ_off0;
                off110000 += half_targ_off0;
                off110001 += half_targ_off0;
                off110010 += half_targ_off0;
                off110011 += half_targ_off0;
                off110100 += half_targ_off0;
                off110101 += half_targ_off0;
                off110111 += half_targ_off0;
                off111000 += half_targ_off0;
                off111001 += half_targ_off0;
                off111010 += half_targ_off0;
                off111011 += half_targ_off0;
                off111100 += half_targ_off0;
                off111101 += half_targ_off0;
                off111110 += half_targ_off0;
            }
            off000001 += half_targ_off1;
            off000010 += half_targ_off1;
            off000011 += half_targ_off1;
            off000100 += half_targ_off1;
            off000101 += half_targ_off1;
            off000110 += half_targ_off1;
            off000111 += half_targ_off1;
            off001000 += half_targ_off1;
            off001010 += half_targ_off1;
            off001011 += half_targ_off1;
            off001100 += half_targ_off1;
            off001101 += half_targ_off1;
            off001110 += half_targ_off1;
            off001111 += half_targ_off1;
            off010000 += half_targ_off1;
            off010001 += half_targ_off1;
            off010011 += half_targ_off1;
            off010100 += half_targ_off1;
            off010101 += half_targ_off1;
            off010110 += half_targ_off1;
            off010111 += half_targ_off1;
            off011000 += half_targ_off1;
            off011001 += half_targ_off1;
            off011010 += half_targ_off1;
            off011100 += half_targ_off1;
            off011101 += half_targ_off1;
            off011110 += half_targ_off1;
            off011111 += half_targ_off1;
            off100000 += half_targ_off1;
            off100001 += half_targ_off1;
            off100010 += half_targ_off1;
            off100011 += half_targ_off1;
            off100101 += half_targ_off1;
            off100110 += half_targ_off1;
            off100111 += half_targ_off1;
            off101000 += half_targ_off1;
            off101001 += half_targ_off1;
            off101010 += half_targ_off1;
            off101011 += half_targ_off1;
            off101100 += half_targ_off1;
            off101110 += half_targ_off1;
            off101111 += half_targ_off1;
            off110000 += half_targ_off1;
            off110001 += half_targ_off1;
            off110010 += half_targ_off1;
            off110011 += half_targ_off1;
            off110100 += half_targ_off1;
            off110101 += half_targ_off1;
            off110111 += half_targ_off1;
            off111000 += half_targ_off1;
            off111001 += half_targ_off1;
            off111010 += half_targ_off1;
            off111011 += half_targ_off1;
            off111100 += half_targ_off1;
            off111101 += half_targ_off1;
            off111110 += half_targ_off1;
        }
        off000001 += half_targ_off2;
        off000010 += half_targ_off2;
        off000011 += half_targ_off2;
        off000100 += half_targ_off2;
        off000101 += half_targ_off2;
        off000110 += half_targ_off2;
        off000111 += half_targ_off2;
        off001000 += half_targ_off2;
        off001010 += half_targ_off2;
        off001011 += half_targ_off2;
        off001100 += half_targ_off2;
        off001101 += half_targ_off2;
        off001110 += half_targ_off2;
        off001111 += half_targ_off2;
        off010000 += half_targ_off2;
        off010001 += half_targ_off2;
        off010011 += half_targ_off2;
        off010100 += half_targ_off2;
        off010101 += half_targ_off2;
        off010110 += half_targ_off2;
        off010111 += half_targ_off2;
        off011000 += half_targ_off2;
        off011001 += half_targ_off2;
        off011010 += half_targ_off2;
        off011100 += half_targ_off2;
        off011101 += half_targ_off2;
        off011110 += half_targ_off2;
        off011111 += half_targ_off2;
        off100000 += half_targ_off2;
        off100001 += half_targ_off2;
        off100010 += half_targ_off2;
        off100011 += half_targ_off2;
        off100101 += half_targ_off2;
        off100110 += half_targ_off2;
        off100111 += half_targ_off2;
        off101000 += half_targ_off2;
        off101001 += half_targ_off2;
        off101010 += half_targ_off2;
        off101011 += half_targ_off2;
        off101100 += half_targ_off2;
        off101110 += half_targ_off2;
        off101111 += half_targ_off2;
        off110000 += half_targ_off2;
        off110001 += half_targ_off2;
        off110010 += half_targ_off2;
        off110011 += half_targ_off2;
        off110100 += half_targ_off2;
        off110101 += half_targ_off2;
        off110111 += half_targ_off2;
        off111000 += half_targ_off2;
        off111001 += half_targ_off2;
        off111010 += half_targ_off2;
        off111011 += half_targ_off2;
        off111100 += half_targ_off2;
        off111101 += half_targ_off2;
        off111110 += half_targ_off2;
    }
}

VSWAP_Gate_4_4::VSWAP_Gate_4_4(vector<int> targ): Gate(targ) {
    type = VSWAP;
    swap_out = vector<int>(targ.begin()  , targ.begin()+4);
    swap_in  = vector<int>(targ.begin()+4, targ.end());

    for (int i = 0; i < 4; i++)
        half_targ_off.push_back(env.qubit_offset[swap_out[i]]);
    
    half_targ_off.push_back(env.chunk_state);
    for (int i = 0; i < 3; i++)
        half_targ_off.push_back(half_targ_off.back()*2);
    
    for (int i = 0; i < 8; i++)
        targ_off.push_back(half_targ_off[i]*2);
    
    for (int i = 0; i < 256; i++){
        int off = 0;
        for (int j = 0; j < 8; j++){
            if (i & (1 << j))
                off += half_targ_off[j];
        }
        init_off[i] = off;
    }

    for (int i = 0; i < 256; i++) {
        if ((i>>4) != (i%16)) {
            int n1 = i;
            int n2 = ((i%16)<<4)|(i>>4);
            if (n1 < n2) {
                swap_list.push_back(n1);
                dest_list.push_back(n2);
            }
        }
    }

    run = bind(&VSWAP_Gate_4_4::run_nonchunk_chunk, this, _1);
    name += "VSWAP_Gate_4_4";
}

void VSWAP_Gate_4_4::run_nonchunk_chunk(vector<complex<double>> &buffer){
    array<complex<double>, 256> q;
    array<int, 256> off;
    off = init_off;

    for (int i = 0; i < env.chunk_state; i += targ_off[3]) {
        for (int j = 0; j < half_targ_off[3]; j += targ_off[2]) {
            for (int k = 0; k < half_targ_off[2]; k += targ_off[1]) {
                for (int m = 0; m < half_targ_off[1]; m += targ_off[0]) {
                    for (int n = 0; n < half_targ_off[0]; n++) {
                        for (auto p : swap_list) {
                            q[p] = buffer[off[p]];
                        }
                        for (int p = 0; p < 120; p++) {
                            complex<double> temp = buffer[swap_list[p]];
                            buffer[swap_list[p]] = buffer[dest_list[p]];
                            buffer[dest_list[p]] = temp;
                        }
                        for (auto p : off) {
                            p++;
                        }
                    }
                    for (auto p : off) {
                        p += half_targ_off[0];
                    }
                }
                for (auto p : off) {
                    p += half_targ_off[1];
                }
            }
            for (auto p : off) {
                p += half_targ_off[2];
            }
        }
        for (auto p : off) {
            p += half_targ_off[3];
        }
    }
}

VSWAP_Gate_6_6::VSWAP_Gate_6_6(vector<int> targ): Gate(targ) {
    type = VSWAP;
    swap_out = vector<int>(targ.begin()  , targ.begin()+6);
    swap_in  = vector<int>(targ.begin()+6, targ.end());

    for (int i = 0; i < 6; i++)
        half_targ_off.push_back(env.qubit_offset[swap_out[i]]);
    
    half_targ_off.push_back(env.chunk_state);
    for (int i = 0; i < 5; i++)
        half_targ_off.push_back(half_targ_off.back()*2);
    
    for (int i = 0; i < 12; i++)
        targ_off.push_back(half_targ_off[i]*2);
    
    for (int i = 0; i < 4096; i++){
        int off = 0;
        for (int j = 0; j < 12; j++){
            if (i & (1 << j))
                off += half_targ_off[j];
        }
        init_off[i] = off;
    }

    for (int i = 0; i < 4096; i++){
        if ((i>>6) != (i%64)) {
            int n1 = i;
            int n2 = ((i%64)<<6)|(i>>6);
            if (n1 < n2) {
                swap_list.push_back(n1);
                dest_list.push_back(n2);
            }
        }
    }

    run = bind(&VSWAP_Gate_6_6::run_nonchunk_chunk, this, _1);
    name += "VSWAP_Gate_6_6";
}

void VSWAP_Gate_6_6::run_nonchunk_chunk(vector<complex<double>> &buffer){
    array<complex<double>, 4096> q;
    array<int, 4096> off;
    off = init_off;

    for (int i = 0; i < env.chunk_state; i += targ_off[5]) {
        for (int j = 0; j < half_targ_off[5]; j += targ_off[4]) {
            for (int k = 0; k < half_targ_off[4]; k += targ_off[3]) {
                for (int m = 0; m < half_targ_off[3]; m += targ_off[2]) {
                    for (int n = 0; n < half_targ_off[2]; n += targ_off[1]) {
                        for (int r = 0; r < half_targ_off[1]; r += targ_off[0]) {
                            for (int s = 0; s < half_targ_off[0]; s++) {
                                for (auto p : swap_list) {
                                    q[p] = buffer[off[p]];
                                }
                                for (int p = 0; p < 2016; p++) { // (2^12 - 2^6) / 2 = 2016
                                    complex<double> temp = buffer[swap_list[p]];
                                    buffer[swap_list[p]] = buffer[dest_list[p]];
                                    buffer[dest_list[p]] = temp;
                                }
                                for (auto p : off) {
                                    p++;
                                }
                            }
                            for (auto p : off) {
                                p += half_targ_off[0];
                            }
                        }
                        for (auto p : off) {
                            p += half_targ_off[1];
                        }
                    }
                    for (auto p : off) {
                        p += half_targ_off[2];
                    }
                }
                for (auto p : off) {
                    p += half_targ_off[3];
                }
            }
            for (auto p : off) {
                p += half_targ_off[4];
            }
        }
        for (auto p : off) {
            p += half_targ_off[5];
        }
    }
}

/*DIO version*/
/*-----ONE_QUBIT_GATE-----*/
H_Gate_DIO::H_Gate_DIO(vector<int> targ): ONE_QUBIT_GATE(targ) {
    name = "H_Gate_DIO";
    if(mpi_count)
    {
        bind_gate_mpi_1_dio(H_Gate_DIO)
    }
    else
    {
        bind_gate_1_dio(H_Gate_DIO)
    }
}

void H_Gate_DIO::run_chunk(complex<double> *buffer){
    chunk_gate(init_off2(int, 0, half_off), update_off2(1), update_off2(half_off), HGATE)
}

void H_Gate_DIO::run_nonchunk(complex<double> *buffer){
    nonchunk_gate(init_off2(int, 0, half_off), update_off2(1), HGATE)
}

void H_Gate_DIO::run_mpi_chunk(complex<double> *buffer1,complex<double> *buffer2,int offset0,int offset1,int round){
    mpi_gate(init_off2(int, offset0 * env.chunk_state, offset1 * env.chunk_state),round, update_off2(1), HGATEMPI)
}

X_Gate_DIO::X_Gate_DIO(vector<int> targ): ONE_QUBIT_GATE(targ) {
    name = "X_Gate_DIO";
    if(mpi_count)
    {
        bind_gate_mpi_1_dio(X_Gate_DIO)
    }
    else
    {
        bind_gate_1_dio(X_Gate_DIO)
    }
}

void X_Gate_DIO::run_chunk(complex<double> *buffer){
    chunk_gate(init_off2(int, 0, half_off), update_off2(1), update_off2(half_off), XGATE)
}

void X_Gate_DIO::run_nonchunk(complex<double> *buffer){
    nonchunk_gate(init_off2(int, 0, half_off), update_off2(1), XGATE)
}

void X_Gate_DIO::run_mpi_chunk(complex<double> *buffer1,complex<double> *buffer2,int offset0,int offset1,int round){
    mpi_gate(init_off2(int, offset0 * env.chunk_state, offset1 * env.chunk_state),round, update_off2(1), XGATEMPI)
}

Y_Gate_DIO::Y_Gate_DIO(vector<int> targ): ONE_QUBIT_GATE(targ) {
    name = "Y_Gate_DIO";
    if(mpi_count)
    {
        bind_gate_mpi_1_dio(Y_Gate_DIO)
    }
    else
    {
        bind_gate_1_dio(Y_Gate_DIO)
    }
}

void Y_Gate_DIO::run_chunk(complex<double> *buffer){
    chunk_gate(init_off2(int, 0, half_off), update_off2(1), update_off2(half_off), YGATE)
}

void Y_Gate_DIO::run_nonchunk(complex<double> *buffer){
    nonchunk_gate(init_off2(int, 0, half_off), update_off2(1), YGATE)
}

void Y_Gate_DIO::run_mpi_chunk(complex<double> *buffer1,complex<double> *buffer2,int offset0,int offset1,int round){
    mpi_gate(init_off2(int, offset0 * env.chunk_state, offset1 * env.chunk_state),round, update_off2(1), YGATEMPI)
}

Z_Gate_DIO::Z_Gate_DIO(vector<int> targ): ONE_QUBIT_GATE(targ) {
    name = "Z_Gate_DIO";
    bind_gate_1_dio(Z_Gate_DIO)
}

void Z_Gate_DIO::run_chunk(complex<double> *buffer){
    chunk_gate(init_off1(int, half_off), update_off1(1), update_off1(half_off), ZGATE)
}

void Z_Gate_DIO::run_nonchunk(complex<double> *buffer){
    nonchunk_gate(init_off1(int, half_off), update_off1(1), ZGATE)
}

U1_Gate_DIO::U1_Gate_DIO(vector<int> targ, vector<complex<double>> coeff): ONE_QUBIT_GATE(targ), coeff(coeff) {
    name = "U1_Gate_DIO";
    if(mpi_count)
    {
        bind_gate_mpi_1_dio(U1_Gate_DIO)
    }
    else
    {
        bind_gate_1_dio(U1_Gate_DIO)
    }
}

void U1_Gate_DIO::run_chunk(complex<double> *buffer){
    chunk_gate(init_off2(int, 0, half_off), update_off2(1), update_off2(half_off), U1GATE)
}

void U1_Gate_DIO::run_nonchunk(complex<double> *buffer){
    nonchunk_gate(init_off2(int, 0, half_off), update_off2(1), U1GATE)
}

void U1_Gate_DIO::run_mpi_chunk(complex<double> *buffer1,complex<double> *buffer2,int offset0,int offset1,int round){
    mpi_gate(init_off2(int, offset0 * env.chunk_state, offset1 * env.chunk_state),round, update_off2(1), U1GATEMPI)
}

Phase_Gate_DIO::Phase_Gate_DIO(vector<int> targ, double phi): ONE_QUBIT_GATE(targ), phi(cos(phi), sin(phi)) {
    name = "Phase_Gate_DIO";
    bind_gate_1_dio(Phase_Gate_DIO)
}

void Phase_Gate_DIO::run_chunk(complex<double> *buffer){
    chunk_gate(init_off1(int, half_off), update_off1(1), update_off1(half_off), PGATE)
}

void Phase_Gate_DIO::run_nonchunk(complex<double> *buffer){
    nonchunk_gate(init_off1(int, half_off), update_off1(1), PGATE)
}

RX_Gate_DIO::RX_Gate_DIO(vector<int> targ, double phi): ONE_QUBIT_GATE(targ), cos_Phi_2(cos(phi/2), 0), i_sin_Phi_2(0, -sin(phi/2)) {
    name = "RX_Gate_DIO";
    if(mpi_count)
    {
        bind_gate_mpi_1_dio(RX_Gate_DIO)
    }
    else
    {
        bind_gate_1_dio(RX_Gate_DIO)
    }
}

void RX_Gate_DIO::run_chunk(complex<double> *buffer){
    chunk_gate(init_off2(int, 0, half_off), update_off2(1), update_off2(half_off), RXGATE)
}

void RX_Gate_DIO::run_nonchunk(complex<double> *buffer){
    nonchunk_gate(init_off2(int, 0, half_off), update_off2(1), RXGATE)
}

void RX_Gate_DIO::run_mpi_chunk(complex<double> *buffer1,complex<double> *buffer2,int offset0,int offset1,int round){
    mpi_gate(init_off2(int, offset0 * env.chunk_state, offset1 * env.chunk_state),round, update_off2(1), RXGATEMPI)
}

RY_Gate_DIO::RY_Gate_DIO(vector<int> targ, double phi): ONE_QUBIT_GATE(targ), cos_Phi_2(cos(phi/2), 0), sin_Phi_2(sin(phi/2), 0) {
    name = "RY_Gate_DIO";
    if(mpi_count)
    {
        bind_gate_mpi_1_dio(RY_Gate_DIO)
    }
    else
    {
        bind_gate_1_dio(RY_Gate_DIO)
    }
}

void RY_Gate_DIO::run_chunk(complex<double> *buffer){
    chunk_gate(init_off2(int, 0, half_off), update_off2(1), update_off2(half_off), RYGATE)
}

void RY_Gate_DIO::run_nonchunk(complex<double> *buffer){
    nonchunk_gate(init_off2(int, 0, half_off), update_off2(1), RYGATE)
}

void RY_Gate_DIO::run_mpi_chunk(complex<double> *buffer1,complex<double> *buffer2,int offset0,int offset1,int round){
    mpi_gate(init_off2(int, offset0 * env.chunk_state, offset1 * env.chunk_state),round, update_off2(1), RYGATEMPI)
}

RZ_Gate_DIO::RZ_Gate_DIO(vector<int> targ, double phi): ONE_QUBIT_GATE(targ), exp_p_iPhi_2(cos(phi/2), sin(phi/2)), exp_n_iPhi_2(cos(-phi/2), sin(-phi/2)) {
    name = "RZ_Gate_DIO";
    bind_gate_1_dio(RZ_Gate_DIO)
}

void RZ_Gate_DIO::run_chunk(complex<double> *buffer){
    chunk_gate(init_off2(int, 0, half_off), update_off2(1), update_off2(half_off), RZGATE)
}

void RZ_Gate_DIO::run_nonchunk(complex<double> *buffer){
    nonchunk_gate(init_off2(int, 0, half_off), update_off2(1), RZGATE)
}

/*-----TWO_QUBIT_GATE-----*/

U2_Gate_DIO::U2_Gate_DIO(vector<int> targ, vector<complex<double>> matrix): TWO_QUBIT_GATE(targ){
    name = "U2_Gate_DIO";

    assert(targ[0] < targ[1]);
    for (int i = 0; i < 16; i++) {
        int res = 0;
        int b1 = i & 1;
        int b2 = (i >> 1) & 1;
        int b3 = (i >> 2) & 1;
        int b4 = (i >> 3) & 1;
        res = (b3 << 3) | (b4 << 2) | (b1 << 1) | b2;
        coeff[i] = matrix[res];
    }
    if(mpi_count == 1 && chunk_count == 1)
    {
        bind_gate_mpi_1_dio(U2_Gate_DIO)
    }
    else
    {
        bind_gate_2_dio(U2_Gate_DIO)
    }
}

void U2_Gate_DIO::run_chunk_chunk(complex<double> *buffer){
    chunk_chunk_gate(init_off4(int, 0, half_off_0, half_off_1, (half_off_0 + half_off_1)),
                    update_off4(1), 
                    update_off4(half_off_0), 
                    update_off4(half_off_1), 
                    U2Gate)
}

void U2_Gate_DIO::run_nonchunk_chunk(complex<double> *buffer){
    nonchunk_chunk_gate(init_off4(int, 0, half_off_0, half_off_1, (half_off_0 + half_off_1)),
                    update_off4(1), 
                    update_off4(half_off_0), 
                    U2Gate)
}

void U2_Gate_DIO::run_nonchunk_nonchunk(complex<double> *buffer){
    nonchunk_nonchunk_gate(init_off4(int, 0, half_off_0, half_off_1, (half_off_0 + half_off_1)),
                    update_off4(1), 
                    U2Gate)
}

void U2_Gate_DIO::run_mpi_chunk(complex<double>*buffer1,complex<double>*buffer2,int offset0,int offset1,int round){
    mpi_nonchunk_chunk(init_off4(int,offset0 * env.chunk_state,offset0 * env.chunk_state + half_off_0,offset1 * env.chunk_state, offset1 * env.chunk_state + half_off_0),
                    round,
                    update_off4(1),
                    update_off4(half_off_0),
                    U2MPICHUNK)
}

SWAP_Gate_DIO::SWAP_Gate_DIO(vector<int> targ): TWO_QUBIT_GATE(targ) {
    name = "SWAP_Gate_DIO";
    if(mpi_count == 1 && chunk_count == 1)
    {
        bind_gate_special_swap_dio(SWAP_Gate_DIO)
    }
    else
    {
        bind_gate_2_dio(SWAP_Gate_DIO)
    }
}

MPI_VSWAP_Gate_1_1_DIO::MPI_VSWAP_Gate_1_1_DIO(vector<int> targ): SWAP_Gate_DIO(targ) {
    type = VSWAP;
    name = "MPI_VSWAP_Gate_1_1_DIO";
}

void SWAP_Gate_DIO::run_chunk_chunk(complex<double> *buffer){
    chunk_chunk_gate(init_off2(int, half_off_0, half_off_1),
                    update_off2(1), 
                    update_off2(half_off_0), 
                    update_off2(half_off_1), 
                    SWAPGate)
}

void SWAP_Gate_DIO::run_nonchunk_chunk(complex<double> *buffer){
    nonchunk_chunk_gate(init_off2(int, half_off_0, half_off_1),
                    update_off2(1), 
                    update_off2(half_off_0), 
                    SWAPGate)
}

void SWAP_Gate_DIO::run_nonchunk_nonchunk(complex<double> *buffer){
    nonchunk_nonchunk_gate(init_off2(int, half_off_0, half_off_1),
                    update_off2(1), 
                    SWAPGate)
}

void SWAP_Gate_DIO::run_mpi_nonchunk_chunk(complex<double>*buffer1,complex<double>*buffer2,int offset0,int offset1,int round){
    mpi_nonchunk_chunk(init_off2(int, offset0 * env.chunk_state + half_off_0, offset1 * env.chunk_state),round, update_off2(1),update_off2(half_off_0),SWAPGateMPIChunk)
}

VSWAP_Gate_1_1_DIO::VSWAP_Gate_1_1_DIO(vector<int> targ): SWAP_Gate_DIO(targ) {
    type = VSWAP;
    name = "VSWAP_Gate_1_1_DIO";
}

CPhase_Gate_DIO::CPhase_Gate_DIO(vector<int> targ, double phi): TWO_QUBIT_GATE(targ) {
    name = "CPhase_Gate_DIO";
    exp_iPhi = complex<double>(cos(phi), sin(phi));
    bind_gate_2_dio(CPhase_Gate_DIO)
}

void CPhase_Gate_DIO::run_chunk_chunk(complex<double> *buffer){
    chunk_chunk_gate(init_off1(int, (half_off_0 + half_off_1)),
                    update_off1(1), 
                    update_off1(half_off_0), 
                    update_off1(half_off_1), 
                    CPGATE)
}

void CPhase_Gate_DIO::run_nonchunk_chunk(complex<double> *buffer){
    nonchunk_chunk_gate(init_off1(int, (half_off_0 + half_off_1)),
                    update_off1(1), 
                    update_off1(half_off_0), 
                    CPGATE)
}

void CPhase_Gate_DIO::run_nonchunk_nonchunk(complex<double> *buffer){
    nonchunk_nonchunk_gate(init_off1(int, (half_off_0 + half_off_1)),
                    update_off1(1), 
                    CPGATE)
}

RZZ_Gate_DIO::RZZ_Gate_DIO(vector<int> targ, double phi): TWO_QUBIT_GATE(targ) {
    name += "RZZ_Gate_DIO";
    exp_p_iPhi_2 = complex<double>(cos(phi/2), sin(phi/2));
    exp_n_iPhi_2 = complex<double>(cos(-phi/2), sin(-phi/2));
    bind_gate_2_dio(RZZ_Gate_DIO)
}

void RZZ_Gate_DIO::run_chunk_chunk(complex<double> *buffer){
    chunk_chunk_gate(init_off4(int, 0, half_off_0, half_off_1, (half_off_0 + half_off_1)),
                    update_off4(1), 
                    update_off4(half_off_0), 
                    update_off4(half_off_1), 
                    RZZGATE)
}

void RZZ_Gate_DIO::run_nonchunk_chunk(complex<double> *buffer){
    nonchunk_chunk_gate(init_off4(int, 0, half_off_0, half_off_1, (half_off_0 + half_off_1)),
                    update_off4(1), 
                    update_off4(half_off_0), 
                    RZZGATE)
}

void RZZ_Gate_DIO::run_nonchunk_nonchunk(complex<double> *buffer){
    nonchunk_nonchunk_gate(init_off4(int, 0, half_off_0, half_off_1, (half_off_0 + half_off_1)),
                    update_off4(1), 
                    RZZGATE)
}

VSWAP_Gate_2_2_DIO::VSWAP_Gate_2_2_DIO(vector<int> targ): Gate(targ) {
    type = VSWAP;
    name = "VSWAP_Gate_2_2_DIO";

    half_off_0 = env.qubit_offset[targ[0]];
    half_off_1 = env.qubit_offset[targ[1]];
    half_off_2 = env.chunk_state;
    half_off_3 = 2 * env.chunk_state;

    off_0 = half_off_0 * 2;
    off_1 = half_off_1 * 2;
    off_2 = half_off_2 * 2;
    off_3 = half_off_3 * 2;
        
    run_dio = bind(&VSWAP_Gate_2_2_DIO::run_nonchunk_chunk, this, _1);
}

void VSWAP_Gate_2_2_DIO::run_nonchunk_chunk(complex<double> *buffer){
    int off0001 = half_off_0;
    int off0010 = half_off_1;
    int off0011 = half_off_0 + half_off_1;
    int off0100 = half_off_2;
    int off0110 = half_off_1 + half_off_2;
    int off0111 = half_off_0 + half_off_1 + half_off_2;
    int off1000 = half_off_3;
    int off1001 = half_off_0 + half_off_3;
    int off1011 = half_off_0 + half_off_1 + half_off_3;
    int off1100 = half_off_2 + half_off_3;
    int off1101 = half_off_0 + half_off_2 + half_off_3;
    int off1110 = half_off_1 + half_off_2 + half_off_3;

    complex<double> q0001;
    complex<double> q0010;
    complex<double> q0011;
    complex<double> q0100;
    complex<double> q0110;
    complex<double> q0111;
    complex<double> q1000;
    complex<double> q1001;
    complex<double> q1011;
    complex<double> q1100;
    complex<double> q1101;
    complex<double> q1110;

    for (int i = 0; i < env.chunk_state; i += off_1) {
        for (int j = 0; j < half_off_1; j += off_0) {
            for (int k = 0; k < half_off_0; k++) {
                q0001 = buffer[off0001];
                q0010 = buffer[off0010];
                q0011 = buffer[off0011];
                q0100 = buffer[off0100];
                q0110 = buffer[off0110];
                q0111 = buffer[off0111];
                q1000 = buffer[off1000];
                q1001 = buffer[off1001];
                q1011 = buffer[off1011];
                q1100 = buffer[off1100];
                q1101 = buffer[off1101];
                q1110 = buffer[off1110];
                buffer[off0001] = q0100;
                buffer[off0010] = q1000;
                buffer[off0011] = q1100;
                buffer[off0100] = q0001;
                buffer[off0110] = q1001;
                buffer[off0111] = q1101;
                buffer[off1000] = q0010;
                buffer[off1001] = q0110;
                buffer[off1011] = q1110;
                buffer[off1100] = q0011;
                buffer[off1101] = q0111;
                buffer[off1110] = q1011;
                off0001++;
                off0010++;
                off0011++;
                off0100++;
                off0110++;
                off0111++;
                off1000++;
                off1001++;
                off1011++;
                off1100++;
                off1101++;
                off1110++;
            }
            off0001 += half_off_0;
            off0010 += half_off_0;
            off0011 += half_off_0;
            off0100 += half_off_0;
            off0110 += half_off_0;
            off0111 += half_off_0;
            off1000 += half_off_0;
            off1001 += half_off_0;
            off1011 += half_off_0;
            off1100 += half_off_0;
            off1101 += half_off_0;
            off1110 += half_off_0;
        }
        off0001 += half_off_1;
        off0010 += half_off_1;
        off0011 += half_off_1;
        off0100 += half_off_1;
        off0110 += half_off_1;
        off0111 += half_off_1;
        off1000 += half_off_1;
        off1001 += half_off_1;
        off1011 += half_off_1;
        off1100 += half_off_1;
        off1101 += half_off_1;
        off1110 += half_off_1;
    }
}

MPI_VSWAP_Gate_2_2_DIO::MPI_VSWAP_Gate_2_2_DIO(vector<int> targ): Gate(targ) {
    type = VSWAP;
    name = "MPI_VSWAP_Gate_2_2_DIO";
    half_off_0 = env.qubit_offset[targ[0]];
    half_off_1 = env.qubit_offset[targ[1]];
    half_off_2 = 0;
    half_off_3 = 0;
        
    off_0 = half_off_0 * 2;
    off_1 = half_off_1 * 2;
    off_2 = half_off_2 * 2;
    off_3 = half_off_3 * 2;
    run_mpi_vswap2_2_dio = bind(&MPI_VSWAP_Gate_2_2_DIO::run_vswap_2_2,this,_1,_2,_3,_4,_5);
}

void MPI_VSWAP_Gate_2_2_DIO::run_vswap_2_2(complex<double> *buffer1,complex<double> *buffer2,complex<double> *buffer3,complex<double> *buffer4,int round)
{
    int off0001 = half_off_0;
    int off0010 = half_off_1;
    int off0011 = half_off_0 + half_off_1;
    int off0100 = half_off_2;
    int off0110 = half_off_1 + half_off_2;
    int off0111 = half_off_0 + half_off_1 + half_off_2;
    int off1000 = half_off_3;
    int off1001 = half_off_0 + half_off_3;
    int off1011 = half_off_0 + half_off_1 + half_off_3;
    int off1100 = half_off_2 + half_off_3;
    int off1101 = half_off_0 + half_off_2 + half_off_3;
    int off1110 = half_off_1 + half_off_2 + half_off_3;

    complex<double> q0001;
    complex<double> q0010;
    complex<double> q0011;
    complex<double> q0100;
    complex<double> q0110;
    complex<double> q0111;
    complex<double> q1000;
    complex<double> q1001;
    complex<double> q1011;
    complex<double> q1100;
    complex<double> q1101;
    complex<double> q1110;

    for (int i = 0; i < round * env.chunk_state; i += off_1) {
        for (int j = 0; j < half_off_1; j += off_0) {
            for (int k = 0; k < half_off_0; k++) {
                q0001 = buffer1[off0001];
                q0010 = buffer1[off0010];
                q0011 = buffer1[off0011];
                q0100 = buffer2[off0100];
                q0110 = buffer2[off0110];
                q0111 = buffer2[off0111];
                q1000 = buffer3[off1000];
                q1001 = buffer3[off1001];
                q1011 = buffer3[off1011];
                q1100 = buffer4[off1100];
                q1101 = buffer4[off1101];
                q1110 = buffer4[off1110];
                buffer1[off0001] = q0100;
                buffer1[off0010] = q1000;
                buffer1[off0011] = q1100;
                buffer2[off0100] = q0001;
                buffer2[off0110] = q1001;
                buffer2[off0111] = q1101;
                buffer3[off1000] = q0010;
                buffer3[off1001] = q0110;
                buffer3[off1011] = q1110;
                buffer4[off1100] = q0011;
                buffer4[off1101] = q0111;
                buffer4[off1110] = q1011;
                off0001++;
                off0010++;
                off0011++;
                off0100++;
                off0110++;
                off0111++;
                off1000++;
                off1001++;
                off1011++;
                off1100++;
                off1101++;
                off1110++;
            }
            off0001 += half_off_0;
            off0010 += half_off_0;
            off0011 += half_off_0;
            off0100 += half_off_0;
            off0110 += half_off_0;
            off0111 += half_off_0;
            off1000 += half_off_0;
            off1001 += half_off_0;
            off1011 += half_off_0;
            off1100 += half_off_0;
            off1101 += half_off_0;
            off1110 += half_off_0;
        }
        off0001 += half_off_1;
        off0010 += half_off_1;
        off0011 += half_off_1;
        off0100 += half_off_1;
        off0110 += half_off_1;
        off0111 += half_off_1;
        off1000 += half_off_1;
        off1001 += half_off_1;
        off1011 += half_off_1;
        off1100 += half_off_1;
        off1101 += half_off_1;
        off1110 += half_off_1;
    }
}

/*MEM version*/
/*-----ONE_QUBIT_GATE-----*/
ONE_QUBIT_GATE_MEM::ONE_QUBIT_GATE_MEM(vector<int> targ): Gate(targ) {
    type = ONE_QUBIT;
    half_off = env.qubit_offset[targ[0]];
    off = half_off * 2;
}

H_Gate_MEM::H_Gate_MEM(vector<int> targ): ONE_QUBIT_GATE_MEM(targ) {
    name = "H_Gate_MEM";
    if(mpi_count)
    {
        bind_gate_mpi_1_mem(H_Gate_MEM)
    }
    else
    {
        bind_gate_1_mem(H_Gate_MEM)
    }
}

void H_Gate_MEM::run_chunk(vector<complex<double>> &buffer, long long idx){
    chunk_gate_mem(init_off2(long long, idx, idx + half_off), 
                update_off2(1), update_off2(half_off), HGATE)
}

void H_Gate_MEM::run_nonchunk(vector<complex<double>> &buffer, long long idx){
    nonchunk_gate_mem(init_off2(long long, idx, idx + half_off), 
                    update_off2(1), HGATE)
}

void H_Gate_MEM::run_mpi_mem(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,long long offset0,long long offset1,long long length){
    mpi_gate_mem(init_off2(long long, offset0, offset1),length, update_off2(1), HGATEMPI)
}

X_Gate_MEM::X_Gate_MEM(vector<int> targ): ONE_QUBIT_GATE_MEM(targ) {
    name = "X_Gate_MEM";
    if(mpi_count)
    {
        bind_gate_mpi_1_mem(X_Gate_MEM)
    }
    else
    {
        bind_gate_1_mem(X_Gate_MEM)
    }
}

void X_Gate_MEM::run_chunk(vector<complex<double>> &buffer, long long idx){
    chunk_gate_mem(init_off2(long long, idx, idx + half_off), 
                update_off2(1), update_off2(half_off), XGATE)
}

void X_Gate_MEM::run_nonchunk(vector<complex<double>> &buffer, long long idx){
    nonchunk_gate_mem(init_off2(long long, idx, idx + half_off), 
                    update_off2(1), XGATE)
}

void X_Gate_MEM::run_mpi_mem(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,long long offset0,long long offset1,long long length){
    mpi_gate_mem(init_off2(long long, offset0, offset1),length, update_off2(1), XGATEMPI)
}

Y_Gate_MEM::Y_Gate_MEM(vector<int> targ): ONE_QUBIT_GATE_MEM(targ) {
    name = "Y_Gate_MEM";
    if(mpi_count)
    {
        bind_gate_mpi_1_mem(Y_Gate_MEM)
    }
    else
    {
        bind_gate_1_mem(Y_Gate_MEM)
    }
}

void Y_Gate_MEM::run_chunk(vector<complex<double>> &buffer, long long idx){
    chunk_gate_mem(init_off2(long long, idx, idx + half_off), 
                update_off2(1), update_off2(half_off), YGATE)
}

void Y_Gate_MEM::run_nonchunk(vector<complex<double>> &buffer, long long idx){
    nonchunk_gate_mem(init_off2(long long, idx, idx + half_off), 
                    update_off2(1), YGATE)
}

void Y_Gate_MEM::run_mpi_mem(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,long long offset0,long long offset1,long long length){
    mpi_gate_mem(init_off2(long long, offset0, offset1),length, update_off2(1), YGATEMPI)
}

Z_Gate_MEM::Z_Gate_MEM(vector<int> targ): ONE_QUBIT_GATE_MEM(targ) {
    name = "Z_Gate_MEM";
    if(mpi_count)
    {
        bind_gate_mpi_1_mem_diagonal(Z_Gate_MEM)
    }
    else
    {
        bind_gate_1_mem(Z_Gate_MEM)
    }
}

void Z_Gate_MEM::run_chunk(vector<complex<double>> &buffer, long long idx){
    chunk_gate_mem(init_off1(long long, idx + half_off), 
                update_off1(1), update_off1(half_off), ZGATE)
}

void Z_Gate_MEM::run_nonchunk(vector<complex<double>> &buffer, long long idx){
    nonchunk_gate_mem(init_off1(long long, idx + half_off), 
                    update_off1(1), ZGATE)
}

void Z_Gate_MEM::run_mpi_mem(vector<complex<double>> &buffer, long long idx,int pos,long long length){
    if(!pos) return;
    mpi_gate_mem(init_off1(long long, idx),length, update_off1(1), ZGATE)
}

Phase_Gate_MEM::Phase_Gate_MEM(vector<int> targ, double phi): ONE_QUBIT_GATE_MEM(targ), phi(cos(phi), sin(phi)) {
    name = "Phase_Gate_MEM";
    if(mpi_count)
    {
        bind_gate_mpi_1_mem_diagonal(Phase_Gate_MEM)
    }
    else
    {
        bind_gate_1_mem(Phase_Gate_MEM)
    }
}

void Phase_Gate_MEM::run_chunk(vector<complex<double>> &buffer, long long idx){
    chunk_gate_mem(init_off1(long long, idx + half_off), 
                update_off1(1), update_off1(half_off), PGATE)
}

void Phase_Gate_MEM::run_nonchunk(vector<complex<double>> &buffer, long long idx){
    nonchunk_gate_mem(init_off1(long long, idx + half_off), 
                    update_off1(1), PGATE)
}
void Phase_Gate_MEM::run_mpi_mem(vector<complex<double>> &buffer, long long idx,int pos,long long length){
    if(!pos) return;
    mpi_gate_mem(init_off1(long long, idx),length, update_off1(1), PGATE)
}
U1_Gate_MEM::U1_Gate_MEM(vector<int> targ, vector<complex<double>> coeff): ONE_QUBIT_GATE_MEM(targ), coeff(coeff) {
    name = "U1_Gate_MEM";
    if(mpi_count)
    {
        bind_gate_mpi_1_mem(U1_Gate_MEM)
    }
    else
    {
        bind_gate_1_mem(U1_Gate_MEM)
    }
}

void U1_Gate_MEM::run_chunk(vector<complex<double>> &buffer, long long idx){
    chunk_gate_mem(init_off2(long long, idx, idx + half_off), 
                update_off2(1), update_off2(half_off), U1GATE)
}

void U1_Gate_MEM::run_nonchunk(vector<complex<double>> &buffer, long long idx){
    nonchunk_gate_mem(init_off2(long long, idx, idx + half_off), 
                    update_off2(1), U1GATE)
}

void U1_Gate_MEM::run_mpi_mem(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,long long offset0,long long offset1,long long length){
    mpi_gate_mem(init_off2(long long, offset0, offset1),length, update_off2(1), U1GATEMPI)
}

RX_Gate_MEM::RX_Gate_MEM(vector<int> targ, double phi): ONE_QUBIT_GATE_MEM(targ) {
    name += "RX_Gate_MEM";
    cos_Phi_2 = complex<double>(cos(phi/2), 0);
    i_sin_Phi_2 = complex<double>(0, sin(-phi/2));
    if(mpi_count)
    {
        bind_gate_mpi_1_mem(RX_Gate_MEM)
    }
    else
    {
        bind_gate_1_mem(RX_Gate_MEM)
    }
}

void RX_Gate_MEM::run_chunk(vector<complex<double>> &buffer, long long idx){
    chunk_gate_mem(init_off2(long long, idx, idx + half_off), 
                update_off2(1), update_off2(half_off), RXGATE)
}

void RX_Gate_MEM::run_nonchunk(vector<complex<double>> &buffer, long long idx){
    nonchunk_gate_mem(init_off2(long long, idx, idx + half_off), 
                    update_off2(1), RXGATE)
}

void RX_Gate_MEM::run_mpi_mem(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,long long offset0,long long offset1,long long length){
    mpi_gate_mem(init_off2(long long, offset0, offset1),length, update_off2(1), RXGATEMPI)
}

RY_Gate_MEM::RY_Gate_MEM(vector<int> targ, double phi): ONE_QUBIT_GATE_MEM(targ) {
    name += "RY_Gate_MEM";
    cos_Phi_2 = complex<double>(cos(phi/2), 0);
    sin_Phi_2 = complex<double>(sin(phi/2), 0);
    if(mpi_count)
    {
        bind_gate_mpi_1_mem(RY_Gate_MEM)
    }
    else
    {
        bind_gate_1_mem(RY_Gate_MEM)
    }
}

void RY_Gate_MEM::run_chunk(vector<complex<double>> &buffer, long long idx){
    chunk_gate_mem(init_off2(long long, idx, idx + half_off), 
                update_off2(1), update_off2(half_off), RYGATE)
}

void RY_Gate_MEM::run_nonchunk(vector<complex<double>> &buffer, long long idx){
    nonchunk_gate_mem(init_off2(long long, idx, idx + half_off), 
                    update_off2(1), RYGATE)
}

void RY_Gate_MEM::run_mpi_mem(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,long long offset0,long long offset1,long long length){
    mpi_gate_mem(init_off2(long long, offset0, offset1),length, update_off2(1), RYGATEMPI)
}

RZ_Gate_MEM::RZ_Gate_MEM(vector<int> targ, double phi): ONE_QUBIT_GATE_MEM(targ) {
    name += "RZ_Gate_MEM";
    exp_p_iPhi_2 = complex<double>(cos(phi/2), sin(phi/2));
    exp_n_iPhi_2 = complex<double>(cos(-phi/2), sin(-phi/2));
    if(mpi_count)
    {
        bind_gate_mpi_1_mem_diagonal(RZ_Gate_MEM)
    }
    else
    {
        bind_gate_1_mem(RZ_Gate_MEM)
    }
}

void RZ_Gate_MEM::run_chunk(vector<complex<double>> &buffer, long long idx){
    chunk_gate_mem(init_off2(long long, idx, idx + half_off), 
                update_off2(1), update_off2(half_off), RZGATE)
}

void RZ_Gate_MEM::run_nonchunk(vector<complex<double>> &buffer, long long idx){
    nonchunk_gate_mem(init_off2(long long, idx, idx + half_off), 
                    update_off2(1), RZGATE)
}

void RZ_Gate_MEM::run_mpi_mem(vector<complex<double>> &buffer, long long idx,int pos,long long length){
    if(!pos)
    {
        mpi_gate_mem(init_off1(long long, idx),length, update_off1(1), RZ0)
    }
    else
    {
        mpi_gate_mem(init_off1(long long, idx),env.chunk_state, update_off1(1), RZ1)
    }
}

/*-----TWO_QUBIT_GATE-----*/
TWO_QUBIT_GATE_MEM::TWO_QUBIT_GATE_MEM(vector<int> targ): Gate(targ) {
    type = TWO_QUBIT;
    sort(targ.begin(), targ.end());
    half_off_0 = env.qubit_offset[targ[0]];
    half_off_1 = env.qubit_offset[targ[1]];
    off_0 = half_off_0 * 2;
    off_1 = half_off_1 * 2;
}

U2_Gate_MEM::U2_Gate_MEM(vector<int> targ, vector<complex<double>> matrix): TWO_QUBIT_GATE_MEM(targ) {
    name = "U2_Gate_MEM";
    assert(targ[0] < targ[1]);
    coeff.resize(16);
    for (int i = 0; i < 16; i++) {
        int res = 0;
        int b1 = i & 1;
        int b2 = (i >> 1) & 1;
        int b3 = (i >> 2) & 1;
        int b4 = (i >> 3) & 1;
        res = (b3 << 3) | (b4 << 2) | (b1 << 1) | b2;
        coeff[i] = matrix[res];
    }
    if(mpi_count)
    {
        if(chunk_count)
        {
            bind_gate_mpi_1_mem(U2_Gate_MEM)
        }
        else
        {
            bind_gate_mpi_mpi_u2_mem(U2_Gate_MEM)
        }
    }
    else
    {
        bind_gate_2_mem(U2_Gate_MEM)
    }
}

void U2_Gate_MEM::run_chunk_chunk(vector<complex<double>> &buffer, long long idx){
    chunk_chunk_gate_mem(init_off4(long long, idx,
                                    idx + half_off_0, 
                                    idx + half_off_1, 
                                    idx + half_off_0 + half_off_1),
                        update_off4(1),
                        update_off4(half_off_0),
                        update_off4(half_off_1),
                        U2Gate)
}

void U2_Gate_MEM::run_nonchunk_chunk(vector<complex<double>> &buffer, long long idx){
    nonchunk_chunk_gate_mem(init_off4(long long, idx,
                                    idx + half_off_0, 
                                    idx + half_off_1, 
                                    idx + half_off_0 + half_off_1),
                        update_off4(1),
                        update_off4(half_off_0),
                        U2Gate)
}

void U2_Gate_MEM::run_nonchunk_nonchunk(vector<complex<double>> &buffer, long long idx){
    nonchunk_nonchunk_gate_mem(init_off4(long long, idx,
                                    idx + half_off_0, 
                                    idx + half_off_1, 
                                    idx + half_off_0 + half_off_1),
                        update_off4(1),
                        U2Gate)
}

void U2_Gate_MEM::run_mpi_nonchunk(vector<complex<double>>& buffer1,vector<complex<double>>& buffer2,vector<complex<double>>& buffer3,vector<complex<double>>& buffer4,long long offset0,long long offset1,long long offset2,long long offset3,long long length){
    mpi_nonchunk_nonchunk_mem(init_off4(long long,offset0,offset1,offset2,offset3),length,update_off4(1),U2GENERAL)
}

void U2_Gate_MEM::run_mpi_mem(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,long long offset0,long long offset1,long long length){
    mpi_nonchunk_chunk_mem(init_off4(long long, offset0,offset0 + half_off_0,offset1,offset1 + half_off_0),length, update_off4(1),update_off4(half_off_0),U2MPICHUNK)
}

SWAP_Gate_MEM::SWAP_Gate_MEM(vector<int> targ): TWO_QUBIT_GATE_MEM(targ) {
    name = "SWAP_Gate_MEM";
    if(mpi_count == 2)
    {
        bind_gate_swap_mpi_mpi(SWAP_Gate_MEM)
    }
    else if(mpi_count)
    {
        if(chunk_count)
        {
            bind_gate_swap_mpi_chunk(SWAP_Gate_MEM)
        }
        else
        {
            bind_gate_swap_mpi_file_middle(SWAP_Gate_MEM)
        }
    }
    else
    {
        bind_gate_2_mem(SWAP_Gate_MEM);
    }
}

void SWAP_Gate_MEM::run_chunk_chunk(vector<complex<double>> &buffer, long long idx){
    chunk_chunk_gate_mem(init_off2(long long, 
                                    idx + half_off_0,
                                    idx + half_off_1),
                        update_off2(1),
                        update_off2(half_off_0),
                        update_off2(half_off_1),
                        SWAPGate)
}

void SWAP_Gate_MEM::run_nonchunk_chunk(vector<complex<double>> &buffer, long long idx){
    nonchunk_chunk_gate_mem(init_off2(long long, 
                                    idx + half_off_0,
                                    idx + half_off_1),
                        update_off2(1),
                        update_off2(half_off_0),
                        SWAPGate)
}

void SWAP_Gate_MEM::run_nonchunk_nonchunk(vector<complex<double>> &buffer, long long idx){
    nonchunk_nonchunk_gate_mem(init_off2(long long, 
                                    idx + half_off_0,
                                    idx + half_off_1),
                        update_off2(1),
                        SWAPGate)
}

void SWAP_Gate_MEM::run_mpi_mpi(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,long long offset0,long long offset1,long long length){
    mpi_gate_mem(init_off2(long long, offset0, offset1),length, update_off2(1), SWAPGateMPIChunk)
}

void SWAP_Gate_MEM::run_mpi_chunk(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,long long offset0,long long offset1,long long length){
    mpi_nonchunk_chunk_mem(init_off2(long long, offset0 + half_off_0, offset1),length, update_off2(1),update_off2(half_off_0),SWAPGateMPIChunk)
}

void SWAP_Gate_MEM::run_mpi_file_middle(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,long long offset0,long long offset1,long long length){
    mpi_nonchunk_nonchunk_mem(init_off2(long long, offset0, offset1),length, update_off2(1),SWAPGateMPIChunk)
}

VSWAP_Gate_1_1_MEM::VSWAP_Gate_1_1_MEM(vector<int> targ): SWAP_Gate_MEM(targ) {
    type = VSWAP;
    name == "VSWAP_Gate_1_1_MEM";
}

MPI_VSWAP_Gate_1_1_MEM::MPI_VSWAP_Gate_1_1_MEM(vector<int>targ):SWAP_Gate_MEM(targ){
    type = VSWAP;
    name == "MPI_VSWAP_Gate_1_1_MEM";
}

CPhase_Gate_MEM::CPhase_Gate_MEM(vector<int> targ, double phi): TWO_QUBIT_GATE_MEM(targ) {
    name = "CPhase_Gate_MEM";
    exp_iPhi = complex<double>(cos(phi), sin(phi));
    if(mpi_count)
    {
        bind_gate_mpi_1_mem_diagonal(CPhase_Gate_MEM)
    }
    else
    {
        bind_gate_2_mem(CPhase_Gate_MEM)
    }
}

void CPhase_Gate_MEM::run_chunk_chunk(vector<complex<double>> &buffer, long long idx){
    chunk_chunk_gate_mem(init_off1(long long, 
                                    idx + half_off_0 + half_off_1),
                        update_off1(1),
                        update_off1(half_off_0),
                        update_off1(half_off_1),
                        CPGATE)
}

void CPhase_Gate_MEM::run_nonchunk_chunk(vector<complex<double>> &buffer, long long idx){
    nonchunk_chunk_gate_mem(init_off1(long long, 
                                    idx + half_off_0 + half_off_1),
                        update_off1(1),
                        update_off1(half_off_0),
                        CPGATE)
}

void CPhase_Gate_MEM::run_nonchunk_nonchunk(vector<complex<double>> &buffer, long long idx){
    nonchunk_nonchunk_gate_mem(init_off1(long long, 
                                    idx + half_off_0 + half_off_1),
                        update_off1(1),
                        CPGATE)
}
void CPhase_Gate_MEM::run_mpi_mem(vector<complex<double>> &buffer, long long idx,int pos,long long length){
    if(pos < 3) return;
    mpi_gate_mem(init_off1(long long, idx),length, update_off1(1), CPGATE)
}

RZZ_Gate_MEM::RZZ_Gate_MEM(vector<int> targ, double phi): TWO_QUBIT_GATE_MEM(targ) {
    name = "RZZ_Gate_MEM";
    exp_p_iPhi_2 = complex<double>(cos(phi/2), sin(phi/2));
    exp_n_iPhi_2 = complex<double>(cos(-phi/2), sin(-phi/2));
    if(mpi_count)
    {
        bind_gate_mpi_1_mem_diagonal(RZZ_Gate_MEM)
    }
    else
    {
        bind_gate_2_mem(RZZ_Gate_MEM)
    }
}

void RZZ_Gate_MEM::run_chunk_chunk(vector<complex<double>> &buffer, long long idx){
    chunk_chunk_gate_mem(init_off4(long long, idx,
                                    idx + half_off_0, 
                                    idx + half_off_1, 
                                    idx + half_off_0 + half_off_1),
                        update_off4(1),
                        update_off4(half_off_0),
                        update_off4(half_off_1),
                        RZZGATE)
}

void RZZ_Gate_MEM::run_nonchunk_chunk(vector<complex<double>> &buffer, long long idx){
    nonchunk_chunk_gate_mem(init_off4(long long, idx,
                                    idx + half_off_0, 
                                    idx + half_off_1, 
                                    idx + half_off_0 + half_off_1),
                        update_off4(1),
                        update_off4(half_off_0),
                        RZZGATE)
}

void RZZ_Gate_MEM::run_nonchunk_nonchunk(vector<complex<double>> &buffer, long long idx){
    nonchunk_nonchunk_gate_mem(init_off4(long long, idx,
                                    idx + half_off_0, 
                                    idx + half_off_1, 
                                    idx + half_off_0 + half_off_1),
                        update_off4(1),
                        RZZGATE)
}

void RZZ_Gate_MEM::run_mpi_mem(vector<complex<double>> &buffer, long long idx,int pos,long long length){
    if(pos == 0 || pos == 3)
    {
        mpi_gate_mem(init_off1(long long, idx),length, update_off1(1), RZZN)
    }
    else
    {
        mpi_gate_mem(init_off1(long long, idx),length, update_off1(1), RZZP)
    }
}

/*-----THREE_QUBIT_GATE-----*/
U3_Gate_MEM::U3_Gate_MEM(vector<int> targ, vector<complex<double>> coeff): Gate(targ), coeff(coeff) {
    type = THREE_QUBIT;
    name += "U3_Gate";
    half_target_offset0 = env.qubit_offset[targ[0]];
    half_target_offset1 = env.qubit_offset[targ[1]];
    half_target_offset2 = env.qubit_offset[targ[2]];
    target_offset0 = half_target_offset0 * 2;
    target_offset1 = half_target_offset1 * 2;
    target_offset2 = half_target_offset2 * 2;

    if(isChunk(targ[2])){
        _run = bind(&U3_Gate_MEM::run_chunk_chunk_chunk, this, _1, _2);
        name += "_chunk_chunk_chunk";
    }
    if(isChunk(targ[1])){
        _run = bind(&U3_Gate_MEM::run_nonchunk_chunk_chunk, this, _1, _2);
        name += "_nonchunk_chunk_chunk";
    }
    else if (isChunk(targ[0])){
        _run = bind(&U3_Gate_MEM::run_nonchunk_nonchunk_chunk, this, _1, _2);
        name += "_nonchunk_nonchunk_chunk";
    }
    else {
        _run = bind(&U3_Gate_MEM::run_nonchunk_nonchunk_nonchunk, this, _1, _2);
        name += "_nonchunk_nonchunk_nonchunk";
    }
    name += "_MEM";
}

void U3_Gate_MEM::run_chunk_chunk_chunk(vector<complex<double>> &buffer, long long idx){
    long long off000 = idx;
    long long off001 = idx + half_target_offset0;
    long long off010 = idx + half_target_offset1;
    long long off011 = idx + half_target_offset1 + half_target_offset0;
    long long off100 = idx + half_target_offset2;
    long long off101 = idx + half_target_offset2 + half_target_offset0;
    long long off110 = idx + half_target_offset2 + half_target_offset1;
    long long off111 = idx + half_target_offset2 + half_target_offset1 + half_target_offset0;

    complex<double> q000;
    complex<double> q001;
    complex<double> q010;
    complex<double> q011;
    complex<double> q100;
    complex<double> q101;
    complex<double> q110;
    complex<double> q111;

    for (int i = 0; i < env.chunk_state; i += target_offset2) {
        for (int j = 0; j < half_target_offset2; j += target_offset1) {
            for (int k = 0; k < half_target_offset1; k += target_offset0) {
                for (int m = 0; m < half_target_offset0; m++) {
                    q000 = buffer[off000];
                    q001 = buffer[off001];
                    q010 = buffer[off010];
                    q011 = buffer[off011];
                    q100 = buffer[off100];
                    q101 = buffer[off101];
                    q110 = buffer[off110];
                    q111 = buffer[off111];
                    
                    buffer[off000] = coeff[ 0] * q000 + coeff[ 1] * q001 + coeff[ 2] * q010 + coeff[ 3] * q011 + coeff[ 4] * q100 + coeff[ 5] * q101 + coeff[ 6] * q110 + coeff[ 7] * q111;
                    buffer[off001] = coeff[ 8] * q000 + coeff[ 9] * q001 + coeff[10] * q010 + coeff[11] * q011 + coeff[12] * q100 + coeff[13] * q101 + coeff[14] * q110 + coeff[15] * q111;
                    buffer[off010] = coeff[16] * q000 + coeff[17] * q001 + coeff[18] * q010 + coeff[19] * q011 + coeff[20] * q100 + coeff[21] * q101 + coeff[22] * q110 + coeff[23] * q111;
                    buffer[off011] = coeff[24] * q000 + coeff[25] * q001 + coeff[26] * q010 + coeff[27] * q011 + coeff[28] * q100 + coeff[29] * q101 + coeff[30] * q110 + coeff[31] * q111;
                    buffer[off100] = coeff[32] * q000 + coeff[33] * q001 + coeff[34] * q010 + coeff[35] * q011 + coeff[36] * q100 + coeff[37] * q101 + coeff[38] * q110 + coeff[39] * q111;
                    buffer[off101] = coeff[40] * q000 + coeff[41] * q001 + coeff[42] * q010 + coeff[43] * q011 + coeff[44] * q100 + coeff[45] * q101 + coeff[46] * q110 + coeff[47] * q111;
                    buffer[off110] = coeff[48] * q000 + coeff[49] * q001 + coeff[50] * q010 + coeff[51] * q011 + coeff[52] * q100 + coeff[53] * q101 + coeff[54] * q110 + coeff[55] * q111;
                    buffer[off111] = coeff[56] * q000 + coeff[57] * q001 + coeff[58] * q010 + coeff[59] * q011 + coeff[60] * q100 + coeff[61] * q101 + coeff[62] * q110 + coeff[63] * q111;
                    
                    off000++;
                    off001++;
                    off010++;
                    off011++;
                    off100++;
                    off101++;
                    off110++;
                    off111++;
                }
                off000 += half_target_offset0;
                off001 += half_target_offset0;
                off010 += half_target_offset0;
                off011 += half_target_offset0;
                off100 += half_target_offset0;
                off101 += half_target_offset0;
                off110 += half_target_offset0;
                off111 += half_target_offset0;
            }
            off000 += half_target_offset1;
            off001 += half_target_offset1;
            off010 += half_target_offset1;
            off011 += half_target_offset1;
            off100 += half_target_offset1;
            off101 += half_target_offset1;
            off110 += half_target_offset1;
            off111 += half_target_offset1;
        }
        off000 += half_target_offset2;
        off001 += half_target_offset2;
        off010 += half_target_offset2;
        off011 += half_target_offset2;
        off100 += half_target_offset2;
        off101 += half_target_offset2;
        off110 += half_target_offset2;
        off111 += half_target_offset2;
    }
}

void U3_Gate_MEM::run_nonchunk_chunk_chunk(vector<complex<double>> &buffer, long long idx){
    long long off000 = idx;
    long long off001 = idx + half_target_offset0;
    long long off010 = idx + half_target_offset1;
    long long off011 = idx + half_target_offset1 + half_target_offset0;
    long long off100 = idx + half_target_offset2;
    long long off101 = idx + half_target_offset2 + half_target_offset0;
    long long off110 = idx + half_target_offset2 + half_target_offset1;
    long long off111 = idx + half_target_offset2 + half_target_offset1 + half_target_offset0;

    complex<double> q000;
    complex<double> q001;
    complex<double> q010;
    complex<double> q011;
    complex<double> q100;
    complex<double> q101;
    complex<double> q110;
    complex<double> q111;

    for (int j = 0; j < env.chunk_state; j += target_offset1) {
        for (int k = 0; k < half_target_offset1; k += target_offset0) {
            for (int m = 0; m < half_target_offset0; m++) {
                q000 = buffer[off000];
                q001 = buffer[off001];
                q010 = buffer[off010];
                q011 = buffer[off011];
                q100 = buffer[off100];
                q101 = buffer[off101];
                q110 = buffer[off110];
                q111 = buffer[off111];
                
                buffer[off000] = coeff[ 0] * q000 + coeff[ 1] * q001 + coeff[ 2] * q010 + coeff[ 3] * q011 + coeff[ 4] * q100 + coeff[ 5] * q101 + coeff[ 6] * q110 + coeff[ 7] * q111;
                buffer[off001] = coeff[ 8] * q000 + coeff[ 9] * q001 + coeff[10] * q010 + coeff[11] * q011 + coeff[12] * q100 + coeff[13] * q101 + coeff[14] * q110 + coeff[15] * q111;
                buffer[off010] = coeff[16] * q000 + coeff[17] * q001 + coeff[18] * q010 + coeff[19] * q011 + coeff[20] * q100 + coeff[21] * q101 + coeff[22] * q110 + coeff[23] * q111;
                buffer[off011] = coeff[24] * q000 + coeff[25] * q001 + coeff[26] * q010 + coeff[27] * q011 + coeff[28] * q100 + coeff[29] * q101 + coeff[30] * q110 + coeff[31] * q111;
                buffer[off100] = coeff[32] * q000 + coeff[33] * q001 + coeff[34] * q010 + coeff[35] * q011 + coeff[36] * q100 + coeff[37] * q101 + coeff[38] * q110 + coeff[39] * q111;
                buffer[off101] = coeff[40] * q000 + coeff[41] * q001 + coeff[42] * q010 + coeff[43] * q011 + coeff[44] * q100 + coeff[45] * q101 + coeff[46] * q110 + coeff[47] * q111;
                buffer[off110] = coeff[48] * q000 + coeff[49] * q001 + coeff[50] * q010 + coeff[51] * q011 + coeff[52] * q100 + coeff[53] * q101 + coeff[54] * q110 + coeff[55] * q111;
                buffer[off111] = coeff[56] * q000 + coeff[57] * q001 + coeff[58] * q010 + coeff[59] * q011 + coeff[60] * q100 + coeff[61] * q101 + coeff[62] * q110 + coeff[63] * q111;
                
                off000++;
                off001++;
                off010++;
                off011++;
                off100++;
                off101++;
                off110++;
                off111++;
            }
            off000 += half_target_offset0;
            off001 += half_target_offset0;
            off010 += half_target_offset0;
            off011 += half_target_offset0;
            off100 += half_target_offset0;
            off101 += half_target_offset0;
            off110 += half_target_offset0;
            off111 += half_target_offset0;
        }
        off000 += half_target_offset1;
        off001 += half_target_offset1;
        off010 += half_target_offset1;
        off011 += half_target_offset1;
        off100 += half_target_offset1;
        off101 += half_target_offset1;
        off110 += half_target_offset1;
        off111 += half_target_offset1;
    }
}

void U3_Gate_MEM::run_nonchunk_nonchunk_chunk(vector<complex<double>> &buffer, long long idx){
    long long off000 = idx;
    long long off001 = idx + half_target_offset0;
    long long off010 = idx + half_target_offset1;
    long long off011 = idx + half_target_offset1 + half_target_offset0;
    long long off100 = idx + half_target_offset2;
    long long off101 = idx + half_target_offset2 + half_target_offset0;
    long long off110 = idx + half_target_offset2 + half_target_offset1;
    long long off111 = idx + half_target_offset2 + half_target_offset1 + half_target_offset0;

    complex<double> q000;
    complex<double> q001;
    complex<double> q010;
    complex<double> q011;
    complex<double> q100;
    complex<double> q101;
    complex<double> q110;
    complex<double> q111;

    for (int k = 0; k < env.chunk_state; k += target_offset0) {
        for (int m = 0; m < half_target_offset0; m++) {
            q000 = buffer[off000];
            q001 = buffer[off001];
            q010 = buffer[off010];
            q011 = buffer[off011];
            q100 = buffer[off100];
            q101 = buffer[off101];
            q110 = buffer[off110];
            q111 = buffer[off111];
            
            buffer[off000] = coeff[ 0] * q000 + coeff[ 1] * q001 + coeff[ 2] * q010 + coeff[ 3] * q011 + coeff[ 4] * q100 + coeff[ 5] * q101 + coeff[ 6] * q110 + coeff[ 7] * q111;
            buffer[off001] = coeff[ 8] * q000 + coeff[ 9] * q001 + coeff[10] * q010 + coeff[11] * q011 + coeff[12] * q100 + coeff[13] * q101 + coeff[14] * q110 + coeff[15] * q111;
            buffer[off010] = coeff[16] * q000 + coeff[17] * q001 + coeff[18] * q010 + coeff[19] * q011 + coeff[20] * q100 + coeff[21] * q101 + coeff[22] * q110 + coeff[23] * q111;
            buffer[off011] = coeff[24] * q000 + coeff[25] * q001 + coeff[26] * q010 + coeff[27] * q011 + coeff[28] * q100 + coeff[29] * q101 + coeff[30] * q110 + coeff[31] * q111;
            buffer[off100] = coeff[32] * q000 + coeff[33] * q001 + coeff[34] * q010 + coeff[35] * q011 + coeff[36] * q100 + coeff[37] * q101 + coeff[38] * q110 + coeff[39] * q111;
            buffer[off101] = coeff[40] * q000 + coeff[41] * q001 + coeff[42] * q010 + coeff[43] * q011 + coeff[44] * q100 + coeff[45] * q101 + coeff[46] * q110 + coeff[47] * q111;
            buffer[off110] = coeff[48] * q000 + coeff[49] * q001 + coeff[50] * q010 + coeff[51] * q011 + coeff[52] * q100 + coeff[53] * q101 + coeff[54] * q110 + coeff[55] * q111;
            buffer[off111] = coeff[56] * q000 + coeff[57] * q001 + coeff[58] * q010 + coeff[59] * q011 + coeff[60] * q100 + coeff[61] * q101 + coeff[62] * q110 + coeff[63] * q111;
            
            off000++;
            off001++;
            off010++;
            off011++;
            off100++;
            off101++;
            off110++;
            off111++;
        }
        off000 += half_target_offset0;
        off001 += half_target_offset0;
        off010 += half_target_offset0;
        off011 += half_target_offset0;
        off100 += half_target_offset0;
        off101 += half_target_offset0;
        off110 += half_target_offset0;
        off111 += half_target_offset0;
    }
}

void U3_Gate_MEM::run_nonchunk_nonchunk_nonchunk(vector<complex<double>> &buffer, long long idx){
    long long off000 = idx;
    long long off001 = idx + half_target_offset0;
    long long off010 = idx + half_target_offset1;
    long long off011 = idx + half_target_offset1 + half_target_offset0;
    long long off100 = idx + half_target_offset2;
    long long off101 = idx + half_target_offset2 + half_target_offset0;
    long long off110 = idx + half_target_offset2 + half_target_offset1;
    long long off111 = idx + half_target_offset2 + half_target_offset1 + half_target_offset0;

    complex<double> q000;
    complex<double> q001;
    complex<double> q010;
    complex<double> q011;
    complex<double> q100;
    complex<double> q101;
    complex<double> q110;
    complex<double> q111;

    for (int m = 0; m < env.chunk_state; m++) {
        q000 = buffer[off000];
        q001 = buffer[off001];
        q010 = buffer[off010];
        q011 = buffer[off011];
        q100 = buffer[off100];
        q101 = buffer[off101];
        q110 = buffer[off110];
        q111 = buffer[off111];
        
        buffer[off000] = coeff[ 0] * q000 + coeff[ 1] * q001 + coeff[ 2] * q010 + coeff[ 3] * q011 + coeff[ 4] * q100 + coeff[ 5] * q101 + coeff[ 6] * q110 + coeff[ 7] * q111;
        buffer[off001] = coeff[ 8] * q000 + coeff[ 9] * q001 + coeff[10] * q010 + coeff[11] * q011 + coeff[12] * q100 + coeff[13] * q101 + coeff[14] * q110 + coeff[15] * q111;
        buffer[off010] = coeff[16] * q000 + coeff[17] * q001 + coeff[18] * q010 + coeff[19] * q011 + coeff[20] * q100 + coeff[21] * q101 + coeff[22] * q110 + coeff[23] * q111;
        buffer[off011] = coeff[24] * q000 + coeff[25] * q001 + coeff[26] * q010 + coeff[27] * q011 + coeff[28] * q100 + coeff[29] * q101 + coeff[30] * q110 + coeff[31] * q111;
        buffer[off100] = coeff[32] * q000 + coeff[33] * q001 + coeff[34] * q010 + coeff[35] * q011 + coeff[36] * q100 + coeff[37] * q101 + coeff[38] * q110 + coeff[39] * q111;
        buffer[off101] = coeff[40] * q000 + coeff[41] * q001 + coeff[42] * q010 + coeff[43] * q011 + coeff[44] * q100 + coeff[45] * q101 + coeff[46] * q110 + coeff[47] * q111;
        buffer[off110] = coeff[48] * q000 + coeff[49] * q001 + coeff[50] * q010 + coeff[51] * q011 + coeff[52] * q100 + coeff[53] * q101 + coeff[54] * q110 + coeff[55] * q111;
        buffer[off111] = coeff[56] * q000 + coeff[57] * q001 + coeff[58] * q010 + coeff[59] * q011 + coeff[60] * q100 + coeff[61] * q101 + coeff[62] * q110 + coeff[63] * q111;
        
        off000++;
        off001++;
        off010++;
        off011++;
        off100++;
        off101++;
        off110++;
        off111++;
    }
}


VSWAP_Gate_2_2_MEM::VSWAP_Gate_2_2_MEM(vector<int> targ): Gate(targ) {
    type = VSWAP;
    swap_out = {targ[0], targ[1]};
    swap_in  = {targ[2], targ[3]};

    half_target_offset0 = env.qubit_offset[targ[0]];
    half_target_offset1 = env.qubit_offset[targ[1]];
    half_target_offset2 = env.qubit_offset[targ[2]];
    half_target_offset3 = env.qubit_offset[targ[3]];
        
    target_offset0 = half_target_offset0 * 2;
    target_offset1 = half_target_offset1 * 2;
    target_offset2 = half_target_offset2 * 2;
    target_offset3 = half_target_offset3 * 2;
        
    _run = bind(&VSWAP_Gate_2_2_MEM::run_nonchunk_chunk, this, _1, _2);
    name += "VSWAP_Gate_2_2_MEM";
}

MPI_VSWAP_Gate_2_2_MEM::MPI_VSWAP_Gate_2_2_MEM(vector<int>targ): Gate(targ){
    type = VSWAP;
    name = "MPI_VSWAP_Gate_2_2_MEM";
    half_off_0 = env.qubit_offset[targ[0]];
    half_off_1 = env.qubit_offset[targ[1]];
    half_off_2 = 0;
    half_off_3 = 0;
    off_0 = half_off_0 * 2;
    off_1 = half_off_1 * 2;
    off_2 = half_off_2 * 2;
    off_3 = half_off_3 * 2;
    run_mpi_vswap2_2_mem = bind(&MPI_VSWAP_Gate_2_2_MEM::run_vswap_2_2,this,_1,_2,_3,_4,_5,_6,_7,_8,_9);
}

void MPI_VSWAP_Gate_2_2_MEM::run_vswap_2_2(vector<complex<double>>&buffer1,vector<complex<double>>&buffer2,vector<complex<double>>&buffer3,vector<complex<double>>&buffer4,long long offset0,long long offset1,long long offset2,long long offset3,long long length)
{
    int off0001 = offset0 +half_off_0;
    int off0010 = offset0 +half_off_1;
    int off0011 = offset0 +half_off_0 + half_off_1;
    int off0100 = offset1 + half_off_2;
    int off0110 = offset1 + half_off_1 + half_off_2;
    int off0111 = offset1 + half_off_0 + half_off_1 + half_off_2;
    int off1000 = offset2 + half_off_3;
    int off1001 = offset2 + half_off_0 + half_off_3;
    int off1011 = offset2 + half_off_0 + half_off_1 + half_off_3;
    int off1100 = offset3 + half_off_2 + half_off_3;
    int off1101 = offset3 + half_off_0 + half_off_2 + half_off_3;
    int off1110 = offset3 + half_off_1 + half_off_2 + half_off_3;
    complex<double> q0001;
    complex<double> q0010;
    complex<double> q0011;
    complex<double> q0100;
    complex<double> q0110;
    complex<double> q0111;
    complex<double> q1000;
    complex<double> q1001;
    complex<double> q1011;
    complex<double> q1100;
    complex<double> q1101;
    complex<double> q1110;
    for (int i = 0; i < length; i += off_1) {
        for (int j = 0; j < half_off_1; j += off_0) {
            for (int k = 0; k < half_off_0; k++) {
                q0001 = buffer1[off0001];
                q0010 = buffer1[off0010];
                q0011 = buffer1[off0011];
                q0100 = buffer2[off0100];
                q0110 = buffer2[off0110];
                q0111 = buffer2[off0111];
                q1000 = buffer3[off1000];
                q1001 = buffer3[off1001];
                q1011 = buffer3[off1011];
                q1100 = buffer4[off1100];
                q1101 = buffer4[off1101];
                q1110 = buffer4[off1110];
                buffer1[off0001] = q0100;
                buffer1[off0010] = q1000;
                buffer1[off0011] = q1100;
                buffer2[off0100] = q0001;
                buffer2[off0110] = q1001;
                buffer2[off0111] = q1101;
                buffer3[off1000] = q0010;
                buffer3[off1001] = q0110;
                buffer3[off1011] = q1110;
                buffer4[off1100] = q0011;
                buffer4[off1101] = q0111;
                buffer4[off1110] = q1011;
                off0001++;
                off0010++;
                off0011++;
                off0100++;
                off0110++;
                off0111++;
                off1000++;
                off1001++;
                off1011++;
                off1100++;
                off1101++;
                off1110++;
            }
            off0001 += half_off_0;
            off0010 += half_off_0;
            off0011 += half_off_0;
            off0100 += half_off_0;
            off0110 += half_off_0;
            off0111 += half_off_0;
            off1000 += half_off_0;
            off1001 += half_off_0;
            off1011 += half_off_0;
            off1100 += half_off_0;
            off1101 += half_off_0;
            off1110 += half_off_0;
        }
        off0001 += half_off_1;
        off0010 += half_off_1;
        off0011 += half_off_1;
        off0100 += half_off_1;
        off0110 += half_off_1;
        off0111 += half_off_1;
        off1000 += half_off_1;
        off1001 += half_off_1;
        off1011 += half_off_1;
        off1100 += half_off_1;
        off1101 += half_off_1;
        off1110 += half_off_1;
    }
}

void VSWAP_Gate_2_2_MEM::run_nonchunk_chunk(vector<complex<double>> &buffer, long long idx){
    long long off0001 = idx + half_target_offset0;
    long long off0010 = idx + half_target_offset1;
    long long off0011 = idx + half_target_offset0 + half_target_offset1;
    long long off0100 = idx + half_target_offset2;
    long long off0110 = idx + half_target_offset1 + half_target_offset2;
    long long off0111 = idx + half_target_offset0 + half_target_offset1 + half_target_offset2;
    long long off1000 = idx + half_target_offset3;
    long long off1001 = idx + half_target_offset0 + half_target_offset3;
    long long off1011 = idx + half_target_offset0 + half_target_offset1 + half_target_offset3;
    long long off1100 = idx + half_target_offset2 + half_target_offset3;
    long long off1101 = idx + half_target_offset0 + half_target_offset2 + half_target_offset3;
    long long off1110 = idx + half_target_offset1 + half_target_offset2 + half_target_offset3;

    complex<double> q0001;
    complex<double> q0010;
    complex<double> q0011;
    complex<double> q0100;
    complex<double> q0110;
    complex<double> q0111;
    complex<double> q1000;
    complex<double> q1001;
    complex<double> q1011;
    complex<double> q1100;
    complex<double> q1101;
    complex<double> q1110;

    for (int i = 0; i < env.chunk_state; i += target_offset1) {
        for (int j = 0; j < half_target_offset1; j += target_offset0) {
            for (int k = 0; k < half_target_offset0; k++) {
                q0001 = buffer[off0001];
                q0100 = buffer[off0100];
                buffer[off0001] = q0100;
                buffer[off0100] = q0001;

                q0010 = buffer[off0010];
                q1000 = buffer[off1000];
                buffer[off0010] = q1000;
                buffer[off1000] = q0010;

                q0011 = buffer[off0011];
                q1100 = buffer[off1100];
                buffer[off0011] = q1100;
                buffer[off1100] = q0011;

                q0110 = buffer[off0110];
                q1001 = buffer[off1001];
                buffer[off0110] = q1001;
                buffer[off1001] = q0110;

                q0111 = buffer[off0111];
                q1101 = buffer[off1101];
                buffer[off0111] = q1101;
                buffer[off1101] = q0111;

                q1011 = buffer[off1011];
                q1110 = buffer[off1110];
                buffer[off1011] = q1110;
                buffer[off1110] = q1011;
                
                off0001++;
                off0010++;
                off0011++;
                off0100++;
                off0110++;
                off0111++;
                off1000++;
                off1001++;
                off1011++;
                off1100++;
                off1101++;
                off1110++;
            }
            off0001 += half_target_offset0;
            off0010 += half_target_offset0;
            off0011 += half_target_offset0;
            off0100 += half_target_offset0;
            off0110 += half_target_offset0;
            off0111 += half_target_offset0;
            off1000 += half_target_offset0;
            off1001 += half_target_offset0;
            off1011 += half_target_offset0;
            off1100 += half_target_offset0;
            off1101 += half_target_offset0;
            off1110 += half_target_offset0;
        }
        off0001 += half_target_offset1;
        off0010 += half_target_offset1;
        off0011 += half_target_offset1;
        off0100 += half_target_offset1;
        off0110 += half_target_offset1;
        off0111 += half_target_offset1;
        off1000 += half_target_offset1;
        off1001 += half_target_offset1;
        off1011 += half_target_offset1;
        off1100 += half_target_offset1;
        off1101 += half_target_offset1;
        off1110 += half_target_offset1;
    }
}

// VSWAP_Gate_3_3_MEM::VSWAP_Gate_3_3_MEM(vector<int> targ): Gate(targ) {
//     swap_out = {targ[0], targ[1], targ[2]};
//     swap_in  = {targ[3], targ[4], targ[5]};

//     half_targ_off0 = env.qubit_offset[targ[0]];
//     half_targ_off1 = env.qubit_offset[targ[1]];
//     half_targ_off2 = env.qubit_offset[targ[2]];
//     half_targ_off3 = env.chunk_state;
//     half_targ_off4 = 2 * env.chunk_state;
//     half_targ_off5 = 4 * env.chunk_state;
        
//     targ_off0 = half_targ_off0 * 2;
//     targ_off1 = half_targ_off1 * 2;
//     targ_off2 = half_targ_off2 * 2;
//     targ_off3 = half_targ_off3 * 2;
//     targ_off4 = half_targ_off4 * 2;
//     targ_off5 = half_targ_off5 * 2;
        
//     _run = bind(&VSWAP_Gate_3_3_MEM::run_nonchunk_chunk, this, _1, _2);
//     name += "VSWAP_Gate_3_3_MEM";
// }

// void VSWAP_Gate_3_3_MEM::run_nonchunk_chunk(vector<complex<double>> &buffer, long long idx){
//     long long off000001 = idx + half_targ_off0;
//     long long off000010 = idx + half_targ_off1;
//     long long off000011 = idx + half_targ_off0 + half_targ_off1;
//     long long off000100 = idx + half_targ_off2;
//     long long off000101 = idx + half_targ_off0 + half_targ_off2;
//     long long off000110 = idx + half_targ_off1 + half_targ_off2;
//     long long off000111 = idx + half_targ_off0 + half_targ_off1 + half_targ_off2;
//     long long off001000 = idx + half_targ_off3;
//     long long off001010 = idx + half_targ_off1 + half_targ_off3;
//     long long off001011 = idx + half_targ_off0 + half_targ_off1 + half_targ_off3;
//     long long off001100 = idx + half_targ_off2 + half_targ_off3;
//     long long off001101 = idx + half_targ_off0 + half_targ_off2 + half_targ_off3;
//     long long off001110 = idx + half_targ_off1 + half_targ_off2 + half_targ_off3;
//     long long off001111 = idx + half_targ_off0 + half_targ_off1 + half_targ_off2 + half_targ_off3;
//     long long off010000 = idx + half_targ_off4;
//     long long off010001 = idx + half_targ_off0 + half_targ_off4;
//     long long off010011 = idx + half_targ_off0 + half_targ_off1 + half_targ_off4;
//     long long off010100 = idx + half_targ_off2 + half_targ_off4;
//     long long off010101 = idx + half_targ_off0 + half_targ_off2 + half_targ_off4;
//     long long off010110 = idx + half_targ_off1 + half_targ_off2 + half_targ_off4;
//     long long off010111 = idx + half_targ_off0 + half_targ_off1 + half_targ_off2 + half_targ_off4;
//     long long off011000 = idx + half_targ_off3 + half_targ_off4;
//     long long off011001 = idx + half_targ_off0 + half_targ_off3 + half_targ_off4;
//     long long off011010 = idx + half_targ_off1 + half_targ_off3 + half_targ_off4;
//     long long off011100 = idx + half_targ_off2 + half_targ_off3 + half_targ_off4;
//     long long off011101 = idx + half_targ_off0 + half_targ_off2 + half_targ_off3 + half_targ_off4;
//     long long off011110 = idx + half_targ_off1 + half_targ_off2 + half_targ_off3 + half_targ_off4;
//     long long off011111 = idx + half_targ_off0 + half_targ_off1 + half_targ_off2 + half_targ_off3 + half_targ_off4;
//     long long off100000 = idx + half_targ_off5;
//     long long off100001 = idx + half_targ_off0 + half_targ_off5;
//     long long off100010 = idx + half_targ_off1 + half_targ_off5;
//     long long off100011 = idx + half_targ_off0 + half_targ_off1 + half_targ_off5;
//     long long off100101 = idx + half_targ_off0 + half_targ_off2 + half_targ_off5;
//     long long off100110 = idx + half_targ_off1 + half_targ_off2 + half_targ_off5;
//     long long off100111 = idx + half_targ_off0 + half_targ_off1 + half_targ_off2 + half_targ_off5;
//     long long off101000 = idx + half_targ_off3 + half_targ_off5;
//     long long off101001 = idx + half_targ_off0 + half_targ_off3 + half_targ_off5;
//     long long off101010 = idx + half_targ_off1 + half_targ_off3 + half_targ_off5;
//     long long off101011 = idx + half_targ_off0 + half_targ_off1 + half_targ_off3 + half_targ_off5;
//     long long off101100 = idx + half_targ_off2 + half_targ_off3 + half_targ_off5;
//     long long off101110 = idx + half_targ_off1 + half_targ_off2 + half_targ_off3 + half_targ_off5;
//     long long off101111 = idx + half_targ_off0 + half_targ_off1 + half_targ_off2 + half_targ_off3 + half_targ_off5;
//     long long off110000 = idx + half_targ_off4 + half_targ_off5;
//     long long off110001 = idx + half_targ_off0 + half_targ_off4 + half_targ_off5;
//     long long off110010 = idx + half_targ_off1 + half_targ_off4 + half_targ_off5;
//     long long off110011 = idx + half_targ_off0 + half_targ_off1 + half_targ_off4 + half_targ_off5;
//     long long off110100 = idx + half_targ_off2 + half_targ_off4 + half_targ_off5;
//     long long off110101 = idx + half_targ_off0 + half_targ_off2 + half_targ_off4 + half_targ_off5;
//     long long off110111 = idx + half_targ_off0 + half_targ_off1 + half_targ_off2 + half_targ_off4 + half_targ_off5;
//     long long off111000 = idx + half_targ_off3 + half_targ_off4 + half_targ_off5;
//     long long off111001 = idx + half_targ_off0 + half_targ_off3 + half_targ_off4 + half_targ_off5;
//     long long off111010 = idx + half_targ_off1 + half_targ_off3 + half_targ_off4 + half_targ_off5;
//     long long off111011 = idx + half_targ_off0 + half_targ_off1 + half_targ_off3 + half_targ_off4 + half_targ_off5;
//     long long off111100 = idx + half_targ_off2 + half_targ_off3 + half_targ_off4 + half_targ_off5;
//     long long off111101 = idx + half_targ_off0 + half_targ_off2 + half_targ_off3 + half_targ_off4 + half_targ_off5;
//     long long off111110 = idx + half_targ_off1 + half_targ_off2 + half_targ_off3 + half_targ_off4 + half_targ_off5;

//     complex<double> q000001;
//     complex<double> q000010;
//     complex<double> q000011;
//     complex<double> q000100;
//     complex<double> q000101;
//     complex<double> q000110;
//     complex<double> q000111;
//     complex<double> q001000;
//     complex<double> q001010;
//     complex<double> q001011;
//     complex<double> q001100;
//     complex<double> q001101;
//     complex<double> q001110;
//     complex<double> q001111;
//     complex<double> q010000;
//     complex<double> q010001;
//     complex<double> q010011;
//     complex<double> q010100;
//     complex<double> q010101;
//     complex<double> q010110;
//     complex<double> q010111;
//     complex<double> q011000;
//     complex<double> q011001;
//     complex<double> q011010;
//     complex<double> q011100;
//     complex<double> q011101;
//     complex<double> q011110;
//     complex<double> q011111;
//     complex<double> q100000;
//     complex<double> q100001;
//     complex<double> q100010;
//     complex<double> q100011;
//     complex<double> q100101;
//     complex<double> q100110;
//     complex<double> q100111;
//     complex<double> q101000;
//     complex<double> q101001;
//     complex<double> q101010;
//     complex<double> q101011;
//     complex<double> q101100;
//     complex<double> q101110;
//     complex<double> q101111;
//     complex<double> q110000;
//     complex<double> q110001;
//     complex<double> q110010;
//     complex<double> q110011;
//     complex<double> q110100;
//     complex<double> q110101;
//     complex<double> q110111;
//     complex<double> q111000;
//     complex<double> q111001;
//     complex<double> q111010;
//     complex<double> q111011;
//     complex<double> q111100;
//     complex<double> q111101;
//     complex<double> q111110;

//     for (int i = 0; i < env.chunk_state; i += targ_off2) {
//         for (int j = 0; j < half_targ_off2; j += targ_off1) {
//             for (int k = 0; k < half_targ_off1; k += half_targ_off0) {
//                 for (int m = 0; m < half_targ_off0; m++) {
//                     q000001 = buffer[off000001];
//                     q000010 = buffer[off000010];
//                     q000011 = buffer[off000011];
//                     q000100 = buffer[off000100];
//                     q000101 = buffer[off000101];
//                     q000110 = buffer[off000110];
//                     q000111 = buffer[off000111];
//                     q001000 = buffer[off001000];
//                     q001010 = buffer[off001010];
//                     q001011 = buffer[off001011];
//                     q001100 = buffer[off001100];
//                     q001101 = buffer[off001101];
//                     q001110 = buffer[off001110];
//                     q001111 = buffer[off001111];
//                     q010000 = buffer[off010000];
//                     q010001 = buffer[off010001];
//                     q010011 = buffer[off010011];
//                     q010100 = buffer[off010100];
//                     q010101 = buffer[off010101];
//                     q010110 = buffer[off010110];
//                     q010111 = buffer[off010111];
//                     q011000 = buffer[off011000];
//                     q011001 = buffer[off011001];
//                     q011010 = buffer[off011010];
//                     q011100 = buffer[off011100];
//                     q011101 = buffer[off011101];
//                     q011110 = buffer[off011110];
//                     q011111 = buffer[off011111];
//                     q100000 = buffer[off100000];
//                     q100001 = buffer[off100001];
//                     q100010 = buffer[off100010];
//                     q100011 = buffer[off100011];
//                     q100101 = buffer[off100101];
//                     q100110 = buffer[off100110];
//                     q100111 = buffer[off100111];
//                     q101000 = buffer[off101000];
//                     q101001 = buffer[off101001];
//                     q101010 = buffer[off101010];
//                     q101011 = buffer[off101011];
//                     q101100 = buffer[off101100];
//                     q101110 = buffer[off101110];
//                     q101111 = buffer[off101111];
//                     q110000 = buffer[off110000];
//                     q110001 = buffer[off110001];
//                     q110010 = buffer[off110010];
//                     q110011 = buffer[off110011];
//                     q110100 = buffer[off110100];
//                     q110101 = buffer[off110101];
//                     q110111 = buffer[off110111];
//                     q111000 = buffer[off111000];
//                     q111001 = buffer[off111001];
//                     q111010 = buffer[off111010];
//                     q111011 = buffer[off111011];
//                     q111100 = buffer[off111100];
//                     q111101 = buffer[off111101];
//                     q111110 = buffer[off111110];

//                     buffer[off000001] = q001000;
//                     buffer[off000010] = q010000;
//                     buffer[off000011] = q011000;
//                     buffer[off000100] = q100000;
//                     buffer[off000101] = q101000;
//                     buffer[off000110] = q110000;
//                     buffer[off000111] = q111000;
//                     buffer[off001000] = q000001;
//                     buffer[off001010] = q010001;
//                     buffer[off001011] = q011001;
//                     buffer[off001100] = q100001;
//                     buffer[off001101] = q101001;
//                     buffer[off001110] = q110001;
//                     buffer[off001111] = q111001;
//                     buffer[off010000] = q000010;
//                     buffer[off010001] = q001010;
//                     buffer[off010011] = q011010;
//                     buffer[off010100] = q100010;
//                     buffer[off010101] = q101010;
//                     buffer[off010110] = q110010;
//                     buffer[off010111] = q111010;
//                     buffer[off011000] = q000011;
//                     buffer[off011001] = q001011;
//                     buffer[off011010] = q010011;
//                     buffer[off011100] = q100011;
//                     buffer[off011101] = q101011;
//                     buffer[off011110] = q110011;
//                     buffer[off011111] = q111011;
//                     buffer[off100000] = q000100;
//                     buffer[off100001] = q001100;
//                     buffer[off100010] = q010100;
//                     buffer[off100011] = q011100;
//                     buffer[off100101] = q101100;
//                     buffer[off100110] = q110100;
//                     buffer[off100111] = q111100;
//                     buffer[off101000] = q000101;
//                     buffer[off101001] = q001101;
//                     buffer[off101010] = q010101;
//                     buffer[off101011] = q011101;
//                     buffer[off101100] = q100101;
//                     buffer[off101110] = q110101;
//                     buffer[off101111] = q111101;
//                     buffer[off110000] = q000110;
//                     buffer[off110001] = q001110;
//                     buffer[off110010] = q010110;
//                     buffer[off110011] = q011110;
//                     buffer[off110100] = q100110;
//                     buffer[off110101] = q101110;
//                     buffer[off110111] = q111110;
//                     buffer[off111000] = q000111;
//                     buffer[off111001] = q001111;
//                     buffer[off111010] = q010111;
//                     buffer[off111011] = q011111;
//                     buffer[off111100] = q100111;
//                     buffer[off111101] = q101111;
//                     buffer[off111110] = q110111;
                    
//                     off000001++;
//                     off000010++;
//                     off000011++;
//                     off000100++;
//                     off000101++;
//                     off000110++;
//                     off000111++;
//                     off001000++;
//                     off001010++;
//                     off001011++;
//                     off001100++;
//                     off001101++;
//                     off001110++;
//                     off001111++;
//                     off010000++;
//                     off010001++;
//                     off010011++;
//                     off010100++;
//                     off010101++;
//                     off010110++;
//                     off010111++;
//                     off011000++;
//                     off011001++;
//                     off011010++;
//                     off011100++;
//                     off011101++;
//                     off011110++;
//                     off011111++;
//                     off100000++;
//                     off100001++;
//                     off100010++;
//                     off100011++;
//                     off100101++;
//                     off100110++;
//                     off100111++;
//                     off101000++;
//                     off101001++;
//                     off101010++;
//                     off101011++;
//                     off101100++;
//                     off101110++;
//                     off101111++;
//                     off110000++;
//                     off110001++;
//                     off110010++;
//                     off110011++;
//                     off110100++;
//                     off110101++;
//                     off110111++;
//                     off111000++;
//                     off111001++;
//                     off111010++;
//                     off111011++;
//                     off111100++;
//                     off111101++;
//                     off111110++;
//                 }
//                 off000001 += half_targ_off0;
//                 off000010 += half_targ_off0;
//                 off000011 += half_targ_off0;
//                 off000100 += half_targ_off0;
//                 off000101 += half_targ_off0;
//                 off000110 += half_targ_off0;
//                 off000111 += half_targ_off0;
//                 off001000 += half_targ_off0;
//                 off001010 += half_targ_off0;
//                 off001011 += half_targ_off0;
//                 off001100 += half_targ_off0;
//                 off001101 += half_targ_off0;
//                 off001110 += half_targ_off0;
//                 off001111 += half_targ_off0;
//                 off010000 += half_targ_off0;
//                 off010001 += half_targ_off0;
//                 off010011 += half_targ_off0;
//                 off010100 += half_targ_off0;
//                 off010101 += half_targ_off0;
//                 off010110 += half_targ_off0;
//                 off010111 += half_targ_off0;
//                 off011000 += half_targ_off0;
//                 off011001 += half_targ_off0;
//                 off011010 += half_targ_off0;
//                 off011100 += half_targ_off0;
//                 off011101 += half_targ_off0;
//                 off011110 += half_targ_off0;
//                 off011111 += half_targ_off0;
//                 off100000 += half_targ_off0;
//                 off100001 += half_targ_off0;
//                 off100010 += half_targ_off0;
//                 off100011 += half_targ_off0;
//                 off100101 += half_targ_off0;
//                 off100110 += half_targ_off0;
//                 off100111 += half_targ_off0;
//                 off101000 += half_targ_off0;
//                 off101001 += half_targ_off0;
//                 off101010 += half_targ_off0;
//                 off101011 += half_targ_off0;
//                 off101100 += half_targ_off0;
//                 off101110 += half_targ_off0;
//                 off101111 += half_targ_off0;
//                 off110000 += half_targ_off0;
//                 off110001 += half_targ_off0;
//                 off110010 += half_targ_off0;
//                 off110011 += half_targ_off0;
//                 off110100 += half_targ_off0;
//                 off110101 += half_targ_off0;
//                 off110111 += half_targ_off0;
//                 off111000 += half_targ_off0;
//                 off111001 += half_targ_off0;
//                 off111010 += half_targ_off0;
//                 off111011 += half_targ_off0;
//                 off111100 += half_targ_off0;
//                 off111101 += half_targ_off0;
//                 off111110 += half_targ_off0;
//             }
//             off000001 += half_targ_off1;
//             off000010 += half_targ_off1;
//             off000011 += half_targ_off1;
//             off000100 += half_targ_off1;
//             off000101 += half_targ_off1;
//             off000110 += half_targ_off1;
//             off000111 += half_targ_off1;
//             off001000 += half_targ_off1;
//             off001010 += half_targ_off1;
//             off001011 += half_targ_off1;
//             off001100 += half_targ_off1;
//             off001101 += half_targ_off1;
//             off001110 += half_targ_off1;
//             off001111 += half_targ_off1;
//             off010000 += half_targ_off1;
//             off010001 += half_targ_off1;
//             off010011 += half_targ_off1;
//             off010100 += half_targ_off1;
//             off010101 += half_targ_off1;
//             off010110 += half_targ_off1;
//             off010111 += half_targ_off1;
//             off011000 += half_targ_off1;
//             off011001 += half_targ_off1;
//             off011010 += half_targ_off1;
//             off011100 += half_targ_off1;
//             off011101 += half_targ_off1;
//             off011110 += half_targ_off1;
//             off011111 += half_targ_off1;
//             off100000 += half_targ_off1;
//             off100001 += half_targ_off1;
//             off100010 += half_targ_off1;
//             off100011 += half_targ_off1;
//             off100101 += half_targ_off1;
//             off100110 += half_targ_off1;
//             off100111 += half_targ_off1;
//             off101000 += half_targ_off1;
//             off101001 += half_targ_off1;
//             off101010 += half_targ_off1;
//             off101011 += half_targ_off1;
//             off101100 += half_targ_off1;
//             off101110 += half_targ_off1;
//             off101111 += half_targ_off1;
//             off110000 += half_targ_off1;
//             off110001 += half_targ_off1;
//             off110010 += half_targ_off1;
//             off110011 += half_targ_off1;
//             off110100 += half_targ_off1;
//             off110101 += half_targ_off1;
//             off110111 += half_targ_off1;
//             off111000 += half_targ_off1;
//             off111001 += half_targ_off1;
//             off111010 += half_targ_off1;
//             off111011 += half_targ_off1;
//             off111100 += half_targ_off1;
//             off111101 += half_targ_off1;
//             off111110 += half_targ_off1;
//         }
//         off000001 += half_targ_off2;
//         off000010 += half_targ_off2;
//         off000011 += half_targ_off2;
//         off000100 += half_targ_off2;
//         off000101 += half_targ_off2;
//         off000110 += half_targ_off2;
//         off000111 += half_targ_off2;
//         off001000 += half_targ_off2;
//         off001010 += half_targ_off2;
//         off001011 += half_targ_off2;
//         off001100 += half_targ_off2;
//         off001101 += half_targ_off2;
//         off001110 += half_targ_off2;
//         off001111 += half_targ_off2;
//         off010000 += half_targ_off2;
//         off010001 += half_targ_off2;
//         off010011 += half_targ_off2;
//         off010100 += half_targ_off2;
//         off010101 += half_targ_off2;
//         off010110 += half_targ_off2;
//         off010111 += half_targ_off2;
//         off011000 += half_targ_off2;
//         off011001 += half_targ_off2;
//         off011010 += half_targ_off2;
//         off011100 += half_targ_off2;
//         off011101 += half_targ_off2;
//         off011110 += half_targ_off2;
//         off011111 += half_targ_off2;
//         off100000 += half_targ_off2;
//         off100001 += half_targ_off2;
//         off100010 += half_targ_off2;
//         off100011 += half_targ_off2;
//         off100101 += half_targ_off2;
//         off100110 += half_targ_off2;
//         off100111 += half_targ_off2;
//         off101000 += half_targ_off2;
//         off101001 += half_targ_off2;
//         off101010 += half_targ_off2;
//         off101011 += half_targ_off2;
//         off101100 += half_targ_off2;
//         off101110 += half_targ_off2;
//         off101111 += half_targ_off2;
//         off110000 += half_targ_off2;
//         off110001 += half_targ_off2;
//         off110010 += half_targ_off2;
//         off110011 += half_targ_off2;
//         off110100 += half_targ_off2;
//         off110101 += half_targ_off2;
//         off110111 += half_targ_off2;
//         off111000 += half_targ_off2;
//         off111001 += half_targ_off2;
//         off111010 += half_targ_off2;
//         off111011 += half_targ_off2;
//         off111100 += half_targ_off2;
//         off111101 += half_targ_off2;
//         off111110 += half_targ_off2;
//     }
// }

// VSWAP_Gate_3_3_MEM::VSWAP_Gate_3_3_MEM(vector<int> targ): Gate(targ) {
//     swap_out = {targ[0], targ[1], targ[2]};
//     swap_in  = {targ[3], targ[4], targ[5]};

//     for (int i = 0; i < 3; i++){
//         half_targ_off[i] = env.qubit_offset[swap_out[i]];
//         targ_off[i] = half_targ_off[i] * 2;
//     }
    
//     half_targ_off[3] = env.chunk_state;
//     half_targ_off[4] = 2 * env.chunk_state;
//     half_targ_off[5] = 4 * env.chunk_state;

//     targ_off[3] = half_targ_off[3] * 2;
//     targ_off[4] = half_targ_off[4] * 2;
//     targ_off[5] = half_targ_off[5] * 2;
    
//     for (int i = 0; i < 64; i++){
//         int off = 0;
//         for (int j = 0; j < 6; j++){
//             if (i & (1 << j))
//                 off += half_targ_off[j];
//         }
//         init_off[i] = off;
//     }

//     int cur = 0;
//     for (int i = 0; i < 64; i++) {
//         if ((i >> 3) != (i % 8)) {
//             int n1 = i;
//             int n2 = ((i % 8) << 3)|(i >> 3);
//             if (n1 < n2) {
//                 swap_list[cur] = n1;
//                 dest_list[cur] = n2;
//                 cur++;
//             }
//         }
//     }

//     _run = bind(&VSWAP_Gate_3_3_MEM::run_nonchunk_chunk, this, _1, _2);
//     name += "VSWAP_Gate_3_3_MEM";
// }

// void VSWAP_Gate_3_3_MEM::run_nonchunk_chunk(vector<complex<double>> &buffer, long long idx){
//     complex<double> temp;
//     array<long long, 64> off = init_off;

//     for (auto &p : off) {
//         p += idx;
//     }

//     for (int i = 0; i < env.chunk_state; i += targ_off[2]) {
//         for (int j = 0; j < half_targ_off[2]; j += targ_off[1]) {
//             for (int k = 0; k < half_targ_off[1]; k += targ_off[0]) {
//                 for (int m = 0; m < half_targ_off[0]; m++) {
//                     for (int p = 0; p < 28; p++) {
//                         temp = buffer[swap_list[p]];
//                         buffer[swap_list[p]] = buffer[dest_list[p]];
//                         buffer[dest_list[p]] = temp;
//                     }
//                     for (auto &p : off) {
//                         p++;
//                     }
//                 }
//                 for (auto &p : off) {
//                     p += half_targ_off[0];
//                 }
//             }
//             for (auto &p : off) {
//                 p += half_targ_off[1];
//             }
//         }
//         for (auto &p : off) {
//             p += half_targ_off[2];
//         }
//     }
// }

VSWAP_Gate_3_3_MEM::VSWAP_Gate_3_3_MEM(vector<int> targ): Gate(targ) {
    type = VSWAP;
    half_targ_off[0] = env.qubit_offset[targ[0]];
    half_targ_off[1] = env.qubit_offset[targ[1]];
    half_targ_off[2] = env.qubit_offset[targ[2]];
    half_targ_off[3] = env.chunk_state;
    half_targ_off[4] = 2 * env.chunk_state;
    half_targ_off[5] = 4 * env.chunk_state;

    targ_off[0] = half_targ_off[0] * 2;
    targ_off[1] = half_targ_off[1] * 2;
    targ_off[2] = half_targ_off[2] * 2;

    int cur = 0;
    for (int i = 0; i < 64; i++) {
        if ((i >> 3) != (i % 8)) {
            int n1 = i;
            int n2 = ((i % 8) << 3)|(i >> 3);
            if (n1 < n2) {
                swap_list[cur] = n1;
                dest_list[cur] = n2;

                int n1_off = 0;
                int n2_off = 0;
                for (int j = 0; j < 6; j++){
                    if (n1 & (1 << j))
                        n1_off += half_targ_off[j];
                    if (n2 & (1 << j))
                        n2_off += half_targ_off[j];
                }
                init_off[cur * 2] = n1_off;
                init_off[cur * 2 + 1] = n2_off;
                cur++;
            }
        }
    }

    _run = bind(&VSWAP_Gate_3_3_MEM::run_nonchunk_chunk, this, _1, _2);
    name += "VSWAP_Gate_3_3_MEM";
}

void VSWAP_Gate_3_3_MEM::run_nonchunk_chunk(vector<complex<double>> &buffer, long long idx){
    complex<double> temp;
    long long off[56];
    #pragma GCC unroll(56)
    for (int p = 0; p < 56; p++) {
        off[p] = init_off[p] + idx;
    }
    
    for (int i = 0; i < env.chunk_state; i += targ_off[2]) {
        for (int j = 0; j < half_targ_off[2]; j += targ_off[1]) {
            for (int k = 0; k < half_targ_off[1]; k += targ_off[0]) {
                for (int m = 0; m < half_targ_off[0]; m++) {

                    #pragma GCC unroll 7
                    // #pragma GCC unroll 28
                    // #pragma omp unroll full
                    for (int p = 0; p < 28; p++) {
                        temp = buffer[swap_list[p]];
                        buffer[swap_list[p]] = buffer[dest_list[p]];
                        buffer[dest_list[p]] = temp;
                    }
                    
                    #pragma GCC unroll 7
                    // #pragma GCC unroll 56
                    // #pragma omp unroll full
                    for (int p = 0; p < 56; p++) {
                        off[p]++;
                    }
                }
                #pragma GCC unroll 7
                // #pragma GCC unroll 56
                // #pragma omp unroll full
                for (int p = 0; p < 56; p++) {
                    off[p] += half_targ_off[0];
                }
            }
            #pragma GCC unroll 7
            // #pragma GCC unroll 56
            // #pragma omp unroll full
            for (int p = 0; p < 56; p++) {
                off[p] += half_targ_off[1];
            }
        }
        #pragma GCC unroll 7
        // #pragma GCC unroll 56
        // #pragma omp unroll full
        for (int p = 0; p < 56; p++) {
            off[p] += half_targ_off[2];
        }
    }
}

VSWAP_Gate_4_4_MEM::VSWAP_Gate_4_4_MEM(vector<int> targ): Gate(targ) {
    type = VSWAP;
    swap_out = {targ[0], targ[1], targ[2], targ[3]};
    swap_in  = {targ[4], targ[5], targ[6], targ[7]};

    for (int i = 0; i < 4; i++){
        half_targ_off[i] = env.qubit_offset[swap_out[i]];
        targ_off[i] = half_targ_off[i] * 2;
    }
    
    half_targ_off[4] = env.chunk_state;
    targ_off[4] = half_targ_off[4] * 2;

    for (int i = 5; i < 8; i++) {
        half_targ_off[i] = half_targ_off[i-1] * 2;
        targ_off[i] = half_targ_off[i] * 2;
    }
    
    for (int i = 0; i < 256; i++){
        int off = 0;
        for (int j = 0; j < 8; j++){
            if (i & (1 << j))
                off += half_targ_off[j];
        }
        init_off[i] = off;
    }

    int cur = 0;
    for (int i = 0; i < 256; i++) {
        if ((i>>4) != (i%16)) {
            int n1 = i;
            int n2 = ((i%16)<<4)|(i>>4);
            if (n1 < n2) {
                swap_list[cur] = n1;
                dest_list[cur] = n2;
                cur++;
            }
        }
    }

    _run = bind(&VSWAP_Gate_4_4_MEM::run_nonchunk_chunk, this, _1, _2);
    name += "VSWAP_Gate_4_4_MEM";
}

void VSWAP_Gate_4_4_MEM::run_nonchunk_chunk(vector<complex<double>> &buffer, long long idx){
    complex<double> temp;
    array<long long, 256> off = init_off;

    for (auto &p : off) {
        p += idx;
    }

    for (int i = 0; i < env.chunk_state; i += targ_off[3]) {
        for (int j = 0; j < half_targ_off[3]; j += targ_off[2]) {
            for (int k = 0; k < half_targ_off[2]; k += targ_off[1]) {
                for (int m = 0; m < half_targ_off[1]; m += targ_off[0]) {
                    for (int n = 0; n < half_targ_off[0]; n++) {
                        for (int p = 0; p < 120; p++) {
                            temp = buffer[swap_list[p]];
                            buffer[swap_list[p]] = buffer[dest_list[p]];
                            buffer[dest_list[p]] = temp;
                        }
                        for (auto &p : off) {
                            p++;
                        }
                    }
                    for (auto &p : off) {
                        p += half_targ_off[0];
                    }
                }
                for (auto &p : off) {
                    p += half_targ_off[1];
                }
            }
            for (auto &p : off) {
                p += half_targ_off[2];
            }
        }
        for (auto &p : off) {
            p += half_targ_off[3];
        }
    }
}

VSWAP_Gate_6_6_MEM::VSWAP_Gate_6_6_MEM(vector<int> targ): Gate(targ) {
    type = VSWAP;
    swap_out = {targ[0], targ[1], targ[2], targ[3], targ[4], targ[5]};
    swap_in  = {targ[6], targ[7], targ[8], targ[9], targ[10], targ[11]};

    for (int i = 0; i < 6; i++)
        half_targ_off[i] = env.qubit_offset[swap_out[i]];
    
    half_targ_off[6] = env.chunk_state;
    for (int i = 7; i < 12; i++)
        half_targ_off[i] = half_targ_off[i-1] * 2;
    
    for (int i = 0; i < 12; i++)
        targ_off[i] = half_targ_off[i] * 2;
    
    for (int i = 0; i < 4096; i++){
        int off = 0;
        for (int j = 0; j < 12; j++){
            if (i & (1 << j))
                off += half_targ_off[j];
        }
        init_off[i] = off;
    }

    int cur = 0;
    for (int i = 0; i < 4096; i++){
        if ((i>>6) != (i%64)) {
            int n1 = i;
            int n2 = ((i%64)<<6)|(i>>6);
            if (n1 < n2) {
                swap_list[cur] = n1;
                dest_list[cur] = n2;
                cur++;
            }
        }
    }

    _run = bind(&VSWAP_Gate_6_6_MEM::run_nonchunk_chunk, this, _1, _2);
    name += "VSWAP_Gate_6_6_MEM";
}

void VSWAP_Gate_6_6_MEM::run_nonchunk_chunk(vector<complex<double>> &buffer, long long idx){
    complex<double> temp;
    array<long long, 4096> off = init_off;

    for (auto & i : off) {
        i += idx;
    }

    for (int i = 0; i < env.chunk_state; i += targ_off[5]) {
        for (int j = 0; j < half_targ_off[5]; j += targ_off[4]) {
            for (int k = 0; k < half_targ_off[4]; k += targ_off[3]) {
                for (int m = 0; m < half_targ_off[3]; m += targ_off[2]) {
                    for (int n = 0; n < half_targ_off[2]; n += targ_off[1]) {
                        for (int r = 0; r < half_targ_off[1]; r += targ_off[0]) {
                            for (int s = 0; s < half_targ_off[0]; s++) {
                                for (int p = 0; p < 2016; p++) { // (2^12 - 2^6) / 2 = 2016
                                    temp = buffer[swap_list[p]];
                                    buffer[swap_list[p]] = buffer[dest_list[p]];
                                    buffer[dest_list[p]] = temp;
                                }
                                for (auto &p : off) {
                                    p++;
                                }
                            }
                            for (auto &p : off) {
                                p += half_targ_off[0];
                            }
                        }
                        for (auto &p : off) {
                            p += half_targ_off[1];
                        }
                    }
                    for (auto &p : off) {
                        p += half_targ_off[2];
                    }
                }
                for (auto &p : off) {
                    p += half_targ_off[3];
                }
            }
            for (auto &p : off) {
                p += half_targ_off[4];
            }
        }
        for (auto &p : off) {
            p += half_targ_off[5];
        }
    }
}
