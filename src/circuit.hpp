#ifndef CIRCUIT_HPP
#define CIRCUIT_HPP
#include <cmath>
#include <vector>
#include <complex>
#include <functional>
#include <algorithm>
using namespace std;

typedef enum gate_type{
    MEASURE = 0,
    ONE_QUBIT = 1,
    TWO_QUBIT = 2,
    THREE_QUBIT = 3,
    VSWAP = 4,
    ERROR = 5,
} GATE_TYPE;

class Gate
{
public:
    GATE_TYPE type;
    int mpi_count = 0;
    int file_count = 0;
    int middle_count = 0;
    int chunk_count = 0;
    int nonchunk_count = 0;
    vector<int> targs;


    string name;
    function<void (vector<complex<double>>& buffer)> run;
    function<void (vector<complex<double>>*buffer1,vector<complex<double>>*buffer2,int off0,int off1,int round)> run_one_qubit_mpi;
    function<void (complex<double> *buffer)> run_dio;  // DIO version
    function<void (vector<complex<double>>& buffer, long long idx)> _run; // MEM version
    Gate(vector<int>);
    virtual ~Gate(){};
};

/*-----ONE_QUBIT_GATE-----*/

class ONE_QUBIT_GATE: public Gate
{
public:
    int half_off;
    int off;
    ONE_QUBIT_GATE(vector<int>);
};

class H_Gate: public ONE_QUBIT_GATE
{
public:
    // static double coeff = 0.70710678118;    // 1/sqrt(2) = 0.70710678118
    H_Gate(vector<int>);
    void run_chunk(vector<complex<double>> &);
    void run_nonchunk(vector<complex<double>> &);
    void run_mpi_chunk(vector<complex<double>>*,vector<complex<double>>*,int,int,int);
};

class X_Gate: public ONE_QUBIT_GATE
{
public:
    X_Gate(vector<int>);
    void run_chunk(vector<complex<double>> &);
    void run_nonchunk(vector<complex<double>> &);
    void run_mpi_chunk(vector<complex<double>>*,vector<complex<double>>*,int,int,int);
};

class Y_Gate: public ONE_QUBIT_GATE
{
public:
    Y_Gate(vector<int>);
    void run_chunk(vector<complex<double>> &);
    void run_nonchunk(vector<complex<double>> &);
    void run_mpi_chunk(vector<complex<double>>*,vector<complex<double>>*,int,int,int);
};

class Z_Gate: public ONE_QUBIT_GATE
{
public:
    Z_Gate(vector<int>);
    void run_chunk(vector<complex<double>> &);
    void run_nonchunk(vector<complex<double>> &);
};

class U1_Gate: public ONE_QUBIT_GATE
{
public:
    vector<complex<double>> coeff;
    U1_Gate(vector<int>, vector<complex<double>>);
    void run_chunk(vector<complex<double>> &);
    void run_nonchunk(vector<complex<double>> &);
    void run_mpi_chunk(vector<complex<double>>*,vector<complex<double>>*,int,int,int);
};

class Phase_Gate: public ONE_QUBIT_GATE
{
public:
    complex<double> phi;
    Phase_Gate(vector<int>, double);
    void run_chunk(vector<complex<double>> &);
    void run_nonchunk(vector<complex<double>> &);
};

class RX_Gate: public ONE_QUBIT_GATE
{
public:
    complex<double> cos_Phi_2;
    complex<double> i_sin_Phi_2;

    RX_Gate(vector<int>, double);
    void run_chunk(vector<complex<double>> &);
    void run_nonchunk(vector<complex<double>> &);
    void run_mpi_chunk(vector<complex<double>>*,vector<complex<double>>*,int,int,int);
};

class RY_Gate: public ONE_QUBIT_GATE
{
public:
    complex<double> cos_Phi_2;
    complex<double> sin_Phi_2;

    RY_Gate(vector<int>, double);
    void run_chunk(vector<complex<double>> &);
    void run_nonchunk(vector<complex<double>> &);
    void run_mpi_chunk(vector<complex<double>>*,vector<complex<double>>*,int,int,int);
};

class RZ_Gate: public ONE_QUBIT_GATE
{
public:
    complex<double> exp_p_iPhi_2;
    complex<double> exp_n_iPhi_2;

    RZ_Gate(vector<int>, double);
    void run_chunk(vector<complex<double>> &);
    void run_nonchunk(vector<complex<double>> &);
};

/*-----TWO_QUBIT_GATE-----*/
class TWO_QUBIT_GATE: public Gate{
public:
    int half_off_0;
    int half_off_1;

    int off_0;
    int off_1;
    TWO_QUBIT_GATE(vector<int>);
};

class U2_Gate: public TWO_QUBIT_GATE
{
public:
    complex<double> coeff[16];
    U2_Gate(vector<int>, vector<complex<double>>);
    void run_chunk_chunk(vector<complex<double>> &);
    void run_nonchunk_chunk(vector<complex<double>> &);
    void run_nonchunk_nonchunk(vector<complex<double>> &);
};

class SWAP_Gate: public TWO_QUBIT_GATE
{
public:
    SWAP_Gate(vector<int>);
    void run_chunk_chunk(vector<complex<double>> &);
    void run_nonchunk_chunk(vector<complex<double>> &);
    void run_nonchunk_nonchunk(vector<complex<double>> &);
    void run_mpi_nonchunk_chunk(vector<complex<double>>*,vector<complex<double>>*,int,int,int);
};

class VSWAP_Gate_1_1: public SWAP_Gate
{
public:
    VSWAP_Gate_1_1(vector<int>);
};

class CPhase_Gate: public TWO_QUBIT_GATE
{
public:
    complex<double> exp_iPhi;

    CPhase_Gate(vector<int>, double);
    void run_chunk_chunk(vector<complex<double>> &);
    void run_nonchunk_chunk(vector<complex<double>> &);
    void run_nonchunk_nonchunk(vector<complex<double>> &);
};

class RZZ_Gate: public TWO_QUBIT_GATE
{
public:
    complex<double> exp_p_iPhi_2;
    complex<double> exp_n_iPhi_2;
    
    RZZ_Gate(vector<int>, double);
    void run_chunk_chunk(vector<complex<double>> &);
    void run_nonchunk_chunk(vector<complex<double>> &);
    void run_nonchunk_nonchunk(vector<complex<double>> &);
};

/*-----THREE_QUBIT_GATE-----*/
class U3_Gate: public Gate
{
public:
    int half_target_offset0;
    int half_target_offset1;
    int half_target_offset2;

    int target_offset0;
    int target_offset1;
    int target_offset2;

    vector<complex<double>> coeff;
    U3_Gate(vector<int>, vector<complex<double>>);
    void run_chunk_chunk_chunk(vector<complex<double>> &);
    void run_nonchunk_chunk_chunk(vector<complex<double>> &);
    void run_nonchunk_nonchunk_chunk(vector<complex<double>> &);
    void run_nonchunk_nonchunk_nonchunk(vector<complex<double>> &);
};

/*-----VSWAP_GATE-----*/
// Limit: SWAP that must be half in chunk and half out of chunk.
// [TODO]: findout the possibility if it can be template class.

class VSWAP_Gate_2_2: public Gate
{
public:
    int half_off_0;
    int half_off_1;
    int half_off_2;
    int half_off_3;

    int off_0;
    int off_1;
    int off_2;
    int off_3;

    VSWAP_Gate_2_2(vector<int>);
    void run_nonchunk_chunk(vector<complex<double>> &);
};

class VSWAP_Gate_3_3: public Gate
{
public:
    int half_targ_off0;
    int half_targ_off1;
    int half_targ_off2;
    int half_targ_off3;
    int half_targ_off4;
    int half_targ_off5;

    int targ_off0;
    int targ_off1;
    int targ_off2;
    int targ_off3;
    int targ_off4;
    int targ_off5;
    
    vector<int> swap_in;
    vector<int> swap_out;

    VSWAP_Gate_3_3(vector<int>);
    void run_nonchunk_chunk(vector<complex<double>> &);
};

class VSWAP_Gate_4_4: public Gate
{
public:
    vector<int> half_targ_off;
    vector<int> targ_off;

    array<int, 256> init_off;
    vector<int> swap_list;
    vector<int> dest_list;
    
    vector<int> swap_in;
    vector<int> swap_out;

    VSWAP_Gate_4_4(vector<int>);
    void run_nonchunk_chunk(vector<complex<double>> &);
};

class VSWAP_Gate_6_6: public Gate
{
public:
    vector<int> half_targ_off;
    vector<int> targ_off;

    array<int, 4096> init_off;
    vector<int> swap_list;
    vector<int> dest_list;
    
    vector<int> swap_in;
    vector<int> swap_out;

    VSWAP_Gate_6_6(vector<int>);
    void run_nonchunk_chunk(vector<complex<double>> &);
};

/*DIO version*/
/*-----ONE_QUBIT_GATE-----*/

class H_Gate_DIO: public ONE_QUBIT_GATE
{
public:
    // static double coeff = 0.70710678118;    // 1/sqrt(2) = 0.70710678118
    H_Gate_DIO(vector<int>);
    void run_chunk(complex<double> *);
    void run_nonchunk(complex<double> *);
};

class X_Gate_DIO: public ONE_QUBIT_GATE
{
public:
    X_Gate_DIO(vector<int>);
    void run_chunk(complex<double> *);
    void run_nonchunk(complex<double> *);
};

class Y_Gate_DIO: public ONE_QUBIT_GATE
{
public:
    Y_Gate_DIO(vector<int>);
    void run_chunk(complex<double> *);
    void run_nonchunk(complex<double> *);
};

class Z_Gate_DIO: public ONE_QUBIT_GATE
{
public:
    Z_Gate_DIO(vector<int>);
    void run_chunk(complex<double> *);
    void run_nonchunk(complex<double> *);
};

class U1_Gate_DIO: public ONE_QUBIT_GATE
{
public:
    vector<complex<double>> coeff;
    U1_Gate_DIO(vector<int>, vector<complex<double>>);
    void run_chunk(complex<double> *);
    void run_nonchunk(complex<double> *);
};

class Phase_Gate_DIO: public ONE_QUBIT_GATE
{
public:
    complex<double> phi;
    Phase_Gate_DIO(vector<int>, double);
    void run_chunk(complex<double> *);
    void run_nonchunk(complex<double> *);
};

class RX_Gate_DIO: public ONE_QUBIT_GATE
{
public:
    complex<double> cos_Phi_2;
    complex<double> i_sin_Phi_2;

    RX_Gate_DIO(vector<int>, double);
    void run_chunk(complex<double> *);
    void run_nonchunk(complex<double> *);
};

class RY_Gate_DIO: public ONE_QUBIT_GATE
{
public:
    complex<double> cos_Phi_2;
    complex<double> sin_Phi_2;

    RY_Gate_DIO(vector<int>, double);
    void run_chunk(complex<double> *);
    void run_nonchunk(complex<double> *);
};

class RZ_Gate_DIO: public ONE_QUBIT_GATE
{
public:
    complex<double> exp_p_iPhi_2;
    complex<double> exp_n_iPhi_2;

    RZ_Gate_DIO(vector<int>, double);
    void run_chunk(complex<double> *);
    void run_nonchunk(complex<double> *);
};

/*-----TWO_QUBIT_GATE-----*/
class U2_Gate_DIO: public TWO_QUBIT_GATE
{
public:
    complex<double> coeff[16];
    U2_Gate_DIO(vector<int>, vector<complex<double>>);
    void run_chunk_chunk(complex<double> *);
    void run_nonchunk_chunk(complex<double> *);
    void run_nonchunk_nonchunk(complex<double> *);
};

class SWAP_Gate_DIO: public TWO_QUBIT_GATE
{
public:
    SWAP_Gate_DIO(vector<int>);
    void run_chunk_chunk(complex<double> *);
    void run_nonchunk_chunk(complex<double> *);
    void run_nonchunk_nonchunk(complex<double> *);
};

class VSWAP_Gate_1_1_DIO: public SWAP_Gate_DIO
{
public:
    VSWAP_Gate_1_1_DIO(vector<int>);
};

class CPhase_Gate_DIO: public TWO_QUBIT_GATE
{
public:
    complex<double> exp_iPhi;

    CPhase_Gate_DIO(vector<int>, double);
    void run_chunk_chunk(complex<double> *);
    void run_nonchunk_chunk(complex<double> *);
    void run_nonchunk_nonchunk(complex<double> *);
};

class RZZ_Gate_DIO: public TWO_QUBIT_GATE
{
public:
    complex<double> exp_p_iPhi_2;
    complex<double> exp_n_iPhi_2;
    
    RZZ_Gate_DIO(vector<int>, double);
    void run_chunk_chunk(complex<double> *);
    void run_nonchunk_chunk(complex<double> *);
    void run_nonchunk_nonchunk(complex<double> *);
};

/*-----THREE_QUBIT_GATE-----*/
class U3_Gate_DIO: public Gate
{
public:
    int half_target_offset0;
    int half_target_offset1;
    int half_target_offset2;

    int target_offset0;
    int target_offset1;
    int target_offset2;

    vector<complex<double>> coeff;
    U3_Gate_DIO(vector<int>, vector<complex<double>>);
    void run_chunk_chunk_chunk(complex<double> *);
    void run_nonchunk_chunk_chunk(complex<double> *);
    void run_nonchunk_nonchunk_chunk(complex<double> *);
    void run_nonchunk_nonchunk_nonchunk(complex<double> *);
};

/*-----VSWAP_GATE-----*/
// Limit: SWAP that must be half in chunk and half out of chunk.
// [TODO]: findout the possibility if it can be template class.

class VSWAP_Gate_2_2_DIO: public Gate
{
public:
    int half_off_0;
    int half_off_1;
    int half_off_2;
    int half_off_3;

    int off_0;
    int off_1;
    int off_2;
    int off_3;

    VSWAP_Gate_2_2_DIO(vector<int>);
    void run_nonchunk_chunk(complex<double> *);
};

/*MEM version*/
/*-----ONE_QUBIT_GATE-----*/
class ONE_QUBIT_GATE_MEM: public Gate
{
public:
    long long half_off;
    long long off;
    ONE_QUBIT_GATE_MEM(vector<int>);
};

class H_Gate_MEM: public ONE_QUBIT_GATE_MEM
{
public:
    H_Gate_MEM(vector<int>);
    void run_chunk(vector<complex<double>> &buffer, long long idx);
    void run_nonchunk(vector<complex<double>> &buffer, long long idx);
};

class X_Gate_MEM: public ONE_QUBIT_GATE_MEM
{
public:
    X_Gate_MEM(vector<int>);
    void run_chunk(vector<complex<double>> &buffer, long long idx);
    void run_nonchunk(vector<complex<double>> &buffer, long long idx);
};

class Y_Gate_MEM: public ONE_QUBIT_GATE_MEM
{
public:
    Y_Gate_MEM(vector<int>);
    void run_chunk(vector<complex<double>> &buffer, long long idx);
    void run_nonchunk(vector<complex<double>> &buffer, long long idx);
};

class Z_Gate_MEM: public ONE_QUBIT_GATE_MEM
{
public:
    Z_Gate_MEM(vector<int>);
    void run_chunk(vector<complex<double>> &buffer, long long idx);
    void run_nonchunk(vector<complex<double>> &buffer, long long idx);
};

class Phase_Gate_MEM: public ONE_QUBIT_GATE_MEM
{
public:
    complex<double> phi;
    Phase_Gate_MEM(vector<int>, double);
    void run_chunk(vector<complex<double>> &buffer, long long idx);
    void run_nonchunk(vector<complex<double>> &buffer, long long idx);
};

class U1_Gate_MEM: public ONE_QUBIT_GATE_MEM
{
public:
    vector<complex<double>> coeff;
    U1_Gate_MEM(vector<int>, vector<complex<double>>);
    void run_chunk(vector<complex<double>> &buffer, long long idx);
    void run_nonchunk(vector<complex<double>> &buffer, long long idx);
};

class RX_Gate_MEM: public ONE_QUBIT_GATE_MEM
{
public:
    complex<double> cos_Phi_2;
    complex<double> i_sin_Phi_2;
    RX_Gate_MEM(vector<int>, double);
    void run_chunk(vector<complex<double>> &, long long);
    void run_nonchunk(vector<complex<double>> &, long long);
};

class RY_Gate_MEM: public ONE_QUBIT_GATE_MEM
{
public:
    complex<double> cos_Phi_2;
    complex<double> sin_Phi_2;
    RY_Gate_MEM(vector<int>, double);
    void run_chunk(vector<complex<double>> &, long long);
    void run_nonchunk(vector<complex<double>> &, long long);
};

class RZ_Gate_MEM: public ONE_QUBIT_GATE_MEM
{
public:
    complex<double> exp_p_iPhi_2;
    complex<double> exp_n_iPhi_2;
    RZ_Gate_MEM(vector<int>, double);
    void run_chunk(vector<complex<double>> &, long long);
    void run_nonchunk(vector<complex<double>> &, long long);
};

/*-----TWO_QUBIT_GATE-----*/
class TWO_QUBIT_GATE_MEM: public Gate
{
public:
    long long half_off_0;
    long long half_off_1;
    long long off_0;
    long long off_1;
    TWO_QUBIT_GATE_MEM(vector<int>);
};

class U2_Gate_MEM: public TWO_QUBIT_GATE_MEM
{
public:
    vector<complex<double>> coeff;
    U2_Gate_MEM(vector<int>, vector<complex<double>>);
    void run_chunk_chunk(vector<complex<double>> &buffer, long long idx);
    void run_nonchunk_chunk(vector<complex<double>> &buffer, long long idx);
    void run_nonchunk_nonchunk(vector<complex<double>> &buffer, long long idx);
};

class SWAP_Gate_MEM: public TWO_QUBIT_GATE_MEM
{
public:
    SWAP_Gate_MEM(vector<int>);
    void run_chunk_chunk(vector<complex<double>> &, long long idx);
    void run_nonchunk_chunk(vector<complex<double>> &, long long idx);
    void run_nonchunk_nonchunk(vector<complex<double>> &, long long idx);
};


class VSWAP_Gate_1_1_MEM: public SWAP_Gate_MEM
{
public:
    VSWAP_Gate_1_1_MEM(vector<int>);
};

class CPhase_Gate_MEM: public TWO_QUBIT_GATE_MEM
{
public:
    complex<double> exp_iPhi;
    CPhase_Gate_MEM(vector<int>, double);
    void run_chunk_chunk(vector<complex<double>> &, long long idx);
    void run_nonchunk_chunk(vector<complex<double>> &, long long idx);
    void run_nonchunk_nonchunk(vector<complex<double>> &, long long idx);
};

class RZZ_Gate_MEM: public TWO_QUBIT_GATE_MEM
{
public:
    complex<double> exp_p_iPhi_2;
    complex<double> exp_n_iPhi_2;
    RZZ_Gate_MEM(vector<int>, double);
    void run_chunk_chunk(vector<complex<double>> &, long long);
    void run_nonchunk_chunk(vector<complex<double>> &, long long);
    void run_nonchunk_nonchunk(vector<complex<double>> &, long long);
};

/*-----THREE_QUBIT_GATE-----*/
class U3_Gate_MEM: public Gate
{
public:
    long long half_target_offset0;
    long long half_target_offset1;
    long long half_target_offset2;
    long long target_offset0;
    long long target_offset1;
    long long target_offset2;

    vector<complex<double>> coeff;
    U3_Gate_MEM(vector<int>, vector<complex<double>>);
    void run_chunk_chunk_chunk(vector<complex<double>> &buffer, long long idx);
    void run_nonchunk_chunk_chunk(vector<complex<double>> &buffer, long long idx);
    void run_nonchunk_nonchunk_chunk(vector<complex<double>> &buffer, long long idx);
    void run_nonchunk_nonchunk_nonchunk(vector<complex<double>> &buffer, long long idx);
};

/*-----VSWAP_GATE-----*/
// Limit: SWAP that must be half in chunk and half out of chunk.
// [TODO]: findout the possibility if it can be template class.

class VSWAP_Gate_2_2_MEM: public Gate
{
public:
    long long half_target_offset0;
    long long half_target_offset1;
    long long half_target_offset2;
    long long half_target_offset3;

    long long target_offset0;
    long long target_offset1;
    long long target_offset2;
    long long target_offset3;
    
    vector<int> swap_in;
    vector<int> swap_out;

    VSWAP_Gate_2_2_MEM(vector<int>);
    void run_nonchunk_chunk(vector<complex<double>> &, long long);
};

// class VSWAP_Gate_3_3_MEM: public Gate
// {
// public:
//     long long half_targ_off0;
//     long long half_targ_off1;
//     long long half_targ_off2;
//     long long half_targ_off3;
//     long long half_targ_off4;
//     long long half_targ_off5;

//     long long targ_off0;
//     long long targ_off1;
//     long long targ_off2;
//     long long targ_off3;
//     long long targ_off4;
//     long long targ_off5;
    
//     array<int, 3> swap_in;
//     array<int, 3> swap_out;

//     VSWAP_Gate_3_3_MEM(vector<int>);
//     void run_nonchunk_chunk(vector<complex<double>> &, long long);
// };

class VSWAP_Gate_3_3_MEM: public Gate
{
public:
    long long half_targ_off[6];
    long long targ_off[3];

    long long init_off[56];
    int swap_list[28];
    int dest_list[28];

    VSWAP_Gate_3_3_MEM(vector<int>);
    void run_nonchunk_chunk(vector<complex<double>> &, long long);
};

class VSWAP_Gate_4_4_MEM: public Gate
{
public:
    array<long long, 8> half_targ_off;
    array<long long, 8> targ_off;

    array<long long, 256> init_off;
    array<int, 120> swap_list;
    array<int, 120> dest_list;
    
    array<int, 4> swap_in;
    array<int, 4> swap_out;

    VSWAP_Gate_4_4_MEM(vector<int>);
    void run_nonchunk_chunk(vector<complex<double>> &, long long);
};

class VSWAP_Gate_6_6_MEM: public Gate
{
public:
    array<long long, 6> half_targ_off;
    array<long long, 6> targ_off;

    array<long long, 4096> init_off;
    array<int, 2016> swap_list;
    array<int, 2016> dest_list;
    
    array<int, 6> swap_in;
    array<int, 6> swap_out;

    VSWAP_Gate_6_6_MEM(vector<int>);
    void run_nonchunk_chunk(vector<complex<double>> &, long long);
};
#endif