from qiskit import Aer, QuantumCircuit, transpile
from qiskit.extensions import UnitaryGate
from ini_generator import *
from math import sqrt
import numpy as np
import sys
import os
import struct
import re
import time

class verifier:
    def __init__(self, args:Args, circuit_path) -> None:
        self.args = args
        self.circuit_path = circuit_path
        if args.runner_type == "IO":
            self.file_paths = self.args.state_paths.split(",")
            self.run = self.run_IO
        #     self.loader = loader_IO
        elif args.runner_type == "MEM":
            self.run = self.run_MEM
        #     self.loader = loader_MEM

    # [Qiskit]: init to 0 state
    def qiskit_circ_init(self):
        circ = QuantumCircuit(self.args.total_qbit)
        initial_state = [1,0]   # Define initial_state as |0>
        circ.initialize(initial_state, 0)
        return circ

    # [Qiskit]: state_vector
    def qiskit_init_state_vector(self, circ):
        simulator = Aer.get_backend('aer_simulator_statevector')
        circ.save_statevector(label=f'save')
        circ = transpile(circ, simulator)

        # start_time = time.perf_counter()
        data = simulator.run(circ).result().data(0)
        # end_time = time.perf_counter()
        # print(f"Time: {(end_time - start_time)*1000000:.0f} (us)\n")

        res = np.array(data['save'])
        return res

    # [Qiskit]: density_matrix
    def qiskit_init_density_matrix(self, circ):
        simulator = Aer.get_backend('aer_simulator_density_matrix')
        circ.save_density_matrix(label=f'save')
        circ = transpile(circ, simulator)

        # start_time = time.perf_counter()
        data = simulator.run(circ).result().data(0)
        # end_time = time.perf_counter()
        # print(f"Time: {(end_time - start_time)*1000000:.0f} (us)\n")

        res = np.array(data['save'].T.reshape(-1))
        return res

    def loader_IO(self):
        N = self.args.total_qbit
        NUMFILE = len(self.file_paths)
        FILESIZE = int((1 << N) / NUMFILE)
        fd_arr = [np.zeros(FILESIZE, dtype=np.complex128) for i in range(NUMFILE)]

        for fd, state_path in enumerate(self.file_paths):
            f = fd_arr[fd]
            with open(state_path, mode="rb") as state_file:
                try:
                    state = state_file.read()
                    k = 0
                    for i in range(FILESIZE):
                        (real, imag) = struct.unpack("dd", state[k:k+16])
                        f[i] = real+imag*1j
                        k += 16
                except Exception as e:
                    print(e)
                    print(f"read from {state_path}")
                    print(f"[ERROR]: error at reading {k}th byte")
                    exit()
        return fd_arr

    def loader_MEM(self):
        FILESIZE = 1 << self.args.total_qbit
        sv = [np.zeros(FILESIZE, dtype=np.complex128)]

        with open(self.args.dump_file, mode="r") as sv_file:
            for i in range(FILESIZE):
                try:
                    real = float(sv_file.readline())
                    imag = float(sv_file.readline())
                    sv[0][i] = real + imag*1j
                except Exception as e:
                    print(e)
                    print(f"read from {self.args.dump_file}")
                    print(f"[ERROR]: error at reading {i}th byte")
                    exit()
        return sv

    def qiskit_set_circuit(self):
        circ = self.qiskit_circ_init()
        ops = []
        with open(self.circuit_path,'r') as f:
            for line in f.readlines():
                s = line.split()
                ops.append(s)

        # for op in ops:
        #     print(op)
        for op in ops:
            if op[0]=="H":
                circ.h(int(op[1]))
            if op[0]=="S":
                circ.s(int(op[1]))
            if op[0]=="T":
                circ.t(int(op[1]))
            if op[0]=="X":
                circ.x(int(op[1]))
            if op[0]=="Y":
                circ.y(int(op[1]))
            if op[0]=="Z":
                circ.z(int(op[1]))
            if op[0]=="RX":
                circ.rx(float(op[2]), int(op[1]))
            if op[0]=="RY":
                circ.ry(float(op[2]), int(op[1]))
            if op[0]=="RZ":
                circ.rz(float(op[2]), int(op[1]))
            if op[0]=="P":
                circ.p(float(op[2]), int(op[1]))
            if op[0]=="U1": # UnitaryGate
                for i in range(2, 10):
                    op[i] = float(op[i])
                gate = UnitaryGate([[op[2]+op[3]*1j, op[4]+op[5]*1j], 
                                    [op[6]+op[7]*1j, op[8]+op[9]*1j]])
                circ.append(gate,[int(op[1])])

            if op[0]=="CX":
                circ.cx(int(op[1]), int(op[2]))
            if op[0]=="CY":
                circ.cy(int(op[1]), int(op[2]))
            if op[0]=="CZ":
                circ.cy(int(op[1]), int(op[2]))
            if op[0]=="CP":
                circ.cp(float(op[3]), int(op[1]), int(op[2]))
            if op[0]=="CU1":
                # control-unitary
                for i in range(3, 11):
                    op[i] = float(op[i])
                gate = UnitaryGate([[1, 0, 0, 0],
                                    [0, 1, 0, 0],
                                    [0, 0, op[3]+op[4]*1j, op[5]+op[ 6]*1j],
                                    [0, 0, op[7]+op[8]*1j, op[9]+op[10]*1j]])
                circ.append(gate,[int(op[1]), int(op[2])])
            if op[0]=="SWAP":
                circ.swap(int(op[1]), int(op[2]))
            if op[0]=="RZZ":
                circ.rzz(float(op[3]), int(op[1]), int(op[2]))
            if op[0]=="TOFFOLI":
                # control1 control2 targert
                circ.toffoli(int(op[1]), int(op[2]), int(op[3]))

            if op[0]=="U2": # 2 qubit UnitaryGate
                for i in range(3, 35):
                    op[i] = float(op[i])
                gate = UnitaryGate([[op[ 3]+op[ 4]*1j,op[ 5]+op[ 6]*1j,op[ 7]+op[ 8]*1j,op[ 9]+op[10]*1j],
                                    [op[11]+op[12]*1j,op[13]+op[14]*1j,op[15]+op[16]*1j,op[17]+op[18]*1j],
                                    [op[19]+op[20]*1j,op[21]+op[22]*1j,op[23]+op[24]*1j,op[25]+op[26]*1j],
                                    [op[27]+op[28]*1j,op[29]+op[30]*1j,op[31]+op[32]*1j,op[33]+op[34]*1j]])
                # Qiskit treat the order in a reverse way.
                # WHAAAT?!!!WHYYYYYY?????
                # circ.append(gate,[reorder(op[4], N), reorder(op[5], N)])
                circ.append(gate,[int(op[2]), int(op[1])])

            if op[0]=="U3":  # 3 qubit UnitaryGate
                for i in range(4, 132):
                    op[i] = float(op[i])
                gate = UnitaryGate([[op[  4]+op[  5]*1j,op[  6]+op[  7]*1j,op[  8]+op[  9]*1j,op[ 10]+op[ 11]*1j,op[ 12]+op[ 13]*1j,op[ 14]+op[ 15]*1j,op[ 16]+op[ 17]*1j,op[ 18]+op[ 19]*1j],
                                    [op[ 20]+op[ 21]*1j,op[ 22]+op[ 23]*1j,op[ 24]+op[ 25]*1j,op[ 26]+op[ 27]*1j,op[ 28]+op[ 29]*1j,op[ 30]+op[ 31]*1j,op[ 32]+op[ 33]*1j,op[ 34]+op[ 35]*1j],
                                    [op[ 36]+op[ 37]*1j,op[ 38]+op[ 39]*1j,op[ 40]+op[ 41]*1j,op[ 42]+op[ 43]*1j,op[ 44]+op[ 45]*1j,op[ 46]+op[ 47]*1j,op[ 48]+op[ 49]*1j,op[ 50]+op[ 51]*1j],
                                    [op[ 52]+op[ 53]*1j,op[ 54]+op[ 55]*1j,op[ 56]+op[ 57]*1j,op[ 58]+op[ 59]*1j,op[ 60]+op[ 61]*1j,op[ 62]+op[ 63]*1j,op[ 64]+op[ 65]*1j,op[ 66]+op[ 67]*1j],
                                    [op[ 68]+op[ 69]*1j,op[ 70]+op[ 71]*1j,op[ 72]+op[ 73]*1j,op[ 74]+op[ 75]*1j,op[ 76]+op[ 77]*1j,op[ 78]+op[ 79]*1j,op[ 80]+op[ 81]*1j,op[ 82]+op[ 83]*1j],
                                    [op[ 84]+op[ 85]*1j,op[ 86]+op[ 87]*1j,op[ 88]+op[ 89]*1j,op[ 90]+op[ 91]*1j,op[ 92]+op[ 93]*1j,op[ 94]+op[ 95]*1j,op[ 96]+op[ 97]*1j,op[ 98]+op[ 99]*1j],
                                    [op[100]+op[101]*1j,op[102]+op[103]*1j,op[104]+op[105]*1j,op[106]+op[107]*1j,op[108]+op[109]*1j,op[110]+op[111]*1j,op[112]+op[113]*1j,op[114]+op[115]*1j],
                                    [op[116]+op[117]*1j,op[118]+op[119]*1j,op[120]+op[121]*1j,op[122]+op[123]*1j,op[124]+op[125]*1j,op[126]+op[127]*1j,op[128]+op[129]*1j,op[130]+op[131]*1j]])
                circ.append(gate,[int(op[3]), int(op[2]), int(op[1])])
        # print(circ)
        return circ

    def verify_IO(self, file_state, qiskit_state):
        N = self.args.total_qbit
        NUMFILE = 1 << self.args.file_qbit
        FILESIZE = int((1 << N)/NUMFILE)
        flag = True
        for i in range(NUMFILE):
            if (not np.alltrue(np.abs(qiskit_state[i*FILESIZE:(i+1)*FILESIZE]-file_state[i]) < 1e-9)):
                flag = False
                break
        return flag

    def verify_MEM(self, file_state, qiskit_state):
        N = self.args.total_qbit
        FILESIZE = 1 << N
        flag = True
        if (not np.alltrue(np.abs(qiskit_state[:FILESIZE]-file_state[0]) < 1e-9)):
            flag = False
        return flag

    def run_IO(self, name:str, quiet:bool = False):
        flag = True
        if(self.args.is_density):
            # [not done yet]
            # fd_state = read_state(state_path, 2*N, NGQB)
            # circ = set_circuit(circuit_path, N)
            # qiskit_state = qiskit_init_density_matrix(circ)
            # if(not check(fd_state, qiskit_state, 2*N, NGQB)):
            #     print("test not pass")
            #     print()
            #     flag = False
            pass
        else:
            file_state = self.loader_IO()
            circ = self.qiskit_set_circuit()
            qiskit_state = self.qiskit_init_state_vector(circ)
            if(not self.verify_IO(file_state, qiskit_state)):
                print("test not pass")
                print()
                flag = False
        if(quiet):
            if(not flag):
                print("[x]", name, "not pass under 1e-9", flush=True)
        else:
            if(flag):
                print("[pass]", name, ": match with qiskit under 1e-9", flush=True)
            else:
                print("[x]", name, "not pass under 1e-9", flush=True)
        return flag

    def run_MEM(self, name:str, quiet:bool = False):
        flag = True
        if(self.args.is_density):
            # [not done yet]
            # fd_state = read_state(state_path, 2*N, NGQB)
            # circ = set_circuit(circuit_path, N)
            # qiskit_state = qiskit_init_density_matrix(circ)
            # if(not check(fd_state, qiskit_state, 2*N, NGQB)):
            #     print("test not pass")
            #     print()
            #     flag = False
            pass
        else:
            file_state = self.loader_MEM()
            circ = self.qiskit_set_circuit()
            qiskit_state = self.qiskit_init_state_vector(circ)
            if(not self.verify_MEM(file_state, qiskit_state)):
                print("test not pass")
                print()
                flag = False
        if(quiet):
            if(not flag):
                print("[x]", name, "not pass under 1e-9", flush=True)
        else:
            if(flag):
                print("[pass]", name, ": match with qiskit under 1e-9", flush=True)
            else:
                print("[x]", name, "not pass under 1e-9", flush=True)
        return flag

def print_header(N, NGQB, NSQB, NLQB, isDensity):
    NUMFILE = 1 << NGQB
    NUMTD = 1 << NSQB
    CHUNKSIZE = 1 << NLQB

    if(isDensity):
        FILESIZE = 1 << (2*N-NGQB)
        print(f"N = {N}")
        print(f"NGQB = {NGQB}, #FILE = {NUMFILE}, FILESIZE = {FILESIZE}")
    else:
        FILESIZE = 1 << (N-NGQB)
        print(f"N = {N}")
        print(f"NGQB = {NGQB}, #FILE = {NUMFILE}, FILESIZE = {FILESIZE}")
        print(f"NSQB = {NSQB}, #Thread = {NUMTD}")
        print(f"NLQB = {NLQB}, CHUNKSIZE = {CHUNKSIZE}")
        print(f"global = {[i for i in range(0, NGQB)]}")
        print(f"thread = {[i for i in range(NGQB, NSQB)]}")
        print(f"middle = {[i for i in range(NSQB, N-NLQB)]}")
        print(f"local  = {[i for i in range(N-NLQB, N)]}")

    print("===========================", flush=True)

# [Qiskit]: init to 0 state
def circ_init(N:int):
    circ = QuantumCircuit(N)
    initial_state = [1,0]   # Define initial_state as |0>
    circ.initialize(initial_state, 0)
    return circ

# [Qiskit]: state_vector
def qiskit_init_state_vector(circ):
    simulator = Aer.get_backend('aer_simulator_statevector')
    circ.save_statevector(label=f'save')
    circ = transpile(circ, simulator)

    # start_time = time.perf_counter()
    data = simulator.run(circ).result().data(0)
    # end_time = time.perf_counter()
    # print(f"Time: {(end_time - start_time)*1000000:.0f} (us)\n")

    return data['save']

# [Qiskit]: density_matrix
def qiskit_init_density_matrix(circ):
    simulator = Aer.get_backend('aer_simulator_density_matrix')
    circ.save_density_matrix(label=f'save')
    circ = transpile(circ, simulator)

    # start_time = time.perf_counter()
    data = simulator.run(circ).result().data(0)
    # end_time = time.perf_counter()
    # print(f"Time: {(end_time - start_time)*1000000:.0f} (us)\n")

    return data['save'].T.reshape(-1)

def loader_IO(state_paths, N, NUMFILE):
    FILESIZE = int((1 << N) / NUMFILE)
    fd_arr = [np.zeros(FILESIZE, dtype=np.complex128) for i in range(NUMFILE)]
    for fd, state_path in enumerate(state_paths.split(',')):
        f = fd_arr[fd]
        with open(state_path, mode="rb") as state_file:
            try:
                state = state_file.read()
                k = 0
                for i in range(FILESIZE):
                    (real, imag) = struct.unpack("dd", state[k:k+16])
                    f[i] = real+imag*1j
                    k += 16
            except Exception as e:
                print(e)
                print(f"read from {state_path}")
                print(f"[ERROR]: error at reading {k}th byte")
                exit()
    return fd_arr

def loader_MEM(path, N):
    FILESIZE = 1 << N
    sv = np.zeros(FILESIZE, dtype=np.complex128)

    with open(path, mode="r") as sv_file:
        for i in range(FILESIZE):
            try:
                real = float(sv_file.readline())
                imag = float(sv_file.readline())
                sv[i] = real + imag*1j
            except Exception as e:
                print(e)
                print(f"read from {path}")
                print(f"[ERROR]: error at reading {i}th byte")
                exit()

    return sv

def set_circuit(path, N):
    circ = circ_init(N)
    ops = []
    with open(path,'r') as f:
        for line in f.readlines():
            s = line.split()
            ops.append(s)

    # for op in ops:
    #     print(op)
    for op in ops:
        if op[0]=="H":
            circ.h(int(op[1]))
        if op[0]=="S":
            circ.s(int(op[1]))
        if op[0]=="T":
            circ.t(int(op[1]))
        if op[0]=="X":
            circ.x(int(op[1]))
        if op[0]=="Y":
            circ.y(int(op[1]))
        if op[0]=="Z":
            circ.z(int(op[1]))
        if op[0]=="P":
            circ.p(float(op[2]), int(op[1]))
        if op[0]=="U1": # UnitaryGate
            for i in range(2, 10):
                op[i] = float(op[i])
            gate = UnitaryGate([[op[2]+op[3]*1j, op[4]+op[5]*1j], 
                                [op[6]+op[7]*1j, op[8]+op[9]*1j]])
            circ.append(gate,[int(op[1])])

        if op[0]=="CX":
            circ.cx(int(op[1]), int(op[2]))
        if op[0]=="CY":
            circ.cy(int(op[1]), int(op[2]))
        if op[0]=="CZ":
            circ.cy(int(op[1]), int(op[2]))
        if op[0]=="CP":
            circ.cp(float(op[3]), int(op[1]), int(op[2]))
        if op[0]=="CU1":
            # control-unitary
            for i in range(3, 11):
                op[i] = float(op[i])
            gate = UnitaryGate([[1, 0, 0, 0],
                                [0, 1, 0, 0],
                                [0, 0, op[3]+op[4]*1j, op[5]+op[ 6]*1j],
                                [0, 0, op[7]+op[8]*1j, op[9]+op[10]*1j]])
            circ.append(gate,[int(op[1]), int(op[2])])
        if op[0]=="SWAP":
            circ.swap(int(op[1]), int(op[2]))
        if op[0]=="TOFFOLI":
            # control1 control2 targert
            circ.toffoli(int(op[1]), int(op[2]), int(op[3]))

        if op[0]=="U2": # 2 qubit UnitaryGate
            for i in range(3, 35):
                op[i] = float(op[i])
            gate = UnitaryGate([[op[ 3]+op[ 4]*1j,op[ 5]+op[ 6]*1j,op[ 7]+op[ 8]*1j,op[ 9]+op[10]*1j],
                                [op[11]+op[12]*1j,op[13]+op[14]*1j,op[15]+op[16]*1j,op[17]+op[18]*1j],
                                [op[19]+op[20]*1j,op[21]+op[22]*1j,op[23]+op[24]*1j,op[25]+op[26]*1j],
                                [op[27]+op[28]*1j,op[29]+op[30]*1j,op[31]+op[32]*1j,op[33]+op[34]*1j]])
            # Qiskit treat the order in a reverse way.
            # WHAAAT?!!!WHYYYYYY?????
            # circ.append(gate,[reorder(op[4], N), reorder(op[5], N)])
            circ.append(gate,[int(op[2]), int(op[1])])

        if op[0]=="U3":  # 3 qubit UnitaryGate
            for i in range(4, 132):
                op[i] = float(op[i])
            gate = UnitaryGate([[op[  4]+op[  5]*1j,op[  6]+op[  7]*1j,op[  8]+op[  9]*1j,op[ 10]+op[ 11]*1j,op[ 12]+op[ 13]*1j,op[ 14]+op[ 15]*1j,op[ 16]+op[ 17]*1j,op[ 18]+op[ 19]*1j],
                                [op[ 20]+op[ 21]*1j,op[ 22]+op[ 23]*1j,op[ 24]+op[ 25]*1j,op[ 26]+op[ 27]*1j,op[ 28]+op[ 29]*1j,op[ 30]+op[ 31]*1j,op[ 32]+op[ 33]*1j,op[ 34]+op[ 35]*1j],
                                [op[ 36]+op[ 37]*1j,op[ 38]+op[ 39]*1j,op[ 40]+op[ 41]*1j,op[ 42]+op[ 43]*1j,op[ 44]+op[ 45]*1j,op[ 46]+op[ 47]*1j,op[ 48]+op[ 49]*1j,op[ 50]+op[ 51]*1j],
                                [op[ 52]+op[ 53]*1j,op[ 54]+op[ 55]*1j,op[ 56]+op[ 57]*1j,op[ 58]+op[ 59]*1j,op[ 60]+op[ 61]*1j,op[ 62]+op[ 63]*1j,op[ 64]+op[ 65]*1j,op[ 66]+op[ 67]*1j],
                                [op[ 68]+op[ 69]*1j,op[ 70]+op[ 71]*1j,op[ 72]+op[ 73]*1j,op[ 74]+op[ 75]*1j,op[ 76]+op[ 77]*1j,op[ 78]+op[ 79]*1j,op[ 80]+op[ 81]*1j,op[ 82]+op[ 83]*1j],
                                [op[ 84]+op[ 85]*1j,op[ 86]+op[ 87]*1j,op[ 88]+op[ 89]*1j,op[ 90]+op[ 91]*1j,op[ 92]+op[ 93]*1j,op[ 94]+op[ 95]*1j,op[ 96]+op[ 97]*1j,op[ 98]+op[ 99]*1j],
                                [op[100]+op[101]*1j,op[102]+op[103]*1j,op[104]+op[105]*1j,op[106]+op[107]*1j,op[108]+op[109]*1j,op[110]+op[111]*1j,op[112]+op[113]*1j,op[114]+op[115]*1j],
                                [op[116]+op[117]*1j,op[118]+op[119]*1j,op[120]+op[121]*1j,op[122]+op[123]*1j,op[124]+op[125]*1j,op[126]+op[127]*1j,op[128]+op[129]*1j,op[130]+op[131]*1j]])
            circ.append(gate,[int(op[3]), int(op[2]), int(op[1])])
    # print(circ)
    return circ

def check(fd_state, qiskit_state, N, NUMFILE):
    FILESIZE = int((1<< N) / NUMFILE)
    qiskit_state = np.array(qiskit_state)
    flag = True
    for i in range(NUMFILE):
        if (not np.alltrue(np.abs(qiskit_state[i*FILESIZE:(i+1)*FILESIZE]-fd_state[i]) < 1e-9)):
            flag = False
            break
    return flag


def check_mem(sv_state, qiskit_state):
    qiskit_state = np.array(qiskit_state)
    tol = 1e-9
    return np.alltrue(np.abs(qiskit_state-sv_state) < tol)

def simple_test(name:str, isDensity:bool, circuit_path, state_path, N, NUMFILE, quiet:bool = False):
    flag = True
    if(isDensity):
        # [not done yet]
        # fd_state = read_state(state_path, 2*N, NGQB)
        # circ = set_circuit(circuit_path, N)
        # qiskit_state = qiskit_init_density_matrix(circ)
        # if(not check(fd_state, qiskit_state, 2*N, NGQB)):
        #     print("test not pass")
        #     print()
        #     flag = False
        pass

    else:
        fd_state = loader_IO(state_path, N, NUMFILE)
        circ = set_circuit(circuit_path, N)
        qiskit_state = qiskit_init_state_vector(circ)
        if(not check(fd_state, qiskit_state, N, NUMFILE)):
            print("test not pass")
            print()
            flag = False
    if(quiet):
        if(not flag):
            print("[x]", name, "not pass under 1e-9", flush=True)
    else:
        if(flag):
            print("[pass]", name, ": match with qiskit under 1e-9", flush=True)
        else:
            print("[x]", name, "not pass under 1e-9", flush=True)

    return flag

def simple_test_mem(name:str, isDensity:bool, circuit_path, state_path, N, quiet:bool = False):
    flag = True
    if(isDensity):
        # [not done yet]
        # fd_state = read_state(state_path, 2*N, NGQB)
        # circ = set_circuit(circuit_path, N)
        # qiskit_state = qiskit_init_density_matrix(circ)
        # if(not check(fd_state, qiskit_state, 2*N, NGQB)):
        #     print("test not pass")
        #     print()
        #     flag = False
        pass

    else:
        fd_state = loader_MEM(state_path, N)
        circ = set_circuit(circuit_path, N)
        qiskit_state = qiskit_init_state_vector(circ)
        if(not check_mem(fd_state, qiskit_state)):
            print("test not pass")
            print()
            flag = False
    if(quiet):
        if(not flag):
            print("[x]", name, "not pass under 1e-9", flush=True)
    else:
        if(flag):
            print("[pass]", name, ": match with qiskit under 1e-9", flush=True)
        else:
            print("[x]", name, "not pass under 1e-9", flush=True)

    return flag