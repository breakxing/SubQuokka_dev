import math

def get_circuit():
    circuit=[]
    return circuit

def H(circuit, qbit):
    circuit.append(f'H {qbit}')
def S(circuit, qbit):
    circuit.append(f'S {qbit}')
def T(circuit, qbit):
    circuit.append(f'T {qbit}')
def X(circuit, qbit):
    circuit.append(f'X {qbit}')
def Y(circuit, qbit):
    circuit.append(f'Y {qbit}')
def Z(circuit, qbit):
    circuit.append(f'Z {qbit}')
def Phase(circuit, qbit, angle):
    circuit.append(f'P {qbit} {angle}')
def U1(circuit, qbit, coeff:list):    
    circuit.append(f"U1 {qbit} {' '.join(map(str, coeff))}")
def RX(circuit, qbit, angle):
    circuit.append(f'RX {qbit} {angle}')

def CX(circuit, control, target):
    circuit.append(f'CX  {control} {target}')
def CY(circuit, control, target):
    circuit.append(f'CY  {control} {target}')
def CZ(circuit, control, target):
    circuit.append(f'CZ  {control} {target}')
def CPhase(circuit, control, target, angle):
    circuit.append(f'CP  {control} {target} {angle}')
def CU1(circuit, control, target, coeff:list):
    circuit.append(f"CU1 {control} {target} {' '.join(map(str, coeff))}")
def RZZ(circuit, target0, target1, angle):
    circuit.append(f'RZZ {target0} {target1} {angle}')


def SWAP(circuit, target1, target2):
    circuit.append(f'SWAP {target1} {target2}')

def measure(circuit, targ, shots, set_src = 0, set_dst = 1):
    circuit.append(f"20 3 1 2 {set_src} {set_dst} {shots} {targ} 0 0 0 0")
def measure_multi(circuit, targ:list, shots, set_src = 0, set_dst = 1):
    circuit.append(f"21 3 0 {len(targ)} {set_src} {set_dst} {shots} {' '.join(map(str, targ))} {' '.join(map(str, targ))}")

def U2(circuit, q0, q1, coeff:list):
    circuit.append(f"U2  {q0} {q1} {' '.join(map(str, coeff))}")

def CCX_true(circuit, ctrl0, ctrl1, targ):
    circuit.append(f"14 2 1 0  {ctrl0} {ctrl1} {targ}")

def CCX(circuit, ctrl0, ctrl1, targ):
    if ctrl0 == ctrl1 or ctrl0 == targ or ctrl1 == targ:
        print(f"c0:{ctrl0}, c1:{ctrl1}, t:{targ} must not be the same qubit.")
        exit()
    imag = []
    for _ in range(64):
        imag.append(0)

    real = []
    if ctrl0 < targ and ctrl1 < targ:
        if ctrl0 < ctrl1:
            q0 = ctrl0
            q1 = ctrl1
        else:
            q0 = ctrl1
            q1 = ctrl0
        q2 = targ
        real = [1, 0, 0, 0, 0, 0, 0, 0,
                0, 1, 0, 0, 0, 0, 0, 0,
                0, 0, 1, 0, 0, 0, 0, 0,
                0, 0, 0, 1, 0, 0, 0, 0,
                0, 0, 0, 0, 1, 0, 0, 0,
                0, 0, 0, 0, 0, 1, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 1,
                0, 0, 0, 0, 0, 0, 1, 0]
        # for i in range(8):
        #     if i == 6:
        #         for j in range(8):
        #             if j == 7:  real.append(1)
        #             else:       real.append(0)
        #     elif i == 7:            
        #         for j in range(8):
        #             if j == 6:  real.append(1)
        #             else:       real.append(0)
        #     else:
        #         for j in range(8):
        #             if i == j:  real.append(1)
        #             else:       real.append(0)
    elif (ctrl0 < targ and targ < ctrl1) or (ctrl1 < targ and targ < ctrl0):
        if ctrl0 < targ and targ < ctrl1:
            q0 = ctrl0
            q1 = targ
            q2 = ctrl1
        else:
            q0 = ctrl1
            q1 = targ
            q2 = ctrl0
        real = [1, 0, 0, 0, 0, 0, 0, 0,
                0, 1, 0, 0, 0, 0, 0, 0,
                0, 0, 1, 0, 0, 0, 0, 0,
                0, 0, 0, 1, 0, 0, 0, 0,
                0, 0, 0, 0, 1, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 1,
                0, 0, 0, 0, 0, 0, 1, 0,
                0, 0, 0, 0, 0, 1, 0, 0]
        # for i in range(8):
        #     if i == 5:
        #         for j in range(8):
        #             if j == 7:  real.append(1)
        #             else:       real.append(0)
        #     elif i == 7:            
        #         for j in range(8):
        #             if j == 5:  real.append(1)
        #             else:       real.append(0)
        #     else:
        #         for j in range(8):
        #             if i == j:  real.append(1)
        #             else:       real.append(0)

    else:
        q0 = targ
        if ctrl0 < ctrl1:
            q1 = ctrl0
            q2 = ctrl1
        else:
            q1 = ctrl1
            q2 = ctrl0
        real = [1, 0, 0, 0, 0, 0, 0, 0,
                0, 1, 0, 0, 0, 0, 0, 0,
                0, 0, 1, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 1,
                0, 0, 0, 0, 1, 0, 0, 0,
                0, 0, 0, 0, 0, 1, 0, 0,
                0, 0, 0, 0, 0, 0, 1, 0,
                0, 0, 0, 1, 0, 0, 0, 0]
        # for i in range(8):
        #     if i == 5:
        #         for j in range(8):
        #             if j == 7:  real.append(1)
        #             else:       real.append(0)
        #     elif i == 7:            
        #         for j in range(8):
        #             if j == 5:  real.append(1)
        #             else:       real.append(0)
        #     else:
        #         for j in range(8):
        #             if i == j:  real.append(1)
        #             else:       real.append(0)
    circuit.append(f"32 0 3 64  {q0} {q1} {q2} {' '.join(map(str, real))} {' '.join(map(str, imag))}")

def U3(circuit, q0, q1, q2, coeff:list):
    circuit.append(f"U3  {q0} {q1} {q2} {' '.join(map(str, coeff))}")

def COPY(circuit, set_src, set_dst):
    circuit.append(f"22 2 0 0 {set_src} {set_dst}")

def inverse_qft_phase_on_work_section(circuit,work_section):
    H(circuit,work_section[0])
    for id in range(1,len(work_section)):
        CPhase(circuit,work_section[id],work_section[0],math.pi/(2**id))

def qft_phase_on_work_section(circuit, work_section):
    for id in range(len(work_section)-1, 0, -1):
        CPhase(circuit, work_section[id], work_section[0], math.pi / (2**id))
    H(circuit,work_section[0])


# def qft(circuit, start_qubit, end_qubit):
#     total_bit = end_qubit - start_qubit
#     for swap_idx in range(int(total_bit / 2)):
#         SWAP(circuit, start_qubit + swap_idx, end_qubit - 1 - swap_idx)
#     for i,times in enumerate(range(end_qubit - start_qubit)):
#         work_section=list(range(end_qubit - 1 - i, end_qubit))
#         qft_phase_on_work_section(circuit, work_section)

def qft(circuit, start_qubit, total_bit):
    end_qubit = start_qubit + total_bit
    for swap_idx in range(int(total_bit / 2)):
        SWAP(circuit, start_qubit + swap_idx, end_qubit - 1 - swap_idx)
    for i,times in enumerate(range(total_bit)):
        work_section = list(range(end_qubit - 1 - i, end_qubit))
        qft_phase_on_work_section(circuit, work_section)

def pseudo_qft_phase_on_work_section(circuit, work_section):
    for id in range(len(work_section)-1, 0, -1):
        q0 = work_section[id]
        q1 = work_section[0]
        if q0 > q1:
            t = q1
            q1 = q0
            q0 = t
        U2(circuit, q0, q1, [0.5, 0,   0, 0,   0, 0,   0, 0,
                               0, 0, 0.5, 0,   0, 0,   0, 0,
                               0, 0,   0, 0, 0.5, 0,   0, 0,
                               0, 0,   0, 0,   0, 0, 0.5, 0])
    H(circuit, work_section[0])

def pseudo_qft(circuit, start_qubit, end_qubit):
    total_bit = end_qubit - start_qubit+1
    for swap_idx in range(int(total_bit/2)):
        SWAP(circuit, start_qubit+swap_idx, end_qubit-swap_idx)
    for i,times in enumerate(range(end_qubit - start_qubit+1)):
        work_section = list(range(end_qubit-i, end_qubit+1))
        pseudo_qft_phase_on_work_section(circuit, work_section)

def inverse_qft(circuit,start_qubit,end_qubit):
    for i,times in enumerate(range(end_qubit-start_qubit+1)):
        work_section=list(range(start_qubit+i,end_qubit+1))
        inverse_qft_phase_on_work_section(circuit,work_section)
    total_bit=end_qubit-start_qubit+1
    for swap_idx in range(int(total_bit/2)):
        SWAP(circuit,start_qubit+swap_idx,end_qubit-swap_idx)

def create_circuit(circuit, path):
    f=open(path,'w')
    for data in circuit:
        print(data,file=f)
    f.close()

def create_subcircuit(circuit, path):
    f=open(path,'w')
    for subcircuit in circuit:
        print(len(subcircuit), file=f)
        for data in subcircuit:
            print(data, file=f)
    f.close()
