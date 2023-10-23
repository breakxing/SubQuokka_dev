from circuit_generator import *
from ini_generator import *
from test_util import verifier
import random
import os
import numpy as np
from scipy.stats import unitary_group

# Test for gate set on SubQuokka

gate_list = ["H U1 X Y Z P RX RY RZ", "RZZ, SWAP, CPhase, VSWAP_2_2"]

# One qubit gate
# N = 3
# H U1 X Y Z P RX RY RZ
n = 3
ini_path = "res.ini"
cir_path = "cir.txt"

args = Args(total_qbit=n, file_qbit=1, chunk_qbit=1, 
                is_subcircuit=0,
                runner_type="MEM", dump_file="res.txt")
ini = Ini(args, ini_path)
ini.out()

# H
for j in range(n):
    cir = get_circuit()
    for i in range(n):
        H(cir, i)
    H(cir, j)
    create_circuit(cir, cir_path)
    
    os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null 2>&1")
    # os.system(f"../Quokka -c {cir_path} -i {ini_path}")
    v = verifier(args, cir_path)
    flag = v.run(f"H gate on qubit {j}")
    if(not flag):
        exit()

print("[PASS]: H gate done.")
print("==================================================")

# U1
for j in range(n):
    random_unitary = unitary_group.rvs(2)
    coeff = []
    for x in random_unitary:
        for y in x:
            coeff.append(y.real)
            coeff.append(y.imag)
    cir = get_circuit()
    for i in range(n):
        H(cir, i)
    U1(cir, j, coeff)
    create_circuit(cir, cir_path)
    
    os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null 2>&1")
    # os.system(f"../Quokka -c {cir_path} -i {ini_path}")
    v = verifier(args, cir_path)
    flag = v.run(f"U1 gate on qubit {j}")
    if(not flag):
        exit()

print("[PASS]: U1 gate done.")
print("==================================================")

# X
for j in range(n):
    cir = get_circuit()
    for i in range(n):
        if i == j:
            X(cir, i)
        else:
            H(cir, i)
    create_circuit(cir, cir_path)
    
    os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null 2>&1")
    # os.system(f"../Quokka -c {cir_path} -i {ini_path}")
    v = verifier(args, cir_path)
    flag = v.run(f"X gate on qubit {j}")
    if(not flag):
        exit()

print("[PASS]: X gate done.")
print("==================================================")

# Y
for j in range(n):
    cir = get_circuit()
    for i in range(n):
        if i == j:
            Y(cir, i)
        else:
            H(cir, i)
    create_circuit(cir, cir_path)
    
    os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null 2>&1")
    # os.system(f"../Quokka -c {cir_path} -i {ini_path}")
    v = verifier(args, cir_path)
    flag = v.run(f"Y gate on qubit {j}")
    if(not flag):
        exit()

print("[PASS]: Y gate done.")
print("==================================================")

# Z
for j in range(n):
    cir = get_circuit()
    for i in range(n):
        if i == j:
            X(cir, i)
            Z(cir, i)
        else:
            H(cir, i)
    create_circuit(cir, cir_path)
    
    os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null 2>&1")
    # os.system(f"../Quokka -c {cir_path} -i {ini_path}")
    v = verifier(args, cir_path)
    flag = v.run(f"Z gate on qubit {j}")
    if(not flag):
        exit()

print("[PASS]: Z gate done.")
print("==================================================")

# P
for j in range(n):
    cir = get_circuit()
    for i in range(n):
        if i == j:
            X(cir, i)
            Phase(cir, i, 2*math.pi * random.random())
        else:
            H(cir, i)
    create_circuit(cir, cir_path)
    
    os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null 2>&1")
    # os.system(f"../Quokka -c {cir_path} -i {ini_path}")
    v = verifier(args, cir_path)
    flag = v.run(f"P gate on qubit {j}")
    if(not flag):
        exit()

print("[PASS]: P gate done.")
print("==================================================")

# RX
for j in range(n):
    cir = get_circuit()
    for i in range(n):
        if i == j:
            X(cir, i)
            RX(cir, i, 2*math.pi * random.random())
        else:
            H(cir, i)
    create_circuit(cir, cir_path)
    
    os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null 2>&1")
    # os.system(f"../Quokka -c {cir_path} -i {ini_path}")
    v = verifier(args, cir_path)
    flag = v.run(f"RX gate on qubit {j}")
    if(not flag):
        exit()

print("[PASS]: RX gate done.")
print("==================================================")

# RY
for j in range(n):
    cir = get_circuit()
    for i in range(n):
        if i == j:
            X(cir, i)
            RY(cir, i, 2*math.pi * random.random())
        else:
            H(cir, i)
    create_circuit(cir, cir_path)
    
    os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null 2>&1")
    # os.system(f"../Quokka -c {cir_path} -i {ini_path}")
    v = verifier(args, cir_path)
    flag = v.run(f"RY gate on qubit {j}")
    if(not flag):
        exit()

print("[PASS]: RY gate done.")
print("==================================================")

# RZ
for j in range(n):
    cir = get_circuit()
    for i in range(n):
        if i == j:
            X(cir, i)
            RZ(cir, i, 2*math.pi * random.random())
        else:
            H(cir, i)
    create_circuit(cir, cir_path)
    
    os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null 2>&1")
    # os.system(f"../Quokka -c {cir_path} -i {ini_path}")
    v = verifier(args, cir_path)
    flag = v.run(f"RZ gate on qubit {j}")
    if(not flag):
        exit()

print("[PASS]: RZ gate done.")
print("==================================================")


print("[PASS]: ONE QUBIT GATE DONE.")
print("==================================================")



# Two qubit gate
# [重要]: qubit傳入的順序 除非是control gate 不然一定要由小排到大
# N = 6
# RZZ SWAP CPhase

n = 6
args = Args(total_qbit=n, file_qbit=2, chunk_qbit=2, 
                is_subcircuit=0,
                runner_type="MEM", dump_file="res.txt")
ini = Ini(args, ini_path)
ini.out()


# RZZ
for k in range(n):
    for j in range(n):
        if k == j:
            continue
        cir = get_circuit()
        for i in range(n):
            H(cir, i)
        RZZ(cir, j, k, 2*math.pi * random.random())
        create_circuit(cir, cir_path)
        
        os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null 2>&1")
        # os.system(f"../Quokka -c {cir_path} -i {ini_path}")
        v = verifier(args, cir_path)
        flag = v.run(f"RZZ gate on qubit {k} {j}")
        if(not flag):
            exit()

print("[PASS]: RZZ gate done.")
print("==================================================")

# SWAP
for k in range(n):
    for j in range(n):
        if k == j:
            continue
        cir = get_circuit()
        for i in range(n):
            if i == j:
                X(cir, i)
            elif i == k:
                Y(cir, i)
            else:
                H(cir, i)
        SWAP(cir, k, j)
        create_circuit(cir, cir_path)
        
        os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null 2>&1")
        # os.system(f"../Quokka -c {cir_path} -i {ini_path}")
        v = verifier(args, cir_path)
        flag = v.run(f"SWAP gate on qubit {k} {j}")
        if(not flag):
            exit()

print("[PASS]: SWAP gate done.")
print("==================================================")

# CPhase
for k in range(n):
    for j in range(n):
        if k == j:
            continue
        cir = get_circuit()
        for i in range(n):
            if i == j or i == k:
                pass
            else:
                H(cir, i)
        CPhase(cir, k, j, 2*math.pi*random.random())
        create_circuit(cir, cir_path)

        os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null 2>&1")
        # os.system(f"../Quokka -c {cir_path} -i {ini_path}")
        v = verifier(args, cir_path)
        flag = v.run(f"CPhase gate on qubit {k} {j}")
        if(not flag):
            exit()

print("[PASS]: CPhase gate done.")
print("==================================================")

print("[PASS]: TWO QUBIT GATE DONE.")
print("==================================================")


n = 7
args = Args(total_qbit=n, file_qbit=2, chunk_qbit=3, 
                is_subcircuit=0,
                runner_type="MEM", dump_file="res.txt")
ini = Ini(args, ini_path)
ini.out()

# VSWAP_2_2
set0 = [[0, 1], [0, 2], [1, 2]]
set1 = [[3, 4], [3, 5], [3, 6], [4, 5], [4, 6], [5, 6]]

for s0 in set0:
    for s1 in set1:
        q0 = s0[0]
        q1 = s0[1]
        q2 = s1[0]
        q3 = s1[1]
        
        cir = get_circuit()
        for i in range(n):
            if i == q0 or i == q1 or i == q2 or i == q3:
                pass
            else:
                H(cir, i)
        X(cir, q0)
        X(cir, q3)
        cir_file = cir.copy()
        cir_verify = cir.copy()
        VSWAP_2_2(cir_file, q0, q1, q2, q3)
        create_circuit(cir_file, cir_path)
        
        os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null 2>&1")
        # os.system(f"../Quokka -c {cir_path} -i {ini_path}")
        SWAP(cir_verify, q0, q2)
        SWAP(cir_verify, q1, q3)
        create_circuit(cir_verify, cir_path)
        v = verifier(args, cir_path)
        flag = v.run(f"VSWAP_2_2 gate on qubit {q0} {q1} {q2} {q3}")
        if(not flag):
            exit()

print("[PASS]: VSWAP_2_2 gate done.")
print("==================================================")




print(f"[PASS]: {gate_list}")
print("==================================================")

exit()





state_paths = "./state/path0,./state/path1,./state/path2,./state/path3,./state/path4,./state/path5,./state/path6,./state/path7,./state/path8,./state/path9,./state/path10,./state/path11,./state/path12,./state/path13,./state/path14,./state/path15,./state/path16,./state/path17,./state/path18,./state/path19,./state/path20,./state/path21,./state/path22,./state/path23,./state/path24,./state/path25,./state/path26,./state/path27,./state/path28,./state/path29,./state/path30,./state/path31,./state/path32,./state/path33,./state/path34,./state/path35,./state/path36,./state/path37,./state/path38,./state/path39,./state/path40,./state/path41,./state/path42,./state/path43,./state/path44,./state/path45,./state/path46,./state/path47,./state/path48,./state/path49,./state/path50,./state/path51,./state/path52,./state/path53,./state/path54,./state/path55,./state/path56,./state/path57,./state/path58,./state/path59,./state/path60,./state/path61,./state/path62,./state/path63"
numFILES = 6
chunk_seg = 12


cir = get_circuit()
for i in range(10):
    H(cir, i)
create_circuit(cir, cir_path)

subcir_path = "subcir.txt"
cir = get_circuit()
subcir = get_circuit()
for _ in range(10):
    H(subcir, 0)
cir.append(subcir)
create_subcircuit(cir, subcir_path)


# # MEM part
# # -circuit
# for n in range(33, 34):
# # for n in range(21, 34):
# # for n in range(33, 34):
#     args = Args(total_qbit=n, file_qbit=6, chunk_qbit=12, 
#                 is_subcircuit=0, dump_file="")
#     ini = Ini(args, ini_path)
#     ini.out()
#     # os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null")
#     os.system(f"../Quokka -c {cir_path} -i {ini_path}")
#     # simple_test_mem(f"[H {n} origin]", False, cir_path, "res.txt", N)

# exit()

# -subcircuit
# for n in range(21, 22):
# # for n in range(21, 34):
# # for n in range(33, 34):
#     args = Args(total_qbit=n, file_qbit=6, chunk_qbit=12, 
#                 is_subcircuit=1, dump_file="")
#     ini = Ini(args, ini_path)
#     ini.out()

#     # os.system(f"../Quokka -c {subcir_path} -i {ini_path} > /dev/null")
#     os.system(f"../Quokka -c {subcir_path} -i {ini_path}")
#     # simple_test_mem(f"[H {n} subcircuit]", False, cir_path, "res.txt", N

# exit()

# IO part
# -circuit
# for n in range(21, 34):
# for n in range(34, 36):
# for n in [32]:
for n in [30]:
    args = Args(total_qbit=n, file_qbit=6, chunk_qbit=12, 
                is_subcircuit=0,
                runner_type="IO", state_paths=state_paths)
    ini = Ini(args, ini_path)
    ini.out()

    # os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null")
    os.system(f"../Quokka_orig -c {cir_path} -i {ini_path}")
    os.system(f"../Quokka -c {cir_path} -i {ini_path}")
    # simple_test(f"[H {n} subcircuit]", False, cir_path, state_paths, N, numFILES)


# -subcircuit
# # for n in range(21, 34):
# # for n in range(34, 36):
# for n in [32]:
for n in [30]:
    args = Args(total_qbit=n, file_qbit=6, chunk_qbit=12, 
                is_subcircuit=1,
                runner_type="IO", state_paths=state_paths)
    ini = Ini(args, ini_path)
    ini.out()

    # os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null")
    os.system(f"../Quokka_orig -c {subcir_path} -i {ini_path}")
    os.system(f"../Quokka -c {subcir_path} -i {ini_path}")
    # simple_test(f"[H {n} subcircuit]", False, cir_path, state_paths, N, numFILES)

# exit()