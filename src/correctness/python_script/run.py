import os

qubits  = list(range(34, 39+1))
#qubits  = list(range(21, 31+1))
options = ['ori', 'sub', 'ideal']
#options = ['ideal']
#cirs = ["qft"]
#cirs = ["qaoa"]

#dirs = list(range(21, 30+1))
#options = ['sub']

'''
qubits  = list(range(37, 39+1))
cirs = ["h"]
for c in cirs:
    for q in qubits:
        for o in options:
            dir_path = f"../circuit_io/{o}/{q}"
            #dir_path = f"../circuit_io/{o}/{q}"
            print(f'../../Quokka -i {dir_path}/cpu.ini -c {dir_path}/{c}{q}.txt > {dir_path}/{c}{q}.log')
            os.system(f'../../Quokka -i {dir_path}/cpu.ini -c {dir_path}/{c}{q}.txt > {dir_path}/{c}{q}.log')


qubits  = list(range(34, 39+1))
cirs = ["rx", "rzz"]
for c in cirs:
    for q in qubits:
        for o in options:
            dir_path = f"../circuit_io/{o}/{q}"
            #dir_path = f"../circuit_io/{o}/{q}"
            print(f'../../Quokka -i {dir_path}/cpu.ini -c {dir_path}/{c}{q}.txt > {dir_path}/{c}{q}.log')
            os.system(f'../../Quokka -i {dir_path}/cpu.ini -c {dir_path}/{c}{q}.txt > {dir_path}/{c}{q}.log')

'''
'''
qubits  = list(range(32, 34))
cirs = ["qft", "qaoa"]
for c in cirs:
    for q in qubits:
        for o in options:
            dir_path = f"../circuit_io/{o}/{q}"
            #dir_path = f"../circuit_io/{o}/{q}"
            print(f'../../Quokka -i {dir_path}/cpu.ini -c {dir_path}/{c}{q}.txt > {dir_path}/{c}{q}.log')
            os.system(f'../../Quokka -i {dir_path}/cpu.ini -c {dir_path}/{c}{q}.txt > {dir_path}/{c}{q}.log')
'''

qubits = [36, 39]
#cirs = ["qft", "qaoa"]
cirs = ["qaoa"]
for c in cirs:
    for q in qubits:
        for o in options:
            dir_path = f"../circuit_io/{o}/{q}"
            #dir_path = f"../circuit_io/{o}/{q}"
            print(f'../../Quokka -i {dir_path}/cpu.ini -c {dir_path}/{c}{q}.txt > {dir_path}/{c}{q}.log')
            os.system(f'../../Quokka -i {dir_path}/cpu.ini -c {dir_path}/{c}{q}.txt > {dir_path}/{c}{q}.log')

