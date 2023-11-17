import os


qubits = 30
os.system(f"g++ ../../../subcircuitFinder/finder.cpp -o ../../../subcircuitFinder/test")
files = [f"h{qubits}.txt",f"qaoa{qubits}.txt",f"qft{qubits}.txt"]

for file in files:
    os.system(f"../../../subcircuitFinder/test ori/{qubits}/{file} 0 ../circuit_io/sub/{qubits}/{file}")