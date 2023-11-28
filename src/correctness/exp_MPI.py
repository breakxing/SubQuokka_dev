from circuit_generator import *
from ini_generator import *
from test_util import *
import os

RDMA_CARD1 = 'mlx5_1:1'
RDMA_CARD2 = 'mlx5_0:1'
RDMA_TCP_CARD1 = 'enp225s0f1np1'
RDMA_TCP_CARD2 = 'enp225s0f0np0'
MORMAL_TCP_CARD1 = 'eno1'
MORMAL_TCP_CARD2 = 'enp2s0'


os.chdir("..")
os.system("make clean")
os.system("make -j")
os.system("scp -r -P 9048 ./*.cpp rdma2:~/SubQuokka_dev/src/")
os.system("scp -r -P 9048 ./*.hpp rdma2:~/SubQuokka_dev/src/")
os.system("scp -r -P 9048 ./*.h rdma2:~/SubQuokka_dev/src/")
input("Go to other computer to make")
os.chdir("correctness")
os.system(f"scp -r -P 9048 ./MPI_circuit/ rdma2:~/SubQuokka_dev/src/correctness/")
file_seg = 5
chunk_seg = 12


protocols = ["IO"]
versions = ["ori","sub"]
qubits = ["21","24","27","30"]
files = ["h","qft","qaoa"]
for protocol in protocols:
    for version in versions:
        for qubit in qubits:
            N = int(qubit)
            ini_path = os.path.join("MPI_circuit",version,qubit,"myini.ini")
            state_paths = [f"/mnt/state{i % 4}/path{i}" for i in range(1 << file_seg)]
            state_paths = ",".join(state_paths)
            args = Args(total_qbit=N, file_qbit=file_seg, chunk_qbit=chunk_seg,mpi_qbit=0,
                is_subcircuit = 1 if version == "sub" else 0,
                runner_type="IO", state_paths=state_paths)
            ini = Ini(args, ini_path)
            ini.out()
            for file in files:
                print(protocol,version,file + qubit,"Start")
                cir_path = os.path.join("circuit_io",version,qubit,file + qubit + ".txt")
                log_path = os.path.join("MPI_circuit",version,qubit,file + qubit + f".{protocol}.log")
                os.system(f"../Quokka -c {cir_path} -i {ini_path} > {log_path}")
            print(f"{protocol}",version,file + qubit,"End")


protocols = ["Pure_TCP"]
versions = ["ori","sub"]
qubits = ["21","24","27","30"]
files = ["h","qft","qaoa"]
for protocol in protocols:
    for version in versions:
        for qubit in qubits:
            N = int(qubit)
            if version == "ori" and N < 27:
                continue
            ini_path = os.path.join("MPI_circuit",version,qubit,"myini.ini")
            state_paths = [f"/mnt/state{i % 4}/path{i}" for i in range(1 << file_seg)]
            state_paths = ",".join(state_paths)
            args = Args(total_qbit=N, file_qbit=file_seg, chunk_qbit=chunk_seg,mpi_qbit=1,
                is_subcircuit = 1 if version == "sub" else 0,
                runner_type="MPI", state_paths=state_paths)
            ini = Ini(args, ini_path)
            ini.out()
            os.system(f"scp -P 9048 {ini_path} rdma2:~/SubQuokka_dev/src/correctness/{ini_path}")
            for file in files:
                print(protocol,version,file + qubit,"Start")
                cir_path = os.path.join("MPI_circuit",version,qubit,file + qubit + ".txt")
                log_path = os.path.join("MPI_circuit",version,qubit,file + qubit + f".{protocol}.log")
                if protocol == "RDMA":
                    os.system(f"mpirun -x UCX_NET_DEVICES={RDMA_CARD1},{RDMA_CARD2} -x LD_LIBRARY_PATH --bind-to none --hostfile ../hf --map-by ppr:1:node \"$(pwd)/../Quokka\" -i {ini_path} -c {cir_path} > {log_path}")
                elif protocol == "RDMA_TCP":
                    os.system(f"UCX_TLS=tcp mpirun -x UCX_NET_DEVICES={RDMA_TCP_CARD1},{RDMA_TCP_CARD2} -x LD_LIBRARY_PATH --bind-to none --hostfile ../hf --map-by ppr:1:node \"$(pwd)/../Quokka\" -i {ini_path} -c {cir_path} > {log_path}")
                elif protocol == "Pure_TCP":
                    os.system(f"mpirun -x UCX_NET_DEVICES={MORMAL_TCP_CARD1},{MORMAL_TCP_CARD2} -x LD_LIBRARY_PATH --bind-to none --hostfile ../hf --map-by ppr:1:node \"$(pwd)/../Quokka\" -i {ini_path} -c {cir_path} > {log_path}")
                else:
                    print("Protocol Error")
                    exit(-1)
            os.system("rm state/path*")
            print(f"{protocol}",version,file + qubit,"End")