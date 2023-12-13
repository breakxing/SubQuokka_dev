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
# os.system("make clean")
os.system("make -j")
os.system("scp -r -P 9048 ./*.cpp rdma2:~/SubQuokka_dev/src/")
os.system("scp -r -P 9048 ./*.hpp rdma2:~/SubQuokka_dev/src/")
os.system("scp -r -P 9048 ./makefile rdma2:~/SubQuokka_dev/src/")
os.system("scp -r -P 9048 ./*.h rdma2:~/SubQuokka_dev/src/")
input("Go to other computer to make")
os.chdir("correctness")
file_seg = 5
chunk_seg = 12

cir_path = "cir.txt"
ini_path = "res.ini"
state_paths = [f"/mnt/state{i % 4}/path{i}" for i in range(1 << file_seg)]
state_paths = ",".join(state_paths)
cir = get_circuit()
N = 33
mpiqubit = 1


H(cir,N - 1)
create_circuit(cir, cir_path)
os.system(f"scp -P 9048 {cir_path} rdma2:~/SubQuokka_dev/src/correctness/{cir_path}")


buffer_sizes = []
for i in range(0,15):
    buffer_sizes.append(1 << i)



version = "mpi_qulacs"

for protocol in ["RDMA","RDMA_TCP","Pure_TCP"]:
    for buffer_size in buffer_sizes:
        args = Args(total_qbit=N, file_qbit=file_seg, chunk_qbit=chunk_seg,mpi_qbit=1,
            is_subcircuit = 0,
            runner_type="MPI", state_paths=state_paths,MPI_buffe_size=buffer_size)
        ini = Ini(args, ini_path)
        ini.out()
        log_path = os.path.join("exp_MPI_buffer_size_logs",version,str(N),f"{protocol}_Buffer_Size_{buffer_size}.log")
        os.system(f"scp -P 9048 {ini_path} rdma2:~/SubQuokka_dev/src/correctness/{ini_path}")
        print(f'{version}_{protocol}_{buffer_size} start')
        if protocol == "RDMA":
            os.system(f"mpirun -x UCX_NET_DEVICES={RDMA_CARD1},{RDMA_CARD2} -x LD_LIBRARY_PATH --bind-to none --hostfile ../hf --map-by ppr:1:node \"$(pwd)/../Quokka\" -i {ini_path} -c {cir_path} > {log_path}")
        elif protocol == "RDMA_TCP":
            os.system(f"UCX_TLS=tcp mpirun -x UCX_NET_DEVICES={RDMA_TCP_CARD1},{RDMA_TCP_CARD2} -x LD_LIBRARY_PATH --bind-to none --hostfile ../hf --map-by ppr:1:node \"$(pwd)/../Quokka\" -i {ini_path} -c {cir_path} > {log_path}")
        elif protocol == "Pure_TCP":
            os.system(f"mpirun -x UCX_NET_DEVICES={MORMAL_TCP_CARD1},{MORMAL_TCP_CARD2} -x LD_LIBRARY_PATH --bind-to none --hostfile ../hf --map-by ppr:1:node \"$(pwd)/../Quokka\" -i {ini_path} -c {cir_path} > {log_path}")
        print(f'{version}_{protocol}_{buffer_size} end')