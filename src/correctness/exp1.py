from circuit_generator import *
from ini_generator import *
from test_util import *
import os
from scipy.stats import unitary_group
import random
file_seg = 3
chunk_seg = 12
RDMA_CARD1 = 'mlx5_1:1'
RDMA_CARD2 = 'mlx5_0:1'
RDMA_TCP_CARD1 = 'enp225s0f1np1'
RDMA_TCP_CARD2 = 'enp225s0f0np0'
MORMAL_TCP_CARD1 = 'eno1'
MORMAL_TCP_CARD2 = 'enp2s0'
ini_path = "res.ini"
# state_paths = [f"/mnt/state{i % 4}/path{i}" for i in range(1 << file_seg)]
state_paths = [f"./state/path{i}" for i in range(1 << file_seg)]
state_paths = ",".join(state_paths)

# cir_path = "../../subcircuitFinder/aaa.txt"
cir_path = "cir.txt"
cir = get_circuit()
N = 21
mpiqubit = 2
for i in range(N):
    H(cir,i)
Phase(cir,N - 1,random.random())
Phase(cir,N - 2,random.random())
SWAP(cir,N-1,N-2)
# random_unitary = unitary_group.rvs(4)
# coeff = []
# for x in random_unitary:
#     for y in x:
#         coeff.append(y.real)
#         coeff.append(y.imag)
# U2(cir,N - 2,N - 1,coeff)
# SWAP(cir,2,N - 1)
# for i in range(10):
    # H(cir,5)
# for i in  range(5):
# CPhase(cir,N-2,N-1,math.pi * random.random())
# for i in range(10)
# for i in range(5):
#     H(cir,29)
# for i in range(10):
#     H(cir,24)
# for i in range(2):
#     H(cir,17)

# H(cir,N-1)
# for i in range(5):
    # H(cir,12)
    # H(cir,N-1)
# for i in range(5):
#     H(cir,N-1)
# for i in range(5):
#     H(cir,N-1)
# H(cir,N - 1)
# H(cir,N - 2)
#     H(cir,29)
# CPhase(cir,12,15,math.pi)
# H(cir,33)
# H(cir,N-1)
# for i in range(N):
#     H(cir,i)
# for i in range(N):
#     H(cir,i)



create_circuit(cir, cir_path)
# exit(0)
os.chdir("..")
# os.system("make clean")
os.system("make -j")
os.chdir("correctness")

for n in [N]:
    args = Args(total_qbit=n, file_qbit=file_seg, chunk_qbit=chunk_seg,mpi_qbit=mpiqubit,
                is_subcircuit=0,MPI_buffe_size=64,MPI_testing=1,dump_file="res.txt",
                runner_type="MPI_IO" if mpiqubit != 0 else "IO", state_paths=state_paths)
    ini = Ini(args, ini_path)
    ini.out()
    if args.mpi_qbit == 0:
        os.system(f"../Quokka -c {cir_path} -i {ini_path}")
        simple_test(f"[H {n} subcircuit]", False, cir_path, state_paths, n, (1 << file_seg))
        # v = verifier(args, cir_path)
        # flag = v.run(f"Exp1")
    else:
        os.chdir("..")
        # os.system("scp -r -P 22 ./*.cpp rdma1:~/SubQuokka_dev/src/")
        # os.system("scp -r -P 22 ./*.hpp rdma1:~/SubQuokka_dev/src/")
        # os.system("scp -r -P 22 ./*.h rdma1:~/SubQuokka_dev/src/")
        # os.system("scp -r -P 22 ./makefile rdma1:~/SubQuokka_dev/src/")
        os.system(f"scp -r -P 22 ./correctness/{ini_path} rdma1:~/SubQuokka_dev/src/correctness/")
        os.system(f"scp -r -P 22 ./correctness/{cir_path} rdma1:~/SubQuokka_dev/src/correctness/")
        # input("Go other computer to make")
        os.chdir("correctness")
        # os.system(f"mpirun -x UCX_NET_DEVICES=mlx5_1:1,mlx5_0:1 -x LD_LIBRARY_PATH --bind-to none --hostfile ../hf --map-by ppr:1:node \"$(pwd)/../Quokka\" -i {ini_path} -c {cir_path}")
        # os.system(f"mpirun -x OMPI_MCA_btl_tcp_links -x UCX_NET_DEVICES={RDMA_TCP_CARD1},{RDMA_TCP_CARD2} -x LD_LIBRARY_PATH --bind-to none --hostfile ../hf --map-by ppr:1:node \"$(pwd)/../Quokka\" -i {ini_path} -c {cir_path}")
        # os.system(f"mpirun -x OMPI_MCA_btl_tcp_links -x UCX_NET_DEVICES={MORMAL_TCP_CARD1},{MORMAL_TCP_CARD2} -x LD_LIBRARY_PATH --bind-to none --hostfile ../hf --map-by ppr:1:node \"$(pwd)/../Quokka\" -i {ini_path} -c {cir_path}")
        os.system(f"mpirun -np {1 << mpiqubit} --bind-to none ../Quokka -i {ini_path} -c {cir_path}")
        for i in range(8,16):
            os.system(f"cp ./state1/path{i - 8} ./state0/path{i}")
            os.system(f"cp ./state2/path{i - 8} ./state0/path{i + 8}")
            os.system(f"cp ./state3/path{i - 8} ./state0/path{i + 16}")
        state_path_compare = [f"./state0/path{i}" for i in range((1 << file_seg) * (1 << args.mpi_qbit))]
        state_path_compare = ",".join(state_path_compare)
        simple_test(f"[H {n} subcircuit]", False, cir_path, state_path_compare, n, (1 << file_seg) * (1 << args.mpi_qbit))\
        # if args.dump_file != "":
        #     with open("res0.txt",'a') as rootfile:
        #         for i in range(1,(1 << mpiqubit)):
        #             filename = f'res{i}.txt'
        #             with open(filename,'r') as file:
        #                 for line in file.readlines():
        #                     rootfile.write(line)
        #             os.system(f"rm {filename}")
        #     os.system("mv res0.txt res.txt")
        #     v = verifier(args, cir_path)
        #     flag = v.run(f"Exp1")
        #     os.system(f"rm res.txt")


# -subcircuit
# # for n in range(21, 34):
# # for n in range(34, 36):
# for n in [32]:
# for n in [30]:
#     args = Args(total_qbit=n, file_qbit=6, chunk_qbit=12, 
#                 is_subcircuit=1,
#                 runner_type="IO", state_paths=state_paths)
#     ini = Ini(args, ini_path)
#     ini.out()

#     # os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null")
#     os.system(f"../Quokka_orig -c {subcir_path} -i {ini_path}")
#     os.system(f"../Quokka -c {subcir_path} -i {ini_path}")
    # simple_test(f"[H {n} subcircuit]", False, cir_path, state_paths, N, numFILES)

# exit()