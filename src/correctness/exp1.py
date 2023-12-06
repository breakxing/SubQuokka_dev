from circuit_generator import *
from ini_generator import *
from test_util import *
import os
from scipy.stats import unitary_group

file_seg = 5
chunk_seg = 12

ini_path = "res.ini"
state_paths = [f"/mnt/state{i % 4}/path{i}" for i in range(1 << file_seg)]
# state_paths = [f"./state/path{i}" for i in range(1 << file_seg)]
state_paths = ",".join(state_paths)

# cir_path = "./MPI_circuit/ori/27/qft27.txt"
cir_path = "cir.txt"
cir = get_circuit()
N = 30
mpiqubit = 1

for i in range(50):
    H(cir,N - 1)



create_circuit(cir, cir_path)
# exit(0)
os.chdir("..")
# os.system("make clean")
os.system("make -j")
os.chdir("correctness")

for n in [N]:
    args = Args(total_qbit=n, file_qbit=file_seg, chunk_qbit=chunk_seg,mpi_qbit=mpiqubit,
                is_subcircuit=0,MPI_testing=0,
                runner_type="MPI" if mpiqubit != 0 else "IO", state_paths=state_paths)
    ini = Ini(args, ini_path)
    ini.out()
    if args.mpi_qbit == 0:
        os.system(f"../Quokka -c {cir_path} -i {ini_path}")
        simple_test(f"[H {n} subcircuit]", False, cir_path, state_paths, n, (1 << file_seg))
    else:
        os.chdir("..")
        os.system("scp -r -P 9048 ./*.cpp rdma2:~/SubQuokka_dev/src/")
        os.system("scp -r -P 9048 ./*.hpp rdma2:~/SubQuokka_dev/src/")
        os.system("scp -r -P 9048 ./*.h rdma2:~/SubQuokka_dev/src/")
        os.system("scp -r -P 9048 ./makefile rdma2:~/SubQuokka_dev/src/")
        os.system(f"scp -r -P 9048 ./correctness/{ini_path} rdma2:~/SubQuokka_dev/src/correctness/")
        os.system(f"scp -r -P 9048 ./correctness/{cir_path} rdma2:~/SubQuokka_dev/src/correctness/")
        input("Go other computer to make")
        os.chdir("correctness")
        os.system(f"mpirun -x UCX_NET_DEVICES=mlx5_1:1,mlx5_0:1 -x LD_LIBRARY_PATH --bind-to none --hostfile ../hf --map-by ppr:1:node \"$(pwd)/../Quokka\" -i {ini_path} -c {cir_path}")
        # os.system(f"mpirun -np {1 << mpiqubit} --bind-to none ../Quokka -i {ini_path} -c {cir_path}")
        # for i in range(8,16):
        #     os.system(f"cp ./state1/path{i - 8} ./state0/path{i}")
        #     os.system(f"cp ./state2/path{i - 8} ./state0/path{i + 8}")
        #     os.system(f"cp ./state3/path{i - 8} ./state0/path{i + 16}")
        # state_path_compare = [f"./state0/path{i}" for i in range((1 << file_seg) * (1 << args.mpi_qbit))]
        # state_path_compare = ",".join(state_path_compare)
        # simple_test(f"[H {n} subcircuit]", False, cir_path, state_path_compare, n, (1 << file_seg) * (1 << args.mpi_qbit))


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