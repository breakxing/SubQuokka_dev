from circuit_generator import *
from ini_generator import *
from test_util import *
import os


file_seg = 3
chunk_seg = 12

ini_path = "res.ini"
state_paths = [f"./state/path{i}" for i in range(1 << file_seg)]
state_paths = ",".join(state_paths)

cir_path = "./circuit/ori/21/qft21.txt"
cir = get_circuit()
N = 21

os.chdir("..")
os.system("make")
os.chdir("correctness")
for n in [N]:
    args = Args(total_qbit=n, file_qbit=file_seg, chunk_qbit=chunk_seg,mpi_qbit=2,
                is_subcircuit=0,MPI_testing=1,
                runner_type="MPI", state_paths=state_paths)
    ini = Ini(args, ini_path)
    ini.out()
    if args.mpi_qbit == 0:
        os.system(f"../Quokka -c ../../subcircuitFinder/mpi_test.txt -i {ini_path}")
        simple_test(f"[H {n} subcircuit]", False, cir_path, state_paths, n, (1 << file_seg))
    else:
        # os.chdir("..")
        # os.system("scp -r -P 9037 ./*.cpp rdma2:~/SubQuokka/src/")
        # os.system("scp -r -P 9037 ./*.hpp rdma2:~/SubQuokka/src/")
        # os.system("scp -r -P 9037 ./*.h rdma2:~/SubQuokka/src/")
        # os.system(f"scp -r -P 9037 ./correctness/{ini_path} rdma2:~/SubQuokka/src/correctness/")
        # os.system(f"scp -r -P 9037 ./correctness/{cir_path} rdma2:~/SubQuokka/src/correctness/")
        # os.system("scp -r -P 9037 ./Quokka rdma2:~/SubQuokka/src/")
        # os.chdir("correctness")
        # os.system(f"mpirun -x LD_LIBRARY_PATH --bind-to none --hostfile ../hf --map-by ppr:1:node \"$(pwd)/../Quokka\" -i {ini_path} -c {cir_path}")
        # os.system(f"mpirun -np {1 << int(args.mpi_qbit)} ../Quokka -i {ini_path} -c ../../subcircuitFinder/mpi_test.txt")
        os.system(f"mpirun -np {1 << int(args.mpi_qbit)} ../Quokka -i {ini_path} -c {cir_path}")
        for i in range(8,16):
            # os.system(f"scp -P 9037 paslab@140.112.90.37:~/SubQuokka/src/correctness/state1/path{i - 8} ./state0/path{i}")
            os.system(f"cp ./state1/path{i - 8} ./state0/path{i}")
            os.system(f"cp ./state2/path{i - 8} ./state0/path{i + 8}")
            os.system(f"cp ./state3/path{i - 8} ./state0/path{i + 16}")
        state_path_compare = [f"./state0/path{i}" for i in range((1 << file_seg) << int(args.mpi_qbit))]
        state_path_compare = ",".join(state_path_compare)
        simple_test(f"[H {n} subcircuit]", False,cir_path, state_path_compare, n, (1 << file_seg) * (1 << args.mpi_qbit))


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