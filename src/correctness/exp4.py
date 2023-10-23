from circuit_generator import *
from ini_generator import *
from test_util import *

# Test for QFT gate on SubQuokka

# N = 21, 24, 27, 30
numFILES = 6
chunk_seg = 12

ini_path = "res.ini"
# cir_path = "../subcircuit.txt"
# cir_path = "../sub2.txt"
# cir_path = "../sub3.txt"
# cir_path = "../sub4.txt"
# cir_path = "../sub5.txt"
state_paths = "./state/path0,./state/path1,./state/path2,./state/path3,./state/path4,./state/path5,./state/path6,./state/path7,./state/path8,./state/path9,./state/path10,./state/path11,./state/path12,./state/path13,./state/path14,./state/path15,./state/path16,./state/path17,./state/path18,./state/path19,./state/path20,./state/path21,./state/path22,./state/path23,./state/path24,./state/path25,./state/path26,./state/path27,./state/path28,./state/path29,./state/path30,./state/path31,./state/path32,./state/path33,./state/path34,./state/path35,./state/path36,./state/path37,./state/path38,./state/path39,./state/path40,./state/path41,./state/path42,./state/path43,./state/path44,./state/path45,./state/path46,./state/path47,./state/path48,./state/path49,./state/path50,./state/path51,./state/path52,./state/path53,./state/path54,./state/path55,./state/path56,./state/path57,./state/path58,./state/path59,./state/path60,./state/path61,./state/path62,./state/path63"

# subcir_path_prefix = "../qft/opt"
# cir = get_circuit()
# subcir = get_circuit()
# for _ in range(10):
#     H(subcir, 0)
# cir.append(subcir)
# create_subcircuit(cir, subcir_path)


# MEM part
# # -circuit
# for n in [24, 24, 24, 24]:

cir_path = "cir_phase_30.txt"
sub_path = "sub_phase.txt"

# for i in range(20):
#     angle = math.pi/2
#     for j in range(i+1, 20):
#         CPhase(cir, i, j, angle)
#         angle /= 2
# for i in range(10):
#     SWAP(cir, i, 20-1-i)
# for i in range(20):
#     H(cir, i)
# for i in range(25):
#     angle = math.pi/2
#     for j in range(i+1, 25):
#         CPhase(cir, i, j, angle)
#         angle /= 2
# for i in range(12):
#     SWAP(cir, i, 25-1-i)
# for i in range(25):
#     H(cir, i)
#for quBit in range(21, 31):3#[21, 24, 27, 30]:
for quBit in [21, 24, 27, 30]:
#for quBit in [30]:
    cir = get_circuit()
    '''
    for i in range(quBit):
        H(cir, i)
    '''
    
    for i in range(quBit):
        H(cir, i)
    for p in range(5):
        for i in range(quBit-1):
            angle = math.pi
            for j in range(i+1, quBit):
                RZZ(cir, i, j, angle)
                #angle /= 2
        for i in range(quBit):
            angle = math.pi
            RX(cir, i, angle)
    
    
    '''
    for i in range(quBit):
        H(cir, i)
        angle = math.pi/2
        for j in range(i+1, quBit):
            CPhase(cir, i, j, angle)
            angle /= 2
    for i in range(int(quBit/2)):
        SWAP(cir, i, quBit-1-i)
    '''

    cir_path = "./mycircuit/ori_qaoa"+ str(quBit) + ".txt"
    create_circuit(cir, cir_path)

    continue
    sub_path = "sub_qft" + str(quBit) + ".txt" 
    os.system(f"../../subcircuitFinder/gpu_ver {cir_path} 0 {sub_path} > /dev/null")
    sub_path = "opt_qft" + str(quBit) + ".txt"
    os.system(f"../../subcircuitFinder/gpu_ver_opt {cir_path} 0 {sub_path} > /dev/null")

    continue
    
    #args = Args(total_qbit=quBit, file_qbit=6, chunk_qbit=12, 
    #            is_subcircuit=1, dump_file="")
    
    '''
    args = Args(total_qbit=quBit, file_qbit=6, chunk_qbit=12, 
                is_subcircuit=0, dump_file="")
    ini = Ini(args, ini_path)
    ini.out()
    os.system(f"../Quokka -c {cir_path} -i {ini_path}")
    '''
    
    sub_path = "sub_qft" + str(quBit) + ".txt"
   
    args = Args(total_qbit=quBit, file_qbit=6, chunk_qbit=12, 
                is_subcircuit=1, dump_file="")
    ini = Ini(args, ini_path)
    ini.out()
    
    
    #os.system(f"../../subcircuitFinder/opt_vs2 {cir_path} 0 {sub_path}")
    #os.system(f"../Quokka -c {sub_path} -i {ini_path}")

    sub_path = "opt_qft" + str(quBit) + ".txt"
    args = Args(total_qbit=quBit, file_qbit=6, chunk_qbit=12, 
                is_subcircuit=1, dump_file="")
    ini = Ini(args, ini_path)
    ini.out()
    os.system(f"../../subcircuitFinder/opt_vs2 {cir_path} 0 {sub_path}")
    #os.system(f"../Quokka -c {sub_path} -i {ini_path}")

    '''
    args = Args(total_qbit=quBit, file_qbit=6, chunk_qbit=12, 
        is_subcircuit=0, runner_type="IO", state_paths=state_paths)
    ini = Ini(args, ini_path)
    ini.out()
    os.system(f"../Quokka -c {cir_path} -i {ini_path}")
    
    
    args = Args(total_qbit=quBit, file_qbit=6, chunk_qbit=12, 
                is_subcircuit=1, runner_type="IO", dump_file="")
    ini = Ini(args, ini_path)
    ini.out()
    os.system(f"../../subcircuitFinder/finder {cir_path} 0 {sub_path}")
    os.system(f"../Quokka -c {sub_path} -i {ini_path}")

    args = Args(total_qbit=quBit, file_qbit=6, chunk_qbit=12, 
                is_subcircuit=1, runner_type="IO", dump_file="")
    ini = Ini(args, ini_path)
    ini.out()
    os.system(f"../../subcircuitFinder/finder_opt {cir_path} 0 {sub_path}")
    os.system(f"../Quokka -c {sub_path} -i {ini_path}")
    '''
    #args = Args(total_qbit=quBit, file_qbit=6, chunk_qbit=12, 
    #            is_subcircuit=0, runner_type="IO", state_paths=state_paths)

    
    #os.system(f"../../subcircuitFinder/finder {cir_path} 0 {sub_path}")
    # os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null")
    
    #os.system(f"../Quokka -c {cir_path} -i {ini_path}")
    #os.system(f"../Quokka -c {sub_path} -i {ini_path}")



    # simple_test_mem(f"[H {n} origin]", False, cir_path, "res.txt", N)
exit()

# # -subcircuit
for n in [21, 24, 27, 30]:
    subcir_path = subcir_path_prefix + str(n) + ".txt"

    args = Args(total_qbit=n, file_qbit=6, chunk_qbit=12, 
                is_subcircuit=1, dump_file="")
    ini = Ini(args, ini_path)
    ini.out()

    # os.system(f"../Quokka -c {subcir_path} -i {ini_path} > /dev/null")
    os.system(f"../Quokka -c {subcir_path} -i {ini_path}")
    # simple_test_mem(f"[H {n} subcircuit]", False, cir_path, "res.txt", N

# # exit()

# IO part
# -circuit
# # for n in [24, 24, 24, 24]:
# for n in [21, 24, 27, 30]:
#     cir_path = cir_path_prefix + str(n) + ".txt"
#     args = Args(total_qbit=n, file_qbit=6, chunk_qbit=12, 
#                 is_subcircuit=0,
#                 runner_type="IO", state_paths=state_paths)
#     ini = Ini(args, ini_path)
#     ini.out()

#     cir = get_circuit()
#     pseudo_qft(cir, 0, n-1)
#     create_circuit(cir, cir_path)

#     # os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null")
#     os.system(f"../Quokka -c {cir_path} -i {ini_path}")
#     # simple_test(f"[H {n} subcircuit]", False, cir_path, state_paths, N, numFILES)

# # -subcircuit
# for n in [21, 24, 27, 30]:
# # for n in [24, 27, 30]:
# # for n in [21, 24, 27]:
# # # for n in [21]:
#     subcir_path = subcir_path_prefix + str(n) + ".txt"
#     args = Args(total_qbit=n, file_qbit=6, chunk_qbit=12, 
#                 is_subcircuit=1,
#                 runner_type="IO", state_paths=state_paths)
#     ini = Ini(args, ini_path)
#     ini.out()

#     # os.system(f"../Quokka -c {cir_path} -i {ini_path} > /dev/null")
#     os.system(f"../Quokka -c {subcir_path} -i {ini_path}")
#     # simple_test(f"[H {n} subcircuit]", False, cir_path, state_paths, N, numFILES)

# exit()