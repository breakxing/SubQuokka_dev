import os
from cir_utils import *

#for quBit in [21, 24, 27, 30]:
for quBit in range(21, 40):
    cir = get_circuit()
    '''
    cir_type = "h"
    for i in range(quBit):
        H(cir, i)
    '''
    '''
    cir_type = "qaoa"
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
    cir_type = "qft"
    for i in range(quBit):
        H(cir, i)
        angle = math.pi/2
        for j in range(i+1, quBit):
            CPhase(cir, i, j, angle)
            angle /= 2
    for i in range(int(quBit/2)):
        SWAP(cir, i, quBit-1-i)
    

    dir_path = "../circuit/"

    cir_path = dir_path + "ori/" + str(quBit) + "/" + cir_type + str(quBit) + ".txt"
    create_circuit(cir, cir_path)
    sub_path = dir_path + "sub/" + str(quBit) + "/" + cir_type + str(quBit) + ".txt"
    opt_path = dir_path + "ideal/" + str(quBit) + "/" + cir_type + str(quBit) + ".txt"

    print(f"../finder_chunk12_2 {cir_path} 0 {sub_path} > /dev/null")
    os.system(f"../finder_chunk12_2 {cir_path} 0 {sub_path} > /dev/null")
    os.system(f"../finder_chunk12_woVswap_2 {cir_path} 0 {opt_path} > /dev/null")
