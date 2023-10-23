import os
from cir_utils import *
import math

Gate = RZZ
cir_type = "rzz"
angle = math.pi

foo = lambda cir, targ1, angle : {'circuit': cir, 
    'target0': 0,
    'target1': targ1,
    'angle': angle
}

for quBit in range(21, 40):

    cir = get_circuit()
    for i in range(1, quBit):
        Gate(cir, 0, i, angle)

    dir_path = "../circuit_io/"

    cir_path = dir_path + "ori/" + str(quBit) + "/" + cir_type + str(quBit) + ".txt"
    print("--->", cir_path)
    create_circuit(cir, cir_path)
    sub_path = dir_path + "sub/" + str(quBit) + "/" + cir_type + str(quBit) + ".txt"

    
    subcircuits = []
    cir = get_circuit()
    for i in range(1, 12):
        Gate(cir, 0, i, angle)
    subcircuits.append(cir)

    for i in range(12, quBit):
        cir = get_circuit()
        Gate(cir, 0, i, angle)
        subcircuits.append(cir)

    create_subcircuit(subcircuits, sub_path)

    #os.system(f"../finder_chunk10 {cir_path} 0 {sub_path} > /dev/null")
