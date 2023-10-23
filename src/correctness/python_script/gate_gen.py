import os
from cir_utils import *
import math

#Gate = RX
#cir_type = "rx"
cir_type = "h"
angle = math.pi

for quBit in range(21, 40):

    cir = get_circuit()
    for i in range(quBit):
        #Gate(cir, i, angle)
        H(cir, i)

    dir_path = "../circuit_io/"

    cir_path = dir_path + "ori/" + str(quBit) + "/" + cir_type + str(quBit) + ".txt"
    print("--->", cir_path)
    create_circuit(cir, cir_path)
    sub_path = dir_path + "sub/" + str(quBit) + "/" + cir_type + str(quBit) + ".txt"

    
    subcircuits = []
    cir = get_circuit()
    for i in range(12):
        #Gate(cir, i, angle)
        H(cir, i)
    subcircuits.append(cir)

    for i in range(12, quBit):
        cir = get_circuit()
        #Gate(cir, i, angle)
        H(cir, i)
        subcircuits.append(cir)

    create_subcircuit(subcircuits, sub_path)

    #os.system(f"../finder_chunk10 {cir_path} 0 {sub_path} > /dev/null")
