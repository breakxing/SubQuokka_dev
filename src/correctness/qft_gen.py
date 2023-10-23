from circuit_generator import *
from ini_generator import *
from test_util import *

# cir_path = "origin_qft21.txt"
# cir = get_circuit()
# # pseudo_qft(cir, 0, 20)
# qft(cir, 0, 21)
# create_circuit(cir, cir_path)


# cir_path = "origin_qft24.txt"
# cir = get_circuit()
# # pseudo_qft(cir, 0, 23)
# qft(cir, 0, 24)
# create_circuit(cir, cir_path)


# cir_path = "origin_qft27.txt"
# cir = get_circuit()
# # pseudo_qft(cir, 0, 26)
# qft(cir, 0, 27)
# create_circuit(cir, cir_path)


# cir_path = "origin_qft30.txt"
# cir = get_circuit()
# # pseudo_qft(cir, 0, 29)
# qft(cir, 0, 30)
# create_circuit(cir, cir_path)


cir_path = "origin_qft32.txt"
cir = get_circuit()
qft(cir, 0, 32)
create_circuit(cir, cir_path)