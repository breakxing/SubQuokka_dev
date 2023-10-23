qubit = 30

h_str = ""
for i in range(qubit):
    h_str = "H " + str(i)
    print(h_str)
    for j in range(i +1, qubit):
        cz_str = "CZ " + str(j) + " " + str(i)
        print(cz_str)

half_qubit = int(qubit / 2)
for i in range(half_qubit):
    cz_str = "CZ " + str(i) + " " + str(qubit - 1 - i)
    print(cz_str)