# def gen_swap(swap_in, swap_out):
#     swap_size = len(swap_in)
#     str2 = "VSWAP_2_2 "
#     str1 = "SWAP"
#     str_ = ""
#     #for i in range(swap_size):
#     i = 0
#     while(i < swap_size):
#         if i + 2 <= swap_size:
#            str_ = str2
#            str_ += " " + str(swap_in[i]) + " " + str(swap_in[i+1]) + " " + str(swap_out[i]) + " " + str(swap_out[i+1])
#            i = i + 1
#         else:
#            str_ = str1
#            str_ += " " + str(swap_in[i]) + " " + str(swap_out[i])
#         print("1")
#         print(str_)
#         i = i + 1

def process(line, mapping):
    res = line.split()
    if res[0] == "H":
        res[1] = str(mapping[int(res[1])])

    elif res[0] == "U2" or res[0] == "SWAP":
        q0 = mapping[int(res[1])]
        q1 = mapping[int(res[2])]
        if q0 > q1:
            tmp = q0
            q0 = q1
            q1 = tmp
        res[1] = str(q0)
        res[2] = str(q1)
    
    return " ".join(res) + "\n"

dict_ = {}

filename = "/home/paslab/Downloads/SubQuokka/src/qft/qft21.txt"
fileInfo = "QQ.txt"
fileWrite = "qqqqqqqq.txt"

with open(filename) as f:
    lines =f.readlines()


with open(fileInfo) as f:
    lineInfos =f.readlines()


idx = 0
for l in lines:
    ltmp = l.split()[0]
    key = ltmp + str(idx)
    idx += 1
    dict_[key] = l

idx = 0
lineInfos_size = len(lineInfos)
with open(fileWrite, 'w+') as f:
    chunk = list(range(40))
    while(idx < lineInfos_size):
        if idx % 3 == 0:
            l_array = lineInfos[idx].split()
            print(l_array)
            f.write(str(len(l_array)) + "\n")
            for ele in l_array:
                line = dict_[ele]
                line = process(line, chunk)
                f.write(line)
            idx += 1
        else:
            l1_array = lineInfos[idx].split(",")
            l1_array = l1_array[:-1]
            print(l1_array)
            l2_array = lineInfos[idx + 1].split(",")
            l2_array = l2_array[:-1]
            swap_size = len(l1_array)
            str2 = "VSWAP_2_2 "
            str1 = "SWAP "
            str_ = ""
            #for i in range(swap_size):
            i = 0
            while(i < swap_size):
                if i + 2 <= swap_size:
                    str_ = str2
                    
                    q0 = chunk[int(l2_array[i])]
                    q1 = chunk[int(l2_array[i + 1])]
                    q2 = chunk[int(l1_array[i])]
                    q3 = chunk[int(l1_array[i + 1])]


                    if q0 > q1:
                        tmp = q0
                        q0 = q1
                        q1 = tmp
                    if q2 > q3:
                        tmp = q2
                        q2 = q3
                        q3 = tmp

                    if q0 > q2:
                        tmp = q0
                        q0 = q2
                        q2 = tmp
              
                    if q1 > q3:
                        tmp = q1
                        q1 = q3
                        q3 = tmp

                    str_ += " " + str(q0) + " " + str(q1) + " " + str(q2) + " " + str(q3)

                    tmp1 = chunk[int(l1_array[i])]
                    chunk[int(l1_array[i])] = chunk[int(l2_array[i])]
                    chunk[int(l2_array[i])] = tmp1

                    tmp2 = chunk[int(l1_array[i+1])]
                    chunk[int(l1_array[i+1])] = chunk[int(l2_array[i+1])]
                    chunk[int(l2_array[i+1])] = tmp2
                    i = i + 1


                else:
                    str_ = str1
                    q0 = chunk[int(l2_array[i])]
                    q1 = chunk[int(l1_array[i])]
                    if q0 > q1:
                        tmp = q0
                        q0 = q1
                        q1 = tmp
                    str_ += " " + str(q0) + " " + str(q1)

                    tmp1 = chunk[int(l1_array[i])]
                    chunk[int(l1_array[i])] = chunk[int(l2_array[i])]
                    chunk[int(l2_array[i])] = tmp1

                #print("1")


                #f.write("1\n")
                #f.write(str_ + "\n")
                i = i + 1
            idx += 2
        

