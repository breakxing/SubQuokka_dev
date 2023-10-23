
def gen_swap(swap_in, swap_out):
    swap_size = len(swap_in)
    str2 = "VSWAP_2_2 "
    str1 = "SWAP"
    str_ = ""
    #for i in range(swap_size):
    i = 0
    while(i < swap_size):
        if i + 2 <= swap_size:
           str_ = str2
           str_ += " " + str(swap_in[i]) + " " + str(swap_in[i+1]) + " " + str(swap_out[i]) + " " + str(swap_out[i+1])
           i = i + 1
        else:
           str_ = str1
           str_ += " " + str(swap_in[i]) + " " + str(swap_out[i])
        print("1")
        print(str_)
        i = i + 1

           

swap_in = [12, 13, 14, 15, 16, 17, ]
swap_out = [6, 7, 8, 9, 10, 11, ]
gen_swap(swap_in, swap_out)
print("-------")



swap_in = [6, 7, 8, 9, 10, 11, ]
swap_out = [0, 1, 2, 3, 4, 5, ]
gen_swap(swap_in, swap_out)
print("-------")

swap_in = [0, 1, 2, 3, 4, 5, 18, 19, 20, 21, 22, 23, ]
swap_out = [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, ]
gen_swap(swap_in, swap_out)
print("-------")

swap_in = [6, 7, 8, 9, 10, 11, ]
swap_out = [0, 1, 2, 3, 4, 5, ]
gen_swap(swap_in, swap_out)
print("-------")

swap_in = [12, 13, 14, 15, 16, 17, ]
swap_out = [6, 7, 8, 9, 10, 11, ]
gen_swap(swap_in, swap_out)
print("-------")


swap_in = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 24, ]
swap_out = [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, ]
gen_swap(swap_in, swap_out)
print("-------")

swap_in = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ]
swap_out = [0, 1, 2, 3, 5, 6, 7, 8, 9, 10, ]
gen_swap(swap_in, swap_out)
print("-------")


swap_in = [0, 1, 2, 3, 9, 10, 21, 22, 23, ]
swap_out = [4, 11, 12, 13, 16, 17, 18, 19, 20, ]
gen_swap(swap_in, swap_out)
print("-------")


swap_in = [5, 6, 7, 8, 16, 17, 18, 19, ]
swap_out = [0, 1, 2, 3, 9, 10, 14, 15, 21, 22, 23, 24, ]
gen_swap(swap_in, swap_out)
print("-------")