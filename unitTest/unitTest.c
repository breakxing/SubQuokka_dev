#include <stdlib.h>
#include <stdio.h>

#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <omp.h>
#include <libgen.h>
#include <dirent.h>
#include <sys/stat.h>
double total = 0.0;

double loop_range = 1ULL << (30 - 12);
double chunk_state = 1ULL << 12;
double chunk_size = (1ULL << 12) * 2 * sizeof(double);
int fd_arr[32];
void *buffer[32];
int fd_offset[32] = {0};

void init(){
    fd_arr[0] = open("state/path0", O_RDWR);
    fd_arr[1] = open("state/path1", O_RDWR);
    fd_arr[2] = open("state/path2", O_RDWR);
    fd_arr[3] = open("state/path3", O_RDWR);
    fd_arr[4] = open("state/path4", O_RDWR);
    fd_arr[5] = open("state/path5", O_RDWR);
    fd_arr[6] = open("state/path6", O_RDWR);
    fd_arr[7] = open("state/path7", O_RDWR);
    fd_arr[8] = open("state/path8", O_RDWR);
    fd_arr[9] = open("state/path9", O_RDWR);
    fd_arr[10] = open("state/path10", O_RDWR);
    fd_arr[11] = open("state/path11", O_RDWR);
    fd_arr[12] = open("state/path12", O_RDWR);
    fd_arr[13] = open("state/path13", O_RDWR);
    fd_arr[14] = open("state/path14", O_RDWR);
    fd_arr[15] = open("state/path15", O_RDWR);
    fd_arr[16] = open("state/path16", O_RDWR);
    fd_arr[17] = open("state/path17", O_RDWR);
    fd_arr[18] = open("state/path18", O_RDWR);
    fd_arr[19] = open("state/path19", O_RDWR);
    fd_arr[20] = open("state/path20", O_RDWR);
    fd_arr[21] = open("state/path21", O_RDWR);
    fd_arr[22] = open("state/path22", O_RDWR);
    fd_arr[23] = open("state/path23", O_RDWR);
    fd_arr[24] = open("state/path24", O_RDWR);
    fd_arr[25] = open("state/path25", O_RDWR);
    fd_arr[26] = open("state/path26", O_RDWR);
    fd_arr[27] = open("state/path27", O_RDWR);
    fd_arr[28] = open("state/path28", O_RDWR);
    fd_arr[29] = open("state/path29", O_RDWR);
    fd_arr[30] = open("state/path30", O_RDWR);
    fd_arr[31] = open("state/path31", O_RDWR);

    for (int i = 0; i < 32; i++){
        buffer[i] = malloc(chunk_size);
        memset(buffer[i], 0, chunk_size);
    }
}
void H_gate(void *buffer){
    int up_off = 0;
    int lo_off = 1;
    std::complex<double> q0;
    std::complex<double> q1;
    double coeff = 1./sqrt(2);

    for (int i = 0; i < chunk_state; i += target_offset) {
        for (int j = 0; j < half_target_offset; j++){
            q0 = buffer[up_off];
            q1 = buffer[lo_off];
            buffer[up_off] = coeff * (q0 + q1);
            buffer[lo_off] = coeff * (q0 - q1);
            up_off += 1;
            lo_off += 1;
        }
        up_off += half_target_offset;
        lo_off += half_target_offset;
    }
}

void run() {
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        #pragma omp barrier

        for (int i = 0; i < loop_range; i += chunk_state)
        {
            pread(fd_arr[tid], buffer[tid], chunk_state, fd_offset[tid]);
            H_gate(&buffer[tid]);
            pwrite(fd_arr[tid], buffer[tid], chunk_state, fd_offset[tid]);
            fd_offset[tid] += chunk_size;
        }
        
        #pragma omp barrier
    }
}

int main(){
    init();

    double st = omp_get_wtime();
    run();
    double ed = omp_get_wtime();
    printf("Time: %lf s\n", ed - st);
    return 0;
    
}