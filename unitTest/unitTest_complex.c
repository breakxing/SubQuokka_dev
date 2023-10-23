#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>

#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <omp.h>
#include <libgen.h>
#include <dirent.h>
#include <sys/stat.h>
double total = 0.0;

int loop_range = 1ULL << (30 - 12 - 5);
int chunk_state = 1ULL << 12;
unsigned long long chunk_size = (1ULL << 12) * sizeof(double complex);
int fd_arr[32];
void *buffer[32];
unsigned long long fd_offset[32] = {0};

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

    omp_set_num_threads(32);
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int fd_off = 0;
        for (int i = 0; i < loop_range; i++)
        {
            pwrite(fd_arr[tid], buffer[tid], chunk_size, fd_off);
            fd_off += chunk_size;
        }
    }
    ((double complex *)(buffer[0]))[0] += 1;
    if(pwrite(fd_arr[0], buffer[0], sizeof(double complex), 0));
}

void H_gate(double complex *buffer){
    int up_off = 0;
    int lo_off = 1;
    double complex q0;
    double complex q1;
    double coeff = 1./sqrt(2);
    for (int i = 0; i < chunk_state; i += 2) {
        for (int j = 0; j < 1; j++){
            q0 = buffer[up_off];
            q1 = buffer[lo_off];
            // printf("结果: %lf + %lfi\n", creal(q0), cimag(q0));
            // printf("结果: %lf + %lfi\n", creal(q1), cimag(q1));
            // // printf("结果: %f + %fi\n", creal(q0+q1), cimag(q0+q1));
            // fflush(stdout);

            buffer[up_off] = coeff * (q0 + q1);
            buffer[lo_off] = coeff * (q0 - q1);
            up_off += 1;
            lo_off += 1;
        }
        up_off += 1;
        lo_off += 1;
    }

}


void run() {
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();

        for (int i = 0; i < loop_range; i++)
        {
            pread(fd_arr[tid], buffer[tid], chunk_size, fd_offset[tid]);
            H_gate((double complex *)(buffer[tid]));
            
    // printf("safe\n");
    // fflush(stdout);
            pwrite(fd_arr[tid], buffer[tid], chunk_size, fd_offset[tid]);
            fd_offset[tid] += chunk_size;
        }
        fd_offset[tid] = 0;
        #pragma omp barrier
    }
}

int main(){
    init();
    double st = omp_get_wtime();
    run();
    run();
    run();
    run();
    run();
    run();
    run();
    run();
    run();
    run();
    run();
    run();
    double ed = omp_get_wtime();
    printf("Time: %lf s\n", ed - st);
    return 0;
    
}