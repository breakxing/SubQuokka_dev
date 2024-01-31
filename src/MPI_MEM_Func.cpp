#include <omp.h>
#include <unistd.h> //pread, pwrite
#include <utility>
#include <complex>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>
#include <algorithm>
#include "circuit.hpp"
#include <functional>
#include "MEM_Runner.hpp"
void MEM_Runner::MPI_gate_scheduler(thread_MEM_task &task,Gate * &g)
{
    // long long mpi_targ_mask_0 = 1 << (g->targs[1] - seg.N);
    long long mpi_targ_mask_2 = 1 << (g->targs[0] - seg.N);
    if(g->targs.size() == 1)
    {
        int partner_rank = env.rank ^ mpi_targ_mask_2;
        long long total_chunk_per_thread = env.thread_state / env.chunk_state;
        long long per_chunk = min(total_chunk_per_thread,(long long)env.MPI_buffer_size);
        task.partner_using = {partner_rank};
        if(g->name == "Z_Gate_MEM" || g->name == "Phase_Gate_MEM" || g->name == "RZ_Gate_MEM")
        {
            MPI_one_qubit_gate_diagonal(task,g);
            return;
        }
        for(int i = task.tid * env.thread_state;i < ((task.tid + 1) * env.thread_state);i+=per_chunk * env.chunk_state)
        {
            _mpi_one_gate_inner(task,g,i);
        }
    }
}
void MEM_Runner::MPI_one_qubit_gate_diagonal(thread_MEM_task &task,Gate* &g)
{
    if(g->name == "Z_Gate_MEM" || g->name == "Phase_Gate_MEM" || g->name == "RZ_Gate_MEM")
    {
        bool isupper = env.rank < task.partner_using[0];
        for(int i = task.tid * env.thread_state;i < ((task.tid + 1) * env.thread_state);i+=env.chunk_state)
        {
            g->run_one_qubit_mpi_mem_diagonal(buffer,i,isupper);
        }
        return;
    }
}
void MEM_Runner::_mpi_one_gate_inner(thread_MEM_task &task,Gate* &g,long long startIdx)
{
    int partner_rank = task.partner_using[0];
    long long per_chunk = min(env.thread_state / env.chunk_state,(long long)env.MPI_buffer_size);
    vector<complex<double>>*buffer1_ptr = &buffer;
    vector<complex<double>>*buffer2_ptr = &task.buffer2;
    int pos1 = startIdx;
    int pos2 = 0;
    if(env.rank > partner_rank)
    {
        swap(buffer1_ptr,buffer2_ptr);
        swap(pos1,pos2);
    }
    // MPI_Sendrecv(&buffer[startIdx],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,&task.buffer2[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Isend(&buffer[startIdx],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request_send[0]);
    MPI_Irecv(&task.buffer2[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request_recv[0]);
    MPI_Wait(&task.request_recv[0],MPI_STATUS_IGNORE);
    MPI_Wait(&task.request_send[0],MPI_STATUS_IGNORE);
    g->run_one_qubit_mpi_mem(*buffer1_ptr,*buffer2_ptr,pos1,pos2,per_chunk * env.chunk_state);
}