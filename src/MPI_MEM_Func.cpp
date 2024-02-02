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
    if(g->targs.size() == 1)
    {
        long long rank_mask_0 = 1 << (g->targs[0] - seg.N);
        int partner_rank = env.rank ^ rank_mask_0;
        long long total_chunk_per_thread = env.thread_state / env.chunk_state;
        long long per_chunk = min(total_chunk_per_thread,(long long)env.MPI_buffer_size);
        task.partner_using = {partner_rank};
        long long startIdx;
        for(long long i = 0;i < env.thread_state;i+=per_chunk * env.chunk_state)
        {
            startIdx = task.tid * env.thread_state + i;
            _mpi_one_gate_inner(task,g,startIdx);
        }
    }
    else//U2
    {
        if(isMpi(g->targs[0]))
        {
            int rank_mask_0 = 1 << (g->targs[0] - seg.N);
            int rank_mask_1 = 1 << (g->targs[1] - seg.N);
            int rank0 = env.rank & (~rank_mask_0) & (~rank_mask_1);
            int rank1 = rank0 | rank_mask_0;
            int rank2 = rank0 | rank_mask_1;
            int rank3 = rank0 | rank_mask_0 | rank_mask_1;
            unordered_map<int,int>rank_mapping{{rank0,0},{rank1,1},{rank2,2},{rank3,3}};
            int rank_order = rank_mapping[env.rank];
            long long per_chunk = min(env.thread_state / env.chunk_state,(long long)env.MPI_buffer_size);
            vector<MPI_Request*>request_send = {&task.request1_send,&task.request2_send,&task.request3_send};
            vector<MPI_Request*>request_recv = {&task.request1_recv,&task.request2_recv,&task.request3_recv};
            vector<vector<complex<double>>*>buffer_recv_using = {&task.buffer2,&task.buffer3,&task.buffer4};
            vector<complex<double>>*buffer1_ptr = &buffer;
            vector<complex<double>>*buffer2_ptr = &task.buffer2;
            vector<complex<double>>*buffer3_ptr = &task.buffer3;
            vector<complex<double>>*buffer4_ptr = &task.buffer4;
            long long pos0 = 0;
            long long pos1 = 0;
            long long pos2 = 0;
            long long pos3 = 0;
            long long *local_pos;
            if(rank_order == 0)
            {
                task.partner_using = {rank1,rank2,rank3};
                local_pos = &pos0;
            }
            else if(rank_order == 1)
            {
                task.partner_using = {rank0,rank2,rank3};
                buffer1_ptr = &task.buffer2;
                buffer2_ptr = &buffer;
                local_pos = &pos1;
            }
            else if(rank_order == 2)
            {
                task.partner_using = {rank0,rank1,rank3};
                buffer1_ptr = &task.buffer2;
                buffer2_ptr = &task.buffer3;
                buffer3_ptr = &buffer;
                local_pos = &pos2;
            }
            else
            {
                task.partner_using = {rank0,rank1,rank2};
                buffer1_ptr = &task.buffer2;
                buffer2_ptr = &task.buffer3;
                buffer3_ptr = &task.buffer4;
                buffer4_ptr = &buffer;
                local_pos = &pos3;
            }
            for(long long i = 0;i <env.thread_state;i+=per_chunk * env.chunk_state)
            {
                *local_pos = task.tid * env.thread_state + i;
                for(int j = 0;j < 3;j++)
                {
                    MPI_Isend(&buffer[*local_pos],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[j],task.tid,MPI_COMM_WORLD,&(*request_send[j]));
                    MPI_Irecv(&(*buffer_recv_using[j])[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[j],task.tid,MPI_COMM_WORLD,&(*request_recv[j]));
                }
                for(int j = 0;j < 3;j++)
                {
                    MPI_Wait(&(*request_send[j]),MPI_STATUS_IGNORE);
                    MPI_Wait(&(*request_recv[j]),MPI_STATUS_IGNORE);
                }
                g->run_mem_u2(*buffer1_ptr,*buffer2_ptr,*buffer3_ptr,*buffer4_ptr,pos0,pos1,pos2,pos3,per_chunk * env.chunk_state);
            }
            // for(int i = task.tid * env.thread_state;i < ((task.tid + 1) * env.thread_state);i+=per_chunk * env.chunk_state)
            // {
            //     for(int j = 0;j < 3;j++)
            //     {
            //         MPI_Isend(&buffer[i],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[j],task.tid,MPI_COMM_WORLD,&(*request_send[j]));
            //         MPI_Irecv(&(*buffer_recv_using[j])[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[j],task.tid,MPI_COMM_WORLD,&(*request_recv[j]));
            //     }
            //     for(int j = 0;j < 3;j++)
            //     {
            //         MPI_Wait(&(*request_send[j]),MPI_STATUS_IGNORE);
            //         MPI_Wait(&(*request_recv[j]),MPI_STATUS_IGNORE);
            //     }
            //     *local_pos = i;
            //     g->run_mem_u2(*buffer1_ptr,*buffer2_ptr,*buffer3_ptr,*buffer4_ptr,pos0,pos1,pos2,pos3,per_chunk * env.chunk_state);
            // }
        }
        else if(isFile(g->targs[0]))
        {
            int rank_mask_1 = 1 << (g->targs[1] - seg.N);
            int thread_mask = 1 << (g->targs[0] - seg.middle - seg.chunk);
            int partner_rank = env.rank ^ rank_mask_1;
            int t0 = task.tid & (~thread_mask);
            int t1 = task.tid | thread_mask;
            bool thread_equal_chunk = (env.thread_state == env.chunk_state);
            long long total_chunk_per_thread = env.thread_state / env.chunk_state / 2;
            long long per_chunk = min(total_chunk_per_thread,(long long)env.MPI_buffer_size);
            long long pos0;
            long long pos1;
            long long pos2;
            long long pos3;
            long long startIdx0;
            long long startIdx1;
            vector<complex<double>>*buffer1_ptr;
            vector<complex<double>>*buffer2_ptr;
            vector<complex<double>>*buffer3_ptr;
            vector<complex<double>>*buffer4_ptr;
            if(thread_equal_chunk)
            {
                if(task.tid == t1) return;
                startIdx0 = t0 * env.thread_state;
                startIdx1 = t1 * env.thread_state;
                pos0 = startIdx0;
                pos1 = startIdx1;
                pos2 = 0;
                pos3 = 0;
                buffer1_ptr = &buffer;
                buffer2_ptr = &buffer;
                buffer3_ptr = &task.buffer2;
                buffer4_ptr = &task.buffer3;
                if(env.rank > partner_rank)
                {
                    swap(pos0,pos2);
                    swap(pos1,pos3);
                    swap(buffer1_ptr,buffer3_ptr);
                    swap(buffer2_ptr,buffer4_ptr);
                }
                MPI_Isend(&buffer[startIdx0],env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request1_send);
                MPI_Isend(&buffer[startIdx1],env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request2_send);
                MPI_Irecv(&task.buffer2[0],env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request1_recv);
                MPI_Irecv(&task.buffer3[0],env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request2_recv);
                MPI_Wait(&task.request1_send,MPI_STATUS_IGNORE);
                MPI_Wait(&task.request2_send,MPI_STATUS_IGNORE);
                MPI_Wait(&task.request1_recv,MPI_STATUS_IGNORE);
                MPI_Wait(&task.request2_recv,MPI_STATUS_IGNORE);
                g->run_mem_u2(*buffer1_ptr,*buffer2_ptr,*buffer3_ptr,*buffer4_ptr,pos0,pos1,pos2,pos3,env.chunk_state);
            }
            else
            {
                startIdx0 = t0 * env.thread_state + (task.tid != t0) * (env.thread_state >> 1);
                startIdx1 = t1 * env.thread_state + (task.tid != t0) * (env.thread_state >> 1);
                pos0 = startIdx0;
                pos1 = startIdx1;
                pos2 = 0;
                pos3 = 0;
                buffer1_ptr = &buffer;
                buffer2_ptr = &buffer;
                buffer3_ptr = &task.buffer2;
                buffer4_ptr = &task.buffer3;
                if(env.rank > partner_rank)
                {
                    swap(buffer1_ptr,buffer3_ptr);
                    swap(buffer2_ptr,buffer4_ptr);
                }
                for(long long i = 0;i < env.thread_state >> 1;i+=per_chunk * env.chunk_state)
                {
                    pos0 = i + startIdx0;
                    pos1 = i + startIdx1;
                    if(env.rank > partner_rank)
                    {
                        swap(pos0,pos2);
                        swap(pos1,pos3);
                    }
                    MPI_Isend(&buffer[i + startIdx0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request1_send);
                    MPI_Isend(&buffer[i + startIdx1],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request2_send);
                    MPI_Irecv(&task.buffer2[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request1_recv);
                    MPI_Irecv(&task.buffer3[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request2_recv);
                    MPI_Wait(&task.request1_send,MPI_STATUS_IGNORE);
                    MPI_Wait(&task.request2_send,MPI_STATUS_IGNORE);
                    MPI_Wait(&task.request1_recv,MPI_STATUS_IGNORE);
                    MPI_Wait(&task.request2_recv,MPI_STATUS_IGNORE);
                    g->run_mem_u2(*buffer1_ptr,*buffer2_ptr,*buffer3_ptr,*buffer4_ptr,pos0,pos1,pos2,pos3,per_chunk * env.chunk_state);
                }
            }
        }
        else if(isMiddle(g->targs[0]))
        {
            int rank_mask_1 = 1 << (g->targs[1] - seg.N);
            int partner_rank = env.rank ^ rank_mask_1;
            long long total_chunk_per_thread = env.qubit_offset[g->targs[0]] / env.chunk_state;
            long long per_chunk = min(total_chunk_per_thread,(long long)env.MPI_buffer_size / 2);
            long long pos0;
            long long pos1;
            long long pos2;
            long long pos3;
            vector<complex<double>>*buffer1_ptr;
            vector<complex<double>>*buffer2_ptr;
            vector<complex<double>>*buffer3_ptr;
            vector<complex<double>>*buffer4_ptr;
            for(long long i = 0;i < env.thread_state;i+=env.qubit_offset[g->targs[0] + 1])
            {
                for(long long j = 0;j < env.qubit_offset[g->targs[0]];j+=per_chunk * env.chunk_state)
                {
                    long long startIdx0 = task.tid * env.thread_state + i + j;
                    long long startIdx1 = task.tid * env.thread_state + i + j + env.qubit_offset[g->targs[0]];
                    pos0 = startIdx0;
                    pos1 = startIdx1;
                    pos2 = 0;
                    pos3 = 0;
                    buffer1_ptr = &buffer;
                    buffer2_ptr = &buffer;
                    buffer3_ptr = &task.buffer2;
                    buffer4_ptr = &task.buffer3;
                    MPI_Isend(&buffer[startIdx0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request1_send);
                    MPI_Isend(&buffer[startIdx1],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request2_send);
                    MPI_Irecv(&task.buffer2[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request1_recv);
                    MPI_Irecv(&task.buffer3[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request2_recv);
                    MPI_Wait(&task.request1_send,MPI_STATUS_IGNORE);
                    MPI_Wait(&task.request2_send,MPI_STATUS_IGNORE);
                    MPI_Wait(&task.request1_recv,MPI_STATUS_IGNORE);
                    MPI_Wait(&task.request2_recv,MPI_STATUS_IGNORE);
                    if(env.rank > partner_rank)
                    {
                        swap(buffer1_ptr,buffer3_ptr);
                        swap(buffer2_ptr,buffer4_ptr);
                        swap(pos0,pos2);
                        swap(pos1,pos3);
                    }
                    g->run_mem_u2(*buffer1_ptr,*buffer2_ptr,*buffer3_ptr,*buffer4_ptr,pos0,pos1,pos2,pos3,per_chunk * env.chunk_state);
                }
            }
        }
        else
        {
            long long total_chunk_per_thread = env.thread_state / env.chunk_state;
            long long per_chunk = min(total_chunk_per_thread,(long long)env.MPI_buffer_size);
            task.partner_using = {env.rank ^ (1 << (g->targs[1] - seg.N))};
            for(long long i = 0;i < env.thread_state;i+=per_chunk * env.chunk_state)
            {
                _mpi_one_gate_inner(task,g,i + task.tid * env.thread_state);
            }
        }
    }
}
void MEM_Runner::MPI_one_qubit_gate_diagonal(Gate* &g)
{
    int rank_mask_0 = 1 << (g->targs[0] - seg.N);
    int rank0 = env.rank & (~rank_mask_0);
    int rank1 = rank0 | rank_mask_0;
    unordered_map<int,int>rank_mapping{{rank0,0},{rank1,1}};
    int rank_order = rank_mapping[env.rank];
    #pragma omp parallel for schedule(static)
    for(long long i = 0;i < (1 << seg.N);i+=env.chunk_state)
    {
        g->run_one_qubit_mpi_mem_diagonal(buffer,i,rank_order,env.chunk_state);
    }
}
void MEM_Runner::MPI_two_qubit_gate_diagonal(Gate* &g)
{
    if(isMpi(g->targs[0]))
    {
        int rank_mask_0 = 1 << (g->targs[0] - seg.N);
        int rank_mask_1 = 1 << (g->targs[1] - seg.N);
        int rank0 = env.rank & (~rank_mask_0) & ~(rank_mask_1);
        int rank1 = rank0 | rank_mask_0;
        int rank2 = rank0 | rank_mask_1;
        int rank3 = rank0 | rank_mask_0 | rank_mask_1;
        unordered_map<int,int>rank_mapping{{rank0,0},{rank1,1},{rank2,2},{rank3,3}};
        int rank_order = rank_mapping[env.rank];
        if(g->name == "CPhase_Gate_MEM" && rank_order != 3) return;
        #pragma omp parallel for schedule(static)
        for(long long i = 0;i < (1 << seg.N);i+=env.chunk_state)
        {
            g->run_one_qubit_mpi_mem_diagonal(buffer,i,rank_order,env.chunk_state);
        }
    }
    else
    {
        long long stride = env.qubit_offset[g->targs[0] + 1];
        if(g->name == "CPhase_Gate_MEM")
        {
            if(env.rank != (env.rank | (1 << (g->targs[1] - seg.N)))) return;
            #pragma omp parallel for schedule(static)
            for(long long i = 0;i < (1 << seg.N);i+=stride)
            {
                g->run_one_qubit_mpi_mem_diagonal(buffer,i + (stride >> 1),3,stride >> 1);
            }
        }
        else
        {
            int rank_off = (env.rank == (env.rank | (1 << (g->targs[1] - seg.N)))) << 1;
            #pragma omp parallel for schedule(static)
            for(long long i = 0;i < (1 << seg.N);i+=stride)
            {
                g->run_one_qubit_mpi_mem_diagonal(buffer,i,rank_off,stride >> 1);
                g->run_one_qubit_mpi_mem_diagonal(buffer,i + (stride >> 1),rank_off | 1,stride >> 1);
            }
        }
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
    MPI_Isend(&buffer[startIdx],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request1_send);
    MPI_Irecv(&task.buffer2[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request1_recv);
    MPI_Wait(&task.request1_recv,MPI_STATUS_IGNORE);
    MPI_Wait(&task.request1_send,MPI_STATUS_IGNORE);
    g->run_one_qubit_mpi_mem(*buffer1_ptr,*buffer2_ptr,pos1,pos2,per_chunk * env.chunk_state);
}
void MEM_Runner::MPI_Swap_1_1(thread_MEM_task &task,Gate* &g)
{
    if(isMpi(g->targs[0]))
    {
        int rank_mask_0 = 1 << (g->targs[0] - seg.N);
        int rank_mask_1 = 1 << (g->targs[1] - seg.N);
        int rank0 = env.rank & (~rank_mask_0) & ~(rank_mask_1);
        int rank1 = rank0 | rank_mask_0;
        int rank2 = rank0 | rank_mask_1;
        int rank3 = rank0 | rank_mask_0 | rank_mask_1;
        unordered_map<int,int>rank_mapping{{rank0,0},{rank1,1},{rank2,2},{rank3,3}};
        int rank_order = rank_mapping[env.rank];
        if(rank_order == 0 || rank_order == 3) return;
        long long total_chunk_per_thread = env.thread_state / env.chunk_state;
        long long per_chunk = min(total_chunk_per_thread,(long long)env.MPI_buffer_size);
        task.partner_using = {env.rank == rank1? rank2 : rank1};
        long long startIdx;
        for(long long i = 0;i < env.thread_state;i+=per_chunk * env.chunk_state)
        {
            startIdx = task.tid * env.thread_state + i;
            _mpi_one_gate_inner(task,g,startIdx);
        }
    }
    else
    {
        long long total_chunk_per_thread = env.thread_state / env.chunk_state;
        long long per_chunk = min(total_chunk_per_thread,(long long)env.MPI_buffer_size);
        task.partner_using = {env.rank ^ (1 << (g->targs[1] - seg.N))};
        if(g->chunk_count)
        {
            long long startIdx;
            for(long long i = 0;i < env.thread_state;i+=per_chunk * env.chunk_state)
            {
                startIdx = task.tid * env.thread_state + i;
                _mpi_one_gate_inner(task,g,startIdx);
            }
        }
        else
        {
            if(isFile(g->targs[0]))
            {
                int partner_rank = env.rank ^ (1 << (g->targs[1] - seg.N));
                int thread_mask = 1 << (g->targs[0] - seg.middle - seg.chunk);
                int t0 = task.tid & (~thread_mask);
                int t1 = task.tid | thread_mask;
                bool thread_equal_chunk = (env.thread_state == env.chunk_state);
                int represent_thread = (env.rank < partner_rank)? t1 : t0;
                if(thread_equal_chunk)
                {
                    int recv_tag = represent_thread == t1? t0 : t1;
                    if(task.tid != represent_thread) return;
                    long long startIdx = represent_thread * env.thread_state;
                    MPI_Isend(&buffer[startIdx],env.thread_state,MPI_DOUBLE_COMPLEX,partner_rank,represent_thread,MPI_COMM_WORLD,&task.request1_send);
                    MPI_Irecv(&task.buffer2[0],env.thread_state,MPI_DOUBLE_COMPLEX,partner_rank,recv_tag,MPI_COMM_WORLD,&task.request1_recv);
                    MPI_Wait(&task.request1_recv,MPI_STATUS_IGNORE);
                    MPI_Wait(&task.request1_send,MPI_STATUS_IGNORE);
                    g->run_one_qubit_mpi_mem(buffer,task.buffer2,startIdx,0,env.thread_state);
                }
                else
                {
                    long long total_chunk_per_thread = (env.thread_state / env.chunk_state) >> 1;
                    long long per_chunk = min(total_chunk_per_thread,(long long)env.MPI_buffer_size);
                    for(long long i = 0;i < env.thread_state >> 1;i+=per_chunk * env.chunk_state)
                    {
                        long long startIdx = i + represent_thread * env.thread_state + (task.tid != represent_thread) * (env.thread_state >> 1);
                        MPI_Isend(&buffer[startIdx],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request1_send);
                        MPI_Irecv(&task.buffer2[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request1_recv);
                        MPI_Wait(&task.request1_recv,MPI_STATUS_IGNORE);
                        MPI_Wait(&task.request1_send,MPI_STATUS_IGNORE);
                        g->run_one_qubit_mpi_mem(buffer,task.buffer2,startIdx,0,per_chunk * env.chunk_state);
                    }
                }
            }
            else
            {
                int partner_rank = env.rank ^ (1 << (g->targs[1] - seg.N));
                // long long block_chunk = env.qubit_offset[g->targs[0]] / env.chunk_state;
                // long long per_chunk = min((long long)env.MPI_buffer_size,block_chunk);
                // for(long long i = 0;i < env.thread_state;i+=env.qubit_offset[g->targs[0] + 1])
                // {
                //     long long base = task.tid * env.thread_state + i + (env.rank < partner_rank) * env.qubit_offset[g->targs[0]];
                //     for(int j = 0;j < block_chunk;j+=per_chunk)
                //     {
                //         MPI_Isend(&buffer[base + j * env.chunk_state],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request1_send);
                //         MPI_Irecv(&task.buffer2[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request1_recv);
                //         MPI_Wait(&task.request1_recv,MPI_STATUS_IGNORE);
                //         MPI_Wait(&task.request1_send,MPI_STATUS_IGNORE);
                //         g->run_one_qubit_mpi_mem(buffer,task.buffer2,base + j * env.chunk_state,0,per_chunk * env.chunk_state);
                //     }
                // }
                // algo faster
                long long total_chunk_per_thread = (env.thread_state >> 1) / env.chunk_state;
                long long per_chunk = min(total_chunk_per_thread,(long long)env.MPI_buffer_size);
                vector<complex<double>>buffer_tmp(per_chunk * env.chunk_state);
                int chunk_cnt = 0;
                unordered_map<int,int>m;
                for(long long i = 0;i < env.thread_state;i+=env.qubit_offset[g->targs[0] + 1])
                {
                    for(long long j = 0;j < env.qubit_offset[g->targs[0]];j+=env.chunk_state)
                    {
                        long long startIdx = task.tid * env.thread_state + i + j + (env.rank < partner_rank) * env.qubit_offset[g->targs[0]];
                        copy(buffer.begin() + startIdx,buffer.begin() + startIdx + env.chunk_state,buffer_tmp.begin() + chunk_cnt * env.chunk_state);
                        m[chunk_cnt++] = startIdx;
                        if(chunk_cnt == per_chunk)
                        {
                            MPI_Isend(&buffer_tmp[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request1_send);
                            MPI_Irecv(&task.buffer2[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request1_recv);
                            MPI_Wait(&task.request1_recv,MPI_STATUS_IGNORE);
                            MPI_Wait(&task.request1_send,MPI_STATUS_IGNORE);
                            if(env.rank < partner_rank)
                            {
                                g->run_one_qubit_mpi_mem(buffer_tmp,task.buffer2,0,0,per_chunk * env.chunk_state);
                            }
                            else
                            {
                                g->run_one_qubit_mpi_mem(task.buffer2,buffer_tmp,0,0,per_chunk * env.chunk_state);
                            }
                            chunk_cnt = 0;
                            for(int q = 0;q < per_chunk;q++)
                            {
                                copy(buffer_tmp.begin() + q * env.chunk_state,buffer_tmp.begin() + (q + 1) * env.chunk_state,buffer.begin() + m[q]);
                            }
                        }
                    }
                }
            }
        }
    }
}
void MEM_Runner::MPI_Swap_2_2(thread_MEM_task &task,Gate* &g)
{
    int rank_mask_2 = 1 << (g->targs[2] - seg.N);
    int rank_mask_3 = 1 << (g->targs[3] - seg.N);
    int rank0 = env.rank & (~rank_mask_2) & (~rank_mask_3);
    int rank1 = rank0 | rank_mask_2;
    int rank2 = rank0 | rank_mask_3;
    int rank3 = rank0 | rank_mask_2 | rank_mask_3;
    unordered_map<int,int>rank_mapping{{rank0,0},{rank1,1},{rank2,2},{rank3,3}};
    int rank_order = rank_mapping[env.rank];
    long long per_chunk = min(env.thread_state / env.chunk_state,(long long)env.MPI_buffer_size);
    vector<MPI_Request*>request_send = {&task.request1_send,&task.request2_send,&task.request3_send};
    vector<MPI_Request*>request_recv = {&task.request1_recv,&task.request2_recv,&task.request3_recv};
    vector<vector<complex<double>>*>buffer_recv_using = {&task.buffer2,&task.buffer3,&task.buffer4};
    vector<complex<double>>*buffer1_ptr = &buffer;
    vector<complex<double>>*buffer2_ptr = &task.buffer2;
    vector<complex<double>>*buffer3_ptr = &task.buffer3;
    vector<complex<double>>*buffer4_ptr = &task.buffer4;
    long long pos0 = 0;
    long long pos1 = 0;
    long long pos2 = 0;
    long long pos3 = 0;
    long long *local_pos;
    if(rank_order == 0)
    {
        task.partner_using = {rank1,rank2,rank3};
        local_pos = &pos0;
    }
    else if(rank_order == 1)
    {
        task.partner_using = {rank0,rank2,rank3};
        buffer1_ptr = &task.buffer2;
        buffer2_ptr = &buffer;
        local_pos = &pos1;
    }
    else if(rank_order == 2)
    {
        task.partner_using = {rank0,rank1,rank3};
        buffer1_ptr = &task.buffer2;
        buffer2_ptr = &task.buffer3;
        buffer3_ptr = &buffer;
        local_pos = &pos2;
    }
    else
    {
        task.partner_using = {rank0,rank1,rank2};
        buffer1_ptr = &task.buffer2;
        buffer2_ptr = &task.buffer3;
        buffer3_ptr = &task.buffer4;
        buffer4_ptr = &buffer;
        local_pos = &pos3;
    }
    for(int i = task.tid * env.thread_state;i < ((task.tid + 1) * env.thread_state);i+=per_chunk * env.chunk_state)
    {
        for(int j = 0;j < 3;j++)
        {
            MPI_Isend(&buffer[i],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[j],task.tid,MPI_COMM_WORLD,&(*request_send[j]));
            MPI_Irecv(&(*buffer_recv_using[j])[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[j],task.tid,MPI_COMM_WORLD,&(*request_recv[j]));
        }
        for(int j = 0;j < 3;j++)
        {
            MPI_Wait(&(*request_send[j]),MPI_STATUS_IGNORE);
            MPI_Wait(&(*request_recv[j]),MPI_STATUS_IGNORE);
        }
        *local_pos = i;
        g->run_mpi_vswap2_2_mem(*buffer1_ptr,*buffer2_ptr,*buffer3_ptr,*buffer4_ptr,pos0,pos1,pos2,pos3,per_chunk * env.chunk_state);
    }
}