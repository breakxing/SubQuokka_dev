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
        int rank_mask_0 = 1 << (g->targs[0] - seg.N);
        int partner_rank = env.rank ^ rank_mask_0;
        int total_chunk_per_thread = env.thread_state / env.chunk_state;
        int per_chunk = min(total_chunk_per_thread,env.MPI_buffer_size);
        task.partner_using = {partner_rank};
        long long startIdx;
        for(long long cur_offset = 0;cur_offset < env.thread_state;cur_offset+=per_chunk * env.chunk_state)
        {
            startIdx = task.tid * env.thread_state + cur_offset;
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
            int per_chunk = min(env.thread_state / env.chunk_state,(long long)(env.MPI_buffer_size >> 2));
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
            for(long long cur_offset = 0;cur_offset <env.thread_state;cur_offset+=per_chunk * env.chunk_state)
            {
                *local_pos = task.tid * env.thread_state + cur_offset;
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
            int total_chunk_per_thread = env.thread_state / env.chunk_state / 2;
            int per_chunk = min(total_chunk_per_thread,env.MPI_buffer_size);
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
                buffer1_ptr = &buffer;
                buffer2_ptr = &buffer;
                buffer3_ptr = &task.buffer2;
                buffer4_ptr = &task.buffer3;
                if(env.rank > partner_rank)
                {
                    swap(buffer1_ptr,buffer3_ptr);
                    swap(buffer2_ptr,buffer4_ptr);
                }
                for(long long cur_offset = 0;cur_offset < env.thread_state >> 1;cur_offset+=per_chunk * env.chunk_state)
                {
                    pos0 = cur_offset + startIdx0;
                    pos1 = cur_offset + startIdx1;
                    pos2 = 0;
                    pos3 = 0;
                    if(env.rank > partner_rank)
                    {
                        swap(pos0,pos2);
                        swap(pos1,pos3);
                    }
                    MPI_Isend(&buffer[cur_offset + startIdx0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request1_send);
                    MPI_Isend(&buffer[cur_offset + startIdx1],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request2_send);
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
            int total_chunk_per_thread = (env.thread_state >> 1) / env.chunk_state;
            int per_chunk = min(total_chunk_per_thread,env.MPI_buffer_size >> 1);
            vector<complex<double>>buffer_tmp0(per_chunk * env.chunk_state);
            vector<complex<double>>buffer_tmp1(per_chunk * env.chunk_state);
            vector<complex<double>>*buffer1_ptr = &buffer_tmp0;
            vector<complex<double>>*buffer2_ptr = &buffer_tmp1;
            vector<complex<double>>*buffer3_ptr = &task.buffer2;
            vector<complex<double>>*buffer4_ptr = &task.buffer3;
            if(env.rank > partner_rank)
            {
                swap(buffer1_ptr,buffer3_ptr);
                swap(buffer2_ptr,buffer4_ptr);
            }
            stack<long long>st0;
            stack<long long>st1;
            for(long long cur_offset = 0;cur_offset < env.thread_state;cur_offset+=env.qubit_offset[g->targs[0]] << 1)
            {
                for(long long j = 0;j < env.qubit_offset[g->targs[0]];j+=env.chunk_state)
                {
                    long long startIdx0 = task.tid * env.thread_state + cur_offset + j;
                    long long startIdx1 = task.tid * env.thread_state + cur_offset + j + env.qubit_offset[g->targs[0]];
                    copy(buffer.begin() + startIdx0,buffer.begin() + startIdx0 + env.chunk_state,buffer_tmp0.begin() + st0.size() * env.chunk_state);
                    copy(buffer.begin() + startIdx1,buffer.begin() + startIdx1 + env.chunk_state,buffer_tmp1.begin() + st1.size() * env.chunk_state);
                    st0.push(startIdx0);
                    st1.push(startIdx1);
                    if(st0.size() == (unsigned long int) per_chunk)
                    {
                        MPI_Isend(&buffer_tmp0[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request1_send);
                        MPI_Isend(&buffer_tmp1[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request2_send);
                        MPI_Irecv(&task.buffer2[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request1_recv);
                        MPI_Irecv(&task.buffer3[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request2_recv);
                        MPI_Wait(&task.request1_send,MPI_STATUS_IGNORE);
                        MPI_Wait(&task.request2_send,MPI_STATUS_IGNORE);
                        MPI_Wait(&task.request1_recv,MPI_STATUS_IGNORE);
                        MPI_Wait(&task.request2_recv,MPI_STATUS_IGNORE);
                        g->run_mem_u2(*buffer1_ptr,*buffer2_ptr,*buffer3_ptr,*buffer4_ptr,0,0,0,0,per_chunk * env.chunk_state);
                        MPI_Swap_restore(buffer_tmp0,st0,env.chunk_state);
                        MPI_Swap_restore(buffer_tmp1,st1,env.chunk_state);
                    }
                }
            }
        }
        else
        {
            int total_chunk_per_thread = env.thread_state / env.chunk_state;
            int per_chunk = min(total_chunk_per_thread,env.MPI_buffer_size);
            task.partner_using = {env.rank ^ (1 << (g->targs[1] - seg.N))};
            for(long long cur_offset = 0;cur_offset < env.thread_state;cur_offset+=per_chunk * env.chunk_state)
            {
                _mpi_one_gate_inner(task,g,cur_offset + task.tid * env.thread_state);
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
    for(long long cur_offset = 0;cur_offset < (1LL << seg.N);cur_offset+=env.chunk_state)
    {
        g->run_one_qubit_mpi_mem_diagonal(buffer,cur_offset,rank_order,env.chunk_state);
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
        for(long long cur_offset = 0;cur_offset < (1LL << seg.N);cur_offset+=env.chunk_state)
        {
            g->run_one_qubit_mpi_mem_diagonal(buffer,cur_offset,rank_order,env.chunk_state);
        }
    }
    else
    {
        long long stride = env.qubit_offset[g->targs[0]] << 1;
        if(g->name == "CPhase_Gate_MEM")
        {
            if(env.rank < (1 <<(g->targs[1] - seg.N))) return;
            #pragma omp parallel for schedule(static)
            for(long long cur_offset = 0;cur_offset < (1LL << seg.N);cur_offset+=stride)
            {
                MPI_CPhase(g,cur_offset);
            }
        }
        else
        {
            #pragma omp parallel for schedule(static)
            for(long long cur_offset = 0;cur_offset < (1LL << seg.N);cur_offset+=stride)
            {
                MPI_RZZ(g,cur_offset);
            }

        }
    }
}
void MEM_Runner::_mpi_one_gate_inner(thread_MEM_task &task,Gate* &g,long long startIdx)
{
    int partner_rank = task.partner_using[0];
    int per_chunk = min(env.thread_state / env.chunk_state,(long long)env.MPI_buffer_size);
    vector<complex<double>>*buffer1_ptr = &buffer;
    vector<complex<double>>*buffer2_ptr = &task.buffer2;
    long long pos1 = startIdx;
    long long pos2 = 0;
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
        int total_chunk_per_thread = env.thread_state / env.chunk_state;
        int per_chunk = min(total_chunk_per_thread,env.MPI_buffer_size);
        task.partner_using = {env.rank == rank1? rank2 : rank1};
        long long startIdx;
        for(long long cur_offset = 0;cur_offset < env.thread_state;cur_offset+=per_chunk * env.chunk_state)
        {
            startIdx = task.tid * env.thread_state + cur_offset;
            _mpi_one_gate_inner(task,g,startIdx);
        }
    }
    else
    {
        int total_chunk_per_thread = env.thread_state / env.chunk_state;
        int per_chunk = min(total_chunk_per_thread,env.MPI_buffer_size);
        task.partner_using = {env.rank ^ (1 << (g->targs[1] - seg.N))};
        if(g->chunk_count)
        {
            long long total_state_per_thread = env.thread_state >> 1;
            long long per_state = min(total_state_per_thread,env.MPI_buffer_size * env.chunk_state);
            long long one_time_state = env.qubit_offset[g->targs[0]];
            vector<complex<double>>buffer_tmp(per_state);
            stack<long long>st;
            for(long long cur_offset = 0;cur_offset < env.thread_state;cur_offset+=env.qubit_offset[g->targs[0]] << 1)
            {
                long long startIdx = task.tid * env.thread_state + cur_offset + (env.rank < task.partner_using[0]) * env.qubit_offset[g->targs[0]];
                copy(buffer.begin() + startIdx,buffer.begin() + startIdx + one_time_state,buffer_tmp.begin() + st.size() * one_time_state);
                st.push(startIdx);
                if(st.size() * one_time_state == (unsigned long int)per_state)
                {
                    MPI_Isend(&buffer_tmp[0],per_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,&task.request1_send);
                    MPI_Irecv(&task.buffer2[0],per_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,&task.request1_recv);
                    MPI_Wait(&task.request1_recv,MPI_STATUS_IGNORE);
                    MPI_Wait(&task.request1_send,MPI_STATUS_IGNORE);
                    g->run_one_qubit_mpi_mem(buffer_tmp,task.buffer2,0,0,per_state);
                    MPI_Swap_restore(buffer_tmp,st,one_time_state);
                }
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
                    total_chunk_per_thread = (env.thread_state / env.chunk_state) >> 1;
                    per_chunk = min(total_chunk_per_thread,env.MPI_buffer_size);
                    for(long long cur_offset = 0;cur_offset < env.thread_state >> 1;cur_offset+=per_chunk * env.chunk_state)
                    {
                        long long startIdx = cur_offset + represent_thread * env.thread_state + (task.tid != represent_thread) * (env.thread_state >> 1);
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
                total_chunk_per_thread = (env.thread_state >> 1) / env.chunk_state;
                per_chunk = min(total_chunk_per_thread,env.MPI_buffer_size);
                vector<complex<double>>buffer_tmp(per_chunk * env.chunk_state);
                stack<long long>st;
                for(long long cur_offset = 0;cur_offset < env.thread_state;cur_offset+=env.qubit_offset[g->targs[0]] << 1)
                {
                    for(long long j = 0;j < env.qubit_offset[g->targs[0]];j+=env.chunk_state)
                    {
                        long long startIdx = task.tid * env.thread_state + cur_offset + j + (env.rank < task.partner_using[0]) * env.qubit_offset[g->targs[0]];
                        copy(buffer.begin() + startIdx,buffer.begin() + startIdx + env.chunk_state,buffer_tmp.begin() + st.size() * env.chunk_state);
                        st.push(startIdx);
                        if(st.size() == (unsigned long int)per_chunk)
                        {
                            MPI_Isend(&buffer_tmp[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,&task.request1_send);
                            MPI_Irecv(&task.buffer2[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,&task.request1_recv);
                            MPI_Wait(&task.request1_recv,MPI_STATUS_IGNORE);
                            MPI_Wait(&task.request1_send,MPI_STATUS_IGNORE);
                            if(env.rank < task.partner_using[0])
                            {
                                g->run_one_qubit_mpi_mem(buffer_tmp,task.buffer2,0,0,per_chunk * env.chunk_state);
                            }
                            else
                            {
                                g->run_one_qubit_mpi_mem(task.buffer2,buffer_tmp,0,0,per_chunk * env.chunk_state);
                            }
                            MPI_Swap_restore(buffer_tmp,st,env.chunk_state);
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
    int per_chunk = min(env.thread_state / env.chunk_state,(long long)env.MPI_buffer_size);
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
    for(long long cur_offset = task.tid * env.thread_state;cur_offset < ((task.tid + 1) * env.thread_state);cur_offset+=per_chunk * env.chunk_state)
    {
        for(int j = 0;j < 3;j++)
        {
            MPI_Isend(&buffer[cur_offset],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[j],task.tid,MPI_COMM_WORLD,&(*request_send[j]));
            MPI_Irecv(&(*buffer_recv_using[j])[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[j],task.tid,MPI_COMM_WORLD,&(*request_recv[j]));
        }
        for(int j = 0;j < 3;j++)
        {
            MPI_Wait(&(*request_send[j]),MPI_STATUS_IGNORE);
            MPI_Wait(&(*request_recv[j]),MPI_STATUS_IGNORE);
        }
        *local_pos = cur_offset;
        g->run_mpi_vswap2_2_mem(*buffer1_ptr,*buffer2_ptr,*buffer3_ptr,*buffer4_ptr,pos0,pos1,pos2,pos3,per_chunk * env.chunk_state);
    }
}
void MEM_Runner::MPI_CPhase(Gate* &g,long long idx)
{
    long long half_stride = env.qubit_offset[g->targs[0]];
    CPhase_Gate_MEM* cphaseGatePtr = static_cast<CPhase_Gate_MEM*>(g);
    complex<double> q0;
    idx += half_stride;
    for (long long cur_offset = 0; cur_offset < half_stride; cur_offset++)
    {
        buffer[idx] *= cphaseGatePtr->exp_iPhi;
        idx += 1;
    }
}
void MEM_Runner::MPI_RZZ(Gate* &g,long long idx)
{
    long long half_stride = env.qubit_offset[g->targs[0]];
    RZZ_Gate_MEM* rzzGatePtr = static_cast<RZZ_Gate_MEM*>(g);
    long long off0 = idx;
    long long off1 = idx + half_stride;
    int rank_order = env.rank == (env.rank | 1 << (g->targs[1] - seg.N));
    for (long long cur_offset = 0; cur_offset < half_stride; cur_offset++)
    {
        if(rank_order == 0)
        {
            buffer[off0] *= rzzGatePtr->exp_n_iPhi_2;
            buffer[off1] *= rzzGatePtr->exp_p_iPhi_2;
        }
        else
        {
            buffer[off0] *= rzzGatePtr->exp_p_iPhi_2;
            buffer[off1] *= rzzGatePtr->exp_n_iPhi_2;
        }
        off0 += 1;
        off1 += 1;
    }
}

void MEM_Runner::MPI_Swap_restore(vector<complex<double>>&buffer_tmp,stack<long long>&st,long long size)
{
    long long buffer_size = st.size() * size;
    buffer_size -= size;
    while(!st.empty())
    {
        long long off = st.top();st.pop();
        copy(buffer_tmp.begin() + buffer_size,buffer_tmp.begin() + buffer_size + size,buffer.begin() + off);
        buffer_size -= size;
    }
}