#include <omp.h>
#include <unistd.h> //pread, pwrite

#include <complex>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>
#include <algorithm>
#include "circuit.hpp"
#include <functional>
#include "MPI_Runner.hpp"



void MPI_Runner::MPI_Swap(thread_MPI_task &task,Gate * &g)
{
    long long mpi_targ_mask_0 = 1 << (g->targs[0] - seg.N);
    long long mpi_targ_mask_1 = 1 << (g->targs[1] - seg.N);
    long long loop_bound = env.thread_size;
    long long loop_stride;
    bool firstround = true;
    task.partner_using = {int((long long)env.rank ^ mpi_targ_mask_1)};
    if(isMpi(g->targs[0]))
    {
        task.fd_using = {env.fd_arr[task.tid],env.fd_arr[task.tid]};
        task.fd_offset_using = {0};
        loop_stride = env.chunk_size;
        if(((env.rank >> (g->targs[1] - seg.N)) & 1) == ((env.rank >> (g->targs[0] - seg.N)) & 1)) return;
        task.partner_using = {int(((long long)env.rank) ^ mpi_targ_mask_0 ^ mpi_targ_mask_1)};
    }
    else if(isFile(g->targs[0]))
    {
        long long file_mask = 1 << (g->targs[0] - seg.middle - seg.chunk);
        int fd0 = task.tid & (~file_mask);
        int fd1 = task.tid | file_mask;
        task.fd_using = {env.fd_arr[fd0],env.fd_arr[fd1]};
        task.fd_offset_using = {0};
        bool thread_equal_chunk = (env.thread_size == env.chunk_size);
        loop_stride = (thread_equal_chunk)? env.chunk_size : env.chunk_size << 1;
        if(thread_equal_chunk && task.tid == fd1) return;
    }
    else if(isMiddle(g->targs[0]))
    {
        task.fd_using = {env.fd_arr[task.tid],env.fd_arr[task.tid]};
        loop_stride = env.qubit_size[g->targs[0] + 1];
    }
    else
    {
        task.fd_using = {env.fd_arr[task.tid]};
        task.fd_offset_using = {0};
        loop_stride = env.chunk_size << 1;
    }
    for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
    {
        if(isMiddle(g->targs[0]))
        {
            task.fd_offset_using = (env.rank < task.partner_using[0])? vector<long long>{cur_offset + env.qubit_size[g->targs[0]]} : vector<long long>{cur_offset};
            for(long long j = 0;j < env.qubit_size[g->targs[0]];j += env.chunk_size)
            {
                _thread_MPI_swap(task,g,firstround);
                update_offset(task,env.chunk_size);
            }
        }
        else if(isChunk(g->targs[0]))
        {
            if(env.thread_state == env.chunk_state && env.rank > task.partner_using[0])
                _thread_no_exec_MPI(task,1);
            else
            {
                for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
                {
                    _thread_read1_recv1(task,g);
                    update_offset(task,loop_stride);
                }
            }
        }
        else
        {
            _thread_MPI_swap(task,g,firstround);
            update_offset(task,loop_stride);
        }
    }
}
void MPI_Runner::_thread_MPI_swap(thread_MPI_task &task,Gate * &g,bool &firstround)
{
    int partner_rank = task.partner_using[0];
    int member_th = env.rank > partner_rank? 1 : 0;
    int thread_member_th = isMpi(g->targs[1]) && isFile(g->targs[0]) && (task.tid & (1 << (g->targs[0] - seg.middle - seg.chunk)))? 1 : 0;//For all thread drive such as D F,otherwise it is 0
    MPI_Irecv(&task.buffer[(!member_th) * env.chunk_state], env.chunk_state, MPI_DOUBLE_COMPLEX, partner_rank,task.tid, MPI_COMM_WORLD,&task.request[!member_th]);
    if(!firstround)
        MPI_Wait(&task.request[member_th],MPI_STATUS_IGNORE);
    firstround = false;
    if(pread(task.fd_using[!member_th],&task.buffer[member_th * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + thread_member_th * env.chunk_size));
    MPI_Isend(&task.buffer[member_th * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request[member_th]);
    MPI_Wait(&task.request[!member_th],MPI_STATUS_IGNORE);
    if(pwrite(task.fd_using[!member_th],&task.buffer[(!member_th) * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + thread_member_th * env.chunk_size));
}
void MPI_Runner::MPI_gate_scheduler(thread_MPI_task &task,Gate * &g)
{
    long long mpi_targ_mask_0 = 1 << (g->targs[1] - seg.N);
    long long mpi_targ_mask_2 = 1 << (g->targs[0] - seg.N);
    long long loop_bound = env.thread_size;
    long long loop_stride;
    long long threshold;
    int rank0;
    int rank1;
    int rank2;
    int rank3;
    if(g->targs.size() == 1)
    {
        task.fd_using = {env.fd_arr[task.tid]};
        task.fd_offset_using = {0};
        task.gate_buffer_using = {1};
        if(g->name == "Z_Gate" || g->name == "Phase_Gate")
        {
            if(env.rank != (env.rank | mpi_targ_mask_2)) return;
            inner_all_thread(task,g,env.thread_size,1);
            return;
        }
        task.partner_using = {int(((long long)env.rank) ^ mpi_targ_mask_2)};
        loop_stride = env.chunk_size << 1;
        if(env.thread_state == env.chunk_state && env.rank > task.partner_using[0])
            _thread_no_exec_MPI(task,1);
        else
        {
            for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
            {
                _thread_read1_recv1(task,g);
                update_offset(task,loop_stride);
            }
        }
    }
    else
    {
        if(isMpi(g->targs[0]))
        {
            rank0 = env.rank & (~(mpi_targ_mask_0 | mpi_targ_mask_2));
            rank1 = rank0 | mpi_targ_mask_2;
            rank2 = rank0 | mpi_targ_mask_0;
            rank3 = rank0 | mpi_targ_mask_0 | mpi_targ_mask_2;
            if(g->name == "CPhase_Gate")
            {
                if(env.rank != rank3) return;
                task.fd_using = {env.fd_arr[task.tid]};
                task.fd_offset_using = {0};
                task.gate_buffer_using = {3};
                inner_all_thread(task,g,env.thread_size,1);
                return;
            }
            task.partner_using = {rank0,rank1,rank2,rank3};
            task.fd_using = {env.fd_arr[task.tid],env.fd_arr[task.tid],env.fd_arr[task.tid],env.fd_arr[task.tid]};
            task.fd_offset_using = {0,env.chunk_size,env.chunk_size << 1,env.chunk_size * 3};
            threshold = min(env.thread_state / env.chunk_state,(long long)4);
            loop_stride = (1 << int(log2(threshold))) * env.chunk_size;
            if(env.rank > task.partner_using[threshold - 1])
            {
                _thread_no_exec_MPI(task,1 << int(log2(threshold)));
            }
            else
            {
                for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
                {
                    _thread_read1_recv3(task,g,1 << int(log2(threshold)));
                    update_offset(task,loop_stride);
                }
            }
        }
        else if(isFile(g->targs[0]))
        {
            long long file_mask = 1 << (g->targs[0] - seg.middle - seg.chunk);
            rank0 = env.rank ^ mpi_targ_mask_0;
            rank1 = env.rank ^ mpi_targ_mask_0;
            if(g->name == "CPhase_Gate")
            {
                if(env.rank < rank1) return;
                int fd1 = task.tid | file_mask;
                bool thread_equal_chunk = (env.thread_size == env.chunk_size);
                if(thread_equal_chunk && task.tid == fd1) return;
                task.fd_using = {env.fd_arr[fd1]};
                task.fd_offset_using = (task.tid != fd1)? vector<long long> {0} : vector<long long> {env.chunk_size};
                task.gate_buffer_using = {3};
                inner_all_thread(task,g,env.thread_size,1);
                return;
            }
            task.partner_using = {rank0,rank1};
            task.fd_using = {env.fd_arr[task.tid & (~file_mask)],env.fd_arr[task.tid | file_mask]};
            task.fd_offset_using = {0,0};
            threshold = min(env.thread_state / env.chunk_state,(long long)4);
            loop_stride = (threshold > 2)? (env.chunk_size << 2) : (env.chunk_size << 1);
            int num_worker = (threshold == 1)? 1 : 2;
            if((threshold <= 2) && (task.tid & ((1 << (g->targs[0] - seg.middle - seg.chunk)))))
                return;
            if(threshold == 1 && env.rank > task.partner_using[0])
            {
                _thread_no_exec_MPI(task,2);
                return;
            }
            for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
            {
                _thread_read2_recv2(task,g,num_worker,env.chunk_size);
                update_offset(task,loop_stride);
            }
        }
        else if(isMiddle(g->targs[0]))
        {
            rank0 = env.rank ^ mpi_targ_mask_0;
            rank1 = env.rank ^ mpi_targ_mask_0;
            if(g->name == "CPhase_Gate")
            {
                if(env.rank < rank1) return;
                task.fd_using = {env.fd_arr[task.tid]};
                task.gate_buffer_using = {3};
                for(long long cur_offset = 0;cur_offset < env.thread_size;cur_offset += env.qubit_size[g->targs[0] + 1])
                {
                    task.fd_offset_using = {cur_offset + env.qubit_size[g->targs[0]]};
                    inner_all_thread(task,g,env.qubit_size[g->targs[0]],1);
                }
                return;
            }
            task.partner_using = {rank0,rank1};
            task.fd_using = {env.fd_arr[task.tid],env.fd_arr[task.tid]};
            loop_bound = (env.thread_state == env.qubit_offset[g->targs[0] + 1])? env.qubit_size[g->targs[0]] : env.thread_size;
            loop_stride = (env.thread_state == env.qubit_offset[g->targs[0] + 1])? (env.chunk_size << 1) : (env.qubit_size[g->targs[0] + 1] << 1);
            if(env.thread_state == env.qubit_offset[g->targs[0] + 1])
            {
                task.fd_offset_using = {0,env.qubit_size[g->targs[0]]};
                if(env.chunk_state == env.qubit_offset[g->targs[0]])
                {
                    if(env.rank > task.partner_using[0])
                        _thread_no_exec_MPI(task,2);
                    else
                        _thread_read2_recv2(task,g,1,0);
                }
                else
                {
                    for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
                    {   
                        _thread_read2_recv2(task,g,2,loop_stride >> 1);
                        update_offset(task,loop_stride);
                    }
                }
            }
            else
            {
                for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
                {
                    task.fd_offset_using = {cur_offset,cur_offset + env.qubit_size[g->targs[0]]};
                    for(long long j = 0;j < env.qubit_size[g->targs[0]];j += env.chunk_size)
                    {
                        _thread_read2_recv2(task,g,2,loop_stride >> 1);
                        update_offset(task,env.chunk_size);
                    }
                }
            }
        }
        else
        {
            rank0 = env.rank ^ mpi_targ_mask_0;
            if(g->name == "CPhase_Gate")
            {
                if(env.rank < rank0) return;
                task.fd_using = {env.fd_arr[task.tid]};
                task.fd_offset_using = {0};
                task.gate_buffer_using = {1};
                inner_all_thread(task,g,env.thread_size,1);
                return;
            }
            task.partner_using = {rank0};
            task.fd_using = {env.fd_arr[task.tid]};
            task.fd_offset_using = {0};
            loop_stride = env.chunk_size << 1;
            if(env.thread_state == env.chunk_state && env.rank > task.partner_using[0])
                _thread_no_exec_MPI(task,1);
            else
            {
                for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
                {
                    _thread_read1_recv1(task,g);
                    update_offset(task,loop_stride);
                }
            }
        }
    }
}
void MPI_Runner::_thread_read1_recv1(thread_MPI_task &task,Gate * &g)
{
    int num_worker = (env.chunk_state == env.thread_state)? 1 : 2;
    int member_th = env.rank > task.partner_using[0]? 1 : 0;
    int partner_rank = task.partner_using[0];

    MPI_Irecv(&task.buffer[(!member_th) * env.chunk_state], env.chunk_state, MPI_DOUBLE_COMPLEX, partner_rank,task.tid, MPI_COMM_WORLD,&task.request[!member_th]);
    if(num_worker == 2)
    {
        if(pread(task.fd_using[0],&task.buffer[env.chunk_state << 1],env.chunk_size,task.fd_offset_using[0] + (!member_th) * env.chunk_size));
        MPI_Isend(&task.buffer[env.chunk_state << 1],env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request[2]);
    }
    if(pread(task.fd_using[0],&task.buffer[member_th * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + member_th * env.chunk_size));
    MPI_Wait(&task.request[!member_th],MPI_STATUS_IGNORE);
    g->run(task.buffer);
    MPI_Isend(&task.buffer[(!member_th) * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request[!member_th]);
    if(pwrite(task.fd_using[0],&task.buffer[member_th * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + member_th * env.chunk_size));
    if(num_worker == 2)
    {
        MPI_Recv(&task.buffer[env.chunk_state << 1],env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        if(pwrite(task.fd_using[0],&task.buffer[env.chunk_state << 1],env.chunk_size,task.fd_offset_using[0] + (!member_th) * env.chunk_size));
    }
}
void MPI_Runner::_thread_read2_recv2(thread_MPI_task &task,Gate * &g,int num_worker,long long stride)
{
    int member_th = env.rank > task.partner_using[0]? 1 : 0;//for rank
    int thread_member_th = isMpi(g->targs[1]) && isFile(g->targs[0]) && (task.tid & (1 << (g->targs[0] - seg.middle - seg.chunk)))? 1 : 0;//For all thread drive such as D F,otherwise it is 0
    MPI_Irecv(&task.buffer[((!member_th) << 1) * env.chunk_state], env.chunk_state, MPI_DOUBLE_COMPLEX, task.partner_using[0],task.tid, MPI_COMM_WORLD,&task.request[(!member_th) << 1]);
    MPI_Irecv(&task.buffer[(((!member_th) << 1) | 1) * env.chunk_state], env.chunk_state, MPI_DOUBLE_COMPLEX, task.partner_using[1],task.tid, MPI_COMM_WORLD,&task.request[((!member_th) << 1) | 1]);
    if(num_worker == 2)
    {
        if(pread(task.fd_using[0],&task.buffer[env.chunk_state << 2],env.chunk_size,task.fd_offset_using[0] + (!member_th) * stride + thread_member_th * 2 * stride));
        MPI_Isend(&task.buffer[env.chunk_state << 2],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,&task.request[4]);
        if(pread(task.fd_using[1],&task.buffer[5 * env.chunk_state],env.chunk_size,task.fd_offset_using[1] + (!member_th) * stride + thread_member_th * 2 * stride));
        MPI_Isend(&task.buffer[5 * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[1],task.tid,MPI_COMM_WORLD,&task.request[5]);
    }
    if(pread(task.fd_using[0],&task.buffer[(member_th << 1) * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + member_th * stride + thread_member_th * 2 * stride));
    if(pread(task.fd_using[1],&task.buffer[((member_th << 1) | 1) * env.chunk_state],env.chunk_size,task.fd_offset_using[1] + member_th * stride + thread_member_th * 2 * stride));
    MPI_Waitall(2,&task.request[(!member_th) << 1],MPI_STATUS_IGNORE);
    g->run(task.buffer);
    MPI_Isend(&task.buffer[((!member_th) << 1) * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,&task.request[(!member_th) << 1]);
    MPI_Isend(&task.buffer[(((!member_th) << 1) | 1) * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[1],task.tid,MPI_COMM_WORLD,&task.request[((!member_th) << 1) | 1]);
    if(num_worker == 2)
    {
        MPI_Irecv(&task.buffer[6 * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,&task.request[6]);
        MPI_Irecv(&task.buffer[7 * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[1],task.tid,MPI_COMM_WORLD,&task.request[7]);
    }
    if(pwrite(task.fd_using[0],&task.buffer[(member_th << 1) * env.chunk_state],env.chunk_size, task.fd_offset_using[0] + member_th * stride + thread_member_th * 2 * stride));
    if(pwrite(task.fd_using[1],&task.buffer[((member_th << 1) | 1) * env.chunk_state],env.chunk_size,task.fd_offset_using[1] + member_th * stride + thread_member_th * 2 * stride));

    if(num_worker == 2)
    {
        vector<bool>used(2,false);
        int next_buffer_idx = Get_Next_Undone_Buffer_index(task.request,used,2,6);
        while(next_buffer_idx != -1)
        {
            if(next_buffer_idx == 6 && pwrite(task.fd_using[0],&task.buffer[6 * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + (!member_th) * stride + thread_member_th * 2 * stride));
            if(next_buffer_idx == 7 && pwrite(task.fd_using[1],&task.buffer[7 * env.chunk_state],env.chunk_size,task.fd_offset_using[1] + (!member_th) * stride + thread_member_th * 2 * stride));
            next_buffer_idx = Get_Next_Undone_Buffer_index(task.request,used,2,6);
        }
    }
}
void MPI_Runner::_thread_read1_recv3(thread_MPI_task &task,Gate * &g,int num_worker)
{
    int member_th;
    for(int i = 0;i < 4;i++)
    {
        if(env.rank == task.partner_using[i])
            member_th = i;
        else
            MPI_Irecv(&task.buffer[i * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[i],task.tid,MPI_COMM_WORLD,&task.request[i]);
    }
    for(int i = 0;i < num_worker;i++)
    {
        if(i != member_th)
        {
            if(pread(task.fd_using[i],&task.buffer[(i + 4) * env.chunk_state],env.chunk_size,task.fd_offset_using[i]));
            MPI_Isend(&task.buffer[(i + 4) * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[i],task.tid,MPI_COMM_WORLD,&task.request[i + 4]);
        }
    }
    if(pread(task.fd_using[member_th],&task.buffer[member_th * env.chunk_state],env.chunk_size,task.fd_offset_using[member_th]));
    for(int i = 0;i < 4;i++)
    {
        if(i != member_th)
            MPI_Wait(&task.request[i],MPI_STATUS_IGNORE);
    }
    g->run(task.buffer);
    for(int i = 0;i < 4;i++)
        if(i != member_th) MPI_Isend(&task.buffer[i * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[i],task.tid,MPI_COMM_WORLD,&task.request[i]);
    if(pwrite(task.fd_using[member_th],&task.buffer[member_th * env.chunk_state],env.chunk_size,task.fd_offset_using[member_th]));
    unordered_map<int,int>m;
    for(int i = 0,rd_off = 8;i < num_worker;i++)
    {
        if(i != member_th)
        {
            MPI_Irecv(&task.buffer[rd_off * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[i],task.tid,MPI_COMM_WORLD,&task.request[rd_off]);
            m[rd_off++] = i;
        }
    }
    vector<bool>used(num_worker - 1,false);
    int next_buffer_idx = Get_Next_Undone_Buffer_index(task.request,used,num_worker - 1,8);
    while(next_buffer_idx != -1)
    {
        if(pwrite(task.fd_using[m[next_buffer_idx]],&task.buffer[next_buffer_idx * env.chunk_state],env.chunk_size,task.fd_offset_using[m[next_buffer_idx]]));
        next_buffer_idx = Get_Next_Undone_Buffer_index(task.request,used,num_worker - 1,8);
    }
}
void MPI_Runner::_thread_no_exec_MPI(thread_MPI_task &task,int round)
{
    vector<bool>used(round,false);
    for(int i = 0;i < round;i++)
    {
        MPI_Irecv(&task.buffer[(round + i) * env.chunk_state],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[i],task.tid,MPI_COMM_WORLD,&task.request[round + i]);
    }
    for(int i = 0;i < round;i++)
    {
        if(pread(task.fd_using[i], &task.buffer[i * env.chunk_state], env.chunk_size, task.fd_offset_using[i]));
        MPI_Isend(&task.buffer[i * env.chunk_state], env.chunk_state, MPI_DOUBLE_COMPLEX, task.partner_using[i], task.tid, MPI_COMM_WORLD,&task.request[i]);
    }
    int next_buffer_idx = Get_Next_Undone_Buffer_index(task.request,used,round,round);
    while(next_buffer_idx != -1)
    {
        if(pwrite(task.fd_using[next_buffer_idx - round], &task.buffer[next_buffer_idx * env.chunk_state], env.chunk_size, task.fd_offset_using[next_buffer_idx - round]));
        next_buffer_idx = Get_Next_Undone_Buffer_index(task.request,used,round,round);
    }
}
int MPI_Runner::Get_Next_Undone_Buffer_index(vector<MPI_Request> &request,vector<bool>&used,int round,int start)
{
    int res = -1;
    int fl;
    int done_cnt = 0;
    while(done_cnt != round)
    {
        done_cnt = 0;
        for(int i = 0;i < round;i++)
        {
            MPI_Test(&request[start + i],&fl,MPI_STATUS_IGNORE);
            if(fl)
            {
                if(!used[i])
                {
                    res = start + i;
                    used[i] = true;
                    return res;
                }
                else
                {
                    done_cnt++;
                }
            }
        }
    }
    return res;
}