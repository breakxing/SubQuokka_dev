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
#include "IO_Runner.hpp"

bool IO_Runner::skip_read_write(Gate * & g,const int &idx)
{
    if((g->name == "Z_Gate" || g->name == "Phase_Gate") && (idx == 0) && (!isChunk(g->targs[0]) && (!isMpi(g->targs[0])))) return true;
    else if((g->name == "SWAP_Gate") && (idx == 0 || idx == 3) && (!isChunk(g->targs[0]))) return true;
    else if(g->name == "CPhase_Gate" && (idx != 1) && (!isChunk(g->targs[1]) && (isChunk(g->targs[0])) && (!isMpi(g->targs[1])))) return true;
    else if(g->name == "CPhase_Gate" && (idx != 3) && (!isChunk(g->targs[1]) && (!isMpi(g->targs[1])))) return true;
    return false;
}
void IO_Runner::all_thread_drive_scheduler(thread_IO_task &task,Gate * &g)
{
    vector<int> targ = g->targs;
    int tid = task.fd_table[0];
    long long func_loop_size;
    if(targ.size() == 1)
    {
        int file_mask = 1 << (g->targs[0] - seg.middle - seg.chunk);
        int fd0 = tid & (~file_mask);
        int fd1 = tid | file_mask;
        if(env.thread_state == env.chunk_state && tid == fd1) return;
        task.fd_using = {env.fd_arr[fd0],env.fd_arr[fd1]};
        func_loop_size = (env.thread_state == env.chunk_state)? env.thread_size : env.thread_size >> 1;
        task.fd_offset_using = (tid == fd0)? vector<long long>{0,0} : vector<long long>{func_loop_size,func_loop_size};
        task.gate_buffer_using = {0,1};
        inner_all_thread(task,g,func_loop_size,2);
    }
    else if(targ.size() == 2)
    {
        int file_mask_left = 1 << (g->targs[1] - seg.middle - seg.chunk);
        int file_mask_right = 1 << (g->targs[0] - seg.middle - seg.chunk);
        if(isFile(targ[0]))
        {
            int fd0 = tid & (~file_mask_left) & (~file_mask_right);
            int fd1 = fd0 | file_mask_right;
            int fd2 = fd0 | file_mask_left;
            int fd3 = fd1 | file_mask_left;
            task.fd_using = {env.fd_arr[fd0],env.fd_arr[fd1],env.fd_arr[fd2],env.fd_arr[fd3]};
            int shift = min(int(log2(env.thread_state / env.chunk_state)),2);
            func_loop_size = env.thread_size >> shift;
            vector<long long>off0 = {0}; off0.resize(4,off0[0]);
            vector<long long>off1 = {func_loop_size}; off1.resize(4,off1[0]);
            vector<long long>off2 = {func_loop_size << 1};off2.resize(4,off2[0]);
            vector<long long>off3 = {func_loop_size * 3};off3.resize(4,off3[0]);
            unordered_map<int,vector<long long>>m = {{fd0,off0},{fd1,off1},{fd2,off2},{fd3,off3}};
            task.fd_offset_using = m[tid];
            if((shift == 0) && (tid != fd0)) return;
            else if((shift == 1) && (tid != fd0 && tid != fd1)) return;
            task.gate_buffer_using = {0,1,2,3};
            inner_all_thread(task,g,func_loop_size,4);
        }
        else if(isMiddle(targ[0]))
        {
            int fd0 = tid & (~file_mask_left);
            int fd1 = tid | file_mask_left;
            task.fd_using = {env.fd_arr[fd0],env.fd_arr[fd0],env.fd_arr[fd1],env.fd_arr[fd1]};
            long long base1 = 0;
            long long base2 = env.qubit_size[targ[0]];
            bool cond1 = env.thread_state == (env.qubit_offset[targ[0]]) << 1;
            bool cond2 = env.chunk_state == env.qubit_offset[targ[0]];
            func_loop_size = (cond1 && !cond2)?(env.qubit_size[targ[0]] >> 1) : env.qubit_size[targ[0]];
            if(cond1 && cond2 && tid == fd1) return;
            task.fd_offset_using = {base1,base2,base1,base2};
            if(cond1 && !cond2 && tid == fd1) task.fd_offset_using = {base1 + (env.qubit_size[targ[0]] >> 1),base2 + (env.qubit_size[targ[0]] >> 1),base1 + (env.qubit_size[targ[0]] >> 1),base2 + (env.qubit_size[targ[0]] >> 1)};
            if(!cond1 && tid == fd1) task.fd_offset_using = {base1 + (env.thread_size >> 1),base2 + (env.thread_size >> 1),base1 + (env.thread_size >> 1),base2 + (env.thread_size >> 1)};
            task.gate_buffer_using = {0,1,2,3};
            if(cond1)
                inner_all_thread(task,g,func_loop_size,4);
            else
            {
                for(long long cur_offset = 0;cur_offset < (env.thread_size >> 1);cur_offset += env.qubit_size[targ[0]] << 1)
                {
                    inner_all_thread(task,g,env.qubit_size[targ[0]],4);
                    update_offset(task,env.qubit_size[targ[0]]);
                }
            }
        }
        else
        {
            int fd0 = tid & (~file_mask_left);
            int fd1 = tid | file_mask_left;
            if(env.thread_state == env.chunk_state && tid == fd1) return;
            task.fd_using = {env.fd_arr[fd0],env.fd_arr[fd1]};
            func_loop_size = (env.thread_state == env.chunk_state)? env.thread_size : env.thread_size >> 1;
            unordered_map<int,vector<long long>>m = {{fd0,{0,0}},{fd1,{func_loop_size,func_loop_size}}};
            task.fd_offset_using = m[tid];
            task.gate_buffer_using = {0,1};
            inner_all_thread(task,g,func_loop_size,2);
        }
    }
}
void IO_Runner::all_thread_drive_vs2_2(thread_IO_task &task,Gate * &g)
{
    vector<int> targ = g->targs;
    int tid = task.fd_table[0];
    long long func_loop_size;
    if(isFile(targ[2]))//L L D D
    {
        int file_mask_left = 1 << (g->targs[3] - seg.middle - seg.chunk);
        int file_mask_right = 1 << (g->targs[2] - seg.middle - seg.chunk);
        int fd0 = tid & (~file_mask_left) & (~file_mask_right);
        int fd1 = fd0 | file_mask_right;
        int fd2 = fd0 | file_mask_left;
        int fd3 = fd1 | file_mask_left;
        task.fd_using = {env.fd_arr[fd0],env.fd_arr[fd1],env.fd_arr[fd2],env.fd_arr[fd3]};
        int shift = min(int(log2(env.thread_state / env.chunk_state)),2);
        func_loop_size = env.thread_size >> shift;
        vector<long long>off0 = {0}; off0.resize(4,off0[0]);
        vector<long long>off1 = {func_loop_size}; off1.resize(4,off1[0]);
        vector<long long>off2 = {func_loop_size << 1};off2.resize(4,off2[0]);
        vector<long long>off3 = {func_loop_size * 3};off3.resize(4,off3[0]);
        unordered_map<int,vector<long long>>m = {{fd0,off0},{fd1,off1},{fd2,off2},{fd3,off3}};
        task.fd_offset_using = m[tid];
        if((shift == 0) && (tid != fd0)) return;
        else if((shift == 1) && (tid != fd0 && tid != fd1)) return;
        task.gate_buffer_using = {0,1,2,3};
        inner_all_thread(task,g,func_loop_size,4);
    }
    else//L L M D
    {
        int file_mask = 1 << (g->targs[3] - seg.middle - seg.chunk);
        int fd0 = tid & (~file_mask);
        int fd1 = tid | file_mask;
        task.fd_using = {env.fd_arr[fd0],env.fd_arr[fd0],env.fd_arr[fd1],env.fd_arr[fd1]};
        long long base1 = 0;
        long long base2 = env.qubit_size[targ[2]];
        bool cond1 = env.thread_state == (env.qubit_offset[targ[2]] << 1);
        bool cond2 = env.chunk_state == env.qubit_offset[targ[2]];
        func_loop_size = (cond1 && !cond2)?(env.qubit_size[targ[2]] >> 1) : env.qubit_size[targ[2]];
        if(cond1 && cond2 && tid == fd1) return;
        task.fd_offset_using = {base1,base2,base1,base2};
        if(cond1 && !cond2 && tid == fd1) task.fd_offset_using = {base1 + (env.qubit_size[targ[2]] >> 1),base2 + (env.qubit_size[targ[2]] >> 1),base1 + (env.qubit_size[targ[2]] >> 1),base2 + (env.qubit_size[targ[2]] >> 1)};
        if(!cond1 && tid == fd1) task.fd_offset_using = {base1 + (env.thread_size >> 1),base2 + (env.thread_size >> 1),base1 + (env.thread_size >> 1),base2 + (env.thread_size >> 1)};
        task.gate_buffer_using = {0,1,2,3};
        if(cond1)
            inner_all_thread(task,g,func_loop_size,4);
        else
        {
            for(long long cur_offset = 0;cur_offset < (env.thread_size >> 1);cur_offset += env.qubit_size[targ[2]] << 1)
            {
                inner_all_thread(task,g,env.qubit_size[targ[2]],4);
                update_offset(task,env.qubit_size[targ[2]]);
            }
        }
    }
}
void IO_Runner::inner_all_thread(thread_IO_task &task,Gate * &g,long long func_loop_size,int round)
{
    for(long long cur_offset = 0;cur_offset < func_loop_size;cur_offset += env.chunk_size)
    {
        for(int j = 0;j < round;j++)
        {
            if(skip_read_write(g,j)) continue;
            if(pread(task.fd_using[j],&task.buffer1[task.gate_buffer_using[j] * env.chunk_state],env.chunk_size,task.fd_offset_using[j]));
        }
        g->run(task.buffer1);
        for(int j = 0;j < round;j++)
        {
            if(skip_read_write(g,j)) continue;
            if(pwrite(task.fd_using[j],&task.buffer1[task.gate_buffer_using[j] * env.chunk_state],env.chunk_size,task.fd_offset_using[j]));
        }
        update_offset(task,env.chunk_size);
    }
}
void IO_Runner::MPI_special_gate_inner(thread_IO_task &task,Gate * &g,long long func_loop_size,int round)
{
    long long stride = env.chunk_size;
    int lowest_idx = 3;
    if(g->name == "CPhase_Gate" && isMpi(g->targs[1]) && isFile(g->targs[0]))
        stride = env.chunk_size << 1;
    for(long long cur_offset = 0;cur_offset < func_loop_size;cur_offset += stride)
    {
        for(int j = 0;j < round;j++)
        {
            if(skip_read_write(g,j)) continue;
            if(pread(task.fd_using[j],&task.buffer1[task.gate_buffer_using[j] * env.chunk_state],env.chunk_size,task.fd_offset_using[j]));
            lowest_idx = min(lowest_idx,task.gate_buffer_using[j]);
        }
        if(g->name == "Z_Gate" || g->name == "Phase_Gate" || g->name == "RZ_Gate" || g->name == "RZZ_Gate" || g->name == "CPhase_Gate")
            g->run(task.buffer1);
        else
            g->run_one_qubit_mpi_io(task.buffer1,task.buffer2,lowest_idx,lowest_idx + 1,1);
        for(int j = 0;j < round;j++)
        {
            if(skip_read_write(g,j)) continue;
            if(pwrite(task.fd_using[j],&task.buffer1[task.gate_buffer_using[j] * env.chunk_state],env.chunk_size,task.fd_offset_using[j]));
        }
        update_offset(task,stride);
    }
}
void IO_Runner::update_offset(thread_IO_task &task,long long &off)
{
    for(auto &x:task.fd_offset_using)
        x+=off;
}
void IO_Runner::MPI_vs2_2(thread_IO_task &task,Gate * &g)
{
    int mpi_targ_mask_2 = 1 << (g->targs[2] - seg.N);
    int mpi_targ_mask_3 = 1 << (g->targs[3] - seg.N);
    int rank0 = env.rank & (~(mpi_targ_mask_3 | mpi_targ_mask_2));
    int rank1 = rank0 | mpi_targ_mask_2;
    int rank2 = rank0 | mpi_targ_mask_3;
    int rank3 = rank0 | mpi_targ_mask_3 | mpi_targ_mask_2;
    long long loop_bound = env.thread_size;
    long long loop_stride = env.chunk_size * min((long long) MPI_buffer_size,env.thread_state / env.chunk_state);
    int one_round_chunk = min((long long) MPI_buffer_size,env.thread_state / env.chunk_state);
    task.partner_using = {rank0,rank1,rank2,rank3};
    task.fd_using = {env.fd_arr[task.tid],env.fd_arr[task.tid],env.fd_arr[task.tid],env.fd_arr[task.tid]};
    task.fd_offset_using = {0,0,0,0};
    int member_th = find(task.partner_using.begin(),task.partner_using.end(),env.rank) - task.partner_using.begin();
    vector<vector<complex<double>>*>buffer_using = {&task.buffer1,&task.buffer2,&task.buffer3,&task.buffer4};
    vector<complex<double>>*rank_buffer = buffer_using[member_th];
    for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
    {
        if(pread(task.fd_using[0],rank_buffer->data(),one_round_chunk * env.chunk_size,task.fd_offset_using[0]));
        for(int i = 0;i < 4;i++)
        {
            if(member_th != i)
            {
                MPI_Isend(rank_buffer->data(),one_round_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[i],task.tid,MPI_COMM_WORLD,&task.request[i]);
                MPI_Irecv(buffer_using[i]->data(),one_round_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[i],task.tid,MPI_COMM_WORLD,&task.request[i + 4]);
            }
        }
        for(int i = 0;i < 4;i++)
        {
            if(member_th != i)
            {
                MPI_Wait(&task.request[i],MPI_STATUS_IGNORE);
                MPI_Wait(&task.request[i + 4],MPI_STATUS_IGNORE);
            }
        }
        g->run_mpi_vswap2_2_io(task.buffer1,task.buffer2,task.buffer3,task.buffer4,one_round_chunk);
        if(pwrite(task.fd_using[0],rank_buffer->data(),one_round_chunk * env.chunk_size,task.fd_offset_using[0]));
        update_offset(task,loop_stride);
    }
}
void IO_Runner::MPI_Swap(thread_IO_task &task,Gate * &g)
{
    int mpi_targ_mask_0 = 1 << (g->targs[0] - seg.N);
    int mpi_targ_mask_1 = 1 << (g->targs[1] - seg.N);
    long long loop_bound = env.thread_size;
    long long loop_stride;
    task.partner_using = {int((long long)env.rank ^ mpi_targ_mask_1)};
    if(isMpi(g->targs[0]))
    {
        int rank0 = env.rank & (~mpi_targ_mask_0) & ~(mpi_targ_mask_1);
        int rank1 = rank0 | mpi_targ_mask_0;
        int rank2 = rank0 | mpi_targ_mask_1;
        int rank3 = rank0 | mpi_targ_mask_0 | mpi_targ_mask_1;
        unordered_map<int,int>rank_mapping{{rank0,0},{rank1,1},{rank2,2},{rank3,3}};
        int rank_order = rank_mapping[env.rank];
        if(rank_order == 0 || rank_order == 3) return;

        task.fd_using = {env.fd_arr[task.tid]};
        task.fd_offset_using = {0};
        loop_stride = env.chunk_size;
        task.partner_using = {env.rank == rank1? rank2 : rank1};

        int one_round_chunk = min((long long) MPI_buffer_size,env.thread_state / env.chunk_state);
        loop_stride = one_round_chunk * env.chunk_size;
        for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
        {
            for(int i = 0;i < one_round_chunk;i++)
            {
                if(pread(task.fd_using[0],&task.buffer1[i * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + i * env.chunk_size)); 
            }
            MPI_Sendrecv(&task.buffer1[0],one_round_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,&task.buffer2[0],one_round_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            for(int i = 0;i < one_round_chunk;i++)
            {
                if(pwrite(task.fd_using[0],&task.buffer2[i * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + i * env.chunk_size)); 
            }
            update_offset(task,loop_stride);
        }
    }
    else if(isFile(g->targs[0]))
    {
        int file_mask = 1 << (g->targs[0] - seg.middle - seg.chunk);
        bool thread_equal_chunk = (env.thread_size == env.chunk_size);
        int t0 = task.tid & (~file_mask);
        int t1 = task.tid | file_mask;
        int represent_thread = (env.rank < task.partner_using[0])? t1 : t0;
        task.fd_using = {env.fd_arr[represent_thread]};
        task.fd_offset_using = {0};
        if(thread_equal_chunk)
        {
            int recv_tag = represent_thread == t1? t0 : t1;
            if(task.tid != represent_thread) return;
            if(pread(task.fd_using[0],&task.buffer1[0],env.chunk_size,task.fd_offset_using[0]));
            MPI_Sendrecv(&task.buffer1[0],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,&task.buffer2[0],env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],recv_tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            if(pwrite(task.fd_using[0],&task.buffer2[0],env.chunk_size,task.fd_offset_using[0]));
        }
        else
        {
            loop_bound = loop_bound >> 1;
            long long one_round_chunk = min((long long) MPI_buffer_size,(env.thread_state / env.chunk_state) >> 1);
            loop_stride = one_round_chunk * env.chunk_size;
            for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
            {
                long long startIdx = cur_offset + (task.tid != represent_thread) * (env.thread_size >> 1);
                for(int j = 0;j < one_round_chunk;j++)
                {
                    if(pread(task.fd_using[0],&task.buffer1[j * env.chunk_state],env.chunk_size,startIdx + j * env.chunk_size));
                }
                MPI_Isend(&task.buffer1[0],one_round_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,&task.request[0]);
                MPI_Irecv(&task.buffer2[0],one_round_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,&task.request[1]);
                MPI_Wait(&task.request[1],MPI_STATUS_IGNORE);
                for(int j = 0;j < one_round_chunk;j++)
                {
                    if(pwrite(task.fd_using[0],&task.buffer2[j * env.chunk_state],env.chunk_size,startIdx + j * env.chunk_size));
                }
                MPI_Wait(&task.request[0],MPI_STATUS_IGNORE);
            }
        }
    }
    else if(isMiddle(g->targs[0]))
    {
        task.fd_using = {env.fd_arr[task.tid]};
        int total_chunk_per_thread = (env.thread_state >> 1) / env.chunk_state;
        int per_chunk = min(total_chunk_per_thread,env.MPI_buffer_size);
        stack<long long>st;
        for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset+=env.qubit_size[g->targs[0]] << 1)
        {
            for(long long j = 0;j < env.qubit_size[g->targs[0]];j+=env.chunk_size)
            {
                long long startIdx = cur_offset + j + (env.rank < task.partner_using[0]) * env.qubit_size[g->targs[0]];
                if(pread(task.fd_using[0],&task.buffer1[st.size() * env.chunk_state],env.chunk_size,startIdx));
                st.push(startIdx);
                if(st.size() == (unsigned long int)per_chunk)
                {
                    MPI_Sendrecv(&task.buffer1[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,&task.buffer2[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    int n = st.size();
                    for(int q = 0;q < per_chunk;q++)
                    {
                        if(pwrite(task.fd_using[0],&task.buffer2[(n - 1 - q) * env.chunk_state],env.chunk_size,st.top()));
                        st.pop();
                    }
                }
            }
        }
    }
    else //Use new algo
    {
        task.fd_using = {env.fd_arr[task.tid]};
        task.fd_offset_using = {0};
        loop_stride = env.chunk_size * min((long long) MPI_buffer_size,env.thread_state / env.chunk_state);
        for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
        {
            _mpi_one_gate_inner(task,g);
            update_offset(task,loop_stride);
        }
    }
}
void IO_Runner::MPI_gate_scheduler(thread_IO_task &task,Gate * &g)
{
    if(g->name == "Z_Gate" || g->name == "Phase_Gate" || g->name == "RZ_Gate")
    {
        MPI_one_qubit_gate_diagonal(task,g);
        return;
    }
    else if(g->name == "CPhase_Gate" || g->name == "RZZ_Gate")
    {
        MPI_two_qubit_gate_diagonal(task,g);
        return;
    }
    int mpi_targ_mask_0 = 1 << (g->targs[1] - seg.N);
    int mpi_targ_mask_2 = 1 << (g->targs[0] - seg.N);
    long long loop_bound = env.thread_size;
    long long loop_stride;
    int rank0;
    int rank1;
    int rank2;
    int rank3;
    if(g->targs.size() == 1)
    {
        task.fd_using = {env.fd_arr[task.tid]};
        task.fd_offset_using = {0};
        task.partner_using = {int(((long long)env.rank) ^ mpi_targ_mask_2)};
        loop_stride = env.chunk_size * min((long long) MPI_buffer_size,env.thread_state / env.chunk_state);
        for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
        {
            _mpi_one_gate_inner(task,g);
            update_offset(task,loop_stride);
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
            task.partner_using = {rank0,rank1,rank2,rank3};
            task.fd_using = {env.fd_arr[task.tid]};
            task.fd_offset_using = {0};
            int total_chunk_per_thread = env.thread_state / env.chunk_state;
            int per_chunk = min(total_chunk_per_thread,env.MPI_buffer_size >> 2);
            loop_stride = per_chunk * env.chunk_size;
            vector<complex<double>>*buffer1_ptr = &task.buffer1;
            vector<complex<double>>*buffer2_ptr = &task.buffer2;
            vector<complex<double>>*buffer3_ptr = &task.buffer3;
            vector<complex<double>>*buffer4_ptr = &task.buffer4;
            if(env.rank == task.partner_using[1])
            {
                buffer1_ptr = &task.buffer2;
                buffer2_ptr = &task.buffer1;
            }
            else if(env.rank == task.partner_using[2])
            {
                buffer1_ptr = &task.buffer2;
                buffer2_ptr = &task.buffer3;
                buffer3_ptr = &task.buffer1;
            }
            else if(env.rank == task.partner_using[3])
            {
                buffer1_ptr = &task.buffer2;
                buffer2_ptr = &task.buffer3;
                buffer3_ptr = &task.buffer4;
                buffer4_ptr = &task.buffer1;
            }
            for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
            {
                for(int i = 0;i < per_chunk;i++)
                {
                    if(pread(task.fd_using[0],&task.buffer1[i * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + i * env.chunk_size));
                }
                if(env.rank == task.partner_using[0])
                {
                    MPI_Isend(&task.buffer1[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[1],task.tid,MPI_COMM_WORLD,&task.request[0]);
                    MPI_Isend(&task.buffer1[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[2],task.tid,MPI_COMM_WORLD,&task.request[1]);
                    MPI_Isend(&task.buffer1[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[3],task.tid,MPI_COMM_WORLD,&task.request[2]);
                    MPI_Irecv(&task.buffer2[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[1],task.tid,MPI_COMM_WORLD,&task.request[3]);
                    MPI_Irecv(&task.buffer3[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[2],task.tid,MPI_COMM_WORLD,&task.request[4]);
                    MPI_Irecv(&task.buffer4[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[3],task.tid,MPI_COMM_WORLD,&task.request[5]);
                }
                else if(env.rank == task.partner_using[1])
                {
                    MPI_Isend(&task.buffer1[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,&task.request[0]);
                    MPI_Isend(&task.buffer1[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[2],task.tid,MPI_COMM_WORLD,&task.request[1]);
                    MPI_Isend(&task.buffer1[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[3],task.tid,MPI_COMM_WORLD,&task.request[2]);
                    MPI_Irecv(&task.buffer2[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,&task.request[3]);
                    MPI_Irecv(&task.buffer3[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[2],task.tid,MPI_COMM_WORLD,&task.request[4]);
                    MPI_Irecv(&task.buffer4[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[3],task.tid,MPI_COMM_WORLD,&task.request[5]);
                }
                else if(env.rank == task.partner_using[2])
                {
                    MPI_Isend(&task.buffer1[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,&task.request[0]);
                    MPI_Isend(&task.buffer1[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[1],task.tid,MPI_COMM_WORLD,&task.request[1]);
                    MPI_Isend(&task.buffer1[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[3],task.tid,MPI_COMM_WORLD,&task.request[2]);
                    MPI_Irecv(&task.buffer2[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,&task.request[3]);
                    MPI_Irecv(&task.buffer3[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[1],task.tid,MPI_COMM_WORLD,&task.request[4]);
                    MPI_Irecv(&task.buffer4[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[3],task.tid,MPI_COMM_WORLD,&task.request[5]);
                }
                else
                {
                    MPI_Isend(&task.buffer1[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,&task.request[0]);
                    MPI_Isend(&task.buffer1[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[1],task.tid,MPI_COMM_WORLD,&task.request[1]);
                    MPI_Isend(&task.buffer1[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[2],task.tid,MPI_COMM_WORLD,&task.request[2]);
                    MPI_Irecv(&task.buffer2[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[0],task.tid,MPI_COMM_WORLD,&task.request[3]);
                    MPI_Irecv(&task.buffer3[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[1],task.tid,MPI_COMM_WORLD,&task.request[4]);
                    MPI_Irecv(&task.buffer4[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,task.partner_using[2],task.tid,MPI_COMM_WORLD,&task.request[5]);
                }
                MPI_Waitall(6,&task.request[0],MPI_STATUSES_IGNORE);
                g->run_mpi_u2_nonchunk_io(*buffer1_ptr,*buffer2_ptr,*buffer3_ptr,*buffer4_ptr,per_chunk);
                for(int i = 0;i < per_chunk;i++)
                {
                    if(pwrite(task.fd_using[0],&task.buffer1[i * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + i * env.chunk_size));
                }
                update_offset(task,loop_stride);
            }
        }
        else if(isFile(g->targs[0]))
        {
            int file_mask = 1 << (g->targs[0] - seg.middle - seg.chunk);
            int f0 = task.tid & (~file_mask);
            int f1 = task.tid | file_mask;
            task.fd_using = {env.fd_arr[f0],env.fd_arr[f1]};
            task.fd_offset_using = {0,env.thread_size >> 1};
            int partner_rank = env.rank ^ mpi_targ_mask_0;
            bool thread_equal_chunk = env.thread_state == env.chunk_state;
            int total_chunk_per_thread = (thread_equal_chunk)? 1 : ((env.thread_state / env.chunk_state) >> 1);
            int per_chunk = min(total_chunk_per_thread,thread_equal_chunk? env.MPI_buffer_size : (env.MPI_buffer_size >> 1));
            loop_bound = (thread_equal_chunk)? env.thread_size : env.thread_size >> 1;
            loop_stride = per_chunk * env.chunk_size;
            if(thread_equal_chunk && task.tid == f1) return;
            vector<complex<double>>*buffer1_ptr = &task.buffer1;
            vector<complex<double>>*buffer2_ptr = &task.buffer2;
            vector<complex<double>>*buffer3_ptr = &task.buffer3;
            vector<complex<double>>*buffer4_ptr = &task.buffer4;
            if(env.rank > partner_rank)
            {
                buffer1_ptr = &task.buffer3;
                buffer2_ptr = &task.buffer4;
                buffer3_ptr = &task.buffer1;
                buffer4_ptr = &task.buffer2;
            }
            for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
            {                
                if(task.tid == f0)
                {
                    for(int i = 0;i < per_chunk;i++)
                    {
                        if(pread(task.fd_using[0],&task.buffer1[i * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + i * env.chunk_size));
                        if(pread(task.fd_using[1],&task.buffer2[i * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + i * env.chunk_size));
                    }
                }
                else
                {
                    for(int i = 0;i < per_chunk;i++)
                    {
                        if(pread(task.fd_using[0],&task.buffer1[i * env.chunk_state],env.chunk_size,task.fd_offset_using[1] + i * env.chunk_size));
                        if(pread(task.fd_using[1],&task.buffer2[i * env.chunk_state],env.chunk_size,task.fd_offset_using[1] + i * env.chunk_size));
                    }
                }
                MPI_Isend(&task.buffer1[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request[0]);
                MPI_Isend(&task.buffer2[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request[1]);
                MPI_Irecv(&task.buffer3[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request[2]);
                MPI_Irecv(&task.buffer4[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request[3]);
                MPI_Waitall(4,&task.request[0],MPI_STATUSES_IGNORE);
                g->run_mpi_u2_nonchunk_io(*buffer1_ptr,*buffer2_ptr,*buffer3_ptr,*buffer4_ptr,per_chunk);
                if(task.tid == f0)
                {
                    for(int i = 0;i < per_chunk;i++)
                    {
                        if(pwrite(task.fd_using[0],&task.buffer1[i * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + i * env.chunk_size));
                        if(pwrite(task.fd_using[1],&task.buffer2[i * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + i * env.chunk_size));
                    }
                }
                else
                {
                    for(int i = 0;i < per_chunk;i++)
                    {
                        if(pwrite(task.fd_using[0],&task.buffer1[i * env.chunk_state],env.chunk_size,task.fd_offset_using[1] + i * env.chunk_size));
                        if(pwrite(task.fd_using[1],&task.buffer2[i * env.chunk_state],env.chunk_size,task.fd_offset_using[1] + i * env.chunk_size));
                    }
                }
                update_offset(task,loop_stride);
            }
        }
        else if(isMiddle(g->targs[0]))
        {
            task.fd_using = {env.fd_arr[task.tid]};
            int partner_rank = env.rank ^ mpi_targ_mask_0;
            int total_chunk_per_thread = (env.thread_state >> 1) / env.chunk_state;
            int per_chunk = min(total_chunk_per_thread,MPI_buffer_size >> 1);
            loop_stride = env.qubit_size[g->targs[0]] << 1;
            stack<long long>st;
            vector<complex<double>>*buffer1_ptr = &task.buffer1;
            vector<complex<double>>*buffer2_ptr = &task.buffer2;
            vector<complex<double>>*buffer3_ptr = &task.buffer3;
            vector<complex<double>>*buffer4_ptr = &task.buffer4;
            if(env.rank > partner_rank)
            {
                swap(buffer1_ptr,buffer3_ptr);
                swap(buffer2_ptr,buffer4_ptr);
            }
            for(long long cur_offset = 0;cur_offset < env.thread_size;cur_offset += loop_stride)
            {
                task.fd_offset_using = {cur_offset,cur_offset + env.qubit_size[g->targs[0]]};
                for(long long chunk_off = 0;chunk_off < env.qubit_size[g->targs[0]];chunk_off += env.chunk_size)
                {
                    if(pread(task.fd_using[0],&task.buffer1[st.size() * env.chunk_state],env.chunk_size,task.fd_offset_using[0]));
                    if(pread(task.fd_using[0],&task.buffer2[st.size() * env.chunk_state],env.chunk_size,task.fd_offset_using[1]));
                    st.push(task.fd_offset_using[0]);
                    update_offset(task,env.chunk_size);
                    if(st.size() == (unsigned long int)per_chunk)
                    {
                        MPI_Isend(&task.buffer1[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request[0]);
                        MPI_Isend(&task.buffer2[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request[1]);
                        MPI_Irecv(&task.buffer3[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request[2]);
                        MPI_Irecv(&task.buffer4[0],per_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,&task.request[3]);
                        MPI_Waitall(4,&task.request[0],MPI_STATUSES_IGNORE);
                        g->run_mpi_u2_nonchunk_io(*buffer1_ptr,*buffer2_ptr,*buffer3_ptr,*buffer4_ptr,per_chunk);
                        int n = st.size();
                        for(int i = 0;i < n;i++)
                        {
                            long long off = st.top();st.pop();
                            if(pwrite(task.fd_using[0],&task.buffer1[(n - i - 1) * env.chunk_state],env.chunk_size,off));
                            if(pwrite(task.fd_using[0],&task.buffer2[(n - i - 1) * env.chunk_state],env.chunk_size,off + env.qubit_size[g->targs[0]]));
                        }
                    }
                }
            }
        }
        else // use new algo
        {
            rank0 = env.rank ^ mpi_targ_mask_0;
            task.partner_using = {rank0};
            task.fd_using = {env.fd_arr[task.tid]};
            task.fd_offset_using = {0};
            loop_stride = env.chunk_size * min((long long) MPI_buffer_size,env.thread_state / env.chunk_state);
            for(long long cur_offset = 0;cur_offset < loop_bound;cur_offset += loop_stride)
            {
                _mpi_one_gate_inner(task,g);
                update_offset(task,loop_stride);
            }
        }
    }
}

void IO_Runner::_mpi_one_gate_inner(thread_IO_task &task,Gate * &g)
{
    int partner_rank = task.partner_using[0];
    int member_th = env.rank > partner_rank;
    int one_round_chunk = min((long long) MPI_buffer_size,env.thread_state / env.chunk_state);
    for(int i = 0;i < one_round_chunk;i++)
    {
        if(pread(task.fd_using[0],&task.buffer1[i * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + i * env.chunk_size)); 
    }
    MPI_Sendrecv(&task.buffer1[0],one_round_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,&task.buffer2[0],one_round_chunk * env.chunk_state,MPI_DOUBLE_COMPLEX,partner_rank,task.tid,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    vector<complex<double>> *buffer1 = &task.buffer1;
    vector<complex<double>> *buffer2 = &task.buffer2;
    if(member_th) swap(buffer1,buffer2);
    g->run_one_qubit_mpi_io(*buffer1,*buffer2,0,0,one_round_chunk);
    for(int i = 0;i < one_round_chunk;i++)
    {
        if(pwrite(task.fd_using[0],&task.buffer1[i * env.chunk_state],env.chunk_size,task.fd_offset_using[0] + i * env.chunk_size)); 
    }
}
void IO_Runner::MPI_one_qubit_gate_diagonal(thread_IO_task &task,Gate * &g)
{
    int mpi_mask = 1 << (g->targs[0] - seg.N);
    task.fd_using = {env.fd_arr[task.tid]};
    task.fd_offset_using = {0};
    if(g->name == "Z_Gate" || g->name == "Phase_Gate")
    {        
        if(env.rank != (env.rank | mpi_mask)) return;
        task.gate_buffer_using = {1};
        MPI_special_gate_inner(task,g,env.thread_size,1);
    }
    else if(g->name == "RZ_Gate")
    {
        task.gate_buffer_using = {env.rank == (env.rank | mpi_mask)};
        MPI_special_gate_inner(task,g,env.thread_size,1);
    }
}
void IO_Runner::MPI_two_qubit_gate_diagonal(thread_IO_task &task,Gate * &g)
{
    int mpi_mask1 = 1 << (g->targs[1] - seg.N);
    int mpi_mask0 = 1 << (g->targs[0] - seg.N);
    if(isMpi(g->targs[0]))
    {
        int rank0 = env.rank & (~(mpi_mask0 | mpi_mask1));
        int rank1 = rank0 | mpi_mask0;
        int rank2 = rank0 | mpi_mask1;
        int rank3 = rank0 | mpi_mask0 | mpi_mask1;
        if(g->name == "CPhase_Gate")
        {
            if(env.rank != rank3) return;
            task.fd_using = {env.fd_arr[task.tid]};
            task.fd_offset_using = {0};
            task.gate_buffer_using = {3};
            MPI_special_gate_inner(task,g,env.thread_size,1);
        }
        else if(g->name == "RZZ_Gate")
        {
            task.fd_using = {env.fd_arr[task.tid]};
            task.fd_offset_using = {0};
            unordered_map<int,int>m = {
                {rank0,0},
                {rank1,1},
                {rank2,2},
                {rank3,3},
            };
            task.gate_buffer_using = {m[env.rank]};
            MPI_special_gate_inner(task,g,env.thread_size,1);
        }
    }
    else if(isFile(g->targs[0]))
    {
        int file_mask = 1 << (g->targs[0] - seg.middle - seg.chunk);
        int partner_rank = env.rank ^ mpi_mask1;
        if(g->name == "CPhase_Gate")
        {
            if(env.rank < partner_rank) return;
            int fd1 = task.tid | file_mask;
            bool thread_equal_chunk = (env.thread_size == env.chunk_size);
            if(thread_equal_chunk && task.tid == fd1) return;
            task.fd_using = {env.fd_arr[fd1]};
            task.fd_offset_using = (task.tid != fd1)? vector<long long> {0} : vector<long long> {env.chunk_size};
            task.gate_buffer_using = {3};
            MPI_special_gate_inner(task,g,env.thread_size,1);
        }
        else if(g->name == "RZZ_Gate")
        {
            int rank_member = (env.rank > partner_rank);
            int file_member = (task.tid == (task.tid | file_mask));
            task.fd_using = {env.fd_arr[task.tid]};
            task.fd_offset_using = {0};
            task.gate_buffer_using = {(rank_member << 1) | file_member};
            MPI_special_gate_inner(task,g,env.thread_size,1);
        }
    }
    else if(isMiddle(g->targs[0]))
    {
        int partner_rank = env.rank ^ mpi_mask1;
        if(g->name == "CPhase_Gate")
        {
            if(env.rank < partner_rank) return;
            task.fd_using = {env.fd_arr[task.tid]};
            task.gate_buffer_using = {3};
            for(long long cur_offset = 0;cur_offset < env.thread_size;cur_offset += env.qubit_size[g->targs[0]] << 1)
            {
                task.fd_offset_using = {cur_offset + env.qubit_size[g->targs[0]]};
                MPI_special_gate_inner(task,g,env.qubit_size[g->targs[0]],1);
            }
        }
        else if(g->name == "RZZ_Gate")
        {
            int rank_member = (env.rank > partner_rank);
            task.fd_using = {env.fd_arr[task.tid],env.fd_arr[task.tid]};
            task.gate_buffer_using = {rank_member << 1,(rank_member << 1) | 1};
            for(long long cur_offset = 0;cur_offset < env.thread_size;cur_offset += env.qubit_size[g->targs[0]] << 1)
            {
                task.fd_offset_using = {cur_offset,cur_offset + env.qubit_size[g->targs[0]]};
                MPI_special_gate_inner(task,g,env.qubit_size[g->targs[0]],2);
            }
        }
    }
    else
    {
        int partner_rank = env.rank ^ mpi_mask1;
        if(g->name == "CPhase_Gate")
        {
            if(env.rank < partner_rank) return;
            task.fd_using = {env.fd_arr[task.tid]};
            task.fd_offset_using = {0};
            task.gate_buffer_using = {1};
            MPI_special_gate_inner(task,g,env.thread_size,1);
        }
        else if(g->name == "RZZ_Gate")
        {
            int rank_member = (env.rank > partner_rank);
            task.fd_using = {env.fd_arr[task.tid]};
            task.fd_offset_using = {0};
            task.gate_buffer_using = {rank_member};
            MPI_special_gate_inner(task,g,env.thread_size,1);
        }
    }
}