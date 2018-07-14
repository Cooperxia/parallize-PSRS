#include<iostream>
#include<fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string>
#include<mpi.h>
#include<omp.h>
#include <cstdio>
using namespace std;

#define BLOCK 1024*1024
#define FLAGS O_WRONLY | O_CREAT | O_TRUNC
#define MODE S_IRWXU | S_IXGRP | S_IROTH | S_IXOTH
int PNUM = 0;
int LOOP = 0;
char* NAME = NULL;

int comm_sz;
int my_rank;
int totalSize;
int blockSize;

void readData(double *block,int i,int flag)
{
	char a[10];
	sprintf(a,"%d",i);
	
	string ss;
	ss.append(a);
	ss.append(".txt");
	
	
	ifstream file(ss.data(), ios_base::binary);
	file.seekg(my_rank*blockSize/2 * sizeof(double));
	file.read((char*)(&block[flag*blockSize/2]), blockSize/2 * sizeof(double));
	
	file.close();
}

//用来存储中间结果的函数，方便调试以及查看运行结果的正确性
void store(double* data, string name, int blockSize)
{
	ofstream file(name, ios_base::binary);
	file.write((char*)data, blockSize * sizeof(double));
	file.close();
}

void check(double* data,int a,int b,int pnum)
{
	int k=a;
	if (b<a)
		return;
	while(data[k]<=data[k+1]&&k!=b)
		k++;
	if (k == b)
		cout << pnum << " is true" << endl;
	else
		cout << pnum << " is false" << endl;
}

void quickSort(double *A, int left, int right)
{
	if (left >= right) return;
	double x = A[(left + right) >> 1];
	int low = left, high = right;
	while (low<high)
	{
		while (A[low]<x)
			low++;
		while (A[high]>x)
			high--;
		if (low <= high)
		{
			double Temp = A[low];
			A[low] = A[high];
			A[high] = Temp;
			low++;
			high--;
		}
	}
	quickSort(A, left, high);
	quickSort(A, low, right);
}

//用于参照主元划分block的函数，输入：temp，size，block，输出：num_blockt
void judge1(double *temp, int size, double *block, int *&num_blockt)
{
	int begint = 0;
	for (int j = 0, k = 0; (j < comm_sz - 1) && (k < size);)
	{
		if (temp[j] < block[k])
		{
			j++;//进入下一个主元判断
			begint = begint + num_blockt[j - 1];
			if (j != comm_sz - 1)
				num_blockt[j] = 0;//下一个主元片段长度初始为0
			else if (k != size)//主元到头，block未到头
				num_blockt[j] = size - begint;
		}
		else
		{
			num_blockt[j]++;
			k++;
			if ((k == size) && (j != comm_sz - 1))
			{
				while (j != comm_sz - 1)
					num_blockt[++j] = 0;
			}
		}
	}
}

//用于并行的处理，judge1的任务，输入：blockSize，temp，block，输出：num_block
void judge2(int blockSize, double *temp, double *block, int *&num_block)
{
	int size = blockSize / PNUM;
	int last = blockSize - (PNUM - 1)*size;
#pragma omp parallel for num_threads(PNUM)
	for (int i = 0; i < PNUM; i++)
	{
		int begint = 0;
		int *num_blockt = new int[comm_sz];
		num_blockt[0] = 0;
		if (i != PNUM - 1)
			judge1(temp, size, &block[i*size], num_blockt);
		else
			judge1(temp, last, &block[i*size], num_blockt);

		//将中间变量num_blockt的值赋给num_block
		for (int j = 0; j < comm_sz; j++)
		{
#pragma omp atomic
			num_block[j] += num_blockt[j];
		}
	}
}

void PSRS(double* block, double*& data)
{
	double begin_time = 0, end_time = 0;


	begin_time = MPI_Wtime();

	double *temp = new double[comm_sz*comm_sz];//从p个处理器中分别提出的p个数据
	int diff = blockSize / (comm_sz + 1);//取主元时的间隔

	quickSort(block, 0, blockSize - 1);//块内快排

	double *tempt = new double[comm_sz];

	for (int i = 0; i < comm_sz; i++)
		tempt[i] = block[(i + 1)*diff];
	MPI_Gather(tempt, comm_sz, MPI_DOUBLE, temp, comm_sz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (my_rank == 0)
	{
		quickSort(temp, 0, comm_sz*comm_sz - 1);
		for (int i = 1; i < comm_sz; i++)//选取p-1个主元
			temp[i - 1] = temp[i*(comm_sz - 1)];//temp的0至comm_sz-2存comm_sz-1个有效主元
												//store(temp, "block主元0.txt", comm_sz - 1);
	}
	MPI_Bcast(temp, comm_sz - 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//全局交换思路：存储每一段的数据个数

	int *num_block = new int[comm_sz];
	int *begin = new int[comm_sz];
	for (int i = 0; i < comm_sz; i++)
	{
		num_block[i] = 0;
		begin[i] = 0;
	}

	judge2(blockSize, temp, block, num_block);

	for (int i = 1; i < comm_sz; i++)
		begin[i] = begin[i - 1] + num_block[i - 1];

	//num是来自comm_sz个进程的，要接收的数据量，例：num[2]是进程2发来的数据量
	int *num = new int[comm_sz];
	MPI_Alltoall(num_block, 1, MPI_INT, num, 1, MPI_INT, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	//total：该进程要接收的数据总量，方便创建动态数组
	int total = 0;
	for (int i = 0; i < comm_sz; i++)
		total += num[i];

	//new_block：该进程接收到的所有数据
	double *new_block = new double[total];
	int *recvDisp = new int[comm_sz];
	recvDisp[0] = 0;

	//recvDisp：该进程每次接收数据时的偏移量
	for (int i = 1; i < comm_sz; i++)
		recvDisp[i] = num[i - 1] + recvDisp[i - 1];

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Alltoallv(block, num_block, begin, MPI_DOUBLE, new_block, num, recvDisp, MPI_DOUBLE, MPI_COMM_WORLD);

	quickSort(new_block, 0, total - 1);

	MPI_Gather(&total, 1, MPI_INT, num, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	if (my_rank == 0)
	{
		recvDisp[0] = 0;
		for (int i = 1; i < comm_sz; i++)
			recvDisp[i] = num[i - 1] + recvDisp[i - 1];
	}

	//主进程收集最终数据
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gatherv(new_block, total, MPI_DOUBLE, data, num, recvDisp, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	end_time = MPI_Wtime();
	// if (my_rank == 0)
		// cout << "time: " << end_time - begin_time << endl;

	if (my_rank == 0)
	{
		int bsize=totalSize/PNUM;
#pragma omp parallel for num_threads(PNUM)
		for(int i=0;i<PNUM;i++)
			check(data,i*bsize,(i+1)*bsize-1,i);
	}
	/*if (my_rank == 0)
		store(data, "total.txt", BLOCK*LOOP);*/

}

//输入参数：进程数PNUM，文件名NAME，循环次数LOOP
int main(int argc, char* argv[])
{

	//路径
	if (argc != 4)
	{
		cout << "Input number wrong" << endl;
		return 0;
	}
	PNUM = atoi(argv[1]);
	LOOP = atoi(argv[2]);
	int num = atoi(argv[3]);
	
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	totalSize = BLOCK*LOOP*2;
	blockSize = totalSize / comm_sz;
	
	for (int i = 2; i <= num;i+=2)
	{
		char a[10];
		sprintf(a,"%d",i/2);
		double *block = NULL;
		double *data = NULL;
		data = (double*)malloc(sizeof(double)*totalSize);
		block = (double *)malloc(sizeof(double)*blockSize);
		readData(block,i-1,0);
		readData(block,i,1);

		if (comm_sz == 1)
		{
			double begin_time, end_time;
			begin_time = MPI_Wtime();
			quickSort(block, 0, blockSize - 1);
			end_time = MPI_Wtime();
			store(block, string("middle_").append(a).append(".txt"), blockSize);
			cout << "time: " << end_time - begin_time << endl;
		}
		else
		{
			PSRS(block, data);
			if(my_rank==0)
				store(data, (string("middle_").append(a).append(".txt")).data(), totalSize);
		}
		
		delete[] data;
		delete[] block;
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	
	MPI_Finalize();

	return 0;
}
