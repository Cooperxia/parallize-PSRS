#include<iostream>
#include<fstream>
#include<time.h>
#include<string>
#include<algorithm>
#include<mpi.h>
#include<omp.h>
using namespace std;

#define BLOCK 1024*1024
int PNUM = 0;
int LOOP = 0;
char* NAME = NULL;

int comm_sz;
int my_rank;
int totalSize;
int blockSize;
int lastSize;

//用来存储中间结果的函数，方便调试以及查看运行结果的正确性
void store(int* data, string name, int blockSize)
{
	ofstream file(name);
	for (int i = 0; i < blockSize; i++)
		file << data[i] << " ";
	file.close();
}

template<typename E>
void quickSort(E *A, int left, int right)
{
	if (left >= right) return;
	int x = A[(left + right) >> 1], low = left, high = right;
	while (low<high)
	{
		while (A[low]<x)
			low++;
		while (A[high]>x)
			high--;
		if (low <= high)
		{
			int Temp = A[low];
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
void judge1(int *temp, int size, int *block, int *&num_blockt)
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
void judge2(int blockSize, int *temp, int *block, int *&num_block)
{
	int size = blockSize / PNUM;
	int last = blockSize - (PNUM - 1)*size;
#pragma omp parallel for num_threads(PNUM)
	for (int i = 0; i < PNUM; i++)
	{
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

void PSRS(int* block, int*& data)
{
	double begin_time = 0, end_time = 0;


	begin_time = MPI_Wtime();

	int *temp = new int[comm_sz*comm_sz];//从p个处理器中分别提出的p个数据
	int diff = blockSize / (comm_sz + 1);//取主元时的间隔

	if (my_rank != comm_sz - 1)
		quickSort<int>(block, 0, blockSize - 1);//块内快排
	else//分配最后一个Block
		quickSort<int>(block, 0, lastSize - 1);//块内快排

	int *tempt = new int[comm_sz];

	for (int i = 0; i < comm_sz; i++)
		tempt[i] = block[(i + 1)*diff];
	MPI_Gather(tempt, comm_sz, MPI_INT, temp, comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
	if (my_rank == 0)
	{
		quickSort<int>(temp, 0, comm_sz*comm_sz - 1);
		for (int i = 1; i < comm_sz; i++)//选取p-1个主元
			temp[i - 1] = temp[i*(comm_sz - 1)];//temp的0至comm_sz-2存comm_sz-1个有效主元
												//store(temp, "block主元0.txt", comm_sz - 1);
	}
	MPI_Bcast(temp, comm_sz - 1, MPI_INT, 0, MPI_COMM_WORLD);

	//全局交换思路：存储每一段的数据个数

	int *num_block = new int[comm_sz];
	int *begin = new int[comm_sz];
	for (int i = 0; i < comm_sz; i++)
	{
		num_block[i] = 0;
		begin[i] = 0;
	}

	if (my_rank != comm_sz - 1)
		judge2(blockSize, temp, block, num_block);
	else
		judge2(lastSize, temp, block, num_block);

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
	int *new_block = new int[total];
	int *recvDisp = new int[comm_sz];
	recvDisp[0] = 0;

	//recvDisp：该进程每次接收数据时的偏移量
	for (int i = 1; i < comm_sz; i++)
		recvDisp[i] = num[i - 1] + recvDisp[i - 1];

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Alltoallv(block, num_block, begin, MPI_INT, new_block, num, recvDisp, MPI_INT, MPI_COMM_WORLD);

	quickSort<int>(new_block, 0, total - 1);

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
	MPI_Gatherv(new_block, total, MPI_INT, data, num, recvDisp, MPI_INT, 0, MPI_COMM_WORLD);

	end_time = MPI_Wtime();
	if (my_rank == 0)
		cout << "time: " << end_time - begin_time << endl;

	if (my_rank == 0)
		store(data, "total.txt", BLOCK*LOOP);

	MPI_Finalize();
}

//输入参数：进程数PNUM，文件名NAME，循环次数LOOP
int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		cout << "参数数目不对" << endl;
		return 0;
	}
	PNUM = atoi(argv[1]);
	NAME = argv[2];
	LOOP = atoi(argv[3]);

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	totalSize = BLOCK*LOOP;
	blockSize = totalSize / comm_sz;
	lastSize = totalSize - ((comm_sz - 1)*blockSize);

	int *block = NULL;
	int *data = NULL;
	data = (int*)malloc(sizeof(int)*totalSize);
	ifstream file(NAME, ios_base::binary);

	file.seekg(my_rank*blockSize * sizeof(int));
	if (my_rank != comm_sz - 1)
	{
		block = (int *)malloc(sizeof(int)*blockSize);
		file.read((char*)block, blockSize * sizeof(int));
	}
	else
	{
		block = (int *)malloc(sizeof(int)*lastSize);
		file.read((char*)block, lastSize * sizeof(int));
	}
	file.close();

	if (comm_sz == 1)
	{
		double begin_time, end_time;
		begin_time = MPI_Wtime();
		quickSort<int>(data, 0, blockSize - 1);
		end_time = MPI_Wtime();
		cout << "time: " << end_time - begin_time << endl;
	}
	else
		PSRS(block, data);

	delete block;
	return 0;
}