#include<iostream>
#include<fstream>
#include<time.h>
#include<string>
#include<algorithm>
#include<mpi.h>
using namespace std;

#define BLOCK 1024*1024
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

	//每个进程选p个主元待选
	for (int i = 0; i < comm_sz; i++)
		tempt[i] = block[(i + 1)*diff];

	//交给主进程
	MPI_Gather(tempt, comm_sz, MPI_INT, temp, comm_sz, MPI_INT, 0, MPI_COMM_WORLD);

	//选取最终的p个主元
	if (my_rank == 0)
	{
		quickSort<int>(temp, 0, comm_sz*comm_sz - 1);
		for (int i = 1; i < comm_sz; i++)//选取p-1个主元
			temp[i - 1] = temp[i*(comm_sz - 1)];//temp的0至comm_sz-2存comm_sz-1个有效主元
												//store(temp, "block主元0.txt", comm_sz - 1);
	}

	//广播主元
	MPI_Bcast(temp, comm_sz - 1, MPI_INT, 0, MPI_COMM_WORLD);

	//全局交换
	int *num_block = new int[comm_sz];//根据主元划分的每一段的长度
	int *begin = new int[comm_sz];//根据主元划分的每一段的起始位置
	num_block[0] = 0;
	begin[0] = 0;

	//根据主元分段
	if (my_rank != comm_sz - 1)
	{
		for (int j = 0, k = 0; (j < comm_sz - 1) && (k < blockSize);)
		{
			if (temp[j] < block[k])
			{
				j++;//进入下一个主元判断
				begin[j] = begin[j - 1] + num_block[j - 1];
				if (j != comm_sz - 1)
					num_block[j] = 0;//下一个主元片段长度初始为0
				else if (k != blockSize)//主元到头，block未到头
					num_block[j] = blockSize - begin[j];
			}
			else
			{
				num_block[j]++;
				k++;
				if ((k == blockSize) && (j != comm_sz - 1))
				{
					while (j != comm_sz - 1)
					{
						num_block[++j] = 0;
						begin[j] = blockSize;
					}
				}
			}
		}
	}
	else
	{
		if (comm_sz == 1)
			num_block[0] = lastSize;
		else
		{
			for (int j = 0, k = 0; (j < comm_sz - 1) && (k < lastSize);)
			{
				if (temp[j] < block[k])
				{
					j++;//进入下一个主元判断
					begin[j] = begin[j - 1] + num_block[j - 1];
					if (j != comm_sz - 1)
						num_block[j] = 0;//下一个主元片段长度初始为0
					else if (k != lastSize)//主元到头，block未到头
						num_block[j] = lastSize - begin[j];
				}
				else
				{
					num_block[j]++;
					k++;
					if ((k == lastSize) && (j != comm_sz - 1))
					{
						while (j != comm_sz - 1)
						{
							num_block[++j] = 0;
							begin[j] = lastSize;
						}
					}
				}
			}
		}
	}

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
	MPI_Gatherv(new_block, total, MPI_INT, data, num, recvDisp, MPI_INT, 0, MPI_COMM_WORLD);

	end_time = MPI_Wtime();
	if (my_rank == 0)
		cout << "time: " << end_time - begin_time << endl;

	if (my_rank == 0)
		store(data, "total.txt", BLOCK*LOOP);

}

//输入参数：进程数PNUM，文件名NAME，循环次数LOOP
int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		cout << "参数数目不对" << endl;
		return 0;
	}
	NAME = argv[1];
	LOOP = atoi(argv[2]);

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

	PSRS(block, data);
	delete block;
	MPI_Finalize();
	return 0;
}