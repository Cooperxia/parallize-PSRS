#include<iostream>
#include<fstream>
#include<fcntl.h>
#include<omp.h>
using namespace std;

#define BLOCK 1024*1024
int PNUM = 0;
int LOOP = 0;
char* NAME = NULL;

//用于调试用的存储程序运行中数据的函数，*data为数据指针，name为存储文件名，size为存储的数据个数
//void store(int* data, string name, int size)
//{
//	ofstream file(name);
//	for (int i = 0; i < size; i++)
//		file << data[i] << " ";
//	file.close();
//}

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

int* PSRS(int* data)
{
	int* temp = new int[PNUM*PNUM];//从p个处理器中分别提出的p个数据
	int size = BLOCK*LOOP / PNUM;//每个块的大小
	int last = BLOCK*LOOP - size*(PNUM - 1);//最后一块大小(将剩余一块多一些的都分给最后一块)
	int diff = size / PNUM;//取主元时的间隔
	omp_set_num_threads(PNUM);
#pragma omp parallel for
	for (int i = 0; i < PNUM; i++)
	{
		if (i != PNUM - 1)//分配PNUM-1个Block
			quickSort<int>(data, i*size, (i + 1)*size - 1);//块内快排
		else//分配最后一个Block
			quickSort<int>(data, i*size, BLOCK*LOOP - 1);//块内快排
		/*string ss;
		ss.append("block");
		ss.append(to_string(i));
		ss.append(".txt");
		store(block[i], ss, size);*/
		for (int j = 0; j < PNUM; j++)//每个块选取p个样本，共p*p个
			temp[i*PNUM + j] = data[i*size + j*diff];
	}
	quickSort<int>(temp, 0, PNUM*PNUM - 1);//p2个样本排序
	for (int i = 1; i < PNUM; i++)//选取p-1个主元
		temp[i - 1] = temp[i*(PNUM - 1)];//temp的0至PNUM-2存PNUM-1个有效主元
	//store(temp, "block主元.txt", PNUM - 1);

	//全局交换思路：存储每一段的数据个数

	int** num_block = new int*[PNUM];
	int** begin = new int*[PNUM];
#pragma omp parallel for
	for (int i = 0; i < PNUM; i++)
	{
		num_block[i] = new int[PNUM];
		begin[i] = new int[PNUM];
		num_block[i][0] = 0;
		begin[i][0] = 0;
		if (i != PNUM - 1)
		{
			for (int j = 0, k = 0; (j < PNUM - 1) && (k < size);)
			{
				if (temp[j] < data[i*size + k])
				{
					j++;//进入下一个主元判断
					begin[i][j] = begin[i][j - 1] + num_block[i][j - 1];
					if (j != PNUM - 1)
						num_block[i][j] = 0;//下一个主元片段长度初始为0
					else if (k != size)//主元到头，block未到头
						num_block[i][j] = size - begin[i][j];
				}
				else
				{
					num_block[i][j]++;
					k++;
					if ((k == size) && (j != PNUM - 1))
					{
						while (j != PNUM - 1)
						{
							num_block[i][++j] = 0;
							begin[i][j] = size;
						}
					}
				}
			}
		}
		else
		{
			for (int j = 0, k = 0; (j < PNUM - 1) && (k < last);)
			{
				if (temp[j] < data[i*size + k])
				{
					j++;//进入下一个主元判断
					begin[i][j] = begin[i][j - 1] + num_block[i][j - 1];
					if (j != PNUM - 1)
						num_block[i][j] = 0;//下一个主元片段长度初始为0
					else if (k != last)//主元到头，block未到头
						num_block[i][j] = last - begin[i][j];
				}
				else
				{
					num_block[i][j]++;
					k++;
					if ((k == last) && (j != PNUM - 1))
					{
						while (j != PNUM - 1)
						{
							num_block[i][++j] = 0;
							begin[i][j] = last;
						}
					}
				}
			}
		}
	}
	/*for (int i = 0; i < PNUM; i++)
	{
	string ss;
	ss.append("block段");
	ss.append(to_string(i));
	ss.append(".txt");
	store(num_block[i], ss, PNUM);
	}
	for (int i = 0; i < PNUM; i++)
	{
	string ss;
	ss.append("block开始");
	ss.append(to_string(i));
	ss.append(".txt");
	store(begin[i], ss, PNUM);
	}*/

	int index = 0;
	int* result = (int *)malloc(sizeof(int)*BLOCK * LOOP);
	int* sum = new int[PNUM];
	int** temp2 = new int*[PNUM];
	//对各自段进行排序
#pragma omp parallel for private(index)
	for (int i = 0; i < PNUM; i++)
	{
		index = 0;
		sum[i] = 0;
		for (int j = 0; j < PNUM; j++)
			sum[i] += num_block[j][i];
		temp2[i] = new int[sum[i]];
		for (int j = 0; j < PNUM; j++)
			for (int k = 0; k < num_block[j][i]; k++)
				temp2[i][index++] = data[j*size + begin[j][i] + k];
		quickSort<int>(temp2[i], 0, sum[i] - 1);
	}
	for (int i = 0, a = 0; i < PNUM; i++)
	{
		for (int j = 0; j < sum[i]; j++)
			result[a + j] = temp2[i][j];
		a += sum[i];
	}
	return result;
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
	int* data = NULL;
	ifstream file(NAME, ios_base::binary);
	data = (int *)malloc(sizeof(int)*BLOCK * LOOP);
	file.read((char*)data, BLOCK * LOOP * sizeof(int));
	double start, finish;
	start = omp_get_wtime();
	data = PSRS(data);
	finish = omp_get_wtime();
	cout << finish - start << endl;
	file.close();
	//store(data, "total.txt", BLOCK*LOOP);
	delete data;
	return 0;
}