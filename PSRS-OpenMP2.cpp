#include<iostream>
#include<fstream>
#include<fcntl.h>
#include<omp.h>
using namespace std;

#define BLOCK 1024*1024
int PNUM = 0;
int LOOP = 0;
char* NAME = NULL;

//���ڵ����õĴ洢�������������ݵĺ�����*dataΪ����ָ�룬nameΪ�洢�ļ�����sizeΪ�洢�����ݸ���
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
	int* temp = new int[PNUM*PNUM];//��p���������зֱ������p������
	int size = BLOCK*LOOP / PNUM;//ÿ����Ĵ�С
	int last = BLOCK*LOOP - size*(PNUM - 1);//���һ���С(��ʣ��һ���һЩ�Ķ��ָ����һ��)
	int diff = size / PNUM;//ȡ��Ԫʱ�ļ��
	omp_set_num_threads(PNUM);
#pragma omp parallel for
	for (int i = 0; i < PNUM; i++)
	{
		if (i != PNUM - 1)//����PNUM-1��Block
			quickSort<int>(data, i*size, (i + 1)*size - 1);//���ڿ���
		else//�������һ��Block
			quickSort<int>(data, i*size, BLOCK*LOOP - 1);//���ڿ���
		/*string ss;
		ss.append("block");
		ss.append(to_string(i));
		ss.append(".txt");
		store(block[i], ss, size);*/
		for (int j = 0; j < PNUM; j++)//ÿ����ѡȡp����������p*p��
			temp[i*PNUM + j] = data[i*size + j*diff];
	}
	quickSort<int>(temp, 0, PNUM*PNUM - 1);//p2����������
	for (int i = 1; i < PNUM; i++)//ѡȡp-1����Ԫ
		temp[i - 1] = temp[i*(PNUM - 1)];//temp��0��PNUM-2��PNUM-1����Ч��Ԫ
	//store(temp, "block��Ԫ.txt", PNUM - 1);

	//ȫ�ֽ���˼·���洢ÿһ�ε����ݸ���

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
					j++;//������һ����Ԫ�ж�
					begin[i][j] = begin[i][j - 1] + num_block[i][j - 1];
					if (j != PNUM - 1)
						num_block[i][j] = 0;//��һ����ԪƬ�γ��ȳ�ʼΪ0
					else if (k != size)//��Ԫ��ͷ��blockδ��ͷ
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
					j++;//������һ����Ԫ�ж�
					begin[i][j] = begin[i][j - 1] + num_block[i][j - 1];
					if (j != PNUM - 1)
						num_block[i][j] = 0;//��һ����ԪƬ�γ��ȳ�ʼΪ0
					else if (k != last)//��Ԫ��ͷ��blockδ��ͷ
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
	ss.append("block��");
	ss.append(to_string(i));
	ss.append(".txt");
	store(num_block[i], ss, PNUM);
	}
	for (int i = 0; i < PNUM; i++)
	{
	string ss;
	ss.append("block��ʼ");
	ss.append(to_string(i));
	ss.append(".txt");
	store(begin[i], ss, PNUM);
	}*/

	int index = 0;
	int* result = (int *)malloc(sizeof(int)*BLOCK * LOOP);
	int* sum = new int[PNUM];
	int** temp2 = new int*[PNUM];
	//�Ը��Զν�������
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

//���������������PNUM���ļ���NAME��ѭ������LOOP
int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		cout << "������Ŀ����" << endl;
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