#include<iostream>
#include<fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "unistd.h"
#include <string>
#include <mpi.h>
#include <omp.h>
#include <cstdio>
#include <limits>
using namespace std;

#define BLOCK 1024
int PNUM = 0;
int LOOP = 0;
char* NAME = NULL;

int comm_sz;
int my_rank;
int totalSize;
int blockSize;

void readData(string path, double *block, int i, int flag);
void store(double* data, string name, int blockSize);
void check(double* data, int a, int b, int pnum);
void quickSort(double *A, int left, int right);
int findsmall(double* block, int begin, double smallvalue);
void merge(string path, int filenum, int fileSize, int eachBlockSize);


//���������������PNUM���ļ���NAME��ѭ������LOOP
int main(int argc, char* argv[])
{

	//·��
	if (argc != 5)
	{
		cout << "Input number wrong" << endl;
		return 0;
	}
	char* path = argv[1];
	int filenum = atoi(argv[2]);
	int fileSize = atoi(argv[3]);
	int eachBlockSize = atoi(argv[4]);

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	merge(path, filenum, fileSize*BLOCK, eachBlockSize*BLOCK);

	MPI_Finalize();

	return 0;
}

//Ϊ�������ļ��ϲ���һ���ļ������Ķ����ݺ���
void readData(string path, double *block, int i, int flag)
{
	char a[10];
	sprintf(a, "%d", i);

	string ss;
	ss.append(path);
	ss.append("/");
	ss.append(a);
	ss.append(".txt");


	ifstream file(ss.data(), ios_base::binary);
	file.seekg(my_rank*blockSize / 2 * sizeof(double));
	file.read((char*)(&block[flag*blockSize / 2]), blockSize / 2 * sizeof(double));

	file.close();
}

//�����洢�м����ĺ�������������Լ��鿴���н������ȷ��
void store(double* data, string name, int blockSize)
{
	ofstream file(name);
	for (int i = 0; i < blockSize; i++)
		file << data[i] << endl;
	file.close();
}

void check(double* data, int a, int b, int pnum)
{
	int k = a;
	if (b<a)
		return;
	while (data[k] <= data[k + 1] && k != b)
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

int findsmall(double* block, int begin, double smallvalue)
{
	int i = 0;
	while (block[begin] > smallvalue&&begin > 0)
	{
		begin--;
		i++;
	}
	return i;
}

int checkFinish(int* a, int size)
{
	for (int i = 0; i < size; i++)
		if (a[i] == 0)
			return 0;
	return 1;
}

void merge(string path, int filenum, int fileSize, int eachBlockSize)
{
	int num_perp = filenum / comm_sz;
	double smallvalue;
	double *data;
	int index = 0;

	int *begin = new int[num_perp];//ÿ���зֵ�����ʱ��ʾ����ʼλ��
	double **block = new double*[num_perp];//���������
	int *localBlockSize = new int[num_perp];//�ҵ���Сֵ�������к�ÿ�����ݿ�Ĵ�С
	int *finishFlag = new int[filenum];//ÿ������iҪ������ļ�j�Ƿ�����ɣ�finishFlag���������ʼ����
	int *localFlag = new int[num_perp];
	double *smallvalueList = new double[comm_sz];//���н��̽����е���Сֵ���б�
	double *pivot = new double[comm_sz - 1];
	int **num_block = new int*[num_perp];
	int **begin_block = new int*[num_perp];
	int **num = new int*[num_perp];
	int *total = new int[num_perp];
	int **recvDisp = new int*[num_perp];
	int *finalBlockSize = new int[comm_sz];
	int *recvDis = new int[comm_sz];
	cout << my_rank << "new����" << endl;
	for (int i = 0; i < filenum; i++)
	{
		finishFlag[i] = 0;
		localFlag[i] = 0;
	}
	for (int i = 0; i < num_perp; i++)
	{
		begin[i] = 0;
		num_block[i] = new int[comm_sz];
		begin_block[i] = new int[comm_sz];
		num[i] = new int[comm_sz];
		recvDisp[i] = new int[comm_sz];
	}
	cout << my_rank << "��ʼ" << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	do
	{
		smallvalue = (numeric_limits<double>::max)();
		for (int i = 0; i < num_perp; i++)
		{

			char a[10];
			sprintf(a, "%d", i + 1 + my_rank*num_perp);

			string ss;
			ss.append(path);
			ss.append("/middle_");
			ss.append(a);
			ss.append(".txt");

			localBlockSize[i] = (fileSize - begin[i] < eachBlockSize) ? fileSize - begin[i] : eachBlockSize;
			block[i] = new double[localBlockSize[i]];
			cout << my_rank << " blocksize of " << ss << " is " << localBlockSize[i] << endl;

			ifstream file(ss.data(), ios_base::binary);
			file.seekg(begin[i] * sizeof(double));
			file.read((char*)block[i], sizeof(double)*localBlockSize[i]);
			cout << my_rank << "������������" << endl;
			if (block[i][localBlockSize[i] - 1] < smallvalue)
				smallvalue = block[i][localBlockSize[i] - 1];
		}
		MPI_Barrier(MPI_COMM_WORLD);
		cout << my_rank << "���ж�������" << endl;
		
		MPI_Allgather(&smallvalue, 1, MPI_DOUBLE, smallvalueList, 1, MPI_DOUBLE, MPI_COMM_WORLD);
		cout << my_rank << "allgather���" << endl;
		MPI_Barrier(MPI_COMM_WORLD);
		
		int localtotal = 0;
		int *localtotalList = new int[comm_sz];
		for (int i = 0; i < num_perp; i++)
			localtotal += localBlockSize[i];
		MPI_Allgather(&localtotal, 1, MPI_INT, localtotalList, 1, MPI_INT, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		localtotal = 0;
		for (int i = 0; i < comm_sz; i++)
			localtotal += localtotalList[i];
		cout << my_rank << "�ܹ�" << localtotal << endl;
		if (localtotal < eachBlockSize*BLOCK*filenum)
			for (int i = 0; i < comm_sz - 1; i++)
				pivot[i] = (numeric_limits<double>::max)();
		else
		{
			cout << my_rank << "ѡ����С" << endl;
			smallvalue = smallvalueList[0];
			for (int i = 1; i < comm_sz; i++)
				if (smallvalueList[i] < smallvalue)
					smallvalue = smallvalueList[i];
			cout << my_rank << "ѡ����Сֵ���" << endl;
			localBlockSize[0] -= findsmall(block[0], localBlockSize[0] - 1, smallvalue);
			cout << my_rank << " blocksize of total0" << " is " << localBlockSize[0] << endl;
			for (int i = 1; i < num_perp; i++)
			{
				localBlockSize[i] -= findsmall(block[i], localBlockSize[i] - 1, smallvalue);
				cout << my_rank << " blocksize of total" << i << " is " << localBlockSize[i] << endl;
			}
			cout << my_rank << "localtotal�������" << endl;
			if (my_rank == 0)
			{
				if (comm_sz != 1)
				{
					int diff = localBlockSize[0] / comm_sz;
					cout << "diff " << diff << endl;
					pivot[0] = block[0][diff];
					cout << "block is " << pivot[0] << endl;
					cout << "pivot" << endl;
					for (int i = 1; i < comm_sz; i++)
					{
						pivot[i - 1] = block[0][i*diff];
						cout << pivot[i - 1] << " ";
					}
					cout << endl;


				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			cout << my_rank << "��Ԫѡ�����" << endl;
		}
		
		if (comm_sz != 1)
			MPI_Bcast(pivot, comm_sz - 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		cout << my_rank << "bcast���" << endl;
		for (int i = 0; i < num_perp; i++)
		{
			if (comm_sz == 1)
				num_block[i][0] = localBlockSize[i];
			else
			{
				num_block[i][comm_sz - 1] = findsmall(block[i], localBlockSize[i] - 1, pivot[comm_sz - 2]);
				num_block[i][0] = num_block[i][comm_sz - 1];
				for (int j = comm_sz - 2; j > 0; j--)
				{
					num_block[i][j] = findsmall(block[i], localBlockSize[i] - 1 - num_block[i][0], pivot[j - 1]);
					num_block[i][0] += num_block[i][j];//num_block[i][0]ͳ������block�������ٱ��������������ֵ
				}
				num_block[i][0] = localBlockSize[i] - num_block[i][0];
			}
		}
		cout << my_rank << "num_block���" << endl;
		for (int i = 0; i < num_perp; i++)
		{
			begin_block[i][0] = 0;
			for (int j = 1; j < comm_sz; j++)
				begin_block[i][j] = begin_block[i][j - 1] + num_block[i][j - 1];

		}
		cout << my_rank << "begin_block���" << endl;
		//����Ҫ���յ���Ŀ
		//num[i][j]:Ҫ����j�̵߳�i�ļ��е�������Ŀ
		
		for (int i = 0; i < num_perp; i++)
			MPI_Alltoall(num_block[i], 1, MPI_INT, num[i], 1, MPI_INT, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		cout << my_rank << "��������alltoall���" << endl;
		int alltotal = 0;
		for (int i = 0; i < num_perp; i++)
		{
			recvDisp[i][0] = 0;
			for (int j = 1; j < comm_sz; j++)
				recvDisp[i][j] = recvDisp[i][j - 1] + num[i][j - 1];
			total[i] = recvDisp[i][comm_sz - 1] + num[i][comm_sz - 1];
			alltotal += total[i];
		}
		double *new_block = new double[alltotal];

		for (int i = 0, localsize = 0; i < num_perp; i++)
		{
			MPI_Alltoallv(block[i], num_block[i], begin_block[i], MPI_DOUBLE, &new_block[localsize], num[i], recvDisp[i], MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			localsize += total[i];
		}
		cout << my_rank << "alltoallv���" << endl;
		quickSort(new_block, 0, alltotal - 1);

		

		MPI_Gather(&alltotal, 1, MPI_INT, finalBlockSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		cout << my_rank << "gather���" << endl;
		int finaltotal = 0;
		if (my_rank == 0)
		{
			recvDis[0] = 0;
			finaltotal = finalBlockSize[0];
			for (int i = 1; i < comm_sz; i++)
			{
				recvDis[i] = finalBlockSize[i - 1] + recvDis[i - 1];
				finaltotal += finalBlockSize[i];
			}
			data = new double[finaltotal];
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gatherv(new_block, alltotal, MPI_DOUBLE, data, finalBlockSize, recvDis, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		cout << my_rank << "gatherv���" << endl;
		char a[10];
		sprintf(a, "%d", index);

		string ss;
		ss.append("result");
		ss.append(a);

		if (my_rank == 0)
			store(data, ss, finaltotal);
		cout << my_rank << "�洢���" << endl;
		for (int i = 0; i < num_perp; i++)
		{
			begin[i] += localBlockSize[i];
			if (begin[i] == fileSize)
				localFlag[i] = 1;
			cout << my_rank << "begin of file " << i << " is " << begin[i] << endl;
		}
		cout << my_rank << "����finishflag��begin���" << endl;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allgather(localFlag, num_perp, MPI_INT, finishFlag, num_perp, MPI_INT, MPI_COMM_WORLD);
		cout << "finishFlag ";
		for (int i = 0; i < filenum; i++)
			cout << finishFlag[i] << " ";
		cout << endl;
		if (my_rank == 0)
			delete[] data;
		for (int i = 0; i < num_perp; i++)
			delete[] block[i];
		delete[] new_block;
		index++;
		cout << my_rank << "��һ�ֿ�ʼ" << endl;
	} while (checkFinish(finishFlag, filenum) != 1);


}