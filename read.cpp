#include<iostream>
#include<fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace std;

int BLOCK;
int loop;
int size;

int main(int argc, char* argv[])
{
	BLOCK = atoi(argv[2]);
	loop = atoi(argv[3]);
	size = atoi(argv[4]);
	double *data;
	data = (double*)malloc(sizeof(double)*loop*BLOCK);
	ifstream file(argv[1], ios_base::binary);
	file.read((char*)data, sizeof(double)*loop*BLOCK);
	int order=0;
	int start=0,end=0;
	cout<<"������ָ����в�����"<<endl;
	cout<<"-1:�����������"<<endl;
	cout<<"-2:������ʼ���յ�λ��"<<endl;
	cout<<"-3:��ֹ����"<<endl;
	cout<<"��������:���ʸ�λ������"<<endl;
	while(cin>>order)
	{
		switch (order){
			case -1:
			{
				cout<<"����������С��"<<loop*BLOCK<<endl;
				break;
			}
			case -2:
			{
				cin>>start>>end;
				for(int i=start;i<=end;i++)
					cout<<"i"<<i<<"��"<<data[i]<<endl;
				break;
			}
			case -3:
			{
				break;
			}
			default:
			{
				cout<<"order "<<order<<"��"<<data[order]<<endl;
			}
		}
		if(order==-3)
			break;
		
	}

	return 0;
}