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
	cout<<"请输入指令进行操作："<<endl;
	cout<<"-1:输出数据总量"<<endl;
	cout<<"-2:输入起始和终点位置"<<endl;
	cout<<"-3:终止输入"<<endl;
	cout<<"其他数字:访问该位置数据"<<endl;
	while(cin>>order)
	{
		switch (order){
			case -1:
			{
				cout<<"数据总量大小："<<loop*BLOCK<<endl;
				break;
			}
			case -2:
			{
				cin>>start>>end;
				for(int i=start;i<=end;i++)
					cout<<"i"<<i<<"："<<data[i]<<endl;
				break;
			}
			case -3:
			{
				break;
			}
			default:
			{
				cout<<"order "<<order<<"："<<data[order]<<endl;
			}
		}
		if(order==-3)
			break;
		
	}

	return 0;
}