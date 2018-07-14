#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctime>
#include "unistd.h"

#define H 1
#define BLOCK 1024
double f()
{
	return double(rand()%2000000);
}
double data[BLOCK];
int main(int argc, char* argv[])
{
	int loop;
	int i, j;
	if (argc != 3)
	{
		printf("Usage:%s start_value loop_times data_file\n", argv[0]);
	}
	loop = atoi(argv[1]);
	int handle;
	printf("loop=%d\n", loop);
	FILE *fid;
	fid = fopen(argv[2], "wb");
	if (fid == NULL)
	{
		printf("写出文件出错");
		return -1;
	}
	srand((unsigned)time(NULL));
	for (i = 0; i<loop; i++)
	{
		for (j = 0; j < BLOCK; j++)
			data[j] = f();
		fwrite(data, sizeof(double), BLOCK, fid);

	}
	fclose(fid);
}
