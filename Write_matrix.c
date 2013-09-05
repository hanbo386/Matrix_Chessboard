#include <stdio.h>

typedef int dtype;

int main(int argc, char *argv[])
{
	int m,n;		     	//j表示输入矩阵行数，k表示输入矩阵列数
	int i,j;	    		//循环控制变量
	FILE *fp;
	int a;
	//m = *argv[2] - 48;
	//n = *argv[3] - 48;
	
	if((fp = fopen(argv[1],"wb")) == NULL)
	{
		printf("Cannot open file!");
		//这里还不知道添加什么退出语句 exit(1);
	}				//打开指定文件

	printf("请输入矩阵行数:");
	scanf("%d", &m);
	if(m > 1)
		fwrite(&m,sizeof(int),1,fp);
	printf("请输入矩阵列数:");
	scanf("%d", &n);
	fwrite(&n,sizeof(int),1,fp);	//读入矩阵行数和列数
				
	for(i = 0; i < m; i++)
	{
		for(j = 0; j < n; j++)
		{
			scanf("%d", &a);
			fwrite(&a,sizeof(int),1,fp);
		}
	}				//用户按行输入矩阵元素，矩阵元素以2进制形式存储
	return 0;	
}
