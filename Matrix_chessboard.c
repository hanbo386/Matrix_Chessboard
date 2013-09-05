#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include "MyMPI.h"
#include "MyMPI.c"

#define MPI_TYPE MPI_INT
//typedef int dtype;		//定义本程序设计到的数据类型

int main(int argc, char *argv[])
{
	int ** a;		//矩阵
	int * b;		//向量
	int * storage;		
	int * buff;
	int * reduce_buff;
	int nprime;		//向量中元素
	int id;			//进程编号
	int p;			//进程数
	int m;			//矩阵行数
	int n;			//矩阵列数
	int size[2];		//网格中每一维度的大小
	int periodic[2];
	MPI_Comm grid_comm;
	MPI_Status status;

	size[0] = size[1] = 0;
	periodic[0] = periodic[1] = 0;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	/*读取列向量，打印至屏幕方便查看*/
	read_replicated_vector(argv[1], (void * ) &b, MPI_TYPE, &nprime, MPI_COMM_WORLD);	//在某一通信域内，所有进程读取相同列向量//其实在这个题目里，仅需要0号进程操作列向量读和打印就够了
	if(!id)
		printf("程序从文件中读取的列向量：\n");
	print_replicated_vector(b, MPI_TYPE, nprime, MPI_COMM_WORLD);				//但是既然有现成的，那就懒得自己写了，这样效率可能不会高
	
	MPI_Dims_create(p, 2, size);								//创建虚拟进程网格
	MPI_Cart_create(MPI_COMM_WORLD, 2, size, periodic, 1, &grid_comm); 			//为上面拓扑建立通信域						
	read_checkerboard_matrix(argv[2], (void *) &a, (void * ) &storage, MPI_TYPE, &m, &n, grid_comm);	//读取矩阵
	if(!id)
		printf("程序从文件中读取的矩阵：\n");
	print_checkerboard_matrix((void ** )a, MPI_TYPE, m, n, grid_comm);					//打印矩阵，方便查看

	/*事实上下面这个函数总共包含3个部分，并不只是函数名所写的那样，它包含了：读取列向量并分发至进程网格；各进程进行子矩阵乘法；打印最终向量至屏幕。。。只是开始那么写了，就依旧延续名字了*/
	scatter_vector(argv[1], (int ** ) a, buff, reduce_buff, &m, &n, MPI_TYPE, status, grid_comm);	
	MPI_Finalize();		//并行计算报错很烦，前面有些函数执行不正确，它报这个函数的错
	return 0;
}
/*程序总结： 	1. 基本上好多函数都是要参与的所有进程都调用的，并不是之前想的那样，一个进程广播一组数据，其它进程就不调用Bcast函数了
           	2. 调试程序之前优先查看各MPI函数参数是否都传递正确
		3. 刚开始懒了，现在只限定到整数类型的矩阵了，不过问题不不大，下次再说吧
		4. 划分的各个通信域之间好像不怎么干扰，开始的时候理解错了
		5. 347行那个指针的移动没弄明白，再回头看看书吧
		6. 进程不要给自己Send数据
		7. 现在对各个进程同时分配的内存空间比较迷惑
		8. 感觉搞并行的人好少，出了问题都不好搜啊
		9. MPI_Comm_free不知道有什么用，没释放好像也没出什么问题
*/

