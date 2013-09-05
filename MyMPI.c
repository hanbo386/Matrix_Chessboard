/*
 *   MyMPI.c -- A library of matrix/vector
 *   input/output/redistribution functions
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 4 September 2002
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "MyMPI.h"


							/***************** 用到的函数 *****************/

/*****************************************************   计算并返回给定类型的大小   **********************************************************************/
int get_size (MPI_Datatype t) {
   if (t == MPI_BYTE) return sizeof(char);
   if (t == MPI_DOUBLE) return sizeof(double);
   if (t == MPI_FLOAT) return sizeof(float);
   if (t == MPI_INT) return sizeof(int);
   printf ("Error: Unrecognized argument to 'get_size'\n");
   fflush (stdout);
   MPI_Abort (MPI_COMM_WORLD, TYPE_ERROR);
}




/**************************************************   输出错误，结束程序    ****************************************************************************/
void terminate (
   int   id,            /* IN - Process rank */
   char *error_message) /* IN - Message to print */
{
   if (!id) {
      printf ("Error: %s\n", error_message);
      fflush (stdout);
   }
   MPI_Finalize();
   exit (-1);
}




/*************************************   在某一通信域内，0号进程读取相同列向量，并使该通信域内每个一 个进程都拥有其副本   ****************************************/
void read_replicated_vector (
   char        *s,      /* IN - File name */
   void       **v,      /* OUT - Vector */
   MPI_Datatype dtype,  /* IN - Vector type */
   int         *n,      /* OUT - Vector length */
   MPI_Comm     comm)   /* IN - Communicator */
{
   int        datum_size; /* Bytes per vector element */
   int        i;
   int        id;         /* Process rank */
   FILE      *infileptr;  /* Input file pointer */
   int        p;          /* Number of processes */

   MPI_Comm_rank (comm, &id);
   MPI_Comm_size (comm, &p);
   datum_size = get_size (dtype);
   if (id == (p-1)) {
      infileptr = fopen (s, "r");
      if (infileptr == NULL) *n = 0;
      else fread (n, sizeof(int), 1, infileptr);
   }
   MPI_Bcast (n, 1, MPI_INT, p-1, MPI_COMM_WORLD);
   if (! *n) terminate (id, "Cannot open vector file");

   *v = (int * )malloc (*n * datum_size);

   if (id == (p-1)) {
      fread (*v, datum_size, *n, infileptr);
      fclose (infileptr);
   }
   MPI_Bcast (*v, *n, dtype, p-1, MPI_COMM_WORLD);
}




/********************************************************   打印一行向量   *************************************************************************/
void print_subvector (
   void        *a,       /* IN - Array pointer */
   MPI_Datatype dtype,   /* IN - Array type */
   int          n)       /* IN - Array size */
{
   int i;
   for (i = 0; i < n; i++) {
      if (dtype == MPI_DOUBLE)
         printf ("%6.3f ", ((double *)a)[i]);
      else {
         if (dtype == MPI_FLOAT)
            printf ("%6.3f ", ((float *)a)[i]);
         else if (dtype == MPI_INT)
            printf ("%6d ", ((int *)a)[i]);
      }
   }
}




/*******************************************************   由0号进程打印向量   *************************************************************/
void print_replicated_vector (
   void        *v,      /* IN - Address of vector */
   MPI_Datatype dtype,  /* IN - Vector element type */
   int          n,      /* IN - Elements in vector */
   MPI_Comm     comm)   /* IN - Communicator */
{
   int id;              /* Process rank */

   MPI_Comm_rank (comm, &id);
   
   if (!id) {
      print_subvector (v, dtype, n);
      printf ("\n\n");
   }
}


/*****************************************************   按照棋盘方式读取矩阵   *************************************************************/ 
void read_checkerboard_matrix (						//从文件读取一个矩阵
   char *s,              /* IN - File name */
   void ***subs,         /* OUT - 2D array */
   void **storage,       /* OUT - Array elements */
   MPI_Datatype dtype,   /* IN - Element type */
   int *m,               /* OUT - Array rows */
   int *n,               /* OUT - Array cols */
   MPI_Comm grid_comm)   /* IN - Communicator */
{
   void      *buffer;         /* File buffer */
   int        coords[2];      /* Coords of proc receiving	//?????????????????????
                                 next row of matrix */
   int        datum_size;     /* Bytes per elements */
   int        dest_id;        /* Rank of receiving proc */
   int        grid_coord[2];  /* Process coords */
   int        grid_id;        /* Process rank */
   int        grid_period[2]; /* Wraparound */
   int        grid_size[2];   /* Dimensions of grid */
   int        i, j, k;
   FILE      *infileptr;      /* Input file pointer */
   void      *laddr;          /* Used when proc 0 gets row */
   int        local_cols;     /* Matrix cols on this proc */
   int        local_rows;     /* Matrix rows on this proc */
   void     **lptr;           /* Pointer into 'subs' */
   int        p;              /* Number of processes */
   void      *raddr;          /* Address of first element
                                 to send */
   void      *rptr;           /* Pointer into 'storage' */
   MPI_Status status;         /* Results of read */

   MPI_Comm_rank (grid_comm, &grid_id);
   MPI_Comm_size (grid_comm, &p);
   datum_size = get_size (dtype);

   /* Process 0 opens file, gets number of rows and
      number of cols, and broadcasts this information
      to the other processes. */

   if (grid_id == 0) {				//0号进程负责打开文件，读取文件
      infileptr = fopen (s, "r");
      if (infileptr == NULL) *m = 0;
      else {
         fread (m, sizeof(int), 1, infileptr);
         fread (n, sizeof(int), 1, infileptr);
      }
   }
   MPI_Bcast (m, 1, MPI_INT, 0, grid_comm);

   if (!(*m)) MPI_Abort (MPI_COMM_WORLD, OPEN_FILE_ERROR);  //这里到底如何保证进程同步的？？？

   MPI_Bcast (n, 1, MPI_INT, 0, grid_comm);

   //printf("fuck %d \n",*m);
   //printf("fuck %d \n",*n); 
   /* Each process determines the size of the submatrix
      it is responsible for. */

   MPI_Cart_get (grid_comm, 2, grid_size, grid_period,		//通信域是笛卡尔拓扑通信域 //grid_size返回的是各维度的进程数
      grid_coord);						//grid_coord返回的是当前进程所在网格内坐标
   local_rows = BLOCK_SIZE(grid_coord[0],grid_size[0],*m);	//id,p,n
   local_cols = BLOCK_SIZE(grid_coord[1],grid_size[1],*n);
   //记录本进程子矩阵的行数和列数

   /* Dynamically allocate two-dimensional matrix 'subs' */
   *storage = (int * )malloc (local_rows * local_cols * datum_size);
   *subs = (void **) malloc (local_rows*PTR_SIZE);
   lptr = (void *) *subs;
   rptr = (void *) *storage;
   for (i = 0; i < local_rows; i++) {
      *(lptr++) = (void *) rptr;
      rptr += local_cols * datum_size;
   }
   
   /* Grid process 0 reads in the matrix one row at a time
      and distributes each row among the MPI processes. */
   if (grid_id == 0)
      buffer = (int * )malloc(*n * datum_size);

   /* For each row of processes in the process grid... 这里一共三重循环，每一行进程组，每一进程行，每一进程行的进程列*/
   for (i = 0; i < grid_size[0]; i++) {				//grid_size[0]是行数，grid_size[1]是列数
      	coords[0] = i;

      	/* For each matrix row controlled by this proc row...*/
      	for (j = 0; j < BLOCK_SIZE(i,grid_size[0],*m); j++) 
	{ 	//每一行进程组，根据其行数，进行循环

        	 /* Read in a row of the matrix */

         	if (grid_id == 0) 
		{
            		fread (buffer, datum_size, *n, infileptr);   				//0号进程读入原始矩阵的一行
         	}

        	 /* Distribute it among process in the grid row */
		/*在当前大行进行分发*/
         	for (k = 0; k < grid_size[1]; k++) 
		{			//每一大行，根据其列数，进行循环
            		coords[1] = k;

            		/* Find address of first element to send */
            		raddr = buffer + BLOCK_LOW(k,grid_size[1],*n) * datum_size;		//获取要分发数据的首地址

            		/* Determine the grid ID of the process getting the subrow */
            		MPI_Cart_rank (grid_comm, coords, &dest_id);				//coords记录了指定进程在虚拟网格中的坐标； dest_id返回该进程的进程号

            		/* Process 0 is responsible for sending...*/
            		if (grid_id == 0) 
			{	//依然由0号进程负责发送

               			/* It is sending (copying) to itself */
               			if (dest_id == 0) 
				{
                  			laddr = (*subs)[j];
                  			memcpy (laddr, raddr, local_cols * datum_size);		//目的地址，源地址，数据大小； 根据该进程列数取数据

               				/* It is sending to another process */
               			} 
				else 
				{
                  			MPI_Send (raddr, BLOCK_SIZE(k,grid_size[1],*n), dtype, dest_id, 0, grid_comm);
               			}

            			/* Process 'dest_id' is responsible for receiving... */
            		} 
			else if (grid_id == dest_id) 
			{									//其它进程等待属于自己的数据
               			MPI_Recv ((*subs)[j], local_cols, dtype, 0, 0, grid_comm,&status);
            		}
         	}
      	}
   }
  if (grid_id == 0) free (buffer);				//最后释放0号进程分配的临时内存
}



/************************************************   本程序的核心函数，功能不描述了，主函数里有说明   *****************************************************************************/
void scatter_vector(char *s, int **a, int *buff, int *reduce_buff, int *m, int *n, MPI_Datatype dtype, MPI_Status status, MPI_Comm grid_comm)
{
	int 	i,j,k,l,r;		//循环控变量
	int 	id;			//进程编号
	int 	p;			//进程数
	int	vector_len;		//记录向量长度
	int     grid_coord[2];  	//记录进程在网格内坐标
   	int     grid_id;        	//
   	int     grid_period[2]; 	//
   	int     grid_size[2];  		//记录各维度进程数
	int 	coords[2]; 		//用于保存当前循环下坐标
	int 	coords_temp[2];		//
	int	coords_self[2];
	int 	current_id;		//记录循环中当前坐标对应id	
	int	temp_id;		//
	int 	datum_size;		//变量大小
	int	local_rows;		//本进程子矩阵行数
	int	local_cols;		//本进程子矩阵列数
	int	*raddr;			//0进程分发向量时用于记录分发位置的指针
	int	*vector_ptr;		//指向列向量的指针
	int	*vector;		//指向每个进程执行矩阵乘法后得到的向量指针
	MPI_Comm reduce_comm;		//进行按行归约时划分的通信域
	MPI_Comm broad_comm;		//进行按列广播时划分的通信域
	MPI_Comm final_comm;		//最后用于打印列向量的通信域
	FILE *fp;			//文件指针
	
	MPI_Comm_rank (grid_comm, &id);
	MPI_Comm_size (grid_comm, &p);
	datum_size = get_size (dtype);


	MPI_Cart_get (grid_comm, 2, grid_size, grid_period, grid_coord);//通信域是笛卡尔拓扑通信域 ; grid_coord返回的是当前进程所在网格内坐标; grid_size返回的是各维度的进程数
	MPI_Cart_coords(grid_comm, id, 2, coords_self);			//获得本进程所在坐标

	local_rows = BLOCK_SIZE(grid_coord[0],grid_size[0],*m);		//id,p,n ;本进程子矩阵行数和列数
	local_cols = BLOCK_SIZE(grid_coord[1],grid_size[1],*n);
	
	vector = (int * )malloc(local_rows * datum_size); 		//保存本进程执行一次矩阵乘法之后的向量
	buff = (int *)malloc(local_cols * datum_size);			//为本进程分配内存,以存储本进程分配的向量部分

	/*读向量文件的部分，先读取向量长度，然后根据此长度一次取出整个向量*/
	if (!id) 							//0号进程读取向量文件
	{
		fp = fopen (s, "r");
	      	if (fp == NULL) vector_len = 0;
	      	else fread (&vector_len, sizeof(int), 1, fp);
  	}
	MPI_Bcast (&vector_len, 1, MPI_INT, 0, grid_comm);	
	if (! vector_len) terminate (id, "Cannot open vector file");
	if(!id)
	{
		vector_ptr = (int *)malloc(datum_size * vector_len);
		fread(vector_ptr, sizeof(int), vector_len, fp);
		fclose(fp);
	}

	/*在网格通信域基础上，创建需要的以列和以行划分的各自的通信域*/
	MPI_Comm_split(grid_comm, coords_self[1], coords_self[0], &broad_comm);		//按列划分广播通信域,以及最后打印结果的通信域
	int broad_id;
	int broad_p;				
	MPI_Comm_rank(broad_comm, &broad_id);
	MPI_Comm_size(broad_comm, &broad_p);

	MPI_Comm_split(grid_comm, coords_self[0], coords_self[1], &reduce_comm);	//按行划分归约通信域
	int 	reduce_p;
	int 	reduce_id;
	MPI_Comm_rank(reduce_comm, &reduce_id);
	MPI_Comm_size(reduce_comm, &reduce_p);
	reduce_buff = (int *)malloc(local_rows * datum_size);

	for(i = 0; i < grid_size[0]; i++)						//每个进程执行此循环//网格行
	{		
		coords[0] = i;								//记录当前横坐标
		for(j = 0; j < grid_size[1]; j++)					//网格列
		{	
			coords[1] = j;							//记录当前纵坐标
			MPI_Cart_rank(grid_comm, coords, &current_id);			//根据记录的坐标计算出

			if(current_id == id)						//第一行内进程接收属于自己的进程，这个条件必须有
			{
				if(id == 0)						//0号进程为本行进程分发向量
				{
					for(k = 0; k < grid_size[1]; k++)		//分发向量的循环过程
					{
						raddr = vector_ptr + BLOCK_LOW(k, grid_size[1], *n);		//根据各块大小移动指针
						coords_temp[0] = 0;
						coords_temp[1] = k;
						MPI_Cart_rank(grid_comm, coords_temp, &temp_id);		//获得对应坐标的id
						if(!temp_id)							//如果是0号进程，则不能用send给自己发送数据，所以单独处理
						{
							memcpy(buff, raddr, local_cols * datum_size);
									
						}
						else								//为0号进程同一行的其它进程发送该属于进程的部分向量
						{	
							MPI_Send (raddr, BLOCK_SIZE(k, grid_size[1], *n), dtype, temp_id, 0, grid_comm);
							
						}
					}
				}
				if(coords[0] == 0 && id)							//处于第一行的非0号进程，接收来自0号进程的数据
					MPI_Recv(buff, BLOCK_SIZE(j, grid_size[1], *n), dtype, 0, 0, grid_comm, &status);

				/*至此第一行的进程都有了自己的部分向量，接下来需要把自己的数据广播至同列的其它进程；这里注意：Bcast和Reduce都需要指定通信域内所有进程调用*/
				MPI_Bcast(buff, BLOCK_SIZE(j, grid_size[1], *n), dtype, 0, broad_comm);		//MPI_Bcast 和 MPI_Reduce 可能会导致新的进程域不能划分，最初出现过
														//上述错误，没有进一步验证，有待考究


				for(l = 0; l < local_rows; l++)							//各进程执行矩阵乘法
				{
					vector[l] = 0;
					for(r = 0; r < local_cols; r++)
					{
						vector[l] = buff[r] * a[l][r] + vector[l];
					}
				}
			}
		}	
	}

	/*以上是列向量分发以及进行各进程内子矩阵运算的过程，下面是数据的收集过程*/
	MPI_Reduce(vector, reduce_buff, local_rows, dtype, MPI_SUM, 0, reduce_comm);//发送消息缓冲区地址；接收消息缓冲区地址；发送消息缓冲区数据个数...
	/*按照进程网格各行进行归约求和。这里可能犯错，一定要弄清楚第三个参数*/
	
	int *print_buff;			//定义一个指向最终结果向量的指针，其大小为原始矩阵行数
	print_buff = (int *)malloc(*m * datum_size);

	replicate_block_vector(reduce_buff, *m, print_buff, dtype, broad_comm);		//作者给出的一个从指定通信域内各进程收集数据并将最终结果复制给该通信域内所有进程的函数
											//这里根据实际应用需要，把其中的MPI_Allgatgerv改为MPI_Gatherv，使得最终只有0号进
											//程保留并打印该最终结果
	
	/*下面由0号进程打印最终结果*/
	if(!id)
	{
		printf("最终计算结果向量为：\n");
		for(i = 0; i < *m; i++)
			printf("%d   ", print_buff[i]);
		printf("\n");
	}
	/*
	MPI_Comm_free(&reduce_comm);	
	MPI_Comm_free(&final_comm);
	MPI_Comm_free(&broad_comm);
	MPI_Comm_free(&grid_comm);*/
	
}



/*********************************************************************   打印棋盘矩阵   ********************************************************************************/
void print_checkerboard_matrix (					//打印一个矩阵
   void       **a,            /* IN -2D matrix */
   MPI_Datatype dtype,        /* IN -Matrix element type */
   int          m,            /* IN -Matrix rows */
   int          n,            /* IN -Matrix columns */
   MPI_Comm     grid_comm)    /* IN - Communicator */
{
   void      *buffer;         /* Room to hold 1 matrix row */
   int        coords[2];      /* Grid coords of process
                                 sending elements */
   int        datum_size;     /* Bytes per matrix element */
   int        els;            /* Elements received */
   int        grid_coords[2]; /* Coords of this process */
   int        grid_id;        /* Process rank in grid */
   int        grid_period[2]; /* Wraparound */
   int        grid_size[2];   /* Dims of process grid */
   int        i, j, k;
   void      *laddr;          /* Where to put subrow */
   int        local_cols;     /* Matrix cols on this proc */
   int        p;              /* Number of processes */
   int        src;            /* ID of proc with subrow */
   MPI_Status status;         /* Result of receive */

   MPI_Comm_rank (grid_comm, &grid_id);
   MPI_Comm_size (grid_comm, &p);
   datum_size = get_size (dtype);

   MPI_Cart_get (grid_comm, 2, grid_size, grid_period,
      grid_coords);
   local_cols = BLOCK_SIZE(grid_coords[1],grid_size[1],n);

  	if (!grid_id)
      	buffer = (int *)malloc(n*datum_size);

   	/* For each row of the process grid */
  for (i = 0; i < grid_size[0]; i++) 
  {
      	coords[0] = i;

      		/* For each matrix row controlled by the process row */
      	for (j = 0; j < BLOCK_SIZE(i,grid_size[0],m); j++) 
	{

         	/* Collect the matrix row on grid process 0 and
            	print it */
         	if (!grid_id) 			//0号进程负责打印一行数据
		{
            		for (k = 0; k < grid_size[1]; k++) 
			{
               			coords[1] = k;
               			MPI_Cart_rank (grid_comm, coords, &src);
               			els = BLOCK_SIZE(k,grid_size[1],n);
               			laddr = buffer + BLOCK_LOW(k,grid_size[1],n) * datum_size;
               			if (src == 0) 
				{
                  			memcpy (laddr, a[j], els * datum_size);
               			} 
				else 
				{
                  			MPI_Recv(laddr, els, dtype, src, 0, grid_comm, &status);
               			}
            		}
           	 	print_subvector (buffer, dtype, n);
            		putchar ('\n');
        	} 
		else if (grid_coords[0] == i) 
		{
            		MPI_Send (a[j], local_cols, dtype, 0, 0, grid_comm);
         	}
      	}
   }
   	if (!grid_id) {
      free (buffer);
      putchar ('\n');
   }
}



/*****************************************   当从通信域内各进程需要收集数量不等的数据时，调用该函数获得每个进程的数据个数和偏移量数组   *****************************************/
void create_mixed_xfer_arrays (
   	int id,       /* IN - Process rank */
   	int p,        /* IN - Number of processes */
	int n,        /* IN - Total number of elements */
   	int **count,  /* OUT - Array of counts */
   	int **disp)   /* OUT - Array of displacements */
{

   	int i;
	
	
	//printf(" p is %d \n", p);
	//printf("n is %d \n", n);
   	*count = (int *)malloc (p * sizeof(int));
   	*disp = (int *)malloc (p * sizeof(int));
   	(*count)[0] = BLOCK_SIZE(0,p,n);
   	(*disp)[0] = 0;
	
   	for (i = 1; i < p; i++) 
	{
      		(*disp)[i] = (*disp)[i-1] + (*count)[i-1];
      		(*count)[i] = BLOCK_SIZE(i,p,n);
   	}
}



/*************************************   从通信域各个进程收集数量不等的数据，并将收集结果复制给所有进程；本程序修改为仅有0号进程保留此结果   ****************************************/
void replicate_block_vector (
   	void        *ablock,  /* IN - Block-distributed vector */
   	int          n,       /* IN - Elements in vector */
   	void        *arep,    /* OUT - Replicated vector */
   	MPI_Datatype dtype,   /* IN - Element type */
   	MPI_Comm     comm)    /* IN - Communicator */
{
   	int *cnt;  /* Elements contributed by each process */
   	int *disp; /* Displacement in concatenated array */
   	int id;    /* Process id */
   	int p;     /* Processes in communicator */

   	MPI_Comm_size (comm, &p);
   	MPI_Comm_rank (comm, &id);
	
   	create_mixed_xfer_arrays (id, p, n, &cnt, &disp);
   	MPI_Gatherv (ablock, cnt[id], dtype, arep, cnt, disp, dtype, 0, comm); 		//这里稍作修改
   	free (cnt);
   	free (disp);
}
