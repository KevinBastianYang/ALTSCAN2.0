#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <libgen.h>
#include <time.h>
#include "ZOE.h"


int produce_i_limit(i_limit)
{
	i_limit+=rand()%10*100000;
	if (i_limit<100000000)
	{
		return i_limit;
	}
	printf("LIMIT REACHED\n");
	return 100000000;
}

int produce_k_limit(k_limit,i_limit)
{
	k_limit+=rand()%10*100000;
	if (k_limit<i_limit)
	{
		return k_limit;
	}
	printf("KKK LIMIT REACHED\n");
	return  i_limit;
}

void InitTreeNode(zTBTreeNode** layer)
{
	int k;
	for(k=0; k < 1000; k++){
		layer[k]=zCreateTBTreeNode();
		layer[k]->score = MIN_SCORE;
	}
}
void FreeLayer(zTBTreeNode** p){
	if (p != NULL)
	{
		zFree(p);
		p = NULL;

	}
}

int main (int argc, char *argv[]) {
	int i=0;
// if (false)
// {
// 	void *p = calloc(10000,1);
//     printf( "end:%p heap:%p rodata:%p data:%p\n", &_end, p, ppp1, ppp0 );
//     sleep(10000); /* sleep to give chance to look at the process memory map */
//     return 0;
// }

	zTBTreeNode **tmp;
	tmp=zCalloc(100000,sizeof(zTBTreeNode*),"zAllocViterbi layer");
	InitTreeNode(tmp);
	sleep(10);
//	tmp=zCreateTBTreeNode();
	for (i = 0; i < 100000; ++i)
	{
		zFreeTreeNode(&tmp[i]);
	}
	sleep(10);
	zFreeLayer(&tmp);
	sleep(5);

	tmp=zCalloc(100000,sizeof(zTBTreeNode*),"zAllocViterbi layer");
	InitTreeNode(tmp);
	sleep(10);
//	tmp=zCreateTBTreeNode();
	for (i = 0; i < 100000; ++i)
	{
		zFreeTreeNode(&tmp[i]);
	}
	zFreeLayer(&tmp);



	if (false)
	{//1Billion apply memory

		srand( (unsigned)time( NULL ) ); 
		int i=0;
		int j=0;
		int k=0;
		int i_limit=0;
		int k_limit=0;
		int mark=0;
		int total_size=0;
		int little_size=0;

		int **a;
		total_size=2e8*sizeof(int*);
		a=zMalloc(total_size,"treetry try");
		printf("%p\n",a );
		printf("%p\n",&a );

		i_limit=10000000;
		little_size=10*sizeof(int);
		printf("apply %d -> %d :\n",i,i_limit,i_limit-i );
		sleep(0.5);
		for (i = 0; i < i_limit; ++i)
		{	
			
			a[i]=zMalloc(little_size,"treetry");
		//	printf("%p\n",a[i] );
		//	printf("%p\n",&(a[i]) );
			for (j = 0; j < 10; ++j)
			{
			//	printf("%p\n",a[i][j] );
			//	printf("%p\n",&(a[i][j]) );
				a[i][j]=rand()%100+1;
			}
			total_size+=little_size;
		}
		printf("TOTAL %d\n", total_size );
		while(i_limit<190000000)
		{
			i_limit=produce_i_limit(i_limit);
			printf("apply %d -> %d :\n",i,i_limit,i_limit-i );
			sleep(0.5);
			for (; i < i_limit ; ++i)
			{
				a[i]=zMalloc(little_size,"treetry");
				printf("%p\n",&(a[i]) );
				for (j = 0; j < 10; ++j)
				{
					a[i][j]=rand()%100+1;
			//		printf("%p\n",&(a[i][j]) );
				}
				total_size+=little_size;
			}
			printf("TOTAL %d\n", total_size );

			k_limit=produce_k_limit(k_limit,i_limit);
			printf("Release %d -> %d :\n",k,k_limit,k_limit-k );
			sleep(0.5);
			for (; k < k_limit; ++k)
			{
				zFree(a[k]);
				total_size-=little_size;
			}
			printf("TOTAL %d\n", total_size );
		}

	printf("----------------------------------\n");

		sleep(3);
		i_limit=200000000;
		printf("apply %d -> %d :\n",i,i_limit,i_limit-i );
		sleep(0.5);
		for (; i < i_limit ; ++i)
		{
			a[i]=zMalloc(little_size,"treetry");
			printf("%p\n",&(a[i]) );
			for (j = 0; j < 10; ++j)
			{
				a[i][j]=rand()%100+1;
			//	printf("%p\n",&(a[i][j]) );
			}
			total_size+=little_size;
		}
		printf("TOTAL %d\n", total_size );

		sleep(3);

		k_limit=100000000;
		printf("Release %d -> %d :\n",k,k_limit,k_limit-k );
		sleep(0.5);
		for (; k < k_limit; ++k)
		{
			zFree(a[k]);
			total_size-=little_size;
		}
		printf("TOTAL %d\n", total_size );

		sleep(3);
		printf("FIN\n");
}



































// 	int i,j,k;
// 	int count=20;
// 	zTBTree *t;
// 	t=malloc(sizeof(zTBTree));
// 	zTBTreeNode *layer[count];
// 	zTBTreeNode*c;
// 	zInitTBTree(t);


// 	int dna=6000;
// 	int size=49*sizeof(zTrellisCell);

// printf("start\n");
// zTrellisCell **cell;
// 	cell = zCalloc(dna, sizeof(zTrellisCell*), 
// 							"zAllocViterbi cells");
// for (i = 0; i < dna; ++i)
// {
// cell[i] = zMalloc(size, "zAllocViterbi cells[i]");





// 	for (j = 0; j < 49; ++j)
// 	{
// 		cell[i][j].layer=zCalloc(top_num,sizeof(zTBTreeNode*),"asdasd");
// 	//	zInitTreeNode(cell[i][j].layer);
// 		for(k=0; k < top_num; k++){
// 			cell[i][j].layer[k]=zCreateTBTreeNode();

// 			cell[i][j].layer[k]->score = MIN_SCORE;
// 		}
// 	}

// }

// printf("%f\n",cell[4][4].layer[4]->score );
// printf("FIN\n");
// zDeleteTBTreeNode(cell[4][4].layer[4]);

// zFree(cell[4][4].layer[4]);
// zFree(cell[4][4].layer);

// printf("%f\n",cell[4][4].layer[4]->score );



	return 0;
}