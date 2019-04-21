/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
  zTrellis.c - part of the ZOE library for genomic analysis
  
  Copyright (C) 2001-2002 Ian F. Korf
  Modified by Zhiqiang Hu in 2012
\******************************************************************************/

#ifndef ZOE_TRELLIS_C
#define ZOE_TRELLIS_C

#include <assert.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>



#include "zTrellis.h"
#include "zFastaFile.h"
//Added by Jim
#include "zTBTree.h"
#include "zGTF.h"
#include <math.h>

/* static const int PADDING = 50; */
static const score_t MIN_INIT_SCORE = -10000;
extern coor_t any_new_cpoint;
//Added by Jim
static const int WINDOW=300;
static const int CHUNK=5000;//default 18000
static const int CPOINT_WINDOW=18000;//origin 18000
static const score_t EPISLON=1e-10;//threshold for handshake checking

static const int LARGE_TOP=1000;

//Added by Jim

void zFreeLayer(zTBTreeNode*** p){
	if (*p != NULL)
	{
	//	printf("%d\n",*p );
		zFree(*p);
	//	printf("%d\n",*p );
	//	printf("FLEEE\n");
//		
		*p = NULL;

	}


}

void zFreeViterbiVarsSingleNew(zTrellis *trellis, coor_t i) {
	int j,k;
//	if (i>600){printf("FreeViterSingle %d\n",i );}

	if (NULL == trellis->cell) return;
	if (trellis->cell[i]==NULL) return;
	
	if (i%100000==0){printf("FreeViterSingleNew %d\n",i );}

	for (j = 0; j < 49; ++j)
	{
		if (trellis->cell[i][j].layer==NULL){continue;}
		
		for (k = 0; k < top_num; ++k)
		{
			if (trellis->cell[i][j].layer[k]!=NULL)
			{
		//		printf("FreeNode %d,%d,%d\n",i,j,k );
				zFreeTreeNode(&trellis->cell[i][j].layer[k]);
			}
		}
	//	printf("FreeLayer %d,%d\n",i,j );
		zFreeLayer(&trellis->cell[i][j].layer);
	//	if (i>600){printf("FreeLayer %d\n",i );}
		

		/*/////////////////kevin
		if((trellis->hmm->state[j].type == EXPLICIT)
					|| (trellis->hmm->state[j].type == INTERNAL))
		{
			if (trellis->cell[i][j].submax!=NULL)
			{
				if (trellis->cell[i][j].submax->elem != NULL){
					zFreeSFVec(trellis->cell[i][j].submax);
					trellis->cell[i][j].submax->elem = NULL;
				}
				zFree(trellis->cell[i][j].submax);
				trellis->cell[i][j].submax = NULL;
			}
		//	if (i>600){printf("FreeSubmax %d\n",i );}
		}*/
	}
//	printf("FreeCell %d\n",i );
	zFree(trellis->cell[i]);
	trellis->cell[i] = NULL;
}

void zFreeViterbiVarsRange(zTrellis *trellis, coor_t lower,coor_t upper) {//delete not include upper
	int i;
	printf("RANGEguy %d - %d \n",lower,upper );
	for (i = lower; i < upper; ++i)
	{
		zFreeViterbiVarsSingleNew(trellis,i);
	}
}




void zFreeViterbiVarsSingle(zTrellis *trellis, coor_t lower,coor_t upper) {
	int i;
	int j,k;
//	printf("PPP: %d-%d\n",lower,upper);
	printf("FreeViterSingle %d-%d\n",lower,upper );
	if(lower==upper){printf("SAME: %d-%d\n",lower,upper);return;}
	if (lower > upper ||  upper >=trellis->dna->length)
	{
		printf("RRR: %d-%d\n",lower,upper);
		zDie("zFreeViterbiVarsSingle range error");
		return;
	}


	if (NULL == trellis->cell) return;
//	if (NULL == trellis->cell[upper-1]) return;
	for (i = lower; i < upper; ++i)
	{
	//	if (NULL == trellis->cell[i]) {printf("cell[%d] already null\n",i );continue;}
		for (j = 0; j < 49; ++j)
		{
			for (k = 0; k < top_num; ++k)
			{
				if (trellis->cell[i][j].layer[k]!=NULL)
				{
		//			printf("OMG %d,%d,%d\n",i,j,k );
					zFreeTreeNode(&trellis->cell[i][j].layer[k]);
	
				}
			}
			if (trellis->cell[i][j].layer!=NULL)
			{
	//			printf("Xuewen %d,%d\n",i,j );
		//		zFree(&trellis->cell[i][j].layer);
			//	 printf("%d\n",trellis->cell[i][j].layer );
				zFreeLayer(&trellis->cell[i][j].layer);
	//			printf("XXsuc\n");
	//			if (trellis->cell[i][j].layer==NULL)printf("CCCOOONNNGGG\n");
	//			else printf("%d\n",trellis->cell[i][j].layer );

			}

		}
		zFree(trellis->cell[i]);
		trellis->cell[i] = NULL;
		// if (trellis->cell[i]==NULL)
		// {
		// 	printf("Jinyuan\n");
		// }
		// if (trellis->cell[i][1].layer==NULL)
		// {
		// 	printf("hahahaha\n");
		// }
	}

}
void zFreeViterbiVarsFinal(zTrellis *trellis) {
	int    i;
	int j;
	int	k;

	if (NULL == trellis->cell) return;
	for (i = trellis->dna->length-CHUNK-CHUNK; i < trellis->dna->length; i++) {
//	for (i = 0; i < trellis->hmm->states; i++) {
		if (NULL == trellis->cell[i]) continue;
 		//Added by Jim
		for (j = 0; j < trellis->hmm->states; j++)
		{
			for (k = 0; k < top_num; k++)
			{
				if (trellis->cell[i][j].layer[k]!=NULL)
				{
					zFreeTreeNode(&trellis->cell[i][j].layer[k]);
				}
							
			}
			if (trellis->cell[i][j].layer!=NULL)
			{
	//			printf("Xuewen %d,%d\n",i,j );
	//			zFree(&trellis->cell[i][j].layer);
			//	 printf("%d\n",trellis->cell[i][j].layer );
				zFreeLayer(&trellis->cell[i][j].layer);
	//			printf("XXsuc\n");
	//			if (trellis->cell[i][j].layer==NULL)printf("CCCOOONNNGGG\n");
	//			else printf("%d\n",trellis->cell[i][j].layer );

			}

		//	zFree(trellis->cell[i][j].layer);
		}

		zFree(trellis->cell[i]);
		trellis->cell[i] = NULL;

		
	}
	for (i = 0; i < trellis->hmm->states; i++) 
	{
		zFreeIVec(&trellis->extpos[i]);
	}
	
	zFree(trellis->cell);
	trellis->cell = NULL;

	zFree(trellis->extpos);
	trellis->extpos = NULL;
}

//Added by Jim
void zTracePartialTrellisSingle(zTrellis *trellis, zTBTreeNode *n,int state, int story,zTBTree *tree)
{
	if (n->sfv!=NULL)
	{
		return;
	}

	zHMM *hmm = trellis->hmm;
	int trace;
	int trace_top;
	zSfeature istate;
	int tmp_state = state;
	int tmp_story = story;
	


	n->sfv=zMalloc(sizeof(zSFVec),"zMalloc sfv[i] in zTracePartialTrellisSingle");
	zInitSFVec(n->sfv,10);

	zTBTreeNode *c;
	c = n;
//	printf("INSIDE\n");
	int mark=0;
	c->left_pos=1;//will be deleted soon only for test specific function
	while (c!=tree->root){
		c->left_pos=1;//will be deleted soon only for test specific function
//		printf("LOOP: %d\n",mark);
		mark++;
	   	trace = c->trace;
		trace_top = c->trace_top;
		

		istate.name   = hmm->state[tmp_state].name;
		istate.start  = c->pos - c->length + 1;

		// printf("FAKE: %d , REAL:%d",istate.start, c->parent->pos+1);
		// if (istate.start-c->parent->pos-1!=0)
		// {
		// 	printf("FAIL\n");
		// }
		// else{printf("\n");}

		istate.end    = c->pos;
		istate.strand = hmm->state[tmp_state].strand;
		zChar2FragFrame(c->frag_frame, &istate);
	
		istate.state  = tmp_state;
		istate.group  = NULL;

		//protection revealed soonnnn!
		if (c->parent==NULL)
		{
			printf("TTTTT9\n");//LETS DO IT
			printf("HANDSHAKE self: %d,%s,%d,%d\n",c->length,zStrIdx2Char(trellis->hmm->state[state].name),istate.end,c->story);
			printf("HANDSHAKE father: %s,%d,%d\n",zStrIdx2Char(trellis->hmm->state[trace].name),istate.start-1,trace_top);
			zHandShake(trellis,c->state,c->pos,c->story);
		}
	
		if (istate.start == PADDING)
		{
			istate.score = c->score
				- trellis->cell[istate.start-1][trace].layer[trace_top]->score;
				- zGetTransitionScore(trellis->hmm,trace,tmp_state,trellis->tiso_group);
		}
		else if (c->parent==tree->root)
		{
			if(c->sfv!=NULL)
			{
				zMergeSFVec(c->sfv,n->sfv);
				printf("C:%d,N:%d\n",c->sfv->size,n->sfv->size );
			}
			return;
		}
		else
		{
			if(trace_top == -1) zDie("trace_top error!\n");
			istate.score = c->score
				- c->parent->score
				- zGetTransitionScore(trellis->hmm,trace,tmp_state,trellis->tiso_group);
		}

		/*EVAN now removing transition scores from output */
		zPushSFVec(n->sfv, &istate);//Notice: here should be n->sfv
//	printf("%d, %s,%d,%d\n",c->length,zStrIdx2Char(trellis->hmm->state[c->state].name) ,c->pos,c->story);
//	printf("\t %d, %s,%d,%d\n",c->parent->length,zStrIdx2Char(trellis->hmm->state[c->parent->state].name) ,c->parent->pos,c->parent->story);
		assert(c->length >  0);
		tmp_state = trace;
		// printf("fake: %d , real:%d",tmp_state, c->parent->state);
		// if (tmp_state!=c->parent->state)
		// {
		// 	printf("fail\n");
		// }
		// else{printf("\n");}
		tmp_story = trace_top;
		c = c->parent;
//Fl	printf("TPTS8\n");

	}

//	zFreeSFVec(n->sfv);
	
}
bool zTBTChildSurpass(zTBTreeNode* n, coor_t chop_pos)
{
	printf("child:%d\n",n->children);
	zTBTreeNode *c;
	c = n->child;
	while(c!=NULL)
	{
		printf("c->pos: %d\n",c->pos);
		if (c->pos > chop_pos)
		{
			return true;

		}
		c=c->rsib;
	}
	
	return false;
}
void zManifoldCut(zTrellis *trellis,zTBTreeNode *n,zTBTree *new_t, coor_t chop_pos,int depth)
{
	int tmp_size;
	zTBTreeNode *c;
	printf("%d:%d,%s,%d\nchop: %d, depth:%d\n",n->pos,n->length,zStrIdx2Char(trellis->hmm->state[n->state].name),n->story,chop_pos,depth );
	if (zTBTChildSurpass(n,chop_pos))
	{
		printf("SURPASS\n");
		if (n->pos > new_t->root->pos)
		{
			new_t->root->pos = n->pos;
		}
		if (new_t->root->left_pos==0 || n->pos < new_t->root->left_pos)//smallest for chop
		{
			new_t->root->left_pos = n->pos;
		}

		zTracePartialTrellisSingle(trellis, n, n->state,n->story,n->tree);
		zTBTreeSetChild(new_t->root,n);
	}
	else{
		printf("not surpass\n");
		c = n->child;
		while(c!=NULL)
		{
			zManifoldCut(trellis,c,new_t, chop_pos,depth+1);
			c=n->child;
		}
		printf("mani del %d,%s,%d\n",n->pos,zStrIdx2Char(trellis->hmm->state[n->state].name),n->story );
		if(n==n->tree->root)
		{
			printf("TREEID:%d\n", n->tree->id);
			printf("TREEpos:%d\n", n->tree->root->pos);
			if (n->children!=0)
			{
				printf("aREUOK\n");
			}
			else{
				printf("IMOK\n");
			}
		}
		else
		{zFreeTreeNode(&trellis->cell[n->pos][n->state].layer[n->story]);
		}//not totally deleted
		new_t->size--;
	}
}

void zChunkChop(zTrellis *trellis,coor_t pos,zTBTree *tree)
{
	coor_t old_pos = tree->root->left_pos;
	coor_t diff = pos - tree->root->pos;
	if (diff <= CHUNK+CHUNK+WINDOW)
	{
	//	printf("small diff:%d\n",diff);
		return;
	}
	printf("diff:%d\n",diff);
//	coor_t chop_pos = tree->root->pos+CHUNK; //which cause one point condensed
	coor_t chop_pos = pos-CHUNK-WINDOW;
	zTBTree *new_t;
	new_t=zMalloc(sizeof(zTBTree),"zRunViterbiAndForward: tree");
	zInitTBTree(new_t);

	new_t->size = tree->size+1;//Prevent the extra minus 1 when deleting t->root
	zManifoldCut(trellis,tree->root,new_t,chop_pos,0);
printf("OLD:%d,%d\n",tree->root->left_pos,tree->root->pos );

	tree->root = new_t->root;// here transport the position
	tree->id = new_t->id;
printf("NEW:%d,%d\n",tree->root->left_pos,tree->root->pos );
printf("CHOP: %d,%d\n",old_pos, tree->root->left_pos);
	zFreeViterbiVarsSingle(trellis,old_pos,tree->root->left_pos);


}


coor_t zCheckDomainSingle(zTrellis *trellis,int term_pos,zTBTreeNode *tmp,zTBTree *tree)
{
	coor_t tmp_cpoint=0;
	coor_t new_cpoint=0;
	zTBTreeNode *p;
	zTBTreeNode *c;
	p=tmp;
	if (p==tree->root)
	{
		printf("domain root\n");
	}
	else
	{
		if (p==NULL)
		{
			printf("ERRORRR\n");
			return -1;
		}
		if (p->pos>term_pos)
		{
			printf("OutBound %d\n",p->pos);
			if (p->parent==tree->root)
			{
				return -1;
			}
			return p->parent->pos;
		}

		printf("domain tmp point: %d,%s,%d\n",tmp->pos,zStrIdx2Char(trellis->hmm->state[tmp->state].name),tmp->story);
	}

	c=p->child;
	while(c!=NULL)
	{
		printf("ininin\n");
		new_cpoint=zCheckDomainSingle(trellis,term_pos,c,tree);
		if (new_cpoint==-1)
		{
			return -1;
		}
		if (tmp_cpoint==0)
		{
			tmp_cpoint=new_cpoint;
		}
		else if (tmp_cpoint!=new_cpoint)
		{
			return -1;
		}
		c=c->rsib;
	}    

	return tmp_cpoint;
	//new added by Jim
	if (p->child==NULL)
	{
		printf("SSSSTATE: NULL\n" );
		printf("%d\n", tmp_cpoint);
		return -1;
	}
//	else if(p->child->state<17){printf("TSSSSTATE:%d\n",p->child->state ); printf("%d\n", tmp_cpoint);return -1;}
	else{printf("SSSSTATE:%d\n",p->child->state );printf("%d\n", p->child->pos); printf("%d\n", tmp_cpoint);return tmp_cpoint;}
	

	return tmp_cpoint;
	

	if (c->state<17)
	{
		return tmp_cpoint;
	}
	else{
		return -1;
	}
}


coor_t zCheckDomain(zTrellis *trellis, int fore_pos,int term_pos,zTBTree *tree)
{
	int j,k;
	coor_t tmp_cpoint=0;
	coor_t new_cpoint=0;
	coor_t return_mark=0;
	if (fore_pos==-1)
	{
		return zCheckDomainSingle(trellis,term_pos,tree->root,tree);
	}
	for (j = 0; j<trellis->hmm->states; ++j)
	{
		if (trellis->cell[fore_pos][j].layer==NULL){continue;}
		for (k = 0; k < top_num; ++k)
		{
			if (trellis->cell[fore_pos][j].layer[k]==NULL){continue;}
			if (trellis->cell[fore_pos][j].layer[k]->children>0)
			{
				printf("domain inside\n");
				new_cpoint=zCheckDomainSingle(trellis,term_pos,trellis->cell[fore_pos][j].layer[k],tree);
				printf("New_cpoint:%d\n",new_cpoint );
				if (new_cpoint==-1)
				{
					return -1;
				}
				if (tmp_cpoint==0)
				{
					tmp_cpoint=new_cpoint;
				}
				else if (tmp_cpoint!=new_cpoint)
				{
					return -1;
				}

			}
			
		}
	}
	return tmp_cpoint;
}



coor_t zCheckDomainSingleNew(zTrellis *trellis,int term_pos,zTBTreeNode *tmp,zTBTree *tree)
{
	coor_t tmp_cpoint=0;
	coor_t new_cpoint=0;
	zTBTreeNode *p;
	zTBTreeNode *c;
	p=tmp;
	if (p==tree->root)
	{
		printf("domain root\n");
	}
	else
	{
		if (p==NULL || p->left_pos==0)
		{
			printf("ERRORRR\n");
			return -1;
		}
		if (p->pos>term_pos)
		{
			printf("OutBound %d\n",p->pos);
			if (p->parent==tree->root)
			{
				return -1;
			}
			return p->parent->pos;
		}

		printf("domain tmp point: %d,%s,%d\n",tmp->pos,zStrIdx2Char(trellis->hmm->state[tmp->state].name),tmp->story);
	}

	c=p->child;
	while(c!=NULL || c->left_pos!=0)
	{
		printf("ininin\n");
		new_cpoint=zCheckDomainSingleNew(trellis,term_pos,c,tree);
		if (new_cpoint==-1)
		{
			return -1;
		}
		if (tmp_cpoint==0)
		{
			tmp_cpoint=new_cpoint;
		}
		else if (tmp_cpoint!=new_cpoint)
		{
			return -1;
		}
		c=c->rsib;
	}    

	return tmp_cpoint;
	//new added by Jim
	if (p->child==NULL || p->child->left_pos==0)
	{
		printf("SSSSTATE: NULL\n" );
		printf("%d\n", tmp_cpoint);
		return -1;
	}
//	else if(p->child->state<17){printf("TSSSSTATE:%d\n",p->child->state ); printf("%d\n", tmp_cpoint);return -1;}
	else{printf("SSSSTATE:%d\n",p->child->state );printf("%d\n", p->child->pos); printf("%d\n", tmp_cpoint);return tmp_cpoint;}
	

	return tmp_cpoint;
	
}


coor_t zCheckDomainNew(zTrellis *trellis, int fore_pos,int term_pos,zTBTree *tree)
{
	int j,k;
	coor_t tmp_cpoint=0;
	coor_t new_cpoint=0;
	coor_t return_mark=0;
	if (fore_pos==-1)
	{
		return zCheckDomainSingle(trellis,term_pos,tree->root,tree);
	}
	for (j = 0; j<trellis->hmm->states; ++j)
	{
		if (trellis->cell[fore_pos][j].layer==NULL){continue;}
		for (k = 0; k < top_num; ++k)
		{
		//	if (trellis->cell[fore_pos][j].layer[k]==NULL || trellis->cell[fore_pos][j].layer[k]->left_pos==0){printf("1check\n");continue;}
			if (trellis->cell[fore_pos][j].layer[k]==NULL){printf("1check\n");continue;}

			if (trellis->cell[fore_pos][j].layer[k]->children>0)
			{
				printf("domain inside\n");
				new_cpoint=zCheckDomainSingle(trellis,term_pos,trellis->cell[fore_pos][j].layer[k],tree);
				printf("New_cpoint:%d\n",new_cpoint );
				if (new_cpoint==-1)
				{
					return -1;
				}
				if (tmp_cpoint==0)
				{
					tmp_cpoint=new_cpoint;
				}
				else if (tmp_cpoint!=new_cpoint)
				{
					return -1;
				}

			}
		//	printf("3check  \n");
			
		}
	}
	return tmp_cpoint;
}
coor_t zFindCPointStupid(zTrellis *trellis, int fore_pos,int term_pos,zTBTree *tree)
{
	coor_t step=100;
	coor_t tmp_pos=term_pos;
	coor_t tmp_cpoint=0;
	while(tmp_pos-fore_pos>step)
	{
//		tmp_cpoint=zCheckDomain(trellis,fore_pos,tmp_pos,tree);
		printf("stupid %d,%d\t",fore_pos,tmp_pos);
		tmp_cpoint=zCheckDomainNew(trellis,fore_pos,tmp_pos,tree);
		if(tmp_cpoint!=-1 && tmp_cpoint!=0)
		{
			return tmp_cpoint;
		}	
		tmp_pos-=step;
	}
	return -1;
}

zSfeature zFindPotentialCPoint(zSFVec** tmp_sfv,int mark)
{
	int i=1;
	int j=1;
	while(i<tmp_sfv[0]->size && j<tmp_sfv[1]->size)
	{
		if (tmp_sfv[0]->elem[i].end > tmp_sfv[1]->elem[j].end){i++;}
		else if(tmp_sfv[0]->elem[i].end < tmp_sfv[1]->elem[j].end){j++;}
	//	else if (tmp_sfv[0]->elem[i].state == tmp_sfv[1]->elem[j].state) //original one
	//	else if (tmp_sfv[0]->elem[i].state == tmp_sfv[1]->elem[j].state && tmp_sfv[0]->elem[i].state < 17) // limit cpoint to the end of internal form
		else if (tmp_sfv[0]->elem[i].state == tmp_sfv[1]->elem[j].state && (tmp_sfv[0]->elem[i].state ==14 ||tmp_sfv[0]->elem[i].state ==15) ) // limit cpoint to the end of internal form
		{
			mark--;
			if (mark==0)
			{
				return tmp_sfv[0]->elem[i];
			}
			i++;j++;
		}
		else{i++;j++;}


	}
	zSfeature tmp;
	tmp.state=-1;
	return tmp;

}


void zFilterTmpSFV(zSFVec **tmp_sfv)
{
	int i,j,k,t;
	k=0;t=0;
	int mark=0;
	int tmp_sfv_size=zCountRealSFV(tmp_sfv,top_num);
	printf("inside filter %d\n",tmp_sfv_size);
	

	zSFVec **tmp_tmp_sfv;
	tmp_tmp_sfv = zCalloc(top_num,sizeof(zSFVec*),"zFilterTmpSFV sfv"); // waste storage soon correct
	int size = sizeof(zSFVec);
	for(i=0;i<top_num;i++){
		tmp_tmp_sfv[i]=zMalloc(size,"zMalloc tmp_tmp_sfv[i]");
		zInitSFVec(tmp_tmp_sfv[i],10);
	}

	for (i = 0; i < tmp_sfv_size; ++i)
	{

		if (tmp_sfv[i]==NULL){continue;}
		for (j = i+1; j < tmp_sfv_size; ++j)
		{
			k=0;t=0;
			if (tmp_sfv[j]==NULL){continue;}
	//		if (tmp_sfv[i]->size==tmp_sfv[j]->size)
	//		{
			mark=1;
			while(k < tmp_sfv[i]->size && t < tmp_sfv[j]->size)
			{
				if (tmp_sfv[i]->elem[k].state<17 || (tmp_sfv[i]->elem[k].state>=43 && tmp_sfv[i]->elem[k].state<=46) ){k++;}
				else if (tmp_sfv[j]->elem[t].state<17 || (tmp_sfv[j]->elem[t].state>=43 && tmp_sfv[j]->elem[t].state<=46)){t++;}					

				// if (tmp_sfv[i]->elem[k].state<17){k++;}
				// else if (tmp_sfv[j]->elem[t].state<17){t++;}
				else
				{
					if (tmp_sfv[i]->elem[k].state==tmp_sfv[j]->elem[t].state
						  && tmp_sfv[i]->elem[k].start==tmp_sfv[j]->elem[t].start
						  && tmp_sfv[i]->elem[k].end==tmp_sfv[j]->elem[t].end)
					{
						k++;t++;
					}
					else{mark=0;break;}
				}
			}
	//		}
			if (mark==1)
			{
				if (k==tmp_sfv[i]->size)
				{
					for (; t < tmp_sfv[j]->size; ++t)
					{
						if (!(tmp_sfv[j]->elem[t].state<17 || (tmp_sfv[j]->elem[t].state>=43 && tmp_sfv[j]->elem[t].state<=46)) ){mark=0;break;}
					}
				}
				else if (t==tmp_sfv[j]->size)
				{
					for (; k < tmp_sfv[i]->size; ++k)
					{
						if (!(tmp_sfv[i]->elem[k].state<17 || (tmp_sfv[i]->elem[k].state>=43 && tmp_sfv[i]->elem[k].state<=46) ) ){mark=0;break;}
					}
				}
			}

			if(mark==1)
			{

				// printf("FIL %d,%d,%d: ",i,j,tmp_sfv[i]->size==tmp_sfv[j]->size);
				// for (k = 0; k < tmp_sfv[i]->size; ++k)
				// {
				// 	printf("%d,%d-%d |\n",tmp_sfv[i]->elem[k].state,tmp_sfv[i]->elem[k].start,tmp_sfv[i]->elem[k].end);
				// }
				// printf("\n");
				
				// for (k = 0; k < tmp_sfv[j]->size; ++k)
				// {
				// 	printf("%d,%d-%d |\n",tmp_sfv[j]->elem[k].state,tmp_sfv[j]->elem[k].start,tmp_sfv[j]->elem[k].end);
				// }
				// printf("\n");

				printf("%d is freed\n",j );
				
				zFreeSFVec(tmp_sfv[j]);
				zFree(tmp_sfv[j]);
				tmp_sfv[j]=NULL;
			}
			mark=0;
		}
	}
	k=0;
	for (i = 0; i < tmp_sfv_size; ++i)
	{
		if (tmp_sfv[i]==NULL){continue;}
		zMergeSFVec(tmp_sfv[i],tmp_tmp_sfv[k]);
		k++;
	}
	for (i = 0; i < top_num; ++i)
	{
		if (tmp_sfv[i]!=NULL)
		{
			zFreeSFVec(tmp_sfv[i]);
			zFree(tmp_sfv[i]);
		}
		tmp_sfv[i]=tmp_tmp_sfv[i];
	}
}


coor_t zFindCPointByTmpSFV(zSFVec **tmp_sfv)
{
	int i,j;
	int mark=1;
	int range=2;
	int flag=0;
	zSfeature istate;
	while(mark < tmp_sfv[0]->size)
	{
		istate=zFindPotentialCPoint(tmp_sfv,mark);
	//	printf("ISTATE %d\n",istate.end);
		if (istate.state==-1)
		{
			return -1;
		}
		for (i = 2; i < top_num; ++i)
		{
			flag=0;
			for (j = 1; j <tmp_sfv[i]->size; j++)
			{
				if (istate.end == tmp_sfv[i]->elem[j].end && istate.state == tmp_sfv[i]->elem[j].state)
				{
					flag=1;
					break;
				}
			}
			if (flag!=1)
			{
				break;
			}
			
		}
		if (flag==1)
		{
			return istate.end;
		}
		mark++;
//		printf("PLUSPLUS\n");
	}
	return -1;

}


void zCheckCPoint(zTrellis *trellis,int pos,zTBTree *tree)
{
	// int i,j,k,mark,tmp_state,nullcount;
	// int outbreak=0;
	// for (i = 50; i < pos; ++i)
	// {
	// 	if (trellis->cell[i]==NULL){continue;}
	// 	outbreak=0;
		
	// 	mark=0;
	// 	tmp_state=-1;
	// 	for (j = 0; j < trellis->hmm->states; ++j)
	// 	{
	// 		nullcount=0;
	// 		if (trellis->cell[i][j].layer==NULL){continue;}
	// 		for (k = 0; k < top_num; ++k)
	// 		{
	// 			if (trellis->cell[i][j].layer[k]==NULL){nullcount+=1;continue;}
	// 			if(tmp_state==-1){tmp_state=j;mark=1;}
	// 			else if(tmp_state==j){mark+=1;}
	// 			else{tmp_state=j;mark=0;outbreak=1;break;}


	// 		}
	// 		if (outbreak==1){outbreak=0;break;}
	// 	}
	// 	if(mark!=0)
	// 	{
	// 		printf("CHECK c point: %d,%s,mark:%d,%d\n",i,zStrIdx2Char(trellis->hmm->state[tmp_state].name),mark,nullcount);
	// 	}



	// }

	int tmp_state=-1;
	int mark=1;
	zTBTreeNode *c;
	zTBTreeNode *z;
	zTBTreeNode *tmp_z;
	c=tree->root;

	while(c!=NULL)
	{
		if (c->children==1)
		{
			printf("CHECKTREE c point: %d,%s,%d\n",c->pos,zStrIdx2Char(trellis->hmm->state[c->state].name),c->story);
			c=c->child;
			continue;
		}
		else if (c->children>1)
		{
			tmp_z=z=c->child;
			mark=1;
			while(z!=NULL)
			{
				if (z->state!=tmp_state)
				{
					if (tmp_state==-1)
					{
						tmp_state=z->state;
					}
					else{
						mark=0;
						break;
					}
				}
				z=z->rsib;
			}
			if (mark==0)
			{
				printf("CHECKTREE MORE NOT SINGLE c point: %d,%s,%d    %d\n",c->pos,zStrIdx2Char(trellis->hmm->state[c->state].name),c->story,c->children);
				break;
			}
			else{
				printf("CHECKTREE MORE SINGLE c point: %d,%s,%d\n",c->pos,zStrIdx2Char(trellis->hmm->state[c->state].name),c->story);
				c=tmp_z;
				continue;
			}
			continue;
		}
		else{printf("TREERROR\n");break;}
		printf("TREERROR\n");break;
	}	

}



void zShowTrellis(zTrellis *trellis);
//Added by Jim
void zPlot3DShowTrellis(zTrellis *trellis);
score_t zScoreInternalStateRange(zTrellis *trellis, int state, coor_t start, coor_t end){
	zHMM* hmm = trellis->hmm;
	score_t score = 0;
	score_t tscore;
	coor_t i;

	if(hmm->state[state].type != INTERNAL &&
	   hmm->state[state].type != GINTERNAL){
		zDie("State %s is not an internal state\n",zStrIdx2Char(trellis->hmm->state[state].name));
	}
	
	if(hmm->state[state].type == INTERNAL){
		tscore = zGetTransitionScore(trellis->hmm,state,state,trellis->tiso_group);
	}
	else{
		tscore = zGetIntergenicContinueScore(trellis);
	}
	
	score = zScoreInternalState(trellis,state,start,true);
	
	for(i = start+1; i <= end; i++) {
		score += tscore;
		score += zScoreInternalState(trellis,state,i,false);
	}
	
	return score;
}

//Added by Jim
void zHandShake(zTrellis *trellis,int ext_state,coor_t ext_pos,int story)
{//when found null in tracetrellis, use handshake to recover the lost ones 
	if (trellis->cell[ext_pos][ext_state].layer==NULL)
	{
		zDie("zHandShake NULL error");
	}
	int int_state,int_story,mid_state,mid_story;	
	coor_t int_pos,mid_pos;
	score_t from_score,to_score;
	int i,k;
	zTBTreeNode *c;

	
	mid_pos = ext_pos-trellis->cell[ext_pos][ext_state].layer[story]->length;
	int_state = mid_state = trellis->cell[ext_pos][ext_state].layer[story]->trace;
	mid_story = trellis->cell[ext_pos][ext_state].layer[story]->trace_top;

	trellis->cell[mid_pos][mid_state].layer[mid_story]=zCreateTBTreeNode();

	c = trellis->cell[mid_pos][mid_state].layer[mid_story];

	to_score = zExternalBackInternalScore(trellis,ext_state,ext_pos,story);

	for (i = mid_pos-1; ; i--)
	{
		for (k = 0; k <= mid_story; ++k)
		{
			if(trellis->cell[i][int_state].layer[k]==NULL){continue;}
			from_score = zInternalStartToSpecificScore(trellis,int_state, i, k,mid_pos);
			if (abs(to_score-from_score) < EPISLON)
			{
				c->score = from_score;
				c->length = trellis->cell[i][int_state].layer[k]->length + mid_pos - i;
				c->trace = trellis->cell[i][int_state].layer[k]->trace;
				c->trace_top = trellis->cell[i][int_state].layer[k]->trace_top;
				c->pos = mid_pos;

printf("checkA\n");
		//		zTBTreeSetChild(trellis->cell[c->trace][i-1].layer[c->trace_top]->parent,c);
				zTBTreeSetChild(trellis->cell[i][int_state].layer[k]->parent,c);//newer
				zTBTreeSetChild(c,trellis->cell[ext_pos][ext_state].layer[story]);
				return;
			}
		}
	}
	printf("WTF?\n");
}


void zTracePartialTrellisNew (zTrellis *trellis, coor_t start, coor_t end, int state, zSFVec **sfv,score_t *max_score, zTBTree *tree){
	zHMM      *hmm = trellis->hmm;
	int        max_state, trace;
	int        top_trace_top[top_num]; 
	coor_t     i;
	zSfeature  istate;
	score_t    max_score_temp[top_num];
	int        m,n1,n2;
	int        top_state[top_num];
	int        top_state_temp[top_num];
	int        top_trace_top_temp[top_num];
	int        trace_top;
	bool       flag_tr;

//Added by Jim
//for protection
	coor_t t;
	int f;//for protection
	int tmp;
	int f_mark;
	int minus_count;
	int matched;

	int former_pos;
	int former_state;
	int former_story;


printf("TTTTT1\n");

	if (state < 0) {// input in runv() is -1
		if(!zOption("top")){
			max_score[0] = MIN_SCORE;
			if (trellis->cell[end] == NULL) zDie("zTracePartialTrellisNew end reverse error");
			for (max_state = 0; max_state < hmm->states; max_state++) {
			//	if (trellis->cell[end][max_state] == NULL) continue;
		
			//	max_score_temp[0] = trellis->cell[max_state][end].score[0];

				max_score_temp[0] = trellis->cell[end][max_state].layer[0]->score;
				if (max_score_temp[0] > max_score[0]) {
					state = max_state; max_score[0]= max_score_temp[0];
				}
			}
		}
		else{
			if (trellis->cell[end] == NULL) zDie("zTracePartialTrellisNew end reverse error");
//printf("TTTTT2\n");
			for(m=0;m<top_num;m++){
//				printf("TTTTT2 %d\n",m);
				max_score[m] = MIN_SCORE;
				top_state[m] = -1;
				top_trace_top[m] = -1;
			}
		//	printf("TTTTT2.1\n");

			for (max_state = 0; max_state < hmm->states; max_state++) {
		//	printf("TTTTT2.2\n");
	
				for(m=0;m<top_num;m++){
					max_score_temp[m] = max_score[m]; 
					top_state_temp[m] = top_state[m];
					top_trace_top_temp[m] = top_trace_top[m];
				}
				n1=0;n2=0;
				for(m=0;m<top_num;m++){
		//			if(trellis->cell[max_state][end].score[n1] != MIN_SCORE){

					if(trellis->cell[end][max_state].layer[n1]==NULL || trellis->cell[end][max_state].layer[n1]->score == MIN_SCORE){
						top_state[m] = top_state_temp[n2];
					   	max_score[m] = max_score_temp[n2];
					   	top_trace_top[m] = top_trace_top_temp[n2];
					   	n2++;

					}
					else{
			//		if(trellis->cell[end][max_state].layer[n1]->score != MIN_SCORE){
		//				printf("TTT2.222\n");
						flag_tr = (max_score_temp[n2] == MIN_SCORE);
						if(!flag_tr)
							flag_tr = (max_score_temp[n2] < trellis->cell[end][max_state].layer[n1]->score);
			//						   trellis->cell[max_state][end].score[n1]);
						if (flag_tr) {
							top_state[m] = max_state; 
					//		max_score[m] = trellis->cell[max_state][end].score[n1];
							max_score[m] = trellis->cell[end][max_state].layer[n1]->score;
							top_trace_top[m] = n1;
							n1++;
						}
						else{
							top_state[m] = top_state_temp[n2];
							max_score[m] = max_score_temp[n2];
							top_trace_top[m] = top_trace_top_temp[n2];
							n2++;
						}
					}
					// else{
					//    	top_state[m] = top_state_temp[n2];
					//    	max_score[m] = max_score_temp[n2];
					//    	top_trace_top[m] = top_trace_top_temp[n2];
					//    	n2++;
					// }
				}
			}

//printf("TTTTT3\n");
		}
	}
//	zShowTrellis(trellis);
//	zDie("show trellis now!\n");
//	printf("TTTTT4\n");
	if(!zOption("top")){
		i = end;

	//protection
//		if(trellis->cell[state][i].layer[0]==NULL){return;}

	//	if (-1 == trellis->cell[state][i].trace[0]) {
		if (-1 == trellis->cell[state][i].layer[0]->trace) {
			return; /* the sequence can't end in this state */
		}
	
		while (i > start) {
		//protection  easy one maybe changed soon
			if(trellis->cell[state][i].layer[0]==NULL)
			{
					if ((trellis->hmm->state[state].type==INTERNAL)
						|| (trellis->hmm->state[state].type==GINTERNAL)) 
					{	
						trellis->cell[state][i].layer[0]=zCreateTBTreeNode();
						t=i-1;
						while(trellis->cell[state][t].layer[0]==NULL)
						{	
							t--;
						}
						trellis->cell[state][i].layer[0]->score=trellis->cell[state][t].layer[0]->score;//should be changed soon
						trellis->cell[state][i].layer[0]->length=trellis->cell[state][t].layer[0]->length+t-i;
						trellis->cell[state][i].layer[0]->trace=trellis->cell[state][t].layer[0]->trace;
						trellis->cell[state][i].layer[0]->trace_top=trellis->cell[state][t].layer[0]->trace_top;


					}
					else
						{zDie("zTracePartialTrellis protection external error.");}
				
			}


		//	trace = trellis->cell[state][i].trace[0];
			trace = trellis->cell[state][i].layer[0]->trace;
			if (-1 == trace) {
				zDie("Impossible to perform traceback from this state.");
			}

			istate.name   = hmm->state[state].name;
		//	istate.start  = i - trellis->cell[state][i].length[0] + 1;
			istate.start  = i - trellis->cell[state][i].layer[0]->length + 1;
			istate.end    = i;
			istate.strand = hmm->state[state].strand;

		//	zChar2FragFrame(trellis->cell[state][i].frag_frame[0], &istate);
			zChar2FragFrame(trellis->cell[state][i].layer[0]->frag_frame, &istate);

			istate.state  = state;
			istate.group  = NULL;

			// istate.score = trellis->cell[state][istate.end].score[0]
			// 	- trellis->cell[trace][istate.start-1].score[0]

			istate.score = trellis->cell[state][istate.end].layer[0]->score
				- trellis->cell[trace][istate.start-1].layer[0]->score				
				- zGetTransitionScore(trellis->hmm,trace,state,trellis->tiso_group);
			/*EVAN now removing transition scores from output */

			zPushSFVec(sfv[0], &istate);

			// assert(trellis->cell[state][i].length[0] >  0);
			// i -= trellis->cell[state][i].length[0];

			assert(trellis->cell[state][i].layer[0]->length >  0);
			i -= trellis->cell[state][i].layer[0]->length;
			state=trace;
		}
	}
	else{
		int ttt;
//printf("TTTTT5\n");
		for(m=0;m<top_num;m++){
			if(max_score[m] != MIN_SCORE){
				if(top_trace_top[m]<0 || top_trace_top[m] >= top_num){
					zDie("top_trace_top out of range!\n");
				}
			}
		}
//printf("TTTTT6\n");		
		//Added by Jim

		


		zTBTreeNode *c;
		for(m=0;m<top_num;m++){
			if(max_score[m] == MIN_SCORE) break;

			i = end;

			state = top_state[m];
			ttt = top_trace_top[m];
			c=trellis->cell[i][state].layer[ttt];

		//	printf("SPPP\n");
			zTracePartialTrellisSingle(trellis, c, state, ttt, tree);
		//				printf("SPPP2\n");

			zMergeSFVec(c->sfv,sfv[m]);

//printf("hahaa\n");
//printf("ccSFV:%d,%f,%d,%d,%d\n",c->pos,c->score,c->length,c->state,c->story);
			zFreeSFVec(c->sfv);//newly added
			zFree(c->sfv);
			c->sfv=NULL;
		}
	}
//	zShowTrellis(trellis);//Jim
}

void zTraceTrellisNew (zTrellis *trellis, int state, zSFVec **sfv, score_t *score,zTBTree *tree) {
    zTracePartialTrellisNew(trellis, PADDING, trellis->dna->length - 1 - PADDING, state, sfv, score, tree);
}
void zTracePartialTrellis (zTrellis *trellis, coor_t start, coor_t end, int state, zSFVec **sfv,score_t *max_score){
	zHMM      *hmm = trellis->hmm;
	int        max_state, trace;
	int        top_trace_top[top_num]; 
	coor_t     i;
	zSfeature  istate;
	score_t    max_score_temp[top_num];
	int        m,n1,n2;
	int        top_state[top_num];
	int        top_state_temp[top_num];
	int        top_trace_top_temp[top_num];
	int        trace_top;
	bool       flag_tr;

//Added by Jim
//for protection
	coor_t t;
	int f;//for protection
	int tmp;
	int f_mark;
	int minus_count;
	int matched;

	int former_pos;
	int former_state;
	int former_story;

printf("TTTTT1\n");

	if (state < 0) {// input in runv() is -1
		if(!zOption("top")){
			max_score[0] = MIN_SCORE;
			for (max_state = 0; max_state < hmm->states; max_state++) {
				if (trellis->cell[max_state] == NULL) continue;
		
			//	max_score_temp[0] = trellis->cell[max_state][end].score[0];

				max_score_temp[0] = trellis->cell[max_state][end].layer[0]->score;
				if (max_score_temp[0] > max_score[0]) {
					state = max_state; max_score[0]= max_score_temp[0];
				}
			}
		}
		else{

//printf("TTTTT2\n");
			for(m=0;m<top_num;m++){
				max_score[m] = MIN_SCORE;
				top_state[m] = -1;
				top_trace_top[m] = -1;
			}
			for (max_state = 0; max_state < hmm->states; max_state++) {
				if (trellis->cell[max_state] == NULL) continue;
				for(m=0;m<top_num;m++){
					max_score_temp[m] = max_score[m]; 
					top_state_temp[m] = top_state[m];
					top_trace_top_temp[m] = top_trace_top[m];
				}
				n1=0;n2=0;
				for(m=0;m<top_num;m++){
		//			if(trellis->cell[max_state][end].score[n1] != MIN_SCORE){
					if(trellis->cell[max_state][end].layer[n1]->score != MIN_SCORE){
						flag_tr = (max_score_temp[n2] == MIN_SCORE);
						if(!flag_tr)
							flag_tr = (max_score_temp[n2] < trellis->cell[max_state][end].layer[n1]->score);
			//						   trellis->cell[max_state][end].score[n1]);
						if (flag_tr) {
							top_state[m] = max_state; 
					//		max_score[m] = trellis->cell[max_state][end].score[n1];
							max_score[m] = trellis->cell[max_state][end].layer[n1]->score;
							top_trace_top[m] = n1;
							n1++;
						}
						else{
							top_state[m] = top_state_temp[n2];
							max_score[m] = max_score_temp[n2];
							top_trace_top[m] = top_trace_top_temp[n2];
							n2++;
						}
					}
					else{
					   	top_state[m] = top_state_temp[n2];
					   	max_score[m] = max_score_temp[n2];
					   	top_trace_top[m] = top_trace_top_temp[n2];
					   	n2++;
					}
				}
			}

//printf("TTTTT3\n");
		}
	}
//	zShowTrellis(trellis);
//	zDie("show trellis now!\n");
//	printf("TTTTT4\n");
	if(!zOption("top")){
		i = end;

	//protection
//		if(trellis->cell[state][i].layer[0]==NULL){return;}

	//	if (-1 == trellis->cell[state][i].trace[0]) {
		if (-1 == trellis->cell[state][i].layer[0]->trace) {
			return; /* the sequence can't end in this state */
		}
	
		while (i > start) {
		//protection  easy one maybe changed soon
			if(trellis->cell[state][i].layer[0]==NULL)
			{
					if ((trellis->hmm->state[state].type==INTERNAL)
						|| (trellis->hmm->state[state].type==GINTERNAL)) 
					{	
						trellis->cell[state][i].layer[0]=zCreateTBTreeNode();
						t=i-1;
						while(trellis->cell[state][t].layer[0]==NULL)
						{	
							t--;
						}
						trellis->cell[state][i].layer[0]->score=trellis->cell[state][t].layer[0]->score;//should be changed soon
						trellis->cell[state][i].layer[0]->length=trellis->cell[state][t].layer[0]->length+t-i;
						trellis->cell[state][i].layer[0]->trace=trellis->cell[state][t].layer[0]->trace;
						trellis->cell[state][i].layer[0]->trace_top=trellis->cell[state][t].layer[0]->trace_top;


					}
					else
						{zDie("zTracePartialTrellis protection external error.");}
				
			}


		//	trace = trellis->cell[state][i].trace[0];
			trace = trellis->cell[state][i].layer[0]->trace;
			if (-1 == trace) {
				zDie("Impossible to perform traceback from this state.");
			}

			istate.name   = hmm->state[state].name;
		//	istate.start  = i - trellis->cell[state][i].length[0] + 1;
			istate.start  = i - trellis->cell[state][i].layer[0]->length + 1;
			istate.end    = i;
			istate.strand = hmm->state[state].strand;

		//	zChar2FragFrame(trellis->cell[state][i].frag_frame[0], &istate);
			zChar2FragFrame(trellis->cell[state][i].layer[0]->frag_frame, &istate);

			istate.state  = state;
			istate.group  = NULL;

			// istate.score = trellis->cell[state][istate.end].score[0]
			// 	- trellis->cell[trace][istate.start-1].score[0]

			istate.score = trellis->cell[state][istate.end].layer[0]->score
				- trellis->cell[trace][istate.start-1].layer[0]->score				
				- zGetTransitionScore(trellis->hmm,trace,state,trellis->tiso_group);
			/*EVAN now removing transition scores from output */

			zPushSFVec(sfv[0], &istate);

			// assert(trellis->cell[state][i].length[0] >  0);
			// i -= trellis->cell[state][i].length[0];

			assert(trellis->cell[state][i].layer[0]->length >  0);
			i -= trellis->cell[state][i].layer[0]->length;
			state=trace;
		}
	}
	else{
		int ttt;
//printf("TTTTT5\n");
		for(m=0;m<top_num;m++){
			if(max_score[m] != MIN_SCORE){
				if(top_trace_top[m]<0 || top_trace_top[m] >= top_num){
					zDie("top_trace_top out of range!\n");
				}
			}
		}
printf("TTTTT6\n");		
		for(m=0;m<top_num;m++){
			if(max_score[m] == MIN_SCORE) break;
			i = end;
			
			/*	if (-1 == trellis->cell[state][i].trace[trace_top]) {
			return; 
			}*/ 
			state = top_state[m];
			ttt = top_trace_top[m];
printf("TTTTT7\n");
			while (i > start) {
	//			printf("TTTTT8: %s,%d,%d\n",zStrIdx2Char(trellis->hmm->state[state].name),i,ttt);
				

				
			 //   	trace = trellis->cell[state][i].trace[ttt];
				// trace_top = trellis->cell[state][i].trace_top[ttt];
			   	trace = trellis->cell[state][i].layer[ttt]->trace;
				trace_top = trellis->cell[state][i].layer[ttt]->trace_top;
				if(trace_top == -1) zDie("trace_top error!\n");

				istate.name   = hmm->state[state].name;
			//	istate.start  = i - trellis->cell[state][i].length[ttt] + 1;
				istate.start  = i - trellis->cell[state][i].layer[ttt]->length + 1;
				istate.end    = i;
				istate.strand = hmm->state[state].strand;

			//	zChar2FragFrame(trellis->cell[state][i].frag_frame[ttt], &istate);
				zChar2FragFrame(trellis->cell[state][i].layer[ttt]->frag_frame, &istate);
			
				istate.state  = state;
				istate.group  = NULL;
				
				//protection revealed soonnnn!
				if (trellis->cell[trace][istate.start-1].layer[trace_top]==NULL)
				{
				//	printf("TTTTT9\n");
					printf("HANDSHAKE: %s,%d,%d\n",zStrIdx2Char(trellis->hmm->state[trace].name),istate.start-1,trace_top);
					zHandShake(trellis,state,i,ttt);
				}
		

			  //  	istate.score = trellis->cell[state][istate.end].score[ttt]
					// - trellis->cell[trace][istate.start-1].score[trace_top]
				istate.score = trellis->cell[state][istate.end].layer[ttt]->score
					- trellis->cell[trace][istate.start-1].layer[trace_top]->score
					- zGetTransitionScore(trellis->hmm,trace,state,trellis->tiso_group);
				/*EVAN now removing transition scores from output */
				zPushSFVec(sfv[m], &istate);

				//Added by Jim
				// former_pos = i;
				// former_state = state;
				// former_story = ttt;


				// assert(trellis->cell[state][i].length[ttt] >  0);
				// i -= trellis->cell[state][i].length[ttt];			
				assert(trellis->cell[state][i].layer[ttt]->length >  0);
				i -= trellis->cell[state][i].layer[ttt]->length;
				state = trace;
			   	ttt = trace_top;
			}
		}
	}
//	zShowTrellis(trellis);//Jim
}

void zTraceTrellis (zTrellis *trellis, int state, zSFVec **sfv, score_t *score) {
    zTracePartialTrellis(trellis, PADDING, trellis->dna->length - 1 - PADDING, state, sfv, score);
}
//Added by Jim
//void zInitTrellisCellLayer(zTrellisCellLayer **layer)
void zInitTreeNode(zTBTreeNode** layer)
{
	int k;
//	layer=zMalloc(top_num*sizeof(zTrellisCellLayer*), "zAllocViterbi layer");
//	printf("pafflo\n");

//	layer = zMalloc(10000*sizeof(zTBTreeNode), "zAllocViterbi layer");


	for(k=0; k < top_num; k++){
		layer[k]=zCreateTBTreeNode();
	//	layer[k]=zMalloc(sizeof(zTBTreeNode), "zInitTrellisCellLayer layer");
//		layer[k]=malloc(sizeof(zTrellisCellLayer));
		// printf("COOLLLL\n");
		// free(layer[k]);
		// layer[k]=NULL;
		// if (layer[k]==NULL)
		// {
		// 	printf("MORRRR\n");
		// }
		layer[k]->score = MIN_SCORE;
	}
//	printf("PPP\n");
//	printf("%f\n",layer[0]->score );
}

static void zAllocViterbiVars(zTrellis *trellis) {
	coor_t    j;
	int       i;
	int       k;
	size_t size;
	if (NULL != trellis->cell) return;
//	trellis->cell = zCalloc(trellis->hmm->states, sizeof(zTrellisCell*), "zAllocViterbi cells");
	trellis->cell = zCalloc(trellis->dna->length, sizeof(zTrellisCell*), 
							"zAllocViterbi cells");
	/* Allocating backward link arrays for all states                 */
	/* But only those referring to external states will be used later */
	trellis->extpos = zMalloc((trellis->hmm->states*sizeof(zIVec)), "zAllocViterbi extpos");

//	size = trellis->dna->length * sizeof(zTrellisCell);
	size = trellis->hmm->states * sizeof(zTrellisCell);
	printf("SSSize %d\n",size);
	for (i = 0; i < trellis->dna->length; i++) {
		trellis->cell[i] = zMalloc(size, "zAllocViterbi cells[i]");
	//	zInitIVec(&trellis->extpos[i], 1);
		//Added by Jim
		for (j = 0; j < trellis->hmm->states; j++){
			if (i==0)
			{
				zInitIVec(&trellis->extpos[j], 1);
			}
		//	trellis->cell[i][j].layer=zMalloc(top_num*sizeof(zTrellisCellLayer*), "zAllocViterbi layer");
		//	zInitTrellisCellLayer(trellis->cell[i][j].layer);
			zInitTreeNode(trellis->cell[i][j].layer);
		}
	}
}
//Added by Jim
//void zAllocViterbiVarsSingle(zTrellis *trellis, coor_t j, int i) 
void zAllocViterbiVarsInitial(zTrellis *trellis) {
	coor_t    j;
	int       i;
	int       k;
	size_t size;

	if (NULL != trellis->cell){return;}

//	trellis->cell = zCalloc(trellis->hmm->states, sizeof(zTrellisCell*), "zAllocViterbi cells");
	trellis->cell = zCalloc(trellis->dna->length, sizeof(zTrellisCell*), "zAllocViterbi cells");


	/* Allocating backward link arrays for all states                 */
	/* But only those referring to external states will be used later */
	trellis->extpos = zMalloc((trellis->hmm->states*sizeof(zIVec)), "zAllocViterbi extpos");

	for (j = 0; j < trellis->hmm->states; j++){
		zInitIVec(&trellis->extpos[j], 1);
		//	trellis->cell[i][j].layer=zMalloc(top_num*sizeof(zTrellisCellLayer*), "zAllocViterbi layer");
		//	zInitTrellisCellLayer(trellis->cell[i][j].layer);
	}

}

void zResetViterbiVarsSingle(zTrellis *trellis, coor_t i,int j) {
//	int j;
	int       k;
	size_t size;
//	if (NULL != trellis->cell[i]) return;

	if (NULL == trellis->cell[i] || trellis->cell[i][j].layer==NULL)
	{
		zDie("zResetViterbiVarsSingle cell error");
	}
	printf("inside RESET\n");
	for (k = 0; k < top_num; ++k)
	{
		printf("RESET %d,%d,%d\n",i,j,k );
		zResetTBTreeNode(trellis->cell[i][j].layer[k]);
	}
}


void zAllocViterbiVarsSingle(zTrellis *trellis, coor_t i,int j) {
//	int j;
	int       k;
	size_t size;
//	if (NULL != trellis->cell[i]) return;

	if (NULL == trellis->cell[i])
	{
		size = trellis->hmm->states * sizeof(zTrellisCell);
	//	printf("alloc cell %d\n",i );
		trellis->cell[i] = zMalloc(size, "zAllocViterbiVarsSingle cells[i]");
	//	trellis->cell[i] = zCalloc(trellis->hmm->states,sizeof(zTrellisCell), "zAllocViterbi cells[i]");

	}
//	printf("alloc layer %d,%d\n",i,j );
	trellis->cell[i][j].layer=zCalloc(top_num,sizeof(zTBTreeNode*),"zAllocViterbiVarsSingle layer");
	zInitTreeNode(trellis->cell[i][j].layer);

// 	size = trellis->hmm->states * sizeof(zTrellisCell);
// //	printf("SSSize %d\n",size);
// 	trellis->cell[i] = zMalloc(size, "zAllocViterbi cells[i]");

// 	//Added by Jim
// //	for (j = 0; j < trellis->hmm->states; j++){
// 		zInitTreeNode(trellis->cell[i][j].layer);
// //	}

}
void zAllocViterbiVarsSingleForCPoint(zTrellis *trellis, coor_t i,int j) {//useless now
	int       k;
	size_t size;
//	if (NULL != trellis->cell[i]) return;
printf("FC:%d,%d\n",i,j );
	if (NULL == trellis->cell[i])
	{printf("in FC:%d,%d\n",i,j );
		size = trellis->hmm->states * sizeof(zTrellisCell);
		trellis->cell[i] = zMalloc(size, "zAllocViterbi cells[i]");
	}
printf("FC2:%d,%d\n",i,j );	
	trellis->cell[i][j].layer=zCalloc(top_num,sizeof(zTBTreeNode*),"zAllocViterbi layer");
	zInitTreeNode(trellis->cell[i][j].layer);
	
printf("FC3:%d,%d\n",i,j );	

	
}
	
void zFreeViterbiVars(zTrellis *trellis) {
	int    i;
	int j;
	int	k;

	if (NULL == trellis->cell) return;
	for (i = 0; i < trellis->dna->length; i++) {
//	for (i = 0; i < trellis->hmm->states; i++) {
		if (NULL == trellis->cell[i]) continue;
 		//Added by Jim
	// 	for (j = 0; j < trellis->dna->length; j++)
	// 	{
	// 		for (k = 0; k < top_num; k++)
	// 		{
	// //			if(trellis->cell[i][j].layer[k]!=NULL)
	// //			{zFree(trellis->cell[i][j].layer[k]);}
				
	// 		}
	// 		zFree(trellis->cell[i][j].layer);
	// 	}

		zFree(trellis->cell[i]);
		trellis->cell[i] = NULL;

		
	}
	for (i = 0; i < trellis->hmm->states; i++) 
	{
		zFreeIVec(&trellis->extpos[i]);
	}
	

	zFree(trellis->cell);
	trellis->cell = NULL;

	zFree(trellis->extpos);
	trellis->extpos = NULL;
}





void zAllocForwardVars(zTrellis *trellis) {
	int i;
	if (NULL != trellis->forward) return;
	trellis->forward = zCalloc(trellis->hmm->states, sizeof(score_t*), 
							   "zAllocForward forward");
	for (i = 0; i < trellis->hmm->states; i++) {
		trellis->forward[i] = zMalloc(trellis->dna->length * sizeof(score_t),
									  "zAllocForward forward[i]");
	}
}

void zFreeForwardVars(zTrellis *trellis) {
	int i;
	if (NULL == trellis->forward) return;
	for (i = 0; i < trellis->hmm->states; i++) {
		if (NULL == trellis->forward[i]) continue;
		zFree(trellis->forward[i]);
		trellis->forward[i] = NULL;
	}
	zFree(trellis->forward);
	trellis->forward = NULL;
}

static void zAllocBackwardVars(zTrellis *trellis) {
	int i;
	if (NULL != trellis->forward) return;
	trellis->backward = zCalloc(trellis->hmm->states, sizeof(score_t*),
								"zAllocBackwardVars backward");
	for (i = 0; i < trellis->hmm->states; i++) {
		trellis->backward[i] = zMalloc(trellis->dna->length * sizeof(score_t),
									   "zAllocBackward forward[i]");
	}
}

void zFreeBackwardVars(zTrellis *trellis) {
	int i;
	if (NULL == trellis->backward) return;
	for (i = 0; i < trellis->hmm->states; i++) {
		if (NULL == trellis->backward[i]) continue;
		zFree(trellis->backward[i]);
		trellis->backward[i] = NULL;
	}
	zFree(trellis->forward);
	trellis->forward = NULL;
}

static void zAllocFactories(zTrellis *trellis, bool conseq_enabled) {
	zFeatureFactory **ffactory;
	zHMM_State        *state;
	zIVec           **fmap5,**fmap3,*jumps;

	zHMM    *hmm = trellis->hmm;
	int      i,j,k,insert;
	int      ns;

	/* create one zStopSeq for each strand */
	trellis->stopseq = zMalloc(sizeof(zStopSeq),"zAllocFactories stopseq");
	trellis->rstopseq = zMalloc(sizeof(zStopSeq),"zAllocFactories rstopseq");
	zInitStopSeq(trellis->stopseq,trellis->unmasked_dna == NULL ? trellis->dna :
				 trellis->unmasked_dna);
	zInitStopSeq(trellis->rstopseq,trellis->unmasked_rdna == NULL ? trellis->rdna :
				 trellis->unmasked_rdna);
	
	
	/* create factories for non-internal states */
	trellis->factory = zCalloc(hmm->feature_count, sizeof(zFeatureFactory*),
							   "zAllocFactories: factory[]");
	trellis->rfactory = zCalloc(hmm->feature_count, sizeof(zFeatureFactory*),
								"zAllocFactories: ractory[]");
	trellis->fmap5   = zCalloc(hmm->feature_count, sizeof(zIVec*),
														 "zAllocFactories: fmap5[]");
	trellis->rfmap5  = zCalloc(hmm->feature_count, sizeof(zIVec*),
														 "zAllocFactories: rfmap5[]");
	trellis->fmap3   = zCalloc(hmm->feature_count, sizeof(zIVec*),
														 "zAllocFactories: fmap3[]");
	trellis->rfmap3  = zCalloc(hmm->feature_count, sizeof(zIVec*),
														 "zAllocFactories: rfmap3[]");
	for (i = 0; i <  hmm->states; i++) {
		state = &hmm->state[i];
		if (NULL == state->ffactory) continue;

		ns = ('-' == state->strand);
		ffactory = (ns) ? &trellis->rfactory[state->model] : 
			&trellis->factory[state->model];
		fmap5 = (ns) ? &trellis->rfmap5[state->model] : &trellis->fmap5[state->model];
		fmap3 = (ns) ? &trellis->rfmap3[state->model] : &trellis->fmap3[state->model];
		if (NULL != *ffactory){
			zPushIVec(*fmap3,i);
			jumps = trellis->hmm->jmap[i];
			for(j = 0;j < jumps->size; j++){
				insert = 1;
				for(k = 0;k < (*fmap5)->size; k++){
					if((*fmap5)->elem[k] == jumps->elem[j]){
						insert = 0;
						break;
					}
				}
				if(insert == 1){
					zPushIVec(*fmap5,jumps->elem[j]);
				}
			}
			continue;
		}
		*fmap5 = zMalloc(sizeof(zIVec),"zAllocFacotries: fmap5");
		*fmap3 = zMalloc(sizeof(zIVec),"zAllocFacotries: fmap3");
		zInitIVec(*fmap5,3);
		zInitIVec(*fmap3,3);
		zPushIVec(*fmap3,i);
		jumps = trellis->hmm->jmap[i];
		for(j = 0;j < jumps->size; j++){
			zPushIVec(*fmap5,jumps->elem[j]);
		}
		*ffactory = zMalloc(sizeof(zFeatureFactory),"zAllocFac: factory[i]");
		state->ffactory(*ffactory,
						(ns) ? trellis->rdna     : trellis->dna, 
						(ns) ? trellis->rscanner : trellis->scanner, 
						(ns) ? trellis->rconseqscanner: trellis->conseqscanner,
							hmm->feature_count, 
						state->strand, 
						conseq_enabled,
							(ns) ? trellis->rstopseq : trellis->stopseq
						);
	
	}
	/* allcate external arrays */
	trellis->external = zMalloc(hmm->feature_count*sizeof(zSFList*),
								"zAllocFactories: external[]");
	trellis->rexternal = zMalloc(hmm->feature_count*sizeof(zSFList*),
								 "zAllocFactories: rexternal[]");
	for (i = 0; i <  hmm->feature_count; i++) {
		if(trellis->factory[i] != NULL){
			trellis->external[i] = zMalloc(sizeof(zSFList),"zAllocFactories: external[i]");
			zInitSFList(trellis->external[i]);
		}
		if(trellis->rfactory[i] != NULL){
			trellis->rexternal[i] = zMalloc(sizeof(zSFList),"zAllocFactories: rexternal[i]");
			zInitSFList(trellis->rexternal[i]);
		}
	}
}

static void zFreeFactories(zTrellis *trellis) {
	int i;

	for (i = 0; i < trellis->hmm->feature_count; i++) {
		if (NULL != trellis->factory[i]) {
			zFreeFeatureFactory(trellis->factory[i]);
			zFree(trellis->factory[i]);
			zFreeSFList(trellis->external[i]);
			zFree(trellis->external[i]);
			zFreeIVec(trellis->fmap5[i]);
			zFreeIVec(trellis->fmap3[i]);
			zFree(trellis->fmap5[i]);
			zFree(trellis->fmap3[i]);
		}
		if (NULL != trellis->rfactory[i]) {
			zFreeFeatureFactory(trellis->rfactory[i]);
			zFree(trellis->rfactory[i]);
			zFreeSFList(trellis->rexternal[i]);
			zFree(trellis->rexternal[i]);
			zFreeIVec(trellis->rfmap5[i]);
			zFreeIVec(trellis->rfmap3[i]);
			zFree(trellis->rfmap5[i]);
			zFree(trellis->rfmap3[i]);
		}
	}
	zFree(trellis->fmap5);
	zFree(trellis->rfmap5);
	zFree(trellis->fmap3);
	zFree(trellis->rfmap3);
	zFree(trellis->factory);
	zFree(trellis->rfactory);
	zFree(trellis->external);
	zFree(trellis->rexternal);
	trellis->external = trellis->rexternal = NULL;
	trellis->factory = trellis->rfactory = NULL;

	zFreeStopSeq(trellis->stopseq);
	zFreeStopSeq(trellis->rstopseq);
	zFree(trellis->stopseq);
	zFree(trellis->rstopseq);
	trellis->stopseq = NULL;
	trellis->rstopseq = NULL;
}

void zInitTrellis (zTrellis *trellis, char* dna_file, zHMM *hmm, char* conseq_file, 
				   char* unmasked_dna_file, char* snp_file, char* p_top_num) {
	int         i;
	zDNA       *dna;
	bool        conseq_enabled;

	/* clear out pointers */
	trellis->dna       = NULL;
	trellis->rdna      = NULL;
	trellis->conseq    = NULL;
	trellis->rconseq    = NULL;
	trellis->hmm       = NULL;
	trellis->fexternal = NULL;
	trellis->cell      = NULL;
	trellis->forward   = NULL;
	trellis->backward  = NULL;
	trellis->scanner   = NULL;
	trellis->factory   = NULL;
	trellis->extpos    = NULL;
	trellis->external  = NULL;
	trellis->rexternal = NULL;

	if(conseq_file != NULL) {
		conseq_enabled = true;
	}
	else {
		conseq_enabled = false;
	}
	
	/* initial setup */


	trellis->dna  = zMalloc(sizeof(zDNA), "zInitTrellis dna");
	trellis->rdna = zMalloc(sizeof(zDNA), "zInitTrellis rdna");
	zInitDNA(trellis->dna);
	zInitDNA(trellis->rdna);
	dna = trellis->dna;
	zLoadDNAFromFasta(dna,dna_file,snp_file);
	zSetDNAPadding(dna,PADDING);
	trellis->padding = PADDING;
	zCopyDNA(trellis->dna, trellis->rdna);
	zAntiDNA(trellis->rdna);
	if(unmasked_dna_file != NULL){
		trellis->unmasked_dna  = zMalloc(sizeof(zDNA), "zInitTrellis unmasked_dna");
		trellis->unmasked_rdna = zMalloc(sizeof(zDNA), "zInitTrellis unmasked_rdna");
		zInitDNA(trellis->unmasked_dna);
		zInitDNA(trellis->unmasked_rdna);
		dna = trellis->dna;
		zLoadDNAFromFasta(trellis->unmasked_dna,unmasked_dna_file,snp_file);
		zSetDNAPadding(trellis->unmasked_dna,PADDING);
		zCopyDNA(trellis->unmasked_dna, trellis->unmasked_rdna);
		zAntiDNA(trellis->unmasked_rdna);
	}
	else{
		trellis->unmasked_dna = NULL;
		trellis->unmasked_rdna = NULL;
	}
	if(p_top_num!=NULL){
	}
	if(conseq_enabled) {
		trellis->conseq = zMalloc(sizeof(zConseq), "zInitTrellis conseq");
		trellis->rconseq = zMalloc(sizeof(zConseq), "zInitTrellis rconseq");
		zInitConseq(trellis->conseq,3);
		zInitConseq(trellis->rconseq,3);
		zLoadConseqFromFasta(trellis->conseq,conseq_file);
		zSetConseqPadding(trellis->conseq,PADDING);
		zCopyConseq(trellis->conseq, trellis->rconseq);
		zReverseConseq(trellis->rconseq);
	}
	
	trellis->hmm  = hmm;

	/* isocore group */
	trellis->tiso_group = -1;
	if(trellis->hmm->iso_transitions > 0) {
		for(i=0;i<trellis->hmm->iso_transitions;i++){
			if(zGetDNAGC(trellis->dna) < trellis->hmm->iso_transition[i]){
				trellis->tiso_group = i;
				break;
			}
		}
	}
	trellis->iiso_group = -1;
	if(hmm->iso_states > 0){
		for(i = 0; i < hmm->iso_states; i++){
			if(zGetDNAGC(trellis->dna) < hmm->iso_state[i]){
				trellis->iiso_group=i;
				break;
			}
		}
	}	

	/* create scanners */
	trellis->scanner  = zCalloc(hmm->feature_count, sizeof(zScanner*), 
								"zInitTrellis: scanner[]");
	trellis->rscanner = zCalloc(hmm->feature_count, sizeof(zScanner*), 
								"zInitTrellis: rscanner[]");
	
	for (i = 0; i < hmm->feature_count; i++) {
		if (hmm->mmap[i] == NULL) continue;
		if (hmm->mmap[i]->seq_type == DNA) {
			trellis->scanner[i] = zMalloc(sizeof(zScanner), "zInitTrellis scan");
			zInitScanner(trellis->scanner[i], trellis->dna->seq, hmm->mmap[i]);
			trellis->rscanner[i] = zMalloc(sizeof(zScanner), "zInitTrellis rscan");
			zInitScanner(trellis->rscanner[i], trellis->rdna->seq, hmm->mmap[i]);

			/* Pre-compute DNA scanners for EXPLICIT    */
			/* states as a "running sum". Required by   */
			/* the implementation of EXPLICIT states    */
			/* (see zScoreScanners in zTransition.c for */
			/* details).                                */
			
			if (zUsedInExplicitState(hmm, i)){
				zPreComputeScanner(trellis->scanner[i]);
				zPreComputeScanner(trellis->rscanner[i]);
			}
		}
		else {
			zDie("non-DNA model in sequence models\n");
		}
	}

	/* create conseq scanners */
	if(conseq_enabled) {
		trellis->conseqscanner = 
			zCalloc(hmm->feature_count, sizeof(zScanner*), 
					"zInitTrellis: conseqscanner[]");
		trellis->rconseqscanner = 
			zCalloc(hmm->feature_count, sizeof(zScanner*), 
					"zInitTrellis: rconseqscanner[]");
	} else {
		trellis->conseqscanner = trellis->rconseqscanner = NULL;
	}
	for (i = 0; i < hmm->feature_count; i++) {
		if (NULL == hmm->cmmap[i]) continue;
		if (hmm->cmmap[i]->seq_type == CONSEQ) {
			if (conseq_enabled) {
				trellis->conseqscanner[i] = zMalloc(sizeof(zScanner), 
													"zInitTrellis conseqscan");
				zInitScanner(trellis->conseqscanner[i], trellis->conseq->seq, 
							 hmm->cmmap[i]);
				trellis->rconseqscanner[i] = zMalloc(sizeof(zScanner), 
													 "zInitTrellis rcscan");
				zInitScanner(trellis->rconseqscanner[i], 
							 trellis->rconseq->seq, hmm->cmmap[i]);

                                /* Pre-compute conseq scanners for EXPLICIT */
                                /* states as a "running sum". Required by   */
                                /* the implementation of EXPLICIT states    */
                                /* (see zScoreScanners in zTransition.c for */
                                /* details).                                */

				/*if (zUsedInExplicitState(hmm, i)){ EVAN removed this for speed*/
				zPreComputeScanner(trellis->conseqscanner[i]);
				zPreComputeScanner(trellis->rconseqscanner[i]);
				/*}*/

			}
		} else {
			zDie("non-conservation model in list of conservation models");
		}
	}

	zAllocFactories(trellis, conseq_enabled);
}

void zFreeTrellis (zTrellis *trellis) {
	int i;

	

	for (i = 0; i < trellis->hmm->feature_count; i++) {
		if (trellis->scanner[i]  != NULL) {
			zFreeScanner(trellis->scanner[i]);
			zFree(trellis->scanner[i]);
		}
		if (trellis->rscanner[i] != NULL) {
			zFreeScanner(trellis->rscanner[i]);
			zFree(trellis->rscanner[i]);
		}

		if (trellis->conseq) {
			if (trellis->conseqscanner[i] != NULL) {
				zFreeScanner(trellis->conseqscanner[i]);
				zFree(trellis->conseqscanner[i]);
			}
			if (trellis->rconseqscanner[i] != NULL) {
				zFreeScanner(trellis->rconseqscanner[i]);
				zFree(trellis->rconseqscanner[i]);
			}
		}
	}

	zFree(trellis->scanner);        trellis->scanner = NULL; 
	zFree(trellis->rscanner);       trellis->rscanner = NULL;
	zFree(trellis->conseqscanner);  trellis->conseqscanner = NULL;
	zFree(trellis->rconseqscanner); trellis->rconseqscanner = NULL;

	
	zFreeFactories(trellis);

//	zFreeViterbiVars(trellis);

//Added by Jim
	zFreeViterbiVarsFinal(trellis);

	zFreeForwardVars(trellis);
	zFreeBackwardVars(trellis);
	
	zFreeDNA(trellis->rdna);
	zFree(trellis->rdna);
	trellis->rdna = NULL;
	zFreeDNA(trellis->dna);
	zFree(trellis->dna);
	trellis->dna = NULL;

	if (trellis->conseq != NULL) {
		zFreeConseq(trellis->conseq);
		zFreeConseq(trellis->rconseq);
		zFree(trellis->conseq);
		zFree(trellis->rconseq);
	}
}

void zRunPartialViterbiAndForward(zTrellis *trellis, coor_t start, coor_t end){
	zHMM         *hmm = trellis->hmm;

	coor_t        i;        /* iterator for sequence */
	int           j;        /* iterator for internal states */
	int           k;        /* iterator for previous states */
	zIVec*        jumps;

	/* induction */
	zTrace2("beginning induction\n");
	for (i = start; i <= end; i++) {
		for (j = 0; j < hmm->states; j++) {
			if (NULL == trellis->cell[j]) continue;
			/* 			if (NULL == trellis->forward[j]) continue; */

			// trellis->cell[j][i].score[0] = MIN_SCORE;
			// trellis->cell[j][i].trace[0] = -1;
			// trellis->cell[j][i].length[0] = 0;

			trellis->cell[j][i].layer[0]->score = MIN_SCORE;
			trellis->cell[j][i].layer[0]->trace = -1;
			trellis->cell[j][i].layer[0]->length = 0;


			/* 			trellis->forward[j][i] = MIN_SCORE; */

			jumps = hmm->jmap[j];
			for(k = 0; k < jumps->size; k++) {
				zGetTransFunc(hmm->state[j].type)
					(trellis, jumps->elem[k], j, i);
			}
		}
	}
}

//Added by Jim
void zTreeterbiCollectNew(zTrellis* trellis,coor_t pos,int state,int story,int window,int depth)
{//want to del pos (already minused window in superior method)

	if (pos<PADDING || pos>trellis->dna->length-1-PADDING)
	{
//		printf("RANGE error %d\n",pos);
		return;
	}
	coor_t from_pos;

//	if (tmp_trace==-1)//temparory
//if(pos>11748)printf("AAA0\n");
	if (trellis->cell[pos][state].layer[story]==NULL)
	{
//		printf("self already %d,%s null %d\n",pos,zStrIdx2Char(trellis->hmm->state[state].name),depth);
		return;
	}
	if (trellis->cell[pos][state].layer[story]->lock==1)
	{
//		printf("self already %d,%s LOCK %d\n",pos,zStrIdx2Char(trellis->hmm->state[state].name),depth);
		return;
	}
	int tmp_length=trellis->cell[pos][state].layer[story]->length;
	int tmp_trace=trellis->cell[pos][state].layer[story]->trace;
	int tmp_trace_top=trellis->cell[pos][state].layer[story]->trace_top;
	if (tmp_trace==-1||tmp_length==0)//No parent no kid
	{		
		zFreeTreeNode(&trellis->cell[pos][state].layer[story]);
//		if (depth>0){	printf("DEPTH: %d\n",depth);}
//		printf("del %d,%s success %d\n",pos,zStrIdx2Char(trellis->hmm->state[state].name),depth);
		return;
	}
//fix for deleting padding in strict INTERNAL parent deletion strategy
	int w;
//	printf("www%d\n",w );

	// if ((trellis->hmm->state[state].type==INTERNAL)
	// 	|| (trellis->hmm->state[state].type==GINTERNAL)) 
	// {
	// 	//compensation or protection for not long enough window
	// 	if(trellis->cell[pos][state].layer[story]->length==1)
	// 	{
	// //		printf("protection\n");
	// 		zTBTreeLockNode(trellis->cell[pos][state].layer[story]);
	// //		printf("ppp %d\n",trellis->cell[pos][state].layer[story]->lock);
	// 		return;
	// 	}
	// }

	zTBTreeNode* layer_temp=trellis->cell[pos][state].layer[story]->parent;

	zFreeTreeNode(&trellis->cell[pos][state].layer[story]);
//	if (depth>0){	printf("DEPTH: %d\n",depth);}
//	printf("del %d,%s success %d\n",pos,zStrIdx2Char(trellis->hmm->state[state].name),depth);


	if(layer_temp==NULL)
	{
		zDie("zTreeterbiCollectNew parent null error");
//		printf("father already %d,%s null %d\n",pos,zStrIdx2Char(trellis->hmm->state[state].name),depth);
		return;
	}
	//Being more strict to Internal states, so annotate this part. Having a direct internal son doesn't count here 
	// if ((trellis->hmm->state[state].type==INTERNAL)
	// || (trellis->hmm->state[state].type==GINTERNAL)) 
	// {from_pos=pos-1;}
	// if (trellis->hmm->state[state].type==EXTERNAL)
	// 	{from_pos=pos-tmp_length;}

	from_pos=pos-tmp_length;
//if(pos>12047)	printf("TTT1\n");
 	if (layer_temp->children <= 0)
	{
 //		printf("start del %d,%s\n",from_pos,zStrIdx2Char(trellis->hmm->state[tmp_trace].name));
		zTreeterbiCollectNew(trellis,from_pos,tmp_trace,tmp_trace_top,window,depth+1);//recursion, could be change into iteration for safety
 	}
//if(pos>12047)	printf("TTT2\n");
}
//Added by Jim
//only for specific function
void zTreeterbiCollectDrawingNew(zTrellis* trellis,coor_t pos,int state,int story,int depth)
{//want to del pos (already minused window in superior method)

	if (pos<PADDING || pos>trellis->dna->length-1-PADDING)
	{
//		printf("RANGE error %d\n",pos);
		return;
	}
	coor_t from_pos;

//	if (tmp_trace==-1)//temparory
//if(pos>11748)printf("AAA0\n");
	if (trellis->cell[pos][state].layer[story]==NULL)
	{
//		printf("self already %d,%s null %d\n",pos,zStrIdx2Char(trellis->hmm->state[state].name),depth);
		return;
	}
	if (trellis->cell[pos][state].layer[story]->left_pos==1)
	{
		return;
	}
	int tmp_length=trellis->cell[pos][state].layer[story]->length;
	int tmp_trace=trellis->cell[pos][state].layer[story]->trace;
	int tmp_trace_top=trellis->cell[pos][state].layer[story]->trace_top;
	if (tmp_trace==-1||tmp_length==0)//No parent no kid
	{		
		if (pos==24264)
		{ 
			printf("NEGONE DELETE %d,%d,%d\n",pos,state,story);
		}
		zFreeTreeNode(&trellis->cell[pos][state].layer[story]);
//		if (depth>0){	printf("DEPTH: %d\n",depth);}
//		printf("del %d,%s success %d\n",pos,zStrIdx2Char(trellis->hmm->state[state].name),depth);
		return;
	}
//fix for deleting padding in strict INTERNAL parent deletion strategy
	int w;

	zTBTreeNode* layer_temp=trellis->cell[pos][state].layer[story]->parent;
//   printf("children-num%d\n", layer_temp->children);



	if (pos==24264)
	{ 
		printf("NORMAL DELETE %d,%d,%d\n",pos,state,story);
	}

	zFreeTreeNode(&trellis->cell[pos][state].layer[story]);
//	if (depth>0){	printf("DEPTH: %d\n",depth);}
//	printf("del %d,%s success %d\n",pos,zStrIdx2Char(trellis->hmm->state[state].name),depth);


	if(layer_temp==NULL)
	{
		
		printf("Father already %d,%s null %d\n",pos,zStrIdx2Char(trellis->hmm->state[state].name),depth);
	//	zDie("zTreeterbiCollectDrawingNew parent null error");
		return;
	}
	from_pos=pos-tmp_length;
//if(pos>12047)	printf("TTT1\n");
//	printf("children-num%d\n", layer_temp->children);
 	if (layer_temp->children <= 0)
	{
// 		printf("Drawing inside start del %d,%s\n",from_pos,zStrIdx2Char(trellis->hmm->state[tmp_trace].name));
		zTreeterbiCollectDrawingNew(trellis,from_pos,tmp_trace,tmp_trace_top,depth+1);//recursion, could be change into iteration for safety
 	}
//if(pos>12047)	printf("TTT2\n");
}
void zTreeterbiCollectDrawingLast(zTrellis* trellis,coor_t pos,int state,int story,int depth)
{//want to del pos (already minused window in superior method)

	if (pos<PADDING || pos>trellis->dna->length-1-PADDING)
	{
//		printf("RANGE error %d\n",pos);
		return;
	}
	coor_t from_pos;

//	if (tmp_trace==-1)//temparory
//if(pos>11748)printf("AAA0\n");
	if (trellis->cell[pos][state].layer[story]==NULL)
	{
//		printf("self already %d,%s null %d\n",pos,zStrIdx2Char(trellis->hmm->state[state].name),depth);
		return;
	}
	if (trellis->cell[pos][state].layer[story]->left_pos==1)
	{
		return;
	}
	int tmp_length=trellis->cell[pos][state].layer[story]->length;
	int tmp_trace=trellis->cell[pos][state].layer[story]->trace;
	int tmp_trace_top=trellis->cell[pos][state].layer[story]->trace_top;
	if (tmp_trace==-1||tmp_length==0)//No parent no kid
	{		
		zFreeTreeNode(&trellis->cell[pos][state].layer[story]);
//		if (depth>0){	printf("DEPTH: %d\n",depth);}
//		printf("del %d,%s success %d\n",pos,zStrIdx2Char(trellis->hmm->state[state].name),depth);
		return;
	}
//fix for deleting padding in strict INTERNAL parent deletion strategy
	int w;

	zTBTreeNode* layer_temp=trellis->cell[pos][state].layer[story]->parent;
//printf("children-num%d\n", layer_temp->children);
	zFreeTreeNode(&trellis->cell[pos][state].layer[story]);
//	if (depth>0){	printf("DEPTH: %d\n",depth);}
//	printf("del %d,%s success %d\n",pos,zStrIdx2Char(trellis->hmm->state[state].name),depth);


	if(layer_temp==NULL)
	{
		zDie("zTreeterbiCollectDrawingLast parent null error");
	//	printf("father already %d,%s null %d\n",pos,zStrIdx2Char(trellis->hmm->state[state].name),depth);
		return;
	}
	from_pos=pos-tmp_length;
//if(pos>12047)	printf("TTT1\n");
//	printf("children-num%d\n", layer_temp->children);
 	if (layer_temp->left_pos != 1)
	{
 //		printf("Last Drawing inside start del %d,%s\n",from_pos,zStrIdx2Char(trellis->hmm->state[tmp_trace].name));
		zTreeterbiCollectDrawingLast(trellis,from_pos,tmp_trace,tmp_trace_top,depth+1);//recursion, could be change into iteration for safety
 	}
//if(pos>12047)	printf("TTT2\n");
}
//void zFreeTrellisCellLayer (zTrellisCellLayer **p) {
 void zFreeTreeNode(zTBTreeNode** p){
	if (*p != NULL)
	{
//		printf("before FREE%d\n", (*p)->state);
//		printf("before FREE%d\n",*p);
	//	printf("DEL:%d,%f,%d,%d,%d\t",(*p)->pos,(*p)->score,(*p)->length,(*p)->state,(*p)->story);
		if((*p)->sfv!=NULL)
		{
			printf("SFV:%d,%f,%d,%d,%d\n",(*p)->pos,(*p)->score,(*p)->length,(*p)->state,(*p)->story);
			if((*p)->sfv->elem==NULL)printf("CONG\n");

			printf("SFV: %d\n",(*p)->sfv->size);
		//	printf("%d\n",(*p)->sfv->elem[0].state );
			zFreeSFVec((*p)->sfv);// here problem
			printf("SFV2\n");
			if((*p)->sfv->elem==NULL)printf("CONG\n");
			zFree((*p)->sfv);//new added
			printf("SFV3\n");
		}

 		zDeleteTBTreeNode((*p));//not real delete only delete relationship
//printf("mid FREE%d\n", *p);
		zFree(*p);
	//	free(p);
	//	*p=NULL	//printf("after FREE%d\n", (*p)->state);
		*p = NULL;

		
	}
 // zDeleteTBTreeNode(*p);
 // zFree(*p);

}

void zTreeterbiNew(zTrellis* trellis, coor_t pos, int window, coor_t root_pos)
{


	int j,k;
	int tmp_pos=pos-window;
	if (window==0)
	{
		return;
	}

//	if ((tmp_pos < PADDING) ||(tmp_pos > trellis->dna->length-PADDING-1))
	if ((tmp_pos < root_pos+window) ||(tmp_pos > trellis->dna->length-PADDING-1))
	{
//		printf("Treeterbi error\n");
		return;
	}

	for (j = 0; j < trellis->hmm->states; j++)
	{
		for (k = 0; k < top_num; k++)
		{
//			printf("TTTT0\n");
			if(trellis->cell[tmp_pos][j].layer[k]->trace==-1 || trellis->cell[tmp_pos][j].layer[k]->trace_top==-1|| trellis->cell[tmp_pos][j].layer[k]->length==0)//temparory replace the real deletion
			{
				// if (trellis->cell[tmp_pos][j].layer[k]->parent==NULL)
				// {
				// 	printf("TTTT0.5\n");
				// }				
		//		printf("tree start EMPTY del %d,%s\n",tmp_pos,zStrIdx2Char(trellis->hmm->state[j].name));
	//			printf("tree start EMPTY del %d,%d,%d\n",tmp_pos,j,k);
				zTreeterbiCollectNew(trellis,tmp_pos,j,k,window,0);
//				if(tmp_pos>13326)printf("EMPTY out\n");
				continue;}
//		printf("TTTT1\n");
			if (trellis->cell[tmp_pos][j].layer[k]->children <= 0)
			{
//				printf("TTTT2\n");
	//			printf("tree start del %d,%s\n",tmp_pos,zStrIdx2Char(trellis->hmm->state[j].name));
	//			printf("tree start del %d,%d,%d\n",tmp_pos,j,k);
				zTreeterbiCollectNew(trellis,tmp_pos,j,k,window,0);//delete the cell and his parents
//				printf("TTTT3\n");
			}
		}
	}
}

//Added by Jim
//only for spoecific function
void zTreeterbiForDrawing(zTrellis* trellis,coor_t term_pos)
{
	int i,j,k;
//	for (i = PADDING; i <= trellis->dna->length-PADDING-1; ++i)//orignially is < not <=
	for (i = PADDING; i <= term_pos; ++i)
	{
		if (trellis->cell[i]==NULL)
		{
			continue;
		}
		for (j = 0; j < trellis->hmm->states; ++j)
		{
			if (trellis->cell[i][j].layer==NULL)
			{
				continue;
			}
			for (k = 0; k < top_num; ++k)
			{
				if(trellis->cell[i][j].layer[k]==NULL){continue;}

				if(trellis->cell[i][j].layer[k]->trace==-1 || trellis->cell[i][j].layer[k]->trace_top==-1|| trellis->cell[i][j].layer[k]->length==0)//temparory replace the real deletion
				{
		//			printf("drawing start EMPTY del %d,%s,%d\n",i,zStrIdx2Char(trellis->hmm->state[j].name),k);
					zTreeterbiCollectDrawingNew(trellis,i,j,k,0);
					continue;}

				//extra
				if(i==term_pos)
				{
					if (trellis->cell[i][j].layer[k]->left_pos!=1)
					{
		//				printf("drawing start END del %d,%s,%d\n",i,zStrIdx2Char(trellis->hmm->state[j].name),k);
						if (trellis->cell[i][j].layer[k]->parent!=NULL)
						{
		//					printf("parent exist\n");
						}
						zTreeterbiCollectDrawingNew(trellis,i,j,k,0);
				//		zTreeterbiCollectDrawingLast(trellis,i,j,k,0);//delete the cell and his parents
					}
					continue;
				}


				if (trellis->cell[i][j].layer[k]->children <= 0)
				{
		//			printf("drawing start del %d,%s,%d\n",i,zStrIdx2Char(trellis->hmm->state[j].name),k);
					zTreeterbiCollectDrawingNew(trellis,i,j,k,0);//delete the cell and his parents
				}
			}
		}
	}

}

void zExtraClean(zTrellis *trellis)//Brutally clean only use for temporary use
{
	int i,j,k;
	for (i = PADDING; i <trellis->dna->length-PADDING-1; ++i)//orignially is < not <=
	{
		if (trellis->cell[i]==NULL)
		{
			continue;
		}
		for (j = 0; j < trellis->hmm->states; ++j)
		{
			if (trellis->cell[i][j].layer==NULL)
			{
				continue;
			}
			for (k = 0; k < top_num; ++k)
			{
				if(trellis->cell[i][j].layer[k]==NULL){continue;}
				if(trellis->cell[i][j].layer[k]->left_pos!=1)
				{
					zFreeTreeNode(&trellis->cell[i][j].layer[k]);
				}
			}
		}
	}

}


void zConnectTBTreeLayers(zTBTree *tree,zTrellis *trellis,coor_t pos,int state, coor_t root_pos)
{
	int k;
	coor_t tmp_pos;
	int tmp_trace,tmp_trace_top;
	for (k = 0; k < top_num; k++)
	{
		if (trellis->cell[pos][state].layer[k]==NULL)
		{
			continue;
		}
		tmp_pos=pos-trellis->cell[pos][state].layer[k]->length;
		tmp_trace=trellis->cell[pos][state].layer[k]->trace;
		tmp_trace_top=trellis->cell[pos][state].layer[k]->trace_top;
		if (tmp_pos==pos || tmp_trace==-1 || tmp_trace_top==-1)
		{
			continue;
		}

		trellis->cell[pos][state].layer[k]->pos = pos;
		trellis->cell[pos][state].layer[k]->state = state;
		trellis->cell[pos][state].layer[k]->story = k;


		if ((trellis->hmm->state[state].type==INTERNAL)
			|| (trellis->hmm->state[state].type==GINTERNAL)) 
		{
			//compensation or protection for not long enough window
			if(trellis->cell[pos][state].layer[k]->length==1)
			{
	//			printf("protection\n");
				zTBTreeLockNode(trellis->cell[pos][state].layer[k]);
		//		printf("ppp %d\n",trellis->cell[pos][state].layer[story]->lock);
			}
		}



//		printf("%d,%d\n",state,k );
//		if (tmp_pos==PADDING-1)
		if(tmp_pos==root_pos)
		{
			zTBTreeSetChild(tree->root,trellis->cell[pos][state].layer[k]);
//			printf("%d,%d,Parent: %d\tChild: %d\n",pos,state,tree->root->id,trellis->cell[state][pos].layer[k]->id );
		}
		else
		{
			if (trellis->cell[tmp_pos][tmp_trace].layer[tmp_trace_top]==NULL)
			{
				printf("Need Hand shake!!!\n");
				printf("HANDSHAKE: %s,%d,%d\n",zStrIdx2Char(trellis->hmm->state[tmp_trace].name),tmp_pos,tmp_trace_top);
				zHandShake(trellis,state,pos,k);

			}
			else{
				zTBTreeSetChild(trellis->cell[tmp_pos][tmp_trace].layer[tmp_trace_top],trellis->cell[pos][state].layer[k]);
			}
			
//			printf("%d,%d,Parent: %d\tChild: %d\n",pos,state,trellis->cell[tmp_trace][tmp_pos].layer[tmp_trace_top]->id,trellis->cell[state][pos].layer[k]->id );
		}
	}
}

void zCountCPoint(zTrellis *trellis, coor_t posi,int window)
{
	int pos=posi-window;
	if (pos<50)
	{
		printf("Cpoint error\n");
		return;
	}
	
	int i,j,k;
	int nullCount=0;
	int inCount=0;
	int exCount=0;
	int cpointCount=0;
	int matched=0;
	int former_pos_set[49];
	int trace_set[49];
	int trace_top_set[49];
	int tmp_pos;
	int tmp_trace;
	int tmp_trace_top;

	for (i = 0; i < 49; ++i)
	{
		former_pos_set[i]=-1;
		trace_set[i]=-1;
		trace_top_set[i]=-1;
	}
	
	for (j = 0; j < 49; ++j)
	{
		i=0;
		matched=0;
	//	printf("AAA:%d,%d\n",j,pos);
		if(trellis->cell[j][pos].layer[0]==NULL){nullCount++;continue;}
//printf("AAA1\n");
		if(trellis->cell[j][pos].layer[0]->trace==-1 || trellis->cell[j][pos].layer[0]->trace_top==-1 || trellis->cell[j][pos].layer[0]->length==0)
		{
			nullCount++;
			continue;
		}
		tmp_pos=pos-trellis->cell[j][pos].layer[0]->length;
		tmp_trace=trellis->cell[j][pos].layer[0]->trace;
		tmp_trace_top=0;
		// if (j<17)
		// {
		// 	tmp_pos=pos-1;
		// 	tmp_trace=j;
		// 	tmp_trace_top=0;
		// }
		// else{
		// 	tmp_pos=pos-trellis->cell[j][pos].layer[0]->length;
		// 	tmp_trace=trellis->cell[j][pos].layer[0]->trace;
		// 	tmp_trace_top=0;
		// }
		if (trace_set[i]<17)
		{
			inCount++;

		}else{
			exCount++;

		}
		
		while(former_pos_set[i]!=-1)
		{
			if((former_pos_set[i]==tmp_pos)&&(trace_set[i]==tmp_trace)&&(trace_top_set[i]==tmp_trace_top))
			{
				matched=1;

				break;
			}
			i++;
		}
		if (matched==0)
		{
			former_pos_set[i]=tmp_pos;
			trace_set[i]=tmp_trace;
			trace_top_set[i]=tmp_trace_top;
			cpointCount++;
		}
	}
	if (cpointCount<10)
	{
		printf("%d: n-%d, i-%d, e-%d, c-%d LOL\n",pos,nullCount,inCount,exCount,cpointCount );
		return;
	}
	printf("%d: n-%d, i-%d, e-%d, c-%d\n",pos,nullCount,inCount,exCount,cpointCount );
}
void zCountCPointNew(zTrellis *trellis, coor_t posi,int window, coor_t after)
{
	int pos=posi-window;
	if (pos<50)
	{
		printf("Cpoint error\n");
		return;
	}
	
	int i,j,k;
	int nullCount=0;
	int surpassInterCount=0;
	int surpassExterCount=0;

	int cpointCount=0;
	int laterCount=0;
	int back_set[trellis->hmm->states][top_num];
	int tmp_pos;
	int tmp_trace;
	int tmp_trace_top;

	for (j = 0; j < trellis->hmm->states; ++j)
	{
		for (k = 0; k < top_num; ++k)
		{
			back_set[j][k]=0;
		}
	}
	for (i = pos+1; i <= pos+after; ++i)
	{
		for (j = 0; j < trellis->hmm->states; ++j)
		{
			for (k = 0; k < top_num; ++k)
			{
			//	printf("AAA:%d,%d\n",j,pos);
				if(trellis->cell[j][i].layer[k]==NULL){nullCount++;continue;}
		//printf("AAA1\n");
				if(trellis->cell[j][i].layer[k]->trace==-1 || trellis->cell[j][i].layer[k]->trace_top==-1 || trellis->cell[j][i].layer[k]->length==0)
				{
					nullCount++;
					continue;
				}
				tmp_pos=i-trellis->cell[j][i].layer[k]->length;
				if (tmp_pos<pos && j<17)
				{
					surpassInterCount++;//problem had, should be changed soon
				}
				if (tmp_pos<pos && j>=17)
				{
					surpassExterCount++;
					continue;
				}
				if(tmp_pos>pos)
				{
					laterCount++;
					continue;
				}
				tmp_trace=trellis->cell[j][i].layer[k]->trace;
				tmp_trace_top=trellis->cell[j][i].layer[k]->trace_top;
				if (back_set[tmp_trace][tmp_trace_top]==0)
				{
					cpointCount++;
				}
				back_set[tmp_trace][tmp_trace_top]++;
			}
		}
	}
	if (cpointCount<=top_num && cpointCount > 0 && surpassExterCount == 0)
	{
		printf("GET %d\n",pos);
		// if (cpointCount==1)
		// {
		// 	printf("POK\n");
		// 	for (j = 0; j < trellis->hmm->states; ++j)
		// 	{
		// 		for (k = 0; k < top_num; ++k)
		// 		{
		// 			if (back_set[j][k]==0) 
		// 			{
		// 				continue;
		// 			}
		// 			if(back_set[j][k]!=0)
		// 			{

		// 				printf("PPP %d,%s,%d thread:%d\n",pos-trellis->cell[j][pos].layer[k]->length,zStrIdx2Char(trellis->hmm->state[trellis->cell[j][pos].layer[k]->trace].name),trellis->cell[j][pos].layer[k]->trace_top,back_set[j][k]);
		// 				break;
		// 			}
		// 		}
		// 	}
		// }
	}
	if (cpointCount > 0)
	{
		printf("%d: n-%d, sI-%d, sE-%d, l-%d, c-%d LOL\n",pos,nullCount,surpassInterCount,surpassExterCount,laterCount,cpointCount );
		return;
	}
	printf("%d: n-%d, s-%d, sE-%d, l-%d, c-%d\n",pos,nullCount,surpassInterCount,surpassExterCount,laterCount,cpointCount );
}

int zPowerInt(int x, int y)
{
	int i;
	int result=1;
	for (i = 0; i < y; ++i)
	{
		result*=x;
	}
	return result;
}
//should be in zSfeature
score_t zScoreSFVec(zTrellis *trellis,zSFVec *vec)
{
	int i;
	score_t t_score;
	if (vec->size==0)
	{
		return MIN_SCORE;
	}
	t_score = vec->elem[0].score;
	for (i = 1; i < vec->size; ++i)
	{
		t_score += zGetTransitionScore(trellis->hmm,vec->elem[i-1].state,vec->elem[i].state,trellis->tiso_group)
			+ vec->elem[i].score;
	}
	return t_score;
}


void zForTmpSFVScoreForSponge(zTrellis *trellis,zSFVec **huge_sfv,zSFVSpongeNode* index, score_t* score_box)
{	
	int i,j;
	score_t t_score;
	int tmp_sfv_size=zCountRealSFV(index->node,index->size);

	for (i = 0; i < tmp_sfv_size; ++i)
	{
		if (index->node[i]->size==0)
		{
			t_score=MIN_SCORE;
		}
		else{
			t_score = zGetTransitionScore(trellis->hmm,huge_sfv[0]->elem[huge_sfv[0]->size-1].state,index->node[i]->elem[0].state,trellis->tiso_group)
				+ index->node[i]->elem[0].score;
		

			for (j = 1; j < index->node[i]->size; ++j)
			{

				t_score += zGetTransitionScore(trellis->hmm,index->node[i]->elem[j-1].state,index->node[i]->elem[j].state,trellis->tiso_group)
						+ index->node[i]->elem[j].score;
			}
		}
		score_box[i]=t_score;
		printf("SCORE_BOX %d,%f\n",i,score_box[i]);
	}
}

void zForTmpSFVScore(zTrellis *trellis,zSFVec **huge_sfv,zSFVec **tmp_sfv, score_t* score_box)
{	
	int i,j;
	score_t t_score;
//	printf("www\n");
	int tmp_sfv_size=zCountRealSFV(tmp_sfv,top_num);
//	printf("yyy %d\n",tmp_sfv_size);
//	for (i = 0; i < top_num; ++i)
	for (i = 0; i < tmp_sfv_size; ++i)
	{
//		printf("zzz\n");
//		printf("%d\n",tmp_sfv[i]->size );
		if (tmp_sfv[i]->size==0)
		{
			t_score=MIN_SCORE;
	//		printf("HUGE_ERROR zForTmpSFVScore\n");
		}
		else{
//			printf("before T_SCORE %d,%f\n",i,t_score);
			t_score = zGetTransitionScore(trellis->hmm,huge_sfv[0]->elem[huge_sfv[0]->size-1].state,tmp_sfv[i]->elem[0].state,trellis->tiso_group)
				+ tmp_sfv[i]->elem[0].score;

//			printf("after T_SCORE %d,%f\n",i,t_score);
			for (j = 1; j < tmp_sfv[i]->size; ++j)
			{
				t_score += zGetTransitionScore(trellis->hmm,tmp_sfv[i]->elem[j-1].state,tmp_sfv[i]->elem[j].state,trellis->tiso_group)
						+ tmp_sfv[i]->elem[j].score;
//				printf("loop T_SCORE %d,%f   %d\n",i,t_score,j);
			}
		}
		score_box[i]=t_score;
		printf("SCORE_BOX %d,%f\n",i,score_box[i]);
	}
}

void zQuickSort(score_t *score, int low, int high,int *mark)
{
	if (low==0)
	{
		printf("Inside zQuickSort\n");
	}
	
	if (low>=high){return;}
	int first = low;
	int last = high;
	score_t key=score[first];
	int loc=mark[first];
	while(first < last)
	{
		while(first<last && score[last]<=key)
		{
			last--;
		}
		score[first]=score[last];
		mark[first]=mark[last];
		while(first<last && score[first]>=key)
		{
			first++;
		}
		score[last]=score[first];
		mark[last]=mark[first];
	}	
	score[first]=key;
	mark[first]=loc;
	zQuickSort(score,low,first-1,mark);
	zQuickSort(score,first+1,high,mark);
}

void zQuickSortSFVec(score_t *score, int low, int high,zSFVec **mark)
{
	if (low==0)
	{
		printf("Inside zQuickSortSFVec\n");
	}
	if (low>=high){return;}
	int first = low;
	int last = high;
	score_t key=score[first];
	zSFVec* loc=mark[first];
	while(first < last)
	{
		while(first<last && score[last]<=key)
		{
			last--;
		}
		score[first]=score[last];
		mark[first]=mark[last];
		while(first<last && score[first]>=key)
		{
			first++;
		}
		score[last]=score[first];
		mark[last]=mark[first];
	}	
	score[first]=key;
	mark[first]=loc;
	zQuickSortSFVec(score,low,first-1,mark);
	zQuickSortSFVec(score,first+1,high,mark);
}

void zSortScore(score_t *score, int *mark,int size)
{
	printf("Inside zSortScore\n");
	int i,j;
	score_t medium;
	int index;
	for (i = 0; i < size; ++i)
	{
		for (j = i+1; j < size; ++j)
		{
			if (score[i]<score[j])
			{
				medium=score[i];
				score[i]=score[j];
				score[j]=medium;
				index=mark[i];
				mark[i]=mark[j];
				mark[j]=index;
			}
		}
	}
}
void zSortScoreSFVec(score_t *score, zSFVec **mark,int size)
{
	printf("Inside zSortScoreSFV\n");
	int i,j;
	score_t medium;
	zSFVec* index;
	for (i = 0; i < size; ++i)
	{
		for (j = i+1; j < size; ++j)
		{
			if (score[i]<score[j])
			{
				medium=score[i];
				score[i]=score[j];
				score[j]=medium;
				index=mark[i];
				mark[i]=mark[j];
				mark[j]=index;
			}
		}
	}
}

void zModifyScoreForSponge(zTrellis *trellis,zSFVec **huge_sfv,zSFVSpongeNode* index,score_t *score_cal,score_t* score_box,int pre_state)
{	
	printf("MS_Pre_state %d\n",pre_state);
	int i,j,k;
	zHMM *hmm=trellis->hmm;
	score_t t_score=0.0;	
	int tmp_sfv_size=zCountRealSFV(index->node,index->size);

	int huge_sfv_size=zCountRealSFV(huge_sfv,LARGE_TOP);

//	printf("GGG1\n");
	for (i = 0; i < huge_sfv_size; ++i)
	{	// printf("GGG: %d\n",i);

		if (pre_state == -1)
		{
			if(hmm->state[huge_sfv[i]->elem[0].state].type == GINTERNAL)
			{
				t_score=huge_sfv[i]->elem[0].score
					+ zGetInitProb(trellis->hmm, huge_sfv[i]->elem[0].state, trellis->iiso_group)
					+ zGetIntergenicContinueScore(trellis);
			}
			else if (hmm->state[huge_sfv[i]->elem[0].state].type == INTERNAL)
			{
				t_score=huge_sfv[i]->elem[0].score
					+ zGetInitProb(trellis->hmm, huge_sfv[i]->elem[0].state, trellis->iiso_group)
					+ zGetTransitionScore(trellis->hmm,huge_sfv[i]->elem[0].state,huge_sfv[i]->elem[0].state,trellis->tiso_group);
			}
			else{t_score=MIN_SCORE;printf("ForSponge Error\n");}
		}
		else{			
				t_score = zGetTransitionScore(trellis->hmm,pre_state,huge_sfv[i]->elem[0].state,trellis->tiso_group)
					+ huge_sfv[i]->elem[0].score;
		}
		
		for (j = 1; j < huge_sfv[i]->size; ++j)
		{
			t_score += zGetTransitionScore(trellis->hmm,huge_sfv[i]->elem[j-1].state,huge_sfv[i]->elem[j].state,trellis->tiso_group)
				+ huge_sfv[i]->elem[j].score;
		}
		for (k = 0; k < tmp_sfv_size; ++k)
		{
			score_cal[i*tmp_sfv_size+k] = t_score+score_box[k];
		}
	}

}
void zModifyScore(zTrellis *trellis,zSFVec **huge_sfv,zSFVec **tmp_sfv,score_t *score_cal,score_t* score_box)
{	
	int i,j,k;
	zHMM *hmm=trellis->hmm;
	score_t t_score=0.0;	
	int tmp_sfv_size=zCountRealSFV(tmp_sfv,top_num);

	int huge_sfv_size=zCountRealSFV(huge_sfv,LARGE_TOP);

//	printf("GGG1\n");
	for (i = 0; i < huge_sfv_size; ++i)
	{	// printf("GGG: %d\n",i);
//printf("%d\n", huge_sfv[i]->elem[0].state);
		if(hmm->state[huge_sfv[i]->elem[0].state].type == GINTERNAL)
		{
			t_score=huge_sfv[i]->elem[0].score
				+ zGetInitProb(trellis->hmm, huge_sfv[i]->elem[0].state, trellis->iiso_group)
				+ zGetIntergenicContinueScore(trellis);
		}
		else if (hmm->state[huge_sfv[i]->elem[0].state].type == INTERNAL)
		{
			t_score=huge_sfv[i]->elem[0].score
				+ zGetInitProb(trellis->hmm, huge_sfv[i]->elem[0].state, trellis->iiso_group)
				+ zGetTransitionScore(trellis->hmm,huge_sfv[i]->elem[0].state,huge_sfv[i]->elem[0].state,trellis->tiso_group);
		}
		else{
	//		printf("SIZE %d\n",huge_sfv[i]->size );
	//		printf("%d\n", huge_sfv[i]->elem[0].state);
			zDie("Modifiying score error");
		}
	//	t_score=qhuge_sfv[i]->elem[0].score+ zGetInitProb(trellis->hmm, huge_sfv[i]->elem[0].state, trellis->iiso_group)+zGetTransitionScore(trellis->hmm,huge_sfv[i]->elem[0].state,huge_sfv[i]->elem[0].state,trellis->tiso_group);
	//	printf("PECK %d: %f,%f,%f\n",i+1,t_score,zGetInitProb(trellis->hmm, huge_sfv[i]->elem[0].state, trellis->iiso_group),zGetIntergenicContinueScore(trellis));//zGetTransitionScore(trellis->hmm,huge_sfv[i]->elem[0].state,huge_sfv[i]->elem[0].state,trellis->tiso_group));
//printf("GGG2   :%d\n",huge_sfv[i]->size);


		for (j = 1; j < huge_sfv[i]->size; ++j)
		{
//		printf("GGG2.1: %d\n",j);
//			printf("\t %d,%d,%d\n",huge_sfv[i]->elem[j-1].state,huge_sfv[i]->elem[j-1].end,huge_sfv[i]->elem[j].score );
	//		printf("\t %d,%d,%d\n",huge_sfv[i]->elem[j].state,huge_sfv[i]->elem[j].end,huge_sfv[i]->elem[j].score );
			t_score += zGetTransitionScore(trellis->hmm,huge_sfv[i]->elem[j-1].state,huge_sfv[i]->elem[j].state,trellis->tiso_group)
				+ huge_sfv[i]->elem[j].score;
		}
//		printf("GGG3\n");
	//	for (k = 0; k < top_num; ++k)
		for (k = 0; k < tmp_sfv_size; ++k)
		{
		//	score_cal[i*top_num+k] = t_score+score_box[k];
			score_cal[i*tmp_sfv_size+k] = t_score+score_box[k];
		//	score_cal[k*huge_sfv_size+i] = t_score+score_box[k];

		}
//printf("GGG4\n");
	}

}

int zCountRealSFV(zSFVec **sfv, int sfv_limit)
{
	int i;
	for (i = 0; i < sfv_limit; ++i)
	{
		if (sfv[i]->size==0)
		{
			break;
		}
	}
	return i;

}

void zModifySortScoreForSponge(zTrellis *trellis,zSFVec **huge_sfv,score_t* path_score,int large_top,int type, int pre_state,int mode);
void zTracePartialTrellisMergeForSponge(zTrellis *trellis,zSFVec **huge_sfv,zSFVSpongeNode *first,zSFVSpongeNode *last,int large_top,score_t* score, int type)
{

	int i,j,k,g;
	int size;
	score_t tmp_score;
	int huge_sfv_size=0;
//	printf("Large_SIZE %d\n",huge_sfv_size );
	int potential_size;
	int total_kinds=1;
	int mark=0;
	int length = last->location - first->location +1;
	zSFVSpongeNode* index;
	index = first;
	

	int mode = 0;

	int pre_state = first->pre_state;

	while(index!=NULL)
	{
		total_kinds*=index->size;
		if (mark!=3)
		{
			potential_size=total_kinds;
		}
		else{potential_size=large_top*index->size;}
		
		if (potential_size < large_top)
		{
			if (mark==0)
			{
				for (i = 0; i < potential_size; ++i)
				{
			//		printf("tmp_sfv[%d] to large_sfv[%d]\n",i,i);
					zMergeSFVec(index->node[i],huge_sfv[i]);	
				}
				huge_sfv_size=potential_size;
				mark=1;
			}
			else
			{


				for (k = huge_sfv_size; k < potential_size ; ++k)
				{
			//		printf("large_sfv[%d] to large_sfv[%d]\n",k-huge_sfv_size,k);
					zMergeSFVec(huge_sfv[k-huge_sfv_size],huge_sfv[k]);	
				}
				for (i = 0; i < index->size; ++i)
				{
					for (k = 0; k < huge_sfv_size; ++k)
					{
			//			printf("tmp_sfv[%d] to large_sfv[%d]\n",i,i*huge_sfv_size+k);
						zMergeSFVec(index->node[i],huge_sfv[i*huge_sfv_size+k]);
					}
				}
				huge_sfv_size=potential_size;
				mark=2;


			}
		}
		else
		{
			mark=3;
			printf("BIG ININ %d\n",potential_size);
			score_t score_box[index->size];
			score_t score_cal[potential_size];
			int score_mark[potential_size];
			score_t tmp;
			zSFVec **tmp_huge_sfv;
			tmp_huge_sfv = zCalloc(large_top,sizeof(zSFVec*),"zRunViterbiAndForward sfv"); // waste storage soon correct
			size = sizeof(zSFVec);
			for(g=0;g<large_top;g++){
				tmp_huge_sfv[g]=zMalloc(size,"zMalloc tmp_huge_sfv[i]");
				zInitSFVec(tmp_huge_sfv[g],10);
			}
			for (i = 0; i < potential_size; ++i)
			{
				score_mark[i]=i;
			}
	//		printf("oops\n");
			zForTmpSFVScoreForSponge(trellis,huge_sfv,index,score_box);

	//		printf("what\n");
			zModifyScoreForSponge(trellis,huge_sfv,index,score_cal,score_box,pre_state);
		//	printf("the\n");
		//	zSortScore(score_cal,score_mark,LARGE_TOP*top_num);
			zQuickSort(score_cal,0,potential_size-1,score_mark);
		//	zSortScore(score_cal,score_mark,potential_size);
	//		printf("Fuck\n");
			for (k = 0; k < large_top; ++k)
			{
				score[k] = score_cal[k];
				i = score_mark[k]/index->size;
				j = score_mark[k]%index->size;
				// printf("index: %d\n",k );
	//		 	printf("large_sfv[%d] and tmp_sfv[%d]\n",i,j);
				zMergeSFVec(huge_sfv[i],tmp_huge_sfv[k]);
				zMergeSFVec(index->node[j],tmp_huge_sfv[k]);
	//			printf("you\n");

			}
	//		printf("doing\n");
		//	huge_sfv=tmp_huge_sfv;
			for (i = 0; i < large_top; ++i)
			{
				zFreeSFVec(huge_sfv[i]);
				zFree(huge_sfv[i]);
				huge_sfv[i]=tmp_huge_sfv[i];
			}
		}
		if (index==last)
		{
	//		if (mark<3)
	//		{
				printf("remember\n");
				zModifySortScoreForSponge(trellis,huge_sfv,score,large_top,type,pre_state,mode);
				printf("me\n");
	//		}
			break;
		}
		index = index->next;
		
	}
	return;	
	

//	printf("merge fin\n");
}

void zTracePartialTrellisMerge(zTrellis *trellis,zSFVec **huge_sfv,zSFVec **tmp_sfv, int mark,int large_top,score_t* score)
{
	int last=0;
	if (mark<0)
	{
		last=1;
		mark=-mark;
	}
	int i,j,k,g;
	int size;
	score_t tmp_score;
	int gap=zPowerInt(top_num,mark-1);
	int tmp_loc;

////////////////////////////////

	// for (i = 0; i < top_num; ++i)
	// {
	// 	printf("SHOW %d: ",i );
	// 	for (j = 0; j < tmp_sfv[i]->size; ++j)
	// 	{
	// 		printf("%d,%d-%d \n",tmp_sfv[i]->elem[j].state,tmp_sfv[i]->elem[j].start,tmp_sfv[i]->elem[j].end );
	// 	}
	// 	printf("\n");
	// }



	int huge_sfv_size=zCountRealSFV(huge_sfv,large_top);
	printf("HUGE_SIZE %d\n",huge_sfv_size );
	int tmp_sfv_size=zCountRealSFV(tmp_sfv,top_num);
	printf("TMP_SIZE %d\n",tmp_sfv_size );
	int potential_size=huge_sfv_size*tmp_sfv_size;
	printf("POTENT_SIZE %d\n",potential_size );

	for (k = 0; k < tmp_sfv_size; ++k)
	{
		zFlipSFVec(tmp_sfv[k]);
	}

	if (huge_sfv_size==0)
	{
		for (k = 0; k < tmp_sfv_size; ++k)
		{
//			printf("tmp_sfv[%d] to huge_sfv[%d]\n",k,k);

			// zMergeSFVec(tmp_sfv[k],huge_sfv[k]);
			// tmp_loc=tmp_sfv[k]->size-1;	
			// score[k]=tmp_sfv[k]->elem[tmp_loc].score;

			zMergeSFVec(tmp_sfv[k],huge_sfv[k]);
			tmp_loc=huge_sfv[k]->size-1;	
			score[k]=huge_sfv[k]->elem[tmp_loc].score;	

		//	score[k]=huge_sfv[k]->elem[tmp_loc].score;
		}
		return;
	}
	if(potential_size <= large_top)
	{
		for (k = huge_sfv_size; k < potential_size ; ++k)
		{
//			printf("huge_sfv[%d] to huge_sfv[%d]\n",k-huge_sfv_size,k);
			zMergeSFVec(huge_sfv[k-huge_sfv_size],huge_sfv[k]);	
			tmp_loc=huge_sfv[k-huge_sfv_size]->size-1;
		//	score[k]+=0.0;
			score[k]=huge_sfv[k-huge_sfv_size]->elem[tmp_loc].score;
		//	zConcatenateSFVec(huge_sfv[k-top_num],huge_sfv[k]);

		}
		for (i = 0; i < tmp_sfv_size; ++i)
		{
			for (k = 0; k < huge_sfv_size; ++k)
			{
	//			printf("tmp_sfv[%d] to huge_sfv[%d]\n",i,i*huge_sfv_size+k);
				zMergeSFVec(tmp_sfv[i],huge_sfv[i*huge_sfv_size+k]);
				tmp_loc=huge_sfv[i*huge_sfv_size+k]->size-1;
				score[k]=huge_sfv[i*huge_sfv_size+k]->elem[tmp_loc].score;
			}
		}
	}
//	if (last==1){return;}  //mark for the last time merge, no need to continue
	if (potential_size > large_top)
	{
		printf("BIG WARN\n");
		score_t score_box[tmp_sfv_size];
		score_t score_cal[potential_size];
		int score_mark[potential_size];
		score_t tmp;
		zSFVec **tmp_huge_sfv;
		tmp_huge_sfv = zCalloc(large_top,sizeof(zSFVec*),"zRunViterbiAndForward sfv"); // waste storage soon correct
		size = sizeof(zSFVec);
		for(g=0;g<large_top;g++){
			tmp_huge_sfv[g]=zMalloc(size,"zMalloc tmp_huge_sfv[i]");
			zInitSFVec(tmp_huge_sfv[g],10);
		}


		for (i = 0; i < potential_size; ++i)
		{
			score_mark[i]=i;
		}
	//	printf("oops\n");
		zForTmpSFVScore(trellis,huge_sfv,tmp_sfv,score_box);

	//	printf("what\n");
		zModifyScore(trellis,huge_sfv,tmp_sfv,score_cal,score_box);
	//	printf("the\n");
	//	zSortScore(score_cal,score_mark,LARGE_TOP*top_num);
		zSortScore(score_cal,score_mark,potential_size);
	//	printf("Fuck\n");
		for (k = 0; k < large_top; ++k)
		{
			score[k] = score_cal[k];
			i = score_mark[k]/tmp_sfv_size;
			j = score_mark[k]%tmp_sfv_size;
	//		printf("index: %d\n",k );
	//		printf("huge_sfv[%d] and tmp_sfv[%d]\n",i,j);
			zMergeSFVec(huge_sfv[i],tmp_huge_sfv[k]);
			zMergeSFVec(tmp_sfv[j],tmp_huge_sfv[k]);
	//		printf("you\n")

		}
	//	printf("doing\n");
	//	huge_sfv=tmp_huge_sfv;
		for (i = 0; i < large_top; ++i)
		{
			zFreeSFVec(huge_sfv[i]);
			zFree(huge_sfv[i]);
			huge_sfv[i]=tmp_huge_sfv[i];
		}
	//	printf("now\n");

		
	}




///////////////////////////////////////

	return;	
	

//	printf("merge fin\n");
}


void zTraceTrellisSort(zSFVec **huge_sfv,zSFVec **sfv,score_t* score)
{
	int i,j;
	score_t tmp_score;
	int mark[LARGE_TOP];
	int tmp_mark;
	for (i = 0; i < LARGE_TOP; ++i)
	{
		mark[i]=i;
	}

	for (i = 0; i < LARGE_TOP; ++i)
	{
		printf("%f\n",score[i]);
	}
	//bubble sort
	for (j = 0; j < LARGE_TOP-1; ++j)
	{
		for (i = 0; i < LARGE_TOP-1; ++i)
		{
			if (score[i]<score[i+1])
			{
				tmp_score=score[i];
				score[i]=score[i+1];
				score[i+1]=tmp_score;

				tmp_mark=mark[i];
				mark[i]=mark[i+1];
				mark[i+1]=tmp_mark;
			}
		}
	}

printf("________________________________________________\n");
	for (i = 0; i < LARGE_TOP; ++i)
	{
		printf("%d-%d\n",i,mark[i] );
		sfv[i]=huge_sfv[mark[i]];
	}	
printf("________________________________________________\n");
	for (i = 0; i < LARGE_TOP; ++i)
	{
		printf("%f\n",score[i]);
	}

} 

void zModifySortScoreForSponge(zTrellis *trellis,zSFVec **huge_sfv,score_t* path_score,int large_top,int type, int pre_state,int mode)
{	

	printf("MSS_Pre_state %d,%d\n",pre_state,mode);
	int i,j;
	zHMM *hmm=trellis->hmm;
	score_t t_score=0.0;	
	int limit = zCountRealSFV(huge_sfv,large_top);

	int head=-1;
	int tail=-1;

	for (i = 0; i < limit; ++i)
	{
		head = -1;
		tail = -1;
		if( ( huge_sfv[i]->elem[0].state==0 || (huge_sfv[i]->elem[0].state>=43 && huge_sfv[i]->elem[0].state<=46) || (huge_sfv[i]->elem[0].state>=13 && huge_sfv[i]->elem[0].state<=16) ) )
		{
			if (mode==1||mode==4){t_score=0;}
			else if(mode==2||mode==5)
			{
				head=0;
				t_score=0;
			}
			else
			{
				if (pre_state == -1)
				{
					if(hmm->state[huge_sfv[i]->elem[0].state].type == GINTERNAL)
					{
						t_score=huge_sfv[i]->elem[0].score
							+ zGetInitProb(trellis->hmm, huge_sfv[i]->elem[0].state, trellis->iiso_group)
							+ zGetIntergenicContinueScore(trellis);
					}
					else if (hmm->state[huge_sfv[i]->elem[0].state].type == INTERNAL)
					{
						t_score=huge_sfv[i]->elem[0].score
							+ zGetInitProb(trellis->hmm, huge_sfv[i]->elem[0].state, trellis->iiso_group)
							+ zGetTransitionScore(trellis->hmm,huge_sfv[i]->elem[0].state,huge_sfv[i]->elem[0].state,trellis->tiso_group);
					}
					else{t_score=MIN_SCORE;printf("ForSponge Error\n");}
				}
				else{			
					t_score = zGetTransitionScore(trellis->hmm,pre_state,huge_sfv[i]->elem[0].state,trellis->tiso_group)
						+ huge_sfv[i]->elem[0].score;
				}
			}
		}
		else{
			head=1;
			if (pre_state == -1)
			{
				if(hmm->state[huge_sfv[i]->elem[0].state].type == GINTERNAL)
				{
					t_score=huge_sfv[i]->elem[0].score
						+ zGetInitProb(trellis->hmm, huge_sfv[i]->elem[0].state, trellis->iiso_group)
						+ zGetIntergenicContinueScore(trellis);
				}
				else if (hmm->state[huge_sfv[i]->elem[0].state].type == INTERNAL)
				{
					t_score=huge_sfv[i]->elem[0].score
						+ zGetInitProb(trellis->hmm, huge_sfv[i]->elem[0].state, trellis->iiso_group)
						+ zGetTransitionScore(trellis->hmm,huge_sfv[i]->elem[0].state,huge_sfv[i]->elem[0].state,trellis->tiso_group);
				}
				else{t_score=MIN_SCORE;printf("ForSponge Error\n");}
			}
			else{			
				t_score = zGetTransitionScore(trellis->hmm,pre_state,huge_sfv[i]->elem[0].state,trellis->tiso_group)
					+ huge_sfv[i]->elem[0].score;
			}
		}

		for (j = 1; j < huge_sfv[i]->size; ++j)
		{
			if( ( huge_sfv[i]->elem[j].state==0 || (huge_sfv[i]->elem[j].state>=43 && huge_sfv[i]->elem[j].state<=46) || (huge_sfv[i]->elem[j].state>=13 && huge_sfv[i]->elem[j].state<=16) ) )
			{
				if (mode==1||mode==4){printf("contnue14\n");continue;}
				else if(mode==2||mode==5)
				{
					if (head==0){printf("contnue25\n");continue;}
					else if(head==1)
					{
						t_score += zGetTransitionScore(trellis->hmm,huge_sfv[i]->elem[j-1].state,huge_sfv[i]->elem[j].state,trellis->tiso_group)
							+ huge_sfv[i]->elem[j].score;
					}

				}
				else{
					t_score += zGetTransitionScore(trellis->hmm,huge_sfv[i]->elem[j-1].state,huge_sfv[i]->elem[j].state,trellis->tiso_group)
						+ huge_sfv[i]->elem[j].score;
				}
			}
			else{
				head=1;
				t_score += zGetTransitionScore(trellis->hmm,huge_sfv[i]->elem[j-1].state,huge_sfv[i]->elem[j].state,trellis->tiso_group)
					+ huge_sfv[i]->elem[j].score;
			}
		}

		if (mode==2||mode==5)
		{
			for (j = huge_sfv[i]->size-1; j > 0; --j)
			{
				if( ( huge_sfv[i]->elem[j].state==0 || (huge_sfv[i]->elem[j].state>=43 && huge_sfv[i]->elem[j].state<=46) || (huge_sfv[i]->elem[j].state>=13 && huge_sfv[i]->elem[j].state<=16) ) )
				{
					tail=0;
					t_score -= zGetTransitionScore(trellis->hmm,huge_sfv[i]->elem[j-1].state,huge_sfv[i]->elem[j].state,trellis->tiso_group)
						+ huge_sfv[i]->elem[j].score;
					printf("contnueTail\n");
				}
				else{
					tail=1;
					break;
				}

			}
		}
		path_score[i] = t_score;
	}
	zSFVec *index;
	score_t medium;
	if (mode==0 || mode==1 || mode==2)
	{
		printf("first to last: %f~%f\n",path_score[0],path_score[limit-1]);
		zQuickSortSFVec(path_score, 0, limit-1,huge_sfv);
	//	zSortScoreSFVec(path_score, huge_sfv, limit);
		// for (i = 0; i < limit; ++i)
		// {
		// 	for (j = i+1; j < limit; ++j)
		// 	{
		// 		if (path_score[i]<path_score[j])
		// 		{
		// 			medium=path_score[i];
		// 			path_score[i]=path_score[j];
		// 			path_score[j]=medium;
		// 			index=huge_sfv[i];
		// 			huge_sfv[i]=huge_sfv[j];
		// 			huge_sfv[j]=index;
		// 		}
		// 	}
		// }
		printf("after first to last: %f~%f\n",path_score[0],path_score[limit-1]);
	}

//	printf("BUBBLE FIN\n");

}

void zModifySortScore(zTrellis *trellis,zSFVec **huge_sfv,score_t* path_score,int large_top)
{	
	int i,j;
	zHMM *hmm=trellis->hmm;
	score_t t_score=0.0;	
	for (i = 0; i < large_top; ++i)
	{
		if(hmm->state[huge_sfv[i]->elem[0].state].type == GINTERNAL)
		{
			t_score=huge_sfv[i]->elem[0].score
				+ zGetInitProb(trellis->hmm, huge_sfv[i]->elem[0].state, trellis->iiso_group)
				+ zGetIntergenicContinueScore(trellis);

	//	printf("T_SCORE %d,%f\t",i,t_score );
	//	printf("  %d,%s,%f\n",i, zStrIdx2Char(hmm->state[huge_sfv[i]->elem[0].state].name),huge_sfv[i]->elem[0].score );
		}
		else if (hmm->state[huge_sfv[i]->elem[0].state].type == INTERNAL)
		{
			t_score=huge_sfv[i]->elem[0].score
				+ zGetInitProb(trellis->hmm, huge_sfv[i]->elem[0].state, trellis->iiso_group)
				+ zGetTransitionScore(trellis->hmm,huge_sfv[i]->elem[0].state,huge_sfv[i]->elem[0].state,trellis->tiso_group);
	//	printf("T_SCORE %d,%f\t",i,t_score );
	//	printf("  %d,%s,%f\n",i, zStrIdx2Char(hmm->state[huge_sfv[i]->elem[0].state].name),huge_sfv[i]->elem[0].score );
		}
		else{
	//		printf("ERROR_SCORE %d,%s,%f\n",i, zStrIdx2Char(hmm->state[huge_sfv[i]->elem[0].state].name),huge_sfv[i]->elem[0].score );
			t_score = MIN_SCORE;
		//	zDie("Modifiying score error");
		}
	//	t_score=huge_sfv[i]->elem[0].score+ zGetInitProb(trellis->hmm, huge_sfv[i]->elem[0].state, trellis->iiso_group)+zGetTransitionScore(trellis->hmm,huge_sfv[i]->elem[0].state,huge_sfv[i]->elem[0].state,trellis->tiso_group);
	//	printf("PECK %d: %f,%f,%f\n",i+1,t_score,zGetInitProb(trellis->hmm, huge_sfv[i]->elem[0].state, trellis->iiso_group),zGetIntergenicContinueScore(trellis));//zGetTransitionScore(trellis->hmm,huge_sfv[i]->elem[0].state,huge_sfv[i]->elem[0].state,trellis->tiso_group));
		for (j = 1; j < huge_sfv[i]->size; ++j)
		{
			t_score += zGetTransitionScore(trellis->hmm,huge_sfv[i]->elem[j-1].state,huge_sfv[i]->elem[j].state,trellis->tiso_group)
				+ huge_sfv[i]->elem[j].score;
	//		printf("%d,%s,%d,%s\n", huge_sfv[i]->elem[j-1].state,zStrIdx2Char(hmm->state[huge_sfv[i]->elem[j-1].state].name),huge_sfv[i]->elem[j].state,zStrIdx2Char(hmm->state[huge_sfv[i]->elem[j].state].name));
	//		printf("t_SCORE %d,%d,%f,%f\n",i,j,zGetTransitionScore(trellis->hmm,huge_sfv[i]->elem[j-1].state,huge_sfv[i]->elem[j].state,trellis->tiso_group),huge_sfv[i]->elem[j].score );
		}
		path_score[i] = t_score;
	//	printf("T_SCORE %d,%f\n",i,t_score );
	}
	zSFVec *index;
	score_t medium;
	for (i = 0; i < large_top; ++i)
	{
		for (j = i+1; j < large_top; ++j)
		{
			if (path_score[i]<path_score[j])
			{
				medium=path_score[i];
				path_score[i]=path_score[j];
				path_score[j]=medium;
				index=huge_sfv[i];
				huge_sfv[i]=huge_sfv[j];
				huge_sfv[j]=index;
			}
		}
	}
	printf("BUBBLE FIN\n");

}
void zTrimTmpSFV(zSFVec** tmp_sfv,coor_t tmp_cpoint)
{
	int i,g;
	int trim;
	for (i = 0; i < top_num; ++i)
	{
		for (g = 0; g < tmp_sfv[i]->size; ++g)
		{
			if (tmp_sfv[i]->elem[g].end==tmp_cpoint)
			{
				printf("OMG\n");
				break;
			}
		}
		trim = g;
		printf("hello %d\n" ,trim);
		zTrimSFVec(tmp_sfv[i], trim);	
		printf("bye\n");
	}
}

void zShowRootChildren(zTBTree *tree)
{

	zTBTreeNode *p;
	p=tree->root->child;
		printf("ROOT %d\n",tree->root->pos);
	printf("%d: %d,%d\n",p->pos,p->state,p->story );
	while(p->rsib!=NULL)
	{
		p=p->rsib;
		printf("%d: %d,%d\n",p->pos,p->state,p->story );

	}


}
int zCountRootRealChildren(zTBTree *tree)
{
	int num=tree->root->children;
	int lockChild=0;
	zTBTreeNode *c;
	c=tree->root->child;
	while(c!=NULL)
	{
//		printf("ininin %d,%d,%d,%d\n",c->state,c->pos,c->length,c->children);
		if (c->children==0)
		{
		//	printf("ininin %d,%d,%d,%d\n",c->state,c->pos,c->length,c->children);
			lockChild++;
			if (c->child!=NULL)
			{
				printf("HUGE ERROR\n");
			}
		}
		if (c->lock==1)
		{
//			printf("ininIN\n");
			lockChild++;
		}
		c=c->rsib;
	}
	return num-lockChild;
}
void zTreeterbiMarkBranch(zTrellis  *trellis, coor_t pos,zTBTree *tree, zTBTreeNode *n)
{
	if (n->pos > pos)
	{
		return;
	}
	if (n->children<=0)
	{
//		printf("Mark Branch Collect %d,%d,%d\n",n->pos,n->state,n->story);
	//	zTreeterbiCollectNew(trellis,n->pos,n->state,n->story,300,0);
		n->lock=-1;
		return;
	}
	zTBTreeNode *c;
	c=n->child;
	while(c!=NULL)
	{	
		c->lock=0;
		zTreeterbiMarkBranch(trellis,pos,tree,c);
		c=c->rsib;
	}



}
void zTreeterbiReverseCPoint(zTrellis  *trellis, coor_t pos,zTBTree *tree)
{	
	if (pos<=tree->root->pos)
	{
		return;
	}

	zTBTreeNode *c;
	int num=tree->root->children;
	printf("RE: %d\n",num );
	c=tree->root->child;
	while(c!=NULL)
	{
//		printf("rerere %d,%d,%d,%d\n",c->state,c->pos,c->length,c->children);
//		if (c->children==0)
//		{	printf("0rerere %d,%d,%d,%d\n",c->state,c->pos,c->length,c->children);}
		c->lock=0;
		zTreeterbiMarkBranch(trellis,pos,tree,c);
		c=c->rsib;
	}
}

void zTreeterbiClean(zTrellis *trellis, coor_t pos,zTBTree *tree)
{
	
	
	if (pos<=tree->root->pos)
	{
		return;
	}



	int i,j,k;
	for (i = tree->root->pos; i < pos; ++i)
	{
		if (trellis->cell[i]==NULL){continue;}
		for (j = 0; j<49; ++j)
		{
			if (trellis->cell[i][j].layer==NULL){continue;}
			for (k = 0; k< top_num; ++k)
			{
				if (trellis->cell[i][j].layer[k]==NULL){continue;}
				if (trellis->cell[i][j].layer[k]->lock==-1)
				{
	//				printf("clean %d,%d,%d,%d,%d\n",i,j,k,trellis->cell[i][j].layer[k]->length,trellis->cell[i][j].layer[k]->children);
					zTreeterbiCollectNew(trellis,i,j,k,300,0);
				}
			}
		}
	}	
}

void zFindCPointSingle(zTrellis *trellis,coor_t pos,zTBTree *tree,zTBTreeNode *n,int depth)
{
	int i;
	zTBTreeNode *c;
	c=n->child;
	while(c!=NULL &&  c->pos < pos)
	{
		for (i = 0; i < depth; ++i)
		{
			printf(" ");
		}
		printf("CP %d,%d,%d,%d,%d\n",c->pos,c->state,c->story,c->length,c->children);
		zFindCPointSingle(trellis,pos,tree,c,depth+1);
		c=c->rsib;
	}
}

int zExpandMulti(zTBTreeNode **multi,int maxLoc,int multi_limit,coor_t pos,int diff)
{
	int i;
	int maxPos=multi[maxLoc]->pos;
	zTBTreeNode *c;
	int ori_multi_limit=multi_limit;
	int count=0;
	for (i = 0; i < ori_multi_limit; ++i)
	{

		if (diff==0)
		{
			c=multi[i]->child;
			if (c==NULL)
			{
				printf("zExpandMulti diff0 first kid ERROR\n");
				return -1;
			//	zDie("zExpandMulti diff0 first kid ERROR\n");
			}
			count=0;
			while(c!=NULL)
			{
				// if (c->pos>=pos)
				// {
				// 	c=c->rsib;
				// 	continue;
				// }
				if (c->pos>pos){return -1;}



				if (count==0)
				{
					multi[i]=c;
				}
				else{
					multi[multi_limit]=c;
					multi_limit++;
				}
				count++;
				c=c->rsib;
			}
			if (count==0)
			{
			//	zDie("zExpandMulti diff0 first kid ERROR\n");
				printf("diff0 first kid ERROR %d   %d,%d,%d\n",multi_limit,multi[i]->pos,multi[i]->state,multi[i]->story );
				return -1;
			}

		}
		else if (multi[i]->pos<maxPos)
		{
			c=multi[i]->child;
			if (c==NULL)
			{
				zDie("zExpandMulti diff1 first kid ERROR\n");
			}
			count=0;
			while(c!=NULL)
			{
				// if (c->pos>=pos)
				// {
				// 	c=c->rsib;
				// 	continue;
				// }
				if (c->pos>pos){return -1;}

				if (count==0)
				{
					multi[i]=c;
				}
				else{
					multi[multi_limit]=c;
					multi_limit++;
				}
				count++;
				c=c->rsib;
			}
			if (count==0)
			{
			//	zDie("zExpandMulti diff1 first kid ERROR\n");
				printf("diff1 first kid ERROR %d   %d,%d,%d\n",multi_limit,multi[i]->pos,multi[i]->state,multi[i]->story );
				return -1;
			}

		}
	}
	return multi_limit;
}



coor_t zFindCPointNew(zTrellis *trellis,coor_t pos,zTBTree *tree,coor_t tmp_cpoint)
{
	zTBTreeNode *c;
	zTBTreeNode *multi[49*top_num];
	int multi_limit=0;
	
	// c=tree->root->child;
	// if(c==NULL){printf("FindCP NULL\n");}
	// while(c!=NULL)
	// {printf("C_POS %d\n",c->pos);c=c->rsib;}

	c=tree->root->child;

	if (c->children>49*top_num)
	{
		printf("zFindCPointNew out of limit\n");
		return -1;
	//	zDie("zFindCPointNew out of limit");
	}


	if (c->pos>pos)
	{
		return -1;
	}





	while(c!=NULL)
	{
		// printf("CP %d,%d,%d,%d,%d\n",c->pos,c->state,c->story,c->length,c->children);
		// zFindCPointSingle(trellis,pos,tree,c,1);
		// if (c->pos>=pos)
		// {
		// 	c=c->rsib;
		// 	continue;
		// }

		multi[multi_limit]=c;
		multi_limit++;
		// printf("m %d,%d,%d,%d,%d\n",multi[multi_limit]->pos,multi[multi_limit]->state,multi[multi_limit]->story,multi[multi_limit]->length,multi[multi_limit]->children);
		// printf("c %d,%d,%d,%d,%d\n",c->pos,c->state,c->story,c->length,c->children);
		c=c->rsib;
		// printf("m %d,%d,%d,%d,%d\n",multi[multi_limit]->pos,multi[multi_limit]->state,multi[multi_limit]->story,multi[multi_limit]->length,multi[multi_limit]->children);
		// printf("c %d,%d,%d,%d,%d\n",c->pos,c->state,c->story,c->length,c->children);
	}
	if (multi_limit==0)
	{
		return -1;
	}

	int diff=1;
	int tmp_state=-1;
	int state_diff=0;
	int i;
	int maxPos=0;
	int maxLoc=0;
	int finalMaxPos=-1;
	while(true)
	{
		maxPos=0;
		maxLoc=0;
		diff=0;
		state_diff=0;
		tmp_state=-1;
//		printf("Multi %d\n", multi_limit );
		for (i = 0; i < multi_limit; ++i)
		{
			// if (i<100)
			// {
			// 	printf("%d:%d,%d\n",i,multi[i]->pos,multi[i]->state );
			// }
			if (multi[i]->state!=tmp_state)
			{
				if (tmp_state!=-1)
				{
					state_diff=1;	
				}
				tmp_state=multi[i]->state;
			}
	//		printf("XXX\n");
			if (multi[i]->pos>maxPos)
			{
	//			printf("Xin %d - %d - %d - %d\n",pos,maxPos,multi[i]->pos,finalMaxPos );
				if (maxPos!=0)
				{diff=1;}
				maxPos=multi[i]->pos;
				maxLoc=i;
				if (maxPos>pos)
				{
		//			printf("ininin %d\n",finalMaxPos);
					return finalMaxPos;
				}
			}
			else if (multi[i]->pos<maxPos)
			{
		//		printf("Yin %d - %d - %d - %d\n",pos,maxPos,multi[i]->pos,finalMaxPos );
				if (maxPos!=0)
				{diff=1;}
			}
		//	printf("YYY\n");
		}
	//	printf("\n");
		if(diff==0)
		{
			if (state_diff==0)
			{
				finalMaxPos=maxPos;
			//	new added   annotated->last_cp used->first_cp
				if (finalMaxPos!=tmp_cpoint)
				{
					if (tmp_state==0 || tmp_state==23 || tmp_state==30 || (tmp_state>=43 && tmp_state<=48) || (tmp_state>=13 && tmp_state<=16) )
					{
						return finalMaxPos;
					}
				//	return finalMaxPos;
				}
			}
			
			printf("NFinalMaxPos-%d: %d\n",state_diff,finalMaxPos);
		}
		multi_limit=zExpandMulti(multi,maxLoc,multi_limit,pos,diff);
		if (multi_limit==-1)
		{
			return finalMaxPos;
		}
	//	printf("Expand fin\n");
	}
	return -1;
}

int zWhetherInter(int state)
{
	int filter[13]={0,13,14,15,16,23,30,43,44,45,46,47,48};
	int i;
	for (i = 0; i < 13; ++i)
	{
		if (state==filter[i])
		{
			return 1;
		}
	}
	return 0;
}

int zWhetherCrossInter(zSFVec* sfv)
{
	int i;
	int mark = -1;//0:all ex 1:all inter 2:ex-in 3:in-ex 4:whole 5:cross
	for (i = 0; i < sfv->size; ++i)
	{
		if(sfv->elem[i].state==0 || (sfv->elem[i].state>=43 && sfv->elem[i].state<=46) || (sfv->elem[i].state>=13 && sfv->elem[i].state<=16))
		{
			if (mark==-1){mark=1;}
			else if(mark==1){mark=1;}
			else if(mark==0){mark=2;}
			else if(mark==2){mark=2;}
			else{printf("Whole ERROR\n");mark=4;}
		}
		else{
			if (mark==-1){mark=0;}
			else if(mark==1){mark=3;}
			else if(mark==0){mark=0;}
			else if(mark==3){mark=3;}
			else{printf("Cross ERROR\n");mark=5;}

		}
	}
	return mark;
}

int zWhetherNodeCrossInter(zSFVSpongeNode *spn)
{
	int i;
	int mark=-1;
	for (i = 0; i < spn->size; ++i)
	{
		if(zWhetherCrossInter(spn->node[i])==5){mark=1;}
		else{mark=0;return 0;}
	}
	return mark;
}

// void zShrinkSFVSponge(zSFVSponge* sp,zSFVSpongeNode* spn)
// {
// 	int limit=spn->size;
// 	coor_t coord[limit];
	
// 	sfv->elem[i].state==0 || (sfv->elem[i].state>=43 && sfv->elem[i].state<=46)

// }




void zShowSFVSponge(zSFVSponge *sp,zTrellis *trellis)
{
	printf("SPN show\n");
	zSFVSpongeNode* spn;
	int i;
	int inter_mark=0;
	int mark=0;
	int kinds = 1;
	int new_kinds=0;
	int gene_start = 50;
	int gene_end = 0;
	zSFVec* gene_sfv;
	gene_sfv=zMalloc(sizeof(zSFVec),"zMalloc gene_sfv");
	zInitSFVec(gene_sfv,10);



	spn = zSFVSpongeMoveFirst(sp);
	while(spn!=NULL)
	{
		inter_mark = zWhetherInter(spn->end_state);

		mark = zWhetherCrossInter(spn->node[0]);


		kinds*=spn->size;

		printf("%d,%d-%d:%s,%d\nMarkCross ",spn->size,spn->start,spn->end,zStrIdx2Char(trellis->hmm->state[spn->end_state].name),inter_mark);


		// for (i = 0; i < spn->node[0]->size; ++i)
		// {
		// 	printf("\t %s\n",zStrIdx2Char(trellis->hmm->state[spn->node[0]->elem[i].state].name) );
		// }

		zMergeSFVec(spn->node[0],gene_sfv);

		if (inter_mark==1)
		{
		//	gene_end = spn->end;
			printf("Kinds %d,%d\n",kinds,zWhetherCrossInter(gene_sfv));
			printf("Gene: %d-%d\n",gene_sfv->elem[0].start,gene_sfv->elem[gene_sfv->size-1].end );
			for (i = 0; i < gene_sfv->size; ++i)
			{
				printf("\t %s\n",zStrIdx2Char(trellis->hmm->state[gene_sfv->elem[i].state].name) );
			}

			zResetSFVec(gene_sfv);
			inter_mark = 0;
			kinds = 1;
		//	gene_start = gene_end+1;
		}
		spn = zSFVSpongeMoveNext(sp);
	}
	zFreeSFVec(gene_sfv);
	zFree(gene_sfv);
	gene_sfv=NULL;



}
void zScoreSFVSponge(zSFVSponge *sp,zTrellis *trellis)
{
	printf("SPN score\n");
	zSFVSpongeNode* spn;
	int i,g,k;
	int inter_mark=0;
//	int mark=0;
	int kinds = 1;
	int gene_start = 50;
	int gene_end = 0;
	int count;
	int type;
	zSFVSpongeNode* spn_index=sp->head->next;
	zGTFVec gtfvec;

	int small_top=50;

	int mark[top_num];
	for (i = 0; i < top_num; ++i)
	{
		mark[i]=-9;
	}

//printf("AAA1\n");
	zSFVec** large_sfv;
	large_sfv=zCalloc(LARGE_TOP,sizeof(zSFVec*),"zMalloc large_sfv");
	int data_size=sizeof(zSFVec);
	for (g = 0; g < LARGE_TOP; ++g)
	{
		large_sfv[g]=zMalloc(data_size,"zMalloc large_sfv");
		zInitSFVec(large_sfv[g],10);
	}
	
//printf("AAA2\n");	

	zSFVec* gene_sfv;
	gene_sfv=zMalloc(sizeof(zSFVec),"zMalloc gene_sfv");
	zInitSFVec(gene_sfv,10);
//printf("AAA3\n");	
	spn = zSFVSpongeMoveFirst(sp);
//printf("AAA4\n");
	while(spn!=NULL)
	{
//printf("AAA4.2\n");

//printf("A1 %d\n",spn->location );
//printf("A2 %d\n",spn->pre_state );
		inter_mark = zWhetherInter(spn->end_state);
//		mark = zWhetherCrossInter(spn->node[0]);

//printf("AAA4.5 %d \n",inter_mark);

		for (i = 0; i < spn->size; ++i)
		{
			mark[i] = zWhetherCrossInter(spn->node[i]);
		// printf("AAA4.8  %d\n",i);
		// 	if (mark[i]==5)
		// 	{
		// 		printf("MARK_FIVE \n");
		// 		for (g = 0; g < spn->node[i]->size; ++g)
		// 		{
		// 			printf("\t\t%d-%d,%s\n",spn->node[i]->elem[g].start,spn->node[i]->elem[g].end,zStrIdx2Char(trellis->hmm->state[spn->node[i]->elem[g].state].name));
		// 		}
		// 	}


		}
//printf("AAA5\n");
		for (; i < top_num; ++i)
		{
			mark[i]=-9;
		}



//intf("AAA6\n");

		kinds*=spn->size;

		printf("%d,%d-%d:%s,%d\nMark",spn->size,spn->start,spn->end,zStrIdx2Char(trellis->hmm->state[spn->end_state].name),inter_mark );

		for (i = 0; i < top_num; ++i)
		{
			if (mark[i]==-9){break;}
			printf("%d,",mark[i] );
		}
		printf("\n");
//printf("AAA7\n");

		for (i = 0; i < spn->size; ++i)
		{
			printf("\t%f\n",zScoreSFVec(trellis,spn->node[i]) );
			// for (g = 0; g < spn->node[i]->size; ++g)
			// {
			// 	printf("\t\t%d-%d,%s\n",spn->node[i]->elem[g].start,spn->node[i]->elem[g].end,zStrIdx2Char(trellis->hmm->state[spn->node[i]->elem[g].state].name));
			// }
		}
		printf("xxxxxxxxxx\n");
		zMergeSFVec(spn->node[0],gene_sfv);

		if (inter_mark==1 || spn->next==sp->tail)
		{
		//	gene_end = spn->end;
			type = zWhetherCrossInter(gene_sfv);
			printf("Kinds %d,%d\n",kinds,type);
			printf("Gene: %d-%d\n",gene_sfv->elem[0].start,gene_sfv->elem[gene_sfv->size-1].end );

			score_t score_list[LARGE_TOP];
			for (g = 0; g < LARGE_TOP; ++g)
			{
				score_list[g]=MIN_SCORE;
			}

			zTracePartialTrellisMergeForSponge(trellis,large_sfv,spn_index,spn,LARGE_TOP,score_list,type);
			printf("ssssss %d-%d,%f\n",spn_index->location,spn->location,score_list[0] );
			count=zCountRealSFV(large_sfv,LARGE_TOP);
			printf("count %d\n",count);





			for (g = 0; g < count; ++g)
			{
				printf("%f\n",score_list[g] );
			}
			printf("allscore\n");
			printf("toGTF\n");

			


			for (g = 0; g < count; ++g){
				printf("gg %d\n", g);
			//	if (g<small_top)
			//	{	
					if (zOption("i"))
					{
						zWriteSFVec(stdout,large_sfv[g]);
					}
					else{
						zInitGTFVec(&gtfvec, large_sfv[g]->size);
					//	zSFVec2GTFVec(trellis->hmm, large_sfv, &gtfvec, filename, t_num);
						zSFVec2GTFVec(trellis->hmm, large_sfv[g], &gtfvec, "", 0);
						zSortGTFVec(&gtfvec);
						zWriteGTFVec(stdout, &gtfvec);
						zFreeGTFVec(&gtfvec);
					}
			//	}
				zResetSFVec(large_sfv[g]);
			}

			for (i = 0; i < gene_sfv->size; ++i)
			{
				printf("\t %s\n",zStrIdx2Char(trellis->hmm->state[gene_sfv->elem[i].state].name) );
			}

			zResetSFVec(gene_sfv);
			inter_mark = 0;
			kinds = 1;
			spn_index=spn->next;
		//	gene_start = gene_end+1;
		}
		spn = zSFVSpongeMoveNext(sp);
	}
printf("AAA7\n");
	zFreeSFVec(gene_sfv);
	zFree(gene_sfv);
	gene_sfv=NULL;

	for(g=0;g<LARGE_TOP;g++)
	{
		zFreeSFVec(large_sfv[g]);
		zFree(large_sfv[g]);
		large_sfv[g]=NULL;
	}
	zFree(large_sfv);




}



zSFVec** zRunViterbiAndForward (zTrellis *trellis, score_t* path_score, zSFVSponge *sponge) {
	zHMM         *hmm = trellis->hmm;
	zDNA         *dna = trellis->dna;

	coor_t        i,l,tmp_start;
	int           j;        /* iterator for internal states */
	int           k;        /* iterator for previous states */
	int           m;
	/*	int           m;    */    /*top mode*/
	int           g;
	int           *nrange, iso_group, u;
	zDurationGroup *group;
	zIVec*        jumps;
	zSFVec        **sfv;

	zSFList       *sfl;
	zSfeature     *f;
	size_t        size;
	//Added by Jim
	int window;
	zTBTree *tree;
	tree=zMalloc(sizeof(zTBTree),"zRunViterbiAndForward: tree");
	zInitTBTree(tree);
	tree->root->pos = PADDING-1;
	tree->root->score = MIN_SCORE;
	zSFVec  **tmp_sfv;
	zSFVec **huge_sfv;
	score_t tmp_path_score[LARGE_TOP];
	coor_t tmp_cpoint=0;
	coor_t old_cpoint=0;
	coor_t older_cpoint=0;
	int show_mark=0;
	zSfeature *tmp_istate;
	int large_top;
	zTrellisCell cell_temp;
	int tmp_cpoint_state;
	int zz=0;
	int tmp_state_length=0;
	int trim=0;
	int cpoint_check_skip=0;
	
	int zzz=0;
	int skip=0;

//	zSFVSponge *sponge;

	if (zOption("a"))
	{
		sponge = zMalloc(sizeof(zSFVSponge),"zRunViterbiAndForward: sponge");
		zInitSFVSponge(sponge);
	}


	// for (j = 0; j < hmm->states; ++j)
	// {
	// 	printf("%d:%s  %d\n",j,zStrIdx2Char(trellis->hmm->state[j].name),trellis->hmm->state[j].type);
	// }




	/*-----------------------------------------------
	  |            [0] [1] [2] [3] [4] [5] [6] [7] ...
	  |(None)  [0]
	  |(Misc)  [1]              X
	  |(Inter) [2]
	  |:                  X = trellis->cell[1][3] 
	  |:                               cell[j][i]
	  ------------------------------------------------*/
	 	 
	/* 	Viterbi and Forward Alg initialization */
	zTrace2("initializing Viterbi vars\n");

	zAllocViterbiVarsInitial(trellis);

	sfl = zMalloc(sizeof(zSFList),"zRunViterbiAndForward: sfl");
	zInitSFList(sfl);

	/* EVAN what is this nrange crap?*/
	nrange = zMalloc((hmm->states*sizeof(int)), "zRunViterbiAndForward nrange");
	/////////printf states
	///////printf("------------------states number\n");
	///////////printf("%d\n",hmm->states);
	for (j = 0; j < hmm->states; j++){
		group = hmm->dmap[hmm->state[j].duration];
		iso_group = zGetDurationIsochoreGroup(group, zGetDNAGC(dna));

		if ((hmm->state[j].type == EXPLICIT) 
			|| (hmm->state[j].type == INTERNAL)){
			nrange[j] = group->duration[iso_group].distributions;
		}

	}

	/* EVAN there must be a better way to do this explicit stuff */

 	zTrace2("starting initilization\n");
 	//Added by Jim
//	large_top=int(pow(float(top_num),float(show_mark)));

	large_top=LARGE_TOP; //FINDDDD
	huge_sfv = zCalloc(large_top,sizeof(zSFVec*),"zRunViterbiAndForward sfv");
	if (!zOption("a"))
	{
		size = sizeof(zSFVec);
		for(g=0;g<large_top;g++){
			huge_sfv[g]=zMalloc(size,"zMalloc huge_sfv[i]");
			zInitSFVec(huge_sfv[g],10);
		}
	}


	tmp_sfv = zCalloc(top_num,sizeof(zSFVec*),"zRunViterbiAndForward sfv");
			printf("VVV2\n");


			///////////exchange i and j loop
		  	//////////for (j = 0; j < hmm->states; j++) {
			for (i = 0; i < PADDING; i++) {
				for (j = 0; j < hmm->states; j++) {
			zAllocViterbiVarsSingle(trellis,i,j);
	//		printf("VVV1\n");
			for(m=0;m<top_num;m++){
				trellis->cell[i][j].layer[m]->left_pos = 0;
				if(m==0){
					trellis->cell[i][j].layer[m]->score = zGetInitProb(trellis->hmm, j, trellis->iiso_group);
		//			printf("JIMJIM %f,%s\n",trellis->cell[i][j].layer[m]->score,zStrIdx2Char(trellis->hmm->state[j].name));
					trellis->cell[i][j].layer[m]->trace = j;
					trellis->cell[i][j].layer[m]->length = 0;
					trellis->cell[i][j].layer[m]->trace_top = 0;
				}
				else{
					trellis->cell[i][j].layer[m]->score = MIN_SCORE;
					trellis->cell[i][j].layer[m]->trace = -1;
					trellis->cell[i][j].layer[m]->length = -1;
					trellis->cell[i][j].layer[m]->trace_top = -1;
				}
				//Added by Jim
			//	trellis->cell[j][i].layer[m]->pos=i;
			}

			/* Allows introns to extend to the start of the */
			/* sequence when modeled as an EXPLICIT state   */
			//if ((trellis->cell[j][i].score[0] > MIN_SCORE)

			if ((trellis->cell[i][j].layer[0]->score > MIN_SCORE)
				&& (hmm->state[j].type == EXTERNAL)
				&& ((i+1) == PADDING)){ /*JIM logically here is impossible to go inside*/
				zPushIVec(&trellis->extpos[j], i);
			}
			
			/////kevin:at layer 0, MIN_INIT_SCORE is -10000, bigger than other's layers init score MIN_SCORE
			if (trellis->cell[i][j].layer[0]->score == MIN_SCORE) {
				trellis->cell[i][j].layer[0]->score = MIN_INIT_SCORE;
			}
			
			/*////////////////kevin
			if (((i+1) == PADDING)
				&& ((hmm->state[j].type == EXPLICIT)
					|| (hmm->state[j].type == INTERNAL))){
				trellis->cell[i][j].submax = zMalloc(sizeof(zSFVec),
													 "zRunViterbiAndForward submax");
				zInitSFVec(trellis->cell[i][j].submax, nrange[j]);
				trellis->cell[i][j].submax->size = nrange[j];
				trellis->cell[i][j].submax->last =
					&trellis->cell[i][j].submax->elem[nrange[j]-1];
				
				for (u = 0; u < nrange[j]; u++){
					trellis->cell[i][j].submax->elem[u].score =
						trellis->cell[i][j].layer[0]->score;
					//	trellis->cell[j][i].score[0];
					trellis->cell[i][j].submax->elem[u].intrinsic_score =
						trellis->cell[i][j].layer[0]->score;
					//	trellis->cell[j][i].score[0];
					trellis->cell[i][j].submax->elem[u].from_state = -1;
					trellis->cell[i][j].submax->elem[u].end = i;
					trellis->cell[i][j].submax->elem[u].start = i;
				}
			}*/
		}
	}





	printf("zRUN1\n");
//	zShowTrellis(trellis);
	printf("%d,%d\n",dna->length,trellis->dna->length );
	/* induction - first pass */
	zTrace2("beginning first induction\n");

	for (i = PADDING; i < dna->length; i++) {
	//	zAllocViterbiVarsSingle(trellis,i);
	//	if (NULL == trellis->cell[i]) continue;//reverse added
		//Added by Jim
//		zInitOriginList(origin, hmm->states, top_num);
//		printf("Initial %d\n",i);
//		zShowOrigin(origin, hmm->states, top_num);
		/* create external links ending at i */
	//	printf("MMMM0:%d   %d\n",i,trellis->hmm->feature_count);
		for (j = 0; j < trellis->hmm->feature_count; j++) {

			/////////printf("---------------hmm->feature_count\n");
			/////////printf("XXXX1:%d\n",trellis->hmm->feature_count);
			l = i;
			if(trellis->factory[j] != NULL){
//				printf("PPPPPPP:%d,%d\n",i,j);
				zResetSFList(trellis->external[j]);		
				trellis->factory[j]->create3(trellis->factory[j], l, sfl);
				f = zSFListMoveFirst(sfl);
				while(f != NULL){
					f->strand = '+';	
					f->state  = j;
				//	printf("ppp%d,%d\n", i,j);
					zSFListInsert(trellis->external[j], f);
					f = zSFListMoveNext(sfl);
				}
				zResetSFList(sfl);
			} 
			// if (trellis->external[j]!=NULL && trellis->external[j]->size>0 )
			// {
			// 	printf("external %d\n",j );
			// 	zPrintSFList(trellis->external[j]);
			// }
	
	//		printf("XXXX2:%d\n",j);
			l = trellis->dna->length - i - 1;
			if(trellis->rfactory[j] != NULL){
			//	printf("QQQQQQQQ:%d,%d\n",i,j);
				zResetSFList(trellis->rexternal[j]);		
				trellis->rfactory[j]->create5(trellis->rfactory[j], l, sfl);
				f = zSFListMoveFirst(sfl);
				while(f != NULL){
					f->strand = '-';
					f->state  = j;
					tmp_start = f->start;
					f->start  = trellis->dna->length - f->end - 1;
					f->end    = trellis->dna->length - tmp_start - 1;
				//	printf("qqq%d,%d\n", i,j);
					zSFListInsert(trellis->rexternal[j], f);
					f = zSFListMoveNext(sfl);
				}
				zResetSFList(sfl);
			}
			// if (trellis->rexternal[j]!=NULL && trellis->rexternal[j]->size>0)
			// {
			// 	printf("rexternal %d\n",j );
			// 	zPrintSFList(trellis->rexternal[j]);
			// }
	
		}	
	//	 printf("external After loop\n" );
		// zPrintSFList(trellis->external[j]);
		// printf("rexternal After loop%d\n",);
		// zPrintSFList(trellis->external[j]);
//	if(show_mark==1)	printf("XXXX3:%d,%d\n",i,j);
		zFreeViterbiVarsSingleNew(trellis,i);
		for (j = 0; j < hmm->states; j++) {
			zAllocViterbiVarsSingle(trellis,i,j);
			if (NULL == trellis->cell[i]) continue;//reverse added

		//	if (NULL == trellis->cell[j]) continue;

			for (k = 0; k < top_num; ++k)
			{
				trellis->cell[i][j].layer[k]->score = MIN_SCORE;
				trellis->cell[i][j].layer[k]->trace = -1;
				trellis->cell[i][j].layer[k]->length = 0;
				trellis->cell[i][j].layer[k]->trace_top = -1;

				if (trellis->cell[i][j].layer[k]->parent!=NULL)//newly added
				{
					printf("aMAZING\n");
					printf("%d, %s,%d,%d\n",trellis->cell[i][j].layer[k]->length,zStrIdx2Char(trellis->hmm->state[trellis->cell[i][j].layer[k]->state].name) ,trellis->cell[i][j].layer[k]->pos,trellis->cell[i][j].layer[k]->story);
					printf("\t %d, %s,%d,%d\n",trellis->cell[i][j].layer[k]->parent->length,zStrIdx2Char(trellis->hmm->state[trellis->cell[i][j].layer[k]->parent->state].name) ,trellis->cell[i][j].layer[k]->parent->pos,trellis->cell[i][j].layer[k]->parent->story);
				}
			}


//if(show_mark==1)printf("XXXX4:%d,%d\n",i,j);
			/*////////////kevin
			if ((hmm->state[j].type == EXPLICIT)
				|| (hmm->state[j].type == INTERNAL)){
				trellis->cell[i][j].submax = zMalloc(sizeof(zSFVec),
													 "zRunViterbiAndForward submax");
				zInitSFVec(trellis->cell[i][j].submax, nrange[j]);
				trellis->cell[i][j].submax->size = nrange[j];
				trellis->cell[i][j].submax->last =
					&trellis->cell[i][j].submax->elem[nrange[j]-1];
//if(show_mark==1)printf("XXXX5\n");
				for (u = 0; u < nrange[j]; u++){
					trellis->cell[i][j].submax->elem[u].score = MIN_SCORE;
					trellis->cell[i][j].submax->elem[u].intrinsic_score = MIN_SCORE;
					trellis->cell[i][j].submax->elem[u].from_state = -1;
					trellis->cell[i][j].submax->elem[u].start = i;
					trellis->cell[i][j].submax->elem[u].end = i;
				}
//if(show_mark==1)		printf("XXXX6\n");
			}*/

			jumps = hmm->jmap[j];

////if(i==2386)	printf("%d,%d:jumpsize %d\n",i,j,jumps->size );
			for(k = 0; k < jumps->size; k++) {
				////if(i==2386)printf("XXXX666:%d\n",k);
////if(i==2386)			printf("type:%d\n",hmm->state[j].type );
////if(i==2386)		printf("elem:%d\n", jumps->elem[k]);
				//////if(hmm->state[j].type == 3){
				//////printf("-------------------hmm->state[j]-----------\n");
				//////printf("%d\n",hmm->state[j].type);}
				zGetTransFunc(hmm->state[j].type) 
					(trellis, jumps->elem[k], j, i);
			}
////if(i==2386)	printf("XXXX6.01:%d,%d\n",i,j);
			//Added by Jim
			zConnectTBTreeLayers(tree,trellis,i,j,PADDING-1);//new jim
////if(i==2386)printf("XXXX6.1:%d,%d--%d,%d\n",i,j,(hmm->state[j].type == EXTERNAL),(trellis->cell[i][j].layer[0]->score > MIN_SCORE));
			/*	if(!zOption("top")){ */
			if ((hmm->state[j].type == EXTERNAL)
				&& (trellis->cell[i][j].layer[0]->score > MIN_SCORE)){ 
			////if(i==2386)printf("XXXX6.15:%d,%d\n",i,j);
	//		   	&& (trellis->cell[j][i].score[0] > MIN_SCORE)){ 
			   	zPushIVec(&trellis->extpos[j], i);
			}
////if(i==2386)printf("XXXX6.2:%d,%d\n",i,j);
			/*	}
				else{
				if(hmm->state[j].type == EXTERNAL){
					for(m=0;m<top_num;m++){
						if(trellis->cell[j][i].score[m] > MIN_SCORE){
							zPushIVec(&trellis->extpos[m][j]);
						}
					}
					}*/




			/*//////////////kevin
			if ((hmm->state[j].type == EXPLICIT)
				|| (hmm->state[j].type == INTERNAL)){
				if (trellis->cell[i-1][j].submax->elem != NULL){
			//	printf("NANI %d,%d\n",trellis->cell[i-1][j].submax->size,trellis->cell[i-1][j].submax->elem[0].state);
					
					zFreeSFVec(trellis->cell[i-1][j].submax);
			//		if (i>614){printf("NANI2 %d\n",trellis->cell[i-1][j].submax->size);}
					trellis->cell[i-1][j].submax->elem = NULL;
				}
			//	printf("NANI3\n");
				zFree(trellis->cell[i-1][j].submax);
				trellis->cell[i-1][j].submax = NULL;
			}*/
//if(i==2386)printf("XXXX6.3:%d,%d\n",i,j);
		}



		window=3000;
		zTreeterbiNew(trellis,i,window,tmp_cpoint);


		if (i>=old_cpoint+CPOINT_WINDOW*(zzz+1)+2000*skip )
		{

			zTreeterbiReverseCPoint(trellis, i-CPOINT_WINDOW,tree);
			zTreeterbiClean(trellis,i-CPOINT_WINDOW,tree);

			size = sizeof(zSFVec);
	    	for(g=0;g<top_num;g++){
				tmp_sfv[g]=zMalloc(size,"zMalloc tmp_sfv[i]");
				zInitSFVec(tmp_sfv[g],10);
			}


			printf("haha %d\n",i);

			tmp_cpoint=zFindCPointNew(trellis,i-CPOINT_WINDOW,tree,tmp_cpoint);

			printf("FIN %d\n",  tmp_cpoint);
			printf("CHNUM:%d\n",i );

	//		if (tmp_cpoint==-1 || tmp_cpoint-old_cpoint < CPOINT_WINDOW*2/5)
				if (tmp_cpoint==-1 || tmp_cpoint-old_cpoint < 1)
			{
				printf("asdasd %d\n",tmp_cpoint);
				tmp_cpoint = old_cpoint;
	//			cpoint_check_skip++;
				skip++;
				continue;
			}
			else{cpoint_check_skip=0;skip=0;}

			
			zTracePartialTrellisNew(trellis, PADDING, tmp_cpoint, -1, tmp_sfv, tmp_path_score, tree);
			zFilterTmpSFV(tmp_sfv);

			if (zOption("a"))
			{
				printf("Start FILL\n");
				zFillSFVSponge(sponge,tmp_sfv,zCountRealSFV(tmp_sfv,top_num));
		//	zShowSFVSponge(sponge,trellis);
		//	zScoreSFVSponge(sponge,trellis);
			

				printf("After FILL\n");
			}


			for (g = 0; g < top_num; ++g)
			{
				if (tmp_sfv[g]->size!=0)
				{
					printf("%d-size:%d\n",g,tmp_sfv[g]->size );
				}
				
			}

			
		//	zFreeViterbiVarsRange(trellis, old_cpoint,tmp_cpoint);
		//	zTreeterbiForDrawing(trellis,i);
			zFreeViterbiVarsRange(trellis, older_cpoint, old_cpoint);
		older_cpoint = old_cpoint;
		old_cpoint = tmp_cpoint;
			// zTracePartialTrellisNew(trellis, PADDING, tmp_cpoint, -1, tmp_sfv, tmp_path_score, tree);
			// for(g=0;g<5;g++)printf("HAHAH%g\n",tmp_path_score[g] );


			// for (g = 0; g< 49; ++g)
			// {
			// 	printf("AA %d ?? %s\n",g,zStrIdx2Char(hmm->state[g].name));
			// }
			
	

			for (g = 0; g < tmp_sfv[0]->size; ++g)
			{
				printf("spring %s:%d\n",zStrIdx2Char(hmm->state[tmp_sfv[0]->elem[g].state].name),tmp_sfv[0]->elem[g].end );
				if (tmp_sfv[0]->elem[g].end==tmp_cpoint)
				{
					break;
				}
			}

		//	zTrimTmpSFV(tmp_sfv,tmp_cpoint);


//			g--;
			tmp_cpoint_state=tmp_sfv[0]->elem[g].state;
//			i=tmp_cpoint-1;
			printf("CHECKDOMAIN %d,%s\n",tmp_cpoint,zStrIdx2Char(hmm->state[tmp_cpoint_state].name));
			i=tmp_cpoint;

			
		// 	if (tmp_cpoint_state<17)
		// //	if (tmp_cpoint_state==14 || tmp_cpoint_state==15)//for utr3- and utr5
		// 	{
		// 	 	i=tmp_cpoint;
		// 	// 	tmp_state_length=tmp_sfv[0]->elem[g].end-tmp_sfv[0]->elem[g].start;//not +1 because we retreated 1 step
		// 		tmp_cpoint_state=tmp_sfv[0]->elem[g].state;
		// 		printf("autumn %d,%s\n",tmp_cpoint,zStrIdx2Char(hmm->state[tmp_cpoint_state].name));
		// 	//		trellis->cell[i][hmm->state[tmp_cpoint_state].name].layer[0]->length);//,hmm->state[trellis->cell[i][tmp_cpoint_state].layer[0]->trace].name);

		// 	}
		// 	else{
		// //		g--;
		// 		tmp_cpoint_state=tmp_sfv[0]->elem[g].state;
		// 		i=tmp_cpoint;
	
		// 		printf("winter %d,%s\n",tmp_cpoint,zStrIdx2Char(hmm->state[tmp_cpoint_state].name));
		// 	//		trellis->cell[i][hmm->state[tmp_cpoint_state].name].layer[0]->length);//,hmm->state[trellis->cell[i][tmp_cpoint_state].layer[0]->trace].name);

		// 		// i=tmp_cpoint-1;
		// 		// tmp_cpoint_state=tmp_sfv[0]->elem[g].state;
		// 		// printf("winter %s\n",zStrIdx2Char(hmm->state[tmp_cpoint_state].name));
		// 	}

			// if (tmp_cpoint_state<17){printf("XIAOYU17\n");}
			// else{printf("DAYU17\n");}

		//	printf("CHECKDOMAIN:%d\n",tmp_cpoint);
		//	printf("Root children%d\n",tree->root->children );
			// if (tree->root->children>0)
			// {
			// 	zTBTreeNode *z;
			// 	z=tree->root->child;
			// 	while(z!=NULL)
			// 	{
			// 		printf("%d\n",z->pos );
			// 		z=z->rsib;
			// 		break;
			// 	}
			// }
			zTBTreeClearSubTree(tree,tree->root);
			// printf("Root new children%d\n",tree->root->children );
			// if (true)
			// {
			// 	zTBTreeNode *z;
			// 	z=tree->root->child;
			// 	while(z!=NULL)
			// 	{
			// 		printf("%d\n",z->pos );
			// 		z=z->rsib;
			// 		break;
			// 	}
			// }
	//		i=tmp_cpoint;
		//	i=tmp_cpoint-1;

			show_mark+=1;
			printf("showmark fin %d\n",show_mark );
			if (!zOption("a"))
			{
				zTracePartialTrellisMerge(trellis,huge_sfv,tmp_sfv,show_mark,large_top,path_score);
			}

		
if(show_mark>1)		printf("VVV1 %d\n",i);


			zFreeViterbiVarsSingleNew(trellis,i);
			for (j = 0; j < hmm->states; j++)
			{
		//		if(show_mark>1)		printf("j:%d\n",j);
			//	zResetViterbiVarsSingle(trellis,i,j);
				zAllocViterbiVarsSingle(trellis,i,j);
		//		if (trellis->cell[i-1]==NULL) printf("HHHHHHH\n");
		//		else printf("AAAAA\n");
			//	zAllocViterbiVarsSingleForCPoint(trellis,i,j);
				for(m=0;m<top_num;m++){
			//		if(show_mark>1)		printf("m:%d\n",m);
					trellis->cell[i][j].layer[m]->score = MIN_SCORE;
					trellis->cell[i][j].layer[m]->trace = -1;
					trellis->cell[i][j].layer[m]->length = -1;
					trellis->cell[i][j].layer[m]->trace_top = -1;

					trellis->cell[i][j].layer[m]->pos=i;
					trellis->cell[i][j].layer[m]->left_pos=0;
					zTBTreeSetChild(tree->root,trellis->cell[i][j].layer[m]);
							
				}
				tree->root->pos=tmp_cpoint;
		//		zShowRootChildren(tree);
				// zTBTreeNode *z;
				// z=tree->root->child;
				// int abc=0;
				// while(z!=NULL)
				// {
				// 	abc++;
				// 	z=z->rsib;
				// }
				// printf("ABC: %d\n",abc);
			//	printf("Root atarashi children%d\n",tree->root->children );
				if (j==tmp_cpoint_state)
				{
					printf("JJJJJ: %d,%d\n",i,j );
				//	trellis->cell[i][tmp_sfv[0]->last->state]=cell_temp;
					trellis->cell[i][j].layer[0]->score = 0.0;
					trellis->cell[i][j].layer[0]->trace = j;
					trellis->cell[i][j].layer[0]->length =tmp_state_length;//new corrected
					trellis->cell[i][j].layer[0]->trace_top = 0;
					
					trellis->cell[i][j].layer[0]->left_pos=1;
				//	zFreeSFVec(trellis->cell[i][j].layer[0]->sfv);
					
				}
		//		printf("BEFORE submax\n");
				


				/*///////////kevin
				if ((hmm->state[j].type == EXPLICIT)
					|| (hmm->state[j].type == INTERNAL)){
					trellis->cell[i][j].submax = zMalloc(sizeof(zSFVec),
														 "zRunViterbiAndForward submax");
					zInitSFVec(trellis->cell[i][j].submax, nrange[j]);
					trellis->cell[i][j].submax->size = nrange[j];
					trellis->cell[i][j].submax->last =
						&trellis->cell[i][j].submax->elem[nrange[j]-1];
					for (u = 0; u < nrange[j]; u++){
						trellis->cell[i][j].submax->elem[u].score = MIN_SCORE;
				//			trellis->cell[i][j].layer[0]->score;
						trellis->cell[i][j].submax->elem[u].intrinsic_score = MIN_SCORE;
				//			trellis->cell[i][j].layer[0]->score;
						trellis->cell[i][j].submax->elem[u].from_state = -1;
						trellis->cell[i][j].submax->elem[u].start = i;
						trellis->cell[i][j].submax->elem[u].end = i;
					}
				}
				*/
//printf("aFTER submax\n");
			//	trellis->cell[i][tmp_sfv[0]->last->state].submax=cell_temp.submax;
				

				/*////////kevin
				if ((hmm->state[j].type == EXPLICIT)
					|| (hmm->state[j].type == INTERNAL)){
					if (trellis->cell[i-1]!=NULL && trellis->cell[i-1][j].submax!=NULL)
					{
						if (trellis->cell[i-1][j].submax->elem != NULL){
						zFreeSFVec(trellis->cell[i-1][j].submax);
						trellis->cell[i-1][j].submax->elem = NULL;
						}
						zFree(trellis->cell[i-1][j].submax);
						trellis->cell[i-1][j].submax = NULL;
					}	
				}
				
				*/

			}

		//	if(show_mark>1)	
			for(g=0;g<top_num;g++){
				zFreeSFVec(tmp_sfv[g]);
				zFree(tmp_sfv[g]);
				tmp_sfv[g]=NULL;
			}

			printf("VVV2\n");

		//	return tmp_sfv;
		}


		
	}

//Jim






	window=window-1-PADDING;
	for(;window>0;window--)
	{
		zTreeterbiNew(trellis,dna->length-1-PADDING,window,0);
	}
	printf("zRUN2\n");

	printf("mark1\n");
//	zPlot3DShowTrellis(trellis);

//	zShowTrellis(trellis);
	printf("mark2\n");
	//Added by Jim
//	zShowTBTree(tree);
	zFreeSFList(sfl);
	zFree(sfl);
	printf("mark3\n");
//printf("ZZZZ1:%d\n",j);
	/* Viterbi trace back */
	zTrace2("traceback\n");
//	sfv = zCalloc(top_num,sizeof(zSFVec*),"zRunViterbiAndForward sfv");
// 	sfv = zCalloc(LARGE_TOP,sizeof(zSFVec*),"zRunViterbiAndForward sfv");
// 	size = sizeof(zSFVec);
// //    for(g=0;g<top_num;g++){
// 	for(g=0;g<LARGE_TOP;g++){
// 		sfv[g]=zMalloc(size,"zMalloc sfv[i]");
// 		zInitSFVec(sfv[g],10);
// 	}
// 	printf("ZZZZ3:%d\n",j);

//zTraceTrellis(trellis, -1, sfv, path_score);
printf("mark4\n");

size = sizeof(zSFVec);
for(g=0;g<top_num;g++){
	tmp_sfv[g]=zMalloc(size,"zMalloc tmp_sfv[i]");
	zInitSFVec(tmp_sfv[g],10);
}
zTraceTrellisNew(trellis, -1, tmp_sfv, path_score, tree);

//			zTreeterbiForDrawing(trellis,i);
			// printf("Root children%d\n",tree->root->children );
			// if (tree->root->children>0)
			// {
			// 	zTBTreeNode *z;
			// 	z=tree->root->child;
			// 	while(z!=NULL)
			// 	{
			// 		printf("%d\n",z->pos );
			// 		z=z->rsib;
			// 		break;
			// 	}
			// }




show_mark+=1;
show_mark=-show_mark;
printf("showmark fin %d\n",show_mark );

printf("Mario\n");
if (!zOption("a"))
{
	for (i = 0; i < huge_sfv[0]->size; ++i)
	{
		printf("%d: %d-%d\n",huge_sfv[0]->elem[i].state,huge_sfv[0]->elem[i].start,huge_sfv[0]->elem[i].end);
	}
}

 zFilterTmpSFV(tmp_sfv);

if (zOption("a"))
{
	zFillSFVSponge(sponge,tmp_sfv,zCountRealSFV(tmp_sfv,top_num));
	zScoreSFVSponge(sponge,trellis);
}

if (!zOption("a"))
{
	zTracePartialTrellisMerge(trellis,huge_sfv,tmp_sfv,show_mark,large_top,path_score);
}

for(g=0;g<top_num;g++){
	zFreeSFVec(tmp_sfv[g]);
	zFree(tmp_sfv[g]);
	tmp_sfv[g]=NULL;
}
zFree(tmp_sfv);




printf("Odyssey\n");
if (!zOption("a"))
{
	for (i = 0; i < huge_sfv[0]->size; ++i)
	{
		printf("%d: %d-%d\n",huge_sfv[0]->elem[i].state,huge_sfv[0]->elem[i].start,huge_sfv[0]->elem[i].end);
	}
	//Added by Jim
	//Modified Score
	zModifySortScore(trellis,huge_sfv,path_score,large_top);

	for(g=0;g<large_top;g++){
		zFlipSFVec(huge_sfv[g]);
		zTranslateSFVec(huge_sfv[g], -PADDING);
	}

}


//zTraceTrellisSort(huge_sfv,sfv,path_score);  
//return sfv;
// printf("SCORELIST:\n");
// for (i = 0; i < huge_sfv[0]->size; ++i)
// //for (i = 0; i < LARGE_TOP; ++i)
// {
// 	printf("%d-%d: %f\n",huge_sfv[0]->elem[i].start,huge_sfv[0]->elem[i].end,huge_sfv[0]->elem[i].score );
// }


zFreeViterbiVarsRange(trellis,older_cpoint,dna->length);


zFree(nrange);
printf("before free nrange\n");
zFreeTreeNode(&tree->root);
zFree(tree);
printf("after free nrange\n");
return huge_sfv;

printf("mark5\n");
zTreeterbiForDrawing(trellis,trellis->dna->length-PADDING-1);
printf("mark5.5\n");
//zExtraClean(trellis);
printf("mark5.6\n");

//zShowTrellis(trellis);
printf("mark5.7\n");
//zCheckCPoint(trellis,dna->length-1-PADDING,tree);
printf("mark6\n");
//zCheckDomain(trellis,-1,2000,tree);

printf("mark7\n");



	zPlot3DShowTrellis(trellis);

printf("ZZZZ4:%d\n",j);

//	printf("ZZZZ5:%d\n",j);
	/*
 *path_score = trellis->cell[sfv->elem[0].state][dna->length-1-PADDING].score[0];
 */
    /*//////kevin
	for (j = 0; j < hmm->states; j++){
		if ((hmm->state[j].type == EXPLICIT)
			|| (hmm->state[j].type == INTERNAL)){
			if (trellis->cell[dna->length-1][j].submax->elem != NULL){
				zFreeSFVec(trellis->cell[dna->length-1][j].submax);
				trellis->cell[dna->length-1][j].submax->elem = NULL;
			}
			zFree(trellis->cell[dna->length-1][j].submax);

			trellis->cell[dna->length-1][j].submax = NULL;

		}
	}*/

	

	return sfv;
}

void zPlot3DShowTrellis(zTrellis *trellis){

	zHMM         *hmm = trellis->hmm;
	/*	zDNA         *dna = trellis->dna; */

	coor_t        i;        /* iterator for sequence */
	int           j, k,m;        /* iterator for internal states */

	
	score_t tmp_score[hmm->states*top_num];
	int mark;
	mark=0;
	score_t to_score,from_score;
	coor_t i1,i2,i3;
	int j1,j2,k1,k2;

//	zTrellisCellLayer* layer_temp;
zTBTreeNode* layer_temp;
	j = 0;
	printf("===================== 3D Show Trellis ===================================\n");
	printf("Top num: %d\nDNA length:%d\nHmm states:%d\nHmm feature count:%d\nHmm models:%d\n",top_num,trellis->dna->length,trellis->hmm->states,trellis->hmm->feature_count,trellis->hmm->models);

 	int nullCount=0;
 	int totalSpace=(trellis->dna->length-100)*top_num*49;
	FILE *fp;
 	fp=fopen("output.txt","w");

	for (i = 50; i < trellis->dna->length-50; ++i)
	{
		if(trellis->cell[i]==NULL){printf("%d GONE\n",i);continue;}
		for (k = 0; k < 49; ++k)  
		{	
			if (trellis->cell[i][k].layer==NULL){printf("%d,%d GONE\n",i,k);continue;}
			for (m = 0; m < top_num; ++m)
			{//	printf("%d,%d,%d\n",i,k,m );
				if (trellis->cell[i][k].layer[m]==NULL)
				{
				//	fprintf(fp, "-13\t-13\t-13\t-13");
				//	printf("null\n");
					nullCount+=1;
					continue;
				}
				else{
					if (trellis->cell[i][k].layer[m]->score==MIN_SCORE)
					{//printf("null1\n");
				//		printf("|z,%d,",trellis->cell[i][k].layer[m]->length);
				//		fprintf(fp, "-13\t%d\t",trellis->cell[k][i].layer[m]->length);
				nullCount++;
				continue;
					}
					else{//printf("null2\n");
					//	printf("|%3.2f,%d,",trellis->cell[i][k].layer[m]->score,trellis->cell[i][k].layer[m]->length);
					//	printf("%3.2f\t%d\t%d\t%d\t", trellis->cell[i][k].layer[m]->score,
					//		trellis->cell[i][k].layer[m]->pos,tris->cell[i][k].layer[m]->state,trellis->cell[i][k].layer[m]->story);
						fprintf(fp, "%3.2f\t%d\t%d\t%d\t", trellis->cell[i][k].layer[m]->score,
							trellis->cell[i][k].layer[m]->pos,trellis->cell[i][k].layer[m]->state,trellis->cell[i][k].layer[m]->story);
					}
					if (trellis->cell[i][k].layer[m]->trace==-1 || trellis->cell[i][k].layer[m]->parent==NULL || trellis->cell[i][k].layer[m]->trace==trellis->cell[i][k].layer[m]->state)
					{//printf("null3\n");
					//	printf("non,%d", trellis->cell[i][k].layer[m]->trace_top);
						fprintf(fp,"-13\t-13\t-13\t");
					}else{//printf("null4\n");
					//	printf("%s,%d",zStrIdx2Char(hmm->state[trellis->cell[i][k].layer[m]->trace].name),trellis->cell[i][k].layer[m]->trace_top);
						fprintf(fp,"%d\t%d\t%d\t",trellis->cell[i][k].layer[m]->length,trellis->cell[i][k].layer[m]->trace,trellis->cell[i][k].layer[m]->trace_top);

					}
					if (trellis->cell[i][k].layer[m]->left_pos==1)//specific function
					{
						fprintf(fp, "1\n" );
					}
					else{fprintf(fp, "0\n" );}
				}

	//			fprintf(fp, "\t");
			}
	//		printf("|\n");
	//		fprintf(fp, "\n");
			
		}
	// 	printf("\n");

	}
	printf("======================= 3D End Trellis =================================\n");
	double saveRatio=(double)nullCount/(double)totalSpace*100;
	printf("NULL count:%d\n",nullCount );
	printf("Total Space:%d\n",totalSpace );
	printf("Saved Space:%f%%\n",saveRatio );
	fclose(fp);
}
void zShowTrellisSingle(zTrellis* trellis, coor_t i)
{
	int k ,j,m;
	for (k = 0; k < 49; ++k)  
	{			
		printf("%d--%s\n", k,zStrIdx2Char(trellis->hmm->state[k].name));

		if(trellis->cell[i]==NULL)continue;
		printf("%d: ",i);//printf("%d:", trellis->cell[i][k].layer[0]);
		for (m = 0; m < top_num; ++m)
		{
			if (trellis->cell[i][k].layer[m]==NULL)
			{
			//	printf("|NULL" );
			//	fprintf(fp, "-13\t-13\t-13\t-13");
			//	nullCount+=1;
			}
			else{
				if (trellis->cell[i][k].layer[m]->score==MIN_SCORE)
				{
					continue;
				}
				else{
					printf("|%3.2f,%d,",trellis->cell[i][k].layer[m]->score,trellis->cell[i][k].layer[m]->length);
				}
				if (trellis->cell[i][k].layer[m]->trace==-1)
				{
					printf("non,%d", trellis->cell[i][k].layer[m]->trace_top);
				}else{
					printf("%s,%d",zStrIdx2Char(trellis->hmm->state[trellis->cell[i][k].layer[m]->trace].name),trellis->cell[i][k].layer[m]->trace_top);

				}
			}
			printf("|\n");			
		}
	 	printf("\n");
	}
}

void zShowTrellis(zTrellis *trellis){
	zHMM         *hmm = trellis->hmm;
	/*	zDNA         *dna = trellis->dna; */

	coor_t        i;        /* iterator for sequence */
	int           j, k,m;        /* iterator for internal states */

	
	score_t tmp_score[hmm->states*top_num];
	int mark;
	mark=0;
	score_t to_score,from_score;
	coor_t i1,i2,i3;
	int j1,j2,k1,k2;

//	zTrellisCellLayer* layer_temp;
zTBTreeNode* layer_temp;
	j = 0;
	printf("===================== Show Trellis ===================================\n");
	printf("Top num: %d\nDNA length:%d\nHmm states:%d\nHmm feature count:%d\nHmm models:%d\n",top_num,trellis->dna->length,trellis->hmm->states,trellis->hmm->feature_count,trellis->hmm->models);
//	zWriteHMM(stdout,hmm);
	// for (k = 0; k < hmm->states; k++) {
	// 	printf("%d\t",k );
	// 	zWriteHMM_State(stdout, &hmm->state[k], hmm->iso_states);
	// }
	// for(k= 0; k <  80; k++) {
	// 	//printf("%s\t", zStrIdx2Char(k));
	// 	printf("%d: %s, ", k,zStrIdx2Char(hmm->state[k].name));
	// //	if(hmm->state[j].type==INTERNAL){}

	//  }
	// printf("\n");
 	int nullCount=0;
 	int totalSpace=36000*top_num*49;
// 	FILE *fp;
 //	fp=fopen("output.txt","w");

	for (k = 0; k < 49; ++k)  
	{			printf("%d--%s\n", k,zStrIdx2Char(hmm->state[k].name));
		for (i = 13200; i < trellis->dna->length-50; ++i)
		{	if(trellis->cell[i]==NULL)continue;
			printf("%d: ",i);//printf("%d:", trellis->cell[i][k].layer[0]);
			for (m = 0; m < top_num; ++m)
			{
				if (trellis->cell[i][k].layer[m]==NULL)
				{
				//	printf("|NULL" );
				//	fprintf(fp, "-13\t-13\t-13\t-13");
					nullCount+=1;
				}
				else{
					if (trellis->cell[i][k].layer[m]->score==MIN_SCORE)
					{
					//	printf("|z,%d,",trellis->cell[i][k].layer[m]->length);
				//		fprintf(fp, "-13\t%d\t",trellis->cell[i][k].layer[m]->length);
						continue;
					}
					else{
						printf("|%3.2f,%d,",trellis->cell[i][k].layer[m]->score,trellis->cell[i][k].layer[m]->length);
					//	fprintf(fp, "%3.2f\t%d\t", trellis->cell[i][k].layer[m]->score,trellis->cell[i][k].layer[m]->length);
					}
					if (trellis->cell[i][k].layer[m]->trace==-1)
					{
						printf("non,%d", trellis->cell[i][k].layer[m]->trace_top);
					//	fprintf(fp,"-13\t%d", trellis->cell[i][k].layer[m]->trace_top);
					}else{
						printf("%s,%d",zStrIdx2Char(hmm->state[trellis->cell[i][k].layer[m]->trace].name),trellis->cell[i][k].layer[m]->trace_top);
					//	fprintf(fp,"%d\t%d",trellis->cell[i][k].layer[m]->trace,trellis->cell[i][k].layer[m]->trace_top);

					}
				}
				// layer_temp=zFromCellLayer(trellis,i,k,m);
				// if (layer_temp==NULL)
				// {
				// 	printf("NULL|");
				// }
				// else {		printf("%d,%s|", layer_temp->length,zStrIdx2Char(hmm->state[layer_temp->trace].name));}
				// layer_temp=zToCellLayer(trellis,i,k,m,300);
				// if (layer_temp==NULL)
				// {
				// 	printf("NULL");
				// }
				// else {		printf("%d,%s", layer_temp->length,zStrIdx2Char(hmm->state[layer_temp->trace].name));}
			//	fprintf(fp, "\t");
			}
			printf("|\n");
		//	fprintf(fp, "\n");
			
		}
	 	printf("\n");

	}
	printf("======================= End Trellis =================================\n");
	double saveRatio=(double)nullCount/(double)totalSpace*100;
	printf("NULL count:%d\n",nullCount );
	printf("Total Space:%d\n",totalSpace );
	printf("Saved Space:%f%%\n",saveRatio );

}
void zRunBackward(zTrellis *trellis) {
	zHMM*         hmm = trellis->hmm;
	coor_t        i;        /* iterator for sequence */
	coor_t        length = trellis->dna->length;
	int           j;        /* iterator for internal states */
	int           k;        /* iterator for previous states */
	
	/* Init */
	zAllocBackwardVars(trellis);
	for (j = 0; j < hmm->states; j++) {
		for (i = length - PADDING; i < length; i++) {
			trellis->backward[j][i] = 0;
		}
	}

	/* Induction */
	for (i = length - PADDING - 1; i > 1; i--) {
		for (j = 0; j < hmm->states; j++) {
			if (NULL == trellis->backward[j]) continue;
			for (k = 0; k < hmm->states; k++) {
				if (MIN_SCORE == zGetTransitionScore(trellis->hmm, j, k,trellis->tiso_group)) continue;
				zGetBackTransFunc(hmm->state[k].type)(trellis, j, k, i);
			}
		}
	} 
}

float  zForwardBackwardProb(zTrellis* trellis, zSfeature* f) {
	/*zHMM    *hmm = trellis->hmm;
	  zIVec   *states;
	  int      i, pre_state, state, post_state;
	  zPhase_t  phs;
	  score_t  parse_score, dscore, fscore, t1score, t2score, phscore;
	  score_t  pre_score, post_score;
	  score_t  result = MIN_SCORE;
	*/
	/* find the correct HMM state for the feature */
	/*	states = &hmm->simap[hmm->somap[f->name]];
	for (i = 0; i < states->size; i++) {
		state = states->elem[i];
		phs   = hmm->state[state].phase;
		if (zCompatibleJumpFromExon(trellis->dna, phs, f)) continue;
	}
	if (i == states->size) {
		zWarn("couldn't find the corresponding HMM state");
		return 0.0;
	}
	
	dscore = zScoreDuration(hmm->dmap[hmm->state[state].duration].score, 
							f->end - f->start + 1);
	fscore = f->score;
	*/
	/* sum over all paths that use this feature */
	/*	for (pre_state = 0; pre_state < hmm->states; pre_state++) {
		t1score = hmm->tmap[pre_state][post_state];
		pre_score = trellis->forward[pre_state][f->start-1];
		if (MIN_SCORE == t1score) continue;
		for (post_state = 0; post_state < hmm->states; post_state++) {
			t2score = hmm->tmap[pre_state][post_state];
			post_score = trellis->backward[post_state][f->end+1];
			phscore = zFhmm->phasepref
			Ack!! get the previous imp. from 8/15/02!
	*/
	trellis = NULL; /*compiler hush*/
	f = NULL;       /*compiler hush*/
	return 0.0;
}

score_t zGetIntergenicContinueScore(zTrellis* trellis) {
	int iso_group = trellis->tiso_group;

	if (iso_group < 0 || iso_group >= trellis->hmm->iso_transitions ) {
		zDie("zGetIntergenicContinueScore: iso_group %d not found\n",
			 iso_group);
	}
	return trellis->hmm->inter_continue[iso_group];
}

void zClearColumn(zTrellis* trellis, int pos) {
	int state;
	for (state=0; state < trellis->hmm->states; ++state) {
		// trellis->cell[state][pos].length[0] = 0;
		// trellis->cell[state][pos].score[0]  = MIN_SCORE;
		// trellis->cell[state][pos].trace[0] = -1;

		trellis->cell[pos][state].layer[0]->length = 0;
		trellis->cell[pos][state].layer[0]->score  = MIN_SCORE;
		trellis->cell[pos][state].layer[0]->trace = -1;
	}
}

score_t zGetPathScore(zTrellis* trellis, const zSFVec* sfv, zSFVec* result) {
	int             cur_state, prev_state;
	int             i, j;
	coor_t          begin, end, pos;
	zHMM           *hmm = trellis->hmm;
	zHMM_StateType  state_type;
	zSFVec         **external;/*, *exter;
							   zFeatureFactory *ff;*/
	score_t  *score;
	/* initialize trellis */
	end = PADDING-1;
	zAllocViterbiVars(trellis);
	for (pos = 0; pos <= end; ++pos) {
		zClearColumn(trellis, pos);
	}
	for (cur_state = 0; cur_state < hmm->states; cur_state++) {
		// trellis->cell[cur_state][end].score[0] = zGetInitProb(trellis->hmm, cur_state, trellis->iiso_group);
		// trellis->cell[cur_state][end].trace[0] = cur_state;
		trellis->cell[end][cur_state].layer[0]->score = zGetInitProb(trellis->hmm, cur_state, trellis->iiso_group);
		trellis->cell[end][cur_state].layer[0]->trace = cur_state;
	}

	/* loop set up */
	if (sfv->size > 0) 
		cur_state = sfv->elem[0].state; /* for prev_state, below */

	/* trace through each feature */
	for (i = 0; i < sfv->size; ++i) {
		zSfeature      *f;
        external = zMalloc(sizeof(zSFVec), "external in zGetPathScore");
        zInitSFVec(external[0], 1);

		/* get vitals on this feature */
		f = &sfv->elem[i];
		prev_state = cur_state;
		cur_state  = f->state;
		begin = end + 1;

		/* run Viterbi if necessary */
		if (begin != f->start) {
			coor_t pos;
			for (pos = begin; pos < f->start; pos++) {
				zClearColumn(trellis, pos);
			}
			zRunPartialViterbiAndForward(trellis, begin, f->start-1);
			prev_state = -1;
		}
		end = f->end;
		
		/* transition from prev_state to cur_state */
		state_type = hmm->state[cur_state].type;
		if (INTERNAL == state_type || GINTERNAL == state_type ) { 
			/*Internal or Ginternal transitions go one base at a time */
			zClearColumn(trellis, begin);
			if (prev_state < 0) {
				zRunPartialViterbiAndForward(trellis, f->start, f->start);
			} else {
				zGetTransFunc(state_type)(trellis,prev_state,cur_state,begin);
			}
			for(pos = f->start+1; pos <= end; ++pos) {
				zClearColumn(trellis, pos);
				zGetTransFunc(state_type)(trellis, cur_state, cur_state, pos);
			//	trellis->cell[cur_state][pos-1].score[0] = MIN_SCORE;
				trellis->cell[pos-1][cur_state].layer[0]->score = MIN_SCORE;
			//	if (trellis->cell[cur_state][pos].score[0] < -10000) {
				if (trellis->cell[pos][cur_state].layer[0]->score < -10000) {
					printf("messed up on feature %d, at %u\n", i, pos);
				}
			}
		} else { /* Other states make jumps */
			for (pos = f->start; pos <= end; ++pos) {
				zClearColumn(trellis, pos);
			}
			/*EVAN I'm guessing that this comment breaks this so I will need to fix it later
			if (EXTERNAL == state_type) { / put only that feature in the vec/
				exter = ('-' == hmm->state[cur_state].strand) ?
					&trellis->rexternal[hmm->state[cur_state].model][end] :
					&trellis->external[hmm->state[cur_state].model][end];
				zFreeSFVec(exter);
				zInitSFVec(exter, 1);
				ff = ('-' == f->strand)
					? trellis->rfactory[hmm->state[cur_state].model]
					: trellis->factory[hmm->state[cur_state].model];
				if ('-' == f->strand) zAntiSfeature(f, trellis->dna->length);
				f->score = ff->scoref(ff, f);
				if ('-' == f->strand) zAntiSfeature(f, trellis->dna->length);
				zPushSFVec(exter, f);
			}
			*/
			if (prev_state < 0) {
				zRunPartialViterbiAndForward(trellis, end, end);
			} else {
				zGetTransFunc(state_type)(trellis, prev_state, cur_state, end);
			}
		}

		/* gather the results of the transitions */
		zTracePartialTrellis(trellis, begin, end, cur_state, external, score);
		if (NULL == external) {
			zWarn("Couldn't trace back from state %s, is the feature legal?",
				  zStrIdx2Char(hmm->state[cur_state].name));
			return MIN_SCORE;
		}
		for (j = external[0]->size; j > 0; j--) {
			zPushSFVec(result, &external[0]->elem[j-1]);
		}
		zFreeSFVec(external[0]);
		zFree(external[0]);
		external[0] = NULL;

		/* clean up so that it's impossible to jump back here */
		for (pos = begin-1; pos < end; pos++) {
			zClearColumn(trellis, pos);
		}

		/* and a sanity check */
	//	if (trellis->cell[cur_state][end].score[0] < -10000) {
		if (trellis->cell[end][cur_state].layer[0]->score < -10000) {
			zWarn("Model prohibits ending in %s, at %d", 
				  zStrIdx2Char(hmm->state[f->state].name), 
				  end - trellis->padding + 1);
			return MIN_SCORE;
		}
	}

	/* run one last viterbi, if needed */
	if (end < trellis->dna->length - PADDING - 1) {
		coor_t pos;
		external = zMalloc(sizeof(zSFVec), "external in zGetPathScore");
		zInitSFVec(external[0], 1); 

		begin = end + 1;
		end = trellis->dna->length - PADDING - 1;
		for (pos = begin; pos < end; pos++) {
			zClearColumn(trellis, pos);
		}
		zRunPartialViterbiAndForward(trellis, begin, end);
		zTracePartialTrellis(trellis, begin, end, -1, external,score);
		for (j = external[0]->size; j > 0; j--) {
			zPushSFVec(result, &external[0]->elem[j-1]);
		}
		zFreeSFVec(external[0]);
		zFree(external[0]);
		external[0] = NULL;
	}
	assert(end == trellis->dna->length - PADDING - 1);

	/* report end score */
//	return trellis->cell[result->last->state][end].score[0];
	return trellis->cell[end][result->last->state].layer[0]->score;
}

frame_t zFragFrame2Char(frame_t lfrag, frame_t rfrag, frame_t frame)
{
	/* Encodes lfrag, rfrag and frame values in a character      */
	/* in bitwise manner, as shown below. All values range       */
	/* between -1 and 2, hence 2 bits are enough to encode each. */

	/* --------------------------------------------------------- */
	/* bit index          :   7 6     5 4     3 2     1 0        */
	/*                                                           */
	/* character (1 byte) :   0 0 --- 0 1 --- 1 0 --- 0 0        */
	/*                        ^ ^     ^ ^     ^ ^     ^ ^        */
	/*                        | |     | |     | |     | |        */
	/* encoded info       :  unused  frame   rfrag   lfrag       */
	/*                                                           */
	/* values             :            1       2       0         */
	/* --------------------------------------------------------- */

	frame_t result = 0; /* frame_t is typedef'ed as char */

	if ((lfrag > 2) || (rfrag > 2) || (frame > 2))
		{
			fprintf(stderr, "lfrag, rfrag or frame is out of bounds: %d %d %d\n",
					lfrag, rfrag, frame);
		}
	else
		{
			/* value of -1 is allowed for promoters */
			/* and PolyAs and is encoded here as 3  */
			if (lfrag < 0) lfrag = 3;
			if (rfrag < 0) rfrag = 3;
			if (frame < 0) frame = 3;

			result |= lfrag;        /* set bits 0 and 1 using lfrag */

			result |= (rfrag << 2); /* set bits 2 and 3 using rfrag */
			/* shifted left by 2 bits       */

			result |= (frame << 4); /* set bits 4 and 5 using frame */
			/* shifted left by 4 bits       */
		}

	return(result);
}

void zChar2FragFrame(frame_t frag_frame, zSfeature *f)
{
	/* Reads lfrag, rfrag and frame values encoded in character frag_frame */
	/* and assigns those to the corresponding members of zSfeature object. */
	/* The bytes are read two at a time by means of bitwise AND operation  */
	/* with a number that has 1's in the bits to be read and 0's elsewhere */

  /* ------------------------------------------------------------------- */
  /* bit index          :   7 6     5 4     3 2     1 0                  */
  /*                                                                     */
  /* character (1 byte) :   0 0 --- 0 1 --- 1 0 --- 0 0                  */
  /*                                                & &                  */
  /*     & 3            :                           1 1 = 0 0 (lfrag=0)  */
  /*                                        & &                          */
  /*     & (3 << 2)     :                   1 1 >> 2    = 1 0 (rfrag=2)  */
  /*                                & &                                  */
  /*     & (3 << 4)     :           1 1 >> 4            = 0 1 (frame=1)  */
  /* ------------------------------------------------------------------- */

	f->lfrag = frag_frame & 3;               /* set lfrag from bits 0 and 1 */
	f->rfrag = (frag_frame & (3 << 2)) >> 2; /* set rfrag from bits 2 and 3 */
	f->frame = (frag_frame & (3 << 4)) >> 4; /* set frame from bits 4 and 5 */

  
	/* value of -1 is allowed for promoters */
	/* and PolyAs and is enco------- */

	f->lfrag = frag_frame & 3;               /* set lfrag from bits 0 and 1 */
	f->rfrag = (frag_frame & (3 << 2)) >> 2; /* set rfrag from bits 2 and 3 */
	f->frame = (frag_frame & (3 << 4)) >> 4; /* set frame from bits 4 and 5 */

  
	/* value of -1 is allowed for promoters */
	/* and PolyAs and is encoded here as 3  */
	if (f->lfrag == 3) f->lfrag = -1;
	if (f->rfrag == 3) f->rfrag = -1;
	if (f->frame == 3) f->frame = -1;
}


/**********************************************************\
  Pin Data Search Code                                   
	EVAN
\**********************************************************/

/* EVAN removed checking code because it is old style, need to fix it and test feature creation 
static void check_for_feature(zSfeature *f,zSFList *sfl,int erase){
	zSfeature *f2;
	f2 = zSFListMoveFirst(sfl);
	while(f2 != NULL){
		if(f->start == f->start && f->end == f2->end && f->score == f2->score && f->lfrag == f2->lfrag && f->rfrag == f2->rfrag && f->frame == f2->frame){
			if(erase != 0){
				zSFListRemoveCurrent(sfl);
			}
			return;
		}
		f2 = zSFListMoveNext(sfl);
	}
	if(erase){
		fprintf(stderr,"Missed %d (%d->%d) %c\n",f->name,f->start,f->end,f->strand);
	}
	else{
		fprintf(stderr,"Created bad %d (%d->%d) %c\n",f->name,f->start,f->end,f->strand);
	}
}

static void check_for_all_features(zTrellis *trellis,zSFList *sfl, coor_t pos, int fnum, strand_t strand){
	coor_t i,end;
	zSFList* sfl2;
	zSfeature *f2;
	
	end = pos + 50000;
	if(end > trellis->dna->length) end = trellis->dna->length -10;
	for(i = pos; i < end;i++){
		if(strand == '+'){
			sfl2 = zGetTrellisExternalValue(trellis,i,fnum);
		}
		else{
			sfl2 = zGetTrellisRExternalValue(trellis,i,fnum);
		}
		f2 = zSFListMoveFirst(sfl2);
		while(f2 != NULL){
			if(f2->start < pos){
				check_for_feature(f2,sfl,1);
			}
			f2 = zSFListMoveNext(sfl2);
		}
	}
}

static void check_inside_create(zTrellis *trellis){
	int j;
	coor_t i;
	zSFList sfl;
	zSfeature *f;

	zInitSFList(&sfl);

	zComputeBackLinks2(trellis);
	dprintf(0,"External vals ready\n");
	for(i = 10;i <= trellis->dna->length-10; i++){
		if(i % 100 == 0){
			dprintf(9,"Checked up to %d\n",i);
		}
		for(j = 0;j < trellis->hmm->feature_count;j++){
			if(trellis->factory[j] != NULL){
				trellis->factory[j]->icreate(trellis->factory[j],i,&sfl);
				f = zSFListMoveFirst(&sfl);
				while(f != NULL){
					f->strand = '+';
					check_for_feature(f,zGetTrellisExternalValue(trellis,f->end,j),0);
					f = zSFListMoveNext(&sfl);
				}
				check_for_all_features(trellis,&sfl,i,j,'+');
				zResetSFList(&sfl);
			}
			if(trellis->rfactory[j] != NULL){
				trellis->factory[j]->icreate(trellis->rfactory[j],i,&sfl);
				f = zSFListMoveFirst(&sfl);
				while(f != NULL){
					f->strand = '-';
					check_for_feature(f,zGetTrellisRExternalValue(trellis,f->end,j),0);
					f = zSFListMoveNext(&sfl);
				}
				check_for_all_features(trellis,&sfl,i,j,'-');
				zResetSFList(&sfl);
			}
		}
	}

static void check_link_dist(zTrellis *trellis) {
	coor_t     i;
	int        j;
	coor_t     tmp_start;
	zSFList   *sfl;
	zSfeature *f;

	int idx;
	coor_t split_size = 5000;
	int    split_count = trellis->dna->length/split_size + 1;
	int*   count = zMalloc(sizeof(int)*split_count,"");
	for(j = 0;j < split_count;j++){
		count[j] = 0;
	}

	if (NULL != trellis->external) return;

	trellis->external = zMalloc(sizeof(zExternal),
															"zComputeBackLinks trellis->external");
	zInitExternal(trellis->external,trellis->factory,1,
								trellis->hmm->feature_count);
	trellis->rexternal = zMalloc(sizeof(zExternal),
															 "zComputeBackLinks trellis->rexternal");
	zInitExternal(trellis->rexternal,trellis->rfactory,1,
								trellis->hmm->feature_count);
	sfl = zMalloc(sizeof(zSFList), "zComputeBackLinks sfl");
	zInitSFList(sfl);
  
	for (j = 0; j < trellis->hmm->feature_count; j++) {
		if(trellis->factory[j] != NULL){
			for (i = 0; i < trellis->dna->length; i++) {
				trellis->factory[j]->create(trellis->factory[j], i, sfl);
				f = zSFListMoveFirst(sfl);
				while(f != NULL){
					f->strand = '+';	

					idx = (int)f->start/split_size;
					while(idx*split_size <= f->end){
						if(idx*split_size >= f->start){
							count[idx]++;
						}
						idx++;
					}
					f = zSFListMoveNext(sfl);
				}
				zResetSFList(sfl);
			}
		}
		if(trellis->rfactory[j] != NULL){
			for (i = 0; i < trellis->dna->length; i++) {		
				trellis->rfactory[j]->create(trellis->rfactory[j], i, sfl);
				f = zSFListMoveFirst(sfl);
				while(f != NULL){
					f->strand = '-';
					tmp_start = f->start;
					f->start  = trellis->dna->length - f->end - 1;
					f->end    = trellis->dna->length - tmp_start - 1;

					idx = (int)f->start/split_size;
					while(idx*split_size <= f->end){
						if(idx*split_size >= f->start){
							count[idx]++;
						}
						idx++;
					}
					f = zSFListMoveNext(sfl);
				}
				zResetSFList(sfl);
			}
		}
	}
  
	zFreeSFList(sfl);
	zFree(sfl);

	for(j = 0;j < split_count;j++){
		dprintf(-1000,"pos\t%d\t%d\n",j*split_size,count[j]);
	}
}
*/

struct zPinPathData{
	score_t score;
	int     pinstate;
};
typedef struct zPinPathData zPinPathData;

static void zInitPinPathData(zPinPathData* p){
	p->score = MIN_SCORE;
	p->pinstate = -1;
}

struct zPinArrayData{
	int    index;
	int    state;
	coor_t pos;
	/*the following are needed only for 5' long distance links*/
	int link_state;
};
typedef struct zPinArrayData zPinArrayData;

static void* zCreatePinArrayData(){
	zPinArrayData* pad = zMalloc(sizeof(zPinArrayData),"zCreateInt pad");
	pad->index = -1;
	pad->state = -1;
	pad->pos = (coor_t)-1;
	return (void*)pad;
}

static void zFreePinArrayData(void* v){
	zFree(v);
}

static void zResetPinArrayData(void* v){
	zPinArrayData* pad = (zPinArrayData*)v;
	pad->index = -1;
	pad->state = -1;
	pad->pos = (coor_t)-1;
}

static void* zCreateInt(){
	int* i = zMalloc(sizeof(int),"zCreateInt i");
	return i;
}

static void zFreeInt(void* v){
	zFree(v);
}

static void zResetInt(void* v){
	int* i = (int*)v;
	*i = 0;
}

static void zInitIList(zList* l){
	zInitList(l,zCreateInt,zFreeInt,zResetInt);
}

static void zIListAddLast(zList* l, int i){
	int *j;
	j = zListAddLast(l);
	*j = i;
}

static void zIListAddFirst(zList* l, int i){
	int *j;
	j = zListAddFirst(l);
	*j = i;
}

static int zIListPopFirst(zList* l){
	int *j;
	j = zListMoveFirst(l);
	zListRemoveFirst(l);
	if(j == NULL){
		zDie("Poped empty zIList");
	}
	return *j;
}

struct zPinSearch{
	zTrellis*          trellis;          /* zTrellis object */
	zPinPathData**     data;             /* matrix holding all path data */
	int                states;           /* number of states */
	int                size;             /* data matrix size */
	zList              live_to_list;     /* list of live path start nodes
											maps path start to zPinSearch.data row
											sorted by decreasing seq pos */
	zList              live_from_list;   /* list of live path end nodes
											maps path end to zPinSearch.data column
											sorted by increasing seq pos */
	zList              free_to_list;     /* list of unused rows in zPinSearch.data */
	zList              free_from_list;   /* list of unused columns in zPinSearch.data */
	zPinPathData*      tmp;               /* temp buffer to hold data when stepping 
											 forward of backward by one pos */
	coor_t             pin_pos;          /* the pos at which we are searching for a pin */
	coor_t             start_pos;        /* the current maximum path start pos */
	coor_t             end_pos;          /* the current minimum path end pos 
											at any time we know all path which can get us 
											from start_pos to end_pos, including any that 
											include a long distance link across either 
											start_pos or end_pos */
	int                pin_state;        /* the pin state at pin_pos, this is -1 until we
											find the correct value */
	zSFList            sfl;              /* a temp zSfeature list to avoid repeated 
											allocation and initialization */ 
};
typedef struct zPinSearch zPinSearch;

static int INIT_PIN_SEARCH_SIZE = 2500;
static const coor_t PIN_SEARCH_STEPS = 5000;

static void zStepPinSearchBack(zPinSearch*);
static void zStepPinSearchForward(zPinSearch*);
static void zPinSearchCheckForPin(zPinSearch*);
static void zStepPinArrayForward(zPinSearch*, zPinArrayData*);
static void zStepPinArrayBack(zPinSearch*, zPinArrayData*);
static void zPinSearchCreateLiveNodes3(zPinSearch*);
static void zPinSearchCreateLiveNodes5(zPinSearch*);
static void zPinSearchAddFeature3(zPinSearch*, int, strand_t);
static void zPinSearchAddFeature5(zPinSearch*, int, strand_t);
static void zPinSearchAddFeatureInit(zPinSearch*, int, strand_t);
static zPinArrayData* zPinSearchGetToPad(zPinSearch*, coor_t, int);
static zPinArrayData* zPinSearchGetFromPad(zPinSearch*, coor_t, int);

/* EVAN currently ignoring the long distance links which cross the pin_pos */
/* EVAN stupidly adding multiple list nodes for the same start/end state-pos combo */

static void zInitPinSearch(zPinSearch* ps, zTrellis* trellis, coor_t pin_pos){
	int i,j;
	zPinArrayData* pad;
	ps->trellis = trellis;
	ps->states = trellis->hmm->states;
	ps->pin_pos = pin_pos;
	ps->start_pos = pin_pos;
	ps->end_pos = pin_pos;
	ps->pin_state = -1;

	/*EVAN no current mechanism for increasing ps->size.  If it isn't big enough
		program will exit with error via zDie()*/
	ps->size = INIT_PIN_SEARCH_SIZE;
	ps->data = zMalloc(sizeof(zPinPathData*)*ps->size,"zFindTrellisPin data");
	ps->tmp = zMalloc(sizeof(zPinPathData)*ps->size,"zFindTrellisPin tmp");
	for(i = 0;i < ps->size;i++){
		zInitPinPathData(&ps->tmp[i]);
		ps->data[i] = zMalloc(sizeof(zPinPathData)*ps->size,"zFindTrellisPin data[i]");
		for(j = 0;j < ps->size;j++){
			zInitPinPathData(&ps->data[i][j]);
		}
	}

	zInitList(&ps->live_to_list,zCreatePinArrayData,zFreePinArrayData,zResetPinArrayData);
	zInitList(&ps->live_from_list,zCreatePinArrayData,zFreePinArrayData,
						zResetPinArrayData);
	zInitIList(&ps->free_to_list);
	zInitIList(&ps->free_from_list);
	
	/* Add row/col to live list for each state, i, at pos i+1 (adjusted for tmp at 0) */
	for(i = 0;i < ps->states;i++){
		pad = zListAddLast(&ps->live_to_list);
		pad->index = i;
		pad->pos = pin_pos;
		pad->state = i;
		pad = zListAddLast(&ps->live_from_list);
		pad->index = i;
		pad->pos = pin_pos;
		pad->state = i;
		
		if(trellis->hmm->state[i].type == INTERNAL ||
		   trellis->hmm->state[i].type == GINTERNAL ){
			ps->data[i][i].pinstate = i;
			ps->data[i][i].score = 0;
		}
		else if(trellis->hmm->state[i].type == EXTERNAL){
			ps->data[i][i].pinstate = -1;
			ps->data[i][i].score = MIN_SCORE;
		}
		else{
			zDie("Pin Search does not support states other than INTERNAL and EXTERNAL");
		}
	}
	/* Add the rest of the rows/cols to the free list */ 
	for(i = ps->states;i < ps->size-1;i++){
		zIListAddLast(&ps->free_to_list,i);
		zIListAddLast(&ps->free_from_list,i);
	}

	zInitSFList(&ps->sfl);

	/* Add long distance links crossing pin_pos */
	for (i = 0; i < ps->trellis->hmm->feature_count; i++) {
		j = ps->pin_pos;
		if(trellis->factory[i] != NULL){
			/*dprintf(0,"Adding Feature %d at pos %d +\n",i,j);*/
			trellis->factory[i]->createi(trellis->factory[i], j, &ps->sfl);
			zPinSearchAddFeatureInit(ps,i,'+');
			zResetSFList(&ps->sfl);
		}
		j = trellis->dna->length - 1 - ps->pin_pos;
		if(trellis->rfactory[i] != NULL){
			/*dprintf(0,"Adding Feature %d at pos %d -\n",i,j);*/
			trellis->rfactory[i]->createi(trellis->rfactory[i], j, &ps->sfl);
			zPinSearchAddFeatureInit(ps,i,'-');
			zResetSFList(&ps->sfl);
		}		
	}	
	fprintf(stderr,"Done adding Features\n");

}

static void zFreePinSearch(zPinSearch* ps){
	int i;
	zFreeSFList(&ps->sfl);
	zFreeList(&ps->live_to_list);
	zFreeList(&ps->live_from_list);
	zFreeList(&ps->free_to_list);
	zFreeList(&ps->free_from_list);
	zFree(ps->tmp);
	ps->tmp = NULL;
	for(i = 0;i < ps->size;i++){
		zFree(ps->data[i]);
		ps->data[i] = NULL;
	}
	zFree(ps->data);
	ps->data = NULL;
}

static void zPinSearchResetTo(zPinSearch* ps, int i){
	zPinArrayData* pad;
	if(i >= ps->size){
		return;
	}
	pad = zListMoveFirst(&ps->live_from_list);
	while(pad != NULL){
		ps->data[pad->index][i].score = MIN_SCORE;
		ps->data[pad->index][i].pinstate = -1;
		pad = zListMoveNext(&ps->live_from_list);
	}
}

static void zPinSearchResetFrom(zPinSearch* ps, int i){
	zPinArrayData* pad;
	if(i >= ps->size){
		return;
	}
	pad = zListMoveFirst(&ps->live_to_list);
	while(pad != NULL){
		ps->data[i][pad->index].score = MIN_SCORE;
		ps->data[i][pad->index].pinstate = -1;
		pad = zListMoveNext(&ps->live_to_list);
	}
}

static void zPinSearchResetTmp(zPinSearch* ps){
	int i;
	for(i = 0;i < ps->states;i++){
		ps->tmp[i].score = MIN_SCORE;
		ps->tmp[i].pinstate = -1;
	}
}

static zPinArrayData* zPinSearchGetToPad(zPinSearch* ps, coor_t pos, int state){
	zPinArrayData* to_pad;

	to_pad = zListGetCurrent(&ps->live_to_list);
	if(to_pad == NULL){
		/* start at end of list (highest coor)*/
		to_pad = zListMoveLast(&ps->live_to_list);
	}
	/* if the list isn't empty */
	if(to_pad != NULL){
		if(to_pad->pos > pos){
			/* move back until we are at a node <= pos */
			while(to_pad != NULL && to_pad->pos > pos){
				to_pad = zListMovePrev(&ps->live_to_list);
			}
			/* while the current node == pos, see if it is the correct state */
			while(to_pad != NULL && to_pad->pos == pos){
				if(to_pad->state == state){
					return to_pad;
				}
				to_pad = zListMovePrev(&ps->live_to_list);
			}
			/* didn't find the node in the list so add a new one */
			to_pad = zListAddNext(&ps->live_to_list);
		}
		else{
			/* move forward until we are at a node >= pos */
			while(to_pad != NULL && to_pad->pos < pos){
				to_pad = zListMoveNext(&ps->live_to_list);
			}
			/* while the current node == pos, see if it is the correct state */
			while(to_pad != NULL && to_pad->pos == pos){
				if(to_pad->state == state){
					return to_pad;
				}
				to_pad = zListMoveNext(&ps->live_to_list);
			}
			/* didn't find the node in the list so add a new one */
			to_pad = zListAddPrev(&ps->live_to_list);
		}
	}
	/* fill in the new nodes info */
	to_pad->index = zIListPopFirst(&ps->free_to_list);
	to_pad->state = state;
	to_pad->pos = pos;
	/*dprintf(0,"ADD TOPAD %d %d\n",pos,state);EVAN*/
	return to_pad;
}

static zPinArrayData* zPinSearchGetFromPad(zPinSearch* ps, coor_t pos, int state){
	zPinArrayData* from_pad;

	from_pad = zListGetCurrent(&ps->live_from_list);
	if(from_pad == NULL){
		/* start at end of list (lowest coor)*/
		from_pad = zListMoveLast(&ps->live_from_list);
	}
	/* if the list isn't empty */
	if(from_pad != NULL){
		if(from_pad->pos < pos){
			/* move back until we are at a node >= pos */
			while(from_pad != NULL && from_pad->pos < pos){
				from_pad = zListMovePrev(&ps->live_from_list);
			}
			/* while the current node == pos, see if it is the correct state */
			while(from_pad != NULL && from_pad->pos == pos){
				if(from_pad->state == state){
					return from_pad;
				}
				from_pad = zListMovePrev(&ps->live_from_list);
			}
			/* didn't find the node in the list so add a new one */
			from_pad = zListAddNext(&ps->live_from_list);
		}
		else{
			/* move forward until we are at a node <= pos */
			while(from_pad != NULL && from_pad->pos > pos){
				from_pad = zListMoveNext(&ps->live_from_list);
			}
			/* while the current node == pos, see if it is the correct state */
			while(from_pad != NULL && from_pad->pos == pos){
				if(from_pad->state == state){
					return from_pad;
				}
				from_pad = zListMoveNext(&ps->live_from_list);
			}
			/* didn't find the node in the list so add a new one */
			from_pad = zListAddPrev(&ps->live_from_list);
		}
	}
	/* fill in the new nodes info */
	from_pad->index = zIListPopFirst(&ps->free_from_list);
	from_pad->state = state;
	from_pad->pos = pos;
	/*dprintf(0,"ADD FROMPAD %d %d\n",pos,state);EVAN*/
	return from_pad;
}

void zCheckForwardSearch(zPinSearch* ps){
	coor_t i;
	int keep_state = 0;
	zPinArrayData *to_pad,*from_pad;
	coor_t lead_in = 100000;

	ps->start_pos -= lead_in;
	ps->end_pos -= lead_in;
	to_pad = zListMoveFirst(&ps->live_to_list);
	while(to_pad != NULL){
		to_pad->pos -= lead_in;
		to_pad = zListMoveNext(&ps->live_to_list);
	}
	from_pad = zListMoveFirst(&ps->live_from_list);
	while(from_pad != NULL){
		from_pad->pos -= lead_in;

		if(from_pad->state == 9999){/*keep_state){*/	
			/*remove long distance link node */
			/*dprintf(0,"\tremoving %d\n",from_pad->state);*/
			zPinSearchResetFrom(ps,from_pad->index);
			zIListAddFirst(&ps->free_from_list,from_pad->index);
			zListRemoveCurrent(&ps->live_from_list);
		}

		from_pad = zListMoveNext(&ps->live_from_list);
	}

	if(ps->start_pos > ps->trellis->dna->length){
		zDie("I suck\n");
	}

	for(i = 0;i < lead_in;i++){
		zStepPinSearchForward(ps);
	}
	fprintf(stderr,"Search forward to %d (%d->%d)\n",ps->end_pos,ps->live_from_list.size,ps->live_to_list.size);
	zPinSearchCheckForPin(ps);

	fprintf(stderr,"Setting pin at pos %d\n",ps->end_pos);
	fprintf(stderr,"trimming to only state %d at start (%d)\n",keep_state,ps->start_pos);
	from_pad = zListMoveFirst(&ps->live_from_list);
	while(from_pad != NULL){
		to_pad = zListMoveFirst(&ps->live_to_list);
		while(to_pad != NULL){
			if(ps->data[from_pad->index][to_pad->index].score != MIN_SCORE){
				/*dprintf(0,"\t\t(%d->%d) (%d->%d) pin = %d->%d %f\n",from_pad->state,to_pad->state,from_pad->pos,to_pad->pos,ps->data[from_pad->index][to_pad->index].pinstate,to_pad->state,ps->data[from_pad->index][to_pad->index].score);			*/
				ps->data[from_pad->index][to_pad->index].pinstate = to_pad->state;			
			}
			to_pad = zListMoveNext(&ps->live_to_list);
		}
		from_pad = zListMoveNext(&ps->live_from_list);
	}
	fprintf(stderr,"Search forward to %d (%d->%d)\n",ps->end_pos,ps->live_from_list.size,ps->live_to_list.size);
	zPinSearchCheckForPin(ps);
	while(ps->pin_state < 0){
		for(i = 0;i < PIN_SEARCH_STEPS;i++){
			zStepPinSearchForward(ps);
		}
		fprintf(stderr,"Search forward to %d (%d->%d)\n",ps->end_pos,ps->live_from_list.size,ps->live_to_list.size);
		zPinSearchCheckForPin(ps);
	}
}


void zCheckBackSearch(zPinSearch* ps){
	coor_t i;
	int keep_state = 0;
	zPinArrayData *to_pad,*from_pad;
	coor_t lead_in = 100000;

	ps->start_pos += lead_in;
	ps->end_pos += lead_in;
	to_pad = zListMoveFirst(&ps->live_to_list);
	while(to_pad != NULL){
		to_pad->pos += lead_in;
		if(to_pad->state == 9999){/*keep_state){*/
			zPinSearchResetTo(ps,to_pad->index);
			zIListAddFirst(&ps->free_to_list,to_pad->index);
			zListRemoveCurrent(&ps->live_to_list);
		}
		to_pad = zListMoveNext(&ps->live_to_list);
	}
	from_pad = zListMoveFirst(&ps->live_from_list);
	while(from_pad != NULL){
		from_pad->pos += lead_in;

		from_pad = zListMoveNext(&ps->live_from_list);
	}

	if(ps->end_pos > ps->trellis->dna->length){
		zDie("I suck\n");
	}
	
	for(i = 0;i < lead_in;i++){
		zStepPinSearchBack(ps);
	}
	fprintf(stderr,"Search back to %d (%d->%d)\n",ps->start_pos,ps->live_from_list.size,ps->live_to_list.size);
	zPinSearchCheckForPin(ps);

	fprintf(stderr,"Setting pin at pos %d\n",ps->start_pos);
	fprintf(stderr,"trimming to only state %d at start (%d)\n",keep_state,ps->end_pos);
	to_pad = zListMoveFirst(&ps->live_to_list);
	while(to_pad != NULL){
		from_pad = zListMoveFirst(&ps->live_from_list);
		while(from_pad != NULL){
			if(ps->data[from_pad->index][to_pad->index].score != MIN_SCORE){
				/*fprintf(stderr,"\t\t(%d->%d) (%d->%d) pin = %d->%d %f\n",from_pad->state,to_pad->state,from_pad->pos,to_pad->pos,ps->data[from_pad->index][to_pad->index].pinstate,to_pad->state,ps->data[from_pad->index][to_pad->index].score);*/
				if(from_pad->pos < ps->start_pos){
					ps->data[from_pad->index][to_pad->index].pinstate = from_pad->link_state;			
					
				}
				else{
					ps->data[from_pad->index][to_pad->index].pinstate = from_pad->state;			
				}
			}
			from_pad = zListMoveNext(&ps->live_from_list);
		}
		to_pad = zListMoveNext(&ps->live_to_list);
	}
	fprintf(stderr,"Search back to %d (%d->%d)\n",ps->start_pos,ps->live_from_list.size,ps->live_to_list.size);
	zPinSearchCheckForPin(ps);
	while(ps->pin_state < 0){
		for(i = 0;i < PIN_SEARCH_STEPS;i++){
			zStepPinSearchBack(ps);
		}
		fprintf(stderr,"Search back to %d (%d->%d)\n",ps->start_pos,ps->live_from_list.size,ps->live_to_list.size);
		zPinSearchCheckForPin(ps);
	}
}



int zFindTrellisPin(zTrellis* trellis, coor_t pin_pos){
	coor_t i;
	zDNA*        dna;
	zPinSearch*  ps;
	int          pin_state = -1;
	
	pin_pos += PADDING - 1;

	dna = trellis->dna;
	if(pin_pos >= dna->length){
		zDie("zFindTrellisPin pin_pos must be less than the sequence length");
	}
	if(pin_pos == 0){
		zDie("zFindTrellisPin pin_pos must be greater than 0");
	}
	
	/* EVAN testing code
	if(pin_pos == 0){
		check_link_dist(trellis);
		check_inside_create(trellis);
		return -1;
	}
	*/

	/* init pin data structure */
	ps = zMalloc(sizeof(zPinSearch),"zFindTrellisPin ps");
	zInitPinSearch(ps,trellis,pin_pos);
	
	/* Until we have found the pin */
	zPinSearchCheckForPin(ps);
	/*EVAN*/
	/*zCheckBackSearch(ps);*/
	while(ps->pin_state < 0){
		for(i = 0;i < PIN_SEARCH_STEPS;i++){
			if(ps->end_pos >= trellis->dna->length - 1){
				zDie("Pin Search past end of sequence");
			}
			zStepPinSearchForward(ps);
		}
		fprintf(stderr,"Search forward to %d (%d->%d)\n",ps->end_pos,ps->live_from_list.size,ps->live_to_list.size);
		zPinSearchCheckForPin(ps);
		for(i = 0;i < PIN_SEARCH_STEPS;i++){
			if(ps->start_pos == 0){
				zDie("Pin Search past start of sequence");
			}
			zStepPinSearchBack(ps);
		}
		fprintf(stderr,"Search back to %d (%d->%d)\n",ps->start_pos,ps->live_from_list.size,ps->live_to_list.size);
		zPinSearchCheckForPin(ps);
	}
	pin_state = ps->pin_state;
	
	/* free pin data */ 
	zFreePinSearch(ps);
	zFree(ps);

	/* free trellis */
	zFreeTrellis(trellis);
	zFree(trellis);

	/* Return the pin state */
	/* EVAN what if ps->pin_pos != pin_pos */
	return pin_state;
}


static void zPinSearchCheckForPin(zPinSearch* ps){
	int i;
	int paths[100];
	int max_path_to_state[100];
	int min_path_to_state[100];
	int max_path_from_state[100];
	int min_path_from_state[100];
	score_t max_path_score[100];
	score_t min_path_score[100];
	int state = -1;
	coor_t last_from = -1;
	coor_t last_to = -1;
	int live_paths = 0;
	zPinArrayData* to_pad = zListMoveFirst(&ps->live_from_list);
	zPinArrayData* from_pad = zListMoveFirst(&ps->live_from_list);

	for(i = 0;i <= ps->states;i++){
		paths[i] = 0;
		max_path_score[i] = MIN_SCORE;
		min_path_score[i] = MAX_SCORE;
		max_path_to_state[i] = -1;
		min_path_to_state[i] = -1;
		max_path_from_state[i] = -1;
		min_path_from_state[i] = -1;
	}

	if(to_pad != NULL && from_pad != NULL){
		state = ps->data[from_pad->index][to_pad->index].pinstate;
		last_from = from_pad->pos;
		last_to = to_pad->pos;
	}
	else{
		zDie("Empty live_to_list or live_from_list");
	}

	if(from_pad->pos != ps->start_pos){
		zDie("Live from list not starting at start pos\n");
	}
	
	while(from_pad != NULL){
		to_pad = zListMoveFirst(&ps->live_to_list);
		if(to_pad->pos != ps->end_pos){
			zDie("Live to list not starting at end pos\n");
		}
		last_to = to_pad->pos;
		while(to_pad != NULL){
			if(ps->data[from_pad->index][to_pad->index].score > MIN_SCORE){
				i = ps->data[from_pad->index][to_pad->index].pinstate;
				live_paths++;
				paths[i]++;
				/*dprintf(2,"(%d,%d)->(%d,%d) = %d (%.4f)\tp[%d]=%d\n",from_pad->pos,from_pad->state,to_pad->pos,to_pad->state, ps->data[from_pad->index][to_pad->index].pinstate,ps->data[from_pad->index][to_pad->index].score,ps->data[from_pad->index][to_pad->index].pinstate,paths[ps->data[from_pad->index][to_pad->index].pinstate]);*/
				if(state != i){
					state = -1;
					ps->pin_state = -1;
					/*return;*/
				}
				if(max_path_score[i] < ps->data[from_pad->index][to_pad->index].score){
					max_path_score[i] = ps->data[from_pad->index][to_pad->index].score;
					max_path_to_state[i] = to_pad->state;
					max_path_from_state[i] = from_pad->state;
				}
				if(min_path_score[i] > ps->data[from_pad->index][to_pad->index].score){
					min_path_score[i] = ps->data[from_pad->index][to_pad->index].score;
					min_path_to_state[i] = to_pad->state;
					min_path_from_state[i] = from_pad->state;
				}
			}
			if(last_to > to_pad->pos){
				zDie("Out of order to");
			}
			last_to = to_pad->pos;
			to_pad = zListMoveNext(&ps->live_to_list);
		}
		if(last_from < from_pad->pos){
			zDie("Out of order from");
		}
		last_from = from_pad->pos;
		from_pad = zListMoveNext(&ps->live_from_list);
	}
	ps->pin_state = state;
	
	fprintf(stderr,"\tCHECK FOR_PIN: %d live paths\n",live_paths);
	for(i = 0;i <= ps->states;i++){
		if(paths[i] > 0){
			fprintf(stderr,"\tstate %d - %d\t%f (%d->%d)\t%f (%d->%d)\n",i,paths[i],min_path_score[i],min_path_from_state[i],min_path_to_state[i],max_path_score[i],max_path_from_state[i],max_path_to_state[i]);
		}
	}
}

static void zCheckPinSearchListOrder(zPinSearch* ps){
	zPinArrayData*  to_pad;
	zPinArrayData*  from_pad;
	coor_t last;
	to_pad = zListMoveFirst(&ps->live_to_list);
	last = to_pad->pos;
	while(to_pad != NULL){
		if(last > to_pad->pos){
			to_pad = zListMoveFirst(&ps->live_to_list);
			last = to_pad->pos;
			while(to_pad != NULL){
				if(last > to_pad->pos){
					fprintf(stderr,"\t**");
				}
				fprintf(stderr,"\t%d\t%d\n",to_pad->pos,to_pad->state);
				last = to_pad->pos;
				to_pad = zListMoveNext(&ps->live_to_list);
			}
			fprintf(stderr,"LAST\n");
			to_pad = zListMoveLast(&ps->live_to_list);
			last = to_pad->pos;
			while(to_pad != NULL){
				if(last < to_pad->pos){
					fprintf(stderr,"\t**");
				}
				fprintf(stderr,"\t%d\t%d\n",to_pad->pos,to_pad->state);
				last = to_pad->pos;
				to_pad = zListMovePrev(&ps->live_to_list);
			}

			
			zDie("WTF to\n");
		}
		last = to_pad->pos;
		to_pad = zListMoveNext(&ps->live_to_list);
	}
	from_pad = zListMoveFirst(&ps->live_from_list);
	last = from_pad->pos;
	while(from_pad != NULL){
		if(last < from_pad->pos){
			from_pad = zListMoveFirst(&ps->live_from_list);
			last = from_pad->pos;
			while(from_pad != NULL){
				if(last < from_pad->pos){
					fprintf(stderr,"\t**");
				}
				fprintf(stderr,"\t%d\t%d\n",from_pad->pos,from_pad->state);
				last = from_pad->pos;
				from_pad = zListMoveNext(&ps->live_from_list);
			}
			fprintf(stderr,"LAST\n");
			from_pad = zListMoveLast(&ps->live_from_list);
			last = from_pad->pos;
			while(from_pad != NULL){
				if(last > from_pad->pos){
					fprintf(stderr,"\t**");
				}
				fprintf(stderr,"\t%d\t%d\n",from_pad->pos,from_pad->state);
				last = from_pad->pos;
				from_pad = zListMovePrev(&ps->live_from_list);
			}

			
			zDie("WTF from\n");
		}
		last = from_pad->pos;
		from_pad = zListMoveNext(&ps->live_from_list);
	}
}

static void zStepPinSearchForward(zPinSearch* ps){
	zPinArrayData*  to_pad;
	zPinArrayData*  from_pad;

	/* Move end forward */
	ps->end_pos++;

	/* create new 3' live nodes (links starting at end_pos)*/
	zPinSearchCreateLiveNodes3(ps);

	/* foreach live from-node update path to each live to-node at time end_pos */
	from_pad = zListMoveFirst(&ps->live_from_list);
	while(from_pad != NULL){
		zStepPinArrayForward(ps,from_pad);
		from_pad = zListMoveNext(&ps->live_from_list);
	}

	/* update live_to_list states pad->pos (not long dist links) */
	to_pad = zListMoveFirst(&ps->live_to_list);
	while(to_pad != NULL && to_pad->pos == ps->end_pos-1){
		to_pad->pos++;
		to_pad = zListMoveNext(&ps->live_to_list);
	}

	/* update and remove 3' live nodes from long dist links (pad->pos == ps->end_pos)*/
	while(to_pad != NULL && to_pad->pos == ps->end_pos){
		/* incorporate long distance live nodes (EXTERNAL state)*/
		from_pad = zListMoveFirst(&ps->live_from_list);
		while(from_pad != NULL){
			if(ps->data[from_pad->index][to_pad->index].score > 
				 ps->data[from_pad->index][to_pad->state].score){				
				ps->data[from_pad->index][to_pad->state].score = 
					ps->data[from_pad->index][to_pad->index].score;
				ps->data[from_pad->index][to_pad->state].pinstate = 
					ps->data[from_pad->index][to_pad->index].pinstate;
			}
			from_pad = zListMoveNext(&ps->live_from_list);
		}
		/*remove long distance link node */
		zPinSearchResetTo(ps,to_pad->index);
		zIListAddFirst(&ps->free_to_list,to_pad->index);
		zListRemoveCurrent(&ps->live_to_list);
		to_pad = zListMoveNext(&ps->live_to_list);
	}
}

static void zStepPinSearchBack(zPinSearch* ps){
	zPinArrayData*  to_pad;
	zPinArrayData*  from_pad;

	/* Move start back */
	ps->start_pos--;

	/* create new 5' live nodes (links starting at start_pos)*/	
	zPinSearchCreateLiveNodes5(ps);

	/* foreach live to-node update path from each live from-node time start_pos */
	to_pad = zListMoveFirst(&ps->live_to_list);
	while(to_pad != NULL){
		zStepPinArrayBack(ps,to_pad);
		to_pad = zListMoveNext(&ps->live_to_list);
	}

	/* update live_from_list states pad->pos (not long dist links) */
	from_pad = zListMoveFirst(&ps->live_from_list);
	while(from_pad != NULL && from_pad->pos == ps->start_pos+1){
		from_pad->pos--;
		from_pad = zListMoveNext(&ps->live_from_list);
	}

	/* update and remove 5' live nodes from long dist links (pad->pos == ps->start_pos)*/
	while(from_pad != NULL && from_pad->pos == ps->start_pos){
		/* incorporate long distance live nodes (EXTERNAL state)*/
		to_pad = zListMoveFirst(&ps->live_to_list);
		while(to_pad != NULL){
			if(ps->data[from_pad->index][to_pad->index].score > 
				 ps->data[from_pad->state][to_pad->index].score){				
				ps->data[from_pad->state][to_pad->index].score = 
					ps->data[from_pad->index][to_pad->index].score;
				ps->data[from_pad->state][to_pad->index].pinstate = 
					ps->data[from_pad->index][to_pad->index].pinstate;
			}
			to_pad = zListMoveNext(&ps->live_to_list);
		}
		/*remove long distance link node */
		zPinSearchResetFrom(ps,from_pad->index);
		zIListAddFirst(&ps->free_from_list,from_pad->index);
		zListRemoveCurrent(&ps->live_from_list);
		from_pad = zListMoveNext(&ps->live_from_list);
	}
}

static void zStepPinArrayForward(zPinSearch* ps,zPinArrayData* from_pad){
	int            i,prestate,state;
	score_t        score;
	zIVec*         jumps;
	zHMM*          hmm = ps->trellis->hmm;

	/* reset tmp_from data */
	zPinSearchResetTmp(ps);
	
	/* Step each to state node forward */
	for(state = 0;state < ps->states;state++){
		if(hmm->state[state].type == INTERNAL ||
		   hmm->state[state].type == GINTERNAL){
			jumps = hmm->jmap[state];
			for(i = 0; i < jumps->size; i++) {
				prestate = jumps->elem[i];
				score = ps->data[from_pad->index][prestate].score;
				if(score == MIN_SCORE) continue;
				if(hmm->state[state].type == INTERNAL){
					score += zGetTransitionScore(ps->trellis->hmm,prestate,state,ps->trellis->tiso_group);
				}
				else{
					score += zGetIntergenicContinueScore(ps->trellis);
				}
				score += zScoreInternalState(ps->trellis,state,ps->end_pos,(state != prestate));
				if(score > ps->tmp[state].score){
					ps->tmp[state].score = score;
					ps->tmp[state].pinstate = ps->data[from_pad->index][prestate].pinstate;
				}
			}
		}
		else if(hmm->state[state].type == EXTERNAL){
			/* dealt with in zStepPinSearchForward */
		}
		else{
			zDie("zPinSearch only supports INTERNAL and EXTERNAL states");
		}
	}
	/* copy info from buffer back to from_pad->index */
	for(state = 0;state < ps->states;state++){
		ps->data[from_pad->index][state].score = ps->tmp[state].score;
		ps->data[from_pad->index][state].pinstate = ps->tmp[state].pinstate;
	}
}

static void zStepPinArrayBack(zPinSearch* ps,zPinArrayData* to_pad){
	int            i,poststate,state;
	score_t        score;
	zIVec*         jumps;
	zHMM*          hmm = ps->trellis->hmm;
	
	/* reset tmp_to data */
	zPinSearchResetTmp(ps);

	/* Step each from state node back */
	for(state = 0;state < ps->states;state++){
		jumps = hmm->fmap[state];
		for(i = 0; i < jumps->size; i++) {
			poststate = jumps->elem[i];
			if(hmm->state[poststate].type == INTERNAL ||
			   hmm->state[poststate].type == GINTERNAL){
				score = ps->data[poststate][to_pad->index].score;
				if(score == MIN_SCORE) continue;
				/*length 1 for duration score. it works properly since length 0 only need special
				  case for first position in sequence. see funtion in zTransition */
				if(state == 0 && poststate == 0){
					poststate = jumps->elem[i];
				}
				if(hmm->state[poststate].type == INTERNAL){
					score += zGetTransitionScore(ps->trellis->hmm,state,poststate,ps->trellis->tiso_group);
				}
				else{
					score += zGetIntergenicContinueScore(ps->trellis);
				}
				score += zScoreInternalState(ps->trellis,state,ps->start_pos,(state != poststate));
				if(score > ps->tmp[state].score){
					ps->tmp[state].score = score;
					ps->tmp[state].pinstate = ps->data[poststate][to_pad->index].pinstate;
				}
			}
			else if(hmm->state[poststate].type == EXTERNAL){
				/* dealt with when creating features (long distance links) 
					 intenrnal states preceding the external state correspinding to the 
					 feature are incorporated in zStepPinSearchBack */
			}
			else{
				zDie("zPinSearch only supports INTERNAL and EXTERNAL states");
			}
		}
	}
	/* copy info from buffer back to from_pad->index */
	for(state = 0;state < ps->states;state++){
		ps->data[state][to_pad->index].score = ps->tmp[state].score;
		ps->data[state][to_pad->index].pinstate = ps->tmp[state].pinstate;
	}
}

static void zPinSearchCreateLiveNodes5(zPinSearch* ps){
	int i;
	coor_t j;

	for (i = 0; i < ps->trellis->hmm->feature_count; i++) {
		j = ps->start_pos;
		if(ps->trellis->factory[i] != NULL){
			ps->trellis->factory[i]->create3(ps->trellis->factory[i], j, &ps->sfl);
			zPinSearchAddFeature5(ps,i,'+');
			zResetSFList(&ps->sfl);
		}
		j = ps->trellis->dna->length - ps->start_pos - 1;
		if(ps->trellis->rfactory[i] != NULL){
			ps->trellis->rfactory[i]->create5(ps->trellis->rfactory[i], j, &ps->sfl);
			zPinSearchAddFeature5(ps,i,'-');
			zResetSFList(&ps->sfl);
		}
	}	
}

static void zPinSearchCreateLiveNodes3(zPinSearch* ps){
	int i;
	coor_t j;

	for (i = 0; i < ps->trellis->hmm->feature_count; i++) {
		j = ps->end_pos;
		if(ps->trellis->factory[i] != NULL){
			ps->trellis->factory[i]->create5(ps->trellis->factory[i], j, &ps->sfl);
			zPinSearchAddFeature3(ps,i,'+');
			zResetSFList(&ps->sfl);
		}
		j = ps->trellis->dna->length - ps->end_pos - 1;
		if(ps->trellis->rfactory[i] != NULL){
			ps->trellis->rfactory[i]->create3(ps->trellis->rfactory[i], j, &ps->sfl);
			zPinSearchAddFeature3(ps,i,'-');
			zResetSFList(&ps->sfl);
		}
	}	
}

void zPinSearchAddFeature5(zPinSearch* ps, int fnum, strand_t strand){
	zTrellis*    trellis = ps->trellis;
	int          poss,i,j,k,state,prestate,poststate;
	score_t      score,base_score;
	zSfeature*   f;
	zIVec*       fmap;
	zIVec*       prejumps;
	zIVec*       postjumps;
	coor_t       tmp_start;
	zPhase_t     left_phase,right_phase;
	zPinArrayData* to_pad;
	zPinArrayData* from_pad;

	poss = ('-' != strand);
	if(poss){
		fmap = trellis->fmap3[fnum];
	}
	else{
		fmap = trellis->rfmap3[fnum];
	}
	
	f = zSFListMoveFirst(&ps->sfl);
	/* add new live node for each possible start state of each feature in sfl */
	while(f != NULL){
		f->strand = strand;	
		if(poss){
			right_phase = zSFEndPhase(trellis->dna, f);
			left_phase = zSFBeginPhase(trellis->dna, f);
		}
		else{
			tmp_start = f->start;
			f->start  = trellis->dna->length - f->end - 1;
			f->end    = trellis->dna->length - tmp_start - 1;			
			right_phase = zSFBeginPhase(trellis->rdna, f);
			left_phase = zSFEndPhase(trellis->rdna, f);
		}
		if(f->end != ps->start_pos){
			zDie("Cannot add external features to 5' end of pin search which don't end at start_pos");
		}
		/* for each (external) state which this feature can end in */
		for(i = 0;i < fmap->size;i++){
			state = fmap->elem[i];
			/* make sure external state (state) is compatible with feature */
			if(right_phase != trellis->hmm->state[state].phase) continue;
			/* for each (internal) state which can transition into this (external) state */ 
			prejumps = trellis->hmm->jmap[state];
			for(j = 0; j < prejumps->size;j++){
				prestate = prejumps->elem[j];
				if (!zCompatPhases(poss, trellis->hmm->state[prestate].phase, left_phase)) 
					continue;
				/* add a new from pad */
				from_pad = zPinSearchGetFromPad(ps,f->start-1,prestate);
				from_pad->link_state = state;
				/* for each (internal) state which this (external) state can transition to */
				postjumps = trellis->hmm->fmap[state];
				for(k = 0; k < postjumps->size;k++){
					poststate = postjumps->elem[k];
					base_score = zScoreInternalState(trellis,poststate,f->end+1,true) +
						zGetTransitionScore(ps->trellis->hmm,state,poststate,ps->trellis->tiso_group)+
						zScoreExternalState(trellis,state,f->end,f->end-f->start+1,f->score) + 
						zGetTransitionScore(ps->trellis->hmm,prestate,state,ps->trellis->tiso_group);
					
					/* select the best route to each to_pad.  each route can take any
						 (internal) prestate -> (external) state -> (internal) poststate */
					/*fprintf(stderr,"LDL %d -> %d (%d) -> %d\n",prestate,state,fnum,poststate);EVAN*/
					to_pad = zListMoveFirst(&ps->live_to_list);
					while(to_pad != NULL){
						score = base_score + ps->data[poststate][to_pad->index].score;
						if(score > ps->data[from_pad->index][to_pad->index].score){
							ps->data[from_pad->index][to_pad->index].score = score;
							ps->data[from_pad->index][to_pad->index].pinstate = 
								ps->data[poststate][to_pad->index].pinstate;
						}
						to_pad = zListMoveNext(&ps->live_to_list);
					}
				}
			}
		}
		f = zSFListMoveNext(&ps->sfl);
	}
}

void zPinSearchAddFeature3(zPinSearch* ps, int fnum, strand_t strand){
	zTrellis*      trellis = ps->trellis;
	int            poss,i,j,state,prestate;
	score_t        score,base_score;
	zSfeature*     f;
	zIVec*         fmap;
	zIVec*         jumps;
	coor_t         tmp_start;
	zPhase_t       left_phase,right_phase;
	zPinArrayData* to_pad;
	zPinArrayData* from_pad;
	
	poss = ('-' != strand);
	if(poss){
		fmap = trellis->fmap3[fnum];
	}
	else{
		fmap = trellis->rfmap3[fnum];
	}

	f = zSFListMoveFirst(&ps->sfl);
	/* add new live node for each state reachable by each feature in sfl */
	while(f != NULL){
		f->strand = strand;	
		if(poss){
			right_phase = zSFEndPhase(trellis->dna, f);
			left_phase = zSFBeginPhase(trellis->dna, f);
		}
		else{
			tmp_start = f->start;
			f->start  = trellis->dna->length - f->end - 1;
			f->end    = trellis->dna->length - tmp_start - 1;			
			right_phase = zSFBeginPhase(trellis->rdna, f);
			left_phase = zSFEndPhase(trellis->rdna, f);
		}			
		if(f->start != ps->end_pos){
			zDie("Cannot add external features to 3' end of pin search which don't begin at end_pos");
		}
		/* for each (external) state which this feature can end in */
		for(i = 0;i < fmap->size;i++){
			state = fmap->elem[i];
			/* make sure extenral state (state) is compatible with feature */
			if(right_phase != trellis->hmm->state[state].phase) continue;					
			to_pad = zPinSearchGetToPad(ps,f->end,state);
			/* for each from_state select the best path to f->end in state to_pad->state.  
				 This path may go through any of the states which transition to to_pad->state */
			jumps = trellis->hmm->jmap[state];
			for(j = 0;j < jumps->size;j++){
				prestate = jumps->elem[j];
				/* make sure internal state (prestate) is compatible with feature */
				if (!zCompatPhases(poss, trellis->hmm->state[prestate].phase, left_phase)) 
					continue;
				/* for each possible state state/pos, if path to end state though prestate 
					 is the best we've seen so far use it */
				base_score = zScoreExternalState(trellis,state,f->end,f->end-f->start+1,
												 f->score) +
					zGetTransitionScore(trellis->hmm,prestate,state,trellis->tiso_group);
				from_pad = zListMoveFirst(&ps->live_from_list);
				while(from_pad != NULL){
					score = base_score + ps->data[from_pad->index][prestate].score;
					if(score > ps->data[from_pad->index][to_pad->index].score){
						ps->data[from_pad->index][to_pad->index].score = score;
						ps->data[from_pad->index][to_pad->index].pinstate = 
							ps->data[from_pad->index][prestate].pinstate;
					}
					from_pad = zListMoveNext(&ps->live_from_list);
				}
			}
		}
		f = zSFListMoveNext(&ps->sfl);
	}
}

void zPinSearchAddFeatureInit(zPinSearch* ps, int fnum, strand_t strand){
	zTrellis*      trellis = ps->trellis;
	int            poss,i,j,state,prestate;
	score_t        score;
	zSfeature*     f;
	zIVec*         fmap;
	zIVec*         jumps;
	coor_t         tmp_start;
	zPhase_t       left_phase,right_phase;
	zPinArrayData* to_pad;
	zPinArrayData* from_pad;
	
	poss = ('-' != strand);
	if(poss){
		fmap = trellis->fmap3[fnum];
	}
	else{
		fmap = trellis->rfmap3[fnum];
	}

	f = zSFListMoveFirst(&ps->sfl);
	/* add new live to and from nodes for each state reachable by each feature in sfl */
	while(f != NULL){

		zCheckPinSearchListOrder(ps);

		/*fprintf(stderr,"Add %d\t%d -> %d %c\n",fnum,f->start,f->end,f->strand);EVAN*/
		f->strand = strand;	
		if(poss){
			right_phase = zSFEndPhase(trellis->dna, f);
			left_phase = zSFBeginPhase(trellis->dna, f);
		}
		else{
			tmp_start = f->start;
			f->start  = trellis->dna->length - f->end - 1;
			f->end    = trellis->dna->length - tmp_start - 1;			
			right_phase = zSFBeginPhase(trellis->rdna, f);
			left_phase = zSFEndPhase(trellis->rdna, f);
		}			
		/* we need at least one pos before f->start for the (internal) prestate
			 this will never be a feature that would be used anyway and no search 
			 will be successful if started this close to the end of the sequence */
		if(f->start == 0) continue;
		/* for each (external) state which this feature can end in */
		for(i = 0;i < fmap->size;i++){
			state = fmap->elem[i];
			/* make sure extenral state (state) is compatible with feature */
			if(right_phase != trellis->hmm->state[state].phase) continue;					
			to_pad = zPinSearchGetToPad(ps,f->end,state);
			
			/* for each state which transitions to to_pad->state */
			jumps = trellis->hmm->jmap[state];
			for(j = 0;j < jumps->size;j++){
				prestate = jumps->elem[j];
				/* make sure internal state (prestate) is compatible with feature */
				if (!zCompatPhases(poss, trellis->hmm->state[prestate].phase, left_phase)) 
					continue;
				/* for each possible state state/pos, if path to end state though prestate 
				   is the best we've seen so far use it */
				score = zScoreInternalState(trellis,prestate,f->start-1,true)+
					zGetTransitionScore(trellis->hmm,prestate,state,trellis->tiso_group) +
					zScoreExternalState(trellis,state,f->end,f->end-f->start+1,
										f->score);
				from_pad = zPinSearchGetFromPad(ps,f->start-1,prestate);
				if(score > ps->data[from_pad->index][to_pad->index].score){
					ps->data[from_pad->index][to_pad->index].score = score;
					ps->data[from_pad->index][to_pad->index].pinstate = state;
				}			
			}
		}
		f = zSFListMoveNext(&ps->sfl);
	}
}

#endif

