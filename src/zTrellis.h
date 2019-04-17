/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zTrellis.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf
 Modified by Zhiqiang Hu in 2012 
\******************************************************************************/

#ifndef ZOE_TRELLIS_H
#define ZOE_TRELLIS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zDNA.h"
#include "zFeatureFactory.h"
#include "zHMM.h"
#include "zModel.h"
#include "zScanner.h"
#include "zSfeature.h"
#include "zTools.h"
#include "zConseq.h"
#include "zStopSeq.h"
//Added by Jim
#include "zTBTree.h"

#define top_num 100
//struct zTBTreeNode;

//Added by Jim
// struct zTrellisCellLayer
// {
// 	score_t score;
// 	int length;
// 	int trace;
// 	int trace_top;
// 	frame_t frag_frame;
// };
// typedef struct zTrellisCellLayer zTrellisCellLayer;

// struct zTBTreeNode
// {
// 	score_t score;
// 	int length;
// 	int trace;
// 	int trace_top;
// 	frame_t frag_frame;
// };
// typedef struct zTBTreeNode zTBTreeNode;



struct zTrellisCell {
//zTrellisCellLayer **layer;


//	zTrellisCellLayer *layer[top_num];

//	zTBTreeNode *(layer[top_num]);
	zTBTreeNode *(*layer);

	// score_t       score[top_num]; 
	// int           length[top_num];
	// int           trace[top_num];  /* 0 state index */
	// int           trace_top[top_num];/*0<=trace_top<topnum*/
//  struct zTBTreeNode*  tb; 
	zSfeature*    jump;	  

//	frame_t       frag_frame[top_num];   /* keep lfrag, rfrag and frame   */


	/* information in bitwise manner */
	
	/////kevin zSFVec        *submax; /* for states of EXPLICIT and INTERNAL   */
	/* types. Keeps highest-scoring features */ 
	/* for each range of length distribution */
	/* see zExplicitTrans function in        */
	/* zTransition.c for details             */

};
typedef struct zTrellisCell zTrellisCell;

/******************************************************************************\
 zTrellis

zTrellis collects many of the zoe components into a single entity used for the
gene prediciton algorithms. A trellis contains a specific zDNA and zHMM.

\******************************************************************************/

struct zTrellis {
	zDNA             *dna;
	zDNA             *rdna; /* reverse complemented dna */
	zDNA             *unmasked_dna;
	zDNA             *unmasked_rdna; /* reverse complemented dna */
	zConseq          *conseq;
	zConseq          *rconseq;
	bool              top_mode;

#ifdef DEBUG
	bool              show_trellis; /* indicates whether --show_trellis */
	                                /* option has been used             */
#endif

	zHMM             *hmm;
	int               tiso_group; /* isocore group for transitions */
	int               iiso_group; /* isocore group for init. prob. */
	score_t           seq_prob;  /* result of forward alg */
	coor_t            padding;
	
	/* indexed by an hmm state */
	zSFVec          **fexternal; /* external forward links */
	zTrellisCell    **cell;      /* cell[state_idx][seq_pos] */
	score_t         **forward;   /* forward algorithm values */
	score_t         **backward;  /* backward algorithm values */
	zIVec           *extpos;     /* backward links to external       */
                                     /* features scoring above MIN_SCORE */
                                     /* This is used to compute explicit */
                                     /* state features                   */
	zSFList        **external;   /* array of lists of external feature on 
									+ strand ending at current pos */
	zSFList        **rexternal;  /* array of lists of external feature on 
									- strand ending at current pos */

	/* indexed by zStrIdx */
	zScanner        **scanner;   /* map hmm models to scanners here */
	zScanner        **rscanner;  /* scanners for reverse dna */
	zScanner        **conseqscanner;
	zScanner        **rconseqscanner;


	zFeatureFactory **factory;   /* external feature factory */
	zFeatureFactory **rfactory;  /* external feature factory */
	zIVec           **fmap5;
	zIVec           **rfmap5;
	zIVec           **fmap3;
	zIVec           **rfmap3;

	zStopSeq         *stopseq;
	zStopSeq         *rstopseq;
	
};
typedef struct zTrellis zTrellis;

void    zFreeViterbiVars (zTrellis*);
void    zFreeForwardVars (zTrellis*);
void    zFreeBackwardVars (zTrellis*);
void    zFreeBackLinks (zTrellis*);
void    zFreeForwardLinks (zTrellis*);
void    zFreeTrellis (zTrellis*);
void    zInitTrellis (zTrellis*, char*, zHMM*, char*, char*, char*, char*);
zSFVec** zRunViterbiAndForward (zTrellis*,score_t*,zSFVSponge*);
void    zRunPartialViterbiAndForward (zTrellis*, coor_t, coor_t);
void    zTraceTrellis (zTrellis*, int, zSFVec**, score_t*);
void    zTracePartialTrellis (zTrellis*, coor_t, coor_t, int, zSFVec**, score_t*);
void    zRunBackward (zTrellis*);
float   zForwardBackwardProb(zTrellis*, zSfeature*);
void    zShowTrellis(zTrellis*);

score_t zGetIntergenicContinueScore(zTrellis*);
score_t zGetPathScore(zTrellis*, const zSFVec*, zSFVec*);

frame_t zFragFrame2Char(frame_t, frame_t, frame_t);
void zChar2FragFrame(frame_t, zSfeature *);


int zFindTrellisPin(zTrellis*,coor_t);

/*zSFVec* zRunViterbi(zTrellis*, score_t*);*/
zPtrList* zRunSNPViterbi(zTrellis*, score_t*);
zSFVec* zRunPinViterbi(zTrellis*, score_t*, coor_t, coor_t);

//void zInitTrellisCellLayer(zTrellisCellLayer**);
void zNewterbi(zTrellis*, coor_t, int);//Just for testing
//void zTreeterbi(zTrellis*, coor_t, int);
//void zTreeterbiCollect(zTrellis*,coor_t,int,int,int,int);
void zInitOriginList(int **,int,int);
void zShowOrigin(int **, int, int);
//zTrellisCellLayer* zFromCellLayer(zTrellis*, coor_t, int, int);
//zTrellisCellLayer* zToCellLayer(zTrellis*, coor_t, int, int,int);
//Added by Jim new
//zTBTreeNode* zFromCellLayer(zTrellis*, coor_t, int, int);
//zTBTreeNode* zToCellLayer(zTrellis*, coor_t, int, int,int);

void zInitTreeNode(zTBTreeNode**);
void zFreeTreeNode(zTBTreeNode** );

//void zInitTreeNode(void*);
//void zFreeTreeNode(void** );

//void zFreeTrellisCellLayer(zTrellisCellLayer** );


void zHandShake(zTrellis *trellis,int ext_state,coor_t ext_pos,int story);
void zFreeLayer(zTBTreeNode*** p);
#endif

