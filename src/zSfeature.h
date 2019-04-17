/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zSfeature.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_SFEATURE_H
#define ZOE_SFEATURE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "zTools.h"
#include "zDNA.h"
//Added by Jim
// #include "zTrellis.h"
// #include "zTransition.h"

/******************************************************************************\
 zSfeature

zSfeature is for sequence features like splice sites or exons. zSfeatureName is
a very important enumerated type.

The zSFVec copies value, so you can pass a pointer to a zSfeature in and it
will allocate storage for it. It is very important that you make sure that if
the group is unspecified, you give it a NULL value. You don't have to do this
for zReadSfeature, but if you're creating zSfeatures on the fly, be aware that
the zCopySfeature command will segfault if the group points randomly.

\******************************************************************************/

struct zSfeature {
	zStrIdx       name;
	int           state;
	coor_t        start;
	coor_t        end;
	strand_t      strand;
	score_t       score;
	frame_t       lfrag;
	frame_t       rfrag;
	frame_t       frame;
	char         *group;

#ifdef DEBUG
    score_t       begin_score;
    score_t       end_score;
    score_t       content_score;
#endif /* DEBUG */
	score_t       cdscons_score;
	score_t       cdsest_score;
	score_t       align_score;

	/* EVAN the following are Mikhails additions */
	score_t       intrinsic_score;
	int           from_state;
	int           trace;
};
typedef struct zSfeature zSfeature;

struct zSFVec {
	zSfeature *elem;   /* array of features */
	int        size;   /* number of features */
	int        limit;  /* number of elements currently allocated */
	zSfeature *last;
};
typedef struct zSFVec zSFVec;

/**********************************************************************\
  zSFList - list of zSfeature objects 
\**********************************************************************/

struct zSFListNode;
typedef struct zSFListNode zSFListNode;
struct zSFListNode{
	zSfeature    data;
	zSFListNode* next;
	zSFListNode* prev;
};

struct zSFList {
	zSFListNode* head;
	zSFListNode* tail;
	zSFListNode* current;
	int          size; 
};
typedef struct zSFList zSFList;

//added by Jim
/**********************************************************************\
  zSFVSponge - sponge structure of zSFV objects 
\**********************************************************************/

struct zSFVSpongeNode;
typedef struct zSFVSpongeNode zSFVSpongeNode;
struct zSFVSpongeNode{
	zSFVec** node;
	coor_t start;
	coor_t end;
	int end_state;
	int pre_state;
	int size;
	int broken;
	zSFVSpongeNode* prev;
	zSFVSpongeNode* next;
	int location;
	score_t* score_box;
};

struct zSFVSponge{
	zSFVSpongeNode* head;
	zSFVSpongeNode* tail;
	zSFVSpongeNode* current;
	int size;
};
typedef struct zSFVSponge zSFVSponge;

void zInitSFVSponge(zSFVSponge*);
void zFillSFVSponge(zSFVSponge*,zSFVec**,int);//,zTrellis*);
void zFreeSFVSponge(zSFVSponge*);
zSFVSpongeNode* zSFVSpongeMoveFirst(zSFVSponge*);
zSFVSpongeNode* zSFVSpongeMoveLast(zSFVSponge*);
zSFVSpongeNode* zSFVSpongeMoveNext(zSFVSponge*);
zSFVSpongeNode* zSFVSpongeMovePrev(zSFVSponge*);
zSFVSpongeNode* zSFVSpongeGetCurrent(zSFVSponge*);

//score_t zScoreSFVec(zTrellis *trellis,zSFVec *vec);

void       zInitSFList(zSFList*);
void       zFreeSFList(zSFList*);
void       zResetSFList(zSFList*);


zSfeature* zSFListMoveNext(zSFList*);
zSfeature* zSFListMovePrev(zSFList*);
zSfeature* zSFListMoveFirst(zSFList*);
zSfeature* zSFListMoveLast(zSFList*);

zSfeature* zSFListGetCurrent(zSFList*);
zSfeature* zSFListInsert(zSFList*,zSfeature*);
zSfeature* zSFListInsertNext(zSFList*,zSfeature*);
zSfeature* zSFListInsertPrev(zSFList*,zSfeature*);
zSfeature* zSFListAppend(zSFList* l);
zSfeature* zSFListPrepend(zSFList* l);
void       zSFListRemoveFirst(zSFList*);
void       zSFListRemoveLast(zSFList*);
void       zSFListRemoveCurrent(zSFList*);
void       zSFList2SFVec(zSFList*,zSFVec*);
//Added by Jim
void zPrintSFList(zSFList*);

/**********************************************************************\
  zSfeature functions
\**********************************************************************/

void zClearSfeature (zSfeature*);
void zFreeSfeature (zSfeature*);
void zWriteSfeature (FILE*, const zSfeature*);
int  zReadSfeature (FILE*, zSfeature*);
void zCopySfeature (const zSfeature*, zSfeature*);
void zAntiSfeature (zSfeature*, coor_t length);
int  zVerifySfeature (const zSfeature*);
int  zSfeatureCmp (const zSfeature*, const zSfeature*);
int  zSfPtrCmp(const void*, const void*);
//Added by Jim
void zInitSfeature (zSfeature*);

void    zInitSFVec (zSFVec*, int);
void    zPushSFVec (zSFVec*, const zSfeature*);
void    zResetSFVec (zSFVec*);
void    zFreeSFVec (zSFVec*);
void    zFlipSFVec (zSFVec*);
void    zWriteSFVec (FILE*, const zSFVec*);
zSFVec* zReadSFVec  (FILE*);
void    zTranslateSFVec (zSFVec*, int);
void    zTranslateSFList (zSFList*, int);
//Added by Jim
void zPrintSfeature(zSfeature *);
void zMergeSFVec (zSFVec *vec, zSFVec *m_vec);
void zConcatenateSFVec (zSFVec *vec, zSFVec *m_vec);





zPhase_t zSFEndPhase(zDNA*, zSfeature*);
zPhase_t zSFBeginPhase(zDNA*, zSfeature*);
int zCompatPhases(int positive_strand, zPhase_t from, zPhase_t to);

#endif
