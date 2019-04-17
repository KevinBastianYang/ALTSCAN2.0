#include "zTools.h"
#include "zSfeature.h"

/**********************************************************************\
  Traceback Tree

  This is a tree structure which can store all live traceback paths 
  from a trellis more efficiently than retaining the full trellis.

\**********************************************************************/
#ifndef TBT_NODE
#define TBT_NODE
struct zTBTreeNode;
typedef struct zTBTreeNode zTBTreeNode;



struct zTBTree {
	zTBTreeNode*  root;   /* this is the tree root */
	zTBTreeNode*  cpoint; /* the cpoint of this tree */
	zTBTreeNode*  dead_nodes;  /* this faux tree holds the dead nodes stored for reuse
								  the nodes are stored in a single linked list using the 
								  child pointer */
	int           size;   /* total number of nodes in tree */
	short         new_cpoint; /*flag set to 1 when cpoint changes */
	short         id;
	short         check_count;
};
typedef struct zTBTree zTBTree;

struct zTBTreeNode {
	
	
	coor_t           cdna_pos;
	
//	score_t          score;
	strand_t         strand;
	frame_t          frame_data;
	zPhase_t         phase;
	zTBTreeNode*     parent;
	zTBTreeNode*     child;
	zTBTreeNode*     lsib;
	zTBTreeNode*     rsib;
	int              children;
	int              id;
	int              lock;
	short            alive;/*EVAN*/
	zTBTree*         tree;


	//Added by Jim
	score_t score;
	int length;
	int trace;
	int trace_top;
	frame_t frag_frame;

	coor_t	pos;
	zSFVec* sfv;
	int state;
	int story;//It's like the layer of this node

	coor_t left_pos;//only for chunk the smallest pos of a root
	/*EVAN
	score_t          dec_score; / score of all decendants for heuristic /
	bool             sim_trim;
	int              ng[11];*/
};
/* each node sees its leftmost child only.  other accessed through the sibling
   pointers.  each child knows its parent. */
#endif

void zInitTBTree(zTBTree *t);
void zFreeTBTree(zTBTree *t);
void zResetTBTree(zTBTree* tree);
bool zTBTreeCheckNewCpoint(zTBTree* tree);
bool zTBTreeClearNewCpoint(zTBTree* tree);
void zReleaseDeadTBTreeNode(zTBTree* t, zTBTreeNode* n);
void zTBTreeLockNode(zTBTreeNode* n);
void zTBTreeClearSubTree(zTBTree* t, zTBTreeNode* n);
void zCopyTBTreeNode(zTBTree* t1, zTBTreeNode* n1, zTBTree* t2, zTBTreeNode* n2, int depth);
void zAppendTBTreeFromNode(zTBTree* t1, zTBTreeNode* n1, zTBTree* t2, zTBTreeNode* n2);
void zCopyTBTreeFromCPoint(zTBTree* t1,zTBTree* t2,coor_t min_pos);
zTBTreeNode* zGetTBTreeNode(zTBTree *t);
void zTBTreeSetChild(zTBTreeNode* p,zTBTreeNode* c);
void zClearTBTreeNodeChildren(zTBTree* t, zTBTreeNode* n);
void zReleaseRedundantTBTreeNode(zTBTree* t,zTBTreeNode* n);

void* zCreateTBTreeNode();
void zFreeTBTreeNode(void* );
void zResetTBTreeNode(void* );
//Added by Jim
void zShowTBTree(zTBTree*);
void zQueueIn(zTBTreeNode**,zTBTreeNode*);
zTBTreeNode* zQueueOut(zTBTreeNode**);

bool zCheckTBTreeNodeChildren(zTBTreeNode*);
void zDeleteTBTreeNode(zTBTreeNode* n);