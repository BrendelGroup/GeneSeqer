/*
   ;/hdr/ ************************************************************************
   ;
   ; Copyright (c) 2001 VisualMetrics Corporation
   ; All Rights Reserved.
   ;
   ; This code contribution to the GeneSeqer Product has been contributed
   ; by VisualMetrics Corporation in part by work supported by the National
   ; Institutes of Health under Grant No. R4HG01850B.
   ;
   ; Licensees may use this code module only under the terms of the
   ; GeneSeqer Product as contained in the Product's accompanying
   : License in the README file and ONLY in combination with this Product.
   ; Licensees may NOT extract, modify, or employ this code module for any
   ; purpose except to adapt the Product for the Licensee's own internal
   ; research purposes.  Uses of this code module or any code herein outside
   ; the scope of GeneSeqer and its license ARE STRICTLY PROHIBITED.
   ;
   ; Please direct all communications related to this code module to:
   ;
   ;   Michael Bergman
   ;   VisualMetrics Corporation
   ;   380 Knowling Drive
   ;   Coralville, IA  52241
   ;   U.S.A.
   ;
   ;   Phone:  (319) 339-0110; (319) 339-0665 (fax)
   ;   Email:  mkb@visualmetrics.com
   ;
   ; Code Module Description:
   ; ------------------------
   ;
   ; name of module: TopsPost.c -- The Operation of Posting Sets
   ;
   ; A posting set is one of the fundamental building blocks around which the
   ; Tops system organizes persistent information. A "posting" is simply an
   ; abstract integer value whose semantics are irrelavent to this module.
   ; When a posting set is created it is provided with a comparison function
   ; which imposes on ordering on the individual posting values. The posting
   ; set then creates a list of the posting values ordered by the comparison
   ; function. The primary use of posting sets is to maintain symbol tables;
   ; however, they can be used where-ever lists of information must be stored
   ; and then retrieved by value.
   ;
   ; global symbols defined:
   ;
   ; LONG* TopsPostWrite    Write a new posting to the set
   ; void  TopsPostCreate   Create a posting stream
   ; void  TopsPostDestroy  Destroy a posting set
   ; void  TopsPostFirst    Position posting set on first posting
   ; LONG  TopsPostNext     Get net posting from posting set
   ;
   ; useful notes and assumptions:
   ;
   ; 01/16/00 written by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : as part of the Promula Version 8 TOPS development
   ;
   ;/hdr/ ************************************************************************
 */

#define PRIOFUNC					/* Gives access to PROMULA I/O functions */
#define RAWMFUNC					/* Gives access to raw memory operations */
#define STRGFUNC					/* Gives access to string manipulation */
#include "platform.h"					/* Define the platform hosting this code */
#include "TopsPage.h"					/* Transfer onto a paging system */
#include "TopsPost.h"					/* The operation of posting sets */

#define MAX_LEVELS          20				/* Maximum number of levels in posting set */

#define BLK_FOLB             0				/* Offset of following block number */
#define BLK_NREC             4				/* Offset of number of records in block */
#define BLK_DATA             8				/* Offset of start of data */

typedef struct {
  int PageProvider;					/* The paging system provider number */
  int PageSize;						/* The size of each page */
  int iLevels;						/* Number of levels in current tree */
  int PostingsPerPage;					/* Number of postings per page */
  int PostLinksPerPage;					/* Number of post-link records per page */
  int CurRecNo[MAX_LEVELS];				/* Number of current record within block */
  ULONG CurBlock[MAX_LEVELS];				/* Number of current block */
  ULONG FirstBlk[MAX_LEVELS];				/* Number of first block in block chain */
  ULONG LastBlk[MAX_LEVELS];				/* Number to last block in block chain */
  ULONG MaxRecNo[MAX_LEVELS];				/* Maximum record number within stream */
  UBYTE *PostBlock;					/* Current terminal block pointer */
  ULONG uBlock;						/* Number of current terminal block */
  LONG *Postings;					/* Current postings pointer */
  int nPosting;						/* Number of postings in current block */
  int iPost;						/* Current posting from current block */
  int (*Compare) (					/* Postings ordering comparison function */
#ifdef FPROTOTYPE
    LONG Posting1, LONG Posting2
#endif
  );
} tPostingSet;						/* Posting set control structure */

#define MAX_POST       5				/* Maximum number of posting sets */

static tPostingSet PostingSet[MAX_POST];		/* Actual posting set information */
static UBYTE PostProviders[MAX_POST + 1];		/* Provider control vector */
static tPostingSet *PostSet = NULL;			/* Current provider information */
static int PostSelected = 0;				/* Number of current provider */

typedef struct {
  ULONG Link;						/* Number of linked block */
  LONG Post;						/* Posting value or number of postings */
} tPostLink;						/* Posting control structure type */

/*
   ;/doc/ +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ;
   ; name of routine: LocGetProvider
   ;
   ; Get a new provider for a given Tops resource within the current system.
   ;
   ; calling parameters: None
   ;
   ; return parameters:
   ;
   ; The sequence number of the provider available or a zero if there are no
   ; more providers available.
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */
#ifdef FPROTOTYPE
static int LocGetProvider(void)
#else
static int LocGetProvider()
#endif
{
  auto int nUsed;
  auto int Provider;

  nUsed = PostProviders[0];
  if (nUsed < MAX_POST) {
    nUsed++;
    Provider = nUsed;
    PostProviders[0] = (UBYTE) (nUsed);
    PostProviders[Provider] = 1;
  }
  else {
    for (Provider = 1; Provider <= MAX_POST; Provider++) {
      if (!PostProviders[Provider]) {
	PostProviders[Provider] = 1;
	break;
      }
    }
  }
  if (Provider <= MAX_POST)
    return Provider;
  return 0;
}

/*
   ;/doc/ +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ;
   ; name of routine: LocTerminalSpace
   ;
   ; This function attempts to make space for a new record in a terminal posting
   ; block at a particular location. Note that his function does not actually
   ; insert a record into the block, it merely makes room for that insertion.
   ;
   ; calling parameters:
   ;
   ; tPostLink* hBlock   Points to the posting block to receive the new record.
   ;
   ; int        iRecord  Can have a value between -1 and the number of records in
   ;                     the block.
   ;                     Value  Meaning
   ;                     -----  -------
   ;                      -1    requests that the record be inserted after the
   ;                            last record in the block.
   ;                       0    requests that it be inserted before the first
   ;                            record.
   ;                      +n    requests that it be inserted after the indicated
   ;                            record.
   ;               
   ; return parameters:
   ;
   ; If the insertion succeeeds, then this function returns a pointer to the start
   ; of the storage area where the record is to be placed; else it returns a NULL.
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */
#ifdef FPROTOTYPE
static LONG *LocTerminalSpace(tPostLink *hBlock, int iRecord)
#else
static LONG *LocTerminalSpace(hBlock, iRecord)
  tPostLink *hBlock;
  int iRecord;
#endif
{
  auto int nMove;					/* Number of bytes that have to be moved */
  auto int nRecord;					/* Number of records in block */
  auto LONG *Record;					/* Location to receive new record */

   /*------------------------------------------------------------------------
   ;
   ; Step 1: Determine the amount of space available in the block and if
   ; there is not enough to insert the record, return a NULL to the user.
   ;
   ;------------------------------------------------------------------------*/

  nRecord = hBlock->Post;
  if (nRecord >= PostSet->PostingsPerPage)
    return NULL;

   /*------------------------------------------------------------------------
   ;
   ; Step 2: There is enough space in the block for the record. If the new
   ; record is not to be stored at the end of the block, shift the block
   ; contents over to make room for the record. Note that this is the time
   ; a critical point in the entire posting set algorithm and it is VERY
   ; sensitive to the value of the page size.
   ;
   ;------------------------------------------------------------------------*/

  if (iRecord < 0 || iRecord >= nRecord) {
    Record = (LONG *) (hBlock + 1) + nRecord;
  }
  else {
    Record = (LONG *) (hBlock + 1) + iRecord;
    nMove = (nRecord - iRecord) * sizeof(LONG);
    cpyovl((UBYTE *) (Record), (UBYTE *) (Record + 1), nMove);
  }

   /*-------------------------------------------------------------------------
   ;
   ; Step 3: Increment the blocks record count and return the location to
   ; receive the new record.
   ;
   ;-------------------------------------------------------------------------*/

  hBlock->Post = nRecord + 1;
  return Record;
}

/*
   ;/doc/ +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ;
   ; name of routine: LocNonTerminalSpace
   ;
   ; This function attempts to make space for a new record in a nonterminal
   ; posting block at a particular location. Note that his function does not
   ; actually insert a record into the block, it merely makes room for that
   ; insertion.
   ;
   ; calling parameters:
   ;
   ; tPostLink* hBlock   Points to the posting block to receive the new record.
   ;
   ; int        iRecord  Can have a value between -1 and the number of records in
   ;                     the block.
   ;                     Value  Meaning
   ;                     -----  -------
   ;                      -1    requests that the record be inserted after the
   ;                            last record in the block.
   ;                       0    requests that it be inserted before the first
   ;                            record.
   ;                      +n    requests that it be inserted after the indicated
   ;                            record.
   ;               
   ; return parameters:
   ;
   ; If the insertion succeeeds, then this function returns a pointer to the start
   ; of the storage area where the record is to be placed; else it returns a NULL.
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */
#ifdef FPROTOTYPE
static tPostLink *LocNonTerminalSpace(tPostLink *hBlock, int iRecord)
#else
static tPostLink *LocNonTerminalSpace(hBlock, iRecord)
  tPostLink *hBlock;
  int iRecord;
#endif
{
  auto int nMove;					/* Number of bytes that have to be moved */
  auto int nRecord;					/* Number of records in block */
  auto tPostLink *Record;				/* Location to receive new record */

   /*------------------------------------------------------------------------
   ;
   ; Step 1: Determine the amount of space available in the block and if
   ; there is not enough to insert the record, return a NULL to the user.
   ;
   ;------------------------------------------------------------------------*/

  nRecord = hBlock->Post;
  if (nRecord >= PostSet->PostLinksPerPage)
    return NULL;

   /*------------------------------------------------------------------------
   ;
   ; Step 2: There is enough space in the block for the record. If the new
   ; record is not the be stored at the end of the block, shift the block
   ; contents over to make room for the record. Note that this is the time
   ; critical point in the entire posting set algorithm.
   ;
   ;------------------------------------------------------------------------*/

  if (iRecord < 0 || iRecord >= nRecord) {
    Record = hBlock + nRecord + 1;
  }
  else {
    Record = hBlock + iRecord + 1;
    nMove = (nRecord - iRecord) * sizeof(tPostLink);
    cpyovl((UBYTE *) (Record), (UBYTE *) (Record + 1), nMove);
  }

   /*-------------------------------------------------------------------------
   ;
   ; Step 3: Increment the blocks record count and return the location to
   ; receive the new record.
   ;
   ;-------------------------------------------------------------------------*/

  hBlock->Post = nRecord + 1;
  return Record;
}

/*
   ;/doc/ +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ;
   ; name of routine: LocFindPostingInBlock
   ;
   ; This function performs a binary search within a posting block for a specified
   ; posting. It assumes that the postings already in the block are sorted in
   ; ascending order in accordance with the comparison function supplied for the
   ; overall posting set.
   ;
   ; calling parameters:
   ;
   ; tPostLink* hBlock       Handle to the block to be searched for the posting
   ; LONG       Posting      The value of the posting to be searched for
   ; int        iLevel       The level of the block: a zero indicates that that
   ;                         this is a terminal posting block; while a nonzero
   ;                         indicates that this is a linkage block.
   ; return parameters:
   ;
   ; The function returns the sequence number within the block, relative to one if
   ; a matching posting was found. If the specified posting is below all of the
   ; postings in the block, then a zero is returned. If the posting does not match
   ; any in the block, but is greater than the first posting in the block then the
   ; negative sequence number of the last posting which is below the specified one
   ; is returned.
   ;
   ; useful notes and assumptions:
   ;
   ; Remember that all relations used above are relative to the virtual compare
   ; function supplied when the posting set was defined.
   ;
   ;/doc/ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */
#ifdef FPROTOTYPE
static int LocFindPostingInBlock(tPostLink *hBlock, LONG Posting, int iLevel)
#else
static int LocFindPostingInBlock(hBlock, Posting, iLevel)
  tPostLink *hBlock;
  LONG Posting;
  int iLevel;
#endif
{
  auto int iBottom;					/* Current bottom of search range */
  auto int iMiddle;					/* Current mid point of search range */
  auto int iRelation;					/* Comparison relation between postings */
  auto int iTop;					/* Current top of search range */
  auto int nRecords;					/* Number of records in block */
  auto LONG *Postings;					/* Vector of terminal block postings */
  auto tPostLink *PostLinks;				/* Vector of nonterminal postings + links */

   /*------------------------------------------------------------------------
   ;
   ; Step 1: Determine the number of records in the block, it there are none
   ; then the specified posting is below all of the postings in the block, so
   ; return a relation code of zero.
   ;
   ;------------------------------------------------------------------------*/

  if ((nRecords = hBlock->Post) == 0)
    return 0;

   /*-------------------------------------------------------------------------
   ;
   ; Step 2: Initialize the range of postings which might contain the
   ; specified one so that it encompasses all postings in the block.
   ;
   ;------------------------------------------------------------------------*/

  iBottom = iMiddle = iRelation = 0;
  iTop = nRecords - 1;

   /*------------------------------------------------------------------------
   ;
   ; Step 3: The search for the posting will continue until either the search
   ; range collapses -- i.e. until the bottom of the range exceeds the top --
   ; or a matching posting is encountered. Note that we use different search
   ; loops depending upon whether we are searching a terminal block of
   ; simple postings or a nonterminal block of postings paired with links.
   ;
   ;------------------------------------------------------------------------*/

  if (!iLevel) {
    Postings = (LONG *) (hBlock + 1);
    while (iBottom <= iTop) {
      iMiddle = (iBottom + iTop) / 2;
      iRelation = PostSet->Compare(Posting, Postings[iMiddle]);
      if (!iRelation)
	return iMiddle + 1;
      if (iRelation < 0)
	iTop = iMiddle - 1;
      else
	iBottom = iMiddle + 1;
    }
  }
  else {
    PostLinks = hBlock + 1;
    while (iBottom <= iTop) {
      iMiddle = (iBottom + iTop) / 2;
      iRelation = PostSet->Compare(Posting, PostLinks[iMiddle].Post);
      if (iRelation == 0)
	iRelation = -1;
      if (!iRelation)
	return iMiddle + 1;
      if (iRelation < 0)
	iTop = iMiddle - 1;
      else
	iBottom = iMiddle + 1;
    }
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 4: An exact match for the specified posting does not exist in the
   ; block. The last examined posting is either that posting in the block
   ; is directly greater that the specified one or directly less. Which can
   ; be determined by the value of the relation code. Return the negative
   ; sequence number of the record directly less than the specified one or
   ; zero if the first posting is greater.
   ;
   ;-------------------------------------------------------------------------*/

  if (iRelation > 0)
    return -(iMiddle + 1);
  else if (!iMiddle)
    return 0;
  return -iMiddle;
}

/*
   ;/doc/ +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ;
   ; name of routine: LocFindPosting: Find a Posting in a PostingSet
   ;
   ; This routine locates a specified posting within a posting set. As it moves
   ; through the posting set it leaves behind its traversal path. This path
   ; simplifies writing new postings to the posting set once this routine has
   ; established that they are not yet present.
   ;
   ; As this routine moves through the posting set from its initial node to one
   ; of its terminal nodes, it records in the two members "CurBlock" and
   ; "CurRecNo" information about the "related posting" at that level. In the
   ; member "CurBlock" is recorded the block number of the related posting, while
   ; in "CurRecNo" is recorded the sequence number and relation. This value is
   ; the return value from routine "LocFindPostingInBlock" when applied to that
   ; level.
   ;
   ; calling parameters:
   ;
   ; LONG  Posting   The value of the posting to be located.
   ;
   ; return parameters:
   ;
   ; The function returns a pointer to the information for the posting in the 
   ; lowest block if the specified posting was located; else it returns a NULL.
   ; Regardless of whether or not the posting was located, the posting set
   ; structure is is modified to define the traversal path as described above.
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */
#ifdef FPROTOTYPE
static UBYTE *LocFindPosting(LONG Posting)
#else
static UBYTE *LocFindPosting(Posting)
  LONG Posting;
#endif
{
  auto ULONG iBlock;					/* Number of current block being searched */
  auto int iLevel;					/* Current level being processed */
  auto int iRelated;					/* Number of related record in block */
  auto tPostLink *PostBlock;				/* Points to the block containing postings */

   /*-------------------------------------------------------------------------
   ;
   ; Step 1: If the posting set is currently empty, then set the terminal
   ; related record information to indicate that the specified posting is
   ; below all postings now stored in the posting set and return a NULL.
   ;
   ;------------------------------------------------------------------------*/

  if ((iLevel = PostSet->iLevels - 1) < 0) {
    PostSet->CurBlock[0] = 0;
    PostSet->CurRecNo[0] = 0;
    return NULL;
  }

   /*------------------------------------------------------------------------
   ;
   ; Step 2: All appears well to begin the search for the posting. The initial
   ; block to be searched is the first block of the highest level. Note
   ; that this level never contains more than one block.
   ;
   ;------------------------------------------------------------------------*/

  iBlock = PostSet->FirstBlk[iLevel];

   /*-------------------------------------------------------------------------
   ;
   ; Step 3: Search the current block for the specified and branch according
   ; to the result. If the iRelation specification is nonnegative, then no
   ; more searching is required. In this case branch to step 5.
   ;
   ;------------------------------------------------------------------------*/

  for (;;) {
    PostSet->CurBlock[iLevel] = iBlock;
    PostBlock = (tPostLink *) (TopsPageAccess(iBlock, 0));
    if ((iRelated = LocFindPostingInBlock(PostBlock, Posting, iLevel)) >= 0) {
      break;
    }

      /*-----------------------------------------------------------------------
      ;
      ; Step 4: The specified posting exceeds one of the postings in the
      ; current block. If we are at the terminal node then simply return a
      ; NULL. Else obtain the information needed to process the next level and
      ; loop back to step 3.
      ;
      ;----------------------------------------------------------------------*/

    PostSet->CurRecNo[iLevel] = iRelated;
    if (iLevel == 0) {
      return NULL;
    }
    iBlock = (PostBlock - iRelated)->Link;
    iLevel--;
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 5: The search of the posting set has terminated either because we
   ; have found a matching posting or because our posting is less than any
   ; posting now in the set. If we have a matching posting, continue to step 6;
   ; else branch to step 7.
   ;
   ;-------------------------------------------------------------------------*/

  if (iRelated) {
      /*-----------------------------------------------------------------------
      ;
      ; Step 6: The specified posting matches a posting in the block. If we
      ; are not now at the terminal node, move down to that node using the fact
      ; that we must match the first record in the related block for each 
      ; lower block. Once we have reached the bottom, return the pointer to
      ; the record associated information.
      ;
      ;---------------------------------------------------------------------*/

    if (iLevel == 0) {
      return (UBYTE *) (PostBlock + 1) + (iRelated - 1) * sizeof(LONG);
    }
    PostBlock += (iRelated - 1);
    while (iLevel > 0) {
      iBlock = (PostBlock + 1)->Link;
      iLevel--;
      PostBlock = (tPostLink *) (TopsPageAccess(iBlock, 0));
    }
    return (UBYTE *) (PostBlock + 1);
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 7: The specified posting is below all of the postings in the
   ; indicated block. This can only happen if we are at the initial block,
   ; because we can not get to a lower block unless the posting is at least
   ; equal to the first one in the block. Move down the posting levels setting
   ; each related block and line to zero. Then return a NULL.
   ;
   ;------------------------------------------------------------------------*/

  while (iLevel >= 0) {
    PostSet->CurBlock[iLevel] = PostSet->FirstBlk[iLevel];
    PostSet->CurRecNo[iLevel] = 0;
    iLevel--;
  }
  return NULL;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: TopsPostWrite -- Write New Posting to Set
   ;
   ; This routine writes a new posting to a posting set.
   ;
   ; calling parameters:
   ;
   ; LONG   Posting    The posting to be written.
   ;
   ; return parameters:
   ;
   ; If the posting is successfully added to the set then this fuction returns
   ; a NULL. If the posting is already present in the set then a pointer to the
   ; matching posting value is returned.
   ;
   ; useful notes and assumptions:
   ;
   ; Step 1 below simply says "Make certain that the posting is not already in
   ; the posting set. If it is or if an error occurs return an error indication."
   ; and then it calls LocFindPosting.
   ;
   ; The crutial point is that the call to LocFindPosting no only performs the
   ; above service, but also establishes a traversal path thought the blocks
   ; of the posting set. This traversal path is then used by this function to
   ; add the new posting, It is the concepts of "traversal path" and "related
   ; posting" that set the approach here apart from traditional Btree approach.
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
LONG *TopsPostWrite(LONG Posting)
#else
LONG *TopsPostWrite(Posting)
  LONG Posting;
#endif
{
  auto int iLevels;					/* Number of levels in posting set */
  auto int iNeed;					/* Room needed in Block */
  auto ULONG iNewBlock;					/* Number of new block to receive posting */
  auto int iRecord;					/* Sequence number of record to be updated */
  auto int iRoom;					/* Room available in block */
  auto int iSplit;					/* Record at which a block split is made */
  auto ULONG iSplitBlock;				/* Number of new split block created */
  auto UBYTE *PostBlock;				/* Pointer to current control block */
  auto UBYTE *PostSplit;				/* Pointer to new split block */
  auto int cLevel;					/* Current level being modified */
  auto ULONG fblock;					/* Temporary following block number */
  auto int nRecord;					/* Temporary record count */
  auto LONG NewPost;					/* Value of a new posting to be made */
  auto LONG RepPost;					/* Value of a replacement posting */
  auto ULONG NewLink;					/* Value of new link information */
  auto int nMove;					/* Number of bytes to move on a split */
  auto ULONG RepLink;					/* Value of replacement link */
  auto LONG *PostAddr;					/* Address to receive a new posting */
  auto tPostLink *PostLink;				/* A posting link entry */

   /*--------------------------------------------------------------------------
   ;
   ; Step 1: Make certain that the posting is not already in the posting set.
   ; If it is or if an error occurs return an error indication.
   ;
   ;------------------------------------------------------------------------*/

  if ((PostAddr = (LONG *) (LocFindPosting(Posting))) != NULL) {
    return PostAddr;
  }

   /*-------------------------------------------------------------------------
   ;
   ; Step 2: If the posting set is currently empty. Simply insert the posting,
   ; initialize the terminal stream to contain the first posting, and end local
   ; processing.
   ;
   ;------------------------------------------------------------------------*/

  if ((iLevels = PostSet->iLevels) == 0) {
    PostSet->iLevels = 1;
    PostBlock = TopsPageNew(&iNewBlock);
    PostSet->FirstBlk[0] = PostSet->LastBlk[0] =
      PostSet->CurBlock[0] = iNewBlock;
    PostSet->CurRecNo[0] = PostSet->MaxRecNo[0] = 1;
    *(LONG *) (PostBlock + BLK_NREC) = 1;
    *(ULONG *) (PostBlock + BLK_FOLB) = 0;
    *(LONG *) (PostBlock + BLK_DATA) = Posting;
    return NULL;
  }

   /*-------------------------------------------------------------------------
   ;
   ; Step 3: The posting set is not empty. Attempt to insert the new posting
   ; in the terminal block of the posting set and branch to step 4 or 7
   ; depending upon whether the block overflows.
   ;
   ;------------------------------------------------------------------------*/

  cLevel = 0;
  NewLink = RepLink = 0;
  PostLink = NULL;
  iNewBlock = PostSet->CurBlock[0];
  iRecord = PostSet->CurRecNo[0];
  if (iRecord < 0)
    iRecord = -iRecord;
  PostBlock = TopsPageAccess(iNewBlock, 1);
  PostAddr = LocTerminalSpace((tPostLink *) (PostBlock), iRecord);
  if (PostAddr != NULL) {
      /*-----------------------------------------------------------------------
      ;
      ; Step 4: Room for posting has been successfully allocated in the
      ; terminal level. First increment the number of posting in the level
      ; and then store the posting in it.
      ;
      ;----------------------------------------------------------------------*/

    PostSet->MaxRecNo[0]++;
    *PostAddr = Posting;

      /*----------------------------------------------------------------------
      ;
      ; Step 5: The posting has been entered. If it was not a new initial
      ; entry in a block or if the posting set still consists of single
      ; terminal level then we are are finished -- end local processing.
      ;
      ;---------------------------------------------------------------------*/

    if (iRecord != 0 || iLevels == 1) {
      return NULL;
    }

      /*----------------------------------------------------------------------
      ;
      ; Step 6: The initial entry of a terminal block of the posting set has
      ; changed in a posting set which has at least two levels. Therefore, the
      ; link to this block in the directly higher level must be replaced. In
      ; the terminology of this approach, we have a replacement posting but not
      ; a new posting to be incorporated at the next level. Set the local
      ; varibles to reflect this and branch to move to the next level
      ;
      ;---------------------------------------------------------------------*/

    RepPost = Posting;
    RepLink = iNewBlock;
    NewPost = -1;
  }
  else {
      /*----------------------------------------------------------------------
      ;
      ; Step 7: There is insufficient room in the terminal block which is to
      ; receive the new terminal posting. Split that block in the middle
      ; unless the posting is be inserted as the last entry of the rightmost
      ; block. In this case begin a new rightmost block. Some more work
      ; could be done with picking a better split point, but this seems
      ; adequate for now.
      ;
      ;---------------------------------------------------------------------*/

    PostSplit = TopsPageNew(&iSplitBlock);
    *(LONG *) (PostSplit + BLK_NREC) = 0;
    fblock = *(ULONG *) (PostBlock + BLK_FOLB);
    nRecord = *(LONG *) (PostBlock + BLK_NREC);
    if (fblock != 0 || iRecord < nRecord) {
      iSplit = nRecord / 2;
      nMove = (nRecord - iSplit) * sizeof(LONG);
      cpymem((PostBlock + BLK_DATA) + iSplit * sizeof(LONG),
	(PostSplit + BLK_DATA), nMove);
      *(LONG *) (PostBlock + BLK_NREC) = iSplit;
      *(LONG *) (PostSplit + BLK_NREC) = nRecord - iSplit;
    }
    if (iNewBlock == PostSet->LastBlk[0])
      PostSet->LastBlk[0] = iSplitBlock;
    fblock = *(ULONG *) (PostBlock + BLK_FOLB);
    *(ULONG *) (PostSplit + BLK_FOLB) = fblock;
    *(ULONG *) (PostBlock + BLK_FOLB) = iSplitBlock;

      /*-----------------------------------------------------------------------
      ;
      ; Step 8: The terminal block has been split. Determine which block is
      ; to receive the posting and store it there. Then, since we have created
      ; a new block we have a new posting to be sent to the higher level. We
      ; have a replacement posting, only if there is only one level, or if we
      ; are changing the first entry of the first of the two blocks. Set the
      ; local variables to reflect these conditions and move to the next level.
      ;
      ;---------------------------------------------------------------------*/

    nRecord = *(LONG *) (PostBlock + BLK_NREC);
    if (iRecord < nRecord) {
      PostAddr = LocTerminalSpace((tPostLink *) (PostBlock), iRecord);
    }
    else {
      PostAddr = LocTerminalSpace((tPostLink *) (PostSplit), iRecord - nRecord);
    }
    PostSet->MaxRecNo[0]++;
    *PostAddr = Posting;
    NewPost = *(LONG *) (PostSplit + BLK_DATA);
    NewLink = iSplitBlock;
    if (iRecord == 0 || iLevels == 1) {
      RepPost = *(LONG *) (PostBlock + BLK_DATA);
      RepLink = iNewBlock;
    }
    else
      RepPost = -1;
  }

MoveToNextLevel:   /*---------------------------------------------------------
   ;
   ; Step 9: We are ready to move to the next level in the set. We have
   ; potentially altered two blocks at the current level -- the block which
   ; was there before and a new block immediately after it, if we had to split
   ; the original block. First increment the current level number, and if
   ; we already have a postings at that level, branch to step 11.
   ;
   ;-------------------------------------------------------------------------*/

  cLevel++;
  if (cLevel >= iLevels) {
      /*-----------------------------------------------------------------------
      ;
      ; Step 10: We are starting a new level. We must initialize a new left-
      ; most block, insert the two postings from the previous level into it
      ; and end local processing.
      ;
      ;---------------------------------------------------------------------*/

    PostSet->iLevels++;
    iLevels++;
    PostBlock = TopsPageNew(&iNewBlock);
    PostSet->FirstBlk[cLevel] = PostSet->LastBlk[cLevel] =
      PostSet->CurBlock[cLevel] = iNewBlock;
    PostSet->CurRecNo[cLevel] = PostSet->MaxRecNo[cLevel] = 2;
    *(LONG *) (PostBlock + BLK_NREC) = 2;
    *(ULONG *) (PostBlock + BLK_FOLB) = 0;
    PostLink = (tPostLink *) (PostBlock + BLK_DATA);
    PostLink->Link = RepLink;
    PostLink->Post = RepPost;
    PostLink++;
    PostLink->Link = NewLink;
    PostLink->Post = NewPost;
    return NULL;
  }

   /*-------------------------------------------------------------------------
   ;
   ; Step 11: We are within the levels of a posting set and we are potentially
   ; replacing one posting and possibly adding a second posting immediately
   ; after the one being replaced, or both. The first question is "How much
   ; room is available in the block to be updated?".  This equals the room
   ; available in the block less the length of the posting being replaced if
   ; there is a replacement. So access the block and make the determination
   ; of available space.
   ;
   ;------------------------------------------------------------------------*/

  iNewBlock = PostSet->CurBlock[cLevel];
  iRecord = PostSet->CurRecNo[cLevel];
  if (iRecord < 0)
    iRecord = -iRecord;
  else
    iRecord = 1;
  PostBlock = TopsPageAccess(iNewBlock, 1);
  nRecord = *(LONG *) (PostBlock + BLK_NREC);
  iRoom = (PostSet->PageSize - BLK_DATA) - (nRecord * sizeof(tPostLink));
  if (RepPost >= 0) {
    PostLink = (tPostLink *) (PostBlock + BLK_DATA + (iRecord - 1) * sizeof(tPostLink));
    iRoom += sizeof(tPostLink);
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 12: Now compute the room we need, and branch to step 14 if we do
   ; not have enouph room to do the block modifications.
   ;
   ;------------------------------------------------------------------------*/

  if (NewPost >= 0)
    iNeed = sizeof(LONG) + sizeof(ULONG);
  else
    iNeed = 0;
  if (RepPost >= 0)
    iNeed += (sizeof(LONG) + sizeof(ULONG));
  if (iNeed <= iRoom) {
      /*-----------------------------------------------------------------------
      ;
      ; Step 13: There is enough room in the block to perform the update. Do
      ; it and if the first posting has not been changed, end local processing.
      ; If the first posting does change and if there is a higher level loop
      ; back to step 9 to update that level accordingingly.
      ;
      ;----------------------------------------------------------------------*/

    if (RepPost >= 0) {
      RepLink = PostLink->Link;
      PostLink->Post = RepPost;
    }
    if (NewPost >= 0) {
      PostSet->MaxRecNo[cLevel]++;
      PostLink = LocNonTerminalSpace((tPostLink *) (PostBlock), iRecord);
      PostLink->Link = NewLink;
      PostLink->Post = NewPost;
    }
    if (RepPost >= 0 && iRecord == 1 && (cLevel + 1) < iLevels) {
      NewPost = -1;
      goto MoveToNextLevel;
    }
    return NULL;
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 14: We must split the block inorder to perform the desired update.
   ; Split the block in the middle unless the update is being performed on
   ; the last entry of the rightmost block. In this case begin a new rightmost
   ; block.
   ;
   ;------------------------------------------------------------------------*/

  PostSplit = TopsPageNew(&iSplitBlock);
  *(LONG *) (PostSplit + BLK_NREC) = 0;
  *(ULONG *) (PostSplit + BLK_FOLB) = 0;
  fblock = *(ULONG *) (PostBlock + BLK_FOLB);
  nRecord = *(LONG *) (PostBlock + BLK_NREC);
  if (fblock != 0 || iRecord < nRecord)
    iSplit = nRecord / 2;
  else
    iSplit = iRecord - 1;
  nMove = (nRecord - iSplit) * sizeof(tPostLink);
  cpymem((PostBlock + BLK_DATA) + iSplit * sizeof(tPostLink),
    (PostSplit + BLK_DATA), nMove);
  *(LONG *) (PostBlock + BLK_NREC) = iSplit;
  *(LONG *) (PostSplit + BLK_NREC) = nRecord - iSplit;
  if (iNewBlock == PostSet->LastBlk[cLevel])
    PostSet->LastBlk[cLevel] = iSplitBlock;
  fblock = *(ULONG *) (PostBlock + BLK_FOLB);
  *(ULONG *) (PostSplit + BLK_FOLB) = fblock;
  *(ULONG *) (PostBlock + BLK_FOLB) = iSplitBlock;

   /*--------------------------------------------------------------------------
   ;
   ; Step 15: The update block has been split. Determine which block is to
   ; to receive the update and do it. Then, since we have created a new
   ; block we have a new posting to be sent to the higher level. We have a
   ; replacement posting, only if we are at the top level or if we are changing
   ; the first record of the first of the two blocks. Set the local
   ; variables to reflect these conditions and loop back step 9. Again notice
   ; how much easier life is given the power of semi-permanent pointers.
   ;
   ;------------------------------------------------------------------------*/

  nRecord = *(LONG *) (PostBlock + BLK_NREC);
  if (iRecord < nRecord) {
    if (RepPost >= 0) {
      PostLink->Link = RepLink;
      PostLink->Post = RepPost;
    }
    if (NewPost >= 0) {
      PostSet->MaxRecNo[cLevel]++;
      PostLink = LocNonTerminalSpace((tPostLink *) (PostBlock), iRecord);
      PostLink->Link = NewLink;
      PostLink->Post = NewPost;
    }
  }
  else {
    if (RepPost >= 0) {
      PostLink->Link = RepLink;
      PostLink->Post = RepPost;
    }
    if (NewPost >= 0) {
      PostSet->MaxRecNo[cLevel]++;
      nRecord = *(LONG *) (PostBlock + BLK_NREC);
      PostLink =
	LocNonTerminalSpace((tPostLink *) (PostSplit), iRecord - nRecord);
      PostLink->Link = NewLink;
      PostLink->Post = NewPost;
    }
  }
  NewPost = *(LONG *) (PostSplit + BLK_DATA + sizeof(ULONG));
  NewLink = iSplitBlock;
  if ((iRecord == 1 && RepPost >= 0) || ((cLevel + 1) == iLevels)) {
    RepPost = *(LONG *) (PostBlock + BLK_DATA + sizeof(ULONG));
    RepLink = iNewBlock;
  }
  else
    RepPost = -1;
  goto MoveToNextLevel;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: TopsPostCreate -- Create a Posting to Set
   ;
   ; This routine creates a new posting set.
   ;
   ; calling parameters:
   ;
   ; int   PageProvider   The service provider number for the page management
   ;                      provider to be used.
   ;
   ; int    PageSize      The size of each page.
   ;
   ; int   (*Compare)()   The comparison function to be used to determine the
   ;                      sort relationship between two postings.
   ;
   ; return parameters:  None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
int TopsPostCreate(int PageProvider, int PageSize, int (*Compare) (LONG, LONG))
#else
int TopsPostCreate(PageProvider, PageSize, Compare)
  int PageProvider;
  int PageSize;
  int (*Compare) ();
#endif
{
  auto int iLevel;
  auto int Provider;

   /*--------------------------------------------------------------------------
   ;
   ; Step 1: As a slave server it is up to this function to find an available
   ; provider service number. If one is not available, then return a zero.
   ;
   ;-------------------------------------------------------------------------*/

  Provider = LocGetProvider();
  if (Provider == 0)
    return 0;
  PostSet = PostingSet + (Provider - 1);

  for (iLevel = 0; iLevel < MAX_LEVELS; iLevel++) {
    PostSet->CurRecNo[iLevel] = 0;
    PostSet->CurBlock[iLevel] = 0;
    PostSet->FirstBlk[iLevel] = 0;
    PostSet->LastBlk[iLevel] = 0;
    PostSet->MaxRecNo[iLevel] = 0;
  }
  PostSet->iLevels = 0;
  PostSet->PageSize = PageSize;
  PostSet->PageProvider = PageProvider;
  PostSet->Compare = Compare;
  PostSet->PostingsPerPage = (PageSize - sizeof(tPostLink)) / sizeof(LONG);
  PostSet->PostLinksPerPage = (PageSize - sizeof(tPostLink)) / sizeof(tPostLink);
  return Provider;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: TopsPostDestroy
   ;
   ; This function destroys any information associated with the posting set.
   ;
   ; calling parameters: None
   ;
   ; return parameters: None
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
void TopsPostDestroy(void)
#else
void TopsPostDestroy()
#endif
{
  PostSet = NULL;
  PostSelected = 0;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: TopsPostFirst
   ;
   ; This function positions the posting set on its initial terminal posting
   ; so that its content may be read in sorted order.
   ;
   ; calling parameters: None
   ;
   ; return parameters: None
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
void TopsPostFirst(void)
#else
void TopsPostFirst()
#endif
{
  PostSet->uBlock = PostSet->FirstBlk[0];
  PostSet->PostBlock = TopsPageAccess(PostSet->uBlock, 0);
  PostSet->nPosting = *(LONG *) (PostSet->PostBlock + BLK_NREC);
  PostSet->Postings = (LONG *) (PostSet->PostBlock + BLK_DATA);
  PostSet->iPost = 0;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: TopsPostNext
   ;
   ; This function positions the posting set on its next terminal posting and
   ; returns its value.
   ;
   ; calling parameters: None
   ;
   ; return parameters:
   ;
   ; The value of the next posting or a -1 if there are no more postings.
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
LONG TopsPostNext(void)
#else
LONG TopsPostNext()
#endif
{
  auto LONG Posting;

  if (PostSet->iPost < PostSet->nPosting) {
    Posting = PostSet->Postings[PostSet->iPost];
    PostSet->iPost++;
    return Posting;
  }
  PostSet->uBlock = *(ULONG *) (PostSet->PostBlock + BLK_FOLB);
  if (!PostSet->uBlock)
    return -1;
  PostSet->PostBlock = TopsPageAccess(PostSet->uBlock, 0);
  PostSet->nPosting = *(LONG *) (PostSet->PostBlock + BLK_NREC);
  PostSet->Postings = (LONG *) (PostSet->PostBlock + BLK_DATA);
  PostSet->iPost = 1;
  return PostSet->Postings[0];
}
