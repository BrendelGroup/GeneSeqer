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
   ; name of module: result.c
   ;
   ; 07/23/00 changed by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : to use the TOPS database/suffix array algorithm
   ; 10/05/01 changed by: Wei Zhu & Volker Brendel ()
   ; 10/12/01 changed by: Volker Brendel ()
*/

#define DATETIME					/* Gives access to date and time functions */
#define MEMOFUNC					/* Gives access to memory allocation */
#define OSYSFUNC					/* Gives access to the operating system */
#define PRIOFUNC					/* Gives access to PROMULA I/O functions */
#define RAWMFUNC					/* Gives access to raw memory operations */
#define STRGFUNC					/* Gives access to string manipulation */
#include "platform.h"					/* Define the platform hosting this code */
#include "TopsPage.h"					/* The operation of page streams */
#include "TopsPost.h"					/* The operation of postings sets */
#include "TopsSeq.h"					/* Tool Organizing Properties Sequentially */
#include "TopsWire.h"					/* The operation if wiring files to memory */
#include "SufArray.h"
#include "EstData.h"
#include "LcpTree.h"
#include "mytype.h"
#include "list.h"
#include "result.h"
#include "minmax.h"
#include <stdlib.h>


extern FILE *outfp;
extern UBYTE *Tstring;
extern int DataSize;
extern char estName[257];
extern char suffixName[257];
extern char lcptreeName[257];
extern char indexName[257];
extern char dataName[257];
static UBYTE *Pstring;
static int Plength;
static int DEBUGINDEX=-1;
extern int MinQualityHSP;
extern int MinQualityCHAIN;

#define BLOCK_SIZE           4096

#define BUF_SIZE             1024

#define HASH_HITS_TABLE_SIZE 10009

typedef struct {
  int start;
  int end;
  int score;
} SegClass;

static theList HitsTable[HASH_HITS_TABLE_SIZE];
static int HitsTableSize;

static char InxFile[257];
static IndexOfEst *frags;
static long NumOfFrags;
static int InxProvider;
static int ListProvider;
static int PairProvider;
static int ShowTime = 1;

static LONG currEst;
static ValueType newItem;


static int iHash = 0;

#define Min(a,b) ((a>=b)?b:a)
#define Max(a,b) ((a>=b)?a:b)

#ifdef FPROTOTYPE
static int sys_secs(void)
#else
static int sys_secs()
#endif
/* Returns system time in seconds. */
{
  return clock() / CLOCKS_PER_SEC;
}

#ifdef FPROTOTYPE
static int getDebugIndex(void)
#else
static int getDebugIndex()
#endif
{
  char *indexString;
  int index;
  if ((indexString = getenv("GSQINDEX")) != NULL) {
     index = atoi(indexString);
     return index;
  }
  else
     return -1;
}

#ifdef FPROTOTYPE
static void prtPair(thePair* pair)
#else
static void prtPair(pair)
  thePair* pair;
#endif
{
  auto ValueType *item;
  item = &(pair->value);
  fprintf(stderr,"Index: %d (quality %d)\n",pair->index,item->quality);
  fprintf(stderr,"EST match: Start %9d\tEnd %9d\n",item->Start_Est_Match,
  	item->End_Est_Match);
  fprintf(stderr,"Genome match: Start %9d\tEnd %9d\n",
  	item->Start_GenomicDNA_Match, item->End_GenomicDNA_Match);
}

#ifdef FPROTOTYPE
static int getAllowedGapLen(int edelta,int quality)
#else
static int getAllowedGapLen(edelta,quality)
  int edelta, quality;
#endif
/* Returns the maximal allowable gap length. */
{
  int m = 100;

  if (edelta < -LowQuality  ||  edelta+quality < 0)
    return 500;
  if (edelta < 0) edelta= -edelta;
  if (edelta  > LowQuality)
    return   500;
  if (quality > LowQuality)
    return  5000;
  else {
    m+= (LowQuality-edelta);
    return (MIN(m*quality-500,5000));
  }
}

/* used by qsort to sort the  array by their last antidiagonal */
#ifdef FPROTOTYPE
static int compByAntiDiag(ptrPair *pp1,ptrPair *pp2)
#else
static int compByAntiDiag(pp1,pp2)
  ptrPair *pp1, *pp2;
#endif
{
  ptrPair p1=(*pp1);
  ptrPair p2=(*pp2);
  int antid1=(p1->value.End_Est_Match) + (p1->value.End_GenomicDNA_Match);
  int antid2=(p2->value.End_Est_Match) + (p2->value.End_GenomicDNA_Match);
/* fprintf(stderr,"\nantid1=%d \t antid2=%d\n",antid1,antid2); */
  if (antid1 > antid2) {
    return 1;
  }
  else if (antid1 < antid2) {
    return -1;
  }
  /* equal */
  return -1;
}

/*
	To give the penalty for the gaps between p1 and p2.
*/
#define SEGGAPOPEN	2
#define ESTGAPEXT	1
#define MININTRON	50
#define MEDINTRON	500
#define BIGINTRON	5000
#define MININTRONPENAL	5
#define MEDINTRONPENAL	15
#define BIGINTRONPENAL	140

#ifdef FPROTOTYPE
static int gap(ptrPair p1,ptrPair p2)
#else
static int gap(p1,p2)
  ptrPair p1,p2;
#endif
{
  int gGapLen=Max(0,p2->value.Start_GenomicDNA_Match - p1->value.End_GenomicDNA_Match);
	/* the length of a potential mismatch */
  int eGapLen=Max(0,p2->value.Start_Est_Match - p1->value.End_Est_Match);
  int il=gGapLen-eGapLen;
	/* the length of potential intron */
  return SEGGAPOPEN+0.01*eGapLen+0.005*il /* +((il>BIGINTRON)?10:0)*/;
	/* 0.001 is to allow intron as long as 8000 bp */

}


/* To give the score of the largest portion of p2 that has no overlap with p1 */
#ifdef FPROTOTYPE
static int tscore(ptrPair p1,ptrPair p2)
#else
static int tscore(p1,p2)
  ptrPair p1,p2;
#endif
{
  int gOverlapLen=Max(0,-p2->value.Start_GenomicDNA_Match + p1->value.End_GenomicDNA_Match);

/* the length of potential mismatch */
  int eOverlapLen=Max(0,-p2->value.Start_Est_Match + p1->value.End_Est_Match);
  return p2->value.quality-Max(gOverlapLen,eOverlapLen);

}

/* Return 1 if p1 and p2 are close; return 0 other wise */
#define D1      10000
#define D2	20

#ifdef FPROTOTYPE
static int tclose(ptrPair p1,ptrPair p2)
#else
static int tclose(p1,p2)
  ptrPair p1,p2;
#endif
{
  int eDist=p2->value.Start_Est_Match-p1->value.End_Est_Match;
  int gDist=p2->value.Start_GenomicDNA_Match-p1->value.End_GenomicDNA_Match;
  if ( eDist + gDist < D1 && eDist > -D2 && gDist > -D2)
    return 1;
  return 0;
}

#ifdef FPROTOTYPE
int ComparePatterns(LONG Posting1, LONG Posting2)
#else
int ComparePatterns(Posting1, Posting2)
  LONG Posting1;
  LONG Posting2;
#endif
{
  auto int iRelation;

  iRelation = cmpstrn(Pstring + Posting1, Pstring + Posting2, Plength);
  return iRelation;
}

#ifdef FPROTOTYPE
static int HitsTableHashValue(LONG index)
#else
static int HitsTableHashValue(index)
  LONG index;
#endif
{
  return (int) (index % HASH_HITS_TABLE_SIZE);
}

#ifdef FPROTOTYPE
static int HitsTableAddEntry(void)
#else
static int HitsTableAddEntry()
#endif
{
  auto int IsBeginning;
  auto int bucket;
  auto CursorType start;
  auto ptrPair oldPair;
  auto ValueType *oldItem;
  auto thePair *newPair;
  extern int MinMatchLen;

   /*--------------------------------------------------------------------------
   ;
   ; Step 1: Determine whether there now is an entry with the same EST
   ; position.
   ;
   ;-------------------------------------------------------------------------*/

  bucket = HitsTableHashValue(currEst);
  start = NULL;
  ListInitCursor(HitsTable + bucket);
  while (!ListIsOffEnd(HitsTable + bucket) && ListGetCurrItem(HitsTable + bucket)->index != currEst) {
    ListAdvanceCursor(HitsTable + bucket);
  }
  if (!ListIsOffEnd(HitsTable + bucket)) {
    IsBeginning = ListAtBeginning(HitsTable + bucket);
    if (!IsBeginning) {
      ListBackUpCursor(HitsTable + bucket);
      start = ListGetCurrCursor(HitsTable + bucket);
      ListAdvanceCursor(HitsTable + bucket);
    }
    do {
      /*
        ... extend genomic DNA to the right ...
        Start_Genomic_Match should be greater than oldItem, because the search
        is from left to right.
      */
      oldPair = ListGetCurrItem(HitsTable + bucket);
      oldItem = &(oldPair->value);
      /*
        ... extend genomic DNA to the left ...
      */
      if (newItem.End_GenomicDNA_Match <= oldItem->End_GenomicDNA_Match) {
        /*
          ... no extension, nothing to be done ...
          In this case, the matching region usually is a low complexity
          sequence.  Consider lowering the quality value of the hit by
          uncommenting the next statement ...
        */
        /*(oldItem->quality)--;*/
        return 0;
      }
      else if (newItem.Start_GenomicDNA_Match <= oldItem->End_GenomicDNA_Match) {
        /*
          ... some overlap ...
        */
        if (newItem.End_Est_Match <= (oldItem->End_Est_Match + MinMatchLen)) {
          if (newItem.Start_Est_Match > oldItem->Start_Est_Match) {
            /*
              ... new start after original start ...
              => extend
            */
            oldItem->End_GenomicDNA_Match = newItem.End_GenomicDNA_Match;
            if (newItem.End_Est_Match > oldItem->End_Est_Match) {
              oldItem->End_Est_Match = newItem.End_Est_Match;
              (oldItem->quality)++;
            }
            return 1;
          }
          else if (newItem.Start_Est_Match>oldItem->Start_Est_Match-MinMatchLen) {
            /*
              ... new start before original start ...
              Usually, this will be caused by repeats.  Consider lowering the
              quality value of the hit by uncommenting the next statement ...
            */
            /* (oldItem->quality)--; */
            oldItem->Start_Est_Match = newItem.Start_Est_Match;
            oldItem->End_GenomicDNA_Match = newItem.End_GenomicDNA_Match;
            return 1;
          }
        }
      }
      else if (newItem.Start_Est_Match >= oldItem->End_Est_Match) {
        /*
          ... possible intron or exon gap longer than MinMatchLen ...
          SMALLGAP is set as the allowable EST gap between two consecutive
          extensions within GDNAJUMP on the genomic DNA.
        */
        if (newItem.Start_Est_Match < oldItem->End_Est_Match + SMALLGAP) {
          if ((newItem.End_GenomicDNA_Match - oldItem->End_GenomicDNA_Match) < BIGGAP) {
            oldItem->End_GenomicDNA_Match = newItem.End_GenomicDNA_Match;
            oldItem->End_Est_Match = newItem.End_Est_Match;
            oldItem->quality += MinMatchLen;
            return 1;
          }
        }
        else {
          /*
            The gap is greater than SMALLGAP.
            In this case, we extend the matching region only when
            both the gap between the EST matching regions and the gap between
            the genomic DNA matching regions are shorter than BIGGAP.
          */
          if (abs(newItem.Start_Est_Match - oldItem->End_Est_Match) < BIGGAP &&
              abs(newItem.Start_GenomicDNA_Match - oldItem->End_GenomicDNA_Match) < BIGGAP) {
            oldItem->End_GenomicDNA_Match = newItem.End_GenomicDNA_Match;
            oldItem->End_Est_Match = newItem.End_Est_Match;
            oldItem->quality += MinMatchLen;
            return 1;
          }
        }
      }
      /*
        ... other cases; search next hits with the same key ...
      */
      ListAdvanceCursor(HitsTable + bucket);
    } while (!ListIsOffEnd(HitsTable + bucket) && ListGetCurrItem(HitsTable + bucket)->index == currEst);
    /*
      ... nothing found ...
      Add as new hit.
    */
    TopsSeqSelect(PairProvider);
    newPair = (thePair *) (TopsSeqAccess(0));
    TopsSeqSelect(ListProvider);
    if (newPair == NULL) {
      fprintf(stderr, "Insufficient memory for new pair\n");
      exit(1);
    }
    newPair->index = currEst;
    newPair->value = newItem;
    if (IsBeginning) {
      /*
        ... insert at the beginning ...
      */
      ListInsertFirst(HitsTable + bucket, newPair);
    }
    else {
      /*
        ... insert in the middle ...
      */
      ListSetCursor(HitsTable + bucket, start);
      ListInsertAfter(HitsTable + bucket, newPair);
    }
    HitsTableSize++;
    return 1;
  }
  else {
    /*
      ... new ...
    */
    TopsSeqSelect(PairProvider);
    newPair = (thePair *) (TopsSeqAccess(0));
    TopsSeqSelect(ListProvider);
    if (newPair == NULL) {
      fprintf(stderr, "Insufficient memory for new pair\n");
      exit(1);
    }
    newPair->index = currEst;
    newPair->value = newItem;
    ListInsertFirst(HitsTable + bucket, newPair);
    HitsTableSize++;
    return 1;
  }
}

#define SCORE(P) (P->value.quality)

/*
preA:   CurrCursor is pointing to the start position
Action: To chain the segments together.
 */
#ifdef FPROTOTYPE
static void HitsTableExtend(int bucket, CursorType start, int ind)
#else
static void HitsTableExtend(bucket, start, ind)
  int bucket;
  CursorType start;
  int ind;
#endif
{
  auto int number=0;
  auto ptrPair *pairArray;
  auto int *Q, *K;	/* to keep the score and the first segment */
  auto SegClass *classSet;
  auto int classNum;
  auto int i=0,j,newScore;
  auto int isBeginning;
  auto int isFirst,gDist;

  /*
    To get start position
  */
  ListSetCursor(HitsTable + bucket, start);
  isBeginning=ListAtBeginning(HitsTable + bucket);

  /* get the number of hits for the target sequence */
  do {
        number++;
	ListAdvanceCursor(HitsTable+bucket);
  } while ((!ListIsOffEnd(HitsTable + bucket)) &&
           (ListGetCurrItem(HitsTable + bucket)->index == ind));

  ListSetCursor(HitsTable + bucket, start);

  if ( number==1 ) {
    if (ListGetCurrItem(HitsTable + bucket)->value.quality < MinQualityCHAIN) {
      ListDeleteCurrItem(HitsTable+bucket);
      HitsTableSize--;
    }
    return;
  }


/* declare the memory to save those segments */
  pairArray=(ptrPair *)malloc(sizeof(ptrPair)*number);
  Q=(int *)malloc(sizeof(int)*number);
  K=(int *)malloc(sizeof(int)*number);
  classSet=(SegClass *)malloc(sizeof(SegClass)*number);
  if (pairArray==NULL || Q==NULL || K==NULL || classSet==NULL) {
    fprintf(stderr,"Error: memory allocation failed. Exit.\n");
    exit(-1);
  }
  /* remove elements from the list and point them by pairArray */
  number =0;
  do{
	if (ListGetCurrItem(HitsTable + bucket)->value.quality >= MinQualityHSP){
		   /* save only HSPs with score at least MinQualityHSP as specified by
  	  	      the default value or the "-y MinQualityHSP" option
		   */
		pairArray[number]=ListGetCurrItem(HitsTable + bucket);
		number++;
	}
	ListDeleteCurrItem(HitsTable + bucket);
	HitsTableSize--;
  }while ((!ListIsOffEnd(HitsTable + bucket)) &&
	  (ListGetCurrItem(HitsTable + bucket)->index == ind));

  ListSetCursor(HitsTable + bucket, start);
  if(!isBeginning){
	ListBackUpCursor(HitsTable + bucket);
  }

  /* sort the pairArray by their last antidiagonals */
  qsort(pairArray,number,sizeof(ptrPair),(int (*)(const void*, const void*)) compByAntiDiag);

  /* to debug */
  if(ind==DEBUGINDEX){
	fprintf(stderr,"Sorting the pairs by the last anti-diagonal...\n\n");
	for(i=0;i<number;i++) prtPair(pairArray[i]);
  }

  /* Class Set is empty */
  classNum=0;

	for(i=0;i<number;i++){
		/* compute Q(Si) and K(Si); */
		K[i]=i; Q[i]=SCORE(pairArray[i]);
		isFirst=1;
		for(j=i-1;j>=0;j--){
			/* to debug */
			if(ind==DEBUGINDEX){
				fprintf(stderr,"\nComparing ...\n\n");
				prtPair(pairArray[j]);
				prtPair(pairArray[i]);
				fprintf(stderr,"Tscore: %d \n",tscore(pairArray[j],pairArray[i]));
				fprintf(stderr,"Tclose: %d\n",tclose(pairArray[j],pairArray[i]));
				fprintf(stderr,"Gap: %d\n",gap(pairArray[j],pairArray[i]));
				fprintf(stderr,"Q[i]=%d \t Q[j]=%d\n\n",Q[i],Q[j]);
			}
			if(tclose(pairArray[j],pairArray[i]) &&
				tscore(pairArray[j],pairArray[i]) > 10 ){
				newScore=Q[j]+tscore(pairArray[j],pairArray[i])-
							gap(pairArray[j],pairArray[i]);

				gDist=pairArray[i]->value.Start_GenomicDNA_Match
						-pairArray[j]->value.End_GenomicDNA_Match;

				if( !isFirst && gDist > 1000 ){
					/* give penalty if there is a long gap from a segment
					   pair which is not closest one */
					newScore -= (gDist/1000 + 20 );
				if(ind==DEBUGINDEX)
					fprintf(stderr,"\nPenalty for long gap (not closest)= %d\n",-(gDist/1000 + 20));
				}
				if(ind==DEBUGINDEX)
					fprintf(stderr,"\nNew score=%d\n",newScore);
				if(newScore>Q[i]){
					Q[i]=newScore;
					K[i]=K[j];
					if(isFirst) isFirst = 0;
					if(ind==DEBUGINDEX)
						fprintf(stderr,"New start=%d\n\n",K[j]);
				}
			}
		} /* compute Q[i] */
		if(ind==DEBUGINDEX)
			fprintf(stderr,"\nWorking with segment %d with score %d.\n",
				i,Q[i]);
		if(Q[i] >= MinQualityCHAIN){

			/* ... is there a class C with start(c)==K(Si)? */
			for(j=0;j<classNum;j++){
				if(classSet[j].start==K[i]) break;
			}
			if(j<classNum){ /* true */
				if(classSet[j].score<Q[i]){
					/* update */
					if(ind==DEBUGINDEX)
						fprintf(stderr,"\nUpdate the existing class.\n\n");
					classSet[j].end=i;
					classSet[j].score=Q[i];
				}else{
					/* ... a better class with the same starting segment already exists */
					if(ind==DEBUGINDEX)
						fprintf(stderr,"\nA better class already exists with score %d starting at %d.\n\n",classSet[j].score,classSet[j].start);

					if(SCORE(pairArray[i]) >= MinQualityCHAIN){
						/* also add as a new class */
						classSet[classNum].start=i;
						classSet[classNum].end=i;
						classSet[classNum].score=SCORE(pairArray[i]);
						classNum++;

					}
				}
			}else{	/* false */
				/* create a new class C with */
				if(ind==DEBUGINDEX)
					fprintf(stderr,"\nAdd as a new class.\n\n");
				classSet[classNum].start=K[i];
				classSet[classNum].end=i;
				classSet[classNum].score=Q[i];
				classNum++;

			}
		} /* process a chain of score > minQuality */
	}/* each i */

 /* Add classSet information back to hashtable */
	if(ind==DEBUGINDEX)
		fprintf(stderr,"\nAfter extending...\n\n");
	for(i=0;i<classNum;i++){
		pairArray[classSet[i].start]->value.End_Est_Match =
			pairArray[classSet[i].end]->value.End_Est_Match;
		pairArray[classSet[i].start]->value.End_GenomicDNA_Match =
			pairArray[classSet[i].end]->value.End_GenomicDNA_Match;
		pairArray[classSet[i].start]->value.quality=classSet[i].score;
		if(ind==DEBUGINDEX) prtPair(pairArray[classSet[i].start]);
		if(isBeginning){
			ListInsertFirst(HitsTable + bucket, pairArray[classSet[i].start]);

		}else{
			ListInsertAfter(HitsTable + bucket, pairArray[classSet[i].start]);
		}
		HitsTableSize++;
	} /* each class */

/* free the memory */
   free(pairArray);
   free(Q);
   free(K);
   free(classSet);
}

#ifdef FPROTOTYPE
static void HitsTableAssembly(void)
#else
static void HitsTableAssembly()
#endif
{
  auto CursorType start;
  auto int ind;
  auto thePair *oldPair;
  auto int i;

  for (i = 0; i < HASH_HITS_TABLE_SIZE; i++) {
    ListInitCursor(HitsTable + i);
    while (!ListIsOffEnd(HitsTable + i)) {
      start = ListGetCurrCursor(HitsTable + i);
      oldPair = ListGetCurrItem(HitsTable + i);
      ind = oldPair->index;
      HitsTableExtend(i, start, ind);
      ListSetCursor(HitsTable + i, start);
      do {
	ListAdvanceCursor(HitsTable + i);
      }
      while ((!ListIsOffEnd(HitsTable + i)) && (ListGetCurrItem(HitsTable + i)->index == ind));
    }
  }
}

#ifdef FPROTOTYPE
void ResultCreate(char *indexfilename)
#else
void ResultCreate(indexfilename)
  char *indexfilename;
#endif
{
  auto int Length;

  Length = strlen(indexfilename);
  cpymem(indexfilename, InxFile, Length + 1);
}

#ifdef FPROTOTYPE
void ResultDestroy(void)
#else
void ResultDestroy()
#endif
{
  TopsSeqSelect(ListProvider);
  TopsSeqDestroy();
  TopsSeqSelect(PairProvider);
  TopsSeqDestroy();
}

/*
   binary search to find the index's location
 */
#ifdef FPROTOTYPE
static int ResultCheck(LONG index, int length, int PosInPattern)
#else
static int ResultCheck(index, length, PosInPattern)
  LONG index;
  int length;
  int PosInPattern;
#endif
{
  auto long low;
  auto long mid;
  auto long high;
  auto long minus;
  auto int seqlength;
  auto int currPos;
  auto int IsFound;

  IsFound = 0;
  low = 0;
  high = NumOfFrags;

  while ((high - low) > 1 && !IsFound) {
    mid = (low + high) / 2;
    minus = index - frags[mid].PosInFile;
    if (minus == 0) {
      low = mid;
      IsFound = 1;
    }
    else if (minus > 0) {
      low = mid;
    }
    else {
      high = mid;
    }
  }

  seqlength = frags[low].LenOfSeq;
  if (index - frags[low].PosInFile + length > seqlength) {
    return 0;
  }

  currEst = frags[low].giIndex;
  /*
     0 is the first position
   */
  currPos = index - frags[low].PosInFile + frags[low].PosInSeq;
  newItem.Start_Est_Match = currPos;
  newItem.End_Est_Match = currPos + length - 1;
  newItem.Start_GenomicDNA_Match = PosInPattern;
  newItem.End_GenomicDNA_Match = PosInPattern + length - 1;
  newItem.OriginalPos = frags[low].OriginalPos;
  newItem.quality = length;
  return 1;
}

#ifdef FPROTOTYPE
static int ResultFilter(char *seq, int len)
#else
static int ResultFilter(seq, len)
  char *seq;
  int len;
#endif
{
  auto int i;
  auto int have[4];
  auto int count;

  count = 0;
  for (i = 0; i < 4; i++)
    have[i] = 0;
  for (i = 0; i < len; i++) {
    if (seq[i] < 4) {
      if (!have[(int)seq[i]]) {
	have[(int)seq[i]] = 1;
	count++;
	if (count > 2)
	  return 1;
      }
    }
  }
  return 0;
}

#ifdef FPROTOTYPE
static void LocDbSearch2(char *pattern, int length, int matchsize)
#else
static void LocDbSearch2(pattern, length, matchsize)
  char *pattern;
  int length;
  int matchsize;
#endif
{
  auto int len;
  auto int i;
  auto LONG Lw;
  auto LONG Rw;
  auto LONG Buffer[BUF_SIZE];
  auto binfile SufFile;
  auto LONG *LwRw;
  auto int PageProvider;
  auto int iPat;
  auto int iSuffix;
  auto int iRelation;
  auto LONG *Posting;
  auto int iMatch;
  auto int nMatch;
  auto LONG *Suffixs;

   /*-------------------------------------------------------------------------
   ;
   ; Step 1: Using a posting set, sort the input pattern. Note that we
   ; allocate an LwRw array with three entries for each pattern position.
   ;
   ;-------------------------------------------------------------------------*/

  PageProvider = TopsPageMemory(BUF_SIZE);
  TopsPageSelect(PageProvider);
  TopsPostCreate(PageProvider, BUF_SIZE, ComparePatterns);
  len = length - matchsize + 1;
  Pstring = (UBYTE *) (pattern);
  Plength = matchsize;
  LwRw = (LONG *) (getmem(len * sizeof(LONG) * 3));
  for (i = iPat = 0; i < len; i++, iPat += 3) {
    LwRw[iPat] = LwRw[iPat + 1] = LwRw[iPat + 2] = -1;
    if (!ResultFilter(pattern + i, matchsize))
      continue;
    if ((Posting = TopsPostWrite(i)) != NULL) {
      LwRw[iPat] = -2;
      LwRw[iPat + 1] = *Posting;
    }
  }

   /*-------------------------------------------------------------------------
   ;
   ; Step 2: Scan the entries in the suffix array sequentially from start
   ; to finish to establish the locations of the pattern components.
   ;
   ;-------------------------------------------------------------------------*/

  SufFile = rdobinf(suffixName);
  iSuffix = BUF_SIZE;
  TopsPostFirst();
  Lw = -1;
  for (;;) {
    if ((i = TopsPostNext()) < 0)
      goto FoundAllPatterns;
    iPat = i * 3;
    for (;;) {
      if (iSuffix == BUF_SIZE) {
	rdbinf(SufFile, Buffer, sizeof(Buffer));
	iSuffix = 0;
      }
      Lw++;
      iRelation = cmpstrn(Tstring + Buffer[iSuffix], pattern + i, matchsize);
      iSuffix++;
      if (iRelation >= 0)
	break;
      if (Lw > (DataSize - matchsize))
	goto FoundAllPatterns;
    }
    if (iRelation > 0)
    {
       if (Lw > (DataSize - matchsize)) goto FoundAllPatterns;
       continue;
    }
    LwRw[iPat] = Lw;
    for (;;) {
      if (iSuffix == BUF_SIZE) {
	rdbinf(SufFile, Buffer, sizeof(Buffer));
	iSuffix = 0;
      }
      Lw++;
      iRelation = cmpstrn(Tstring + Buffer[iSuffix], pattern + i, matchsize);
      iSuffix++;
      if (iRelation != 0)
	break;
    }
    Lw--;
    iSuffix--;
    LwRw[iPat + 1] = Lw;
  }

FoundAllPatterns:  /*--------------------------------------------------------
   ;
   ; Step 3: We have found all of out matchs. We can destroy the binary string.
   ; Then we must determine the starting position in the matchs vector of each
   ; pattern's matching suffixs.
   ;
   ;-------------------------------------------------------------------------*/

  fprintf(outfp,"\nUseDB: Found all patterns, elapsed seconds = %d\n", sys_secs() - ShowTime);
  EstDestroy();
  nMatch = 0;
  for (i = iPat = 0; i < len; i++, iPat += 3) {
    Lw = LwRw[iPat];
    if (Lw >= 0) {
      Rw = LwRw[iPat + 1];
      LwRw[iPat + 2] = nMatch;
      nMatch += (Rw - Lw + 1);
    }
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 4: We know how many matchs we have and where to put them. We must
   ; now initialize a matches vector and copy the actual suffix values into it.
   ;
   ;-------------------------------------------------------------------------*/

  Suffixs = (LONG *) (getmem(nMatch * sizeof(LONG)));
  TopsPostFirst();
  for (;;) {
    if ((i = TopsPostNext()) < 0)
      break;
    iPat = i * 3;
    Lw = LwRw[iPat];
    if (Lw >= 0) {
      Rw = LwRw[iPat + 1];
      iMatch = LwRw[iPat + 2];
      putbinf(SufFile, Lw * sizeof(LONG));
      rdbinf(SufFile, Suffixs + iMatch, (Rw - Lw + 1) * sizeof(LONG));
    }
  }
  fprintf(outfp,"\nUseDB: Found all matches, elapsed seconds = %d\n", sys_secs() - ShowTime);
  fflush(stdout);
  clsbinf(SufFile);
  TopsPostDestroy();
  TopsPageDestroy();

   /*--------------------------------------------------------------------------
   ;
   ; Step 5: Wire the index file into memory and initialize the property
   ; streams needed to contain the pairs and the list nodes.
   ;
   ;-------------------------------------------------------------------------*/

  InxProvider = TopsWireOpen(InxFile);
  if (InxProvider == 0) {
    fprintf(stderr,"\nUnable to open Index File %s\n", InxFile);
    exit(1);
  }
  NumOfFrags = TopsWireSize(InxProvider) / sizeof(IndexOfEst);
  frags = (IndexOfEst *) (TopsWireAddress(InxProvider));
  PairProvider = TopsSeqCreate(BLOCK_SIZE, sizeof(thePair));
  ListProvider = TopsSeqCreate(BLOCK_SIZE, sizeof(ListNode));
  TopsSeqSelect(ListProvider);

   /*------------------------------------------------------------------------
   ;
   ; Step 6: Move through the matches found above and index them using the
   ; Est data index file.
   ;
   ;------------------------------------------------------------------------*/

  for (i = iPat = 0; i < len; i++, iPat += 3) {
    Lw = LwRw[iPat];
    Rw = LwRw[iPat + 1];
    iMatch = LwRw[iPat + 2];
    if (Lw < 0) {
      if (Lw == -1)
	continue;
      Rw *= 3;
      Lw = LwRw[Rw];
      iMatch = LwRw[Rw + 2];
      Rw = LwRw[Rw + 1];
    }
    while (Lw <= Rw) {
      if (ResultCheck(Suffixs[iMatch], matchsize, i))
	HitsTableAddEntry();
      Lw++;
      iMatch++;
    }
  }
  free(LwRw);
  free(Suffixs);
  TopsWireClose(InxProvider);
  frags = NULL;
}

#ifdef FPROTOTYPE
int ResultSearch2(char *pattern, int length, int matchsize)
#else
int ResultSearch2(pattern, length, matchsize)
  char *pattern;
  int length;
  int matchsize;
#endif
{
  auto int len;
  auto int i;
  auto LONG Lw;
  auto LONG Rw;
  auto int nMatch;
  auto int iMatch;
  auto int Scratch;
  auto LONG LwRw[2];
  auto int *Buffer;
  auto int UseDb = 0;

   /*--------------------------------------------------------------------------
   ;
   ; Step 1: Initialize the hits table.
   ;
   ;-------------------------------------------------------------------------*/

  ShowTime = sys_secs();
  HitsTableSize = 0;
  filmem(HitsTable, sizeof(HitsTable), '\0');

  DEBUGINDEX = getDebugIndex();
   /*--------------------------------------------------------------------------
   ;
   ; Step 2: Read in the binary EST datafile. We need this for all approaches.
   ;
   ;-------------------------------------------------------------------------*/

  DataSize = EstReadDatFile(dataName, &Tstring);

   /*-------------------------------------------------------------------------
   ;
   ; Step 3: If we are using the database approach do so and then branch to
   ; assemble the final result.
   ;
   ;-------------------------------------------------------------------------*/

  if (UseDb) {
    LocDbSearch2(pattern, length, matchsize);
    goto AssembleResult;
  }

   /*-------------------------------------------------------------------------
   ;
   ; Step 4: We are using binary search approach. Restore the entire suffix
   ; array and the longest common prefix tree and then move through the
   ; pattern, find all the matches and write them to a scratch file.
   ;
   ;------------------------------------------------------------------------*/

  SuffixArrayCreate(DataSize);
  SuffixArrayRestore(suffixName, lcptreeName, DataSize);
  Scratch = TopsSeqCreate(BLOCK_SIZE, 2 * sizeof(int));
  nMatch = 0;
  len = length;
  for (i = 0; i <= len - matchsize; i++) {
    if (!ResultFilter(pattern + i, matchsize))
      continue;
    if (SuffixArraySearchLwRw(pattern + i, matchsize, LwRw) != 0) {
      Lw = LwRw[0];
      Rw = LwRw[1];
      while (Lw <= Rw) {
	nMatch++;
	Buffer = (int *) (TopsSeqAccess(0));
	*Buffer = SuffixArray[Lw];
	*(Buffer + 1) = i;
	Lw++;
      }
    }
  }

  fprintf(outfp,"\n... found all matches, elapsed seconds = %d\n", sys_secs() - ShowTime);

   /*--------------------------------------------------------------------------
   ;
   ; Step 5: Release the binary EST datafile, the suffix array and the Lcp tree
   ; since they are no longer needed.
   ;
   ;-------------------------------------------------------------------------*/

  SuffixArrayDestroy();
  EstDestroy();

   /*-------------------------------------------------------------------------
   ;
   ; Step 5: Wire the index file into memory and initialize the property
   ; streams needed to contain the pairs and the list nodes.
   ;
   ;-------------------------------------------------------------------------*/

  InxProvider = TopsWireOpen(InxFile);
  if (InxProvider == 0) {
    fprintf(stderr,"\nUnable to open Index File %s\n", InxFile);
    exit(1);
  }
  NumOfFrags = TopsWireSize(InxProvider) / sizeof(IndexOfEst);
  frags = (IndexOfEst *) (TopsWireAddress(InxProvider));
  PairProvider = TopsSeqCreate(BLOCK_SIZE, sizeof(thePair));
  ListProvider = TopsSeqCreate(BLOCK_SIZE, sizeof(ListNode));

   /*------------------------------------------------------------------------
   ;
   ; Step 6: Move through the matches found above and index them using the
   ; EST data index file.
   ;
   ;------------------------------------------------------------------------*/

  for (iMatch = 1; iMatch <= nMatch; iMatch++) {
    TopsSeqSelect(Scratch);
    Buffer = (int *) (TopsSeqAccess(iMatch));
    TopsSeqSelect(ListProvider);
    Lw = *Buffer;
    i = *(Buffer + 1);
    if (ResultCheck(Lw, matchsize, i))
      HitsTableAddEntry();
  }
  TopsWireClose(InxProvider);
  TopsSeqSelect(Scratch);
  TopsSeqDestroy();
  TopsSeqSelect(ListProvider);
  frags = NULL;

AssembleResult: /*----------------------------------------------------------
   ;
   ; Step 7: Assemble the final result set.
   ;
   ;-------------------------------------------------------------------------*/

  HitsTableAssembly();

  fprintf(outfp,"... matches indexed, elapsed seconds = %d HitsTableSize = %d\n", sys_secs() - ShowTime, HitsTableSize);
  return HitsTableSize;
}

#ifdef FPROTOTYPE
void ResultFirst(void)
#else
void ResultFirst()
#endif
{
  TopsSeqSelect(ListProvider);
  iHash = 0;
  ListInitCursor(HitsTable + iHash);
}

#ifdef FPROTOTYPE
thePair *ResultNext(void)
#else
thePair *ResultNext()
#endif
{
  auto thePair *Result;

  while (ListIsOffEnd(HitsTable + iHash)) {
    iHash++;
    if (iHash >= HASH_HITS_TABLE_SIZE)
      return NULL;
    ListInitCursor(HitsTable + iHash);
  }
  Result = ListGetCurrItem(HitsTable + iHash);
  ListAdvanceCursor(HitsTable + iHash);
  return Result;
}
