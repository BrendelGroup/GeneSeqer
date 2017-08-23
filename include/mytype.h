#ifndef _MYTYPE_H_
#define _MYTYPE_H_

typedef long int ItemType;
typedef char pType;

typedef int IndexType;

#ifdef HP_UX
typedef char bool;
#define true 1
#define false 0
#endif


typedef struct EstHits
  {
    int Start_Est_Match;
    int End_Est_Match;
    int Start_GenomicDNA_Match;
    int End_GenomicDNA_Match;
    IndexType OriginalPos;	
    int quality;
  }
ValueType;

typedef struct
  {
    IndexType index;
    ValueType value;
  }
thePair;
typedef thePair *ptrPair;


/*The following definitions are for result.c*/
#define SMALLGAP                6 /* Specify the small allowed gap between two
				     consecutive matches in the EST sequence. */
#define BIGGAP                 20 /* Specify the big allowed gap between two
				     consecutive matches in the EST sequence. */
#define GDNAJUMP              800 /* Specify the maximum allowed gap bwteen two
				     consecutive matches in the genomic sequence. */
#define NoExtendQuality        70  
#define NoExtendPenalty         4
#define MaxExtension        10000
#define LowQuality             35


/*The following definitions are for GeneSeqer.c */
#define MinWidthOfEST          50
#define MinWidthOfGDNA         50
#define ExdFactor              10
#define ExdLength             600

#endif
