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
   ; name of module: DbSuffix.c - Database Suffix Array Build
   ;
   ; The techique of building a suffix array using both a position and a lexical
   ; ordering vector and a binary tree all in memory at the same time works
   ; well until the platform runs out of physical memory and begin trashing. This
   ; module builds the suffix array and Lcp binary tree using "data management"
   ; techiques -- in particular posting sets which are designed for producing
   ; ordered listing of very large sets of information.
   ;
   ; global symbols defined: 
   ;
   ; void DbSuffixArray
   ;
   ; useful notes and assumptions:
   ;
   ; This module uses a purely memory based paging system. For instances where
   ; the Tstring becomes extremely large relative to platform memory a file-based
   ; paging system can be substituted.
   ;
   ; 01/16/00 written by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : as part of the Promula Version 8 TOPS development
   ;
   ;/hdr/ ************************************************************************
 */

#define OSYSFUNC					/* Gives access to the operating system */
#define PRIOFUNC					/* Gives access to PROMULA I/O functions */
#define RAWMFUNC					/* Gives access to raw memory operations */
#define STRGFUNC					/* Gives access to string manipulation */
#include "platform.h"					/* Define the platform hosting this code */
#include "LcpTree.h"
#include "TopsPage.h"					/* The operation of page streams */
#include "TopsPost.h"					/* The operation of postings sets */

extern UBYTE *Tstring;
static LONG DataSize;
static int PageProvider;

#ifndef PAGE_SIZE
#define PAGE_SIZE         1024				/* Size of each memory page */
#endif

#ifdef FPROTOTYPE
int ComparePostings(LONG Posting1, LONG Posting2)
#else
int ComparePostings(Posting1, Posting2)
  LONG Posting1;
  LONG Posting2;
#endif
{
  auto int iRelation;

  iRelation = cmpstrn(Tstring + Posting1, Tstring + Posting2, DataSize - Posting1);
  if (iRelation == 0)
    iRelation = -1;
  return iRelation;
}

#ifdef FPROTOTYPE
static int SuffixComputeLcp(LONG Posting1, LONG Posting2)
#else
static int SuffixComputeLcp(Posting1, Posting2)
  LONG Posting1;
  LONG Posting2;
#endif
{
  auto int Lcp;
  auto int Length;
  auto UBYTE *Str1;
  auto UBYTE *Str2;

  Str1 = Tstring + Posting1;
  Str2 = Tstring + Posting2;
  if (*Str1 != *Str2)
    return 0;
  if (Posting1 < Posting2)
    Length = DataSize - Posting2;
  else
    Length = DataSize - Posting1;
  if (Length > 255)
    Length = 255;
  for (Lcp = 1; Lcp < Length; Lcp++)
    if (Str1[Lcp] != Str2[Lcp])
      break;
  return Lcp;
}

#ifdef FPROTOTYPE
void DbSuffixArray(char *FileName, char *LcpName, LONG nString, int PageSize)
#else
void DbSuffixArray(FileName, LcpName, nString, PageSize)
  char *FileName;
  char *LcpName;
  LONG nString;
  int PageSize;
#endif
{
  auto LONG Position;
  auto int iCount;
  auto ULONG Buffer[1024];
  auto binfile fpSuffix;
  auto int Lcp;
  auto LONG Previous;

   /*-------------------------------------------------------------------------
   ;
   ; Step 1: Create a paging system to be used to support the posting set.
   ; Note that we are using a purely memory based paging system here.
   ;
   ;-------------------------------------------------------------------------*/

  if (PageSize == 0)
    PageSize = PAGE_SIZE;
  DataSize = nString;
  PageProvider = TopsPageMemory(PageSize);
  TopsPageSelect(PageProvider);
  TopsPostCreate(PageProvider, PageSize, ComparePostings);

   /*--------------------------------------------------------------------------
   ;
   ; Step 2: Just brute force sort the suffixs in the Tstring using the posting
   ; set and the "ComparePostings" function above. The terminal level in the
   ; posting set will literally be the suffix array.
   ;
   ;-------------------------------------------------------------------------*/

  for (Position = 0; Position < DataSize; Position++) {
    if (TopsPostWrite(Position) != NULL) {
      printf("Error encountered at Position %d\n", Position);
      exit(2);
    }
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 3: We now have our suffix array computed in the termimal level of the
   ; posting set. Simply buffer it out to the suffix array file and close
   ; the posting set and its supporting page system.
   ;
   ;-------------------------------------------------------------------------*/

  fpSuffix = inibinf(FileName);
  if (binopner(fpSuffix)) {
    printf("Unable to create Suffix Array File %s\n", FileName);
    exit(1);
  }
  iCount = 0;
  TopsPostFirst();
  for (;;) {
    if ((Position = TopsPostNext()) < 0)
      break;
    Buffer[iCount++] = Position;
    if (iCount == 1024) {
      wrbinf(fpSuffix, Buffer, 1024 * sizeof(LONG));
      iCount = 0;
    }
  }
  if (iCount != 0)
    wrbinf(fpSuffix, Buffer, iCount * sizeof(LONG));
  TopsPostDestroy();
  TopsPageDestroy();

   /*-------------------------------------------------------------------------
   ;
   ; Step 4: We have now created the suffix array, but we still need to create
   ; the Longest Common Prefix Tree. We use brute force again.
   ;
   ;-------------------------------------------------------------------------*/

  LcpCreate(DataSize);
  rewbinf(fpSuffix);
  iCount = 1;
  rdbinf(fpSuffix, Buffer, 1024 * sizeof(LONG));
  Previous = Buffer[0];
  for (Position = 1; Position < DataSize; Position++) {
    if (iCount == 1024) {
      rdbinf(fpSuffix, Buffer, 1024 * sizeof(LONG));
      iCount = 0;
    }
    Lcp = SuffixComputeLcp(Previous, Buffer[iCount]);
    LcpSet(Position, Lcp);
    Previous = Buffer[iCount++];
  }
  clsbinf(fpSuffix);
  LcpSaveTree(LcpName);
}
