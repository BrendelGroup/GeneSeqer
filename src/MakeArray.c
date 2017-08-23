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
   ; name of module: MakeDat.c -- Make Suffix Array dat file
   ;
   ; This module makes a suffix aray from a text file containing gene sequences.
   ;
   ; To use the program enter the nane of the EST data file with no extensions.
   ; In addition the following switches are available. Note they MUST be in lower
   ; case.
   ;
   ; bits   Specifies that the bucket control information should be stored in a
   ;        bit vector as opposed to a byte vector. Bits requires less memory but
   ;        more computational time.
   ;
   ; bytes  Specifies that the bucket control information should be stored in a
   ;        byte vector as opposed to a bit vector. Bytes requires more memory
   ;        but less computational time.
   ;
   ; usedat Specifies that an existing binary sequence file should be used rather
   ;        than creating a new one.
   ;
   ; loops  Specifies that multiple loops should be used to reduce thrashing
   ;
   ; useful notes and assumptions:
   ;
   ; By definition given an m-character string T, a suffix array for T is an array
   ; of integers in the range 1 to m, specifiying the lexicographic order of the
   ; m suffixs of string T. Consider the string
   ;
   ; "FredGoodman" of length 11. The suffixes then are
   ;
   ; Order  Suffix
   ;  1     fredgoodman
   ;  2     redgoodman
   ;  3     edgoodman
   ;  4     dgoodman
   ;  5     goodman
   ;  6     oodman
   ;  7     odman
   ;  8     dman
   ;  9     man
   ; 10     an
   ; 11     n
   ;
   ; Putting these in lexicographic order we get
   ;
   ; Order  Suffix
   ; 10     an
   ;  4     dgoodman
   ;  8     dman
   ;  3     edgoodman
   ;  1     fredgoodman
   ;  5     goodman
   ;  9     man
   ; 11     n
   ;  7     odman
   ;  6     oodman
   ;  2     redgoodman
   ;
   ; Therefore the suffix array for fredgoodman would be
   ;
   ; 10,4,8,3,1,5,9,11,7,6,2
   ;
   ; 07/23/00 changed by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : to use the TOPS database/suffix array algorithm
   ;
   ;/hdr/ ************************************************************************
 */

#define DATETIME					/* Gives access to date and time functions */
#define MEMOFUNC					/* Gives access to memory allocation */
#define OSYSFUNC					/* Gives access to the operating system */
#define PRIOFUNC					/* Gives access to PROMULA I/O functions */
#define RAWMFUNC					/* Gives access to raw memory operations */
#define STRGFUNC					/* Gives access to string manipulation */
#include "platform.h"					/* Define the platform hosting this code */
#include "EstData.h"
#include "SufArray.h"

#ifdef FPROTOTYPE
static int sys_secs(void)
#else
static int sys_secs()
#endif
/* Seconds timer */
{
  return clock() / CLOCKS_PER_SEC;
}

#ifdef FPROTOTYPE
void DbSuffixArray(char *FileName, char *LcpName, LONG nString, int PageSize);
#else
extern void DbSuffixArray();
#endif

UBYTE *Tstring = NULL;

#define USING_BH_BITS              0
#define USING_BHBYTES              1

static int CurrentApproach = USING_BH_BITS;
static int UseDat = 0;
static int SeparateLoops = 0;
static int UseDB = 1;
static int PackDat = 0;
static int DatOnly = 0;
static int ShowTime = 0;

#ifdef FPROTOTYPE
int main(int argc, char **argv)
#else
int main(argc, argv)
  int argc;
  char **argv;
#endif
{
  auto char IndexName[257];
  auto char ArrayName[257];
  auto char LcpName[257];
  auto char DataName[257];
  auto int nChar;
  auto int iArg;
  auto LONG DataSize;
  auto int PageSize;

   /*-------------------------------------------------------------------------
   ;
   ; Step 1: If the filename of the source file is not specified on the
   ; command line, then issue a usage mesage and exit. Else check for any
   ; additional processing requests.
   ;
   ;-------------------------------------------------------------------------*/

  if (argc < 2) {
    printf("Usage:\n\n    %s dbest\n\n  dbest :   EST sequence data file (FASTA-format).\n\n", argv[0]);
    return 1;
  }
  PageSize = 0;
  for (iArg = 2; iArg < argc; iArg++) {
    if (cmpstrn(argv[iArg], "bits", 5) == 0)
      CurrentApproach = USING_BH_BITS;
    else if (cmpstrn(argv[iArg], "bytes", 6) == 0)
      CurrentApproach = USING_BHBYTES;
    else if (cmpstrn(argv[iArg], "usedat", 7) == 0)
      UseDat = 1;
    else if (cmpstrn(argv[iArg], "loops", 6) == 0)
      SeparateLoops = 1;
    else if (cmpstrn(argv[iArg], "uselcp", 7) == 0)
      UseDB = 0;
    else if (cmpstrn(argv[iArg], "packdat", 8) == 0)
      PackDat = 1;
    else if (cmpstrn(argv[iArg], "datonly", 8) == 0)
      DatOnly = 1;
    else if (cmpstrn(argv[iArg], "time", 5) == 0)
      ShowTime = 1;
  }
  nChar = strlen(argv[1]);
  cpymem(argv[1], IndexName, nChar);
  cpymem(".ind", IndexName + nChar, 5);
  cpymem(argv[1], ArrayName, nChar);
  cpymem(".suf", ArrayName + nChar, 5);
  cpymem(argv[1], LcpName, nChar);
  cpymem(".tre", LcpName + nChar, 5);
  cpymem(argv[1], DataName, nChar);
  cpymem(".dat", DataName + nChar, 5);

  ShowTime = sys_secs();

   /*-------------------------------------------------------------------------
   ;
   ; Step 2: If we are to create the data and index files from the EST
   ; database name, then do so now.
   ;
   ;------------------------------------------------------------------------*/

  if (!UseDat)
    EstCreateDatIndFiles(argv[1], DataName, IndexName);
  if (PackDat)
    EstPackDatFile(DataName);
  if (DatOnly)
    return 0;

   /*-------------------------------------------------------------------------
   ;
   ; Step 3: Obtain the length of the sequence string and read it in.
   ;
   ;-------------------------------------------------------------------------*/

  DataSize = EstReadDatFile(DataName, &Tstring);
    if(DataSize == 0){
        /* release memory and delete file */
        printf("Sorry. This EST database is empty.\n\n");
        delfile(DataName);
        delfile(IndexName);
        delfile(ArrayName);
        return -1;
    }

   /*-------------------------------------------------------------------------
   ;
   ; Step 4: Create the basic structures needed to build the longest common
   ; prefix tree and then initialize the storage areas need to build the
   ; suffix array. 
   ;
   ;------------------------------------------------------------------------*/

  if (UseDB) {
    DbSuffixArray(ArrayName, LcpName, DataSize, PageSize);
  }
  else {
    SuffixArrayCreate(DataSize);
    SuffixArraySort();
    SuffixArraySave(ArrayName, LcpName);
    SuffixArrayDestroy();
  }
  printf("Elapsed seconds = %d\n", sys_secs() - ShowTime);
  return 0;
}
