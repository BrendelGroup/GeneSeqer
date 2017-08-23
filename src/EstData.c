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
   ; name of module: EstData.c -- Process Est Data, Implementation
   ;
   ; global symbols defined: None
   ;
   ; global functions defined:
   ;
   ; void EstCreateDatIndFile
   ; LONG EstReadDatFile
   ;
   ; useful notes and assumptions:
   ;
   ; The binary data vector consists of a series of byte in the range 0-3
   ;
   ;  T = 0, C = 1, A = 2, G = 3
   ;
   ; The original code states that the binary data also contains end-of-sequence
   ; codes of 4; however, no such codes are ever written to the dat file.
   ;
   ; 07/23/00 changed by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : to use the TOPS database/suffix array algorithm
   ;
   ;/hdr/ ************************************************************************
 */

#define OSYSFUNC					/* Gives access to the operating system */
#define MEMOFUNC					/* Gives access to memory allocation */
#define PRIOFUNC					/* Gives access to PROMULA I/O functions */
#define RAWMFUNC					/* Gives access to raw memory operations */
#define STRGFUNC					/* Gives access to string manipulation */
#include "platform.h"					/* Define the platform hosting this code */
#include "EstData.h"
#include "RawTextFile.h"				/* Header, raw text file services */
#include "TopsWire.h"					/* The operation if wiring files to memory */

#define MinSeqLen               12
#define BufferSize       	20000
#define EST_FILE                1

#define NA 124
#define EL 125
#define ES 126


/* T=U=0,C=1,A=2,G=3, NSYWRKM >3, others > 15 */	
int nascii[]={
/*       0  1  2  3  5  6  7  8  9 10 11 12 13 14 15 15
         @  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O
         P  Q  R  S  T  U  V  W  X  Y  Z                */
        EL,NA,NA,NA,NA,NA,NA,NA,NA,NA,EL,NA,NA,EL,NA,NA,
        NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
        NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,16,NA,NA,
        NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
        ES, 2,14, 1,11,NA,NA, 3,12,NA,NA,10,NA, 7,15,NA,
         5, 6, 5, 9, 0, 0,13, 8,15, 6,NA,NA,NA,NA,NA,NA,
        NA, 2,14, 1,11,NA,NA, 3,12,NA,NA,10,NA, 7,15,NA,
         5, 6, 5, 9, 0, 0,13, 8,15, 6,NA,NA,NA,NA,NA,NA};


static int EstProvider = 0;



/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: EstCreateDatIndFiles
   ;
   ; Given the name of a sequence file in EST format this function creates a
   ; binary data file and index file which can be used to build and use a
   ; suffix array based on the sequence string.
   ;
   ; calling parameters:
   ;
   ; char*   EST_Name    The name of the input file which is assumed to be in
   ;                     EST format.
   ;
   ; char*   DataName    The name of the binary data file to be created
   ;
   ; char*   IndexName   The name of the index file to be created
   ;
   ; return parameters:
   ;
   ; None -- if there is a problem this function displays an informative message
   ; and exits to the operating system.
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
void EstCreateDatIndFiles(char *EST_Name, char *DataName, char *IndexName)
#else
void EstCreateDatIndFiles(EST_Name, DataName, IndexName)
  char *EST_Name;
  char *DataName;
  char *IndexName;
#endif
{
  auto binfile DatHandle;
  auto binfile IndHandle;
  auto IndexOfEst theIndex;
  auto char NewLine[8192];
  auto int iChar;
  auto char NewChar;
  auto int estCount;
  auto int LenOfString;
  auto int PosOfString;
  auto UBYTE newEST[BufferSize];
  auto int PosOfLine;
  auto int giIndex;
  auto int PosInDataFile;
  auto LONG LengthOfDat;
  auto int p;
  auto int c;
  auto LONG totalLen;
  auto int minLen=0,maxLen=0,currLen;
  auto int stat[11],s,i;
  
  for(i=0;i<11;i++) stat[i]=0;
   /*-------------------------------------------------------------------------
   ;
   ; Step 2: Open the input file and create the two output files. If there is
   ; a problem, issue a message and exit.
   ;
   ;-------------------------------------------------------------------------*/

  if (openRawTextFile(EST_FILE, EST_Name, 1) != 0) {
    printf("Unable to open EST file %s\n", EST_Name);
    exit(1);
  }
  DatHandle = inibinf(DataName);
  if (binopner(DatHandle)) {
    printf("Unable to create DataFile %s\n", DataName);
    exit(1);
  }
  IndHandle = inibinf(IndexName);
  if (binopner(IndHandle)) {
    printf("Unable to create IndexFile %s\n", IndexName);
    exit(1);
  }

   /*-------------------------------------------------------------------------
   ;
   ; Step 3: Read records from the input file until the ">" is encountered.
   ;
   ;-------------------------------------------------------------------------*/
  totalLen=0;
  estCount = 0;
  LengthOfDat = 0;
  PosOfLine = 0;
  NewLine[0] = '\0';
  for (;;) {
    PosOfLine = getPositionRawTextFile(EST_FILE);
    if (readRawTextFile(EST_FILE, NewLine, sizeof(NewLine)) == NULL) {
      printf("The EstsFile %s contains no sequences\n", EST_Name);
      exit(1);
    }
    for (iChar = 0; (NewChar = NewLine[iChar]) != '\0'; iChar++) {
      if (NewChar == '>')
	goto ReadSequence;
    }
  }

ReadSequence:	/*-----------------------------------------------------------
   ;
   ; Step 4: We are ready to read the current sequence. If an end of file
   ;
   ;------------------------------------------------------------------------*/

  theIndex.OriginalPos = PosOfLine + iChar + 1;
  estCount++;
  giIndex = estCount;
  PosInDataFile = posbinf(DatHandle);
  LenOfString = 0;
  PosOfString = 0;
  p=0;
  currLen=0;
  for (;;) {
    PosOfLine = getPositionRawTextFile(EST_FILE);
    if (readRawTextFile(EST_FILE, NewLine, sizeof(NewLine)) == NULL) {
      if (strlen(NewLine)==0) {
        NewChar = -1;
        goto FoundTheEnd;
      }
    }
    for (iChar = 0; (NewChar = NewLine[iChar]) != '\0'; iChar++) {
      if (NewChar=='>') goto FoundTheEnd;
      c=nascii[(int)NewChar];
      if(c<=3){
	if (p >BufferSize-2){
	    /* dump the buffer */
            wrbinf(DatHandle, newEST, p);
            LenOfString+=p;
            p=0;
        }
	newEST[p++] = c;
        currLen++;
      }else if(c<=15){
        currLen++;
        LenOfString+=p;
	if (LenOfString >= MinSeqLen) {
	    theIndex.giIndex = giIndex;
	    theIndex.LenOfSeq = LenOfString;
	    theIndex.PosInSeq = PosOfString;
	    theIndex.PosInFile = PosInDataFile;
	    wrbinf(IndHandle, &theIndex, sizeof(theIndex));
	    if(p>0) wrbinf(DatHandle, newEST, p);
	    LengthOfDat += LenOfString;
	    PosInDataFile = posbinf(DatHandle);
	}
        PosOfString += (LenOfString + 1);
	LenOfString = 0;
        p=0;
       }/* end of else if */
    }/* end of for(iChar) */
	
  }/* end of for(;;) */
FoundTheEnd:  /*-------------------------------------------------------------
   ;
   ; Step 5: We have found the end of the current sequence either via an '>'
   ; or via and end-of-file.
   ;
   ;-------------------------------------------------------------------------*/

  LenOfString+=p;
  
  if (LenOfString >= MinSeqLen) {
    theIndex.giIndex = giIndex;
    theIndex.LenOfSeq = LenOfString;
    theIndex.PosInSeq = PosOfString;
    theIndex.PosInFile = PosInDataFile;
    wrbinf(IndHandle, &theIndex, sizeof(theIndex));
    if(p>0)wrbinf(DatHandle, newEST, p);
    LengthOfDat += LenOfString;
  }
  
  
  if(estCount==1){
    minLen=maxLen=currLen;
  }else{
    if(currLen>maxLen) maxLen=currLen;
    if(currLen<minLen) minLen=currLen;
  }
  totalLen+=currLen;
  
  s=currLen/100;
  if(s>10) s=10;
  stat[s]++;
  
  if (NewChar == '>')
    goto ReadSequence;

   /*--------------------------------------------------------------------------
   ;
   ; Step 6: We have processed the entire set of sequence data. Close any
   ; handles and return.
   ;
   ;-------------------------------------------------------------------------*/

  printf("Total number of ESTs: %d\tTotal sequence length: %d\n",estCount,totalLen);
  printf("Minimum sequence length: %d\tMaximum sequence length: %d \n",minLen,maxLen);
  printf("\nLength distribution (number of sequences of specified length):\n");
  for(i=0;i<10;i++){
    printf("<%5d:\t%d\n",(i+1)*100,stat[i]);
  }
  printf(">=1000:\t%d\n\n",stat[10]);
  
  closeRawTextFile(EST_FILE);
  clsbinf(IndHandle);
  clsbinf(DatHandle);
  return;

}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: EstReadDatFile
   ;
   ; Given the name of a sequence file in EST format this function reads its
   ; corresponding binary data file.
   ;
   ; calling parameters:
   ;
   ; char*   DataName    The name of the binary data file name.
   ;
   ; UBYTE** Tstring     Returns a pointer to the memory area that contains the
   ;                     entire binary data file.
   ;
   ; return parameters:
   ;
   ; If all goes well this function returns the length of the sequence string.
   ; If there is a problem this function displays an informative message
   ; and exits to the operating system.
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
LONG EstReadDatFile(char *DataName, UBYTE **Tstring)
#else
LONG EstReadDatFile(DataName, Tstring)
  char *DataName;
  UBYTE **Tstring;
#endif
{
  auto LONG DataSize;

  EstProvider = TopsWireOpen(DataName);
  if (EstProvider == 0) {
    printf("Unable to open DataFile %s\n", DataName);
    exit(1);
  }
  DataSize = TopsWireSize(EstProvider);
  *Tstring = TopsWireAddress(EstProvider);
  return DataSize;
}

#ifdef FPROTOTYPE
void EstDestroy(void)
#else
void EstDestroy()
#endif
{
  if (EstProvider != 0)
    TopsWireClose(EstProvider);
  EstProvider = 0;
}

#ifdef FPROTOTYPE
void EstPackDatFile(char *DataName)
#else
void EstPackDatFile(DataName)
  char *DataName;
#endif
{
  auto binfile DatFile;
  auto binfile PacFile;
  auto UBYTE DatBuffer[4096];
  auto UBYTE PacBuffer[1024];
  auto LONG DataSize;
  auto char PackName[257];
  auto int nName;
  auto int nRead;
  auto int iDat;
  auto int iPac;

   /*-------------------------------------------------------------------------
   ;
   ; Step 1: Open the dat file and obtain its size.
   ;
   ;-------------------------------------------------------------------------*/

  DatFile = rdobinf(DataName);
  if (binopner(DatFile)) {
    printf("Unable to open EST DataFile %s\n", DataName);
    exit(1);
  }
  sizbinf(DatFile);
  DataSize = posbinf(DatFile);
  rewbinf(DatFile);

   /*-------------------------------------------------------------------------
   ;
   ; Step 2: Create the packed data file.
   ;
   ;------------------------------------------------------------------------*/

  nName = strlen(DataName);
  cpymem(DataName, PackName, nName + 1);
  PackName[nName - 3] = 'p';
  PackName[nName - 2] = 'a';
  PackName[nName - 1] = 'c';
  PacFile = inibinf(PackName);
  if (binopner(PacFile)) {
    printf("Unable to create packed EST DataFile %s\n", PackName);
    exit(1);
  }
   /*--------------------------------------------------------------------------
   ;
   ; Step 3: Do the compaction.
   ;
   ;-------------------------------------------------------------------------*/

  nRead = 4096;
  for (nRead = 4096; DataSize > 0; DataSize -= nRead) {
    if (nRead > DataSize)
      nRead = DataSize;
    rdbinf(DatFile, DatBuffer, nRead);
    for (iDat = iPac = 0; iDat < nRead; iDat += 4, iPac++) {
      PacBuffer[iPac] = (UBYTE) (DatBuffer[iDat] * 64 + DatBuffer[iDat + 1] * 16 +
	DatBuffer[iDat + 2] * 8 + DatBuffer[iDat + 3]);
    }
    wrbinf(PacFile, PacBuffer, iPac);
  }
  clsbinf(DatFile);
  clsbinf(PacFile);
}
