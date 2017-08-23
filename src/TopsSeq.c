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
   ; name of module: TopsSeq.c -- Tool for Organizing Properties Sequentially
   ;
   ; In this module, TOPS is viewed as referring to a "Tool for Organizing
   ; Properties Sequentially". No matter how carefully an information system is
   ; laid out to be persistent, there must ultimately be some set of property
   ; information which must be stored in memory -- this can be symbol table
   ; information, virtual page control information, or whatever other set of
   ; control information is needed. To be consistent with the TOPS approaches
   ; this information is always organized into fixed length records that are
   ; accessed via a simple sequence number.
   ;
   ; global symbols defined:
   ;
   ; UBYTE* TopsSeqAccess        Access a property sequence
   ; int    TopsSeqCreate        Create a property sequential
   ; void   TopsSeqDestroy       Destroy a property sequence
   ; ULONG  TopsSeqLength        Return length of property sequence
   ; int    TopsSeqSelect        Select a property sequence to be the current one
   ;
   ; useful notes and assumptions:
   ;
   ; Sequential streams are viewed as streams consisting of sequences of
   ; relatively short fixed length records numbered from 1 to n, where n is the
   ; number of records in the stream. They are completely memory based and have
   ; no upper limit other than that imposed by the platform limits. Conceptually,
   ; the records of a sequential stream are viewed as records within a linear
   ; vector with no upper bound. The only problem is that since multiple
   ; sequential streams must share the same linear addressing scheme, each stream
   ; must be blocked and the addresses of those blocks must be retained in higher
   ; level blocks.
   ;
   ; NOTE THAT these routines are provided only for cases where maximum speed is
   ; required. They do no well-formedness checking of their calling parameters.
   ; If this capability is to be used for application programs directly, they
   ; should use the higher level fixed stream functions.
   ;
   ; 01/16/00 written by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : as part of the Promula Version 8 TOPS development
   ;
   ;/hdr/ ************************************************************************
 */

#define MEMOFUNC					/* Gives access to memory allocation */
#define RAWMFUNC					/* Gives access to raw memory operations */
#include "platform.h"					/* Define the platform hosting this code */
#include "TopsSeq.h"					/* Tool Organizing Properties Sequentially */

#define MAX_LEVELS          20				/* Maximum number of levels in mem stream */
#define MAX_SEQ              5				/* Maximum number of seq stream providers */

typedef struct {
  UBYTE **root;						/* Address of root block */
  ULONG maxseq;						/* Maximum record number */
  ULONG modulo;						/* Modulus of root block */
  ULONG bsize;						/* Block size */
  ULONG rsize;						/* Record size */
} tSeqStream;						/* Sequence stream structure type */

static tSeqStream SeqInfo[MAX_SEQ];			/* Actual stream information */
static UBYTE SeqProviders[MAX_SEQ + 1];			/* Providers control vector */
static int SeqSelected;					/* Currently selected provider */
static tSeqStream *SeqStream;				/* Currently selected information */

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
   ; providers available.
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

  nUsed = SeqProviders[0];
  if (nUsed < MAX_SEQ) {
    nUsed++;
    Provider = nUsed;
    SeqProviders[0] = (UBYTE) (nUsed);
    SeqProviders[Provider] = 1;
  }
  else {
    for (Provider = 1; Provider <= MAX_SEQ; Provider++) {
      if (!SeqProviders[Provider]) {
	SeqProviders[Provider] = 1;
	break;
      }
    }
  }
  if (Provider <= MAX_SEQ)
    return Provider;
  return 0;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: TopsSeqSelect
   ;
   ; This function selects a specified property sequence provider and returns the
   ; provider selection number of the previously selected provider.
   ;
   ; calling parameters:
   ;
   ; int    Provider       The provider selection number for the desired page
   ;                       stream. A value of zero selects the currently 
   ;                       selected provider.
   ; return parameters:
   ;
   ; This function returns the prior provider selection number. This may be zero
   ; if there was no prior provider selection number.
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
int TopsSeqSelect(int Provider)
#else
int TopsSeqSelect(Provider)
  int Provider;
#endif
{
  auto int OldProvider;

  OldProvider = SeqSelected;
  if (Provider) {
    SeqSelected = Provider;
    SeqStream = SeqInfo + (Provider - 1);
  }
  return OldProvider;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: TopsSeqAccess - Access a Sequential Stream
   ;
   ; This routine returns a pointer to a slot in memory which contains a specific
   ; fixed length property vector, specified via its subscript or sequence number
   ; (relative to one). This routine makes no attempt to check the specified
   ; sequence number for reasonableness. It assumes that no record is accessed
   ; before it has been created. If the next slot in sequence is desired, then
   ; pass a zero sequence number to this routine. The TopsSeqLength function can
   ; be used to obtain the current length of the sequence.
   ;
   ; The only possible error condition that can be encountered by this routine
   ; is insufficient memory to create a new block.
   ;
   ; calling parameters:
   ;
   ; ULONG uSequence     Sequence number of the desired record or zero if a new
   ;                     record is to be allocated to the stream.
   ;
   ; return parameters:
   ;
   ; The routine returns a pointer to the record whose sequence number was
   ; specified, or NULL if there was insufficient memory for a new block.
   ;
   ; useful notes and assumptions:
   ;
   ; For compatability with the information stream codes, this function actually
   ; accepts arbitrary record numbers. If a record number exceeds the current
   ; number of records, this function will create empty records to fill in any
   ; gaps. This feature will not be supported indefinitely as it complicates
   ; the logic, so always add new records one at a time using a zero sequence
   ; number as described above.
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
UBYTE *TopsSeqAccess(ULONG uSequence)
#else
UBYTE *TopsSeqAccess(uSequence)
  ULONG uSequence;
#endif
{
  auto ULONG i;						/* Dummy counter */
  auto ULONG isa;					/* Offset within current block */
  auto ULONG nbyte;					/* Number of bytes in record block */
  auto UBYTE **ioiad;					/* Address of current block */
  auto ULONG modulo;					/* Modulus of current block */
  auto UBYTE *ioa;					/* Address of a record block */
  auto ULONG recperblk;					/* Number of records per block */
  auto ULONG ptrperblk;					/* Number of pointers per block */

   /*------------------------------------------------------------------------
   ;
   ; Step 1: If the specified sequence number is within the current range
   ; of slots allocated for the sequence stream, branch directly to step 3
   ; to find its address.
   ;
   ;------------------------------------------------------------------------*/

  if (uSequence == 0 || uSequence > SeqStream->maxseq) {
      /*-----------------------------------------------------------------------
      ;
      ; Step 2: If the root sequence number block does not encompass the
      ; specified slot, then at least one new root block needs to be
      ; created. Note that the initial root block created is a record block
      ; while subsequent root blocks are pointer blocks. In the case where
      ; a wildly higher sequence number is specified, a HIGHLY dubious
      ; event, multiple levels of blocks may be created. Note the special
      ; case associated with the "zero" sequence number.
      ;
      ;---------------------------------------------------------------------*/

    if (uSequence == 0)
      uSequence = SeqStream->maxseq + 1;
    SeqStream->maxseq = uSequence;
    nbyte = SeqStream->bsize;
    recperblk = nbyte / SeqStream->rsize;
    ptrperblk = nbyte / sizeof(UBYTE **);
    for (;;) {
      modulo = SeqStream->modulo;
      if (modulo == 1) {
	if (uSequence <= recperblk)
	  break;
      }
      else if (uSequence <= (ptrperblk * modulo))
	break;
      if (SeqStream->root == NULL) {
	ioa = (UBYTE *) (getmem(nbyte));
	if (ioa == NULL)
	  return NULL;
	filmem(ioa, nbyte, '\0');
	SeqStream->root = (UBYTE **) (ioa);
	SeqStream->modulo = 1;
      }
      else {
	ioiad = (UBYTE **) (getmem(nbyte));
	if (ioiad == NULL)
	  return NULL;
	*ioiad = (UBYTE *) (SeqStream->root);
	for (i = 1; i < ptrperblk; i++)
	  *(ioiad + i) = NULL;
	SeqStream->root = ioiad;
	if (SeqStream->modulo == 1)
	  SeqStream->modulo = recperblk;
	else
	  SeqStream->modulo *= ptrperblk;
      }
    }
  }
   /*--------------------------------------------------------------------------
   ;
   ; Step 3: The initial sequence number block encompasses the sequence
   ; number. Initialize the variables needed to control the search for the
   ; specified slot. If the sequence stream consists of a single terminal
   ; block, simply return the address of the record from within that block.
   ;
   ;------------------------------------------------------------------------*/

  ioiad = SeqStream->root;
  modulo = SeqStream->modulo;
  if (modulo == 1)
    return ((UBYTE *) (ioiad)) + (uSequence - 1) * SeqStream->rsize;
  nbyte = SeqStream->bsize;
  recperblk = nbyte / SeqStream->rsize;
  ptrperblk = nbyte / sizeof(UBYTE **);

   /*--------------------------------------------------------------------------
   ;
   ; Step 4: Seek the specified slot. If it is found, return its address.
   ; Else if it is not yet in memory, continue to step 5. Notice that as we
   ; move down the tree the value of the sequence number is reduced by the
   ; modulus of the blocks.
   ;
   ;------------------------------------------------------------------------*/

  for (;;) {
    if (modulo == 1)
      return ((UBYTE *) (ioiad)) + (uSequence - 1) * SeqStream->rsize;
    isa = (uSequence - 1) / modulo;
    uSequence -= (isa * modulo);
    ioiad += isa;
    if (*ioiad == NULL)
      break;
    ioiad = (UBYTE **) (*ioiad);
    if (modulo == recperblk)
      modulo = 1;
    else
      modulo /= ptrperblk;
  }
   /*--------------------------------------------------------------------------
   ;
   ; Step 5: The record block that contains the specified slot is not yet
   ; in memory. Create any pointer blocks that are needed to dominate the
   ; record block.
   ;
   ;------------------------------------------------------------------------*/

  for (;;) {
    if (modulo == recperblk)
      break;
    ioa = (UBYTE *) (getmem(nbyte));
    if (ioa == NULL)
      return NULL;
    *ioiad = ioa;
    ioiad = (UBYTE **) (ioa);
    *ioiad = NULL;
    modulo /= ptrperblk;
    for (i = 1; i < ptrperblk; i++)
      *(ioiad + i) = NULL;
    isa = (uSequence - 1) / modulo;
    uSequence -= (isa * modulo);
    ioiad += isa;
  }
   /*--------------------------------------------------------------------------
   ;
   ; Step 6: Create the record block that contains the slot and return a
   ; pointer to the slot.
   ;
   ;------------------------------------------------------------------------*/

  ioa = (UBYTE *) (getmem(nbyte));
  if (ioa == NULL)
    return NULL;
  filmem(ioa, nbyte, 0);
  *ioiad = ioa;
  return ioa + (uSequence - 1) * SeqStream->rsize;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: TopsSeqCreate -- Create a Sequential Stream
   ;
   ; This function creates a memory-based sequential stream -- i.e. a memory-bound
   ; stream consisting of a sequence of relatively short fixed length records
   ; numbered from 1 to n, where n is the number of records in the stream.
   ;
   ; calling parameters:
   ;
   ; ULONG iBlockSize     Size of each memory block to be used in the stream.
   ;
   ; ULONG iRecordSize    Size of each record stored in the stream.
   ;
   ; return parameters:
   ;
   ; This function returns the provider number for the sequential stream or a
   ; zero if no provider can be supplied.
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
int TopsSeqCreate(ULONG iBlockSize, ULONG iRecordSize)
#else
int TopsSeqCreate(iBlockSize, iRecordSize)
  ULONG iBlockSize;					/* Size of each memory block in stream */
  ULONG iRecordSize;					/* Size of each record in stream */
#endif
{
  auto int Provider;

   /*--------------------------------------------------------------------------
   ;
   ; Step 1: Obtain a provider number for the memory stream. If none is
   ; available, return a 0 to the calling function.
   ;
   ;-------------------------------------------------------------------------*/

  if ((Provider = LocGetProvider()) == 0)
    return 0;
  SeqStream = SeqInfo + (Provider - 1);

   /*--------------------------------------------------------------------------
   ;
   ; Step 2: Initialize the stream control information and return its provider
   ; number.
   ;
   ;------------------------------------------------------------------------*/

  SeqStream->root = NULL;
  SeqStream->maxseq = SeqStream->modulo = 0;
  SeqStream->bsize = iBlockSize;
  SeqStream->rsize = iRecordSize;
  return Provider;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: TopsSeqDestroy - Destroy Sequential Stream
   ;
   ; This routine frees all memory associated a sequential stream. Note that this
   ; function assumes that no sequential stream will have more than MAX_LEVELS
   ; levels.
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
void TopsSeqDestroy(void)
#else
void TopsSeqDestroy()
#endif
{
  auto ULONG i;						/* Dummy counter */
  auto int nlevel;					/* Number of levels in the stream */
  auto UBYTE **ioiad;					/* Address of current block */
  auto ULONG modulo;					/* Modulus of current block */
  auto UBYTE **ioa[MAX_LEVELS];				/* Address of left-most block at level i */
  auto ULONG seq[MAX_LEVELS];				/* Sequence number of left-most block */
  auto int ilev;					/* Current block level being processed */
  auto ULONG nbyte;					/* Number of bytes in record block */
  auto ULONG recperblk;					/* Number of records per block */
  auto ULONG ptrperblk;					/* Number of pointers per block */

   /*------------------------------------------------------------------------
   ;
   ; Step 1: If there is no root block associated with this sequence stream,
   ; then the stream is already empty. Branch to step 8 to complete processing.
   ;
   ;------------------------------------------------------------------------*/

  if (SeqStream->root == NULL)
    goto Step08;

Step02:   /*-----------------------------------------------------------------
   ;
   ; Step 2: If the modulus of the root block is one, then this stream
   ; consists of a single record block. Free that block and branch to
   ; step 8 to complete processing.
   ;
   ;------------------------------------------------------------------------*/

  if (SeqStream->modulo == 1) {
    free(SeqStream->root);
    goto Step08;
  }
   /*--------------------------------------------------------------------------
   ;
   ; Step 3: The sequence stream contains at least one level of pointer
   ; blocks. Find the path to the current left-bottommost pointer block.
   ;
   ;------------------------------------------------------------------------*/

  nlevel = 0;
  ioiad = SeqStream->root;
  modulo = SeqStream->modulo;
  nbyte = SeqStream->bsize;
  recperblk = nbyte / SeqStream->rsize;
  ptrperblk = nbyte / sizeof(UBYTE **);
  while (modulo > recperblk) {
    ioa[nlevel] = ioiad;
    for (i = 0; i < ptrperblk; i++)
      if (*(ioiad + i) != NULL)
	break;
    seq[nlevel] = i;
    ioiad = (UBYTE **) (*(ioiad + i));
    modulo /= ptrperblk;
    nlevel++;
  }
Step04:   /*-----------------------------------------------------------------
   ;
   ; Step 4: We are now pointing to the leftmost bottommost block. Free all
   ; the blocks that it points to.
   ;
   ;------------------------------------------------------------------------*/

  for (i = 0; i < ptrperblk; i++) {
    if (*(ioiad + i) != NULL) {
      free(*(ioiad + i));
      *(ioiad + i) = NULL;
    }
  }
   /*------------------------------------------------------------------------
   ;
   ; Step 5: If there are no levels then all blocks have been freed. Branch
   ; to step 8 to complete processing.
   ;
   ;------------------------------------------------------------------------*/

  if (nlevel == 0) {
    free(ioiad);
    goto Step08;
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 6: Attempt to find the sister of the just freed block in the
   ; current terminal pointer block and loop back to step 4 if there is one.
   ;
   ;------------------------------------------------------------------------*/

  for (;;) {
    ilev = nlevel - 1;
    ioiad = ioa[ilev];
    for (i = seq[ilev] + 1; i < ptrperblk; i++)
      if (*(ioiad + i) != NULL) {
	seq[ilev] = i;
	ioiad = (UBYTE **) (*(ioiad + i));
	goto Step04;
      }
      /*---------------------------------------------------------------------
      ;
      ; Step 7: The just freed block has no sister in the current terminal
      ; pointer block. Move up the terminal pointer block chain until a
      ; sister is found. If no sister can be found, branch to step 8.
      ;
      ;---------------------------------------------------------------------*/

  
      if (ilev == 0)
	break;
      ilev--;
      ioiad = ioa[ilev];
      for (i = seq[ilev] + 1; i < ptrperblk; i++) {
	if (ioiad[i] != NULL) {
	  seq[ilev] = i;
	  ioiad = (UBYTE **) (ioiad[i]);
	  ilev++;
	  while (ilev < nlevel) {
	    ioa[ilev] = ioiad;
	    for (i = 0; i < ptrperblk; i++)
	      if (*(ioiad + i) != NULL)
		break;
	    seq[ilev] = i;
	    ioiad = (UBYTE **) (*(ioiad + i));
	    ilev++;
	  }
	  goto Step04;
	}
      }
    }
    SeqStream->modulo /= ptrperblk;
    goto Step02;
    
Step08:   /*-----------------------------------------------------------------
   ;
   ; Step 8: All blocks associated with the memory stream have been freed.
   ; Free the memor; stream structure itself and end local processing.
   ;
   ;-------------------------------------------------------------------------*/

  SeqStream->bsize = 0;
  SeqProviders[SeqSelected] = 0;
  return;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: TopsSeqLength - Return Length of Property Sequence
   ;
   ; This routine returns the current length of the property sequence.
   ;
   ; calling parameters: None
   ;
   ; return parameters:
   ;
   ; The number of property vectors currently stored in the sequence.
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
ULONG TopsSeqLength(void)
#else
ULONG TopsSeqLength()
#endif
{
  return SeqStream->maxseq;
}
