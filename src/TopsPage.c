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
   ; name of module: TopsPage.c -- Transfer Onto Paging System
   ;
   ; In this module TOPS is viewed as referring to the "Transfer onto a Paging
   ; System". It is used to store persistent blocks from a file in memory where
   ; the information in them can be used by applications. In the area of
   ; information management persistent storage must be viewed as containing a
   ; wide variety of different structures; indexes, lists, fixed records, long
   ; unformatted records, etc. However, ultimately all information must be
   ; physically stored in relatively large fixed length blocks, called "pages".
   ; This provider manages the pages.
   ;
   ; global symbols defined:
   ;
   ; UBYTE* TopsPageAccess   Access an existing block from a page stream
   ; int    TopsPageCreate   Create a file-based page stream
   ; void   TopsPageDestroy  Destroy a page stream
   ; ULONG  TopsPageFileSize Get current size of file in bytes
   ; int    TopsPageMemory   Create a memory-based page stream
   ; UBYTE* TopsPageNew      Create a new block in a page stream
   ; int    TopsPageOpen     Open a file-based page stream
   ; int    TopsPageSelect   Select the current page stream to be used
   ;
   ; useful notes and assumptions:
   ;
   ; There are two fundamentally different approaches made available by this 
   ; module:
   ;
   ;   1) Memory Mapping 
   ;   2) Virtual Paging
   ;
   ; Memory mapping is available on all contemporary platforms. It can be used
   ; to wire the file directly into the virtual memory system of the application;
   ; thus, allowing the file to be addressed directly as memory without actually
   ; reading it in. The two limitations of this approach are that file sizes are
   ; limited by the amount of addressable memory and that the size of the file
   ; must be known when it is opened. As an alternative, memory blocks can be
   ; managed via a sequential property stream. This approach avoids the need
   ; to know the size in advance, but pays a price in block location time.
   ;
   ; Virtual Paging has none of the above limitations, but it requires that 
   ; blocks be actually read from and written to the file. The major problem
   ; associated with virtual paging is that time required to determine whether
   ; a given block of the file is presently in memory. The implementation here
   ; uses three different schemes:
   ;
   ;   1) SearchForBlock -- Here a relatively short block, dimensioned by the
   ;         maximum number of pages to be allowed in the memory at one time,
   ;         is used which contains the block numbers of the blocks curently
   ;         in memory. This is the most general approach and can be used for
   ;         all files no matter how big; but since a search for the block
   ;         number must be made it is by far the slowest.
   ;
   ;   2) UseBlockVector -- Here the maximum number of blocks is specified and
   ;         a vector of ubytes is allocated which contains the memory page
   ;         number of the block. The advantage of this approach is that no
   ;         searching is required. The disadvantage is that a potentially
   ;         very long vector must be allocated.
   ;
   ;   3) UseBlockStream -- Here a sequential property stream is used to keep
   ;         track of which pages are in memory. This approach is the cleanest,
   ;         but is slower and has unpredictable memory usage behavior since
   ;         the page information vector itself uses dynamic memory.
   ;
   ; 01/16/00 written by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : as part of the Promula Version 8 TOPS development
   ;
   ;/hdr/ ************************************************************************
 */

#define MEMOFUNC					/* Gives access to memory allocation */
#define PRIOFUNC					/* Gives access to PROMULA I/O functions */
#define RAWMFUNC					/* Gives access to raw memory operations */
#include "platform.h"					/* Define the platform hosting this code */
#include "TopsSeq.h"					/* Tool organizing properties sequentially */
#include "TopsWire.h"					/* File to memory wiring operations */
#include "TopsPage.h"					/* Transfer onto a paging system */

typedef struct {
  int HasPageChanged;					/* Indicates if sheet has been changed */
  int LessRecentlyUsed;					/* Linkage for least recently used chain */
  int MoreRecentlyUsed;					/* Linkage for most recently used chain */
  ULONG BlockNumber;					/* Block number in stream */
  UBYTE *DataBlock;					/* Actual block of data */
} tPageInfo;						/* Virtual block control structure */

typedef struct {
  int MemProvider;					/* The memory page provider */
  binfile FileHandle;					/* The platform handle for the file */
  int Algorithm;					/* The algorithm used for finding blocks */
  UBYTE *PageNumbers;					/* The page numbers sorted by block number */
  int PagesInMemory;					/* Number of pages now in memory */
  int HasAnyPageChanged;				/* Has any page been changed */
  tPageInfo *Pages;					/* The actual pages in memory */
  int MaxNumberOfPages;					/* Maximum number of blocks */
  int PageMemoryExhausted;				/* Memory exhausted flag */
  int LeastRecentlyUsed;				/* Least-recently used record */
  int MostRecentlyUsed;					/* Most-recently user record */
  int PageSize;						/* Size of each page */
  ULONG MaxBlockNumber;					/* Maximum record number */
  ULONG BlockCount;					/* Number of blocks actually used */
  ULONG FileSize;					/* Total size of file */
} tPageSys;						/* Virtual memory system information */

#define MAX_PAGE_SYS    10				/* Maximum number of page systems */

static tPageSys PageInfo[MAX_PAGE_SYS];			/* Actual page system info */
static UBYTE PageProviders[MAX_PAGE_SYS + 1];		/* Provider control vector */
static tPageSys *Page;					/* Current page system information */
static int PageSelected = 0;				/* Number of current page system provider */

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

  nUsed = PageProviders[0];
  if (nUsed < MAX_PAGE_SYS) {
    nUsed++;
    Provider = nUsed;
    PageProviders[0] = (UBYTE) (nUsed);
    PageProviders[Provider] = 1;
  }
  else {
    for (Provider = 1; Provider <= MAX_PAGE_SYS; Provider++) {
      if (!PageProviders[Provider]) {
	PageProviders[Provider] = 1;
	break;
      }
    }
  }
  if (Provider <= MAX_PAGE_SYS)
    return Provider;
  return 0;
}

/*
   ;/doc/ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ;
   ; name of routine: LocRemovePage
   ;
   ; This function is called when a page is to be removed from its current
   ; position in the paging system. This removal occurs most frequently when the
   ; page in question is about to be elevated to the most recently used page. It
   ; is also used when a page is to be swapped out of memory.
   ;
   ; calling parameters:
   ;
   ; int   CurrrentPage   The number of the page to be removed.
   ;
   ; return parameters: None
   ;
   ; useful notes and assumptions: None.
   ;
   ;/doc/ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */
#ifdef FPROTOTYPE
static void LocRemovePage(int CurrentPage)
#else
static void LocRemovePage(CurrentPage)
  int CurrentPage;
#endif
{
  auto tPageInfo *ActualPage;
  auto tPageInfo *LinkedPage;

   /*-------------------------------------------------------------------------
   ;
   ; Step 1: If the page is the current most recently used page, then the
   ; page immediately preceding it becomes the most recently used one.
   ;
   ;------------------------------------------------------------------------*/

  ActualPage = Page->Pages + CurrentPage;
  if (CurrentPage == Page->MostRecentlyUsed) {
    if ((Page->MostRecentlyUsed = ActualPage->LessRecentlyUsed) != 0) {
      LinkedPage = Page->Pages + Page->MostRecentlyUsed;
      LinkedPage->MoreRecentlyUsed = 0;
    }
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 2: If the page is the currently least recently used page, then the
   ; page immediately following it becomes be least recently used page.
   ;
   ;-------------------------------------------------------------------------*/

  else if (CurrentPage == Page->LeastRecentlyUsed) {
    if ((Page->LeastRecentlyUsed = ActualPage->MoreRecentlyUsed) != 0) {
      LinkedPage = Page->Pages + Page->LeastRecentlyUsed;
      LinkedPage->LessRecentlyUsed = 0;
    }
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 3: If the page is neither the least recently nor most recently used
   ; page the pluck it out of both chains.
   ;
   ;------------------------------------------------------------------------*/

  else {
    if (ActualPage->MoreRecentlyUsed) {
      LinkedPage = Page->Pages + ActualPage->MoreRecentlyUsed;
      LinkedPage->LessRecentlyUsed = ActualPage->LessRecentlyUsed;
    }
    if (ActualPage->LessRecentlyUsed) {
      LinkedPage = Page->Pages + ActualPage->LessRecentlyUsed;
      LinkedPage->MoreRecentlyUsed = ActualPage->MoreRecentlyUsed;
    }
  }
}

/*
   ;/doc/++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ;
   ; name of routine: LocFindPage -- Find Page Associated with Block
   ;
   ; Given a block number this function attempts to find a page which represents
   ; that block currently in memory.
   ;
   ; calling parameters:
   ;
   ; ULONG CurrentBlock   The number of the blovck whose page is to be found.
   ;
   ; return parameters:
   ;
   ; The function returns the subscript within the PageNumbers vector of the page
   ; represents the specified block. It returns a -1 if no page is present.
   ;
   ; useful notes and assumptions:
   ;
   ; A standard binary search is performed. We use the relationship between the
   ; desired block and the current block in the list to determine where to look
   ; next. Given that the length of the vector is no more than 255 we will rarely
   ; get more that 7 comparisons and will usually get fewer. Though this is
   ; a bit slower than using a block number vector for small files, as the size of
   ; the files increase this approach might well become equivalent in terms of
   ; time requirements.
   ;
   ;/doc/+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */
#ifdef FPROTOTYPE
static int LocFindPage(ULONG CurrentBlock)
#else
static int LocFindPage(CurrentBlock)
  ULONG CurrentBlock;
#endif
{
  auto int iBottom;					/* Current bottom of search range */
  auto int iMiddle;					/* Current mid point of search range */
  auto int iRelation;					/* Comparison relation between keys */
  auto int iTop;					/* Current top of search range */
  auto tPageInfo *ActualPage;				/* Actual page being considered */
  auto UBYTE *PageNumbers;				/* Page numbers sorted by block number */
  auto tPageInfo *Pages;				/* The actual pages in memory */

  Pages = Page->Pages;
  PageNumbers = Page->PageNumbers;
  iBottom = iMiddle = 0;
  iTop = Page->PagesInMemory - 1;
  while (iBottom <= iTop) {
    iMiddle = (iBottom + iTop) / 2;
    ActualPage = Pages + (int) (PageNumbers[iMiddle]);
    iRelation = (int) (CurrentBlock - ActualPage->BlockNumber);
    if (!iRelation)
      return iMiddle;
    if (iRelation < 0)
      iTop = iMiddle - 1;
    else
      iBottom = iMiddle + 1;
  }
  return -1;
}

/*
   ;/doc/++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ;
   ; name of routine: LocRemovePageNumber - Remove Page Number
   ;
   ; The paging system maintains a sorted vector of page numbers which it uses to
   ; determine if a given block is currently in memory. When a page is swapped out
   ; of memory then its number must be removed from this list
   ;
   ; calling parameters:
   ;
   ; tPageInfo* ActualPage   The actual page whose number is to be removed.
   ;
   ; return parameters: None
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */
#ifdef FPROTOTYPE
static void LocRemovePageNumber(tPageInfo *ActualPage)
#else
static void LocRemovePageNumber(ActualPage)
  tPageInfo *ActualPage;
#endif
{
  auto UBYTE *PageNumbers;
  auto int CurrentPage;
  auto int PagesInMemory;

   /*-------------------------------------------------------------------------
   ;
   ; Step 1: Find the entry controlling the page. We know it is there because
   ; it is still registered with the paging memory system. Nonetheless we will
   ; check and if it is not there we will simply return.
   ;
   ;-------------------------------------------------------------------------*/

  CurrentPage = LocFindPage(ActualPage->BlockNumber);
  if (CurrentPage < 0)
    return;

   /*--------------------------------------------------------------------------
   ;
   ; Step 2: Remove the reference to this page from the page number vector
   ;
   ;------------------------------------------------------------------------*/

  PageNumbers = Page->PageNumbers;
  PagesInMemory = Page->PagesInMemory;
  PagesInMemory--;
  if (CurrentPage < PagesInMemory) {
    cpymem(PageNumbers + CurrentPage + 1, PageNumbers + CurrentPage,
      (PagesInMemory - CurrentPage));
  }
  Page->PagesInMemory = PagesInMemory;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: TopsPageAccess -- Access a Page
   ;
   ; This routine returns a pointer to the start of a page in the paging system
   ; which is identified by its sequential block number in the overall stream.
   ;
   ; calling parameters:
   ;
   ; ULONG BlockNumber  The number of the desired block, where blocks are numbered
   ;                    starting at 1, for the desired page.
   ;
   ; int   delta        A flag, which if non-zero, indicates that the calling
   ;                    function intends to change the information within the
   ;                    block. Only pages accessed via a non-zero delta are
   ;                    written back to their stream when the block is removed
   ;                    from memory.
   ;
   ; return parameters:
   ;
   ; A semi-permanent pointer to the start of the block or a NULL if something
   ; went wrong.
   ;
   ; useful notes and assumptions:
   ;
   ; The notion of "semi-permanent" pointers pervades the TOPS approach to
   ; information management. The basic idea is that a page pointer is a valid
   ; memory pointer until the block which it contains is swapped out of memory.
   ; Thus these pointers can be used safely within local blocks of code where
   ; efficiency is needed, even though they cannot be saved and reused over
   ; longer stretchs of code.
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
UBYTE *TopsPageAccess(ULONG BlockNumber, int Delta)
#else
UBYTE *TopsPageAccess(BlockNumber, Delta)
  ULONG BlockNumber;
  int Delta;
#endif
{
  auto UBYTE **blkadr;					/* Points to an address of a block address */
  auto tPageInfo *ActualPage;				/* An actual page being examimed */
  auto UBYTE *DataBlock;				/* Points to a block to contain stream data */
  auto UBYTE *PageNumbers;				/* The page number control vector */
  auto tPageInfo *Pages;				/* The actual pages in memory */
  auto int iRelation;					/* Comparison relation between keys */
  auto int CurrentPage;					/* Number of current page being examined */
  auto tPageInfo *MruBlock;				/* A more recently accessed page */
  auto int iBottom;					/* Current bottom of search range */
  auto int iMiddle;					/* Current mid point of search range */
  auto int iTop;					/* Current top of search range */
  auto int PagesInMemory;				/* Number of pages currently in memory */
  auto UBYTE *OldNumber;				/* Number of previous page in slot */

   /*-------------------------------------------------------------------------
   ;
   ; Step 1: If the pages are all memory based, then simply determine the
   ; address of the block in memory and return its pointer. If the pages are
   ; file based, then first check to see if we are accessing the currently
   ; most recently used page. If so we can gets its information directly
   ; and end local processing.
   ;
   ;-------------------------------------------------------------------------*/

  switch (Page->Algorithm) {
  case TOPSPAGE_MEMORYSTREAM:
    blkadr = (UBYTE **) (TopsSeqAccess(BlockNumber));
    return *blkadr;

  case TOPSPAGE_MEMORYAREA:
    return (UBYTE *) (Page->Pages) + (BlockNumber - 1) * Page->PageSize;

  default:
    if ((CurrentPage = Page->MostRecentlyUsed) != 0) {
      ActualPage = Page->Pages + CurrentPage;
      if (ActualPage->BlockNumber == BlockNumber) {
	if (Delta) {
	  Page->HasAnyPageChanged = 1;
	  ActualPage->HasPageChanged = 1;
	}
	return ActualPage->DataBlock;
      }
    }
  }

   /*-------------------------------------------------------------------------
   ;
   ; Step 2: We are using file-base pages and the requested block is not the
   ; most recently accessed one in the memory system. Determine if there
   ; currently is a page representing it in memory. If so, simply remove it
   ; from the order of reference chains and then branch so that it can become
   ; the most recently used block.
   ;
   ;-------------------------------------------------------------------------*/

  Pages = Page->Pages;
  switch (Page->Algorithm) {
  case TOPSPAGE_SEARCHFORBLOCK:
    PageNumbers = Page->PageNumbers;
    iRelation = LocFindPage(BlockNumber);
    if (iRelation >= 0) {
      CurrentPage = PageNumbers[iRelation];
      LocRemovePage(CurrentPage);
      goto UpdateChains;
    }
    break;

  case TOPSPAGE_USEBLOCKVECTOR:
    PageNumbers = Page->PageNumbers + (BlockNumber - 1);
    if ((CurrentPage = *PageNumbers) != 0) {
      LocRemovePage(CurrentPage);
      goto UpdateChains;
    }
    break;

  case TOPSPAGE_USEBLOCKSTREAM:
  default:
    PageNumbers = (UBYTE *) (TopsSeqAccess(BlockNumber));
    if ((CurrentPage = *PageNumbers) != 0) {
      LocRemovePage(CurrentPage);
      goto UpdateChains;
    }
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 3: The requested file block is not now in memory. If available
   ; memory has been exhausted, branch to reuse a page already in memory.
   ; Otherwise, allocate memory for a new page and initialize its information.
   ;
   ;------------------------------------------------------------------------*/

  if (!Page->PageMemoryExhausted) {
    DataBlock = (UBYTE *) (getmem(Page->PageSize));
    if (DataBlock != NULL) {
      CurrentPage = Page->PagesInMemory + 1;
      if (Page->Algorithm != TOPSPAGE_SEARCHFORBLOCK)
	Page->PagesInMemory += 1;
      ActualPage = Pages + CurrentPage;
      ActualPage->DataBlock = DataBlock;
      if (CurrentPage == Page->MaxNumberOfPages)
	Page->PageMemoryExhausted = 1;
      goto PageAvailable;
    }
    Page->PageMemoryExhausted = 1;
  }

   /*-------------------------------------------------------------------------
   ;
   ; Step 4: Page memory has been exhausted; therefore a page already in
   ; memory must be reused to support the requested block. By convention this
   ; is the currently least-recently used page. Obtain its number and remove it
   ; from the page reference chains.
   ;
   ;------------------------------------------------------------------------*/

  CurrentPage = Page->LeastRecentlyUsed;
  LocRemovePage(CurrentPage);
  ActualPage = Pages + CurrentPage;

  switch (Page->Algorithm) {
  case TOPSPAGE_SEARCHFORBLOCK:
    LocRemovePageNumber(ActualPage);
    break;

  case TOPSPAGE_USEBLOCKVECTOR:
    *(Page->PageNumbers + ActualPage->BlockNumber - 1) = 0;
    break;

  case TOPSPAGE_USEBLOCKSTREAM:
    OldNumber = (UBYTE *) (TopsSeqAccess(ActualPage->BlockNumber));
    *OldNumber = 0;
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 5: If the content of the datablock associated with the page to be
   ; reused has been changed then we must write its values back to the file
   ; before we use it to contain a different block.
   ;
   ;------------------------------------------------------------------------*/

  if (ActualPage->HasPageChanged) {
    putbinf(Page->FileHandle, (ActualPage->BlockNumber - 1) * Page->PageSize);
    wrbinf(Page->FileHandle, ActualPage->DataBlock, Page->PageSize);
    ActualPage->HasPageChanged = 0;
  }

PageAvailable:	  /*---------------------------------------------------------
   ;
   ; Step 6: A page is available for the block. First we must read then content
   ; of the block from the file into the memory page.
   ;
   ;-------------------------------------------------------------------------*/

  ActualPage->BlockNumber = BlockNumber;
  putbinf(Page->FileHandle, (BlockNumber - 1) * Page->PageSize);
  if (BlockNumber <= Page->BlockCount) {
    rdbinf(Page->FileHandle, ActualPage->DataBlock, Page->PageSize);
  }
  else {
    filmem(ActualPage->DataBlock, Page->PageSize, '\0');
    wrbinf(Page->FileHandle, ActualPage->DataBlock, Page->PageSize);
    Page->BlockCount += 1;
  }
  ActualPage->HasPageChanged = 0;

   /*-------------------------------------------------------------------------
   ;
   ; Step 7: We have the content of the block in memory. We must now update
   ; the block page number information. This is simple except in the case 
   ; where we are searching for page numbers. Here we must keep the vector
   ; of page numbers sorted by their block number.
   ;
   ;-------------------------------------------------------------------------*/

  switch (Page->Algorithm) {
  case TOPSPAGE_SEARCHFORBLOCK:
    PagesInMemory = Page->PagesInMemory;
    if (PagesInMemory == 0) {
      *PageNumbers = (UBYTE) (CurrentPage);
      Page->PagesInMemory = 1;
      goto UpdateChains;
    }
    ActualPage = Pages + (int) (PageNumbers[PagesInMemory - 1]);
    if (BlockNumber > ActualPage->BlockNumber) {
      PageNumbers[PagesInMemory++] = (UBYTE) (CurrentPage);
      Page->PagesInMemory = PagesInMemory;
      goto UpdateChains;
    }
    iBottom = iMiddle = iRelation = 0;
    iTop = PagesInMemory - 1;
    while (iBottom <= iTop) {
      iMiddle = (iBottom + iTop) / 2;
      ActualPage = Pages + (int) (PageNumbers[iMiddle]);
      iRelation = (int) (BlockNumber - ActualPage->BlockNumber);
      if (!iRelation)
	goto UpdateChains;
      if (iRelation < 0)
	iTop = iMiddle - 1;
      else
	iBottom = iMiddle + 1;
    }
    if (iRelation > 0)
      iMiddle++;
    PagesInMemory++;
    cpyovl(PageNumbers + iMiddle, PageNumbers + iMiddle + 1, PagesInMemory - iMiddle);
    PageNumbers[iMiddle] = (UBYTE) (CurrentPage);
    Page->PagesInMemory = PagesInMemory;
    goto UpdateChains;

  case TOPSPAGE_USEBLOCKVECTOR:
    *PageNumbers = (UBYTE) (CurrentPage);
    goto UpdateChains;

  case TOPSPAGE_USEBLOCKSTREAM:
    *PageNumbers = (UBYTE) (CurrentPage);
    goto UpdateChains;
  }

UpdateChains:	 /*-----------------------------------------------------------
   ;
   ; Step 8: If there is now a most recently block then this block becomes
   ; more recently used. Update the previous blocks information accordingly.
   ;
   ;------------------------------------------------------------------------*/

  ActualPage = Pages + CurrentPage;
  if (Page->MostRecentlyUsed) {
    MruBlock = Pages + Page->MostRecentlyUsed;
    MruBlock->MoreRecentlyUsed = CurrentPage;
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 9: Set the current block's information to make it the most recently
   ; used block.
   ;
   ;------------------------------------------------------------------------*/

  ActualPage->LessRecentlyUsed = Page->MostRecentlyUsed;
  ActualPage->MoreRecentlyUsed = 0;

   /*-------------------------------------------------------------------------
   ;
   ; Step 10: Make the current page the most recently used one in the page
   ; group, and if this is the first page in the page group then this page
   ; also becomes the current least recently used block.
   ;
   ;------------------------------------------------------------------------*/

  Page->MostRecentlyUsed = CurrentPage;
  if (!Page->LeastRecentlyUsed)
    Page->LeastRecentlyUsed = CurrentPage;

   /*-------------------------------------------------------------------------
   ;
   ; Step 11: We are ready to end processing as the requested stream block is
   ; now in memory. If the block is to be changed, mark its change indicator,
   ; and return a pointer to the actual block data area.
   ;
   ;-------------------------------------------------------------------------*/

  if (Delta) {
    Page->HasAnyPageChanged = 1;
    ActualPage->HasPageChanged = 1;
  }
  return ActualPage->DataBlock;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: TopsPageNew -- Obtain a New Page
   ;
   ; Obtains a pointer to and the number of a new block within a page stream
   ; that may be used to contain new information.
   ;
   ; calling parameters:
   ;
   ; ULONG*  iblock   Returns the block number of the new block
   ;
   ; return parameters:
   ;
   ; A semi-permanent pointer to the start of the memory page representing the
   ; block, or a NULL if something went wrong.
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
UBYTE *TopsPageNew(ULONG *iblock)
#else
UBYTE *TopsPageNew(iblock)
  ULONG *iblock;					/* Returns number of block */
#endif
{
  auto UBYTE *block;
  auto UBYTE **blkadr;
  auto ULONG BlockNumber;

  switch (Page->Algorithm) {
  case TOPSPAGE_MEMORYSTREAM:
    blkadr = (UBYTE **) (TopsSeqAccess(0));
    *iblock = TopsSeqLength();
    block = (UBYTE *) (getmem(Page->PageSize));
    *blkadr = block;
    return block;

  case TOPSPAGE_MEMORYAREA:
    if ((BlockNumber = Page->BlockCount) >= Page->MaxBlockNumber)
      return NULL;
    *iblock = BlockNumber + 1;
    return (UBYTE *) (Page->Pages) + BlockNumber * Page->PageSize;

  default:
    BlockNumber = Page->BlockCount + 1;
    if ((Page->MaxBlockNumber != 0) && BlockNumber > Page->MaxBlockNumber)
      return NULL;
    *iblock = BlockNumber;
    return TopsPageAccess(BlockNumber, 1);
  }
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: TopsPageDestroy -- Destroy a page System
   ;
   ; Destroys a page system by returning any memory occupied by it to the
   ; operating system. Before doing so, if the page system is backed by a
   ; persistent stream, it writes any of its blocks that have been changed in
   ; memory to the file hosting the stream and then physically closes that
   ; file.
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
void TopsPageDestroy(void)
#else
void TopsPageDestroy()
#endif
{
  auto ULONG i;
  auto UBYTE *mblock;
  auto ULONG Length;
  auto int CurrentPage;					/* Number of current page being examined */
  auto tPageInfo *ActualPage;				/* An actual page being examimed */
  auto int PagesInMemory;				/* Number of pages currently in memory */
  auto tPageInfo *Pages;				/* The actual pages in memory */

  PageProviders[PageSelected] = 0;
  if (Page->Algorithm == TOPSPAGE_MEMORYSTREAM) {
    Length = TopsSeqLength();
    for (i = 1; i <= Length; i++) {
      mblock = *(UBYTE **) (TopsSeqAccess(i));
      if (mblock != NULL)
	free(mblock);
    }
    TopsSeqDestroy();
    return;
  }
  if (Page->Algorithm == TOPSPAGE_MEMORYAREA) {
    TopsWireClose(Page->MemProvider);
    return;
  }
  if (Page->Algorithm == TOPSPAGE_USEBLOCKSTREAM)
    TopsSeqDestroy();
  else
    free(Page->PageNumbers);
  Pages = Page->Pages;
  PagesInMemory = Page->PagesInMemory;
  for (CurrentPage = 1; CurrentPage <= PagesInMemory; CurrentPage++) {
    ActualPage = Pages + CurrentPage;
    if (ActualPage->HasPageChanged) {
      putbinf(Page->FileHandle, (ActualPage->BlockNumber - 1) * Page->PageSize);
      wrbinf(Page->FileHandle, ActualPage->DataBlock, Page->PageSize);
    }
    free(ActualPage->DataBlock);
  }
  clsbinf(Page->FileHandle);
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: TopsPageMemory
   ;
   ; This function creates a purely memory based page stream. The only limit on
   ; the number of pages that can be allocated to the stream is the physical
   ; memory limit imposed by the platform.
   ;
   ; calling parameters:
   ;
   ; int    PageSize       The size of each page in the stream
   ;
   ; return parameters:
   ;
   ; If all goes well this function returns the provider number for the memory
   ; based page stream; else it returns zero.
   ;
   ; useful notes and assumptions:
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
int TopsPageMemory(int PageSize)
#else
int TopsPageMemory(PageSize)
  int PageSize;
#endif
{
  auto int Provider;
  auto tPageSys *Page;

   /*--------------------------------------------------------------------------
   ;
   ; Step 1: As a slave server it is up to this function to find an available
   ; provider service number. If one is not available, then return a zero.
   ;
   ;-------------------------------------------------------------------------*/

  Provider = LocGetProvider();
  if (Provider == 0)
    return 0;
  Page = PageInfo + (Provider - 1);

   /*--------------------------------------------------------------------------
   ;
   ; Step 2: There is an available provider. Obtain a provider for the property
   ; stream to be used to keep track of the memory pages. If no such provider
   ; can be supplied, return the page stream provider and return a zero.
   ;
   ;-------------------------------------------------------------------------*/

  if ((Page->MemProvider = TopsSeqCreate(PageSize, sizeof(UBYTE **))) == 0) {
    PageProviders[Provider] = 0;
    return 0;
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 3: All is well. Initialize the remaining parts of the paging 
   ; information and return the provider number.
   ;
   ;-------------------------------------------------------------------------*/

  Page->Algorithm = TOPSPAGE_MEMORYSTREAM;
  Page->PagesInMemory = 0;
  Page->HasAnyPageChanged = 0;
  Page->MaxNumberOfPages = 0;
  Page->PageMemoryExhausted = 0;
  Page->LeastRecentlyUsed = 0;
  Page->MostRecentlyUsed = 0;
  Page->PageSize = PageSize;
  Page->MaxBlockNumber = 0;
  return Provider;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: TopsPageSelect
   ;
   ; This function selects a specified page provider and returns the provider
   ; selection number of the previously selected provider.
   ;
   ; calling parameters:
   ;
   ; int    Provider       The provider selection number for the desired page
   ;                       stream. A value of zero selects the currently 
   ;                       selected provider.
   ; return parameters:
   ;
   ; If all goes well this function returns the prior provider selection number.
   ; This may be zero if there was no prior provider selection number.
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */

#ifdef FPROTOTYPE
int TopsPageSelect(int Provider)
#else
int TopsPageSelect(Provider)
  int Provider;
#endif
{
  auto int OldProvider;

  if (PageSelected == Provider)
    return Provider;
  OldProvider = PageSelected;
  if (Provider) {
    PageSelected = Provider;
    Page = PageInfo + (Provider - 1);
    if (Page->MemProvider)
      TopsSeqSelect(Page->MemProvider);
  }
  return OldProvider;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: TopsPageOpen
   ;
   ; This function opens a specified file and links it into a page provider.
   ;
   ; calling parameters:
   ;
   ; char*  FileName       The name of the file to be opened.
   ;
   ; int    Algorithm      The algorithm to be used in processing the file. It
   ;                       can have one of the following values:
   ; 
   ;                       TOPSPAGE_MEMORYAREA, TOPSPAGE_SEARCHFORBLOCK,
   ;                       TOPSPAGE_USEBLOCKVECTOR, TOPSPAGE_USEBLOCKSTREAM
   ;
   ; int    PageSize       The pagesize to be used in processing the file.
   ;
   ; int    MaxPages       The maximum number pages to be allowed in memory.
   ;
   ; return parameters:
   ;
   ; If all goes well this function returns the provider number for the memory
   ; based page stream; else it returns zero.
   ;
   ; useful notes and assumptions:
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
int TopsPageOpen(char *FileName, int Algorithm, int PageSize, int MaxPages)
#else
int TopsPageOpen(FileName, Algorithm, PageSize, MaxPages)
  char *FileName;
  int Algorithm;
  int PageSize;
  int MaxPages;
#endif
{
  auto int Provider;
  auto tPageSys *Page;

   /*--------------------------------------------------------------------------
   ;
   ; Step 1: As a slave server it is up to this function to find an available
   ; provider service number. If one is not available, then return a zero.
   ;
   ;-------------------------------------------------------------------------*/

  Provider = LocGetProvider();
  if (Provider == 0)
    return 0;
  Page = PageInfo + (Provider - 1);

   /*--------------------------------------------------------------------------
   ;
   ; Step 2: Initialize the scalar entrties in the virtual memory control
   ; table and branch depending upon the algorithm being used.
   ;
   ;-------------------------------------------------------------------------*/

  Page->MemProvider = 0;
  Page->FileHandle = nulbinf;
  Page->Algorithm = Algorithm;
  Page->PageNumbers = NULL;
  Page->PagesInMemory = 0;
  Page->HasAnyPageChanged = 0;
  Page->Pages = NULL;
  Page->MaxNumberOfPages = MaxPages;
  Page->PageMemoryExhausted = 0;
  Page->LeastRecentlyUsed = 0;
  Page->MostRecentlyUsed = 0;
  Page->PageSize = PageSize;
  Page->MaxBlockNumber = 0;
  Page->BlockCount = 0;
  Page->FileSize = 0;

   /*-------------------------------------------------------------------------
   ;
   ; Step 3: If the file is being memory mapped then we use the memory wiring
   ; operation rather than the standard platform binary I/O operations.
   ;
   ;-------------------------------------------------------------------------*/

  if (Algorithm == TOPSPAGE_MEMORYAREA) {
    if ((Page->MemProvider = TopsWireOpen(FileName)) == 0) {
      PageProviders[Provider] = 0;
      return 0;
    }
    Page->Pages = (tPageInfo *) (TopsWireAddress(Page->MemProvider));
    Page->FileSize = TopsWireSize(Page->MemProvider);
    Page->MaxBlockNumber = (Page->FileSize - 1) / Page->PageSize + 1;
    Page->BlockCount = Page->MaxBlockNumber;
    return Provider;
  }

   /*-------------------------------------------------------------------------
   ;
   ; Step 4: The blocks on the file are going to be managed independently of
   ; each other via memory pages. First open the file.
   ;
   ;------------------------------------------------------------------------*/

  Page->FileHandle = rdobinf(FileName);
  if (binopner(Page->FileHandle)) {
    PageProviders[Provider] = 0;
    return 0;
  }

   /*------------------------------------------------------------------------
   ;
   ; Step 5: We need a storage area that contains each page in the pages to
   ; be allocated. Note that we do 1-base subscripting so we need and extra
   ; page.
   ;
   ;------------------------------------------------------------------------*/

  Page->Pages = (tPageInfo *) (getmem((MaxPages + 1) * sizeof(tPageInfo)));
  if (Page->Pages == NULL) {
    PageProviders[Provider] = 0;
    return 0;
  }

   /*-------------------------------------------------------------------------
   ;
   ; Step 6: We are now ready to initialize the page number control vector.
   ; How we do this depends on the paging algorithm selected.
   ;
   ;-------------------------------------------------------------------------*/

  sizbinf(Page->FileHandle);
  Page->FileSize = posbinf(Page->FileHandle);
  Page->MaxBlockNumber = (Page->FileSize - 1) / Page->PageSize + 1;
  Page->BlockCount = Page->MaxBlockNumber;
  switch (Algorithm) {
  case TOPSPAGE_SEARCHFORBLOCK:
    Page->PageNumbers = (UBYTE *) (getmem(MaxPages + 1));
    if (Page->PageNumbers == NULL) {
      PageProviders[Provider] = 0;
      return 0;
    }
    filmem(Page->PageNumbers, MaxPages, '\0');
    break;

  case TOPSPAGE_USEBLOCKVECTOR:
    Page->PageNumbers = (UBYTE *) (getmem(Page->MaxBlockNumber));
    if (Page->PageNumbers == NULL) {
      PageProviders[Provider] = 0;
      return 0;
    }
    filmem(Page->PageNumbers, Page->MaxBlockNumber, '\0');
    break;

  case TOPSPAGE_USEBLOCKSTREAM:
    if ((Page->MemProvider = TopsSeqCreate(PageSize, sizeof(UBYTE))) == 0) {
      PageProviders[Provider] = 0;
      return 0;
    }
  }
  return Provider;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: TopsPageFileSize -- Get Size of Paged File
   ;
   ; Once a paged file has been created or opened, this funtion returns its
   ; total size.
   ;
   ; calling parameters: None
   ;
   ; return parameters:
   ;
   ; The function returns the total size of the file in bytes.
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
ULONG TopsPageFileSize(void)
#else
ULONG TopsPageFileSize()
#endif
{
  return Page->FileSize;
}
