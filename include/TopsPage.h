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
   ; 01/16/00 written by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : as part of the Promula Version 8 TOPS development
   ;
   ;/hdr/ ************************************************************************
 */

#define TOPSPAGE_MEMORYSTREAM        0
#define TOPSPAGE_MEMORYAREA          1
#define TOPSPAGE_SEARCHFORBLOCK      2
#define TOPSPAGE_USEBLOCKVECTOR      3
#define TOPSPAGE_USEBLOCKSTREAM      4

#ifdef FPROTOTYPE
UBYTE *TopsPageAccess (ULONG BlockNumber, int Delta);
int TopsPageCreate (int PageSize, int MaxPages, int Algorithm);
void TopsPageDestroy (void);
int TopsPageMemory (int PageSize);
UBYTE *TopsPageNew (ULONG * iblock);
int TopsPageOpen (char *FileName, int Algorithm, int PageSize, int MaxPages);
int TopsPageSelect (int Provider);
ULONG TopsPageFileSize (void);
#else
extern UBYTE *TopsPageAccess ();
extern int TopsPageCreate ();
extern void TopsPageDestroy ();
extern int TopsPageMemory ();
extern UBYTE *TopsPageNew ();
extern int TopsPageOpen ();
extern int TopsPageSelect ();
extern ULONG TopsPageFileSize ();
#endif
