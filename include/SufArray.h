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
   ; name of module: SufArray.h -- Suffix Arrays, Header
   ;
   ; This module provides utilities for working with Suffix arrays.
   ;
   ; useful notes and assumptions:
   ;
   ; 07/23/00 changed by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : to use the TOPS database/suffix array algorithm
   ;
   ;/hdr/ ************************************************************************
 */

#ifdef SUFARRAY_DECLARE
LONG *SuffixArray = NULL;
#else
extern LONG *SuffixArray;
#endif

#ifdef FPROTOTYPE
void SuffixArrayCreate (LONG size);
void SuffixArrayDestroy (void);
void SuffixArrayRestore (char *ArrayName, char *LcpName, LONG DataSize);
void SuffixArraySave (char *ArrayName, char *LcpName);
LONG SuffixArraySearchLw (char *pattern, LONG p, int *HaveMatch);
LONG SuffixArraySearchRw (char *pattern, LONG p);
int SuffixArraySearchLwRw (char *pattern, LONG p, LONG * LwRw);
void SuffixArraySort (void);
LONG *SuffixArrayVector (LONG index, LONG * Length);
#else
extern void SuffixArrayCreate ();
extern void SuffixArrayDestroy ();
extern void SuffixArrayRestore ();
extern void SuffixArraySave ();
extern LONG SuffixArraySearchLw ();
extern LONG SuffixArraySearchRw ();
extern int SuffixArraySearchLwRw ();
extern void SuffixArraySort ();
extern LONG *SuffixArrayVector ();
#endif
