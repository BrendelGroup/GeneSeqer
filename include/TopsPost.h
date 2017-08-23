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
   ; and then retreived by value.
   ;
   ; global symbols defined:
   ;
   ; LONG* TopsPostWrite    Write a new posting to the set
   ; void  TopsPostCreate   Create a posting stream
   ; void  TopsPostDestroy  Destroy the posting set
   ; void  TopsPostFirst    Position at first posting
   ; LONG  TopsPostNext     Get Next Posting
   ;
   ; useful notes and assumptions:
   ;
   ; 01/16/00 written by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : as part of the Promula Version 8 TOPS development
   ;
   ;/hdr/ ************************************************************************
 */

#ifdef FPROTOTYPE
LONG *TopsPostWrite (LONG Posting);
int TopsPostCreate (int PageProvider, int PageSize, int (*Compare) (LONG, LONG));
void TopsPostDestroy (void);
void TopsPostFirst (void);
LONG TopsPostNext (void);
#else
extern LONG *TopsPostWrite ();
extern int TopsPostCreate ();
extern void TopsPostDestroy ();
extern void TopsPostFirst ();
extern LONG TopsPostNext ();
#endif
