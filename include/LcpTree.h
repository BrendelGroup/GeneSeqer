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
   ; name of module: LcpTree.h -- Longest Common Prefix Tree, Header
   ;
   ; This module provides a Longest Common Prefix tree as support for the use
   ; and construction of suffix arrays.
   ;
   ; global symbols defined: None
   ;
   ; global functions defined:
   ;
   ; LONG   LcpCreate
   ; UBYTE  LcpLcp
   ; LONG   LcpLeftchild
   ; LONG   LcpMid
   ; UBYTE  LcpMinHeight
   ; LONG   LcpParent
   ; UBYTE  LcpRcp
   ; void   LcpRestoreTree
   ; LONG   LcpRightchild
   ; void   LcpSaveTree
   ; void   LcpSet
   ;
   ; useful notes and assumptions: None
   ;
   ;/hdr/ ************************************************************************
 */

#define LcpLeftchild(current)  (2*current+1)
#define LcpRightchild(current) (2*current+2)
#define LcpParent(current)     ((current-1)/2)

#ifdef FPROTOTYPE
LONG LcpCreate (LONG DataSize);
UBYTE LcpLcp (LONG Position);
LONG LcpLeafToIndex (LONG ith);
LONG LcpMid (LONG Left, LONG Right);
UBYTE LcpMinHeight (LONG left_leaf, LONG right_leaf);
UBYTE LcpRcp (LONG Position);
void LcpRestoreTree (char *TreeName);
void LcpDestroyTree (void);
void LcpSaveTree (char *TreeName);
void LcpSet (LONG ith, LONG value);
LONG LcpSearchLw (LONG Lw, int matchsize);
LONG LcpSearchRw (LONG Lw, int matchsize);
#else
extern LONG LcpCreate ();
extern UBYTE LcpLcp ();
extern LONG LcpLeafToIndex ();
extern LONG LcpMid ();
extern UBYTE LcpMinHeight ();
extern UBYTE LcpRcp ();
extern void LcpRestoreTree ();
extern void LcpDestroyTree ();
extern void LcpSaveTree ();
extern void LcpSet ();
extern LONG LcpSearchLw ();
extern LONG LcpSearchRw ();
#endif
