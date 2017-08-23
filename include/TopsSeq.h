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
   ; name of module: TopsSeq.h -- Tool for Organizing Properties Sequentially
   ;
   ; In this module, TOPS is viewed as referring to a "Tool for Organizing
   ; Properties Sequentially". No matter how carefully an information system is
   ; layed out to be persistent, there must ultimately be some set of property
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
   ; 01/16/00 written by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : as part of the Promula Version 8 TOPS development
   ;
   ;/hdr/ ************************************************************************
 */

#ifdef FPROTOTYPE
UBYTE *TopsSeqAccess (ULONG uSequence);
int TopsSeqCreate (ULONG iBlockSize, ULONG iRecordSize);
void TopsSeqDestroy ();
ULONG TopsSeqLength ();
int TopsSeqSelect (int Provider);
#else
extern UBYTE *TopsSeqAccess ();
extern int TopsSeqCreate ();
extern void TopsSeqDestroy ();
extern ULONG TopsSeqLength ();
extern int TopsSeqSelect ();
#endif
