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
   ; name of module: TOPS.11 RawTextFile.h: Header, Raw Text File Services
   ;
   ; This header file contains the C-style property type names and service
   ; function prototypes needed to use the RawTextFile operational entity from C.
   ;
   ; See the file RawTextFile.c for a detailed description of this operational
   ; entity.
   ;
   ; 01/16/00 written by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : as part of the Promula Version 8 TOPS development
   ;
   ;/hdr/ ************************************************************************
 */

#define RAW_TEXT_EOF      -1

#ifdef FPROTOTYPE
void backupRawTextFile (int Selection, int Offset);
void closeRawTextFile (int Selection);
int getCharacterRawTextFile (int Selection);
int getPositionRawTextFile (int Selection);
int openRawTextFile (int Selection, char *FileName, int ReadOnly);
char *readRawTextFile (int Selection, char *Record, int nRecord);
void setPositionRawTextFile (int Selection, int Offset);
#else
extern void backupRawTextFile ();
extern void closeRawTextFile ();
extern int getCharacterRawTextFile ();
extern int getPositionRawTextFile ();
extern int openRawTextFile ();
extern char *readRawTextFile ();
extern void setPositionRawTextFile ();
#endif
