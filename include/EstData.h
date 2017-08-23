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
   ; name of module: EstData.h -- Process Est Data, Header
   ;
   ; global symbols defined: None
   ;
   ; global functions defined:
   ;
   ; void EstCreateDatIndFile
   ; LONG EstReadDatFile
   ;
   ; useful notes and assumptions: None
   ;
   ;/hdr/ ************************************************************************
 */

typedef struct
  {
    LONG giIndex;
    LONG PosInFile;
    int PosInSeq;
    LONG LenOfSeq;
    LONG OriginalPos;
  }
IndexOfEst;

#ifdef FPROTOTYPE
void EstCreateDatIndFiles (char *EST_Name, char *DataName, char *IndexName);
LONG EstReadDatFile (char *DataName, UBYTE ** Tstring);
void EstPackDatFile (char *DataName);
void EstDestroy (void);
#else
extern void EstCreateDatIndFiles ();
extern LONG EstReadDatFile ();
extern void EstPackDatFile ();
extern void EstDestroy ();
#endif
