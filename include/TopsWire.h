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
   ; name of module: TopsWire.h -- Wire File into Memory
   ;
   ; This module is part of the PROMULA TOPS code base. The operations provided
   ; via wiring allow the user to assign file storage to memory storage and visa
   ; versa. Wiring is used for two primary purposes
   ;
   ;  1) Memory sharing -- to provide independent processes simultaneous access
   ;     to the same information in memory.
   ;
   ;  2) Swap space control -- to provide a process access to virtual memory
   ;     areas without consuming swap space. This is needed especially on
   ;     WIN32 platforms where swap space tends to filled with user interace
   ;     stuff.
   ;
   ; Note that wiring is considered to be a utility operation; therefore, though
   ; it uses a provider approach, there is no currently selected wiring provider.
   ;
   ; global symbols defined: 
   ;
   ; int    TopsWireOpen        Wire an existing file into memory
   ; int    TopsWireCreate      Wire a new file into memory
   ; int    TopsWireClose       Close a wired file
   ; UBYTE* TopsWireAddress     Get memory address of wired file
   ; ULONG  TopsWireSize        Get size of wired file.
   ;
   ; 01/16/00 written by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : as part of the Promula Version 8 TOPS development
   ;
   ;/hdr/ ************************************************************************
 */

#ifdef FPROTOTYPE
UBYTE *TopsWireAddress (int Provider);
int TopsWireClose (int Provider);
int TopsWireCreate (char *FileName, ULONG FileSize);
int TopsWireOpen (char *FileName);
ULONG TopsWireSize (int Provider);
#else
extern UBYTE *TopsWireAddress ();
extern int TopsWireClose ();
extern int TopsWireCreate ();
extern int TopsWireOpen ();
extern ULONG TopsWireSize ();
#endif
