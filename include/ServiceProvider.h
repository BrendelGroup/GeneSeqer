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
   ; name of module: TOPS.2 ServiceProvider.h: Header, Service Provider Services
   ;
   ; This header file contains the C-style property type names and service
   ; function prototypes needed to use the ServiceProvider operational entity
   ; from C.
   ;
   ; See the file ServiceProvider.c for a detailed description of this
   ; operational entity.
   ;
   ; 01/16/00 written by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : as part of the Promula Version 8 TOPS development
   ;
   ;/hdr/ ************************************************************************
 */

#define NUMBER_OF_SERVICES              8

#define PROPERTY_VECTOR                 1
#define BINARY_FILE                     2
#define WIRED_FILE                      3
#define FORMATTED_INPUT                 4
#define FORMATTED_OUTPUT                5
#define MEMORY_BUFFER                   6
#define POSTING_SET                     7
#define STORAGE_AREA                    8

#define NUMBER_OF_PROPERTY_VECTORS     10
#define NUMBER_OF_BINARY_FILES         10
#define NUMBER_OF_WIRED_FILES           5
#define NUMBER_OF_FORMATTED_INPUTS      5
#define NUMBER_OF_FORMATTED_OUTPUTS     5
#define NUMBER_OF_MEMORY_BUFFERS        5
#define NUMBER_OF_POSTING_SETS          5
#define NUMBER_OF_STORAGE_AREAS         5

#define ERR_NO_MEMORY                   1
#define ERR_NO_PROVIDERS                2
#define ERR_NO_MEM_MAPPING              3
#define ERR_NO_SHARED_MEM               4
#define ERR_NOT_FINISHED                5
#define ERR_NO_ROOM                     6
#define ERR_NO_BLOCKS                   7
#define ERR_CREATING_FILE               8
#define ERR_FILE_NOT_EXIST              9
#define ERR_NOT_TOPS_AREA              10
#define ERR_BAD_FILE                   11

#ifdef FPROTOTYPE
int selectServiceProvider (int Service);
void reserveServiceProvider (int Service, int nReserve);
void returnServiceProvider (int Service, int Provider);
void fatalErrorInServiceProvider (int Service, int ErrorCode);
#else
extern int selectServiceProvider ();
extern void reserveServiceProvider ();
extern void returnServiceProvider ();
extern void fatalErrorInServiceProvider ();
#endif
