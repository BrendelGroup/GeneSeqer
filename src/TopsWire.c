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
   ; name of module: TopsWire.c -- Wire File into Memory
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
   ; useful notes and assumptions:
   ;
   ; This particular implementation of the of the wiring application programming
   ; interface is written for the Win32 API as distributed by the Microsoft
   ; Corporation. It has been tested with Windows NT, Version 4.0.
   ;
   ; In addition there is a platform independent implementation which can
   ; be used for single process applications, but which does not support
   ; either shared access or swap space control.
   ;
   ; 01/16/00 written by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : as part of the Promula Version 8 TOPS development
   ;
   ;/hdr/ **********************************************************************
*/

#ifdef MSCPLAT /*---------------------------------------------------------------
;
; Note that the MicroSoft headers generate many compiler warnings so we must
; turn these off and then on again if we wish to check our own code for
; warnings.
;
;----------------------------------------------------------------------------*/
#pragma warning(disable:4081)
#pragma warning(disable:4100)
#pragma warning(disable:4115)
#pragma warning(disable:4214)
#pragma warning(disable:4201)
#pragma warning(disable:4514)
#pragma warning(disable:4701)
#define STRICT
#include <windows.h>
#pragma warning(default:4081)
#pragma warning(default:4100)
#pragma warning(default:4115)
#pragma warning(default:4214)
#pragma warning(default:4201)
#pragma warning(default:4701)
#endif

#define MEMOFUNC                                        /* Gives access to memory allocation */
#define PRIOFUNC                                        /* Gives access to PROMULA I/O functions */
#include "platform.h"                                   /* Define the platform hosting this code */
#include "TopsWire.h"                                   /* File to memory wiring operations */

typedef struct {
   UBYTE *WireArea;                                     /* Memory address of file content */
   ULONG WireFileSize;                                  /* Total size of file in bytes */
#ifdef MSCPLAT
   HANDLE WireFile;                                     /* WIN32 file handle */
   HANDLE WireMapping;                                  /* WIN32 mapping handle */
   HANDLE WireView;                                     /* WIN32 view handle */
#else
   binfile WireFile;                                    /* Handle of file */
   int WriteStatus;                                     /* Should memory be written on close */
#endif
} tWireInfo;                                            /* Wiring provider control structure type */

#define WIRE_PROVIDERS    10                            /* Maximum number of wiring providers */

static tWireInfo WireInfo[WIRE_PROVIDERS];              /* Actual wiring information */
static UBYTE WireProviders[WIRE_PROVIDERS + 1];         /* Providers control vector */

/*
   ;/doc/ **********************************************************************
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
   ; providers available.
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ **********************************************************************
*/
#ifdef FPROTOTYPE
static int LocGetProvider(void)
#else
static int LocGetProvider()
#endif
{
   auto int nUsed;
   auto int Provider;

   nUsed = WireProviders[0];
   if (nUsed < WIRE_PROVIDERS)
   {
      nUsed++;
      Provider = nUsed;
      WireProviders[0] = (UBYTE) (nUsed);
      WireProviders[Provider] = 1;
   }
   else
   {
      for (Provider = 1; Provider <= WIRE_PROVIDERS; Provider++)
      {
         if (!WireProviders[Provider])
         {
            WireProviders[Provider] = 1;
            break;
         }
      }
   }
   if (Provider <= WIRE_PROVIDERS)
      return Provider;
   return 0;
}

/*
   ;/doc/ **********************************************************************
   ;
   ; name of routine: TopsWireOpen
   ;
   ; This function opens the file backing up the memory area and prepares it for
   ; use.
   ;
   ; calling parameters:
   ;
   ; char* FileName   The name of file whose content is to be wired into memory.
   ;
   ; return parameters:
   ;
   ; The function returns the number of the provider supplied for this instance;
   ; or a zero if there is a problem.
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ***********************************************************************
*/
#ifdef FPROTOTYPE
int TopsWireOpen(char *FileName)
#else
int TopsWireOpen(FileName)
   char *FileName;
#endif
{
   auto tWireInfo *This;
   auto int Provider;

#ifdef MSCPLAT /*=============================================================
   ;
   ; MSCPLAT code. The approach is straightforward unded WIN32. Open the
   ; file, get its size, create a file mapping, and finally map that into
   ; memory.
   ;
   ;=========================================================================*/

   auto unsigned long HighWord;

   if ((Provider = LocGetProvider()) == 0)
      return 0;
   This = WireInfo + (Provider - 1);
   This->WireFile = CreateFile(FileName,                /* Full pathname of file */
     GENERIC_READ | GENERIC_WRITE,                      /* Allow read and write access */
     FILE_SHARE_READ | FILE_SHARE_WRITE,                /* No sharing by others */
     NULL,                                              /* No security atributes */
     OPEN_EXISTING,                                     /* File must already exist */
     FILE_ATTRIBUTE_NORMAL,                             /* No special attributes */
     0                                                  /* No template file */
     );
   if (This->WireFile == INVALID_HANDLE_VALUE)
      return 1;
   This->WireFileSize = GetFileSize(This->WireFile, &HighWord);
   This->WireMapping = CreateFileMapping(This->WireFile, /* Handle of file to be mapped */
     NULL,                                              /* No security attributes */
     PAGE_READWRITE,                                    /* Allow read and write access */
     0, This->WireFileSize,                             /* 64-bit size of mapping area */
     NULL                                               /* No name for mapping */
     );
   if (This->WireMapping == NULL)
   {
      WireProviders[Provider] = 0;
      return 0;
   }
   This->WireArea = MapViewOfFile(This->WireMapping,    /* Handle of mapping object */
     FILE_MAP_WRITE,                                    /* Allow read and write access */
     0, 0,                                              /* 64-bit offset into file */
     This->WireFileSize                                 /* Size of view */
     );
   if (This->WireArea == NULL)
   {
      WireProviders[Provider] = 0;
      return 0;
   }
   return Provider;

#else /*=======================================================================
   ;
   ; INDEPENDENT: In the independent version, we open the file, get its size,
   ; allocate memory, and then read in the file.
   ;
   ;=========================================================================*/

   if ((Provider = LocGetProvider()) == 0)
      return 0;
   This = WireInfo + (Provider - 1);
   This->WireFile = rdobinf(FileName);
   if (binopner(This->WireFile))
   {
      WireProviders[Provider] = 0;
      return 0;
   }
   sizbinf(This->WireFile);
   This->WireFileSize = posbinf(This->WireFile);
   This->WireArea = (UBYTE *) (getmem(This->WireFileSize+10));
   if (This->WireArea == NULL)
   {
      WireProviders[Provider] = 0;
      return 0;
   }
   rewbinf(This->WireFile);
   rdbinf(This->WireFile, This->WireArea, This->WireFileSize);
   This->WriteStatus = 0;
   return Provider;
#endif
}

/*
;/doc/ ***********************************************************************
;
; name of routine: TopsWireCreate -- Create a Wired File
;
; This function creates a file backing up a memory area and prepares it for
; use. To use this function the final size of the file must be know before
; this function can be called to create it.
;
; calling parameters:
;
; char* FileName     The name of file whose content is to be created.
;
; ULONG FileSize     The size of the file and memory area to be created.
;
; return parameters:
;
; The function returns the number of the provider supplied for this instance;
; or a zero if there is a problem.
;
; useful notes and assumptions: None
;
;/doc/ ************************************************************************
*/
#ifdef FPROTOTYPE
int TopsWireCreate(char *FileName, ULONG FileSize)
#else
int TopsWireCreate(FileName, FileSize)
   char *FileName;
   ULONG FileSize;
#endif
{
   auto int Provider;
   auto tWireInfo *This;

#ifdef MSCPLAT
   if ((Provider = LocGetProvider()) == 0)
      return 0;
   This = WireInfo + (Provider - 1);
   This->WireFile = CreateFile(FileName,                /* Full pathname of file */
     GENERIC_READ | GENERIC_WRITE,                      /* Allow read and write access */
     FILE_SHARE_READ | FILE_SHARE_WRITE,                /* No sharing by others */
     NULL,                                              /* No security atributes */
     CREATE_ALWAYS,                                     /* Create regardless of existance */
     FILE_ATTRIBUTE_NORMAL,                             /* No special attributes */
     0                                                  /* No template file */
     );
   if (This->WireFile == INVALID_HANDLE_VALUE)
   {
      WireProviders[Provider] = 0;
      return 0;
   }
   This->WireFileSize = FileSize;
   This->WireMapping = CreateFileMapping(This->WireFile, /* Handle of file to be mapped */
     NULL,                                              /* No security attributes */
     PAGE_READWRITE,                                    /* Allow read and write access */
     0, This->WireFileSize,                             /* 64-bit size of mapping area */
     NULL                                               /* No name for mapping */
     );
   if (This->WireMapping == NULL)
   {
      WireProviders[Provider] = 0;
      return 0;
   }
   This->WireArea = MapViewOfFile(This->WireMapping,    /* Handle of mapping object */
     FILE_MAP_WRITE,                                    /* Allow read and write access */
     0, 0,                                              /* 64-bit offset into file */
     This->WireFileSize                                 /* Size of view */
     );
   if (This->WireArea == NULL)
   {
      WireProviders[Provider] = 0;
      return 0;
   }
   return Provider;
#else
   if ((Provider = LocGetProvider()) == 0)
      return 0;
   This = WireInfo + (Provider - 1);
   This->WireFile = inibinf(FileName);
   if (binopner(This->WireFile))
   {
      WireProviders[Provider] = 0;
      return 0;
   }
   This->WireFileSize = FileSize;
   This->WireArea = (UBYTE *) (getmem(This->WireFileSize));
   if (This->WireArea == NULL)
   {
      WireProviders[Provider] = 0;
      return 0;
   }
   This->WriteStatus = 1;
   return Provider;
#endif
}

/*
;/doc/ ***********************************************************************
;
; name of routine: TopsWireClose -- Close a Wired File
;
; This function closes the file backing up the wired memory area and returns
; any resources used by it to the operating system. 
;
; calling parameters:
;
; int  Provider    The number of the wiring provider as returned by either
;                  TopsWireOpen or TopsWireCreate.
;
; return parameters:
;
; The function returns a zero if all went well; else it returns a nonzero
; error code.
;
; useful notes and assumptions: None
;
;/doc/ ************************************************************************
*/
#ifdef FPROTOTYPE
int TopsWireClose(int Provider)
#else
int TopsWireClose(Provider)
   int Provider;
#endif
{
   auto tWireInfo *This;

#ifdef MSCPLAT
   This = WireInfo + (Provider - 1);
   WireProviders[Provider] = 0;
   if (This->WireFile == INVALID_HANDLE_VALUE)
      return 1;
   UnmapViewOfFile(This->WireArea);
   CloseHandle(This->WireMapping);
   CloseHandle(This->WireFile);
   This->WireFile = INVALID_HANDLE_VALUE;
   return 0;
#else
   This = WireInfo + (Provider - 1);
   WireProviders[Provider] = 0;
   if (This->WireFile == nulbinf)
      return 1;
   if (This->WriteStatus)
   {
      rewbinf(This->WireFile);
      wrbinf(This->WireFile, This->WireArea, This->WireFileSize);
   }
   clsbinf(This->WireFile);
   free(This->WireArea);
   This->WireFile = nulbinf;
   return 0;
#endif
}

/*
;/doc/ ***********************************************************************
;
; name of routine: TopsWireAddress -- Get Address of Wired Memory Area
;
; Once a wired memory area has been created and opened, this funtion returns
; a pointer to the actual process local memory area. 
;
; calling parameters:
;
; int  Provider    The number of the wiring provider as returned by either
;                  TopsWireOpen or TopsWireCreate.
;
; return parameters:
;
; The function returns a pointer the memory area assigned to the file.
;
; useful notes and assumptions: None
;
;/doc/ ************************************************************************
*/
#ifdef FPROTOTYPE
UBYTE *TopsWireAddress(int Provider)
#else
UBYTE *TopsWireAddress(Provider)
   int Provider;
#endif
{
   auto tWireInfo *This;

   This = WireInfo + (Provider - 1);
   return This->WireArea;
}

/*
;/doc/ ***********************************************************************
;
; name of routine: TopsWireSize -- Get Size of Wired File
;
; Once a wired  memory area has been created and opened, this funtion returns
; the maximum size available for an actual memory area. 
;
; calling parameters:
;
; int  Provider    The number of the wiring provider as returned by either
;                  TopsWireOpen or TopsWireCreate.
;
; return parameters:
;
; The function returns the total size of the memory area.
;
; useful notes and assumptions: None
;
;/doc/ ************************************************************************
*/
#ifdef FPROTOTYPE
ULONG TopsWireSize(int Provider)
#else
ULONG TopsWireSize(Provider)
   int Provider;
#endif
{
   auto tWireInfo *This;

   This = WireInfo + (Provider - 1);
   return This->WireFileSize;
}
