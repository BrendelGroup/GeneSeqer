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
   ; name of module: Platform.h -- Define platforms Hosting Code
   ;
   ; useful notes and assumptions:
   ;
   ; This header file contains all the symbols and typedefs needed to work with the
   ; the PROMULA systems on a variety of platforms. The assumption in the design of
   ; this header is that each module written within the PROMULA system includes
   ; this header using the notation:
   ;
   ;    #include <platform.h>
   ;
   ; Preceeding this include various symbols are defined specifying the particular
   ; types of functions needed by the module. Using this technique all platform
   ; dependent information can be inserted into this header file. The function
   ; codes never change when the code moves from one platform to another. In
   ; particular, platform specific "ifdefs" are never placed in the code.
   ; The particular symbols recognized by this header are as follows:
   ;
   ; Symbol     Description of use
   ; ------     ------------------
   ; VIQASSRT   Gives access to assert and verify functionality.
   ; CHARTYPE   Gives access to the character type flags
   ; DATETIME   Gives access to date and time functions
   ; ECVTFUNC   Gives access to E-format style double precision conversion
   ; LONGMEMO   Gives access to the PROMULA long memory functions.
   ; MATHFUNC   Gives access to the mathematical functions.
   ; MEMOFUNC   Gives access to the memory allocation function
   ; OSYSFUNC   Gives access to operating system functions
   ; PACKFUNC   Gives access to the fixed point packing functions
   ; PRIOFUNC   Gives access to the PROMULA I/O functions
   ; RAWMFUNC   Gives access to the raw memory functions
   ; STRGFUNC   Gives access to the string functions
   ; THRDFUNC   Gives access to thread related functions
   ; UBIOFUNC   Gives access to raw unbuffered I/O functions
   ; VARARGUS   Gives access to variable arguments
   ;
   ; When compiling a symbol identifying the target platform must be defined on
   ; the command for the C compiler being used. The platform identifiers
   ; recognized by this version are as follows:
   ;
   ; 01/16/00 written by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : as part of the Promula Version 8 TOPS development
   ;
   ;/hdr/ ************************************************************************
 */

#ifdef MSCWIND			/* 32 bit MicroSoft C using windows */
#define MSCPLAT
#define USEGUITXT		/* Use Gui for text I/O */
#endif

#ifdef MSCPLAT			/* 32 bit MicroSoft C under WIN32S */
#define A86PLAT
#endif /* MSCPLAT */

#ifdef MPEPLAT
#define AN4PLAT
#endif

#ifdef ZORPLAT			/* Zortech C under DOS extented */
#define OSMSDOS			/* This platform operates under MS-DOS */
#define ANSPLAT			/* This platform supports ANSI C */
#define FPROTOTYPE		/* This platform uses function prototypes */
#define LITTLENDER   1		/* This platform does have little-ender sex */
#define CRUNTYPE     1		/* Runtime host C is ANSI4 */
#define NEEDSPACKF		/* This platform needs packing functions */
#define HASLNGFN		/* Uses long memory as explicit functions */
typedef int LONG;		/* Signed 4-byte integers */
typedef unsigned int ULONG;	/* Unsigned 4-byte integers */
typedef unsigned char UBYTE;	/* Unsigned 1-byte integer */
#endif /* ZORPLAT */

#ifdef UNXPLAT			/* A generic UNIX platform */
#define OSUNIX			/* This platform operates under UNIX */
#define NEEDMMOVE		/* This platform needs an overlaying move */
#define LITTLENDER   0		/* This platform does not have little-ender sex */
#define NEEDSPACKF		/* This platform needs packing functions */
#define VARUNIXV		/* Use Unix variable number of parameters */
#define CRUNTYPE     5		/* Runtime host C is KR4 */
#define LONG    int		/* Signed 4-byte integers */
#define ULONG   unsigned int	/* Unsigned 4-byte integers */
#define UBYTE unsigned char	/* Unsigned 1-byte integer */
#endif /* UNXPLAT */

#ifdef AIXPLAT			/* An IBM Workstation with  AIX and xlc ANSI C */
#define OSUNIX			/* This platform operates under UNIX */
#define POSIX			/* This platform is Posix compliant */
#define _POSIX_SOURCE
#define _XOPEN_SOURCE
#define FPROTOTYPE		/* This platform uses function prototypes */
#define LITTLENDER   0		/* This platform does not have little-ender sex */
#define NEEDSPACKF		/* This platform needs packing functions */
#define ANSPLAT			/* This platform supports ANSI C */
#define CRUNTYPE      1		/* Runtime host C is ANSI4 */
typedef int LONG;		/* Signed 4-byte integers */
typedef unsigned int ULONG;	/* Unsigned 4-byte integers */
typedef unsigned char UBYTE;	/* Unsigned 1-byte integer */
#endif /* AIXPLAT */

#ifdef WTCPLAT			/* Watcom C under DOS extended */
#define A86PLAT
#endif

#ifdef BO2PLAT			/* Borland C under OS/2 */
#define A86PLAT
#endif

#ifdef A86PLAT			/* ANSI C on a DOS extended platform */
#define OSMSDOS			/* This platform operates under MS-DOS */
#define FPROTOTYPE		/* This platform uses function prototypes */
#define LITTLENDER   1		/* This platform does have little-ender sex */
#define ANSPLAT			/* This platform supports ANSI C */
#define CRUNTYPE     1		/* Runtime host C is ANSI4 */
#define LONG int		/* Signed 4-byte integers */
#define ULONG unsigned int	/* Unsigned 4-byte integers */
#define UBYTE unsigned char	/* Unsigned 1-byte integer */
#endif /* A86PLAT */

#ifdef KRCPLAT			/* 32 bit MicroSoft C WIN32S with K+R */
#pragma warning(disable:4131)
#define OSMSDOS			/* This platform operates under MS-DOS */
#define LITTLENDER   1		/* This platform does have little-ender sex */
#define CRUNTYPE     5		/* Runtime host C is ANSI4 */
#define LONG int		/* Signed 4-byte integers */
#define ULONG unsigned int	/* Unsigned 4-byte integers */
#define UBYTE unsigned char	/* Unsigned 1-byte integer */
#define A86PLAT
#endif /* KRCPLAT */

#ifdef AZTPLAT			/* An Apple Macintosh using Aztec C */
#define OSMAC			/* This platform operates under the MacIntosh */
#define SHORTMEMO		/* This platform has a short memory */
#define LITTLENDER   0		/* This platform does not have little-ender sex */
#define NEEDSPACKF		/* This platform needs packing functions */
#define VARUNIXV		/* Use Unix variable number of parameters */
#define CRUNTYPE     5		/* Runtime host C is KR4 */
#define LONG    int		/* Signed 4-byte integers */
#define ULONG   unsigned int	/* Unsigned 4-byte integers */
#define UBYTE unsigned char	/* Unsigned 1-byte integer */
#endif /* AZTPLAT */

#ifdef MPWPLAT			/* An Apple Macintosh the MPW shell with MPW C */
#define OSMAC			/* This platform operates under the MacIntosh */
#define SHORTMEMO		/* This platform has a short memory */
#define LITTLENDER   0		/* This platform does not have little-ender sex */
#define NEEDSPACKF		/* This platform needs packing functions */
#define CRUNTYPE     5		/* Runtime host C is KR4 */
#define LONG    int		/* Signed 4-byte integers */
#define ULONG   unsigned int	/* Unsigned 4-byte integers */
#define UBYTE unsigned char	/* Unsigned 1-byte integer */
#endif /* MPWPLAT */

#ifdef TCCPLAT			/* Turbo C for the IBM PC DOS 2+. */
#define OSMSDOS			/* This platform operates under MS-DOS */
#define SHORTMEMO		/* This platform has a short memory */
#define FPROTOTYPE		/* This platform uses function prototypes */
#define LITTLENDER   1		/* This platform does have little-ender sex */
#define CRUNTYPE     0		/* Runtime host C is ANSI2 */
typedef long LONG;		/* Signed 4-byte integers */
typedef unsigned long ULONG;	/* Unsigned 4-byte integers */
typedef unsigned char UBYTE;	/* Unsigned 1-byte integer */
#endif /* TCCPLAT */

#ifdef THCPLAT			/* Think C for the Apple Macintosh. */
#define OSMAC			/* This platform operates under the MacIntosh */
#define SHORTMEMO		/* This platform has a short memory */
#define LITTLENDER   0		/* This platform does have little-ender sex */
#define NEEDSPACKF		/* This platform needs packing functions */
#define CRUNTYPE     5		/* Runtime host C is KR4 */
#define LONG    int		/* Signed 4-byte integers */
#define ULONG   unsigned int	/* Unsigned 4-byte integers */
#define UBYTE unsigned char	/* Unsigned 1-byte integer */
#endif /* THCPLAT */

#ifdef TSOPLAT			/* IBM mainframe under TSO using SAS C */
#define OSMSDOS			/* This platform operates under MS-DOS */
#define LITTLENDER   1		/* This platform does have little-ender sex */
#define CRUNTYPE     1		/* Runtime host C is ANSI4 */
typedef int LONG;		/* Signed 4-byte integers */
typedef unsigned int ULONG;	/* Unsigned 4-byte integers */
typedef unsigned char UBYTE;	/* Unsigned 1-byte integer */
#endif /* TSOPLAT */

#ifdef VMSPLAT			/* A VAX under VMS using standard VAX C */
#define OSVMS			/* This platform operates under VMS */
#define LITTLENDER   1		/* This platform does have little-ender sex */
#define NUCLEUS      1		/* This platform supports NUCLEUS dbs */
#define VARUNIXV		/* Use Unix variable number of parameters */
#define CRUNTYPE     3		/* Runtime host C is VMS */
typedef int LONG;		/* Signed 4-byte integers */
typedef unsigned int ULONG;	/* Unsigned 4-byte integers */
typedef unsigned char UBYTE;	/* Unsigned 1-byte integer */
#endif /* VMSPLAT */

#ifdef SKYPLAT			/* A SKYBOLT board with Sky High C */
#define OSUNIX			/* This platform operates under UNIX */
#define FPROTOTYPE		/* This platform uses function prototypes */
#define LITTLENDER   0		/* This platform does not have little-ender sex */
#define NEEDSPACKF		/* This platform needs packing functions */
#define CRUNTYPE     6		/* Runtime host C is I860 */
typedef int LONG;		/* Signed 4-byte integers */
typedef unsigned int ULONG;	/* Unsigned 4-byte integers */
typedef unsigned char UBYTE;	/* Unsigned 1-byte integer */
#endif /*SKYPLAT */

#ifdef SN4PLAT			/* A SUN 4 or SUN Sparcstation under SUNOS */
#define OSUNIX			/* This platform operates under UNIX */
#define NEEDMMOVE		/* This platform needs an overlaying move */
#define LITTLENDER   0		/* This platform does not have little-ender sex */
#define NEEDSPACKF		/* This platform needs packing functions */
#define VARUNIXV		/* Use Unix variable number of parameters */
#define UNXPLAT			/* This platform is generic UNIX otherwise */
#define TERMCAP			/* This platform uses the termcap facility */
#define CRUNTYPE     5		/* Runtime host C is KR4 */
#define LONG    int		/* Signed 4-byte integers */
#define ULONG   unsigned int	/* Unsigned 4-byte integers */
#define UBYTE unsigned char	/* Unsigned 1-byte integer */
#endif /*SN4PLAT */

#ifdef FUJPLAT			/* A Fujitsu DS station */
#define OSUNIX			/* This platform operates under UNIX */
#define FPROTOTYPE		/* This platform uses function prototypes */
#define LITTLENDER   0		/* This platform does not have little-ender sex */
#define NEEDSPACKF		/* This platform needs packing functions */
#define ANSPLAT			/* This platform supports ANSI C */
#define CRUNTYPE      1		/* Runtime host C is ANSI4 */
typedef int LONG;		/* Signed 4-byte integers */
typedef unsigned int ULONG;	/* Unsigned 4-byte integers */
typedef unsigned char UBYTE;	/* Unsigned 1-byte integer */
#endif /*FUJPLAT */

#ifdef AN6PLAT			/* A ANSI C compiler on a 386 UNIX platform */
#define OSUNIX			/* This platform operates under UNIX */
#define FPROTOTYPE		/* This platform uses function prototypes */
#define LITTLENDER   1		/* This platform does have little-ender sex */
#define ANSPLAT			/* This platform supports ANSI C */
#define CRUNTYPE      1		/* Runtime host C is ANSI4 */
#define POSIX			/* This platform is Posix compliant */
typedef int LONG;		/* Signed 4-byte integers */
typedef unsigned int ULONG;	/* Unsigned 4-byte integers */
typedef unsigned char UBYTE;	/* Unsigned 1-byte integer */
#define UNXPLAT			/* This platform is generic UNIX otherwise */
#endif /*AN6PLAT */

#ifdef UCSPLAT			/* A Unisys ANSI */
#define FPROTOTYPE		/* This platform uses function prototypes */
#define LITTLENDER   1		/* This platform does have little-ender sex */
#define ANSPLAT			/* This platform supports ANSI C */
#define NEEDSPACKF		/* This platform needs packing functions */
#define CRUNTYPE      1		/* Runtime host C is ANSI4 */
typedef int LONG;		/* Signed 4-byte integers */
typedef unsigned int ULONG;	/* Unsigned 4-byte integers */
typedef unsigned char UBYTE;	/* Unsigned 1-byte integer */
#endif /*UCSPLAT */

#ifdef ALFPLAT			/* DEC Alpha 64-bit computer */
#define AN4PLAT			/* Do it just like AN4PLAT */
#endif /* ALFPLAT */

#ifdef AN4PLAT			/* A ANSI C compiler on a sun 4 */
#define OSUNIX			/* This platform operates under UNIX */
#define LITTLENDER   0		/* This platform does not have little-ender sex */
#define NEEDSPACKF		/* This platform needs packing functions */
#define FPROTOTYPE		/* This platform uses function prototypes */
#define ANSPLAT			/* This platform supports ANSI C */
#define CRUNTYPE     1		/* Runtime host C is ANSI4 */
typedef int LONG;		/* Signed 4-byte integers */
typedef unsigned int ULONG;	/* Unsigned 4-byte integers */
typedef unsigned char UBYTE;	/* Unsigned 1-byte integer */
#endif /*AN4PLAT */

#ifdef HUXPLAT			/* HP 9000 under Unix 5.3 */
#define OSUNIX			/* This platform operates under UNIX */
#define LITTLENDER   0		/* This platform does not have little-ender sex */
#define NEEDSPACKF		/* This platform needs packing functions */
#define VARUNIXV		/* Use Unix variable number of parameters */
#define UNXPLAT			/* This platform is generic UNIX otherwise */
#define CRUNTYPE     5		/* Runtime host C is KR4 */
#define LONG    int		/* Signed 4-byte integers */
#define ULONG   unsigned int	/* Unsigned 4-byte integers */
#define UBYTE unsigned char	/* Unsigned 1-byte integer */
#endif /*HUXPLAT */

#ifdef U86PLAT			/* Unix system an a 386 workstation */
#define OSUNIX			/* This platform operates under UNIX */
#define LITTLENDER   1		/* This platform does have little-ender sex */
#define NEEDMMOVE		/* This platform needs an overlaying move */
#define VARUNIXV		/* Use Unix variable number of parameters */
#define UNXPLAT			/* This platform is generic UNIX otherwise */
#define CRUNTYPE     5		/* Runtime host C is KR4 */
#define LONG    int		/* Signed 4-byte integers */
#define ULONG   unsigned int	/* Unsigned 4-byte integers */
#define UBYTE unsigned char	/* Unsigned 1-byte integer */
#endif /*U86PLAT */

#ifdef DGEPLAT			/* A Data General with an ANSI C compiler */
#define NEEDMMOVE		/* This platform needs an overlaying move */
#define NEEDSPACKF		/* This platform needs packing functions */
#define ANSI
#define FPROTOTYPE		/* This platform uses function prototypes */
#define LITTLENDER   0		/* This platform does not have little-ender sex */
#define ANSPLAT			/* This platform supports ANSI C */
#define CRUNTYPE      1		/* Runtime host C is ANSI4 */
typedef int LONG;		/* Signed 4-byte integers */
typedef unsigned int ULONG;	/* Unsigned 4-byte integers */
typedef unsigned char UBYTE;	/* Unsigned 1-byte integer */
#endif /* DGEPLAT */


#ifdef LCCPLAT			/* Lattice C for the IBM PC DOS 2+ */
#define OSMSDOS			/* This platform operates under MS-DOS */
#define SHORTMEMO		/* This platform has a short memory */
#define FPROTOTYPE		/* This platform uses function prototypes */
#define LITTLENDER   1		/* This platform does have little-ender sex */
#define CRUNTYPE     0		/* Runtime host C is ANSI2 */
typedef long LONG;		/* Signed 4-byte integers */
typedef unsigned long ULONG;	/* Unsigned 4-byte integers */
typedef unsigned char UBYTE;	/* Unsigned 1-byte integer */
#endif /* LCCPLAT */

#ifdef COHPLAT			/* COHERENT 386 from Mark Williams Company */
#define OSUNIX			/* This platform operates under UNIX */
#define LITTLENDER   0		/* This platform does not have little-ender sex */
#define NEEDMMOVE		/* This platform needs an overlaying move */
#define VARUNIXV		/* Use Unix variable number of parameters */
#define UNXPLAT			/* This platform is generic UNIX otherwise */
#define CRUNTYPE     5		/* Runtime host C is KR4 */
#endif /*COHPLAT */

#include <stdio.h>

#ifndef TRUE
#define TRUE   1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#define EPSILON_FUZZ 1.0e-35	/* Epsilon fuzz value */

#ifdef FPROTOTYPE
#define VOID void
#define HANDLE void*
#define CONST const
#else
#define VOID char
#define HANDLE char*
#define CONST
#endif

#ifdef CHARTYPE
#define CLETTER         1
#define CUPPERC         2
#define CLOWERC         4
#define CIDENTC         8
#define CNUMBER        16
#define CSYMBOL        32
#define CWSPACE        64
#define CQUOTEC       128	/* Indicates a quote character type */
#define NULLCHAR      128
#define BACKSLASH     129	/* Metacharacter code for a backslash */
#endif /* CHARTYPE */

#ifdef DATETIME
#ifdef AZTPLAT
#include <utime.h>
#else
#include <time.h>
#ifdef UNXPLAT
#include <sys/types.h>
#include <sys/timeb.h>
#endif /* UNXPLAT */
#endif /* AZTPLAT */
#endif /* DATETIME */

#ifdef ECVTFUNC
#ifdef DGEPLAT
#define DIGPREC     18		/* max # of significant digits */
#endif
#ifdef UNXPLAT
#define DIGPREC     18		/* max # of significant digits */
#endif /* UNXPLAT */
#ifdef AIXPLAT
#define HIGHEXPO		/* has high precision exponent */
#define DIGPREC     18		/* max # of significant digits */
#endif /* AIXPLAT */
#ifdef FUJPLAT
#define HIGHEXPO		/* has high precision exponent */
#define DIGPREC     18		/* max # of significant digits */
#endif /* FUJPLAT */
#ifdef A86PLAT
#define HIGHEXPO		/* has high precision exponent */
#define DIGPREC     18		/* max # of significant digits */
#endif /* A86PLAT */
#ifdef AN4PLAT
#define HIGHEXPO		/* has high precision exponent */
#define DIGPREC     18		/* max # of significant digits */
#endif /* AN4PLAT */
#ifdef UCSPLAT
#define HIGHEXPO		/* has high precision exponent */
#define DIGPREC     18		/* max # of significant digits */
#endif /* UCSPLAT */
#ifdef AZTPLAT
#define DIGPREC     17		/* max # of significant digits */
#endif /* AZTPLAT */
#ifdef MPWPLAT
#define DIGPREC     17		/* max # of significant digits */
#endif /* MPWPLAT */
#ifdef TCCPLAT
#include <string.h>		/* declarations for "memset+memcpy" */
#include <stdlib.h>		/* declaration for "ecvt" */
#define HASECVTF		/* platform has ecvt */
#endif /* TCCPLAT */
#ifdef THCPLAT
#define DIGPREC     17		/* max # of significant digits */
#endif /* THCPLAT */
#ifdef TSOPLAT
#define DIGPREC     17		/* max # of significant digits */
#endif /* TSOPLAT */
#ifdef VMSPLAT
#include <string.h>		/* declarations for "memset+memcpy" */
#include <stdlib.h>		/* declaration for "ecvt" */
#define HASECVTF		/* platform has ecvt */
#endif /* VMSPLAT */
#ifdef SKYPLAT
#define HIGHEXPO		/* has high precision exponent */
#define DIGPREC     18		/* max # of significant digits */
#endif /* SKYPLAT */
#ifdef LCCPLAT
#include <string.h>		/* declarations for "memset+memcpy" */
#include <math.h>		/* declaration for "ecvt" */
#define HASECVTF		/* platform has ecvt */
#endif /* LCCPLAT */
#endif /* ECVTFUNC */

#ifdef LONGMEMO
#ifdef HASLNGFN
LONG *defptr (unsigned int);
LONG *prcptr (unsigned int);
LONG *valptr (unsigned int);
unsigned int defadr (LONG *);
unsigned int prcadr (LONG *);
#else
extern LONG *defn;		/* Definition storage area */
extern LONG *proc;		/* Procedure storage area */
extern LONG *valu;		/* Value storage area */
#define defptr(x) (defn + x)
#define prcptr(x) (proc + x)
#define valptr(x) (valu + x)
#ifdef THCPLAT
#define defadr(x) ((unsigned int)(x - defn))
#define prcadr(x) ((unsigned int)(x - proc))
#else
#define defadr(x) (x - defn)
#define prcadr(x) (x - proc)
#endif /* THCPLAT */
#endif /* HASLNGFN */
#endif /* LONGMEMO */

#ifdef MATHFUNC
#ifdef MPWPLAT
#define LONGDBL
#endif /* MPWPLAT */
#include <math.h>
#ifdef DGEPLAT
#define BIGDBLE   1.0e75
#else
#ifdef VMSPLAT
#define BIGDBLE   1.0e38
#else
#define BIGDBLE   1.0e300
#endif /* VMSPLAT */
#endif /* DGEPLAT */
#endif /* MATHFUNC */

#ifdef OSYSFUNC
#ifdef UNXPLAT
extern void exit ();
extern char *getenv ();
#else
#ifdef AZTPLAT
extern void exit ();
extern char *getenv ();
#else
#ifdef VMSPLAT
extern void exit ();
extern char *getenv ();
#else
#include <setjmp.h>
#include <stdlib.h>
#endif /* VMSPLAT */
#endif /* AZTPLAT */
#endif /* UNXPLAT */
#endif /* OSYSFUNC */

#ifdef STRGFUNC
#include <string.h>
#ifdef THCPLAT
#define cmpstrn(x,y,z) memcmp(x,y,(size_t)z)
#else
#define cmpstrn  memcmp
#endif
#endif /* STRGFUNC */

#ifdef RAWMFUNC
#include <string.h>
#ifdef THCPLAT
#define cpymem(x,y,z) memcpy(y,x,(size_t)z)
#define cpyovl(x,y,z) memmove(y,x,(size_t)z)
#define filmem(x,y,z) memset(x,z,(size_t)y)
#define zeromem(p,x)  memset(p,0,(size_t)x)
#else
#define filmem(x,y,z) memset(x,z,y)
#define cpymem(x,y,z) memcpy(y,x,z)
#ifdef NEEDMMOVE
extern void cpyovl ();
#else
#define cpyovl(x,y,z) memmove(y,x,z)
#define zeromem(p,n)  memset(p,0,n)
#endif /* NEEDMMOVE */
#endif /* THCPLAT */
#endif /* RAWMFUNC */

#ifdef PACKFUNC
#ifdef MPWPLAT
#define longdbl   4		/* Number of longs per double */
#else
#define longdbl   2		/* Number of longs per double */
#endif /* MPWPLAT */
#ifdef NEEDSPACKF
#define getuns16(x,y) get16u(&(x),y)
#define getsng16(x,y) get16s(&(x),y)
#define getsng32(x,y) get32(&(x),y)
#define getuns32(x,y) get32u(&(x),y)
#define getdble(x,y) get64(&(x),y)
#ifdef FPROTOTYPE
void get16u (unsigned int *x, unsigned char *y);
void get16s (int *x, unsigned char *y);
void putuns16 (unsigned char *x, unsigned int y);
void putsng16 (unsigned char *x, int y);
void get32 (LONG * x, unsigned char *y);
void get32u (ULONG * x, unsigned char *y);
void putsng32 (unsigned char *x, LONG y);
void putuns32 (unsigned char *x, unsigned int y);
void get64 (double *x, unsigned char *y);
void putdble (unsigned char *x, double y);
#else
extern void get16u ();
extern void get16s ();
extern void putuns16 ();
extern void putsng16 ();
extern void get32 ();
extern void get32u ();
extern void putsng32 ();
extern void putuns32 ();
extern void get64 ();
extern void putdble ();
#endif /* FPROTOTYPE */
#else
#define getuns16(x,y)  x = *(unsigned short*)(y)
#define putuns16(x,y)  *(unsigned short*)(x) = (unsigned short) y
#define getsng16(x,y)  x = *(short*)(y)
#define putsng16(x,y)  *(short*)(x) = (short) y
#define getsng32(x,y)  x = *(LONG*)(y)
#define getuns32(x,y)  x = *(ULONG*)(y)
#define putsng32(x,y)  *(LONG*)(x) = y
#define putuns32(x,y)  *(ULONG*)(x) = (ULONG) y
#define getdble(x,y)   x = *(double*)(y)
#define putdble(x,y)   *(double*)(x) = y
#define getptr(x,y)   x = *(void**)(y)
#define putptr(x,y)   *(void**)(x) = y
#define valuns16(y)  *(unsigned short*)(y)
#define valsng16(y)  *(short*)(y)
#define valsng32(y)  *(LONG*)(y)
#define valdble(y)   *(double*)(y)
#endif /* NEEDSPACKF */
#endif /* PACKFUNC */

#ifdef MEMOFUNC
#ifdef UNXPLAT
#ifndef COHPLAT
#include <malloc.h>
#endif
#else
#ifdef AZTPLAT
extern char *malloc ();
#else
#ifdef COHPLAT
extern char *malloc ();
#else
#include <stdlib.h>
#endif /* COHPLAT */
#endif /* AZTPLAT */
#endif /* UNXPLAT */
#ifdef THCPLAT
#define getmem(x) malloc((size_t)x)
#define xgetmem(x) malloc((size_t)x)
#define xfree(x) free(x)
#else
#define getmem(x) malloc(x)
#define xgetmem(x) malloc(x)
#define xfree(x) free(x)
#endif /* THCPLAT */
#endif /* MEMOFUNC */

#ifdef PRIOFUNC
#ifdef WTCPLAT
#include <io.h>
#include <fcntl.h>
#include <sys\stat.h>
#define EOLNCHAR 2		/* End-of-line number of characters */
#define setbufsz(x,y,z)  setvbuf(x,y,_IOFBF,z)
#define fileptr void*		/* Longest basic type needed for handle */
#define txtfile FILE*		/* The type of a text file */
#define binfile int		/* The type of raw binary file handles */
#define SINPUT  stdin		/* The standard input file */
#define SOUTPUT stdout		/* The standard output file */
#define SCONSOL stderr		/* The console file */
#define nultxtf NULL		/* Indicates an unassigned text file */
#define initxtf(x) fopen(x,"w+")	/* Initialize and open a text file */
#define apptxtf(x) fopen(x,"a")	/* Initialize or append to a text file */
#define opntxtf(x) fopen(x,"r+")	/* Open an existing text file */
#define rdotxtf(x) fopen(x,"r")	/* Open a readonly text file */
#ifdef USEGUITXT
extern void wrtxtf ();
extern void wrteol ();
extern void wrtffd ();
extern void *rdtxtf ();
#else
#define wrtxtf(f,r,n) fwrite(r,1,n,f)	/* Write a coded record */
#define wrteol(x)  fputc('\n',x)	/* Write an end-of-line */
#define wrtffd(x)  fputc('\014',x)	/* Write a formfeed */
#define rdtxtf(x,y,z) fgets(y,z,x)	/* Read a text file */
#endif
#define clstxtf(x) fclose(x)	/* Close a text file */
#define puttxtf(x,y) fseek(x,y,SEEK_SET)	/* Position a text file */
#define rewtxtf(x) fseek(x,0L,SEEK_SET)		/* Rewind a text file */
#define siztxtf(x) fseek(x,0L,SEEK_END)		/* Position to end (size) of a text file */
#define postxtf(x) ftell(x)	/* Return the position of a text file */
#define txtopner(x) (x==NULL)	/* Text open error ? */
#define nulbinf    0		/* Indicates an unassigned binary file */
#define inibinf(x) open(x,O_CREAT|O_TRUNC|O_RDWR|O_BINARY,S_IREAD|S_IWRITE)
#define opnbinf(x) open(x,O_RDWR|O_BINARY)
#define rdobinf(x) open(x,O_RDONLY|O_BINARY)
#define rdbinf(f,r,n) read(f,r,n)
#define wrbinf(f,r,n) write(f,r,n)
#define rdlong(f,r,n) read(f,r,n*4)
#define wrlong(f,r,n) write(f,r,n*4)
#define clsbinf(x) close(x)
#define rewbinf(x) lseek(x,0L,SEEK_SET)
#define putbinf(x,y) lseek(x,y,SEEK_SET)
#define posbinf(x) lseek(x,0L,SEEK_CUR)
#define sizbinf(x) lseek(x,0L,SEEK_END)
#define binopner(x) (x<0)
#define delfile(x) remove(x)
#else
#ifdef BO2PLAT
#include <io.h>
#include <fcntl.h>
#include <sys\stat.h>
#define EOLNCHAR 2		/* End-of-line number of characters */
#define setbufsz(x,y,z)  setvbuf(x,y,_IOFBF,z)
#define fileptr void*		/* Longest basic type needed for handle */
#define txtfile FILE*		/* The type of a text file */
#define binfile int		/* The type of raw binary file handles */
#define SINPUT  stdin		/* The standard input file */
#define SOUTPUT stdout		/* The standard output file */
#define SCONSOL stderr		/* The console file */
#define nultxtf NULL		/* Indicates an unassigned text file */
#define initxtf(x) fopen(x,"w+")	/* Initialize and open a text file */
#define apptxtf(x) fopen(x,"a")	/* Initialize or append to a text file */
#define opntxtf(x) fopen(x,"r+")	/* Open an existing text file */
#define rdotxtf(x) fopen(x,"r")	/* Open a readonly text file */
#ifdef USEGUITXT
extern void wrtxtf ();
extern void wrteol ();
extern void wrtffd ();
extern void *rdtxtf ();
#else
#define wrtxtf(f,r,n) fwrite(r,1,n,f)	/* Write a coded record */
#define wrteol(x)  fputc('\n',x)	/* Write an end-of-line */
#define wrtffd(x)  fputc('\014',x)	/* Write a formfeed */
#define rdtxtf(x,y,z) fgets(y,z,x)	/* Read a text file */
#endif
#define clstxtf(x) fclose(x)	/* Close a text file */
#define puttxtf(x,y) fseek(x,y,SEEK_SET)	/* Position a text file */
#define rewtxtf(x) fseek(x,0L,SEEK_SET)		/* Rewind a text file */
#define siztxtf(x) fseek(x,0L,SEEK_END)		/* Position to end (size) of a text file */
#define postxtf(x) ftell(x)	/* Return the position of a text file */
#define txtopner(x) (x==NULL)	/* Text open error ? */
#define nulbinf    0		/* Indicates an unassigned binary file */
#define inibinf(x) open(x,O_CREAT|O_TRUNC|O_RDWR|O_BINARY,S_IREAD|S_IWRITE)
#define opnbinf(x) open(x,O_RDWR|O_BINARY)
#define rdobinf(x) open(x,O_RDONLY|O_BINARY)
#define rdbinf(f,r,n) read(f,r,n)
#define wrbinf(f,r,n) write(f,r,n)
#define rdlong(f,r,n) read(f,r,n*4)
#define wrlong(f,r,n) write(f,r,n*4)
#define clsbinf(x) close(x)
#define rewbinf(x) lseek(x,0L,SEEK_SET)
#define putbinf(x,y) lseek(x,y,SEEK_SET)
#define posbinf(x) lseek(x,0L,SEEK_CUR)
#define sizbinf(x) lseek(x,0L,SEEK_END)
#define binopner(x) (x<0)
#define delfile(x) remove(x)
#else
#ifdef UNXPLAT
#define EOLNCHAR 1		/* End-of-line number of characters */
#define SEEK_SET  0
#define SEEK_END  2
#endif /* UNXPLAT */
#ifdef AIXPLAT
#define EOLNCHAR 1		/* End-of-line number of characters */
#define setbufsz(x,y,z)  setvbuf(x,y,_IOFBF,z)
#endif /* AIXPLAT */
#ifdef FUJPLAT
#define EOLNCHAR 1		/* End-of-line number of characters */
#define setbufsz(x,y,z)  setvbuf(x,y,_IOFBF,z)
#endif /* FUJPLAT */
#ifdef AN4PLAT
#define EOLNCHAR 1		/* End-of-line number of characters */
#define setbufsz(x,y,z)  setvbuf(x,y,_IOFBF,z)
#endif /* AN4PLAT */
#ifdef A86PLAT
#define EOLNCHAR 2		/* End-of-line number of characters */
#define setbufsz(x,y,z)  setvbuf(x,y,_IOFBF,z)
#endif /* A86PLAT */
#ifdef AZTPLAT
#define EOLNCHAR 1		/* End-of-line number of characters */
#define SEEK_SET  0
#define SEEK_END  2
#endif /* AZTPLAT */
#ifdef MPWPLAT
#define EOLNCHAR 2		/* End-of-line number of characters */
#define SEEK_SET  0
#define SEEK_END  2
#endif /* MPWPLAT */
#ifdef TCCPLAT
#define EOLNCHAR 2		/* End-of-line number of characters */
#endif /* TCCPLAT */
#ifdef THCPLAT
#define EOLNCHAR 2		/* End-of-line number of characters */
#endif /* THCPLAT */
#ifdef TSOPLAT
#define EOLNCHAR 2		/* End-of-line number of characters */
#define SEEK_SET  0
#define SEEK_END  2
#endif /* TSOPLAT */
#ifdef VMSPLAT
#define EOLNCHAR 1		/* End-of-line number of characters */
#endif /* VMSPLAT */
#ifdef SKYPLAT
#define EOLNCHAR 1		/* End-of-line number of characters */
#endif /* SKYPLAT */
#ifdef LCCPLAT
#define EOLNCHAR 2		/* End-of-line number of characters */
#endif /* LCCPLAT */
#define fileptr void*		/* Longest basic type needed for handle */
#define txtfile FILE*		/* The type of a text file */
#define binfile FILE*		/* The type of raw binary file handles */
#define SINPUT  stdin		/* The standard input file */
#define SOUTPUT stdout		/* The standard output file */
#define SCONSOL stderr		/* The console file */
#define nultxtf NULL		/* Indicates an unassigned text file */
#define initxtf(x) fopen(x,"w+")	/* Initialize and open a text file */
#define apptxtf(x) fopen(x,"a")	/* Initialize or append to a text file */
#define opntxtf(x) fopen(x,"r+")	/* Open an existing text file */
#define rdotxtf(x) fopen(x,"r")	/* Open a readonly text file */
#ifdef USEGUITXT
extern void wrtxtf ();
extern void wrteol ();
extern void wrtffd ();
extern void *rdtxtf ();
#else
#define rdtxtf(x,y,z) fgets(y,z,x)	/* Read a text file */
#define wrteol(x)  fputc('\n',x)	/* Write an end-of-line */
#define wrtffd(x)  fputc('\014',x)	/* Write a formfeed */
#ifdef THCPLAT
#define wrtxtf(f,r,n) fwrite(r,1L,(size_t)n,f)
#else
#ifdef VMSPLAT
#define wrtxtf(f,r,n) fwrite(r,n,1,f)	/* Write a coded record */
#else
#define wrtxtf(f,r,n) fwrite(r,1,n,f)	/* Write a coded record */
#endif /* VMSPLAT */
#endif /* THCPLAT */
#endif
#define clstxtf(x) fclose(x)	/* Close a text file */
#define puttxtf(x,y) fseek(x,y,SEEK_SET)	/* Position a text file */
#define rewtxtf(x) fseek(x,0L,SEEK_SET)		/* Rewind a text file */
#define siztxtf(x) fseek(x,0L,SEEK_END)		/* Position to end (size) of a text file */
#define postxtf(x) ftell(x)	/* Return the position of a text file */
#define txtopner(x) (x==NULL)	/* Text open error ? */
#define nulbinf    NULL		/* Indicates an unassigned binary file */
#define inibinf(x) fopen(x,"w+b")	/* Initialize and open a binary file */
#define opnbinf(x) fopen(x,"r+b")	/* Open an existing binary file */
#define rdobinf(x) fopen(x,"rb")	/* Open a readonly binary file */
#ifdef THCPLAT
#define rdbinf(f,r,n) fread(r,1L,(size_t)n,f)
#define wrbinf(f,r,n) fwrite(r,1L,(size_t)n,f)
#define rdlong(f,r,n) fread(r,(size_t)n,4L,f)
#define wrlong(f,r,n) fwrite(r,(size_t)n,4L,f)
#else
#define rdbinf(f,r,n) fread(r,1,n,f)	/* Read a binary file */
#define wrbinf(f,r,n) fwrite(r,1,n,f)	/* Write a binary file */
#define rdlong(f,r,n) fread(r,n,4,f)	/* Read longs from a binary file */
#define wrlong(f,r,n) fwrite(r,n,4,f)	/* Write longs to a binary file */
#endif /* THCPLAT */
#define clsbinf(x) fclose(x)	/* Close a binary file */
#define rewbinf(x) fseek(x,0L,SEEK_SET)		/* Rewind a binary file */
#define putbinf(x,y) fseek(x,y,SEEK_SET)	/* Position a binary file */
#define posbinf(x) ftell(x)	/* Position of a binary file */
#define sizbinf(x) fseek(x,0L,SEEK_END)		/* Size of a binary file */
#define binopner(x) (x==NULL)	/* Binary open error ? */
#ifdef VMSPLAT
#define delfile(x) delete(x)	/* Delete a file */
#else
#ifdef SKYPLAT
#define delfile(x) remove(x)	/* Delete a file */
#else
#ifdef AN4PLAT
#define delfile(x) remove(x)	/* Delete a file */
#else
#ifdef A86PLAT
#define delfile(x) remove(x)	/* Delete a file */
#else
#ifdef ZORPLAT
#define delfile(x) remove(x)	/* Delete a file */
#else
#ifdef ANSPLAT
#define delfile(x) remove(x)	/* Delete a file */
#else
#define delfile(x) unlink(x)	/* Delete a file */
#endif /* ANSPLAT */
#endif /* ZORPLAT */
#endif /* A86PLAT */
#endif /* AN4PLAT */
#endif /* SKYPLAT */
#endif /* VMSPLAT */
#endif /* BO2PLAT */
#endif /* WTCPLAT */
#endif /* PRIOFUNC */

#ifdef UBIOFUNC
#ifdef UNXPLAT
#define rawfile FILE*		/* The type of raw binary file handles */
#define nulrawf    NULL		/* Indicates an unassigned binary file */
#define inirawf(x) fopen(x,"w+b")	/* Initialize and open a binary file */
#define opnrawf(x) fopen(x,"r+b")	/* Open an existing binary file */
#define rdorawf(x) fopen(x,"rb")	/* Open an existing binary file */
#define rdrawf(f,r,n) fread(r,1,n,f)	/* Read a binary file */
#define wrrawf(f,r,n) fwrite(r,1,n,f)	/* Write a binary file */
#define clsrawf(x) fclose(x)	/* Close a binary file */
#define rewrawf(x) fseek(x,0L,0)	/* Rewind a binary file */
#define putrawf(x,y) fseek(x,y,0)	/* Position a binary file */
#define posrawf(x) ftell(x)	/* Position of a binary file */
#define sizrawf(x) fseek(x,0L,2)	/* Size of a binary file */
#define rawopner(x) (x==NULL)	/* Binary open error ? */
#endif /* UNXPLAT */
#ifdef ANSPLAT
#define rawfile FILE*		/* The type of raw binary file handles */
#define nulrawf    NULL		/* Indicates an unassigned binary file */
#define inirawf(x) fopen(x,"w+b")	/* Initialize and open a binary file */
#define opnrawf(x) fopen(x,"r+b")	/* Open an existing binary file */
#define rdorawf(x) fopen(x,"rb")	/* Open an existing binary file */
#define rdrawf(f,r,n) fread(r,1,n,f)	/* Read a binary file */
#define wrrawf(f,r,n) fwrite(r,1,n,f)	/* Write a binary file */
#define clsrawf(x) fclose(x)	/* Close a binary file */
#define rewrawf(x) fseek(x,0L,0)	/* Rewind a binary file */
#define putrawf(x,y) fseek(x,y,0)	/* Position a binary file */
#define posrawf(x) ftell(x)	/* Position of a binary file */
#define sizrawf(x) fseek(x,0L,2)	/* Size of a binary file */
#define rawopner(x) (x==NULL)	/* Binary open error ? */
#endif /* ANSPLAT */
#ifdef AZTPLAT
#include <fcntl.h>
extern long lseek ();
#define rawfile int		/* The type of raw binary file handles */
#define nulrawf 0		/* Indicates an unassigned binary file */
#define inirawf(x) open(x,O_CREAT|O_TRUNC|O_RDWR)	/* Initialize and open a binary file */
#define opnrawf(x) open(x,O_RDWR)	/* Open an existing binary file */
#define rdorawf(x) open(x,O_RDONLY)	/* Open an existing binary file */
#define rdrawf(f,r,n) read(f,r,n)	/* Read a binary file */
#define wrrawf(f,r,n) write(f,r,n)	/* Write a binary file */
#define clsrawf(x) close(x)	/* Close a binary file */
#define rewrawf(x) lseek(x,0L,0)	/* Rewind a binary file */
#define putrawf(x,y) lseek(x,y,0)	/* Position a binary file */
#define posrawf(x) lseek(x,0L,1)	/* Position of a binary file */
#define sizrawf(x) lseek(x,0L,2)	/* Size of a binary file */
#define rawopner(x) (x<0)	/* Binary open error ? */
#endif /* AZTPLAT */
#ifdef MPWPLAT
#define rawfile FILE*		/* The type of raw binary file handles */
#define nulrawf    NULL		/* Indicates an unassigned binary file */
#define inirawf(x) fopen(x,"w+b")	/* Initialize and open a binary file */
#define opnrawf(x) fopen(x,"r+b")	/* Open an existing binary file */
#define rdorawf(x) fopen(x,"rb")	/* Open an existing binary file */
#define rdrawf(f,r,n) fread(r,1,n,f)	/* Read a binary file */
#define wrrawf(f,r,n) fwrite(r,1,n,f)	/* Write a binary file */
#define clsrawf(x) fclose(x)	/* Close a binary file */
#define rewrawf(x) fseek(x,0L,0)	/* Rewind a binary file */
#define putrawf(x,y) fseek(x,y,0)	/* Position a binary file */
#define posrawf(x) ftell(x)	/* Position of a binary file */
#define sizrawf(x) fseek(x,0L,2)	/* Size of a binary file */
#define rawopner(x) (x==NULL)	/* Binary open error ? */
#endif /* MPWPLAT */
#ifdef TCCPLAT
#include <io.h>
#include <fcntl.h>		/* File control options for open() */
#include <sys\stat.h>		/* Defines structures for stat */
#define rawfile int		/* The type of raw binary file handles */
#define nulrawf    -1		/* Indicates an unassigned binary file */
#define inirawf(x) open(x,O_CREAT | O_TRUNC | O_RDWR | O_BINARY,S_IWRITE | S_IREAD)
#define opnrawf(x) open(x,O_RDWR | O_BINARY)
#define rdorawf(x) open(x,O_RDONLY | O_BINARY)
#define rdrawf(x,y,z)  read(x,y,z)	/* Read from a binary file */
#define wrrawf(x,y,z) write(x,y,z)	/* Write to a binary file */
#define clsrawf(x)  close(x)	/* Close a binary file */
#define rewrawf(x) lseek(x,0L,0)	/* Rewind a binary file */
#define putrawf(x,y) lseek(x,y,0)	/* Position a binary file */
#define posrawf(x) lseek(x,0L,1)	/* Return position of a binary file */
#define sizrawf(x) lseek(x,0L,2)	/* Size of a binary file */
#define rawopner(x) (x<0)	/* Binary open error ? */
#endif /* TCCPLAT */
#ifdef THCPLAT
#include <unix.h>
#include <fcntl.h>
#define rawfile int		/* The type of raw binary file handles */
#define nulrawf    -1		/* Indicates an unassigned binary file */
#define inirawf(x) open(x,O_CREAT|O_TRUNC|O_RDWR|O_BINARY)
#define opnrawf(x) open(x,O_RDWR|O_BINARY)
#define rdorawf(x) open(x,O_RDONLY|O_BINARY)
#define rdrawf(x,y,z)  read(x,y,z)	/* Read from a binary file */
#define wrrawf(x,y,z) write(x,y,z)	/* Write to a binary file */
#define clsrawf(x)  close(x)	/* Close a binary file */
#define rewrawf(x) lseek(x,0L,SEEK_SET)		/* Rewind a binary file */
#define putrawf(x,y) lseek(x,y,SEEK_SET)	/* Position a binary file */
#define posrawf(x) lseek(x,0L,SEEK_CUR)		/* Return position of a binary file */
#define sizrawf(x) lseek(x,0L,SEEK_END)		/* Size of a binary file */
#define rawopner(x) (x<0)	/* Binary open error ? */
#endif /* THCPLAT */
#ifdef TSOPLAT
#define rawfile FILE*		/* The type of raw binary file handles */
#define nulrawf    NULL		/* Indicates an unassigned binary file */
#define inirawf(x) fopen(x,"w+b")	/* Initialize and open a binary file */
#define opnrawf(x) fopen(x,"r+b")	/* Open an existing binary file */
#define rdorawf(x) fopen(x,"rb")	/* Open an existing binary file */
#define rdrawf(f,r,n) fread(r,1,n,f)	/* Read a binary file */
#define wrrawf(f,r,n) fwrite(r,1,n,f)	/* Write a binary file */
#define clsrawf(x) fclose(x)	/* Close a binary file */
#define rewrawf(x) fseek(x,0L,0)	/* Rewind a binary file */
#define putrawf(x,y) fseek(x,y,0)	/* Position a binary file */
#define posrawf(x) ftell(x)	/* Position of a binary file */
#define sizrawf(x) fseek(x,0L,2)	/* Size of a binary file */
#define rawopner(x) (x==NULL)	/* Binary open error ? */
#endif /* TSOPLAT */
#ifdef VMSPLAT
#define rawfile FILE*		/* The type of raw binary file handles */
#define nulrawf    NULL
#define inirawf(x) fopen(x,"w+b")	/* Initialize a binary file */
#define opnrawf(x) fopen(x,"r+b")	/* Open an existing binary file */
#define rdorawf(x) fopen(x,"rb")	/* Open an existing binary file */
#define rdrawf(f,r,n) fread(r,1,n,f)	/* Read a binary file */
#define wrrawf(f,r,n) fwrite(r,1,n,f)	/* Write a binary file */
#define clsrawf(x) fclose(x)	/* Close a binary file */
#define rewrawf(x) fseek(x,0L,SEEK_SET)		/* Rewind a binary file */
#define putrawf(x,y) fseek(x,y,SEEK_SET)	/* Position a binary file */
#define posrawf(x) ftell(x)	/* Position of a binary file */
#define sizrawf(x) fseek(x,0L,SEEK_END)		/* Size of a binary file */
#define rawopner(x) (x==NULL)	/* Binary open error ? */
#endif /* VMSPLAT */
#ifdef SKYPLAT
#include <system.h>
#include <fcntl.h>		/* File control options for open() */
#include <sys/stdtypes.h>	/* System file information types */
#include <sys/stat.h>		/* Defines structures for stat */
#define rawfile int		/* The type of raw binary file handles */
#define nulrawf    -1		/* Indicates an unassigned binary file */
#define inirawf(x) open(x,O_CREAT|O_TRUNC|O_RDWR,S_IWRITE|S_IREAD)
#define opnrawf(x) open(x,O_RDWR)
#define rdorawf(x) open(x,O_RDONLY)
#define rdrawf(x,y,z)  read(x,y,z)	/* Read from a binary file */
#define wrrawf(x,y,z)  write(x,y,z)	/* Write to a binary file */
#define clsrawf(x)     close(x)	/* Close a binary file */
#define rewrawf(x)     lseek(x,0L,SEEK_SET)	/* Rewind a binary file */
#define putrawf(x,y)   lseek(x,y,SEEK_SET)	/* Position a binary file */
#define posrawf(x) lseek(x,0L,SEEK_CUR)		/* Return position of a binary file */
#define sizrawf(x)     lseek(x,0L,SEEK_END)	/* Size of a binary file */
#define rawopner(x) (x<0)	/* Binary open error ? */
#endif /* SKYPLAT */
#ifdef LCCPLAT
#ifdef UBIOFUNC
#include "fcntl.h"		/* File control options for open() */
#define rawfile int		/* The type of raw binary file handles */
#define nulrawf    -1		/* Indicates an unassigned binary file */
#define inirawf(x) open(x,O_CREAT | O_TRUNC | O_RDWR | O_RAW,S_IWRITE | S_IREAD)
#define opnrawf(x) open(x,O_RDWR | O_RAW)
#define rdorawf(x) open(x,O_RDONLY | O_RAW)
#define rdrawf(x,y,z)  read(x,y,z)	/* Read from a binary file */
#define wrrawf(x,y,z) write(x,y,z)	/* Write to a binary file */
#define clsrawf(x)  close(x)	/* Close a binary file */
#define rewrawf(x) lseek(x,0L,0)	/* Rewind a binary file */
#define putrawf(x,y) lseek(x,y,0)	/* Position a binary file */
#define posrawf(x) lseek(x,0L,1)	/* Return position of a binary file */
#define sizrawf(x) lseek(x,0L,2)	/* Size of a binary file */
#define rawopner(x) (x<0)	/* Binary open error ? */
#endif
#endif /* LCCPLAT */
#endif /* UBIOFUNC */

#ifdef VARARGUS
#ifdef VARUNIXV
#include <varargs.h>
#else
#include <stdarg.h>
#endif /* VARUNIXV */
#endif /* VARARGUS */

#ifdef THRDFUNC
#ifdef MSCPLAT
#define TLS _declspec( thread )
#else
#define TLS
#endif /* MSCPLAT */
#else
#define TLS
#endif /* THRDFUNC */

#ifdef VIQASSRT
#ifdef ANSPLAT
#include <assert.h>
#define  DEBUG_ASSERT assert
#define  ALWAYS_ASSERT 0
#else
#define  DEBUG_ASSERT ((void)0)
#define  ALWAYS_ASSERT
#endif /* ANSPLAT */
#endif /* VIQASSRT */
