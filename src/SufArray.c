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
   ; name of module: SufArray.c -- Suffix Arrays -- Implementation
   ;
   ; This module provides utilities for working with Suffix Arrays.
   ;
   ; By definition, given an m-character string T, a suffix array for T is an array
   ; of integers in the range 1 to m, specifying the lexicographic order of the
   ; m suffixes of string T. Consider the string
   ;
   ; "FredGoodman" of length 11. The suffixes then are
   ;
   ; Order  Suffix
   ;  1     fredgoodman
   ;  2     redgoodman
   ;  3     edgoodman
   ;  4     dgoodman
   ;  5     goodman
   ;  6     oodman
   ;  7     odman
   ;  8     dman
   ;  9     man
   ; 10     an
   ; 11     n
   ;
   ; Putting these in lexicographic order we get
   ;
   ; Order  Suffix
   ; 10     an
   ;  4     dgoodman
   ;  8     dman
   ;  3     edgoodman
   ;  1     fredgoodman
   ;  5     goodman
   ;  9     man
   ; 11     n
   ;  7     odman
   ;  6     oodman
   ;  2     redgoodman
   ;
   ; Therefore the suffix array for fredgoodman would be
   ;
   ; 10,4,8,3,1,5,9,11,7,6,2
   ;
   ; useful notes and assumptions:
   ;
   ; 07/23/00 changed by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : to use the TOPS database/suffix array algorithm
   ;
   ;/hdr/ **********************************************************************
 */

#define MEMOFUNC					/* Gives access to memory allocation */
#define OSYSFUNC					/* Gives access to the operating system */
#define PRIOFUNC					/* Gives access to PROMULA I/O functions */
#define RAWMFUNC					/* Gives access to raw memory operations */
#define STRGFUNC					/* Gives access to string manipulation */
#include "platform.h"					/* Define the platform hosting this code */
#include "EstData.h"
#include "LcpTree.h"
#include "TopsWire.h"					/* The operation if wiring files to memory */
#define SUFARRAY_DECLARE
#include "SufArray.h"

extern char *Tstring;

#define USING_BH_BITS              0
#define USING_BHBYTES              1
#define ALPHABET                   5			/* ATGC and a delimiter */
#define BH_TRUE                    1
#define B2H_TRUE                   2

static int ArraySize = 10000;
static UBYTE *BH = NULL;
static UBYTE *B2H = NULL;
static LONG *SuffixOrder = NULL;
static LONG *Count = NULL;
static int CurrentApproach = USING_BH_BITS;
static int SeparateLoops = 0;
static LONG leftmost = 0;
static LONG H;
static LONG nTrueBH = 0;
static int SufProvider = 0;

static UBYTE true8[8] =					/* Contains the masks needed to determine */
{							/* whether a given bit in a byte is True  */
				/*  0 */ 0x80,
				/* i.e, set. It is also used to set given */
				/*  1 */ 0x40,
				/* bits and thus to make them false.      */
   /*  2 */ 0x20,
   /*  3 */ 0x10,
   /*  4 */ 0x08,
   /*  5 */ 0x04,
   /*  6 */ 0x02,
   /*  7 */ 0x01
};


static UBYTE mask8[8] =
{
   /*  0 */ 0xff,
   /*  1 */ 0x7f,
   /*  2 */ 0x3f,
   /*  3 */ 0x1f,
   /*  4 */ 0x0f,
   /*  5 */ 0x07,
   /*  6 */ 0x03,
   /*  7 */ 0x01
};

static UBYTE zero8[8] =
{
   /*  0 */ 0x00,
   /*  1 */ 0x80,
   /*  2 */ 0xc0,
   /*  3 */ 0xe0,
   /*  4 */ 0xf0,
   /*  5 */ 0xf8,
   /*  6 */ 0xfc,
   /*  7 */ 0xfe
};

static int LeftMostBit[256] =
{
   /* 00 00000000 */ 0,
   /* 01 00000001 */ 1,
   /* 02 00000010 */ 2,
   /* 03 00000011 */ 2,
   /* 04 00000100 */ 3,
   /* 05 00000101 */ 3,
   /* 06 00000110 */ 3,
   /* 07 00000111 */ 3,
   /* 08 00001000 */ 4,
   /* 09 00001001 */ 4,
   /* 0a 00001010 */ 4,
   /* 0b 00001011 */ 4,
   /* 0c 00001100 */ 4,
   /* 0d 00001101 */ 4,
   /* 0e 00001110 */ 4,
   /* 0f 00001111 */ 4,
   /* 10 00010000 */ 5,
   /* 11 00010001 */ 5,
   /* 12 00010010 */ 5,
   /* 13 00010011 */ 5,
   /* 14 00010100 */ 5,
   /* 15 00010101 */ 5,
   /* 16 00010110 */ 5,
   /* 17 00010111 */ 5,
   /* 18 00011000 */ 5,
   /* 19 00011001 */ 5,
   /* 1a 00011010 */ 5,
   /* 1b 00011011 */ 5,
   /* 1c 00011100 */ 5,
   /* 1d 00011101 */ 5,
   /* 1e 00011110 */ 5,
   /* 1f 00011111 */ 5,
   /* 20 00100000 */ 6,
   /* 21 00100001 */ 6,
   /* 22 00100010 */ 6,
   /* 23 00100011 */ 6,
   /* 24 00100100 */ 6,
   /* 25 00100101 */ 6,
   /* 26 00100110 */ 6,
   /* 27 00100111 */ 6,
   /* 28 00101000 */ 6,
   /* 29 00101001 */ 6,
   /* 2a 00101010 */ 6,
   /* 2b 00101011 */ 6,
   /* 2c 00101100 */ 6,
   /* 2d 00101101 */ 6,
   /* 2e 00101110 */ 6,
   /* 2f 00101111 */ 6,
   /* 30 00110000 */ 6,
   /* 31 00110001 */ 6,
   /* 32 00110010 */ 6,
   /* 33 00110011 */ 6,
   /* 34 00110100 */ 6,
   /* 35 00110101 */ 6,
   /* 36 00110110 */ 6,
   /* 37 00110111 */ 6,
   /* 38 00111000 */ 6,
   /* 39 00111001 */ 6,
   /* 3a 00111010 */ 6,
   /* 3b 00111011 */ 6,
   /* 3c 00111100 */ 6,
   /* 3d 00111101 */ 6,
   /* 3e 00111110 */ 6,
   /* 3f 00111111 */ 6,
   /* 40 01000000 */ 7,
   /* 41 01000001 */ 7,
   /* 42 01000010 */ 7,
   /* 43 01000011 */ 7,
   /* 44 01000100 */ 7,
   /* 45 01000101 */ 7,
   /* 46 01000110 */ 7,
   /* 47 01000111 */ 7,
   /* 48 01001000 */ 7,
   /* 49 01001001 */ 7,
   /* 4a 01001010 */ 7,
   /* 4b 01001011 */ 7,
   /* 4c 01001100 */ 7,
   /* 4d 01001101 */ 7,
   /* 4e 01001110 */ 7,
   /* 4f 01001111 */ 7,
   /* 50 01010000 */ 7,
   /* 51 01010001 */ 7,
   /* 52 01010010 */ 7,
   /* 53 01010011 */ 7,
   /* 54 01010100 */ 7,
   /* 55 01010101 */ 7,
   /* 56 01010110 */ 7,
   /* 57 01010111 */ 7,
   /* 58 01011000 */ 7,
   /* 59 01011001 */ 7,
   /* 5a 01011010 */ 7,
   /* 5b 01011011 */ 7,
   /* 5c 01011100 */ 7,
   /* 5d 01011101 */ 7,
   /* 5e 01011110 */ 7,
   /* 5f 01011111 */ 7,
   /* 60 01100000 */ 7,
   /* 61 01100001 */ 7,
   /* 62 01100010 */ 7,
   /* 63 01100011 */ 7,
   /* 64 01100100 */ 7,
   /* 65 01100101 */ 7,
   /* 66 01100110 */ 7,
   /* 67 01100111 */ 7,
   /* 68 01101000 */ 7,
   /* 69 01101001 */ 7,
   /* 6a 01101010 */ 7,
   /* 6b 01101011 */ 7,
   /* 6c 01101100 */ 7,
   /* 6d 01101101 */ 7,
   /* 6e 01101110 */ 7,
   /* 6f 01101111 */ 7,
   /* 70 01110000 */ 7,
   /* 71 01110001 */ 7,
   /* 72 01110010 */ 7,
   /* 73 01110011 */ 7,
   /* 74 01110100 */ 7,
   /* 75 01110101 */ 7,
   /* 76 01110110 */ 7,
   /* 77 01110111 */ 7,
   /* 78 01111000 */ 7,
   /* 79 01111001 */ 7,
   /* 7a 01111010 */ 7,
   /* 7b 01111011 */ 7,
   /* 7c 01111100 */ 7,
   /* 7d 01111101 */ 7,
   /* 7e 01111110 */ 7,
   /* 7f 01111111 */ 7,
   /* 80 10000000 */ 8,
   /* 81 10000001 */ 8,
   /* 82 10000010 */ 8,
   /* 83 10000011 */ 8,
   /* 84 10000100 */ 8,
   /* 85 10000101 */ 8,
   /* 86 10000110 */ 8,
   /* 87 10000111 */ 8,
   /* 88 10001000 */ 8,
   /* 89 10001001 */ 8,
   /* 8a 10001010 */ 8,
   /* 8b 10001011 */ 8,
   /* 8c 10001100 */ 8,
   /* 8d 10001101 */ 8,
   /* 8e 10001110 */ 8,
   /* 8f 10001111 */ 8,
   /* 90 10010000 */ 8,
   /* 91 10010001 */ 8,
   /* 92 10010010 */ 8,
   /* 93 10010011 */ 8,
   /* 94 10010100 */ 8,
   /* 95 10010101 */ 8,
   /* 96 10010110 */ 8,
   /* 97 10010111 */ 8,
   /* 98 10011000 */ 8,
   /* 99 10011001 */ 8,
   /* 9a 10011010 */ 8,
   /* 9b 10011011 */ 8,
   /* 9c 10011100 */ 8,
   /* 9d 10011101 */ 8,
   /* 9e 10011110 */ 8,
   /* 9f 10011111 */ 8,
   /* a0 10100000 */ 8,
   /* a1 10100001 */ 8,
   /* a2 10100010 */ 8,
   /* a3 10100011 */ 8,
   /* a4 10100100 */ 8,
   /* a5 10100101 */ 8,
   /* a6 10100110 */ 8,
   /* a7 10100111 */ 8,
   /* a8 10101000 */ 8,
   /* a9 10101001 */ 8,
   /* aa 10101010 */ 8,
   /* ab 10101011 */ 8,
   /* ac 10101100 */ 8,
   /* ad 10101101 */ 8,
   /* ae 10101110 */ 8,
   /* af 10101111 */ 8,
   /* b0 10110000 */ 8,
   /* b1 10110001 */ 8,
   /* b2 10110010 */ 8,
   /* b3 10110011 */ 8,
   /* b4 10110100 */ 8,
   /* b5 10110101 */ 8,
   /* b6 10110110 */ 8,
   /* b7 10110111 */ 8,
   /* b8 10111000 */ 8,
   /* b9 10111001 */ 8,
   /* ba 10111010 */ 8,
   /* bb 10111011 */ 8,
   /* bc 10111100 */ 8,
   /* bd 10111101 */ 8,
   /* be 10111110 */ 8,
   /* bf 10111111 */ 8,
   /* c0 11000000 */ 8,
   /* c1 11000001 */ 8,
   /* c2 11000010 */ 8,
   /* c3 11000011 */ 8,
   /* c4 11000100 */ 8,
   /* c5 11000101 */ 8,
   /* c6 11000110 */ 8,
   /* c7 11000111 */ 8,
   /* c8 11001000 */ 8,
   /* c9 11001001 */ 8,
   /* ca 11001010 */ 8,
   /* cb 11001011 */ 8,
   /* cc 11001100 */ 8,
   /* cd 11001101 */ 8,
   /* ce 11001110 */ 8,
   /* cf 11001111 */ 8,
   /* d0 11010000 */ 8,
   /* d1 11010001 */ 8,
   /* d2 11010010 */ 8,
   /* d3 11010011 */ 8,
   /* d4 11010100 */ 8,
   /* d5 11010101 */ 8,
   /* d6 11010110 */ 8,
   /* d7 11010111 */ 8,
   /* d8 11011000 */ 8,
   /* d9 11011001 */ 8,
   /* da 11011010 */ 8,
   /* db 11011011 */ 8,
   /* dc 11011100 */ 8,
   /* dd 11011101 */ 8,
   /* de 11011110 */ 8,
   /* df 11011111 */ 8,
   /* e0 11100000 */ 8,
   /* e1 11100001 */ 8,
   /* e2 11100010 */ 8,
   /* e3 11100011 */ 8,
   /* e4 11100100 */ 8,
   /* e5 11100101 */ 8,
   /* e6 11100110 */ 8,
   /* e7 11100111 */ 8,
   /* e8 11101000 */ 8,
   /* e9 11101001 */ 8,
   /* ea 11101010 */ 8,
   /* eb 11101011 */ 8,
   /* ec 11101100 */ 8,
   /* ed 11101101 */ 8,
   /* ee 11101110 */ 8,
   /* ef 11101111 */ 8,
   /* f0 11110000 */ 8,
   /* f1 11110001 */ 8,
   /* f2 11110010 */ 8,
   /* f3 11110011 */ 8,
   /* f4 11110100 */ 8,
   /* f5 11110101 */ 8,
   /* f6 11110110 */ 8,
   /* f7 11110111 */ 8,
   /* f8 11111000 */ 8,
   /* f9 11111001 */ 8,
   /* fa 11111010 */ 8,
   /* fb 11111011 */ 8,
   /* fc 11111100 */ 8,
   /* fd 11111101 */ 8,
   /* fe 11111110 */ 8,
   /* ff 11111111 */ 8
};

#ifdef FPROTOTYPE
static void SetBHtrue(LONG index)
#else
static void SetBHtrue(index)
  LONG index;
#endif
{
  auto LONG iByte;
  auto int iBit;

  iByte = index / 8;
  iBit = index % 8;
  nTrueBH++;
  BH[iByte] |= true8[iBit];
}

#ifdef FPROTOTYPE
static void SetB2Htrue(LONG index)
#else
static void SetB2Htrue(index)
  LONG index;
#endif
{
  auto LONG iByte;
  auto LONG iBit;

  iByte = index / 8;
  iBit = index % 8;
  B2H[iByte] |= true8[iBit];
}

#ifdef FPROTOTYPE
static LONG FindTrueBH(LONG index)
#else
static LONG FindTrueBH(index)
  LONG index;
#endif
{
  auto LONG iByte;
  auto LONG iBit;
  auto UBYTE Bits;

  iByte = index / 8;
  iBit = index % 8;
  Bits = (UBYTE) (BH[iByte] & mask8[iBit]);
  while (!Bits)
    Bits = BH[++iByte];
  return iByte * 8 + (8 - LeftMostBit[Bits]);
}

#ifdef FPROTOTYPE
static int TestB2H(LONG index)
#else
static int TestB2H(index)
  LONG index;
#endif
{
  return B2H[index / 8] & true8[index % 8];
}

#ifdef FPROTOTYPE
static void ZeroBHorNotB2H(LONG index)
#else
static void ZeroBHorNotB2H(index)
  LONG index;
#endif
{
  auto LONG iByte;
  auto LONG iBit;
  auto UBYTE Bits;
  auto int lBit;

  iByte = index / 8;
  iBit = index % 8;
  Bits = (UBYTE) ((BH[iByte] | ~B2H[iByte]) & mask8[iBit]);
  if (!Bits) {
    B2H[iByte] &= zero8[iBit];
    for (;;) {
      iByte++;
      Bits = (UBYTE) (BH[iByte] | ~B2H[iByte]);
      if (Bits)
	break;
      B2H[iByte] = 0;
    }
    iBit = 8 - LeftMostBit[Bits];
    if (iBit)
      B2H[iByte] &= mask8[iBit];
  }
  else {
    lBit = 8 - LeftMostBit[Bits];
    lBit = 8 - LeftMostBit[Bits];
    if (iBit < lBit)
      B2H[iByte] &= (zero8[iBit] | mask8[lBit]);
  }
}

#ifdef FPROTOTYPE
static LONG FindNotBHandB2H(LONG index)
#else
static LONG FindNotBHandB2H(index)
  LONG index;
#endif
{
  auto LONG iByte;
  auto LONG iBit;
  auto UBYTE Bits;

  iByte = index / 8;
  iBit = index % 8;
  Bits = (UBYTE) ((~BH[iByte] & B2H[iByte]) & mask8[iBit]);
  while (!Bits) {
    iByte++;
    Bits = (UBYTE) (~BH[iByte] & B2H[iByte]);
  }
  return iByte * 8 + (8 - LeftMostBit[Bits]);
}

#ifdef FPROTOTYPE
static void InitializeBHandB2H(LONG DataSize)
#else
static void InitializeBHandB2H(DataSize)
  LONG DataSize;
#endif
{
  auto LONG BitBytes;

   /*-------------------------------------------------------------------------
   ;
   ; Step 1: Determine the number of bytes needed to contain DataSize+2 bits.
   ; We need two extra bits to supply stop conditions for FindTrueBH and for
   ; FindNotBHandB2H.
   ; 
   ;------------------------------------------------------------------------*/

  BitBytes = (DataSize + 1) / 8 + 1;

   /*--------------------------------------------------------------------------
   ;
   ; Step 2: Allocate the two memory memory areas. If there is a problem,
   ; simply exit to the operationg system with a message.
   ;
   ;-------------------------------------------------------------------------*/

  BH = (UBYTE *) (getmem(BitBytes));
  if (BH == NULL) {
    printf("Unable to allocate BH bits array of %d bytes\n", BitBytes);
    exit(1);
  }
  B2H = (UBYTE *) (getmem(BitBytes));
  if (B2H == NULL) {
    printf("Unable to allocate B2H bits array of %d bytes\n", BitBytes);
    exit(1);
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 3: Initialize the bit vectors at all 0 (false) and set to termination
   ; conditions for the find function:
   ;
   ; FindBHtrue         BH[DataSize] = true;
   ; FindNotBHandB2H    BH[DataSize+1] = false; B2H[DataSize+1] = true.
   ; 
   ; Also initialize the value of nTrueBH at zero, to indicate that no Buckets
   ; have yet been started.
   ;
   ;-------------------------------------------------------------------------*/

  filmem(BH, BitBytes, 0);
  filmem(B2H, BitBytes, 0);
  SetBHtrue(DataSize);
  SetB2Htrue(DataSize + 1);
  nTrueBH = 0;
}

#ifdef FPROTOTYPE
static void CloseBHandB2H(void)
#else
static void CloseBHandB2H()
#endif
{
  if (BH != NULL)
    free(BH);
  BH = NULL;
  if (B2H != NULL)
    free(B2H);
  B2H = NULL;
}

#ifdef FPROTOTYPE
static int SuffixCompare(LONG suffix, char *str, LONG p)
#else
static int SuffixCompare(suffix, str, p)
  LONG suffix;
  char *str;
  LONG p;
#endif
{
  auto int i;
  auto int maxi;
  auto char *SufString;

  SufString = Tstring + suffix;
  maxi = (int) (ArraySize - suffix);

  for (i = 0; (i < p) && (i < maxi) && (SufString[i] == str[i]); i++);
  return i;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: InitializePositions
   ;
   ; We have initialized the storage areas that we need to build the suffix array
   ; and its associated longest common prefix tree. We are now going to start the
   ; sorting of the suffixes by building a linked list of the character values as
   ; they appear in the string. When we complete the linked list, we can
   ; initialize the four main vectors controlling the construction of the suffix
   ; array: BH, SuffixArray, LcpTree, and SuffixOrder. It is important to note
   ; that the characters are linked from the back of the string to the front. 
   ;
   ; calling parameters: None
   ;
   ; return parameters: None
   ;
   ; useful notes and assumptions:
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
static void InitializePositions(void)
#else
static void InitializePositions()
#endif
{
  auto LONG *PrevPosition;
  auto LONG LastPosition[ALPHABET];
  auto LONG iPos;
  auto LONG iLex;
  auto LONG lPos;
  auto int iCharValue;

   /*--------------------------------------------------------------------------
   ;
   ; Step 1: Initialize memory for storing the suffix array, the lexical
   ; order information, and the linked list. If there is a problem, print a
   ; message and exit to the operating system.
   ;
   ;-------------------------------------------------------------------------*/

  SuffixArray = (LONG *) (getmem(ArraySize * sizeof(LONG)));
  if (SuffixArray == NULL) {
    printf("Not enough memory for Suffix array!\n");
    exit(1);
  }
  SuffixOrder = (LONG *) (getmem(ArraySize * sizeof(LONG)));
  if (SuffixOrder == NULL) {
    printf("Not enough memory for Suffix Order array!\n");
    exit(1);
  }
  PrevPosition = (LONG *) (getmem(ArraySize * sizeof(LONG)));
  if (PrevPosition == NULL) {
    printf("Not enough memory for Link List array!\n");
    exit(1);
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 2: Start the sorting of the suffixes by building a linked list
   ; of the character values as they appear in the string. When we complete
   ; this step LastPosition[iCharValue] will contain the position in the string
   ; of the last occurrence of a character of that value. PrevPosition[iPos]
   ; will contain the position in the string of the previous occurrence of that
   ; character. A -1 is used to mark the end of the list.
   ;
   ;-------------------------------------------------------------------------*/

  for (iCharValue = 0; iCharValue < ALPHABET; iCharValue++) {
    LastPosition[iCharValue] = -1;
  }
  for (iPos = 0; iPos < ArraySize; iPos++) {
    iCharValue = Tstring[iPos];
    lPos = LastPosition[iCharValue];
    LastPosition[iCharValue] = iPos;
    PrevPosition[iPos] = lPos;
  }

   /*-------------------------------------------------------------------------
   ;
   ; Step 3: We no longer need the original Tstring information so we can
   ; release the memory allocated to it.
   ;
   ;-------------------------------------------------------------------------*/

  EstDestroy();
  Tstring = NULL;

   /*--------------------------------------------------------------------------
   ;
   ; Step 4: Initialize the bucket control vectors. This depends upon the
   ; initial approach to be taken to retain this information.
   ;
   ;-------------------------------------------------------------------------*/

  if (CurrentApproach == USING_BH_BITS)
    InitializeBHandB2H(ArraySize);
  else {
    BH = (UBYTE *) (getmem(ArraySize + 2));
    if (BH == NULL) {
      printf("Not enouph memory for Bucket control vector\n");
      exit(1);
    }
    filmem(BH, ArraySize, 0);
    BH[ArraySize] = BH_TRUE;
    BH[ArraySize + 1] = B2H_TRUE;
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 5: We are now ready initialize the four main vectors controlling
   ; the construction of the suffix array: BH, SuffixArray, LcpTree,
   ; and SuffixOrder. It is important to note that the characters are linked
   ; from the back of the string to the front. 
   ;
   ;-------------------------------------------------------------------------*/

  iLex = 0;
  for (iCharValue = 0; iCharValue < ALPHABET; iCharValue++) {
    if ((iPos = LastPosition[iCharValue]) != -1) {
      if (CurrentApproach == USING_BH_BITS)
	SetBHtrue(iLex);
      else {
	BH[iLex] = BH_TRUE;
	nTrueBH++;
      }
      if (iCharValue != 0)
	LcpSet(iLex, 0);
      do {
	SuffixOrder[iPos] = iLex;
	SuffixArray[iLex] = iPos;
	iLex++;
	iPos = PrevPosition[iPos];
      }
      while (iPos != -1);
    }
  }

   /*-------------------------------------------------------------------------
   ;
   ; Step 5: Free the linked list storage area.
   ;
   ;-------------------------------------------------------------------------*/

  free(PrevPosition);
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: BeginNextSortPass
   ;
   ; At this point we are ready to begin the next pass through the information.
   ; We must first set the SuffixOrder vector so that it
   ; points to the leftmost lexical order position for each set of
   ; unresolved suffixs. Note the term "bucket" is used by the original
   ; author to refer to these sets. The original code was as follows:
   ;
   ;  for(iLex = 0; i < ArraySize; iLex++) {
   ;     if(BH[iLex]) {
   ;        Count[iLex] = 0;
   ;        leftmost = iLex;
   ;     }
   ;     SuffixOrder[SuffixArray[iLex]] = leftmost;
   ;  }
   ;  upperH = ArraySize - H;
   ;  order_upperH = SuffixOrder[upperH];
   ;  SuffixOrder[upperH] = order_upperH + Count[order_upperH];
   ;  Count[order_upperH]++;
   ;  B2H[SuffixOrder[upperH]]=true;
   ;
   ; The following changes are made. First the BH boolean vector has been
   ; changed to a bit vector so that we can more rapidly find true bits.
   ; Note that during the first few passes this bit vector is very sparse.
   ;
   ; Next we notice that leftmost is changed many times in the inner loop;
   ; we can simply set it equal to the position of the leftmost member
   ; of the bucket as returned by FindBHTrue.
   ;
   ; Finally the zeroing of the Count matrix can be done outside the
   ; loop using a system fill function.
   ; 
   ; calling parameters:
   ;
   ; return parameters:
   ;
   ; useful notes and assumptions:
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
static void BeginNextSortPass(void)
#else
static void BeginNextSortPass()
#endif
{
  auto LONG iStart;
  auto LONG lStart;
  auto LONG iLex;
  auto LONG upperH;
  auto LONG order_upperH;

  switch (CurrentApproach) {
  case USING_BH_BITS:

    lStart = -1;
    for (;;) {
      iStart = lStart + 1;
      lStart = FindTrueBH(iStart);
      for (iLex = iStart; iLex < lStart; iLex++) {
	SuffixOrder[SuffixArray[iLex]] = leftmost;
      }
      if (lStart >= ArraySize)
	break;
      leftmost = lStart;
    }
    filmem(Count, ArraySize * sizeof(LONG), 0);
    upperH = ArraySize - H;
    order_upperH = SuffixOrder[upperH];
    order_upperH = order_upperH + Count[order_upperH]++;
    SuffixOrder[upperH] = order_upperH;
    SetB2Htrue(order_upperH);
    return;

  case USING_BHBYTES:

    lStart = -1;
    for (;;) {
      iStart = lStart + 1;
      for (lStart = iStart; !(BH[lStart] & BH_TRUE); lStart++);
      for (iLex = iStart; iLex < lStart; iLex++) {
	SuffixOrder[SuffixArray[iLex]] = leftmost;
      }
      if (lStart >= ArraySize)
	break;
      leftmost = lStart;
    }
    filmem(Count, ArraySize * sizeof(LONG), 0);
    upperH = ArraySize - H;
    order_upperH = SuffixOrder[upperH];
    order_upperH = order_upperH + Count[order_upperH]++;
    SuffixOrder[upperH] = order_upperH;
    BH[order_upperH] |= B2H_TRUE;
    return;
  }
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: UpdateSortBuckets
   ;
   ; We must now update the Suffix Order and B2H information through the entire
   ; set of sort buckets. In so far as trashing is concerned this block of code
   ; causes the must problems. The revision attempted to separate the references
   ; to the various arrays as much as possible. This was only partially effective
   ; (see comments below). The original approach was as follows:
   ;
   ;  for(i = 0; i < ArraySize; i++){
   ;     leftmost = i;
   ;     do{
   ;         upperH = SuffixArray[i] - H;
   ;         if(upperH >= 0){
   ;            order_upperH=Prm[upperH];
   ;            Prm[upperH]=order_upperH+Count[order_upperH];
   ;            Count[order_upperH]++;
   ;            B2H[Prm[upperH]]=true;
   ;         }
   ;         i++;
   ;     }while(BH[i]==false&&i<ArraySize);
   ;     i = leftmost;
   ;     do {
   ;        upperH = SuffixArray[i]-H;
   ;        if(upperH >=0 && B2H[Prm[upperH]])
   ;        {
   ;           j=Prm[upperH]+1;
   ;           while((!BH[j])&&B2H[j]){j++;}
   ;           rightmost2H=j;
   ;           for(j=Prm[upperH]+1;j<rightmost2H;j++)
   ;           B2H[j]=false;
   ;        }
   ;        i++;
   ;     }while(BH[i]==false&&i<ArraySize);
   ;     i--;
   ;  }
   ;
   ; The most important observation here is that the buckets are independent of
   ; themselves. The do's over "BH[i] == false && i < ArraySize" are in truth
   ; simply partitioning the operations into the same buckets that were
   ; initialized in BeginNextSortPass above and we can use the same FindBHtrue
   ; logic that we used there; thus greatly simplifying the structure of the code
   ;
   ; calling parameters: None
   ;
   ; return parameters: None
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */

#ifdef FPROTOTYPE
static void UpdateSortBuckets(void)
#else
static void UpdateSortBuckets()
#endif
{
  auto LONG rightmost;
  auto LONG iLex;
  auto LONG upperH;
  auto LONG order_upperH;
  auto LONG Temp;

  if (SeparateLoops)
    goto UseSeparateLoops;

  switch (CurrentApproach) {
  case USING_BH_BITS:
    rightmost = 0;
    for (;;) {
      leftmost = rightmost;
      rightmost = FindTrueBH(leftmost + 1);
      for (iLex = leftmost; iLex < rightmost; iLex++) {
	upperH = SuffixArray[iLex] - H;
	if (upperH >= 0) {
	  order_upperH = SuffixOrder[upperH];
	  SuffixOrder[upperH] = order_upperH + Count[order_upperH]++;
	  SetB2Htrue(SuffixOrder[upperH]);
	}
      }
      for (iLex = leftmost; iLex < rightmost; iLex++) {
	if ((upperH = SuffixArray[iLex] - H) >= 0) {
	  order_upperH = SuffixOrder[upperH];
	  if (TestB2H(order_upperH))
	    ZeroBHorNotB2H(order_upperH + 1);
	}
      }
      if (rightmost >= ArraySize)
	break;
    }
    return;

  case USING_BHBYTES:
    rightmost = 0;
    for (;;) {
      leftmost = rightmost;
      for (rightmost = leftmost + 1; !(BH[rightmost] & BH_TRUE); rightmost++);
      for (iLex = leftmost; iLex < rightmost; iLex++) {
	upperH = SuffixArray[iLex] - H;
	if (upperH >= 0) {
	  order_upperH = SuffixOrder[upperH];
	  SuffixOrder[upperH] = order_upperH + Count[order_upperH]++;
	  BH[SuffixOrder[upperH]] |= B2H_TRUE;
	}
      }
      for (iLex = leftmost; iLex < rightmost; iLex++) {
	if ((upperH = SuffixArray[iLex] - H) >= 0) {
	  order_upperH = SuffixOrder[upperH];
	  if (BH[order_upperH] & B2H_TRUE) {
	    order_upperH++;
	    while (BH[order_upperH] == B2H_TRUE)
	      BH[order_upperH++] = 0;
	  }
	}
      }
      if (rightmost >= ArraySize)
	break;
    }
    return;
  }

UseSeparateLoops:

  rightmost = 0;
  for (;;) {
    leftmost = rightmost;
    rightmost = FindTrueBH(leftmost + 1);
    for (iLex = leftmost; iLex < rightmost; iLex++) {
      upperH = SuffixArray[iLex] - H;
      if (upperH >= 0) {
	order_upperH = SuffixOrder[upperH];
	Temp = Count[order_upperH]++;
	order_upperH += Temp;
	SuffixArray[iLex] = order_upperH;
	SuffixOrder[upperH] = order_upperH;
      }
      else
	SuffixArray[iLex] = -1;
    }
    for (iLex = leftmost; iLex < rightmost; iLex++) {
      if ((order_upperH = SuffixArray[iLex]) >= 0) {
	SetB2Htrue(order_upperH);
      }
    }
    for (iLex = leftmost; iLex < rightmost; iLex++) {
      if ((order_upperH = SuffixArray[iLex]) >= 0) {
	if (TestB2H(order_upperH))
	  ZeroBHorNotB2H(order_upperH + 1);
      }
    }
    if (rightmost >= ArraySize)
      break;
  }
  return;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: RefreshSuffixArray
   ;
   ; There are two basic classification schemes bein manipulated by this program
   ; both of which have ArraySize entries:
   ;
   ;  (1)  LexicographicOrder (Note LexicalOrder seems like a better term)
   ;  (2)  PositionInString
   ;
   ; The vector SuffixArray classified by LexicalOrder contains PositionInString;
   ; while SuffixOrder classified by PositionInString contains LexicalOrder. Thus,
   ; SuffixArray[iLex] contains the position in the string of the indicated
   ; suffix; and SuffixOrder[iPos] contains the lexical order of the indicated
   ; suffix. In UpdateSortBuckets above the values of SuffixOrder are updated.
   ; These updates must now be reflected in SuffixArray.
   ;
   ; calling parameters: None
   ;
   ; return parameters: None
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
static void RefreshSuffixArray(void)
#else
static void RefreshSuffixArray()
#endif
{
  auto LONG iPos;

  for (iPos = 0; iPos < ArraySize; iPos++) {
    SuffixArray[SuffixOrder[iPos]] = iPos;
  }
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: UpdateLongestCommonPrefix
   ;
   ; calling parameters: None
   ;
   ; return parameters: None
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
static void UpdateLongestCommonPrefix(void)
#else
static void UpdateLongestCommonPrefix()
#endif
{
  auto LONG iLex;
  auto LONG aPos;
  auto LONG bPos;
  auto LONG aOrder;
  auto LONG bOrder;
  auto LONG value;

  switch (CurrentApproach) {
  case USING_BH_BITS:

    for (iLex = FindNotBHandB2H(1); iLex < ArraySize; iLex = FindNotBHandB2H(iLex + 1)) {
      SetBHtrue(iLex);
      aPos = SuffixArray[iLex - 1] + H;
      bPos = SuffixArray[iLex] + H;
      if ((aPos < ArraySize) && (bPos < ArraySize)) {
	aOrder = SuffixOrder[aPos];
	bOrder = SuffixOrder[bPos];
	if (aOrder > bOrder)
	  value = H + LcpMinHeight(bOrder + 1, aOrder);
	else
	  value = H + LcpMinHeight(aOrder + 1, bOrder);
      }
      else
	value = H;
      LcpSet(iLex, value);
    }
    return;

  case USING_BHBYTES:

    iLex = 0;
    for (;;) {
      iLex++;
      while (BH[iLex] != B2H_TRUE)
	iLex++;
      if (iLex >= ArraySize)
	break;
      BH[iLex] |= BH_TRUE;
      nTrueBH++;
      aPos = SuffixArray[iLex - 1] + H;
      bPos = SuffixArray[iLex] + H;
      if ((aPos < ArraySize) && (bPos < ArraySize)) {
	aOrder = SuffixOrder[aPos];
	bOrder = SuffixOrder[bPos];
	if (aOrder > bOrder)
	  value = H + LcpMinHeight(bOrder + 1, aOrder);
	else
	  value = H + LcpMinHeight(aOrder + 1, bOrder);
      }
      else
	value = H;
      LcpSet(iLex, value);
    }
    return;
  }
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: SuffixArrayCreate
   ;
   ; This function is used to initialize the information to be used to work with
   ; the suffix array. 
   ;
   ; calling parameters:
   ;
   ; LONG   Size       The length of the Tstring that the suffix array is to
   ;                   support.
   ;
   ; return parameters: None
   ;
   ; useful notes and assumptions:
   ;
   ; This function does nothing but save the information passed in.
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
void SuffixArrayCreate(LONG size)
#else
void SuffixArrayCreate(size)
  LONG size;
#endif
{
  ArraySize = size;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: SuffixArrayDestroy
   ;
   ; This function is used to release all resources associated with the currently
   ; active suffix array.
   ;
   ; calling parameters: None
   ;
   ; return parameters: None
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
void SuffixArrayDestroy(void)
#else
void SuffixArrayDestroy()
#endif
{
  if (SuffixArray != NULL) {
    if (SufProvider != 0)
      TopsWireClose(SufProvider);
    else
      free(SuffixArray);
    SuffixArray = NULL;
    SufProvider = 0;
  }
  LcpDestroyTree();
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: SuffixArraySave
   ;
   ; This functions saves the actual suffix array itself along with the longest
   ; common prefix tree used to support searching the suffix array.
   ;
   ; calling parameters:
   ;
   ; char*  ArrayName   The name of the file to contain the suffix array.
   ;
   ; char*  LcpName     The name of the file to contain the LcpTree.
   ;
   ; return parameters: None
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
void SuffixArraySave(char *ArrayName, char *LcpName)
#else
void SuffixArraySave(ArrayName, LcpName)
  char *ArrayName;
  char *LcpName;
#endif
{
  auto binfile SufHandle;

  SufHandle = inibinf(ArrayName);
  if (binopner(SufHandle)) {
    printf("Unable to create Suffix Array File %s\n", ArrayName);
    exit(5);
  }
  wrbinf(SufHandle, SuffixArray, ArraySize * sizeof(LONG));
  clsbinf(SufHandle);
  LcpSaveTree(LcpName);
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: SuffixArrayRestore
   ;
   ; This functions restores a previously saved suffix array along with the
   ; longest common prefix tree used to support searching the suffix array.
   ;
   ; calling parameters:
   ;
   ; char*  ArrayName   The name of the file containing the suffix array.
   ;
   ; char*  LcpName     The name of the file containing the LcpTree.
   ;
   ; return parameters: None
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
void SuffixArrayRestore(char *ArrayName, char *LcpName, LONG DataSize)
#else
void SuffixArrayRestore(ArrayName, LcpName, DataSize)
  char *ArrayName;
  char *LcpName;
  LONG DataSize;
#endif
{
  LcpRestoreTree(LcpName);
  ArraySize = DataSize;
  SufProvider = TopsWireOpen(ArrayName);
  if (SufProvider == 0) {
    printf("Unable to open Suffix Array File %s\n", ArrayName);
    exit(6);
  }
  SuffixArray = (LONG *) (TopsWireAddress(SufProvider));
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: SuffixArraySort
   ;
   ; This function performs the actual sort of the Tstring into a suffix array.
   ;
   ; This functions restores a previously saved suffix array along with the
   ; longest common prefix tree used to support searching the suffix array.
   ;
   ; calling parameters: None
   ;
   ; return parameters: None
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
void SuffixArraySort(void)
#else
void SuffixArraySort()
#endif
{
  auto int Level;

   /*-------------------------------------------------------------------------
   ;
   ; Step 1: Create the basic structures needed to build the longest common
   ; prefix tree and then initialize the storage areas need to build the
   ; suffix array. 
   ;
   ;------------------------------------------------------------------------*/

  LcpCreate(ArraySize);
  InitializePositions();
  Count = (LONG *) (getmem(ArraySize * sizeof(LONG)));
  if (Count == NULL) {
    printf("Not enough memory for Count array!\n");
    exit(1);
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 2: Sort the Suffix array. Note that there are two possible
   ; convergence criteria:
   ;
   ;  (1) When we make completed examining every possible location
   ;        H >= ArraySize
   ;
   ;  (2) When every position in the lexical order has been updated.
   ;        nTrue >= ArraySize.
   ;
   ; In the original formulation only the first criteria was checked. In
   ; examining the behavior of the code, the second criteria appears to be
   ; reached much earlier in all cases.
   ;
   ;-------------------------------------------------------------------------*/

  leftmost = 0;
  Level = 0;
  for (H = 1; nTrueBH < ArraySize && H < ArraySize; H *= 2) {
    Level++;
    BeginNextSortPass();
    UpdateSortBuckets();
    RefreshSuffixArray();
    UpdateLongestCommonPrefix();
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 3: Return the resources that were needed while to support the sorting
   ; of the suffix array.
   ;
   ;-------------------------------------------------------------------------*/

  if (CurrentApproach == USING_BH_BITS)
    CloseBHandB2H();
  else {
    free(BH);
    BH = NULL;
  }
  free(SuffixOrder);
  SuffixOrder = NULL;
  free(Count);
  Count = NULL;
}

#ifdef FPROTOTYPE
LONG *SuffixArrayVector(LONG index, LONG *Length)
#else
LONG *SuffixArrayVector(index, Length)
  LONG index;
  LONG *Length;
#endif
{
  *Length = ArraySize - index;
  return SuffixArray + index;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: SuffixArraySearchLw
   ;
   ; calling parameters:
   ;
   ; char*  pattern     The pattern to be searched for
   ;
   ; LONG   p           The length of the pattern
   ;
   ; return parameters: None
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
LONG SuffixArraySearchLw(char *pattern, LONG p, int *HaveMatch)
#else
LONG SuffixArraySearchLw(pattern, p, HaveMatch)
  char *pattern;
  LONG p;
  int *HaveMatch;
#endif
{
  auto LONG L;
  auto LONG R;
  auto LONG M;
  auto UBYTE m;
  auto LONG curr;
  auto UBYTE Rlcp;
  auto UBYTE Llcp;
  auto LONG l;
  auto LONG r;
  auto LONG Lpos;
  auto LONG Rpos;
  auto LONG Mpos;

  *HaveMatch = 0;
  Lpos = SuffixArray[0];
  l = SuffixCompare(Lpos, pattern, p);
  if ((l >= p) || (pattern[l] <= Tstring[Lpos + l])) {
    if (l >= p)
      *HaveMatch = 1;
    return 0;
  }
  Rpos = SuffixArray[ArraySize - 1];
  r = SuffixCompare(Rpos, pattern, p);
  if ((r < p) && (pattern[r] > Tstring[Rpos + r])) {
    return ArraySize;
  }

  curr = 0;
  L = 0;
  R = ArraySize - 1;
  while ((R - L) > 1) {
    M = LcpMid(L, R);
    Mpos = SuffixArray[M];
    if (l >= r) {
      if ((Llcp = LcpLcp(curr)) >= l)
	m = (UBYTE) (l + SuffixCompare(Mpos + l, pattern + l, p));
      else
	m = Llcp;
    }
    else {
      if ((Rlcp = LcpRcp(curr)) >= r)
	m = (UBYTE) (r + SuffixCompare(Mpos + r, pattern + r, p));
      else
	m = Rlcp;
    }
    if ((m >= p) || (pattern[m] <= Tstring[Mpos + m])) {
      /*
         choose the left half
       */
      if (m >= p)
	*HaveMatch = 1;
      R = M;
      r = m;
      curr = LcpLeftchild(curr);
    }
    else {
      /*
         choose the right half
       */
      L = M;
      l = m;
      curr = LcpRightchild(curr);
    }
  }
  return R;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: SuffixArraySearchRw
   ;
   ; calling parameters:
   ;
   ; char*  pattern     The pattern to be searched for
   ;
   ; LONG   p           The length of the pattern
   ;
   ; return parameters: None
   ;
   ; useful notes and assumptions: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
LONG SuffixArraySearchRw(char *pattern, LONG p)
#else
LONG SuffixArraySearchRw(pattern, p)
  char *pattern;
  LONG p;
#endif
{
  auto LONG L;
  auto LONG R;
  auto LONG M;
  auto LONG curr;
  auto UBYTE Rlcp;
  auto UBYTE Llcp;
  auto UBYTE m;
  auto LONG l;
  auto LONG r;
  auto LONG Lpos;
  auto LONG Rpos;
  auto LONG Mpos;

  Lpos = SuffixArray[0];
  Rpos = SuffixArray[ArraySize - 1];
  l = SuffixCompare(Lpos, pattern, p);
  r = SuffixCompare(Rpos, pattern, p);
  if ((r >= p) || (pattern[r] > Tstring[Rpos + r])) {
    return ArraySize - 1;
  }
  else if ((l < p) && (pattern[l] < Tstring[Lpos + l])) {
    return -1;
  }
  curr = 0;
  L = 0;
  R = ArraySize - 1;
  while ((R - L) > 1) {
    M = LcpMid(L, R);
    Mpos = SuffixArray[M];
    if (l >= r) {
      if ((Llcp = LcpLcp(curr)) >= l) {
	m = (UBYTE) (l + SuffixCompare(Mpos + l, pattern + l, p));
      }
      else
	m = Llcp;
    }
    else {
      if ((Rlcp = LcpRcp(curr)) >= r) {
	m = (UBYTE) (r + SuffixCompare(Mpos + r, pattern + r, p));
      }
      else
	m = Rlcp;
    }
    if ((m >= p) || (pattern[m] >= Tstring[Mpos + m])) {
      /*
         choose the right half
       */
      L = M;
      l = m;
      curr = LcpRightchild(curr);
    }
    else {
      /*
         bigger or equal, choose the left half
       */
      R = M;
      r = m;
      curr = LcpLeftchild(curr);
    }
  }
  return L;
}

#ifdef FPROTOTYPE
int SuffixArraySearchLwRw(char *pattern, LONG p, LONG *LwRw)
#else
int SuffixArraySearchLwRw(pattern, p, LwRw)
  char *pattern;
  LONG p;
  LONG *LwRw;
#endif
{
  auto LONG L;
  auto LONG R;
  auto LONG M;
  auto UBYTE m;
  auto LONG curr;
  auto UBYTE Rlcp;
  auto UBYTE Llcp;
  auto LONG l;
  auto LONG r;
  auto LONG Lpos;
  auto LONG Rpos;
  auto LONG Mpos;

  Tstring[ArraySize]= -1;	/* ... set the end of the text to the lexicographically least character */
  Lpos = SuffixArray[0];
  l = SuffixCompare(Lpos, pattern, p);
  if (l == p) {
    *LwRw = 0;
    *(LwRw + 1) = LcpSearchRw(0,p);
    return 1;
  }
  if (pattern[l] <= Tstring[Lpos + l])
    goto NoMatch;
  Rpos = SuffixArray[ArraySize - 1];
  r = SuffixCompare(Rpos, pattern, p);
  if (r == p) {
    *LwRw = LcpSearchLw(ArraySize-1,p);
    *(LwRw + 1) = ArraySize-1;
    return 1;
  }
  if (pattern[r] > Tstring[Rpos + r])
    goto NoMatch;

  curr = 0;
  L = 0;
  R = ArraySize - 1;
  while ((R - L) > 1) {
    M = LcpMid(L, R);
    Mpos = SuffixArray[M];
    if (l >= r) {
      if ((Llcp = LcpLcp(curr)) >= l)
	m = (UBYTE) (l + SuffixCompare(Mpos + l, pattern + l, p - l));
      else
	m = Llcp;
    }
    else {
      if ((Rlcp = LcpRcp(curr)) >= r)
	m = (UBYTE) (r + SuffixCompare(Mpos + r, pattern + r, p - r));
      else
	m = Rlcp;
    }
    if (m >= p)
      goto HaveMatch;
    if (pattern[m] <= Tstring[Mpos + m]) {
      /*
         choose the left half
       */
      if (R - L == 2  &&  l >= p) {
	M = L;
	goto HaveMatch;
      }
      R = M;
      r = m;
      curr = LcpLeftchild(curr);
    }
    else {
      /*
         choose the right half
       */
      if (R - L == 2  &&  r >= p) {
	M = R;
	goto HaveMatch;
      }
      L = M;
      l = m;
      curr = LcpRightchild(curr);
    }
  }
NoMatch:
  *LwRw = *(LwRw + 1) = -1;
  return 0;
HaveMatch:
  *LwRw = LcpSearchLw(M, p);
  *(LwRw + 1) = LcpSearchRw(M, p);
  return 1;
}
