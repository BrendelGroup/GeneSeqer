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
   ; name of module: LcpTree.c -- Longest Common Prefix Tree, Implementation
   ;
   ; This module provides a Longest Common Prefix tree as support for the use
   ; and construction of suffix arrays.
   ;
   ; global symbols defined: None
   ;
   ; global functions defined:
   ;
   ; void   LcpCompare
   ; LONG   LcpCreate
   ; UBYTE  LcpLcp
   ; LONG   LcpMid
   ; UBYTE  LcpMinHeight
   ; UBYTE  LcpRcp
   ; void   LcpRestoreTree
   ; void   LcpSaveTree
   ; void   LcpSet
   ;
   ; useful notes and assumptions: None
   ;
   ; 07/23/00 changed by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : to use the TOPS database/suffix array algorithm
   ;
   ;/hdr/ ************************************************************************
 */

#define MEMOFUNC					/* Gives access to memory allocation */
#define OSYSFUNC					/* Gives access to the operating system */
#define PRIOFUNC					/* Gives access to PROMULA I/O functions */
#define RAWMFUNC					/* Gives access to raw memory operations */
#include "platform.h"					/* Define the platform hosting this code */
#include "TopsWire.h"					/* The operation if wiring files to memory */
#include "LcpTree.h"

UBYTE *treearray = NULL;				/*  2 x character */

static LONG Y = 0;
static LONG X2 = 0;
static LONG X2Path = 0;
static LONG X = 0;
static int Height = 0;
static LONG TreeSize = 0;
static int TreProvider = 0;

#ifdef MSCPLAT
static unsigned _int64 height[32] =
#else
static unsigned long long int height[32] =
#endif
{
   /*  0 */ 0,
   /*  1 */ 1,
   /*  2 */ 3,
   /*  3 */ 7,
   /*  4 */ 15,
   /*  5 */ 31,
   /*  6 */ 63,
   /*  7 */ 127,
   /*  8 */ 255,
   /*  9 */ 511,
   /* 10 */ 1023,
   /* 11 */ 2047,
   /* 12 */ 4095,
   /* 13 */ 8191,
   /* 14 */ 16383,
   /* 15 */ 32767,
   /* 16 */ 65535,
   /* 17 */ 131071,
   /* 18 */ 262143,
   /* 19 */ 524287,
   /* 20 */ 1048575,
   /* 21 */ 2097151,
   /* 22 */ 4194303,
   /* 23 */ 8388607,
   /* 24 */ 16777215,
   /* 25 */ 33554431,
   /* 26 */ 67108863,
   /* 27 */ 134217727,
   /* 28 */ 268435455,
   /* 29 */ 536870911,
   /* 30 */ 1073741823,
   /* 31 */ 2147483647
};

#ifdef MSCPLAT
static unsigned _int64 Pow2[] =
#else
static unsigned long long int Pow2[] =
#endif
{
   /*  0 */ 1,
   /*  1 */ 2,
   /*  2 */ 4,
   /*  3 */ 8,
   /*  4 */ 16,
   /*  5 */ 32,
   /*  6 */ 64,
   /*  7 */ 128,
   /*  8 */ 256,
   /*  9 */ 512,
   /* 10 */ 1024,
   /* 11 */ 2048,
   /* 12 */ 4096,
   /* 13 */ 8192,
   /* 14 */ 16384,
   /* 15 */ 32768,
   /* 16 */ 65536,
   /* 17 */ 131072,
   /* 18 */ 262144,
   /* 19 */ 524288,
   /* 20 */ 1048576,
   /* 21 */ 2097152,
   /* 22 */ 4194304,
   /* 23 */ 8388608,
   /* 24 */ 16777216,
   /* 25 */ 33554432,
   /* 26 */ 67108864,
   /* 27 */ 134217728,
   /* 28 */ 268435456,
   /* 29 */ 536870912,
   /* 30 */ 1073741824
   /* 31  2147483648 */
};

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine:
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
LONG LcpCreate(LONG DataSize)
#else
LONG LcpCreate(DataSize)
  LONG DataSize;
#endif
{
  auto int i;
  auto LONG temp;

  temp = DataSize;
  for (i = 0; temp >= 2; i++)
    temp /= 2;
  Height = i;
  X = DataSize - Pow2[Height];
  X2 = 2 * X;
  X2Path = X2 * 2 + 2;
  Y = height[Height] - X;
  if (X > 0)
    Height++;

  TreeSize = 2 * DataSize - 1;
  treearray = (UBYTE *) (getmem(TreeSize));
  if (treearray == NULL) {
    printf("Not enough memory to initialize the InternalTree!\n");
    exit(1);
  }
  filmem(treearray, TreeSize, 255);

  return TreeSize;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine:
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
void LcpSet(LONG ith, LONG value)
#else
void LcpSet(ith, value)
  LONG ith;
  LONG value;
#endif
{
  auto LONG current;
  auto LONG theParent;

  if (value > 255)
    value = 255;

   /*--------------------------------------------------------------------------
   ;
   ; the leaf is at the Height-1 level and the (ith-2x+x)th node
   ; while the index of the leftmost at the Height-1 level node is
   ; (2**(Height-1)-1). so, it should be 2**(Height-1)-1+(ith-x)
   ; to make it easier, set X2=2*X, and set Y=2**(Height-1)-1-x;
   ;
   ; else
   ;
   ; the leaf is at the Height level and the ith node.
   ; the corresponding index is 2**(Height)-1+ith;
   ;
   ;-------------------------------------------------------------------------*/

  if (ith >= X2)
    current = ith + Y;
  else
    current = ith + height[Height];
  treearray[current] = (UBYTE) (value);
  for (;;) {
    theParent = (current - 1) / 2;
    if (value >= treearray[theParent])
      break;
    treearray[theParent] = (UBYTE) (value);
    if (theParent == 0)
      break;
    current = theParent;
  }
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine: LcpPathNumber
   ;
   ; This function computes the pathnumber corresponding to a specified leaf
   ; number. The path number is the number one would get from numbering the
   ; nodes of the binary tree where the leftmost terminal node is one, its
   ; parent is two its next child is 3, the next higer parent is 4, and so on.
   ;
   ;                  ----------------P10---------------
   ;        --------P08---------                    ----P24----
   ;    ---P04---          ---P12---         ---P20---         ---P28---
   ;   P02      P06      P10      P14      P18      P22      P26      P30
   ; P01 P03  P05 P07  P09 P11  P13 P15  P17 P19  P21 P23  P25 P27  P29 P31
   ; L00 L01  L02 L03  L04 L05  L06 L07  L08 L09  L10 L11  L12 L13  L14 L15
   ;
   ; Thus if the number of leaves is a power of two then the pathnumber is simply
   ; the leaf number times 2 plus one. In most cases the number of leaves is not
   ; a power of two. In this case we have a slightly unbalanced tree, where the
   ; point of inbalance occurs at the value of X2 as computed below in LcpCreate.
   ; Note see also LcpLeafToIndex where X2 is used in its original form.
   ;
   ; The calulation of X2 is as follows where Pow2[Height] gives us the largest
   ; power of two less than or equal to DataSize (the number of leaves here).
   ;
   ;   X = DataSize - Pow2[Height];
   ;   X2 = 2*X;
   ;
   ; Thus assume DataSize is 11 in this case. Then Height is 3, Pow2[Height] is 8,
   ; X is 3 and X2 is 6.
   ;
   ; The binary tree the looks as follows:
   ;
   ;                  ----------------P10---------------
   ;        --------P08---------                    ----P24----
   ;    ---P04---          ---P12---         ---P20---         ---P28---
   ;   P02      P06      P10      P14      P18      P22      P26      P30
   ; P01 P03  P05 P07  P09 P11    L06      L07      L08      L09      L10
   ; L00 L01  L02 L03  L04 L05
   ;
   ; This of course is because the trick we are using is that we are using
   ; the second level node storage for our terminal nodes to pack the tree
   ; slightly.
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
static LONG LcpPathNumber(LONG Leaf)
#else
static LONG LcpPathNumber(Leaf)
  LONG Leaf;
#endif
{
  if (X2 == 0 || Leaf < X2)
    return Leaf * 2 + 1;
  return X2Path + (Leaf - X2) * 4;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine:
   ;
   ; calling parameters:
   ;
   ; return parameters:
   ;
   ; useful notes and assumptions:
   ;
   ; the leaf is at the Height-1 level and the (ith-2x+x)th node
   ; while the index of the leftmost at the Height-1 level node is
   ; (2**(Height-1)-1). so, it should be 2**(Height-1)-1+(ith-x)
   ; to make it easier, set X2=2*X, and set Y=2**(Height-1)-1-x;
   ;
   ; else
   ;
   ; the leaf is at the Height level and the ith node.
   ; the corresponding index is 2**(Height)-1+ith;
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
LONG LcpLeafToIndex(LONG ith)
#else
LONG LcpLeafToIndex(ith)
  LONG ith;
#endif
{
  if (ith >= X2)
    return ith + Y;
  else
    return ith + height[Height];
}

#ifdef FPROTOTYPE
static LONG LcpNca(LONG left_path, LONG right_path)
#else
static LONG LcpNca(left_path, right_path)
  LONG left_path;
  LONG right_path;
#endif
{
  auto ULONG theNca;
  auto int theHeight;
  auto LONG thePos;
  auto int Msb;
  auto int Lsb;

  theNca = left_path ^ right_path;

   /*--------------------------------------------------------------------------
   ;
   ; find the postion of leftmost 1-bit in theNca
   ; use ">>": if the leftmost 1-bit is being shifted to the right
   ; theNca==0, the "for" cycle is stopped. So, i keep the postion
   ; of the leftmost 1-bit (but it counts from the right)
   ;
   ;-------------------------------------------------------------------------*/

  for (Msb = 0; theNca; Msb++) {
    theNca = theNca >> 1;
  }

   /*-------------------------------------------------------------------------
   ;
   ; get the common path or common ancestor
   ; get rid the rightmost Msb bit information
   ;
   ;------------------------------------------------------------------------*/

  theNca = left_path >> (Msb - 1);

   /*-------------------------------------------------------------------------
   ;
   ; set the last bit equal to 1
   ;
   ;-------------------------------------------------------------------------*/

  theNca = theNca | 1;
  theNca = theNca << (Msb - 1);

   /*--------------------------------------------------------------------------
   ;
   ; Now, theNca contains the path number of the common ancestor
   ; then, transfer the path number to the index of tree_array
   ; first, decide the position of the node (row, col) or (height,pos)
   ;
   ; find the rightmost 1-bit
   ;
   ;-------------------------------------------------------------------------*/

  Lsb = 0;
  while ((theNca & 1) == 0) {
    theNca = theNca >> 1;
    Lsb++;
  }

  theHeight = Height - Lsb;
  thePos = (theNca - 1) >> 1;				/* keep the rightmost (Lsb-1) bit */
  return (thePos + height[theHeight]);
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine:
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
UBYTE LcpMinHeight(LONG left_leaf, LONG right_leaf)
#else
UBYTE LcpMinHeight(left_leaf, right_leaf)
  LONG left_leaf;
  LONG right_leaf;
#endif
{
  auto LONG theParent;
  auto LONG curr;
  auto LONG right_child;
  auto LONG left_child;
  auto LONG theNca;
  auto UBYTE Rlcp;
  auto LONG left_index;
  auto LONG right_index;
  auto UBYTE Llcp;

   /*-------------------------------------------------------------------------
   ;
   ; get the path number for both left and right path
   ;
   ;------------------------------------------------------------------------*/

  left_index = LcpLeafToIndex(left_leaf);
  right_index = LcpLeafToIndex(right_leaf);
  if (left_index == right_index)
    return treearray[left_index];

   /*-------------------------------------------------------------------------
   ;
   ; get the nca of the two path
   ;
   ;------------------------------------------------------------------------*/

  theNca = LcpNca(LcpPathNumber(left_leaf), LcpPathNumber(right_leaf));

   /*--------------------------------------------------------------------------
   ;
   ; move up to the nca along the left path
   ;
   ;-------------------------------------------------------------------------*/

  curr = left_index;
  Llcp = treearray[curr];
  while ((theParent = LcpParent(curr)) != theNca) {
    if ((curr & 1) == 1) {				/* this is a left child */
      if (Llcp == treearray[curr]) {
	    /*---------------------------------------------------------------
            ;
            ; need not compare with right child
            ; move to parent directly
            ;
            ;----------------------------------------------------------------*/

	curr = theParent;
	Llcp = treearray[curr];
      }
      else {
	curr = theParent;
	right_child = LcpRightchild(theParent);
	if (Llcp > treearray[right_child])
	  Llcp = treearray[right_child];
      }
    }
    else {
	 /*------------------------------------------------------------------
         ;
         ; this is a right child
         ; so, move up directly without changing Llcp
         ;
         ;-------------------------------------------------------------------*/

      curr = theParent;
    }
  }
   /*--------------------------------------------------------------------------
   ;
   ; Similar, get the Rlcp
   ;
   ;-------------------------------------------------------------------------*/

  curr = right_index;
  Rlcp = treearray[curr];
  while ((theParent = LcpParent(curr)) != theNca) {
    if ((curr & 1) == 0) {				/* this is a right child */
      if (Rlcp == treearray[curr]) {
	    /*-----------------------------------------------------------------
            ;
            ; need not compare with right child
            ; move to parent directly
            ;
            ;----------------------------------------------------------------*/

	curr = theParent;
	Rlcp = treearray[curr];
      }
      else {
	curr = theParent;
	left_child = LcpLeftchild(theParent);
	if (Rlcp > treearray[left_child])
	  Rlcp = treearray[left_child];
      }
    }
    else {
	 /*--------------------------------------------------------------------
         ;
         ; this is a left child
         ; so, move up directly without changing Rlcp
         ;
         ;------------------------------------------------------------------*/

      curr = theParent;
    }
  }
  if (Llcp <= Rlcp)
    return Llcp;
  return Rlcp;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine:
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
void LcpSaveTree(char *TreeName)
#else
void LcpSaveTree(TreeName)
  char *TreeName;
#endif
{
  auto binfile TreHandle;

  if (treearray != NULL) {
    TreHandle = inibinf(TreeName);
    if (binopner(TreHandle)) {
      printf("Unable to create Tree File %s\n", TreeName);
      exit(1);
    }
    wrbinf(TreHandle, treearray, TreeSize);
    clsbinf(TreHandle);
    free(treearray);
    treearray = NULL;
  }
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine:
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
void LcpRestoreTree(char *TreeName)
#else
void LcpRestoreTree(TreeName)
  char *TreeName;
#endif
{
  auto LONG DataSize;
  auto int i;
  auto LONG temp;

  TreProvider = TopsWireOpen(TreeName);
  if (TreProvider == 0) {
    printf("Unable to open Tree File %s\n", TreeName);
    exit(1);
  }
  TreeSize = TopsWireSize(TreProvider);
  treearray = TopsWireAddress(TreProvider);
  DataSize = (TreeSize + 1) / 2;
  temp = DataSize;
  for (i = 0; temp >= 2; i++)
    temp /= 2;
  Height = i;
  X = DataSize - Pow2[Height];
  X2 = 2 * X;
  X2Path = X2 * 2 + 2;
  Y = height[Height] - X;
  if (X > 0)
    Height++;

}

#ifdef FPROTOTYPE
void LcpDestroyTree(void)
#else
void LcpDestroyTree()
#endif
{
  if (treearray != NULL) {
    if (TreProvider != 0)
      TopsWireClose(TreProvider);
    else
      free(treearray);
    treearray = NULL;
    TreProvider = 0;
  }
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine:
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
LONG LcpMid(LONG Left, LONG Right)
#else
LONG LcpMid(Left, Right)
  LONG Left;
  LONG Right;
#endif
{
  auto LONG result;

  if ((Left < X2) && (Right >= X2)) {
    result = (Left + 1) / 2 + Right - X;
    if (result < X2)
      return result;
    else
      return (result / 2 + X);
  }
  return (Left + Right) / 2;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine:
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
UBYTE LcpLcp(LONG Position)
#else
UBYTE LcpLcp(Position)
  LONG Position;
#endif
{
  auto LONG LeftPos;

  LeftPos = LcpLeftchild(Position);
  return treearray[LeftPos];
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; name of routine:
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
UBYTE LcpRcp(LONG Position)
#else
UBYTE LcpRcp(Position)
  LONG Position;
#endif
{
  auto LONG RightPos;

  RightPos = LcpRightchild(Position);
  return treearray[RightPos];
}

#ifdef FPROTOTYPE
LONG LcpSearchRw(LONG Lw, int matchsize)
#else
LONG LcpSearchRw(Lw, matchsize)
  LONG Lw;
  int matchsize;
#endif
{
  auto LONG current;
  auto int Value;

  for (;;) {
    Lw++;
    if (Lw >= X2)
      current = Lw + Y;
    else
      current = Lw + height[Height];
    Value = treearray[current];
    if (Value < matchsize)
      return Lw - 1;
  }
}
#ifdef FPROTOTYPE
LONG LcpSearchLw(LONG Lw, int matchsize)
#else
LONG LcpSearchLw(Lw, matchsize)
  LONG Lw;
  int matchsize;
#endif
{
  auto LONG current;
  auto int Value;

  for (;;) {
    if (Lw >= X2)
      current = Lw + Y;
    else
      current = Lw + height[Height];
    Value = treearray[current];
    if (Value < matchsize)
      return Lw;
    Lw--;
  }
}
