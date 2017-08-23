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
   ; name of module: list.c
   ;
   ; 07/23/00 changed by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : to use the TOPS database/suffix array algorithm
*/

#include "platform.h"					/* Define the platform hosting this code */
#include "TopsSeq.h"					/* Tool Organizing Properties Sequentially */
#include "mytype.h"
#include "list.h"

#ifdef FPROTOTYPE
void ListAdvanceCursor(theList *List)
#else
void ListAdvanceCursor(List)
  theList *List;
#endif
{
  List->currCursor = (List->currCursor)->next;
}

#ifdef FPROTOTYPE
int ListAtBeginning(theList *List)
#else
int ListAtBeginning(List)
  theList *List;
#endif
{
  return (List->currCursor == List->FirstNode && List->FirstNode != NULL);
}

#ifdef FPROTOTYPE
void ListBackUpCursor(theList *List)
#else
void ListBackUpCursor(List)
  theList *List;
#endif
{
  List->currCursor = (List->currCursor)->prev;
}

#ifdef FPROTOTYPE
void ListDeleteCurrItem(theList *List)
#else
void ListDeleteCurrItem(List)
  theList *List;
#endif
{
  if (List->FirstNode == List->currCursor) {
    List->FirstNode = (List->FirstNode)->next;
    List->currCursor = List->FirstNode;
    if (List->FirstNode != NULL) {
      (List->FirstNode)->prev = NULL;
    }
  }
  else {
    ((List->currCursor)->prev)->next = (List->currCursor)->next;
    if ((List->currCursor)->next != NULL) {
      ((List->currCursor)->next)->prev = (List->currCursor)->prev;
    }
    List->currCursor = (List->currCursor)->next;
  }
  (List->length)--;
}

#ifdef FPROTOTYPE
thePair *ListGetCurrItem(theList *List)
#else
thePair *ListGetCurrItem(List)
  theList *List;
#endif
{
  return (List->currCursor)->value;
}

#ifdef FPROTOTYPE
CursorType ListGetCurrCursor(theList *List)
#else
CursorType ListGetCurrCursor(List)
  theList *List;
#endif
{
  return List->currCursor;
}

#ifdef FPROTOTYPE
void ListInitCursor(theList *List)
#else
void ListInitCursor(List)
  theList *List;
#endif
{
  List->currCursor = List->FirstNode;
}

#ifdef FPROTOTYPE
void ListInsertAfter(theList *List, thePair *newPair)
#else
void ListInsertAfter(List, newPair)
  theList *List;
  thePair *newPair;
#endif
{
  auto ListNode *nodeAfter;

  nodeAfter = (List->currCursor)->next;
  ((List->currCursor)->next) = (ListNode *) (TopsSeqAccess(0));
  ((List->currCursor)->next)->value = newPair;
  ((List->currCursor)->next)->prev = List->currCursor;
  ((List->currCursor)->next)->next = nodeAfter;
  if (nodeAfter != NULL) {
    nodeAfter->prev = (List->currCursor)->next;
  }
  (List->length)++;
}

#ifdef FPROTOTYPE
void ListInsertFirst(theList *List, thePair *newPair)
#else
void ListInsertFirst(List, newPair)
  theList *List;
  thePair *newPair;
#endif
{
  auto ListNode *oldHead;

  oldHead = List->FirstNode;
  List->FirstNode = (ListNode *) (TopsSeqAccess(0));
  (List->FirstNode)->value = newPair;
  (List->FirstNode)->next = oldHead;
  if (oldHead != NULL) {
    oldHead->prev = List->FirstNode;
  }
  (List->FirstNode)->prev = NULL;
  (List->length)++;
}

#ifdef FPROTOTYPE
int ListIsOffEnd(theList *List)
#else
int ListIsOffEnd(List)
  theList *List;
#endif
{
  return (List->currCursor == NULL);
}

#ifdef FPROTOTYPE
void ListSetCursor(theList *List, CursorType start)
#else
void ListSetCursor(List, start)
  theList *List;
  CursorType start;
#endif
{
  List->currCursor = start;
}
