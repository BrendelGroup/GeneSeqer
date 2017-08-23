#define ListNode struct sListNode

struct sListNode
  {
    thePair *value;
    ListNode *next;
    ListNode *prev;
  };

typedef void *CursorType;

typedef struct
  {
    ListNode *FirstNode;
    ListNode *currCursor;
    int length;
  }
theList;

#ifdef FPROTOTYPE
void ListAdvanceCursor (theList * List);
int ListAtBeginning (theList * List);
void ListBackUpCursor (theList * List);
void ListDeleteCurrItem (theList * List);
thePair *ListGetCurrItem (theList * List);
CursorType ListGetCurrCursor (theList * List);
void ListInitCursor (theList * List);
void ListInsertAfter (theList * List, thePair * newPair);
void ListInsertFirst (theList * List, thePair * newPair);
int ListIsOffEnd (theList * List);
void ListSetCursor (theList * List, CursorType start);
#else
extern void ListAdvanceCursor ();
extern int ListAtBeginning ();
extern void ListBackUpCursor ();
extern void ListDeleteCurrItem ();
extern thePair *ListGetCurrItem ();
extern CursorType ListGetCurrCursor ();
extern void ListInitCursor ();
extern void ListInsertAfter ();
extern void ListInsertFirst ();
extern int ListIsOffEnd ();
extern void ListSetCursor ();
#endif
