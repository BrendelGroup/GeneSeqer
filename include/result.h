#ifndef _RESULT_H_
#define _RESULT_H_

#ifdef FPROTOTYPE
void ResultCreate (char *indexfilename);
int ResultSearch2 (char *pattern, int length, int matchsize);
void ResultDestroy (void);
void ResultFirst (void);
thePair *ResultNext (void);
#else
extern void ResultCreate ();
extern int ResultSearch2 ();
extern void ResultDestroy ();
extern void ResultFirst ();
extern thePair *ResultNext ();
#endif

#endif
