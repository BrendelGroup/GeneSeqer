#ifndef STRAPP_H
#define STRAPP_H

#include "space.h"
#define INC     115000;

#define CHECKALIGNLEN   { if(nextFree+500>=allocated){ \
                             allocated+=INC; \
                             ALLOCSPACE(global_algnmnt,char,allocated); \
                           }\
                         }   

static unsigned int nextFree=0; 
static unsigned int allocated=0;

extern char *global_algnmnt;
static void AppendChar(char X);



static void AppendChar(char X){
    if(nextFree+1>=allocated){
        allocated+=INC; 
        ALLOCSPACE(global_algnmnt,char,allocated); 
    } 
    global_algnmnt[nextFree++]=X; 
}

#endif
