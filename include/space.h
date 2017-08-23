#ifndef _SPACE_H_
#define _SPACE_H_
#include <stdio.h>
#include <stdlib.h>

#define ALLOCSPACE(P,T,S) \
        if((P = realloc(P,(size_t)(sizeof(T)*S))) == NULL ){ \
            fprintf(stderr,"not enough memory!\n"); \
        }
            

#endif
