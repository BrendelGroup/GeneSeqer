#ifndef _abc_h_
#define _abc_h_

char AAUC[] = "LAGSVKETDIRPNFQYHMCWBZX-=", AALC[] =
"lagsvketdirpnfqyhmcwbzx-=";
int AtoA[] =
{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
 21, 22, 23};
char CHPN[] = "00000+-0-0+000000000000/=";
char CHCZ[] = "00000**0*0*000000000000/=";
char HPIC[] = "i000i**0*i*00i000i00000/=";
char H0UC[] = "0000000000000000H000000/=";

char NAUC[] = "TCAGYRSWKMN", NALC[] = "tcagyrswkmn", NAUCCP[] = "AGTCYRSWKMN";
int NtoN[]  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
char NARY[] = "RY", NARYlc[] = "ry";
int NtoRY[] = {1, 1, 0, 0};
char NASW[] = "SW", NASWlc[] = "sw";
int NtoSW[] = {1, 0, 1, 0};
char NAKM[] = "KM", NAKMlc[] = "km";
int NtoKM[] = {0, 1, 1, 0};

char NAVW[] = "Tcagyrswkmn";

int codtoaa[4][4][4] = {
 { {13, 13,  0,  0}, { 3,  3,  3,  3}, {15, 15, 23, 23}, {18, 18, 23, 19} },
 { { 0,  0,  0,  0}, {11, 11, 11, 11}, {16, 16, 14, 14}, {10, 10, 10, 10} },
 { { 9,  9,  9, 17}, { 7,  7,  7,  7}, {12, 12,  5,  5}, { 3,  3, 10, 10} },
 { { 4,  4,  4,  4}, { 1,  1,  1,  1}, { 8,  8,  6,  6}, { 2,  2,  2,  2} } };

#endif
