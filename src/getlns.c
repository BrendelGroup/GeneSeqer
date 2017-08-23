/* GETLNS.C;                                      Last update: April 3, 2013. */
/*   - a subroutine to read a nucleic acid sequence in library file format.   */
/* Dependencies:                                                              */
/* Bugs:                                                                      */

/*   Volker Brendel, Department of Biology                                    */
/*   Indiana University, Bloomington, IN 47405                                */
/*   vbrendel@indiana.edu                                                     */
/*   
   ; 07/23/00 changed by: (FKG) Fred Goodman, Promula Development Corporation
   ;                    : to use the TOPS database/suffix array algorithm
*/

#define OSYSFUNC					/* Gives access to the operating system */
#define STRGFUNC					/* Gives access to string manipulation */
#include "platform.h"
#include "RawTextFile.h"				/* Header, raw text file services */
#include "def.h"


int getlns(int LibraryFile, char detsz, char noffset, char *sfname, char *seq)
{
  auto char buf[LINELGTH];
  auto int n = 0, j, l = 0, nsflag = 0;
  auto long offset;
  int ch;

  for (j = 0; j < 257; ++j)
    sfname[j] = '\0';

  offset = getPositionRawTextFile(LibraryFile);
  if (readRawTextFile(LibraryFile, buf, LINELGTH) == NULL)
    return (0);
  while (buf[0] == '\n') {
    if (readRawTextFile(LibraryFile, buf, LINELGTH) == NULL)
      return (0);
  }

  if (buf[0] == '>') {
    for (j = noffset; buf[j] != '\n' && buf[j] != ':'; ++j) {
	if (buf[j] != ' '  &&  buf[j] != '	'  &&  buf[j] != '|') {
	  sfname[l++] = buf[j];
	  nsflag = 1;
	}
	else if (nsflag)
	  break;
    }
  }
  else {
    fprintf(stderr, "\nLibrary file format error:\n");
    fprintf(stderr, "%s", buf);
    return (0);
  }

  if (strlen(buf) == LINELGTH-1 && buf[LINELGTH-2] != '\n') {
    while ((ch = getCharacterRawTextFile(LibraryFile)) != '\n' && ch != RAW_TEXT_EOF) continue;
  }

  if (detsz) {
    while ((ch = getCharacterRawTextFile(LibraryFile)) != '>' && ch != RAW_TEXT_EOF) {
      if (ch == 'T' || ch == 'U' || ch == 't' || ch == 'u') {
	++n;
	continue;
      }
      if (ch == 'C' || ch == 'c') {
	++n;
	continue;
      }
      if (ch == 'A' || ch == 'a') {
	++n;
	continue;
      }
      if (ch == 'G' || ch == 'g') {
	++n;
	continue;
      }
      if (ch == 'Y' || ch == 'y') {
	++n;
	continue;
      }
      if (ch == 'R' || ch == 'r') {
	++n;
	continue;
      }
      if (ch == 'S' || ch == 's') {
	++n;
	continue;
      }
      if (ch == 'W' || ch == 'w') {
	++n;
	continue;
      }
      if (ch == 'K' || ch == 'k') {
	++n;
	continue;
      }
      if (ch == 'M' || ch == 'm') {
	++n;
	continue;
      }
      if (ch == 'B' || ch == 'b') {
	++n;
	continue;
      }
      if (ch == 'D' || ch == 'd') {
	++n;
	continue;
      }
      if (ch == 'H' || ch == 'h') {
	++n;
	continue;
      }
      if (ch == 'I' || ch == 'i') {
	++n;
	continue;
      }
      if (ch == 'Q' || ch == 'q') {
	++n;
	continue;
      }
      if (ch == 'V' || ch == 'v') {
	++n;
	continue;
      }
      if (ch == 'X' || ch == 'x') {
	++n;
	continue;
      }
      if (ch == 'N' || ch == 'n') {
	++n;
	continue;
      }
    }
    if (n == 0) {
      backupRawTextFile(LibraryFile, 1);
      return (-1);
    }
    else {
      setPositionRawTextFile(LibraryFile, offset);
      return (n);
    }
  }

  while ((ch = getCharacterRawTextFile(LibraryFile)) != '>' && ch != RAW_TEXT_EOF) {
    if (ch == 'T' || ch == 'U' || ch == 't' || ch == 'u') {
      seq[n++] = 0;
      continue;
    }
    if (ch == 'C' || ch == 'c') {
      seq[n++] = 1;
      continue;
    }
    if (ch == 'A' || ch == 'a') {
      seq[n++] = 2;
      continue;
    }
    if (ch == 'G' || ch == 'g') {
      seq[n++] = 3;
      continue;
    }
    if (ch == 'Y' || ch == 'y') {
      seq[n++] = 10;
      continue;
    }
    if (ch == 'R' || ch == 'r') {
      seq[n++] = 10;
      continue;
    }
    if (ch == 'S' || ch == 's') {
      seq[n++] = 10;
      continue;
    }
    if (ch == 'W' || ch == 'w') {
      seq[n++] = 10;
      continue;
    }
    if (ch == 'K' || ch == 'k') {
      seq[n++] = 10;
      continue;
    }
    if (ch == 'M' || ch == 'm') {
      seq[n++] = 10;
      continue;
    }
    if (ch == 'B' || ch == 'b') {
      seq[n++] = 10;
      continue;
    }
    if (ch == 'D' || ch == 'd') {
      seq[n++] = 10;
      continue;
    }
    if (ch == 'H' || ch == 'h') {
      seq[n++] = 10;
      continue;
    }
    if (ch == 'I' || ch == 'i') {
      seq[n++] = 10;
      continue;
    }
    if (ch == 'Q' || ch == 'q') {
      seq[n++] = 10;
      continue;
    }
    if (ch == 'V' || ch == 'v') {
      seq[n++] = 10;
      continue;
    }
    if (ch == 'X' || ch == 'x') {
      seq[n++] = 10;
      continue;
    }
    if (ch == 'N' || ch == 'n') {
      seq[n++] = 10;
      continue;
    }
  }
  if (ch == '>')
    backupRawTextFile(LibraryFile, 1);
  if (n == 0)
    return (-1);
  else
    return (n);

}							/* end getlns() */
