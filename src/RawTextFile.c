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
   ; name of module: TOPS.11 RawTextFile.c: Raw Text File Operational Entity
   ;
   ; The RawTextFile operational entity is used in those instances when text files
   ; imported from other platforms need to be processed in some special way,
   ; need to have their end-of-line conventions interrogated, or need to be 
   ; positioned across their end-of-line markers.
   ;
   ; The standard C runtime system allows two types of file: binary and text.
   ; Text files differ from binary files in that the runtime system checks
   ; for end-of-line marks and treats then as record boundaries. In binary files
   ; no such checking is performed. The problem is that not only do different
   ; platforms expect different text end-of-line characters, they have different
   ; numbers of end-of-line characters.
   ;
   ; Simply reading files with a foreign convention will usually work, but as soon
   ; as more complicated operations are performed, such as repositioning, working
   ; with foreign files breaks down. See the discussion of the service
   ; getPositionRawTextFile() for a detailed discussion of this problem.
   ;
   ; The RawTextFile operational entity allows the user to process text files that
   ; conform to any end-of-line convention; cr, crlf, lt, lfcr on any platform.
   ; When this operational entity opens an existing file it scans for the first
   ; end-of-line and then uses that convention through-out.
   ;
   ; Services Available from Operational Entity:
   ;
   ;  1. backupRawTextFile          Backup position in raw text file
   ;  2. closeRawTextFile           Close a raw text file
   ;  3. getCharacterRawTextFile    Get next character from raw text file
   ;  4. getPositionRawTextFile     Get position within raw text file
   ;  5. openRawTextFile            Open an existing raw text file
   ;  6. readRawTextFile            Read record from raw text file
   ;  7. setPositionRawTextFile     Set position within raw text file
   ;
   ;/hdr/ ************************************************************************
 */

#define MEMOFUNC					/* Gives access to memory allocation */
#define PRIOFUNC					/* Gives access to PROMULA I/O functions */
#include "platform.h"					/* Define the platform hosting this code */
#include "ServiceProvider.h"				/* Header, service provider services */
#include "RawTextFile.h"				/* Header, raw text file services */

#define BUFFER_SIZE        4096				/* Size of each file buffer */

#define END_OF_LINE_CR        0
#define END_OF_LINE_CRLF      1
#define END_OF_LINE_LF        2
#define END_OF_LINE_LFCR      3

#define LF_SYMBOL          0x0A
#define CR_SYMBOL          0x0D

#define NUMBER_OF_RAW_TEXT_FILES   10

typedef struct {
  binfile Handle;					/* Host platform handle to file */
  char *Buffer;						/* Points to file buffer */
  int EndOfLine;					/* End-of-line convention */
  int Position;						/* Current file position of start of buffer */
  int iBuffer;						/* Current position in buffer */
  int nBuffer;						/* Number of characters in current buffer */
} tRawTextFile;

static tRawTextFile RawTextFiles[NUMBER_OF_RAW_TEXT_FILES];

/*
   ;/doc/ ***********************************************************************
   ;
   ; TOPS.11.1 backupRawTextFile: Backup position in raw text file
   ;
   ; Synopsis of Service:
   ;
   ; #include "RawTextFile.h"             Header, raw text file services
   ;
   ; void backupRawTextFile               Backup position in raw text file
   ; (
   ;    int   Selection                   Selected provider to control file
   ;    int   Offset                      Number of positions to backup
   ; )
   ;
   ; Description of Service:
   ;
   ; This service repositions a raw text file by backing up its position within
   ; the current character record. The intent of this service is to allow the user
   ; to read individual characters from within a raw text file using the service
   ; getCharacterRawTextFile() in which a few character look-aheads might be
   ; necessary to simplify parsing. Once a look-ahead has occurred, this service
   ; can be used to backup the position indicator so that the next characters
   ; read will be the just read ones.
   ;
   ; This service will only work properly within records of the raw text file.
   ; Any attempt to cross a record boundary -- an end-of-line sequence --
   ; can cause the file position to become invalid. The problem is that
   ; text-records as read and written are logical structures and their character
   ; lengths do not necessarily correspond directly with their physical lengths.
   ;
   ; Properties of Service:
   ;
   ; Selection   The selected provider controlling the file being processed. It
   ;             must have been previously opened or created by a service of the
   ;             RawTextFile entity.
   ;
   ; Offset      The number of positions to move the current position back. It is
   ;             the callling services responcibility to ensure that this value
   ;             will not cause the traversal of and end-of-line marker.
   ;
   ; Return Value from Service: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
void backupRawTextFile(int Selection, int Offset)
#else
void backupRawTextFile(Selection, Offset)
  int Selection;
  int Offset;
#endif
{
  auto int Position;

   /*--------------------------------------------------------------------------
   ;
   ; Step 1: Obtain the current position in the raw text file, decrement it,
   ; and then reset the position. 
   ;
   ;-------------------------------------------------------------------------*/

  Position = getPositionRawTextFile(Selection);
  setPositionRawTextFile(Selection, Position - Offset);
  return;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; TOPS.11.2 closeRawTextFile: Close a raw text file
   ;
   ; Synopsis of Service:
   ;
   ; #include "RawTextFile.h"             Header, raw text file services
   ;
   ; void closeRawTextFile                Close a raw text file
   ; (
   ;    int   Selection                   Selected provider to control file
   ; )
   ;
   ; Description of Service:
   ;
   ; This service closes a raw text by returning any resources allocated to it
   ; to the host platform.
   ;
   ; Properties of Service:
   ;
   ; Selection   The selected provider controlling the file being processed. It
   ;             must have been previously opened or created by a service of the
   ;             RawTextFile entity.
   ;
   ; Return Value from Service: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
void closeRawTextFile(int Selection)
#else
void closeRawTextFile(Selection)
  int Selection;
#endif
{
  auto tRawTextFile *RawTextFile;

  RawTextFile = RawTextFiles + Selection - 1;
  clsbinf(RawTextFile->Handle);
  RawTextFile->Handle = nulbinf;
  free(RawTextFile->Buffer);
  RawTextFile->Buffer = NULL;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; TOPS.11.3 getCharacterRawTextFile: Get next character from raw text file
   ;
   ; Synopsis of Service:
   ;
   ; #include "RawTextFile.h"             Header, raw text file services
   ;
   ; int getCharacterRawTextFile          Get next character from raw text file
   ; (
   ;    int   Selection                   Selected provider to control file
   ; )
   ;
   ; Description of Service:
   ;
   ; This service reads a character from the raw text file provider indicated.
   ; Each character is read in initially as an unsigned character and is then
   ; converted to a non-negative integer before it is returned to the calling
   ; service. After reading the character this service advances the file
   ; position indicator. Note that regardless of the end-of-record conventions
   ; being used within the raw text file, and end-of-record is always returned
   ; as the value '\n'.
   ;
   ; As implemented this service is merely a simplified interface to the 
   ; readRawTextFile() service.
   ;
   ; This service is intended as a replacement for the C-function "fgetc()" and
   ; has the equivalent calling properties and return value.
   ;
   ; Properties of Service:
   ;
   ; Selection   The selected provider controlling the file being read. It must
   ;             have been previously opened or created by a service of the
   ;             RawTextFile entity.
   ;
   ; Return Value from Service:
   ;
   ; This service returns the character read from the raw text servive provider.
   ; If the file-position indicator is beyond the end of the file, then this
   ; service returns the negative constant RAW_TEXT_EOF.
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
int getCharacterRawTextFile(int Selection)
#else
int getCharacterRawTextFile(Selection)
  int Selection;
#endif
{
  auto UBYTE NextCharacter[2];
  auto char *ReadReturn;

  ReadReturn = readRawTextFile(Selection, (char *) (NextCharacter), 2);
  if (ReadReturn == NULL)
    return RAW_TEXT_EOF;
  return (int) (NextCharacter[0]);
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; TOPS.11.4 getPositionRawTextFile: Get position within raw text file
   ;
   ; Synopsis of Service:
   ;
   ; #include "RawTextFile.h"             Header, raw text file services
   ;
   ; int getPositionRawTextFile           Get position within raw text file
   ; (
   ;    int   Selection                   Selected provider to control file
   ; )
   ;
   ; Description of Service:
   ;
   ; As the raw text file is processed the service provider keeps track of its
   ; current position within the physical file -- i.e. the number of characters
   ; that current position is from the beginning. This service merely returns
   ; that position. In this context, all end-of-record marks are treated as having
   ; their actual lengths. Thus, for raw text files the sum of the number of
   ; characters read does not necessarily equal this value, since two-character
   ; end-of-line sequences are reduced to a single '\n' by the readRawTextFile()
   ; service.
   ;
   ; The intent of this service is to allow the user to obtain a position within
   ; a raw text file so that the file can be repositioned there at a later
   ; point. The return value of this service is compatible with the service
   ; setPositionRawTextFile() which does the repositioning.
   ;
   ; This service is intended as a replacement for the C-function "ftell()" and
   ; has the equivalent calling properties and return value. It is in fact the
   ; characteristics of the C-function ftell() that has necessitated the creation
   ; of the RawTextFile operational entity. As the Mark Williams Company says
   ; in their "ANSI C A Lexical Guide" (Prentice Hall 1988) about ftell(FILE* fp)
   ;
   ;  "If fp has been opened into text mode, however, ftell returns an 
   ;   implementation-defined number. For example, in UNIX-style environments,
   ;   ftell returns the number of characters the current position is from the
   ;   beginning; whereas under MS-DOS, where lines are terminated by a carriage
   ;   return-newline pair, ftell counts each carriage return and each newline
   ;   as a character in its return value."
   ;
   ; Obviously, the above makes porting codes that use ftell()-fseek() C-functions
   ; impossible. What the WIN32 platform (MS-DOS) does to implement the above is
   ; to subtract the number of lines read from the actual character offset
   ; position (UNIX ftell() value); thus, compensating for the extra character
   ; in each end-of-line. Unfortunately they do this even if the actual
   ; end-of-lines (which they recognize with no problem or warning) consist of
   ; a single character.
   ;
   ; Properties of Service:
   ;
   ; Selection   The selected provider controlling the file being processed. It
   ;             must have been previously opened or created by a service of the
   ;             RawTextFile entity.
   ;
   ; Return Value from Service:
   ;
   ; As described above, this service returns the UNIX-style position of the 
   ; raw text file.
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
int getPositionRawTextFile(int Selection)
#else
int getPositionRawTextFile(Selection)
  int Selection;
#endif
{
  auto tRawTextFile *RawTextFile;

   /*--------------------------------------------------------------------------
   ;
   ; Step 1: The "Position" member of control information contains the offset
   ; of the start of the buffer in the physical file and the "iBuffer" member
   ; contains the relative position within the buffer. The desired return value
   ; is simply the sum of these two values.
   ;
   ;------------------------------------------------------------------------*/

  RawTextFile = RawTextFiles + Selection - 1;
  return RawTextFile->Position + RawTextFile->iBuffer;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; TOPS.11.5 openRawTextFile: Open an existing raw text file
   ;
   ; Synopsis of Service:
   ;
   ; #include "RawTextFile.h"             Header, raw text file services
   ;
   ; int openRawTextFile                  Open an existing raw text file
   ; (
   ;    int   Selection                   Selected provider to control file
   ;    char* FileName                    Full name of file
   ;    int   ReadOnly                    Is the access to be read-only?
   ; )
   ;
   ; Description of Service:
   ;
   ; This service opens a text file potentially created on a foreign platform
   ; whose end-of-line conventions differ from the host platform that this service
   ; is running on.
   ;
   ; Properties of Service:
   ;
   ; Selection   The number of the selected RawTextFile provider. Note that it is
   ;             the responsibility of the calling service to determine or assign
   ;             this number. See the selectServiceProvider service for
   ;             information on this topic.
   ;
   ; FileName    The name of the file to be opened. It is the responsibility
   ;             of the calling service to ensure that this name corresponds to
   ;             the conventions of the host platform. This service merely passes
   ;             it on.
   ;
   ; ReadOnly    Indicates whether or not the file is to be opened for read
   ;             only access. A nonzero value requests read only; while a zero
   ;             value requests read-write.
   ;
   ; Return Value from Service:
   ;
   ; If all goes well, this function returns a zero. If the file cannot be
   ; opened or processed, then on the error codes ERR_NO_MEMORY,
   ; ERR_FILE_NOT_EXIST, or ERR_FILE_WRONG_TYPE is returned.
   ; 
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
int openRawTextFile(int Selection, char *FileName, int ReadOnly)
#else
int openRawTextFile(Selection, FileName, ReadOnly)
  int Selection;
  char *FileName;
  int ReadOnly;
#endif
{
  auto char *Buffer;
  auto int EndOfLine;
  auto binfile FileHandle;
  auto int iBuffer;
  auto int nRead;
  auto tRawTextFile *RawTextFile;

   /*--------------------------------------------------------------------------
   ;
   ; Step 1: Attempt to open the file with the appropriate access type. If this
   ; fails, return an ERR_FILE_NOT_EXIST to the calling service.
   ;
   ;-------------------------------------------------------------------------*/

  if (ReadOnly)
    FileHandle = rdobinf(FileName);
  else
    FileHandle = opnbinf(FileName);
  if (binopner(FileHandle))
    return ERR_FILE_NOT_EXIST;

   /*-------------------------------------------------------------------------
   ;
   ; Step 2: Allocate the buffer to be used to process the input from the file.
   ; If there is a problem, close the file and return an ERR_NO_MEMORY to the
   ; calling service.
   ;
   ;-------------------------------------------------------------------------*/

  Buffer = (char *) (getmem(BUFFER_SIZE));
  if (Buffer == NULL) {
    clsbinf(FileHandle);
    return ERR_NO_MEMORY;
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 3: Read in the first buffer of information from the file and then
   ; scan forward to find the first end-of-line marker. If the marker cannot
   ; be found in the first buffer, return all resources and return an
   ; ERR_BAD_FILE to the calling service.
   ;
   ;-------------------------------------------------------------------------*/

  nRead = rdbinf(FileHandle, Buffer, BUFFER_SIZE);
  for (iBuffer = 0; iBuffer < nRead; iBuffer++) {
    if (Buffer[iBuffer] == CR_SYMBOL || Buffer[iBuffer] == LF_SYMBOL)
      break;
  }
  if (iBuffer >= nRead) {
    free(Buffer);
    clsbinf(FileHandle);
    return ERR_BAD_FILE;
  }

   /*--------------------------------------------------------------------------
   ;
   ; Step 3: The first end-of-line has been found, classify it. Then initialize
   ; the control information for the raw text file and return a zero, meaning
   ; that all went well.
   ;
   ;-------------------------------------------------------------------------*/

  if (Buffer[iBuffer] == CR_SYMBOL) {
    if (Buffer[iBuffer + 1] == LF_SYMBOL)
      EndOfLine = END_OF_LINE_CRLF;
    else
      EndOfLine = END_OF_LINE_CR;
  }
  else if (Buffer[iBuffer + 1] == CR_SYMBOL)
    EndOfLine = END_OF_LINE_LFCR;
  else
    EndOfLine = END_OF_LINE_LF;
  RawTextFile = RawTextFiles + Selection - 1;
  RawTextFile->Handle = FileHandle;
  RawTextFile->Buffer = Buffer;
  RawTextFile->EndOfLine = EndOfLine;
  RawTextFile->Position = 0;
  RawTextFile->iBuffer = 0;
  RawTextFile->nBuffer = nRead;
  return 0;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; TOPS.11.6 readRawTextFile: Read record from raw text file
   ;
   ; Synopsis of Service:
   ;
   ; #include "RawTextFile.h"             Header, raw text file services
   ;
   ; char* readRawTextFile                Read record from raw text file
   ; (
   ;    int   Selection                   Selected provider to control file
   ;    char* Record                      Receives the record read
   ;    int   nRecord                     length of the record storage area
   ; )
   ;
   ; Description of Service:
   ;
   ; This service reads characters from the raw text file provider into the 
   ; return area "Record" until either "nrecord"-1 characters have been read,
   ; a newline of the appropriate type is read, or an end-of-file is encountered.
   ; A C-standard '\n' character is returned in the record at the position of
   ; the end-of-line regardless of what the representation of that end-of-line
   ; is in the raw text file. At the back of the record a null-character is
   ; appended.
   ;
   ; This service is intended as a replacement for the C-function "fgets()" and
   ; has the equivalent calling properties and return value.
   ;
   ; Properties of Service:
   ;
   ; Selection   The selected provider controlling the file being read. It must
   ;             have been previously opened or created by a service of the
   ;             RawTextFile entity.
   ;
   ; Record      Returns the actual record read with a '\n' in the position of
   ;             the end-of-line indicator. There is always a null-character 
   ;             inserted behind any information entered. If an EOF is encountered
   ;             immediately, then a null-character is placed in the initial
   ;             position.
   ;
   ; nRecord     The total size of the return area. Including the null-character,
   ;             this is the maximum number of characters that will be entered.
   ;
   ; Return Value from Service:
   ;
   ; This service returns the pointer "Record" if the read was completed 
   ; successfully. It returns a NULL if an end-of-file is encountered. This
   ; convention is selected to match the "fgets()" convention, since that is the
   ; C-function being replaced here.
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
char *readRawTextFile(int Selection, char *Record, int nRecord)
#else
char *readRawTextFile(Selection, Record, nRecord)
  int Selection;
  char *Record;
  int nRecord;
#endif
{
  auto char *Buffer;
  auto int Eol1;
  auto int Eol2;
  auto binfile FileHandle;
  auto int iBuffer;
  auto int iRecord;
  auto int nBuffer;
  auto int nRead;
  auto tRawTextFile *RawTextFile;

   /*-------------------------------------------------------------------------
   ;
   ; Step 1: Obtain the basic information needed from the control structure.
   ;
   ;------------------------------------------------------------------------*/

  RawTextFile = RawTextFiles + Selection - 1;
  FileHandle = RawTextFile->Handle;
  Buffer = RawTextFile->Buffer;
  nBuffer = RawTextFile->nBuffer;
  iBuffer = RawTextFile->iBuffer;

   /*-------------------------------------------------------------------------
   ;
   ; Step 3: Establish the actual end-of-line convention. Note that there 
   ; might be one or two characters. These could be stored in the control
   ; structure, but this appears to be equivalent.
   ;
   ;-------------------------------------------------------------------------*/

  switch (RawTextFile->EndOfLine) {
  case END_OF_LINE_CR:
    Eol1 = CR_SYMBOL;
    Eol2 = 0;
    break;
  case END_OF_LINE_CRLF:
    Eol1 = CR_SYMBOL;
    Eol2 = LF_SYMBOL;
    break;
  case END_OF_LINE_LF:
    Eol1 = LF_SYMBOL;
    Eol2 = 0;
    break;
  case END_OF_LINE_LFCR:
  default:
    Eol1 = LF_SYMBOL;
    Eol2 = CR_SYMBOL;
    break;
  }

   /*-------------------------------------------------------------------------
   ;
   ; Step 4: Initialize the outer loop which is controled by the record length
   ; minus 1. We must be certain that there will be room for the terminating
   ; null-character. If the record fills up, we will end reading before the
   ; end-of-line is encountered.
   ;
   ;-------------------------------------------------------------------------*/


  nRecord--;
  for (iRecord = 0; iRecord < nRecord; iRecord++) {

      /*-----------------------------------------------------------------------
      ;
      ; Step 5: We are ready to examine the next character in the buffer. Note
      ; that the first buffer-full of characters is read by the open service.
      ; If there are no more unexamined characters, terminate the record and
      ; return a NULL which is the end-of-file indicator for this service.
      ;
      ;----------------------------------------------------------------------*/

    if (iBuffer >= nBuffer) {
      RawTextFile->Position += nBuffer;
      nRead = rdbinf(FileHandle, Buffer, BUFFER_SIZE);
      if (nRead <= 0) {
	Record[iRecord] = '\0';
	RawTextFile->iBuffer = iBuffer;
	RawTextFile->nBuffer = 0;
	return NULL;
      }
      iBuffer = 0;
      nBuffer = nRead;
      RawTextFile->nBuffer = nRead;
    }

      /*-----------------------------------------------------------------------
      ;
      ; Step 6: If the character does not match the first end-of-line
      ; character, simply copy into the return record and continue to process
      ; the next character.
      ;
      ;----------------------------------------------------------------------*/

    if (Buffer[iBuffer] != Eol1) {
      Record[iRecord] = Buffer[iBuffer++];
      continue;
    }
    iBuffer++;

      /*----------------------------------------------------------------------
      ;
      ; Step 7: We have encountered the first end-or-line character. If there
      ; is no second end-of-line character, then we have completed a record.
      ; In this break out of loop to end local processing.
      ;
      ;----------------------------------------------------------------------*/

    if (Eol2 == 0) {
      Record[iRecord++] = '\n';
      break;
    }

      /*-----------------------------------------------------------------------
      ;
      ; Step 8: We must now to examine the next character in the buffer to
      ; see it it is the second end-of-line character. If there are no more
      ; unexamined characters, terminate the record and return a NULL
      ; which is the end-of-file indicator for this service.
      ;
      ;----------------------------------------------------------------------*/

    if (iBuffer >= nBuffer) {
      RawTextFile->Position += nBuffer;
      nRead = rdbinf(FileHandle, Buffer, BUFFER_SIZE);
      if (nRead <= 0) {
	Record[iRecord++] = (char) (Eol1);
	Record[iRecord] = '\0';
	RawTextFile->iBuffer = iBuffer;
	RawTextFile->nBuffer = 0;
	return NULL;
      }
      iBuffer = 0;
      nBuffer = nRead;
      RawTextFile->nBuffer = nRead;
    }

      /*----------------------------------------------------------------------
      ;
      ; Step 9: At this point we might have the case where we have a two
      ; character end-of-line, we have encountered the first character, but
      ; the next character is not the second one. If this happens then we
      ; try to store both characters in the return record, and then we 
      ; continue to process characters. Alternatively if we have matched both
      ; characters store a single '\n' in the return record and break out of
      ; the loop to end local processing.
      ;
      ;----------------------------------------------------------------------*/

    if (Buffer[iBuffer] != Eol2) {
      Record[iRecord++] = (char) (Eol1);
      if (iRecord >= nRecord)
	break;
      Record[iRecord] = Buffer[iBuffer++];
      continue;
    }

    Record[iRecord++] = '\n';
    iBuffer++;
    break;
  }
   /*--------------------------------------------------------------------------
   ;
   ; Step 10: Either we have filled our return record or we have found an
   ; end-of-line. Strore a terminating null character and return the non-NULL
   ; pointer to the start of the return record.
   ;
   ;-------------------------------------------------------------------------*/

  Record[iRecord] = '\0';
  RawTextFile->iBuffer = iBuffer;
  return Record;
}

/*
   ;/doc/ ***********************************************************************
   ;
   ; TOPS.11.7 setPositionRawTextFile: Set position within raw textFile
   ;
   ; Synopsis of Service:
   ;
   ; #include "RawTextFile.h"             Header, raw text file services
   ;
   ; void setPositionRawTextFile          Set position within raw text file
   ; (
   ;    int   Selection                   Selected provider to control file
   ;    int   Offset                      Relative offset of new position
   ; )
   ;
   ; Description of Service:
   ;
   ; This service repositions a raw text file on a specified offset, which is
   ; defined as the number of characters that the desired position is from the
   ; beginning of the raw text file. If at all possible, this value should be
   ; a value returned by the getPositionRawTextFile() service. The problem is
   ; that text-records as read and written are logical structures and their
   ; character lengths do not necessarily correspond directly with their
   ; physical lengths.
   ;
   ; The intent of this service is to allow the user to obtain a position within
   ; a raw text file so that the file can be repositioned there at a later
   ; point. This service is compatible with the service getPositionRawTextFile()
   ; which gets the position. 
   ;
   ; Properties of Service:
   ;
   ; Selection   The selected provider controlling the file being processed. It
   ;             must have been previously opened or created by a service of the
   ;             RawTextFile entity.
   ;
   ; Offset      An offset value returned by the service getPositionRawTextFile().
   ;             Do not do arithmetic on these values.
   ;
   ; Return Value from Service: None
   ;
   ;/doc/ ************************************************************************
 */
#ifdef FPROTOTYPE
void setPositionRawTextFile(int Selection, int Offset)
#else
void setPositionRawTextFile(Selection, Offset)
  int Selection;
  int Offset;
#endif
{
  auto tRawTextFile *RawTextFile;
  auto binfile FileHandle;
  auto int iBlock;
  auto int iOffset;

  RawTextFile = RawTextFiles + Selection - 1;
  iBlock = (Offset / BUFFER_SIZE) * BUFFER_SIZE;
  iOffset = Offset % BUFFER_SIZE;
  if (iBlock != RawTextFile->Position) {
    RawTextFile->Position = iBlock;
    FileHandle = RawTextFile->Handle;
    putbinf(FileHandle, iBlock);
    RawTextFile->nBuffer = rdbinf(FileHandle, RawTextFile->Buffer,
      BUFFER_SIZE);
  }
  RawTextFile->iBuffer = iOffset;
  return;
}
