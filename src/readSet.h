/*
Copyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)

    This file is part of Velvet.

    Velvet is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Velvet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Velvet; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/
#ifndef _READSET_H_
#define _READSET_H_

#include "globals.h"
#include "tightString.h"
#include "readSet.h"
#include "utility.h"
#include "binarySequences.h"
#include "autoOpen.h"
#include "kseq.h"


struct readSet_st {
	char **sequences;
	TightString *tSequences;
	char **labels;
	char *tSeqMem;
	Quality **confidenceScores;
	Probability **kmerProbabilities;
	IDnum *mateReads;
	Category *categories;
	unsigned char *secondInPair;
	IDnum readCount;
};

struct mask_st {
	Coordinate start;
	Coordinate finish;
	Mask* next;
}; 

ReadSet *newReadSet();

ShortLength *getSequenceLengths(ReadSet * reads, int wordLength);

void convertSequences(ReadSet * rs);

ReadSet *importReadSet(char *filename);

// The overall argument parser and file reader for the hash function
void parseDataAndReadFiles(char * filename, int argc, char **argv, boolean * double_strand, boolean * noHash);

void logInstructions(int argc, char **argv, char *directory);

// Read pairing info
void createReadPairingArray(ReadSet* reads);
int pairedCategories(ReadSet * reads);
boolean isSecondInPair(ReadSet * reads, IDnum index);
void detachDubiousReads(ReadSet * reads, boolean * dubiousReads);

void destroyReadSet(ReadSet * reads);

inline boolean isCreateBinary();
void setCreateBinary(boolean val);



typedef struct referenceCoordinate_st ReferenceCoordinate;
static Coordinate reference_coordinate_double_strand = true;

struct referenceCoordinate_st {
	char * name;
	Coordinate start;
	Coordinate finish;
	IDnum referenceID;
	IDnum counter;
	boolean positive_strand;
}  ATTRIBUTE_PACKED;


typedef struct referenceCoordinateTable_st ReferenceCoordinateTable;

struct referenceCoordinateTable_st {
	ReferenceCoordinate * array;
	IDnum arrayLength;
}  ATTRIBUTE_PACKED;


ReferenceCoordinateTable * newReferenceCoordinateTable();
void destroyReferenceCoordinateTable(ReferenceCoordinateTable * table);
void readFastXFile(int fileType, SequencesWriter *seqWriteInfo, char *filename, Category cat, IDnum * sequenceIndex, ReferenceCoordinateTable * refCoords);
void readFastXPair(int fileType, SequencesWriter *seqWriteInfo, char *filename1, char *filename2, Category cat, IDnum * sequenceIndex);
void readSAMFile(SequencesWriter *seqWriteInfo, char *filename, Category cat, IDnum *sequenceIndex, ReferenceCoordinateTable * refCoords);
void readBAMFile(SequencesWriter *seqWriteInfo, char *filename, Category cat, IDnum *sequenceIndex, ReferenceCoordinateTable * refCoords);
void readRawFile(SequencesWriter *seqWriteInfo, char *filename, Category cat, IDnum * sequenceIndex);
void readRawGZFile(SequencesWriter *seqWriteInfo, char *filename, Category cat, IDnum *sequenceIndex);

#endif
