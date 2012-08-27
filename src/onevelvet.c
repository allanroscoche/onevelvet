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
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#if defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
#include <uce-dirent.h>
#define Arc v_Arc
#else
#include <dirent.h>
#endif
#include <stdio.h>
#include <unistd.h>

#include "run.h"
#include "binarySequences.h"
#include "globals.h"
#include "readSet.h"


#define FASTQ 1
#define FASTA 2
#define FASTA_GZ 5
#define FASTQ_GZ 6
#define SAM 8
#define BAM 9
#define RAW 10
#define RAW_GZ 11
#define AUTO 12


static void printUsage()
{
	puts("Usage:");
	puts("./velveth directory hash_length {[-file_format][-read_type][-separate|-interleaved] filename1 [filename2 ...]} {...} [options]");
	puts("");
	puts("\tdirectory\t: directory name for output files");
	printf("\thash_length\t: EITHER an odd integer (if even, it will be decremented) <= %i (if above, will be reduced)\n", MAXKMERLENGTH);
	printf("\t\t\t: OR: m,M,s where m and M are odd integers (if not, they will be decremented) with m < M <= %i (if above, will be reduced)\n", MAXKMERLENGTH);
	puts("\t\t\t\tand s is a step (even number). Velvet will then hash from k=m to k=M with a step of s");
	puts("\tfilename\t: path to sequence file or - for standard input");	
	puts("");
	puts("File format options:");
	puts("\t-fasta\t-fastq\t-raw\t-fasta.gz\t-fastq.gz\t-raw.gz\t-sam\t-bam\t-fmtAuto");
        puts("\t(Note: -fmtAuto will detect fasta or fastq, and will try the following programs for decompression : gunzip, pbunzip2, bunzip2");
	puts("");
        puts("File layout options for paired reads (only for fasta and fastq formats):");
        puts("\t-interleaved\t: File contains paired reads interleaved in the one file (default)");
        puts("\t-separate\t: Read 2 separate files for paired reads");
        puts("");
	puts("Read type options:");
	puts("\t-short\t-shortPaired");
#if CATEGORIES <= 5
	Category cat; 
	for (cat = 2; cat <= CATEGORIES; cat++)
	    printf("\t-short%i\t-shortPaired%i\n", cat, cat);
#else
	puts("\t...");
	printf("\t-short%i\t-shortPaired%i\n", CATEGORIES - 1, CATEGORIES - 1);
	printf("\t-short%i\t-shortPaired%i\n", CATEGORIES, CATEGORIES);
#endif
	puts("\t-long\t-longPaired");
	puts("\t-reference");
	puts("");
	puts("Options:");
	puts("\t-strand_specific\t: for strand specific transcriptome sequencing data (default: off)");
	puts("\t-reuse_Sequences\t: reuse Sequences file (or link) already in directory (no need to provide original filenames in this case (default: off)");
	puts("\t-noHash\t\t\t: simply prepare Sequences file, do not hash reads or prepare Roadmaps file (default: off)");
	puts("\t-create_binary  \t: create binary CnyUnifiedSeq file (default: off)");
	puts("");
	puts("Synopsis:");
	puts("");
	puts("- Short single end reads:");
	puts("\tvelveth Assem 29 -short -fastq s_1_sequence.txt");
	puts("");
	puts("- Paired-end short reads (remember to interleave paired reads):");
	puts("\tvelveth Assem 31 -shortPaired -fasta interleaved.fna");
	puts("");
	puts("- Paired-end short reads (using separate files for the paired reads)");
	puts("\tvelveth Assem 31 -shortPaired -fasta -separate left.fa right.fa");
	puts("");
	puts("- Two channels and some long reads:");
	puts("\tvelveth Assem 43 -short -fastq unmapped.fna -longPaired -fasta SangerReads.fasta");
	puts("");
	puts("- Three channels:");
	puts("\tvelveth Assem 35 -shortPaired -fasta pe_lib1.fasta -shortPaired2 pe_lib2.fasta -short3 se_lib1.fa");
	puts("");
	puts("Output:");
	puts("\tdirectory/Roadmaps");
	puts("\tdirectory/Sequences");
	puts("\t\t[Both files are picked up by graph, so please leave them there]");
}


int main(int argc, char **argv)
{
	ReadSet *allSequences = NULL;
	SplayTable *splayTable;
	int hashLength, hashLengthStep, hashLengthMax, h;
	char *directory, *filename, *seqFilename, *baseSeqName, *buf;
	char * token;
	boolean double_strand = true;
	boolean noHash = false;
	boolean multiple_kmers = false;
	DIR *dir;

	//int argIndex = 1;
	int filetype = FASTA;
	Category cat = 0;
	IDnum sequenceIndex = 1;
	short short_var;
	ReferenceCoordinateTable * refCoords = newReferenceCoordinateTable();
	boolean reuseSequences = false;
	boolean separate_pair_files = false;


	
	ReadSet *sequences = NULL;
	RoadMapArray *rdmaps;
	PreGraph *preGraph;
	Graph *graph;
	char *graphFilename, *connectedGraphFilename,
	    *preGraphFilename, *roadmapFilename,
	    *lowCovContigsFilename, *highCovContigsFilename;
	double coverageCutoff = -1;
	double longCoverageCutoff = -1;
	double maxCoverageCutoff = -1;
	double expectedCoverage = -1;
	Coordinate minContigLength = -1;
	Coordinate minContigKmerLength;
	boolean *dubious = NULL;
	Coordinate insertLength[CATEGORIES];
	Coordinate insertLengthLong = -1;
	Coordinate std_dev[CATEGORIES];
	Coordinate std_dev_long = -1;
	short int accelerationBits = 24;
	boolean readTracking = false;
	boolean exportAssembly = false;
	boolean unusedReads = false;
	boolean estimateCoverage = false;
	boolean estimateCutoff = false;
	boolean exportAlignments = false;
	FILE *file;
	int arg_index, arg_int;
	double arg_double;
	char *arg;
	ShortLength *sequenceLengths = NULL;
	//Category cat;
	boolean scaffolding = true;
	int pebbleRounds = 1;
	long long longlong_var;
	//short int short_var;
	boolean exportFilteredNodes = false;
	int clean = 0;
	boolean conserveLong = false;
	boolean shadows[CATEGORIES];
	int coverageMask = 1;
	SequencesReader *seqReadInfo = NULL;


	setProgramName("onevelvet");

	if (argc < 4) {
		printf("velvet\n");
		printf("Version %i.%i.%2.2i\n", VERSION_NUMBER,
		       RELEASE_NUMBER, UPDATE_NUMBER);
		printf("\nCopyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)\n");
		printf("This is free software; see the source for copying conditions.  There is NO\n");
		printf("warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n");
		printf("Compilation settings:\n");
		printf("CATEGORIES = %i\n", CATEGORIES);
		printf("MAXKMERLENGTH = %i\n", MAXKMERLENGTH);
#ifdef _OPENMP
		puts("OPENMP");
#endif
#ifdef LONGSEQUENCES
		puts("LONGSEQUENCES");
#endif
#ifdef BIGASSEMBLY
		puts("BIGASSEMBLY");
#endif
#ifdef COLOR
		puts("COLOR");
#endif
#ifdef DEBUG
		puts("DEBUG");
#endif
		printf("\n");
		printUsage();
		return 0;
	}

	token = strtok(argv[2], ",");
	hashLength = atoi(token);
	token = strtok(NULL, ",");
	if (token == NULL) {
		multiple_kmers = false;
		hashLengthMax = hashLength + 1;
	} else {
		multiple_kmers = true;
		hashLengthMax = atoi(token);
	}
	token = strtok(NULL, ",");
	if (token == NULL) {
		hashLengthStep = 2;
	} else {
		hashLengthStep = atoi(token);
	}

	if (hashLength > MAXKMERLENGTH) {
		velvetLog
		    ("Velvet can't handle k-mers as long as %i! We'll stick to %i if you don't mind.\n",
		     hashLength, MAXKMERLENGTH);
		hashLength = MAXKMERLENGTH;
	} 
	if (hashLength <= 0) {
		velvetLog("Invalid hash length: %s\n", argv[2]);
		printUsage();
		return 0;
	} 
	if (hashLength % 2 == 0) {
		velvetLog
		    ("Velvet can't work with even length k-mers, such as %i. We'll use %i instead, if you don't mind.\n",
		     hashLength, hashLength - 1);
		hashLength--;
	} 

	if (multiple_kmers) {
		if (hashLengthMax > MAXKMERLENGTH + 1) {
			velvetLog
			    ("Velvet can't handle k-mers as long as %i! We'll stick to %i if you don't mind.\n",
			     hashLengthMax, MAXKMERLENGTH + 1);
			hashLengthMax = MAXKMERLENGTH + 1;
		} 
		if (hashLengthMax <= hashLength) {
			velvetLog("hashLengthMin < hashLengthMax is required %s", argv[2]);
			printUsage();
			return 0;
		} 

		if (hashLengthStep <= 0) {
			velvetLog("Non-positive hash length! Setting it to 2\n");
			hashLengthStep = 2;
		}
		if (hashLengthStep % 2 == 1) {
			velvetLog
			    ("Velvet can't work with an odd length k-mer step, such as %i. We'll use %i instead, if you don't mind.\n",
			     hashLengthStep, hashLengthStep + 1);
			hashLengthStep++;
		}
	}
	
	// check if binary sequences should be used
	int argIndex;
	for (argIndex = 3; argIndex < argc; argIndex++)
		if (strcmp(argv[argIndex], "-create_binary") == 0)
			setCreateBinary(true);

	for (h = hashLength; h < hashLengthMax; h += hashLengthStep) {

		resetWordFilter(h);

		buf = mallocOrExit(2 * strlen(argv[1]) + 500, char);

		if ( multiple_kmers ) {
			sprintf(buf,"%s_%d",argv[1],h);
			directory = mallocOrExit(strlen(buf) + 100, char);
			strcpy(directory,buf);
		} else 
			directory = argv[1];

		filename = mallocOrExit(strlen(directory) + 100, char);
		seqFilename = mallocOrExit(strlen(directory) + 100, char);
		baseSeqName = mallocOrExit(100, char);

		dir = opendir(directory);

		if (dir == NULL)
			mkdir(directory, 0777);
		else {
			sprintf(buf, "%s/PreGraph", directory);
			remove(buf);
			sprintf(buf, "%s/Graph", directory);
			remove(buf);
			sprintf(buf, "%s/Graph2", directory);
			remove(buf);
			sprintf(buf, "%s/Graph3", directory);
			remove(buf);
			sprintf(buf, "%s/Graph4", directory);
			remove(buf);
		}

		logInstructions(argc, argv, directory);

		strcpy(seqFilename, directory);
		if (isCreateBinary()) {
			// use the CNY unified seq writer
			strcpy(baseSeqName, "/CnyUnifiedSeq");
			// remove other style sequences file
			sprintf(buf, "%s/Sequences", directory);
			remove(buf);
		} else {
			strcpy(baseSeqName, "/Sequences");
			// remove other style sequences file
			sprintf(buf, "%s/CnyUnifiedSeq", directory);
			remove(buf);
			sprintf(buf, "%s/CnyUnifiedSeq.names", directory);
			remove(buf);
		}
		strcat(seqFilename, baseSeqName);

		if ( h == hashLength ) {
		  //parseDataAndReadFiles(seqFilename, argc - 2, &(argv[2]), &double_strand, &noHash);
	
	if (argc < 4) {
		printUsage();
#ifdef DEBUG 
		abort();
#endif 
		exit(1);
	}
	char * filename_t = seqFilename;
	char ** argv_t = &(argv[2]);
	int argc_t = argc -2 ;
	for (argIndex = 1; argIndex < argc_t; argIndex++) {
		if (strcmp(argv_t[argIndex], "-strand_specific") == 0) {
			double_strand = false;
			reference_coordinate_double_strand = false;
		} else if (strcmp(argv_t[argIndex], "-reuse_Sequences") == 0) {
			reuseSequences = true;
		} else if (strcmp(argv_t[argIndex], "-noHash") == 0) {
			noHash = true;
		}
	}

	//if (reuseSequences) 
	//	return;

	SequencesWriter * seqWriteInfo = NULL;
	if (isCreateBinary()) {
		seqWriteInfo = openCnySeqForWrite(filename_t);
		seqWriteInfo->m_unifiedSeqFileHeader.m_bDoubleStrand = double_strand;
		// file is already open
	} else {
		seqWriteInfo = callocOrExit(1, SequencesWriter);
		seqWriteInfo->m_pFile = fopen(filename_t, "w");
	}
	int argIndex_t;
	for (argIndex_t = 1; argIndex_t < argc_t; argIndex_t++) {
		if (argv_t[argIndex_t][0] == '-' && strlen(argv_t[argIndex_t]) > 1) {

			if (strcmp(argv_t[argIndex_t], "-fastq") == 0)
				filetype = FASTQ;
			else if (strcmp(argv_t[argIndex_t], "-fasta") == 0)
				filetype = FASTA;
			else if (strcmp(argv_t[argIndex_t], "-fastq.gz") == 0)
				filetype = FASTQ_GZ;
			else if (strcmp(argv_t[argIndex_t], "-fasta.gz") == 0)
				filetype = FASTA_GZ;
			else if (strcmp(argv_t[argIndex_t], "-sam") == 0)
				filetype = SAM;
			else if (strcmp(argv_t[argIndex_t], "-bam") == 0)
				filetype = BAM;
			else if (strcmp(argv_t[argIndex_t], "-raw") == 0)
				filetype = RAW;
			else if (strcmp(argv_t[argIndex_t], "-raw.gz") == 0)
				filetype = RAW_GZ;
			else if (strcmp(argv_t[argIndex_t], "-fmtAuto") == 0)
				filetype = AUTO;
			else if (strcmp(argv_t[argIndex_t], "-short") == 0)
				cat = 0;
			else if (strcmp(argv_t[argIndex_t], "-shortPaired") ==
				 0)
				cat = 1;
			else if (strncmp
				 (argv_t[argIndex_t], "-shortPaired",
				  12) == 0) {
				sscanf(argv_t[argIndex_t], "-shortPaired%hd", &short_var);
				cat = (Category) short_var;
				if (cat < 1 || cat > CATEGORIES) {
					velvetLog("Unknown option: %s\n",
					       argv_t[argIndex_t]);
#ifdef DEBUG 
					abort();
#endif 
					exit(1);
				}
				cat--;
				cat *= 2;
				cat++;
			} else if (strncmp(argv_t[argIndex_t], "-short", 6) ==
				   0) {
				sscanf(argv_t[argIndex_t], "-short%hd", &short_var);
				cat = (Category) short_var;
				if (cat < 1 || cat > CATEGORIES) {
					velvetLog("Unknown option: %s\n",
					       argv_t[argIndex_t]);
#ifdef DEBUG 
					abort();
#endif 
					exit(1);
				}
				cat--;
				cat *= 2;
			} else if (strcmp(argv_t[argIndex_t], "-long") == 0)
				cat = LONG;		// CATEGORIES * 2;
			else if (strcmp(argv_t[argIndex_t], "-longPaired") == 0)
				cat = LONG_PAIRED;	// CATEGORIES * 2 + 1;
			else if (strcmp(argv_t[argIndex_t], "-reference") == 0)
				cat = REFERENCE;	// CATEGORIES * 2 + 2
			else if (strcmp(argv_t[argIndex_t], "-strand_specific") == 0) {
				double_strand = false;
				reference_coordinate_double_strand = false;
			} else if (strcmp(argv_t[argIndex_t], "-noHash") == 0) {
				;
			} else if (strcmp(argv_t[argIndex_t], "-create_binary") == 0) {
				;
			} else if (strcmp(argv_t[argIndex_t], "-interleaved") == 0) {
				separate_pair_files = false;
			} else if (strcmp(argv_t[argIndex_t], "-separate") == 0) {
				separate_pair_files = true;
			}
			else {
				velvetLog("velveth: Unknown option: %s\n",
				       argv_t[argIndex_t]);
#ifdef DEBUG 
				abort();
#endif 
				exit(1);
			}

			continue;
		}

		if (cat == -1)
			continue;

		
		switch (filetype) {
		case FASTA:
		case FASTQ:
		case FASTA_GZ:
		case FASTQ_GZ:
		case AUTO:
			// Separate files for paired reads?  Note odd categories used for paired read type
			
			if (separate_pair_files && cat%2==1) {
			  argIndex_t++;
			  if (argIndex_t>=argc_t)
			    exitErrorf(EXIT_FAILURE, false, "Require left & right filename for -separate mode");
			  readFastXPair(filetype, seqWriteInfo, argv_t[argIndex_t-1], argv_t[argIndex_t], cat, &sequenceIndex);
			} else {
			  readFastXFile(filetype, seqWriteInfo, argv_t[argIndex_t], cat, &sequenceIndex, refCoords);
			}
			break;
		case RAW:
			if (separate_pair_files && cat%2==1) {
                        	exitErrorf(EXIT_FAILURE, false, "Currently do not support -separate mode for RAW");
                        }
			readRawFile(seqWriteInfo, argv_t[argIndex_t], cat, &sequenceIndex);
			break;
		case RAW_GZ:
			if (separate_pair_files && cat%2==1) {
                        	exitErrorf(EXIT_FAILURE, false, "Currently do not support -separate mode for RAW");
                        }
			readRawGZFile(seqWriteInfo, argv_t[argIndex_t], cat, &sequenceIndex);
			break;
		case SAM:
			readSAMFile(seqWriteInfo, argv_t[argIndex_t], cat, &sequenceIndex, refCoords);
			break;
		case BAM:
			readBAMFile(seqWriteInfo, argv_t[argIndex_t], cat, &sequenceIndex, refCoords);
			break;
		default:
			velvetLog("Screw up in parser... exiting\n");
#ifdef DEBUG 
			abort();
#endif 
			exit(1);
		}

	}
	destroyReferenceCoordinateTable(refCoords);
	if (isCreateBinary()) {
		closeCnySeqForWrite(seqWriteInfo);
	} else {
		fclose(seqWriteInfo->m_pFile);
}
	if (seqWriteInfo) {
	    free(seqWriteInfo);
	}
		} else {
			sprintf(buf,"rm -f %s",seqFilename);
			if (system(buf)) {
				velvetLog("Command failed!\n");
				velvetLog("%s\n", buf);
#ifdef DEBUG
				abort();
#endif
				exit(1);
			}
			if (argv[1][0] == '/')
				sprintf(buf,"ln -s %s_%d%s %s",argv[1],hashLength,baseSeqName,seqFilename);
			else
				sprintf(buf,"ln -s `pwd`/%s_%d%s %s",argv[1],hashLength,baseSeqName,seqFilename);
			if (system(buf)) {
				velvetLog("Command failed!\n");
				velvetLog("%s\n", buf);
#ifdef DEBUG
				abort();
#endif
				exit(1);
			}
		}
		velvetLog("Até aqui!!!\n");


		if (noHash)
			continue;

		splayTable = newSplayTable(h, double_strand);
		if (isCreateBinary()) {
			allSequences = importCnyReadSet(seqFilename);
		} else {
			allSequences = importReadSet(seqFilename);
		}
		velvetLog("%li sequences in total.\n", (long) allSequences->readCount);

		strcpy(filename, directory);
		strcat(filename, "/Roadmaps");
		inputSequenceArrayIntoSplayTableAndArchive(allSequences,
							   splayTable, filename, seqFilename);

		destroySplayTable(splayTable);
		if (dir)
			closedir(dir);
		if (directory != argv[1])
			free(directory);
		free(filename);
		free(seqFilename);
		free(baseSeqName);
		free(buf);
		if (allSequences) {
			destroyReadSet(allSequences);
		}
	}

	
	/* --------------------------------------
	 *             velvetg
	 *  ------------------------------------*/
	 


	for (cat = 0; cat < CATEGORIES; cat++) {
		insertLength[cat] = -1;
		std_dev[cat] = -1;
		shadows[cat] = false;
	}

#ifdef _OPENMP
		puts("OPENMP");
#endif
#ifdef LONGSEQUENCES
		puts("LONGSEQUENCES");
#endif
#ifdef BIGASSEMBLY
		puts("BIGASSEMBLY");
#endif
#ifdef COLOR
		puts("COLOR");
#endif
#ifdef DEBUG
		puts("DEBUG");
#endif
	if (strcmp(argv[1], "--help") == 0) {
		printUsage();
		return 0;
	}

	// Memory allocation 
	directory = argv[1];
	graphFilename = mallocOrExit(strlen(directory) + 100, char);
	connectedGraphFilename = mallocOrExit(strlen(directory) + 100, char);
	preGraphFilename =
	    mallocOrExit(strlen(directory) + 100, char);
	roadmapFilename = mallocOrExit(strlen(directory) + 100, char);
	seqFilename = mallocOrExit(strlen(directory) + 100, char);
	lowCovContigsFilename = mallocOrExit(strlen(directory) + 100, char);
	highCovContigsFilename = mallocOrExit(strlen(directory) + 100, char);
	
	// Argument parsing
	for (arg_index = 4; arg_index < argc; arg_index++) {
		arg = argv[arg_index++];
		if (arg_index >= argc) {
			velvetLog("Unusual number of arguments!\n");
			//printUsage();
			printf("velvetg error\n");
#ifdef DEBUG 
			abort();
#endif 
			exit(1);
		}

		if (strcmp(arg, "-cov_cutoff") == 0) {
			if (strcmp(argv[arg_index], "auto") == 0) {
				estimateCutoff = true;
			} else {
				sscanf(argv[arg_index], "%lf", &coverageCutoff);
			}
		} else if (strcmp(arg, "-long_cov_cutoff") == 0) {
			sscanf(argv[arg_index], "%lf", &longCoverageCutoff);
		} else if (strcmp(arg, "-exp_cov") == 0) {
			if (strcmp(argv[arg_index], "auto") == 0) {
				estimateCoverage = true;
				readTracking = true;
			} else {
				sscanf(argv[arg_index], "%lf", &expectedCoverage);
				if (expectedCoverage > 0)
					readTracking = true;
			}
		} else if (strcmp(arg, "-ins_length") == 0) {
			sscanf(argv[arg_index], "%lli", &longlong_var);
			insertLength[0] = (Coordinate) longlong_var;
			if (insertLength[0] < 0) {
				velvetLog("Invalid insert length: %lli\n",
				       (long long) insertLength[0]);
#ifdef DEBUG 
				abort();
#endif 
				exit(1);
			}
		} else if (strcmp(arg, "-ins_length_sd") == 0) {
			sscanf(argv[arg_index], "%lli", &longlong_var);
			std_dev[0] = (Coordinate) longlong_var;
			if (std_dev[0] < 0) {
				velvetLog("Invalid std deviation: %lli\n",
				       (long long) std_dev[0]);
#ifdef DEBUG 
				abort();
#endif 
				exit(1);
			}
		} else if (strcmp(arg, "-ins_length_long") == 0) {
			sscanf(argv[arg_index], "%lli", &longlong_var);
			insertLengthLong = (Coordinate) longlong_var;
		} else if (strcmp(arg, "-ins_length_long_sd") == 0) {
			sscanf(argv[arg_index], "%lli", &longlong_var);
			std_dev_long = (Coordinate) longlong_var;
		} else if (strncmp(arg, "-ins_length", 11) == 0
			   && strchr(arg, 'd') == NULL) {
			sscanf(arg, "-ins_length%hi", &short_var);
			cat = (Category) short_var;
			if (cat < 1 || cat > CATEGORIES) {
				velvetLog("Unknown option: %s\n", arg);
#ifdef DEBUG 
				abort();
#endif 
				exit(1);
			}
			sscanf(argv[arg_index], "%lli", &longlong_var);
			insertLength[cat - 1] = (Coordinate) longlong_var;
			if (insertLength[cat - 1] < 0) {
				velvetLog("Invalid insert length: %lli\n",
				       (long long) insertLength[cat - 1]);
#ifdef DEBUG 
				abort();
#endif 
				exit(1);
			}
		} else if (strncmp(arg, "-ins_length", 11) == 0) {
			sscanf(arg, "-ins_length%hi_sd", &short_var);
			cat = (Category) short_var;
			if (cat < 1 || cat > CATEGORIES) {
				velvetLog("Unknown option: %s\n", arg);
#ifdef DEBUG 
				abort();
#endif 
				exit(1);
			}
			sscanf(argv[arg_index], "%lli", &longlong_var);
			std_dev[cat - 1] = (Coordinate) longlong_var;
			if (std_dev[cat - 1] < 0) {
				velvetLog("Invalid std deviation: %lli\n",
				       (long long) std_dev[cat - 1]);
#ifdef DEBUG 
				abort();
#endif 
				exit(1);
			}
		} else if (strcmp(arg, "-read_trkg") == 0) {
			readTracking =
			    (strcmp(argv[arg_index], "yes") == 0);
		} else if (strcmp(arg, "-scaffolding") == 0) {
			scaffolding =
			    (strcmp(argv[arg_index], "yes") == 0);
		} else if (strcmp(arg, "-exportFiltered") == 0) {
			exportFilteredNodes =
			    (strcmp(argv[arg_index], "yes") == 0);
		} else if (strcmp(arg, "-amos_file") == 0) {
			exportAssembly =
			    (strcmp(argv[arg_index], "yes") == 0);
		} else if (strcmp(arg, "-alignments") == 0) {
			exportAlignments =
			    (strcmp(argv[arg_index], "yes") == 0);
		} else if (strcmp(arg, "-min_contig_lgth") == 0) {
			sscanf(argv[arg_index], "%lli", &longlong_var);
			minContigLength = (Coordinate) longlong_var;
		} else if (strcmp(arg, "-coverage_mask") == 0) {
			sscanf(argv[arg_index], "%lli", &longlong_var);
			coverageMask = (IDnum) longlong_var;
		} else if (strcmp(arg, "-accel_bits") == 0) {
			sscanf(argv[arg_index], "%hi", &accelerationBits);
			if (accelerationBits < 0) {
				velvetLog
				    ("Illegal acceleration parameter: %s\n",
				     argv[arg_index]);
				printUsage();
				return -1;
			}
		} else if (strcmp(arg, "-max_branch_length") == 0) {
			sscanf(argv[arg_index], "%i", &arg_int);
			setMaxReadLength(arg_int);
			setLocalMaxReadLength(arg_int);
		} else if (strcmp(arg, "-max_divergence") == 0) {
			sscanf(argv[arg_index], "%lf", &arg_double);
			setMaxDivergence(arg_double);
			setLocalMaxDivergence(arg_double);
		} else if (strcmp(arg, "-max_gap_count") == 0) {
			sscanf(argv[arg_index], "%i", &arg_int);
			setMaxGaps(arg_int);
			setLocalMaxGaps(arg_int);
		} else if (strcmp(arg, "-min_pair_count") == 0) {
			sscanf(argv[arg_index], "%i", &arg_int);
			setUnreliableConnectionCutoff(arg_int);
		} else if (strcmp(arg, "-max_coverage") == 0) {
			sscanf(argv[arg_index], "%lf", &maxCoverageCutoff);
		} else if (strcmp(arg, "-long_mult_cutoff") == 0) {
			sscanf(argv[arg_index], "%i", &arg_int);
			setMultiplicityCutoff(arg_int);
		} else if (strcmp(arg, "-paired_exp_fraction") == 0) {
			sscanf(argv[arg_index], "%lf", &arg_double);
			setPairedExpFraction(arg_double);
		} else if (strcmp(arg, "-clean") == 0) {
			if (strcmp(argv[arg_index], "yes") == 0)
				clean = 1;
		} else if (strcmp(arg, "-very_clean") == 0) {
			if (strcmp(argv[arg_index], "yes") == 0)
				clean = 2;
		} else if (strcmp(arg, "-conserveLong") == 0) {
			if (strcmp(argv[arg_index], "yes") == 0)
				conserveLong = 2;
		} else if (strcmp(arg, "-unused_reads") == 0) {
			unusedReads =
			    (strcmp(argv[arg_index], "yes") == 0);
			if (unusedReads)
				readTracking = true;
		} else if (strcmp(arg, "-shortMatePaired") == 0) {
			shadows[0] = (strcmp(argv[arg_index], "yes") == 0);
		} else if (strncmp(arg, "-shortMatePaired", 16) == 0) {
			sscanf(arg, "-shortMatePaired%hi", &short_var);
			cat = (Category) short_var;
			if (cat < 1 || cat > CATEGORIES) {
				velvetLog("Unknown option: %s\n", arg);
#ifdef DEBUG
				abort();
#endif
				exit(1);
			}
			shadows[cat - 1] = (strcmp(argv[arg_index], "yes") == 0);
		} else if (strcmp(arg, "--help") == 0) {
			printUsage();
			return 0;
		}
		if ((strcmp(argv[argIndex], "-fastq") == 0) ||
		    (strcmp(argv[argIndex], "-fasta") == 0) ||
		    (strcmp(argv[argIndex], "-fastq.gz") == 0) ||
		    (strcmp(argv[argIndex], "-fasta.gz") == 0) ||
		    (strcmp(argv[argIndex], "-sam") == 0) ||
		    (strcmp(argv[argIndex], "-bam") == 0) ||
		    (strcmp(argv[argIndex], "-raw") == 0) ||
		    (strcmp(argv[argIndex], "-raw.gz") == 0) ||
		    (strcmp(argv[argIndex], "-fmtAuto") == 0) ||
		    (strcmp(argv[argIndex], "-short") == 0) ||
		    (strcmp(argv[argIndex], "-shortPaired") == 0) ||
		    (strncmp(argv[argIndex], "-short", 6) ==  0) ||
		    (strcmp(argv[argIndex], "-long") == 0) ||
		    (strcmp(argv[argIndex], "-longPaired") == 0) ||
		    (strcmp(argv[argIndex], "-reference") == 0) ||
		    (strcmp(argv[argIndex], "-strand_specific") == 0) ||
		    (strcmp(argv[argIndex], "-noHash") == 0) ||
		    (strcmp(argv[argIndex], "-create_binary") == 0) ||
		    (strcmp(argv[argIndex], "-interleaved") == 0) ||
		    (strcmp(argv[argIndex], "-separate") == 0) ){
		      ;
		    } else {
		    velvetLog("velvetg: Unknown option: %s;\n", arg);
		    //printUsage();
			return 1;
		  }
	}
	
	// Bookkeeping
	logInstructions(argc, argv, directory);

	seqReadInfo = callocOrExit(1, SequencesReader);
	strcpy(seqFilename, directory);
	// if binary CnyUnifiedSeq exists, use it.  Otherwise try Sequences
	strcat(seqFilename, "/CnyUnifiedSeq");
	if (access(seqFilename, R_OK) == 0) {
		seqReadInfo->m_bIsBinary = true;
	} else {
		seqReadInfo->m_bIsBinary = false;
		strcpy(seqFilename, directory);
	strcat(seqFilename, "/Sequences");
	}
	seqReadInfo->m_seqFilename = seqFilename;
	strcpy(roadmapFilename, directory);
	strcat(roadmapFilename, "/Roadmaps");

	strcpy(preGraphFilename, directory);
	strcat(preGraphFilename, "/PreGraph");

	strcpy(connectedGraphFilename, directory);
	strcat(connectedGraphFilename, "/ConnectedGraph");

	if (!readTracking) {
		strcpy(graphFilename, directory);
		strcat(graphFilename, "/Graph");
	} else {
		strcpy(graphFilename, directory);
		strcat(graphFilename, "/Graph2");
	}

	strcpy(lowCovContigsFilename, directory);
	strcat(lowCovContigsFilename, "/lowCoverageContigs.fa");

	strcpy(highCovContigsFilename, directory);
	strcat(highCovContigsFilename, "/highCoverageContigs.fa");

	// Graph uploading or creation
	if ((file = fopen(graphFilename, "r")) != NULL) {
		fclose(file);

		graph = importGraph(graphFilename);

	} else if ((file = fopen(connectedGraphFilename, "r")) != NULL) {
		fclose(file);
		if (seqReadInfo->m_bIsBinary) {

			sequences = importCnyReadSet(seqFilename);

#if 0
			// compare to velvet's version of a seq
			ReadSet *compareSequences = NULL;
			compareSeqFilename = mallocOrExit(strlen(directory) + 100, char);
			strcpy(compareSeqFilename, directory);
			strcat(compareSeqFilename, "/Sequences");
			compareSequences = importReadSet(compareSeqFilename);
			convertSequences(compareSequences);
			if (sequences->readCount != compareSequences->readCount) {
				printf("read count mismatch\n");
				exit(1);
			}
			int i;
			for (i = 0; i < sequences->readCount; i++) {
				TightString *tString = getTightStringInArray(sequences->tSequences, i);
				TightString *tStringCmp = getTightStringInArray(compareSequences->tSequences, i);
				if (getLength(tString) != getLength(tStringCmp)) {
					printf("sequence %d len mismatch\n", i);
					exit(1);
				}
				if (strcmp(readTightString(tString), readTightString(tStringCmp)) != 0) {
					printf("sequence %d cmp mismatch\n", i);
					printf("seq %s != cmp %s\n", readTightString(tString), readTightString(tStringCmp));
					exit(1);
				}
			}
#endif
		} else {
			sequences = importReadSet(seqFilename);
			convertSequences(sequences);
		}
		seqReadInfo->m_sequences = sequences;

		graph =
		    importConnectedGraph(connectedGraphFilename, sequences,
				   roadmapFilename, readTracking, accelerationBits);

		sequenceLengths =
		    getSequenceLengths(sequences, getWordLength(graph));
		correctGraph(graph, sequenceLengths, sequences->categories, conserveLong);
		exportGraph(graphFilename, graph, sequences->tSequences);
	} else if ((file = fopen(preGraphFilename, "r")) != NULL) {
		fclose(file);
		if (seqReadInfo->m_bIsBinary) {
			sequences = importCnyReadSet(seqFilename);
		} else {
		sequences = importReadSet(seqFilename);
		convertSequences(sequences);
		}
		seqReadInfo->m_sequences = sequences;
		graph =
		    importPreGraph(preGraphFilename, sequences,
				   roadmapFilename, readTracking, accelerationBits);
		sequenceLengths =
		    getSequenceLengths(sequences, getWordLength(graph));
		correctGraph(graph, sequenceLengths, sequences->categories, conserveLong);
		exportGraph(graphFilename, graph, sequences->tSequences);
	} else if ((file = fopen(roadmapFilename, "r")) != NULL) {
		fclose(file);

		rdmaps = importRoadMapArray(roadmapFilename);
		if (seqReadInfo->m_bIsBinary) {
			// pull in sequences first and use in preGraph
			sequences = importCnyReadSet(seqFilename);
			seqReadInfo->m_sequences = sequences;
#if 0
			// compare to velvet's version of a seq
			ReadSet *compareSequences = NULL;
			char *compareSeqFilename = mallocOrExit(strlen(directory) + 100, char);
			strcpy(compareSeqFilename, directory);
			strcat(compareSeqFilename, "/Sequences");
			compareSequences = importReadSet(compareSeqFilename);
			convertSequences(compareSequences);
			if (sequences->readCount != compareSequences->readCount) {
				printf("read count mismatch\n");
				exit(1);
			}
			int i;
			for (i = 0; i < sequences->readCount; i++) {
				TightString *tString = getTightStringInArray(sequences->tSequences, i);
				TightString *tStringCmp = getTightStringInArray(compareSequences->tSequences, i);
				if (getLength(tString) != getLength(tStringCmp)) {
					printf("sequence %d len mismatch\n", i);
					exit(1);
				}
				if (strcmp(readTightString(tString), readTightString(tStringCmp)) != 0) {
					printf("sequence %d cmp mismatch\n", i);
					printf("seq %s != cmp %s\n", readTightString(tString), readTightString(tStringCmp));
					exit(1);
				}
			}
			printf("sequence files match!\n");
#endif
		}
		preGraph = newPreGraph_pg(rdmaps, seqReadInfo);
		concatenatePreGraph_pg(preGraph);
		if (!conserveLong)
		    clipTips_pg(preGraph);
		exportPreGraph_pg(preGraphFilename, preGraph);
		destroyPreGraph_pg(preGraph);
		if (!seqReadInfo->m_bIsBinary) {
		sequences = importReadSet(seqFilename);
		convertSequences(sequences);
			seqReadInfo->m_sequences = sequences;
		}
		graph =
		    importPreGraph(preGraphFilename, sequences,
				   roadmapFilename, readTracking, accelerationBits);
		sequenceLengths =
		    getSequenceLengths(sequences, getWordLength(graph));
		correctGraph(graph, sequenceLengths, sequences->categories, conserveLong);
		exportGraph(graphFilename, graph, sequences->tSequences);
	} else {
		velvetLog("No Roadmap file to build upon! Please run velveth (see manual)\n");
#ifdef DEBUG 
		abort();
#endif 
		exit(1);
	}

	// Set insert lengths and their standard deviations
	for (cat = 0; cat < CATEGORIES; cat++) {
		if (insertLength[cat] > -1 && std_dev[cat] < 0)
			std_dev[cat] = insertLength[cat] / 10;
		setInsertLengths(graph, cat,
				 insertLength[cat], std_dev[cat]);
	}

	if (insertLengthLong > -1 && std_dev_long < 0)
		std_dev_long = insertLengthLong / 10;
	setInsertLengths(graph, CATEGORIES,
			 insertLengthLong, std_dev_long);

	// Coverage cutoff
	if (expectedCoverage < 0 && estimateCoverage == true) {
		expectedCoverage = estimated_cov(graph, directory);
		if (coverageCutoff < 0) {
			coverageCutoff = expectedCoverage / 2;
			estimateCutoff = true;
		}
	} else { 
		estimateCoverage = false;
		if (coverageCutoff < 0 && estimateCutoff) 
			coverageCutoff = estimated_cov(graph, directory) / 2;
		else 
			estimateCutoff = false;
	}

	if (coverageCutoff < 0) {
		velvetLog("WARNING: NO COVERAGE CUTOFF PROVIDED\n");
		velvetLog("Velvet will probably leave behind many detectable errors\n");
		velvetLog("See manual for instructions on how to set the coverage cutoff parameter\n");
	}

	if (sequences == NULL) {
		if (seqReadInfo->m_bIsBinary) {
			sequences = importCnyReadSet(seqFilename);
		} else {
		sequences = importReadSet(seqFilename);
		convertSequences(sequences);
	}
		seqReadInfo->m_sequences = sequences;
	}

	if (minContigLength < 2 * getWordLength(graph))
		minContigKmerLength = getWordLength(graph);
	else
		minContigKmerLength = minContigLength - getWordLength(graph) + 1;		

	dubious =
	    removeLowCoverageNodesAndDenounceDubiousReads(graph,
							  coverageCutoff,
							  sequences,
							  exportFilteredNodes,
							  minContigKmerLength,
							  lowCovContigsFilename);

	removeLowLongCoverageNodesAndDenounceDubiousReads(graph,
							  longCoverageCutoff,
							  sequences,
							  dubious,
							  exportFilteredNodes,
							  minContigKmerLength,
							  lowCovContigsFilename);

	removeHighCoverageNodes(graph, maxCoverageCutoff, exportFilteredNodes, minContigKmerLength, highCovContigsFilename);
	clipTipsHard(graph, conserveLong);

	if (sequences->readCount > 0 && sequences->categories[0] == REFERENCE)
		removeLowArcs(graph, coverageCutoff);

	if (expectedCoverage > 0) {

		// Mixed length sequencing
		readCoherentGraph(graph, isUniqueSolexa, expectedCoverage,
				  sequences);

		// Paired end resolution
		createReadPairingArray(sequences);
		pebbleRounds += pairedCategories(sequences);
		detachDubiousReads(sequences, dubious);
		activateGapMarkers(graph);

		for ( ;pebbleRounds > 0; pebbleRounds--)
			exploitShortReadPairs(graph, sequences, dubious, shadows, scaffolding);

	} else {
		velvetLog("WARNING: NO EXPECTED COVERAGE PROVIDED\n");
		velvetLog("Velvet will be unable to resolve any repeats\n");
		velvetLog("See manual for instructions on how to set the expected coverage parameter\n");
	}

	if (dubious)
		free(dubious);

	concatenateGraph(graph);

	removeLowCoverageReferenceNodes(graph, coverageCutoff, longCoverageCutoff, sequences);

	strcpy(graphFilename, directory);
	strcat(graphFilename, "/contigs.fa");
	sequenceLengths = getSequenceLengths(sequences, getWordLength(graph));
	exportLongNodeSequences(graphFilename, graph, minContigKmerLength, sequences, sequenceLengths, coverageMask); 

	if (exportAlignments) {
		strcpy(graphFilename, directory);
		strcat(graphFilename, "/contig-alignments.psa");
		exportLongNodeMappings(graphFilename, graph, sequences,
					     minContigKmerLength, seqReadInfo);
	}

	strcpy(graphFilename, directory);
	strcat(graphFilename, "/stats.txt");
	displayGeneralStatistics(graph, graphFilename, sequences);

	if (clean == 0) {
		strcpy(graphFilename, directory);
		strcat(graphFilename, "/LastGraph");
		exportGraph(graphFilename, graph, sequences->tSequences);
	}

	if (exportAssembly) {
		strcpy(graphFilename, directory);
		strcat(graphFilename, "/velvet_asm.afg");
		exportAMOSContigs(graphFilename, graph, minContigKmerLength, sequences);
	}

	if (unusedReads)
		exportUnusedReads(graph, sequences, minContigKmerLength, directory);

	if (estimateCoverage) 
		velvetLog("Estimated Coverage = %f\n", expectedCoverage);
	if (estimateCutoff) 
		velvetLog("Estimated Coverage cutoff = %f\n", coverageCutoff);

	logFinalStats(graph, minContigKmerLength, directory);

	if (clean > 0) {
		strcpy(graphFilename, directory);
		strcat(graphFilename, "/Roadmaps");
		remove(graphFilename);	

		strcpy(graphFilename, directory);
		strcat(graphFilename, "/LastGraph");
		remove(graphFilename);	
	} 

	if (clean > 1) {
		strcpy(graphFilename, directory);
		strcat(graphFilename, "/Sequences");
		remove(graphFilename);	

		strcpy(graphFilename, directory);
		strcat(graphFilename, "/Graph2");
		remove(graphFilename);	

		strcpy(graphFilename, directory);
		strcat(graphFilename, "/Graph");
		remove(graphFilename);	
	}

	free(sequenceLengths);
	destroyGraph(graph);
	free(graphFilename);
	free(connectedGraphFilename);
	free(preGraphFilename);
	free(seqFilename);
	free(roadmapFilename);
	free(lowCovContigsFilename);
	free(highCovContigsFilename);
	destroyReadSet(sequences);
	if (seqReadInfo) {
		free(seqReadInfo);
	}
	
	return 0;
}
