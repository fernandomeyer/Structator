/* Modified by Benjamin Albrecht, Jan. 2012 */

/*Copyright (C) 2010  Fernando Meyer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>

#include <limits.h>

#include "af.h"

int main(int argc, char *argv[]) {
	AffixArray affixArray;
	Alphabet   alphabet;

	char *argFastaFile = NULL;
    char *argAlphabetFile  = NULL;
    int predefAlphabet = 0;
	bool argSuf     = 0;
	bool argLcp     = 0;
	bool argSkp     = 0;
	bool argRpref   = 0;
	bool argRlcp    = 0;
	bool argRskp    = 0;
	bool argAfsuf   = 0;
	bool argAfrpref = 0;

	/* Modified by Benjamin Albrecht, Jan. 2012  */
	bool argLcpTree = 0;
	bool argRlcpTree = 0;

	bool argSaveTransfSeq = 1;
	bool argShowTimes = 0;
	char *argIndex    = NULL;
	int argOScreen  = 0;
	char *argTextFile = NULL;

	int opt, longIndex = 0;
	extern int optind;
    static const char *optString = "z:drpa12345678s:xct:";

	static const struct option longOpts[] = {
		{"alph", required_argument, NULL, 'z'},
		{"dna", no_argument, NULL, 'd'},
		{"rna", no_argument, NULL, 'r'},
		{"protein", no_argument, NULL, 'p'},
		{"a", no_argument, NULL, 'a'},
		{"suf", no_argument, NULL, '1'},
		{"lcp", no_argument, NULL, '2'},
		{"skp", no_argument, NULL, '3'},
		{"aflk", no_argument, NULL, '4'},
		{"sufr", no_argument, NULL, '5'},
		{"lcpr", no_argument, NULL, '6'},
		{"skpr", no_argument, NULL, '7'},
		{"aflkr", no_argument, NULL, '8'},

		/* Modified by Benjamin Albrecht, Jan. 2012  */
		{"lcpTree", no_argument, NULL, '9'},
		{"lcprTree", no_argument, NULL, '0'},

		{"s", required_argument, NULL, 's'},
		{"x", no_argument, NULL, 'x'},
		{"c", no_argument, NULL, 'c'},
		{"t", required_argument, NULL, 't'},
		{"time", no_argument, NULL, 'b'},
		{NULL, no_argument, NULL, 0}
	};

	while((opt = getopt_long_only(argc, argv, optString, longOpts, &longIndex)) != -1) {
	    switch(opt) {
			case 'z':
				argAlphabetFile = optarg;
				break;
			case 'd':
				predefAlphabet = 0;
				break;
			case 'r':
				predefAlphabet = 1;
				break;
			case 'p':
				predefAlphabet = 2;
				break;
			case 'a':
				argSuf     = 1;
				argLcp     = 1;
				argSkp     = 1;
				argRpref   = 1;
				argRlcp    = 1;
				argRskp    = 1;
				argAfsuf   = 1;
				argAfrpref = 1;
				break;
			case '1':
				argSuf = 1;
				break;
			case '2':
				argLcp = 1;
				break;
			case '3':
				argSkp = 1;
				break;
			case '4':
				argAfsuf = 1;
				break;
			case '5':
				argRpref = 1;
				break;
			case '6':
				argRlcp = 1;
				break;
			case '7':
				argRskp = 1;
				break;
			case '8':
				argAfrpref = 1;
				break;

			/* Modified by Benjamin Albrecht, Jan. 2012  */
			case '9':
				argLcpTree = 1;
				break;
			case '0':
				argRlcpTree = 1;
				break;

			case 's':
				argIndex = optarg;
				break;
			case 'x':
				argSaveTransfSeq = 0;
				break;
			case 'c':
				argOScreen = 1;
				break;
			case 't':
				argTextFile = optarg;
				break;
			case 'b':
				argShowTimes = 1;
				break;
			case '?':
				exit(1);
		}
	}

    if(argc < optind || argc == 1) {
    	fprintf(stderr,"Usage: %s [OPTION]\n", argv[0]);
	    fprintf(stderr,"  <file>	 Load FASTA file\n");
	    fprintf(stderr,"  -alph <file>	 Use alphabet defined in file\n");
	    fprintf(stderr,"  -dna		 Use 4-letter DNA alphabet (default)\n");
	    fprintf(stderr,"  -rna		 Use 4-letter RNA alphabet\n");
	    fprintf(stderr,"  -protein	 Use 20-letter protein alphabet\n");
	    fprintf(stderr,"  -a    	 Construct all tables\n");
	    fprintf(stderr,"  -suf    	 Construct suf table\n");
	    fprintf(stderr,"  -lcp    	 Construct lcp table\n");

	    /* Modified by Benjamin Albrecht, Jan. 2012  */
	    fprintf(stderr,"  -lcpTree    	 Construct lcpTree table\n");

	    fprintf(stderr,"  -skp    	 Construct skp table\n");
	    fprintf(stderr,"  -aflk  	 Construct aflk table\n");
	    fprintf(stderr,"  -sufr   	 Construct sufr table\n");
	    fprintf(stderr,"  -lcpr    	 Construct lcpr table\n");

	    /* Modified by Benjamin Albrecht, Jan. 2012  */
	   	    fprintf(stderr,"  -lcprTree    	 Construct lcprTree table\n");

	    fprintf(stderr,"  -skpr    	 Construct skpr table\n");
	    fprintf(stderr,"  -aflkr	 Construct aflkr table\n");
	    fprintf(stderr,"  -s <index>	 Save constructed structures to given index name\n");
	    fprintf(stderr,"  -x		 Do not save alphabetically transformed sequence\n");
	    fprintf(stderr,"  -c		 Output constructed structures to screen\n");
	    fprintf(stderr,"  -t <file>	 Output constructed structures to text file\n");
	    fprintf(stderr,"  -time		 Display elapsed times\n");
	    return 1;
    }
	if(optind < argc)
		argFastaFile = argv[optind];

	 /* Modified by Benjamin Albrecht, Jan. 2012  */
    init(&affixArray, &alphabet, argAlphabetFile, predefAlphabet,
    		&argSuf, &argLcp, &argLcpTree, &argSkp, &argAfsuf,
			&argRpref, &argRlcp, &argRlcpTree, &argRskp, &argAfrpref, true);

    if(argFastaFile == NULL) {
		fprintf(stderr,"Please specify a FASTA file.\n");
		return 1;
    } else if(loadFastaFile(&affixArray, argFastaFile, true))
		return 1;

	if(affixArray.multiSeq == NULL)
		return 1;

	if(argFastaFile != NULL)
		printf("FASTA file:          %s\n", argFastaFile);
	printf("Number of sequences: %i\n", affixArray.multiSeq->numSeqs);
	printf("Total length:        %i\n", affixArray.length - affixArray.multiSeq->numSeqs);

	if(convertToAlphabet(affixArray.multiSeq->sequences,
						affixArray.multiSeq->convSequences,
						affixArray.length,
						false,
						&alphabet))
		return 1;

	if(argIndex != NULL && saveAffFiles(&affixArray, argIndex, argSaveTransfSeq))
		return 1;

	free(affixArray.multiSeq->sequences);

	/* Modified by Benjamin Albrecht, Jan. 2012  */
	if(affixArray.multiSeq != NULL && affixArray.esa->xarray == NULL &&
			(argSuf || argLcp || argLcpTree || argSkp || argAfsuf || argRpref || argRlcp || argRlcpTree || argRskp || argAfrpref)) {
		computeAffixArray(&affixArray,
							argSuf, argLcp, argLcpTree, argSkp, argAfsuf,
    						argRpref, argRlcp, argRlcpTree, argRskp, argAfrpref,
    						argOScreen, argIndex, argShowTimes);
	} else {
		printf("Please select the tables to be constructed.\n");
		exit(EXIT_FAILURE);
	}

	if(argOScreen) {
		printf("\n"); printEArray(&affixArray, 1);
		printf("\n"); printEArray(&affixArray, 0);
	}

	if(argTextFile != NULL && (saveTextFile(&affixArray, argTextFile, 1) || (saveTextFile(&affixArray, argTextFile, 0))))
		return 1;

	return 0;
}
