/* Modified by Benjamin Albrecht, Jan. 2012 */

/*Copyright (C) 2011  Fernando Meyer

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
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/time.h>
#include <limits.h>

#include "af.h"

int main(int argc, char *argv[]) {
	AffixArray   affixArray;
	Alphabet     alphabet;
	MultiPattern *multiPattern = NULL;
	bool         **complementarityRules = NULL;
	SearchParam  searchParam;

	int opt, longIndex = 0;

	bool argSuf     = false;
	bool argLcp     = false;
	bool argSkp     = false;
	bool argRpref   = false;
	bool argRlcp    = false;
	bool argRskp    = false;
	bool argAfsuf   = false;
	bool argAfrpref = false;

	/* Modified by Benjamin Albrecht, Jan. 2012  */
	bool argLcpTree = 0;
	bool argRlcpTree = 0;

	bool argShowTimes = false;
	bool argSilent1  = false;
	bool argSilent2  = false;
	bool anyTableSelected;

	int  predefAlphabet = 0;
	char *argAlphabetFile  = NULL;
    char *argLoadFile = NULL;
    char *argStruct   = NULL;
    char *argCompFile = NULL;

	searchParam.bSearchForwardString = false;
	searchParam.bSearchReverseString = false;
	searchParam.bBED                 = false;
	searchParam.bReportAllMatches    = false;
	searchParam.bIncludeSeqDesc      = false;
	searchParam.uiMinDiffMatches 	 = 0;
	searchParam.cTextFile            = NULL;

    // Chaining options
	Chainparam chainParam;
	chainParam.reportfile   = NULL;
	chainParam.isactive     = false;
	chainParam.islocal      = false;
	chainParam.globalparm   = NULL;
	chainParam.localparm    = NULL;
	chainParam.maxgap       = 0;
	chainParam.weightfactor = 1.0;
	chainParam.show         = false;
	chainParam.minchainlength = 0;
	chainParam.minchainscore = 0;
	chainParam.topkscoring   = UINT_MAX;

    extern int optind;
    static const char *optString = "h:kpqs:fbc:n:a124568mt:glw:x:e:r:odiju:v:@#z";

	static const struct option longOpts[] = {
		{"alph", required_argument, NULL, 'h'},
		{"dna", no_argument, NULL, 'k'},
		{"rna", no_argument, NULL, 'q'},
		{"protein", no_argument, NULL, 'p'},
		{"pat", required_argument, NULL, 's'},
		{"for", no_argument, NULL, 'f'},
		{"rev", no_argument, NULL, 'b'},
		{"comp", required_argument, NULL, 'c'},
		{"match", required_argument, NULL, 'n'},
		{"a", no_argument, NULL, 'a'},
		{"suf", no_argument, NULL, '1'},
		{"lcp", no_argument, NULL, '2'},
		{"aflk", no_argument, NULL, '4'},
		{"sufr", no_argument, NULL, '5'},
		{"lcpr", no_argument, NULL, '6'},
		{"aflkr", no_argument, NULL, '8'},

		/* Modified by Benjamin Albrecht, Jan. 2012  */
		{"lcpTree", no_argument, NULL, '9'},
		{"lcprTree", no_argument, NULL, '0'},

		{"bed", no_argument, NULL, '@'},
		{"allm", no_argument, NULL, 'm'},
		{"t", required_argument, NULL, 't'},
		{"seqdesc", no_argument, NULL, 'd'},
		{"time", no_argument, NULL, '#'},
		{"silent1", no_argument, NULL, 'i'},
		{"silent2", no_argument, NULL, 'j'},

		{"global", no_argument, NULL, 'g'},
		{"local", no_argument, NULL, 'l'},
		{"allglobal", no_argument, NULL, 'z'},
		{"wf", required_argument, NULL, 'w'},
		{"maxgap", required_argument, NULL, 'x'},
		{"minscore", required_argument, NULL, 'u'},
		{"minlen", required_argument, NULL, 'e'},
		{"top", required_argument, NULL, 'v'},
		{"chainrep", required_argument, NULL, 'r'},
		{"show", no_argument, NULL, 'o'},
		{NULL, no_argument, NULL, 0}
	};

	while((opt = getopt_long_only(argc, argv, optString, longOpts, &longIndex)) != -1) {
		switch(opt) {
			case 'h':
				argAlphabetFile = optarg;
				break;
			case 'k':
				predefAlphabet = 0;
				break;
			case 'q':
				predefAlphabet = 1;
				break;
			case 'p':
				predefAlphabet = 2;
				break;
			case 's':
				argStruct = optarg;
				break;
			case 'f':
				searchParam.bSearchForwardString = true;
				break;
			case 'b':
				searchParam.bSearchReverseString = true;
				break;
			case 'c':
				argCompFile = optarg;
				break;
			case 'n':
				searchParam.uiMinDiffMatches = atoi(optarg);
				break;
			case 'a':
				argSuf     = 1;
				argLcp     = 1;
				argRpref   = 1;
				argRlcp    = 1;
				argAfsuf   = 1;
				argAfrpref = 1;
				break;
			case '1':
				argSuf = 1;
				break;
			case '2':
				argLcp = 1;
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

			case '@':
				searchParam.bBED = true;
				break;
			case 'm':
				searchParam.bReportAllMatches = true;
				break;
			case 't':
				searchParam.cTextFile = optarg;
				break;
			case 'd':
				searchParam.bIncludeSeqDesc = true;
				break;
			case '#':
				argShowTimes = true;
				break;
			case 'i':
				argSilent1 = true;
				break;
			case 'j':
				argSilent2 = true;
				argSilent1 = true;
				break;
			case 'g':
				chainParam.islocal  = false;
				chainParam.isactive = true;
				break;
			case 'l':
				chainParam.islocal  = true;
				chainParam.isactive = true;
				break;
			case 'z':
				chainParam.globalparm = (char *) malloc(sizeof(char) * 10);
				strcpy(chainParam.globalparm, "all");
				break;
			case 'w':
				chainParam.weightfactor = strtod(optarg, NULL);
				break;
			case 'x':
				chainParam.maxgap = atoi(optarg);
				break;
			case 'u':
				chainParam.minchainscore = atoi(optarg);
				break;
			case 'e':
				chainParam.minchainlength = atoi(optarg);
				break;
			case 'v':
				chainParam.topkscoring = atoi(optarg);
				break;
			case 'r':
				chainParam.reportfile = optarg;
				chainParam.isactive = true;
				break;
			case 'o':
				chainParam.show = true;
				break;
			case '?':
				exit(1);
		}
	}

    if(argc < optind || argc == 1) {
	    fprintf(stderr,"Usage: %s [OPTION]\n", argv[0]);
	    fprintf(stderr,"  <data>                    Index name or FASTA file\n");
	    fprintf(stderr,"  -alph <file>              Use alphabet defined by file (option applies only to FASTA file)\n");
		fprintf(stderr,"  -dna                      Use 4-letter DNA alphabet (default) (option applies only to FASTA file)\n");
		fprintf(stderr,"  -rna                      Use 4-letter RNA alphabet (option applies only to FASTA file)\n");
		fprintf(stderr,"  -protein                  Use 20-letter protein alphabet (option applies only to FASTA file)\n");
	    fprintf(stderr,"  -pat <file>               Search for (structural) patterns\n");
	    fprintf(stderr,"  -for                      Search in the forward sequence (default)\n");
	    fprintf(stderr,"  -rev                      Search in the reverse complement sequence. For searching in the forward sequence as well, combine it with -for\n");
	    fprintf(stderr,"  -comp <file>              Load base pair complementarity rules from file\n");
	    fprintf(stderr,"  -a                        Map all index tables\n");
	    fprintf(stderr,"  -suf                      Map suf table\n");
	    fprintf(stderr,"  -lcp                      Map lcp table\n");

	    /* Modified by Benjamin Albrecht, Jan. 2012  */
	    fprintf(stderr,"  -lcpTree                  Map lcpTree table\n");

	    fprintf(stderr,"  -aflk                     Map aflk table\n");
	    fprintf(stderr,"  -sufr                     Map sufr table\n");
	    fprintf(stderr,"  -lcpr                     Map lcpr table\n");

	    /* Modified by Benjamin Albrecht, Jan. 2012  */
	    fprintf(stderr,"  -lcprTree                 Map lcprTree table\n");

	    fprintf(stderr,"  -aflkr                    Map aflkr table\n");
	    fprintf(stderr,"  -bed                      Output matches in BED format\n");
	    fprintf(stderr,"  -allm                     Report all matches of variable length patterns, i.e. not only the longest ones\n");
	    fprintf(stderr,"  -match <k>                Report only sequences matching at least k different patterns\n");
	    fprintf(stderr,"  -t <file>                 Write matches to text file instead of to screen\n");
	    fprintf(stderr,"  -seqdesc                  Include sequence description in the results, otherwise tag each pattern match with the sequence id\n");
	    fprintf(stderr,"  -time                     Display elapsed times\n");
	    fprintf(stderr,"  -silent1                  Do not output matches\n");
	    fprintf(stderr,"  -silent2                  Do not output anything\n");
	    fprintf(stderr,"                          \nChaining options:\n");
	    fprintf(stderr,"  -global                   Perform global chaining\n");
	    fprintf(stderr,"  -local                    Perform local chaining\n");
	    fprintf(stderr,"  -wf <wf>                  Apply weight factor > 0.0 to fragments\n");
	    fprintf(stderr,"  -maxgap <width>           Allow chain gaps with up to the specified width\n");
	    fprintf(stderr,"  -minscore <score>         Report only chains with at least the specified score\n");
	    fprintf(stderr,"  -minlen <length>          Report only chains with number of fragments >= length\n");
	    fprintf(stderr,"  -top <k>                  Report only top k scoring chains of each sequence\n");
	    fprintf(stderr,"  -allglobal                Report for each sequence all global chains satisfying above criteria\n");
	    fprintf(stderr,"  -chainrep <file>          Write chaining report to text file instead of to screen\n");
	    fprintf(stderr,"  -show                     Show chains in the report\n");
	    return 1;
    }
	if(optind < argc)
		argLoadFile = argv[optind];

	anyTableSelected = argSuf || argLcp || argSkp || argAfsuf || argRpref || argRlcp || argRskp || argAfrpref;

	/* Modified by Benjamin Albrecht, Jan. 2012 */
	init(&affixArray, &alphabet, argAlphabetFile, anyTableSelected ? -1 : predefAlphabet,
			&argSuf, &argLcp, &argLcpTree, &argSkp, &argAfsuf,
			&argRpref, &argRlcp, &argRlcpTree, &argRskp, &argAfrpref, false);

	if(argLoadFile == NULL) {
		fprintf(stderr,"Please specify an index or FASTA file.\n");
		return 1;
	}

	if(argSuf) {
		if (loadAffFiles(&affixArray, argLoadFile,
							argSuf, argLcp, argLcpTree, argSkp, argAfsuf,
							argRpref, argRlcp, argRlcpTree, argRskp, argAfrpref)) {
			return 1;
		}
	} else if(loadFastaFile(&affixArray, argLoadFile, false)) {
		return 1;
	}

	if(affixArray.multiSeq == NULL)
		return 1;

	if (!argSilent2) {
		printf("%cNumber of sequences: %i\n", LINESYMBOL, affixArray.multiSeq->numSeqs);
		printf("%cTotal length:        %i\n", LINESYMBOL, affixArray.length - affixArray.multiSeq->numSeqs);
	}

	if (affixArray.multiSeq->convSequences == NULL) {
		allocConvSequences(&affixArray);

		if (!argSilent2) {
			printf("\n%cPerforming alphabet conversion... ", LINESYMBOL);
			fflush(stdout);
		}

		if (convertToAlphabet(affixArray.multiSeq->sequences,
							affixArray.multiSeq->convSequences,
							affixArray.length,
							false,
							affixArray.alphabet))
			return 1;

		if (!argSilent2) printf("done\n");
	}

	if(argStruct != NULL) {
		multiPattern = loadPatternFile(argStruct, &affixArray);
		setSearchPatterns(multiPattern, searchParam.bSearchForwardString, searchParam.bSearchReverseString, affixArray.alphabet);

		if(argCompFile != NULL)
			complementarityRules = loadComplementarityFile(argCompFile, &affixArray);
		else
			complementarityRules = loadDefaultComplementarityRules(argSilent2, predefAlphabet, &affixArray);

		if(affixArray.esa == NULL)
			matchSlow(multiPattern,
					&searchParam,
					&chainParam, /*chaining*/
					complementarityRules,
					argShowTimes,
					argSilent1,
					argSilent2,
					&affixArray);
		else
			match(multiPattern,
					&searchParam,
					&chainParam, /*chaining*/
					complementarityRules,
					argShowTimes,
					argSilent1,
					argSilent2,
					&affixArray);
	}

	freeAll(&affixArray, multiPattern, complementarityRules);

	return 0;
}
