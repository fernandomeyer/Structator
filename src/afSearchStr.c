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
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <sys/time.h>

#include "af.h"
#include "redblack/redblack.h"

#ifndef DEBUG_ST_SEARCH
#define DEBUG_ST_SEARCH 0
#endif

int iMatches=0, iMatchesInTree=0;
unsigned int *uiNumMatchesSeq, *uiNumMatchesRSeq;
unsigned char *patternSequence, *patternStructure, *cBracketAlphabet = (unsigned char*) "()*";
int iPatternId, patternLength, hLoopStart, hLoopEnd, minhLoopStart, maxhLoopEnd, loopLength, minDepth, maxDepth, maxMispair=0;
bool **compRules;
bool bVariableLengthPattern, bSearchingFSeq, bReportAllMatches;
struct rbtree *rbF, *rbR;
Match **matchesVarLengthArray=NULL;

void setSearchPatterns(MultiPattern *multiPattern, bool argSearchForward, bool argSearchReverse, Alphabet *alphabet) {
	unsigned int i, j, strLength;
	char cA=1, cC=2, cG=3, cT=4;

	if (!argSearchForward && !argSearchReverse) {
		for (i = 0; i < multiPattern->numPatterns; i++) {
			multiPattern->pattern[i].seqR = NULL;
			multiPattern->pattern[i].structureR = NULL;
		}
		return;  // search in the forward seq (default)
	}

	if (argSearchReverse) {
		for (i = 0; i < alphabet->numClasses; i++) {
			switch (toupper(alphabet->classRepresentative[i])) {
				case 'A':
					cA = i + 1;
					break;
				case 'C':
					cC = i + 1;
					break;
				case 'G':
					cG = i + 1;
					break;
				case 'T':
					cT = i + 1;
					break;
				case 'U':
					cT = i + 1;
					break;
			}
		}

		for (i = 0; i < multiPattern->numPatterns; i++) {
			strLength = strlen((char*) multiPattern->pattern[i].seq);

			multiPattern->pattern[i].seqR = reverseStringNewPointer(multiPattern->pattern[i].seq, strLength);
			multiPattern->pattern[i].structureR = reverseStringNewPointer(multiPattern->pattern[i].structure, strLength);

			for (j = 0; j < strLength; j++) {
				switch (toupper(multiPattern->pattern[i].seqR[j])) {
					case 'R':
						multiPattern->pattern[i].seqR[j] = 'Y';
						break;
					case 'Y':
						multiPattern->pattern[i].seqR[j] = 'R';
						break;
					case 'S':
						//multiPattern->pattern[i].seqR[j] = 'S';
						break;
					case 'W':
						//multiPattern->pattern[i].seqR[j] = 'W';
						break;
					case 'K':
						multiPattern->pattern[i].seqR[j] = 'M';
						break;
					case 'M':
						multiPattern->pattern[i].seqR[j] = 'K';
						break;
					case 'B':
						multiPattern->pattern[i].seqR[j] = 'V';
						break;
					case 'V':
						multiPattern->pattern[i].seqR[j] = 'B';
						break;
					case 'D':
						multiPattern->pattern[i].seqR[j] = 'H';
						break;
					case 'H':
						multiPattern->pattern[i].seqR[j] = 'D';
						break;
					case 'N':
						multiPattern->pattern[i].seqR[j] = 'N';
						break;
					case '*':
						multiPattern->pattern[i].seqR[j] = 'N';
						break;
					default:
						if(multiPattern->pattern[i].seqR[j] == cA)
							multiPattern->pattern[i].seqR[j] = cT;
						else if(multiPattern->pattern[i].seqR[j] == cC)
							multiPattern->pattern[i].seqR[j] = cG;
						else if(multiPattern->pattern[i].seqR[j] == cG)
							multiPattern->pattern[i].seqR[j] = cC;
						else if(multiPattern->pattern[i].seqR[j] == cT)
							multiPattern->pattern[i].seqR[j] = cA;
						else {
							fprintf(stderr, "Error. Pattern %d contains a non-IUPAC character.\n", i);
							exit(1);
						}
				}

				if (multiPattern->pattern[i].structureR[j] == ')')
					multiPattern->pattern[i].structureR[j] = '(';
				else if (multiPattern->pattern[i].structureR[j] == '(')
					multiPattern->pattern[i].structureR[j] = ')';
			}
		}
	} else {
		for (i = 0; i < multiPattern->numPatterns; i++) {
			multiPattern->pattern[i].seqR = NULL;
			multiPattern->pattern[i].structureR = NULL;
		}
	}

	if (!argSearchForward) {
		for (i = 0; i < multiPattern->numPatterns; i++) {
			free(multiPattern->pattern[i].seq);
			free(multiPattern->pattern[i].structure);
			multiPattern->pattern[i].seq = NULL;
			multiPattern->pattern[i].structure = NULL;
		}
	}

}

#define swap(x, a, b) { tmpmatch = match[a]; \
			 match[a] = match[b]; match[b] = tmpmatch; }


void sortMatchesByPos(Match **match, int left, int right) {
	int i, last;
	Match *tmpmatch;

	if(left >= right)
		return;

	i = (left + right) >> 1; //i = (left + right)/2;
	swap(match, left, i);
	last = left;
	for(i = left+1; i <= right; i++){
		if(match[i]->pos < match[left]->pos) {
			++last;
			swap(array, last, i);
		}
	}
	swap(match, left, last);
	sortMatchesByPos(match, left, last-1);
	sortMatchesByPos(match, last+1, right);
}

/*Eliminate matches embedded in other matches. The matches must be
 * pre-sorted by matching position*/
unsigned int removeDuplicateMatches(Match **match) {
	unsigned int i, j;

	if (iMatches == 0)
		return 0;
	if (iMatches == 1)
		return 1;

	j = 0;
	for (i = 1; i < iMatches; i++) {
		if (match[i]->endPos < match[j]->endPos) {
			free(match[i]);
			match[i] = NULL;

		} else if (match[i]->endPos > match[j]->endPos) {
			if (match[i]->pos == match[j]->pos) {
				free(match[j]);
				while (j > 0 && match[i]->pos <= match[j - 1]->pos && match[i]->endPos > match[j - 1]->endPos) {
					j--;
					free(match[j]);
				}
				match[j] = match[i];
				match[i] = NULL;

			} else { // match[i]->pos > match[j]->pos
				j++;
				if (i > j) {
					match[j] = match[i];
					match[i] = NULL;
				}
			}

		} else { //match[i]->endPos == match[j]->endPos
			if (match[i]->pos == match[j]->pos) {
				j++;
				if (i > j) {
					match[j] = match[i];
					match[i] = NULL;
				}
			} else {
				free(match[i]);
				match[i] = NULL;
			}
		}
	}

	return j + 1; // Number of remaining matches
}

char* processStructureString(StructureInfo *structureInfo) {
	unsigned int ui, r, l, numLeftExtensions=0;
	char *structureString;

	if (structureInfo->isPairing != NULL)
		structureString = (char*) calloc(structureInfo->length + 1, sizeof(char));
	else
		return NULL;

	for (ui = 0; ui < structureInfo->length; ui++) {
		if (!structureInfo->isLeftRightExt[ui])
			numLeftExtensions++;
	}
	r = numLeftExtensions;
	l = r - 1;
	for (ui = 0; ui < structureInfo->length; ui++) {
		if(structureInfo->isLeftRightExt[ui])
			structureString[r++] = structureInfo->isPairing[ui] ? ')' : '.';
		else
			structureString[l--] = structureInfo->isPairing[ui] ? '(' : '.';
	}

	return structureString;
}


void reportMatch(VirtualNode *node, AffixArray *affixArray) {
	unsigned int ui, uiMatchPos, uiSeqId;
	Match *match;
	const void *rbval;
	char *structureString = processStructureString(&node->structureInfo);

	if(node->type == 0) {
		for(ui = node->i; ui <= node->j; ui++) {
			uiMatchPos    = affixArray->esa->xarray[ui] + node->context;
			uiSeqId = getSeqNumber(affixArray->multiSeq->seqEndPos, affixArray->multiSeq->numSeqs, uiMatchPos);

			match = (Match *) malloc(sizeof(Match));
			match->seqId     = uiSeqId;
			match->pos       = uiMatchPos;
			match->patternId = iPatternId;
			match->endPos    = uiMatchPos + (node->d - node->context) - 1;
			match->structureString = structureString;

			if (!bReportAllMatches) {
				if((iMatches + 1) % BUFFER1 == 0 && (matchesVarLengthArray = realloc(matchesVarLengthArray, (iMatches + BUFFER1 + 1) * sizeof(Match*))) == NULL) {
					fprintf(stderr,"Memory allocation failed - %s %d.\n", __FILE__, __LINE__);
					exit(1);
				}
				matchesVarLengthArray[iMatches] = match;
			} else {
				if (bSearchingFSeq) {
					rbval = rbsearch((void *) match, rbF);
					if (rbval != NULL) {
						uiNumMatchesSeq[uiSeqId]++;
						iMatchesInTree++;
					}
				} else {
					rbval = rbsearch((void *) match, rbR);
					if (rbval != NULL) {
						uiNumMatchesRSeq[uiSeqId]++;
						iMatchesInTree++;
					}
				}
			}
			iMatches++;
		}
	} else if(node->type == 1) {
		for(ui = node->i; ui <= node->j; ui++) {
			uiMatchPos    = (affixArray->length - 1) - affixArray->erpa->xarray[ui] - node->d;
			uiSeqId = getSeqNumber(affixArray->multiSeq->seqEndPos, affixArray->multiSeq->numSeqs, uiMatchPos);

			match = (Match *) malloc(sizeof(Match));
			match->seqId     = uiSeqId;
			match->pos       = uiMatchPos;
			match->patternId = iPatternId;
			match->endPos    = uiMatchPos + (node->d - node->context) - 1;
			match->structureString = structureString;

			if (!bReportAllMatches) {
				if((iMatches + 1) % BUFFER1 == 0 && (matchesVarLengthArray = realloc(matchesVarLengthArray, (iMatches + BUFFER1 + 1) * sizeof(Match*))) == NULL) {
					fprintf(stderr,"Memory allocation failed - %s %d.\n", __FILE__, __LINE__);
					exit(1);
				}
				matchesVarLengthArray[iMatches] = match;
			} else {
				if (bSearchingFSeq) {
					rbval = rbsearch((void *) match, rbF);
					if (rbval != NULL) {
						uiNumMatchesSeq[uiSeqId]++;
						iMatchesInTree++;
					}
				} else {
					rbval = rbsearch((void *) match, rbR);
					if (rbval != NULL) {
						uiNumMatchesRSeq[uiSeqId]++;
						iMatchesInTree++;
					}
				}
			}
			iMatches++;
		}
	} else if(node->type == 2 || node->type == 3) {
		uiSeqId = getSeqNumber(affixArray->multiSeq->seqEndPos, affixArray->multiSeq->numSeqs, node->i);

		match = (Match *) malloc(sizeof(Match));
		match->seqId     = uiSeqId;
		match->pos       = node->i;
		match->patternId = iPatternId;
		match->endPos    = node->j;
		match->structureString = structureString;

		if (!bReportAllMatches) {
			if((iMatches + 1) % BUFFER1 == 0 && (matchesVarLengthArray = realloc(matchesVarLengthArray, (iMatches + BUFFER1 + 1) * sizeof(Match*))) == NULL) {
				fprintf(stderr,"Memory allocation failed - %s %d.\n", __FILE__, __LINE__);
				exit(1);
			}
			matchesVarLengthArray[iMatches] = match;
		} else {
			if (bSearchingFSeq) {
				rbval = rbsearch((void *) match, rbF);
				if (rbval != NULL) {
					uiNumMatchesSeq[uiSeqId]++;
					iMatchesInTree++;
				}
			} else {
				rbval = rbsearch((void *) match, rbR);
				if (rbval != NULL) {
					uiNumMatchesRSeq[uiSeqId]++;
					iMatchesInTree++;
				}
			}

		}
		iMatches++;
	}

}

VirtualNode* childNode(VirtualNode *node, unsigned char *sExtension, int numExpansions, bool backwards, bool isPairingPos, AffixArray *affixArray) {
	VirtualNode *cNode = (VirtualNode *) malloc(sizeof(VirtualNode));
	unsigned int ui;

	cNode->context = node->context;
	cNode->d       = node->d;
	cNode->i       = node->i;
	cNode->j       = node->j;
	cNode->type    = node->type;
	cNode->structureInfo = node->structureInfo;

	if(!backwards) {
		if(extendRight(sExtension, numExpansions, cNode, affixArray)) {
			free(cNode);
			return NULL;
		}
	} else {
		if(extendLeft(sExtension, numExpansions, cNode, affixArray)) {
			free(cNode);
			return NULL;
		}
	}

	if (cNode->structureInfo.isPairing != NULL) {
		for (ui = 0; ui < numExpansions; ui++) {
			cNode->structureInfo.isPairing[cNode->structureInfo.length]      = isPairingPos;
			cNode->structureInfo.isLeftRightExt[cNode->structureInfo.length] = !backwards;
			cNode->structureInfo.length++;
		}
	}

	return cNode;
}

void match(MultiPattern *multiPattern, SearchParam *searchParam, Chainparam *chainParam,
		bool **complementarityRules, bool argShowTimes, bool bSilent1, bool bSilent2, AffixArray *affixArray) {

	VirtualNode node;
	Pattern *patternP = multiPattern->pattern;
	int numPatterns = multiPattern->numPatterns;
	int i, j, iOpeningBracket, k, iTotalMatches=0;
	unsigned int maxPatternLength;
	struct timeval start, end;
	struct timeval start_chain, end_chain;
	const void *rbval;
	struct rbtree *rbCR = NULL;
	bool bSeqSearched, bSeqRSearched, **revComplementarityRules = NULL;

	bReportAllMatches = searchParam->bReportAllMatches;

	if (!searchParam->bSearchForwardString && !searchParam->bSearchReverseString)
		searchParam->bSearchForwardString = true;

	if (chainParam->isactive) {
		if ((rbCR = rbinit(compareRBCR, NULL)) == NULL) { // Tree for chaining report
			fprintf(stderr, "Insufficient memory from rbinit() - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
	}

	if (searchParam->bSearchForwardString) {
		if ((uiNumMatchesSeq = (unsigned int *) calloc(affixArray->multiSeq->numSeqs, sizeof(unsigned int))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"uiNumMatchesSeq\" - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}

		for (i = 0; i < affixArray->multiSeq->numSeqs; i++) {
			uiNumMatchesSeq[i] = 0;
		}

		if ((rbF = rbinit(compareRB, NULL)) == NULL) {
			fprintf(stderr, "Insufficient memory from rbinit() - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
	} else {
		uiNumMatchesSeq = NULL;
		rbF = NULL;
	}
	if (searchParam->bSearchReverseString) {
		revComplementarityRules = loadReverseComplementarityFile(complementarityRules, affixArray);

		if ((uiNumMatchesRSeq = (unsigned int *) calloc(affixArray->multiSeq->numSeqs, sizeof(unsigned int))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"uiNumMatchesRSeq\" - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}

		for (i = 0; i < affixArray->multiSeq->numSeqs; i++) {
			uiNumMatchesRSeq[i] = 0;
		}

		if ((rbR = rbinit(compareRB, NULL)) == NULL) {
			fprintf(stderr, "Insufficient memory from rbinit() - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
	} else {
		uiNumMatchesRSeq = NULL;
		rbR = NULL;
	}

	bSeqSearched = bSeqRSearched = false;
	for(i = 0; i < numPatterns;) {
		if (patternP[i].seq != NULL && !bSeqSearched) {
			if (!bSilent2) printf("\n%cSearching for pattern %s in the forward sequence(s)... ", LINESYMBOL, patternP[i].desc);
			patternSequence  = patternP[i].seq;
			patternStructure = patternP[i].structure;

			compRules = complementarityRules;

			bSeqSearched = true;
			bSearchingFSeq = true;
		} else if (patternP[i].seqR != NULL && !bSeqRSearched) {
			if (!bSilent2) printf("\n%cSearching for pattern %s in the reverse complement sequence(s)... ", LINESYMBOL, patternP[i].desc);
			patternSequence  = patternP[i].seqR;
			patternStructure = patternP[i].structureR;

			compRules = revComplementarityRules;

			bSeqRSearched = true;
			bSearchingFSeq = false;
		} else {
			i++;
			bSeqSearched = bSeqRSearched = false;
			continue;
		}

		iPatternId = i;
		maxMispair = patternP[i].maxmispair;

		if (!bSilent2) fflush(stdout);

		if (argShowTimes) gettimeofday(&start, NULL);

		if (!bReportAllMatches) {
			if ((matchesVarLengthArray = (Match **) calloc(BUFFER1, sizeof(Match*))) == NULL) {
				fprintf(stderr,"Memory allocation failed for \"matchesVarLengthArray\" - %s %d.\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}

		node.context = 0;
		node.i = 0;
		node.j = affixArray->length - 1;
		node.d = 0;
		node.type = 0;

		iOpeningBracket = -1;
		patternLength = strlen((char*) patternStructure);

		j = 0;
		while(j < patternLength && patternStructure[j] != ')') {
			if(patternStructure[j] == '(')
				iOpeningBracket = j;
			j++;
		}

		if (iOpeningBracket != -1)
			hLoopStart = iOpeningBracket + 1;
		else
			hLoopStart = 0;
		hLoopEnd   = j - 1;

		/* Modified by Benjamin Albrecht, Jan. 2012 */
		/********* Find k = pos. of left-most non-wildcard pattern character in loop or loop start *********/
		j = hLoopStart;
		while(patternStructure[j] != ')' && (patternSequence[j] == 'N' || patternSequence[j] == '*') && j < patternLength) {
			j++;
		}
		k = j;
		if(k == (hLoopEnd + 1))
			k = hLoopStart;
		/********* End Find k = pos. of left-most non-wildcard pattern character in loop or loop start *********/

		/* Modified by Benjamin Albrecht, Jan. 2012 */
		/* Now a search can be performed without aflk tables! */
		if (k != 0  && (affixArray->esa->xarray == NULL || affixArray->esa->xlcp == NULL || //affixArray->esa->affixLink == NULL ||
				affixArray->erpa->xarray == NULL || affixArray->erpa->xlcp == NULL)){ //|| affixArray->erpa->affixLink == NULL)) {
			//fprintf(stderr, "\nPlease map tables suf, lcp, aflk, sufr, lcpr, aflkr.\n");
			fprintf(stderr, "\nPlease map tables suf, lcp, sufr, lcpr.\n");
			exit(1);
		} else if (affixArray->esa->xarray == NULL) {
			fprintf(stderr, "\nPlease map table suf.\n");
			exit(1);
		}

		loopLength = (hLoopEnd - hLoopStart + 1);

		if (bSearchingFSeq) {
			minhLoopStart = hLoopStart - patternP[i].maxleftloopextent;
			maxhLoopEnd   = hLoopEnd + patternP[i].maxrightloopextent;
		} else {
			minhLoopStart = hLoopStart - patternP[i].maxrightloopextent;
			maxhLoopEnd   = hLoopEnd + patternP[i].maxleftloopextent;
		}

		minDepth = 0;
		for (j = hLoopStart + 1; j < patternLength; j++) {
			if (patternStructure[j] == ')')
				minDepth++;
		}
		minDepth *= 2;
		if (patternP[i].maxstemlength > 0) {
			maxDepth = patternP[i].maxstemlength * 2;
			if (maxDepth < minDepth) {
				fprintf(stderr, "\nThe maximum stem length of pattern %s cannot be shorter than the minimum length\n", patternP[i].desc);
				exit(1);
			}
		} else {
			maxDepth = minDepth;
		}

		if (!chainParam->show && bSilent1) {
			node.structureInfo.isPairing      = NULL;
			node.structureInfo.isLeftRightExt = NULL;
		} else {
			maxPatternLength                  = maxhLoopEnd - minhLoopStart + maxDepth + patternLength;

			if ((node.structureInfo.isPairing = (bool*) calloc(maxPatternLength + 1, sizeof(bool))) == NULL) {
				fprintf(stderr,"Memory allocation failed for \"node.structureInfo.isPairing\" - %s %d.\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			if ((node.structureInfo.isLeftRightExt = (bool*) calloc(maxPatternLength + 1, sizeof(bool))) == NULL) {
				fprintf(stderr,"Memory allocation failed for \"node.structureInfo.isLeftRightExt\" - %s %d.\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}
		node.structureInfo.length = 0;

		matchHairpinLoop(k - 1,
						k,
						true,
						&node,
						affixArray);

		if (!bReportAllMatches) {
			sortMatchesByPos(matchesVarLengthArray, 0, iMatches - 1);
			iMatches = removeDuplicateMatches(matchesVarLengthArray);

			for (j = 0; j < iMatches; j++) {
				if (bSearchingFSeq) {
					rbval = rbsearch((void *) matchesVarLengthArray[j], rbF);
					if (rbval != NULL) {
						uiNumMatchesSeq[matchesVarLengthArray[j]->seqId]++;
						iMatchesInTree++;
					}
				} else {
					rbval = rbsearch((void *) matchesVarLengthArray[j], rbR);
					if (rbval != NULL) {
						uiNumMatchesRSeq[matchesVarLengthArray[j]->seqId]++;
						iMatchesInTree++;
					}
				}
			}
			free(matchesVarLengthArray);
		}

		if (node.structureInfo.isPairing != NULL)
			free(node.structureInfo.isPairing);
		if (node.structureInfo.isLeftRightExt != NULL)
			free(node.structureInfo.isLeftRightExt);

		if (argShowTimes) gettimeofday(&end, NULL);

		if (!bSilent2) {
			printf("done\n");
			if (argShowTimes) printf("%cTime:     %.4f ms\n", LINESYMBOL, (double)( (end.tv_sec - start.tv_sec)*1000 + (end.tv_usec - start.tv_usec)/1000.0) );
			printf("%c#Matches: %d\n", LINESYMBOL, iMatchesInTree);
		}

		iTotalMatches += iMatchesInTree;
		iMatches = iMatchesInTree = 0;
	}

	if (!bSilent2) printf("\n%c#Total matches: %d\n", LINESYMBOL, iTotalMatches);

	if(searchParam->uiMinDiffMatches > 0) {
		if (searchParam->bSearchForwardString)
			iTotalMatches -= removeNonsignificantMatches(rbF, uiNumMatchesSeq, searchParam->uiMinDiffMatches, numPatterns);
		if (searchParam->bSearchReverseString)
			iTotalMatches -= removeNonsignificantMatches(rbR, uiNumMatchesRSeq, searchParam->uiMinDiffMatches, numPatterns);
		if (!bSilent2) printf("%c#Total matches after removing unqualified sequences: %d\n", LINESYMBOL, iTotalMatches);
	}

	if (!bSilent1) saveSearchResultsRB(rbF, rbR, searchParam->cTextFile, searchParam->bBED,
			searchParam->bIncludeSeqDesc, patternP, numPatterns, affixArray, chainParam->isactive ? 0 : 1);

	if(chainParam->isactive) {
		if (!bSilent2) {
			printf("\n%cChaining matches... ", LINESYMBOL);
			fflush(stdout);

			if (argShowTimes)
				gettimeofday(&start_chain, NULL);
		}

		if (searchParam->bSearchForwardString)
			chainAll(rbF, rbCR, uiNumMatchesSeq, chainParam, multiPattern, affixArray, false, true);
		if (searchParam->bSearchReverseString)
			chainAll(rbR, rbCR, uiNumMatchesRSeq, chainParam, multiPattern, affixArray, false, false);

		if (!bSilent2) {
			printf("done\n");

			if (argShowTimes) {
				gettimeofday(&end_chain, NULL);
				printf("%cTime:     %.4f ms\n", LINESYMBOL, (double)( (end_chain.tv_sec - start_chain.tv_sec)*1000 + (end_chain.tv_usec - start_chain.tv_usec)/1000.0) );
			}

			printf("\n");
			reportChainingResults(rbCR, chainParam, affixArray);
		}
	}

	if (uiNumMatchesSeq != NULL)
		free(uiNumMatchesSeq);
	if (uiNumMatchesRSeq != NULL)
		free(uiNumMatchesRSeq);
	if (revComplementarityRules != NULL)
		free(revComplementarityRules);

	if (rbCR != NULL)
		rbdestroy(rbCR);
	if (rbF != NULL)
		rbdestroy(rbF);
	if (rbR != NULL)
		rbdestroy(rbR);
}


void matchHairpinLoop(int const pointerLeft, int const pointerRight, bool bAllowRightExt, VirtualNode *node, AffixArray *affixArray) {
	VirtualNode *cNode;
	int iAlphClass, numExpansions = 1;
	bool backwards, **iupacTable = affixArray->alphabet->iupacTable;
	unsigned char *currentSeq, *classRep = affixArray->alphabet->classRepresentative;
	unsigned char cSingleExtension[1];

	if (pointerRight <= maxhLoopEnd && bAllowRightExt) { // Right extension
		backwards = 0;

		if (pointerRight > hLoopEnd) {
			currentSeq = cBracketAlphabet + 2;
		} else {
			currentSeq = patternSequence + pointerRight;

			while ((numExpansions + pointerRight <= hLoopEnd) && (iupacTable[(int) *(currentSeq + numExpansions)] == NULL)) {
				numExpansions++;
			}
		}

		if (iupacTable[(int) *currentSeq] != NULL) {
			for(iAlphClass = 0; iAlphClass < affixArray->alphabet->numClasses; iAlphClass++) {
				if (affixArray->alphabet->isWildCard[iAlphClass])
					continue;

				if (!iupacTable[(int) *currentSeq][(int) classRep[iAlphClass]])
					continue;

				cSingleExtension[0] = iAlphClass + 1;
				if((cNode = childNode(node, cSingleExtension, 1, backwards, false, affixArray)) != NULL) {

					if(pointerLeft < hLoopStart && pointerRight + 1 > hLoopEnd) {
						if (hLoopEnd + 1 == patternLength) { //if pattern has no secondary structures
							reportMatch(cNode, affixArray);
						} else {
							matchStem(0, hLoopStart - 1, hLoopEnd + 1, 0, cNode, affixArray);
						}
					}

					if (pointerRight + 1 <= maxhLoopEnd)
						matchHairpinLoop(pointerLeft, pointerRight + 1, true, cNode, affixArray);
					if (pointerRight + 1 > hLoopEnd)
						matchHairpinLoop(pointerLeft, pointerRight + 1, false, cNode, affixArray);
					free(cNode);
				}
			}
		} else {
			if((cNode = childNode(node, currentSeq, numExpansions, backwards, false, affixArray)) != NULL) {

				if(pointerLeft < hLoopStart && pointerRight + numExpansions > hLoopEnd) {

					if (hLoopEnd + 1 == patternLength) { //if pattern has no secondary structures
						reportMatch(cNode, affixArray);
					} else {
						matchStem(0, hLoopStart - 1, hLoopEnd + 1, 0, cNode, affixArray);
					}
				}

				if (pointerRight + numExpansions <= maxhLoopEnd)
					matchHairpinLoop(pointerLeft, pointerRight + numExpansions, true, cNode, affixArray);
				if (pointerRight + numExpansions > hLoopEnd)
					matchHairpinLoop(pointerLeft, pointerRight + numExpansions, false, cNode, affixArray);
				free(cNode);
			}
		}
	} else if(pointerLeft >= minhLoopStart) { // Left extension
		backwards = 1;

		if (pointerLeft < hLoopStart)
			currentSeq = cBracketAlphabet + 2;
		else
			currentSeq = patternSequence + pointerLeft;

		if (iupacTable[(int) *currentSeq] != NULL) {
			for(iAlphClass = 0; iAlphClass < affixArray->alphabet->numClasses; iAlphClass++) {
				if (affixArray->alphabet->isWildCard[iAlphClass])
					continue;

				if (!iupacTable[(int) *currentSeq][(int) classRep[iAlphClass]])
					continue;

				cSingleExtension[0] = iAlphClass + 1;
				if((cNode = childNode(node, cSingleExtension, 1, backwards, false, affixArray)) != NULL) {

					if(pointerLeft - 1 < hLoopStart && pointerRight > hLoopEnd) {

						if (hLoopEnd + 1 == patternLength) { //if pattern has no secondary structures
							reportMatch(cNode, affixArray);
						} else {
							matchStem(0, hLoopStart - 1, hLoopEnd + 1, 0, cNode, affixArray);
						}
					}

					matchHairpinLoop(pointerLeft - 1, pointerRight, bAllowRightExt, cNode, affixArray);
					free(cNode);
				}
			}
		} else {
			if((cNode = childNode(node, currentSeq, numExpansions, backwards, false, affixArray)) != NULL) {

				if(pointerLeft - 1 < hLoopStart && pointerRight > hLoopEnd) {
					if (hLoopEnd + 1 == patternLength) { //if pattern has no secondary structures
						reportMatch(cNode, affixArray);
					} else {
						matchStem(0, hLoopStart - 1, hLoopEnd + 1, 0, cNode, affixArray);
					}
				}

				matchHairpinLoop(pointerLeft - numExpansions, pointerRight, bAllowRightExt, cNode, affixArray);
				free(cNode);
			}
		}

	}

}

void matchStem(int currentDepth, int stemPointerLeft, int stemPointerRight, int mispair, VirtualNode *node, AffixArray *affixArray) {
	int iAlphClass;
	bool backwards = ((currentDepth + 1) / 2) % 2, **iupacTable = affixArray->alphabet->iupacTable;
	VirtualNode *cNode;
	unsigned char c=0, cSingleExtension[1], *seq, *str, *classRep = affixArray->alphabet->classRepresentative;

	if (currentDepth > maxDepth) {
		return;
	} else if(currentDepth % 2 == 0 && currentDepth >= minDepth) {
		reportMatch(node, affixArray);

		if (currentDepth == maxDepth)
			return;
	}

	if (currentDepth >= minDepth) {
		seq = cBracketAlphabet + 2; // *
		str = cBracketAlphabet; // (
	} else if (backwards) {
		seq = patternSequence + stemPointerLeft;
		str = patternStructure + stemPointerLeft;
	} else {
		seq = patternSequence + stemPointerRight;
		str = patternStructure + stemPointerRight;
	}

	if(*str == '(' || *str == ')') {
		if(currentDepth % 2 == 0) {
			if (iupacTable[(int) *seq] != NULL) { //if(*seq == '*') {
				for(iAlphClass = 0; iAlphClass < affixArray->alphabet->numClasses; iAlphClass++) {
					if (affixArray->alphabet->isWildCard[iAlphClass])
						continue;

					if (!iupacTable[(int) *seq][(int) classRep[iAlphClass]])
						continue;

					cSingleExtension[0] = iAlphClass + 1;
					if(affixArray->alphabet->classRepresentative[iAlphClass] != '*' && (cNode = childNode(node, cSingleExtension, 1, backwards, true, affixArray)) != NULL) {
						if(backwards)
							matchStem(currentDepth + 1, stemPointerLeft - 1, stemPointerRight, mispair, cNode, affixArray);
						else
							matchStem(currentDepth + 1, stemPointerLeft, stemPointerRight + 1, mispair, cNode, affixArray);
						free(cNode);
					}
				}
			} else {
				if((cNode = childNode(node, seq, 1, backwards, true, affixArray)) != NULL) {
					if(backwards)
						matchStem(currentDepth + 1, stemPointerLeft - 1, stemPointerRight, mispair, cNode, affixArray);
					else
						matchStem(currentDepth + 1, stemPointerLeft, stemPointerRight + 1, mispair, cNode, affixArray);
					free(cNode);
				}
			}
		} else {
			/* Get character to check complementarity */
			if(node->type == 0) {
				if(!backwards)
					c = *(affixArray->multiSeq->convSequences + affixArray->esa->xarray[node->i] + node->context);
				else
					c = *(affixArray->multiSeq->convSequences + affixArray->esa->xarray[node->i] + node->d - 1);
			} else if(node->type == 1) {
				if(backwards)
					c = *(affixArray->multiSeq->convSequences + (affixArray->length - 2) - affixArray->erpa->xarray[node->i] - node->context);
				else
					c = *(affixArray->multiSeq->convSequences + (affixArray->length - 2) - affixArray->erpa->xarray[node->i] - (node->d - 1));
			} else if(node->type == 2) {
				if(!backwards)
					c = *(affixArray->multiSeq->convSequences + node->i);
				else
					c = *(affixArray->multiSeq->convSequences + node->j);
			} else if(node->type == 3) {
				if(backwards)
					c = *(affixArray->multiSeq->convSequences + node->j);
				else
					c = *(affixArray->multiSeq->convSequences + node->i);
			}
			/* End Get character to check complementarity */

			if (iupacTable[(int) *seq] != NULL) { //if(*seq == '*') {
				for(iAlphClass = 0; iAlphClass < affixArray->alphabet->numClasses; iAlphClass++)
					if(affixArray->alphabet->classRepresentative[iAlphClass] != '*') {
						if (!compRules[iAlphClass][c - 1] && mispair >= maxMispair)
							continue;

						if (affixArray->alphabet->isWildCard[iAlphClass])
							continue;

						if (!iupacTable[(int) *seq][(int) classRep[iAlphClass]])
							continue;

						cSingleExtension[0] = iAlphClass + 1;
						if((cNode = childNode(node, cSingleExtension, 1, backwards, true, affixArray)) != NULL) {
							if(backwards)
								matchStem(currentDepth + 1, stemPointerLeft - 1, stemPointerRight, compRules[iAlphClass][c - 1] ? mispair : mispair + 1, cNode, affixArray);
							else
								matchStem(currentDepth + 1, stemPointerLeft, stemPointerRight + 1, compRules[iAlphClass][c - 1] ? mispair : mispair + 1, cNode, affixArray);
							free(cNode);
						}
					}
			} else {
				if(compRules[*seq - 1][c - 1] || (!compRules[*seq - 1][c - 1] && mispair < maxMispair)) {
					if((cNode = childNode(node, seq, 1, backwards, true, affixArray)) != NULL) {
						if(backwards)
							matchStem(currentDepth + 1, stemPointerLeft - 1, stemPointerRight, compRules[*seq - 1][c - 1] ? mispair : mispair + 1, cNode, affixArray);
						else
							matchStem(currentDepth + 1, stemPointerLeft, stemPointerRight + 1, compRules[*seq - 1][c - 1] ? mispair : mispair + 1, cNode, affixArray);
						free(cNode);
					}
				}
			}
		}
	} else {
		if (iupacTable[(int) *seq] != NULL) { //if(*seq == '*') {
			for(iAlphClass = 0; iAlphClass < affixArray->alphabet->numClasses; iAlphClass++) {
				if (affixArray->alphabet->isWildCard[iAlphClass])
					continue;

				if (!iupacTable[(int) *seq][(int) classRep[iAlphClass]])
					continue;

				cSingleExtension[0] = iAlphClass + 1;
				if(/*affixArray->alphabet->classRepresentative[iAlphClass] != '*' &&*/ (cNode = childNode(node, cSingleExtension, 1, backwards, false, affixArray)) != NULL) {
					if(backwards)
						matchStem(currentDepth, stemPointerLeft - 1, stemPointerRight, mispair, cNode, affixArray);
					else
						matchStem(currentDepth, stemPointerLeft, stemPointerRight + 1, mispair, cNode, affixArray);
					free(cNode);
				}
			}
		} else {
			if((cNode = childNode(node, seq, 1, backwards, false, affixArray)) != NULL) {
				if(backwards)
					matchStem(currentDepth, stemPointerLeft - 1, stemPointerRight, mispair, cNode, affixArray);
				else
					matchStem(currentDepth, stemPointerLeft, stemPointerRight + 1, mispair, cNode, affixArray);
				free(cNode);
			}
		}
	}
}

void matchSlow(MultiPattern *multiPattern, SearchParam *searchParam, Chainparam *chainParam,
		bool **complementarityRules, bool argShowTimes, bool bSilent1, bool bSilent2, AffixArray *affixArray) {

	int i, j, iPattern, searchLimit, iOpeningBracket = 0, iTotalMatches=0;
	int pointerLeft, pointerRight, pointerLeftBak, pointerRightBak, pointerPLeft, pointerPRight, depth;
	struct timeval start, end;
	struct timeval start_chain, end_chain;
	Pattern *patternP = multiPattern->pattern;
	int numPatterns = multiPattern->numPatterns;
	char seqChar, extCharR, extCharL, extStrR, extStrL;
	unsigned int mispair = 0;
	bool doBreak, bSeqSearched, bSeqRSearched, **revComplementarityRules = NULL;
	Match *match;
	const void *rbval;
	struct rbtree *rbCR = NULL;
	bool **iupacTable = affixArray->alphabet->iupacTable;
	unsigned char *classRep = affixArray->alphabet->classRepresentative;

	unsigned int maxPatternLength;
	StructureInfo structureInfo;

	bReportAllMatches = searchParam->bReportAllMatches;

	if (!searchParam->bSearchForwardString && !searchParam->bSearchReverseString)
		searchParam->bSearchForwardString = true;

	if (searchParam->bSearchForwardString) {
		if ((uiNumMatchesSeq = (unsigned int *) calloc(affixArray->multiSeq->numSeqs, sizeof(unsigned int))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"uiNumMatchesSeq\" - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}

		for (i = 0; i < affixArray->multiSeq->numSeqs; i++) {
			uiNumMatchesSeq[i] = 0;
		}

		if ((rbF = rbinit(compareRB, NULL)) == NULL) {
			fprintf(stderr, "Insufficient memory from rbinit() - %s %d.\n", __FILE__, __LINE__);
			exit(1);
		}
	} else {
		uiNumMatchesSeq = NULL;
		rbF = NULL;
	}
	if (searchParam->bSearchReverseString) {
		revComplementarityRules = loadReverseComplementarityFile(complementarityRules, affixArray);

		if ((uiNumMatchesRSeq = (unsigned int *) calloc(affixArray->multiSeq->numSeqs, sizeof(unsigned int))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"uiNumMatchesRSeq\" - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}

		for (i = 0; i < affixArray->multiSeq->numSeqs; i++) {
			uiNumMatchesRSeq[i] = 0;
		}

		if ((rbR = rbinit(compareRB, NULL)) == NULL) {
			fprintf(stderr, "Insufficient memory from rbinit() - %s %d.\n", __FILE__, __LINE__);
			exit(1);
		}
	} else {
		uiNumMatchesRSeq = NULL;
		rbR = NULL;
	}

	int seqLength = affixArray->length;
	unsigned char *seq = affixArray->multiSeq->convSequences;

	if (chainParam->isactive) {
		if ((rbCR = rbinit(compareRBCR, NULL)) == NULL) { // Tree for chaining report
			fprintf(stderr, "Insufficient memory from rbinit() - %s %d.\n", __FILE__, __LINE__);
			exit(1);
		}
	}

	bSeqSearched = bSeqRSearched = false;
	for(iPattern = 0; iPattern < numPatterns;) {
		if (patternP[iPattern].seq != NULL && !bSeqSearched) {
			if (!bSilent2) printf("\n%cSearching for pattern %s in the forward sequence(s)... ", LINESYMBOL, patternP[iPattern].desc);
			patternSequence  = patternP[iPattern].seq;
			patternStructure = patternP[iPattern].structure;

			compRules = complementarityRules;

			bSeqSearched = true;
			bSearchingFSeq = true;
		} else if (patternP[iPattern].seqR != NULL && !bSeqRSearched) {
			if (!bSilent2) printf("\n%cSearching for pattern %s in the reverse complement sequence(s)... ", LINESYMBOL, patternP[iPattern].desc);
			patternSequence  = patternP[iPattern].seqR;
			patternStructure = patternP[iPattern].structureR;

			compRules = revComplementarityRules;

			bSeqRSearched = true;
			bSearchingFSeq = false;
		} else {
			iPattern++;
			bSeqSearched = bSeqRSearched = false;
			continue;
		}

		fflush(stdout);

		if (argShowTimes) gettimeofday(&start, NULL);

		if (!bReportAllMatches) {
			if ((matchesVarLengthArray = (Match **) calloc(BUFFER1, sizeof(Match*))) == NULL) {
				fprintf(stderr,"Memory allocation failed for \"matchesVarLengthArray\" - %s %d.\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}

		iOpeningBracket = -1;
		patternLength = strlen((char*) patternStructure);

		j = 0;
		while(j < patternLength && patternStructure[j] != ')') {
			if(patternStructure[j] == '(')
				iOpeningBracket = j;
			j++;
		}

		if (iOpeningBracket != -1)
			hLoopStart = iOpeningBracket + 1;
		else
			hLoopStart = 0;
		hLoopEnd   = j - 1;

		loopLength = (hLoopEnd - hLoopStart + 1);

		if (bSearchingFSeq) {
			minhLoopStart = hLoopStart - patternP[iPattern].maxleftloopextent;
			maxhLoopEnd   = hLoopEnd + patternP[iPattern].maxrightloopextent;
		} else {
			minhLoopStart = hLoopStart - patternP[iPattern].maxrightloopextent;
			maxhLoopEnd   = hLoopEnd + patternP[iPattern].maxleftloopextent;
		}

		minDepth = 0;
		for (j = hLoopStart + 1; j < patternLength; j++) {
			if (patternStructure[j] == ')')
				minDepth++;
		}
		minDepth *= 2;
		if (patternP[iPattern].maxstemlength > 0) {
			maxDepth = patternP[iPattern].maxstemlength * 2;
			if (maxDepth < minDepth) {
				fprintf(stderr, "\nThe maximum stem length of pattern %s cannot be shorter than the minimum length\n", patternP[iPattern].desc);
				exit(1);
			}
		} else {
			maxDepth = minDepth;
		}

		maxPatternLength             = maxhLoopEnd - minhLoopStart + maxDepth + patternLength;

		if ((structureInfo.isPairing = (bool*) calloc(maxPatternLength + 1, sizeof(bool))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"structureInfo.isPairing\" - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		if ((structureInfo.isLeftRightExt = (bool*) calloc(maxPatternLength + 1, sizeof(bool))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"structureInfo.isLeftRightExt\" - %s %d.\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}

		searchLimit = seqLength - patternLength + 1;

		for(i = 0; i < searchLimit; i++) {
			mispair              = 0;
			structureInfo.length = 0;

			for (pointerRight = hLoopStart; pointerRight <= maxhLoopEnd && pointerRight + i < seqLength; pointerRight++) {
				structureInfo.length = pointerRight - hLoopStart;

				seqChar = seq[i + pointerRight];

				if (pointerRight > hLoopEnd)
					extCharR = '*';
				else
					extCharR = patternSequence[pointerRight];


				if (seqChar == $ || affixArray->alphabet->isWildCard[seqChar - 1]) {
					break;
				}

				if (iupacTable[(int) extCharR] == NULL) {
					if (extCharR != seqChar)
						break;
				} else if (!iupacTable[(int) extCharR][(int) classRep[seqChar - 1]]) {
					break;
				}

				structureInfo.isPairing[structureInfo.length]      = false;
				structureInfo.isLeftRightExt[structureInfo.length] = true;
				structureInfo.length++;

				if (pointerRight >= hLoopEnd) {
					// LOOP MATCHED!
					for (pointerLeft = hLoopStart; pointerLeft >= minhLoopStart;) {

						// Make sure the last two extensions are unpairing
						structureInfo.length = pointerRight - pointerLeft + 1;
						if (structureInfo.length - 1 > 0) {
							structureInfo.isPairing[structureInfo.length - 1]      = false;
							structureInfo.isLeftRightExt[structureInfo.length - 1] = true;
						}
						structureInfo.isPairing[structureInfo.length]      = false;
						structureInfo.isLeftRightExt[structureInfo.length] = false;

						if (hLoopEnd + 1 == patternLength) {
							// MATCH!
							if ((match = (Match *) malloc(sizeof(Match))) == NULL) {
								fprintf(stderr,"Memory allocation failed for \"match\" - %s %d.\n", __FILE__, __LINE__);
								exit(EXIT_FAILURE);
							}
							match->seqId     = getSeqNumber(affixArray->multiSeq->seqEndPos, affixArray->multiSeq->numSeqs, i + pointerLeft);
							match->pos       = i + pointerLeft;
							match->patternId = iPattern;
							match->endPos    = i + pointerRight;
							match->structureString = processStructureString(&structureInfo);

							if (!bReportAllMatches) {
								if((iMatches + 1) % BUFFER1 == 0 && (matchesVarLengthArray = realloc(matchesVarLengthArray, (iMatches + BUFFER1 + 1) * sizeof(Match*))) == NULL) {
									fprintf(stderr,"Memory allocation failed - %s %d.\n", __FILE__, __LINE__);
									exit(1);
								}
								matchesVarLengthArray[iMatches] = match;
							} else {
								if (bSearchingFSeq) {
									rbval = rbsearch((void *) match, rbF);
									if (rbval != NULL) {
										uiNumMatchesSeq[match->seqId]++;
										iMatchesInTree++;
									}
								} else {
									rbval = rbsearch((void *) match, rbR);
									if (rbval != NULL) {
										uiNumMatchesRSeq[match->seqId]++;
										iMatchesInTree++;
									}
								}
							}
							iMatches++;

						} else {
							mispair = 0;

							// pointers for seq char extraction
							pointerLeftBak  = pointerLeft - 1;
							pointerRightBak = pointerRight + 1;
							if (i + pointerLeftBak < 0 || i + pointerRightBak >= seqLength)
								break;
							if(seq[i + pointerRightBak] > affixArray->alphabet->numClasses || affixArray->alphabet->isWildCard[seq[i + pointerRightBak] -1])
								break;

							// pointers for pattern char/str extraction
							pointerPLeft    = hLoopStart - 1;
							pointerPRight   = hLoopEnd + 1;

							depth = 0;
							while (depth < maxDepth) {
								if (depth >= minDepth) {
									extStrR = ')';
									extCharR = '*';
								} else {
									extStrR = patternStructure[pointerPRight];
									extCharR = patternSequence[pointerPRight];
								}

								if (iupacTable[(int) extCharR] == NULL) {
									if (extCharR != seq[i + pointerRightBak])
										break;
								} else if (!iupacTable[(int) extCharR][(int) classRep[seq[i + pointerRightBak] - 1]]) {
									break;
								}


								if (extStrR == ')') {

									if (depth >= minDepth) {
										extStrL = '(';
										extCharL = '*';
									} else {
										extStrL = patternStructure[pointerPLeft];
										extCharL = patternSequence[pointerPLeft];
									}
									doBreak = false;
									while (extStrL != '(') {
										if(seq[i + pointerLeftBak] > affixArray->alphabet->numClasses || affixArray->alphabet->isWildCard[seq[i + pointerLeftBak] -1]) {
											//i += iOpeningBracket;
											doBreak = true;
											break;
										}

										if (iupacTable[(int) extCharL] == NULL) {
											if (extCharL != seq[i + pointerLeftBak])
												break;
										} else if (!iupacTable[(int) extCharL][(int) classRep[seq[i + pointerLeftBak] - 1]]) {
											break;
										}

										structureInfo.isPairing[structureInfo.length]      = false;
										structureInfo.isLeftRightExt[structureInfo.length] = false;
										structureInfo.length++;

										pointerPLeft--; pointerLeftBak--;
										extStrL = patternStructure[pointerPLeft];
										extCharL = patternSequence[pointerPLeft];
									}
									if(doBreak || i + pointerLeftBak < 0 || seq[i + pointerLeftBak] > affixArray->alphabet->numClasses || affixArray->alphabet->isWildCard[seq[i + pointerLeftBak] -1]) {
										break;
									}

									if (iupacTable[(int) extCharL] == NULL) {
										if (extCharL != seq[i + pointerLeftBak])
											break;
									} else if (!iupacTable[(int) extCharL][(int) classRep[seq[i + pointerLeftBak] - 1]]) {
										break;
									}

									if(affixArray->alphabet->isWildCard[seq[i + pointerLeftBak] -1]) {
										break;
									}

									if(!compRules[seq[i + pointerRightBak] - 1][seq[i + pointerLeftBak] - 1]) {
										if (++mispair > patternP[iPattern].maxmispair)
											break;
									}

									structureInfo.isPairing[structureInfo.length]      = true;
									structureInfo.isLeftRightExt[structureInfo.length] = true;
									structureInfo.length++;

									structureInfo.isPairing[structureInfo.length]      = true;
									structureInfo.isLeftRightExt[structureInfo.length] = false;
									structureInfo.length++;

									depth += 2;

									if (depth >= minDepth) {
										// MATCH!
										if ((match = (Match *) malloc(sizeof(Match))) == NULL) {
											fprintf(stderr,"Memory allocation failed for \"match\" - %s %d.\n", __FILE__, __LINE__);
											exit(EXIT_FAILURE);
										}
										match->seqId     = getSeqNumber(affixArray->multiSeq->seqEndPos, affixArray->multiSeq->numSeqs, i + pointerLeftBak);
										match->pos       = i + pointerLeftBak;
										match->patternId = iPattern;
										match->endPos    = i + pointerRightBak;
										match->structureString  = processStructureString(&structureInfo);

										if (!bReportAllMatches) {
											if((iMatches + 1) % BUFFER1 == 0 && (matchesVarLengthArray = realloc(matchesVarLengthArray, (iMatches + BUFFER1 + 1) * sizeof(Match*))) == NULL) {
												fprintf(stderr,"Memory allocation failed - %s %d.\n", __FILE__, __LINE__);
												exit(1);
											}
											matchesVarLengthArray[iMatches] = match;
										} else {
											if (bSearchingFSeq) {
												rbval = rbsearch((void *) match, rbF);
												if (rbval != NULL) {
													uiNumMatchesSeq[match->seqId]++;
													iMatchesInTree++;
												}
											} else {
												rbval = rbsearch((void *) match, rbR);
												if (rbval != NULL) {
													uiNumMatchesRSeq[match->seqId]++;
													iMatchesInTree++;
												}
											}
										}
										iMatches++;
									}

									pointerPLeft--; pointerLeftBak--;

								} else {
									if(seq[i + pointerRightBak] > affixArray->alphabet->numClasses || affixArray->alphabet->isWildCard[seq[i + pointerRightBak] -1]) {
										break;
									}

									if (iupacTable[(int) extCharR] == NULL) {
										if (extCharR != seq[i + pointerRightBak])
											break;
									} else if (!iupacTable[(int) extCharR][(int) classRep[seq[i + pointerRightBak] - 1]]) {
										break;
									}

									structureInfo.isPairing[structureInfo.length]      = false;
									structureInfo.isLeftRightExt[structureInfo.length] = true;
									structureInfo.length++;
								}

								pointerPRight++; pointerRightBak++;
								if (i + pointerRightBak >= seqLength || seq[i + pointerRightBak] > affixArray->alphabet->numClasses || affixArray->alphabet->isWildCard[seq[i + pointerRightBak] -1])
									break;
							}


						}
						pointerLeft--;
						if (pointerLeft < minhLoopStart || i + pointerLeft < 0)
							break;

						seqChar = seq[i + pointerLeft];

						if (seqChar == $ || affixArray->alphabet->isWildCard[seqChar - 1]) {
							break;
						}

					}

				}

			}

		}

		if (!bReportAllMatches) {
			sortMatchesByPos(matchesVarLengthArray, 0, iMatches - 1);
			iMatches = removeDuplicateMatches(matchesVarLengthArray);

			for (j = 0; j < iMatches; j++) {
				if (bSearchingFSeq) {
					rbval = rbsearch((void *) matchesVarLengthArray[j], rbF);
					if (rbval != NULL) {
						uiNumMatchesSeq[matchesVarLengthArray[j]->seqId]++;
						iMatchesInTree++;
					}
				} else {
					rbval = rbsearch((void *) matchesVarLengthArray[j], rbR);
					if (rbval != NULL) {
						uiNumMatchesRSeq[matchesVarLengthArray[j]->seqId]++;
						iMatchesInTree++;
					}
				}
			}
			free(matchesVarLengthArray);
		}

		free(structureInfo.isPairing);
		free(structureInfo.isLeftRightExt);

		if (argShowTimes) gettimeofday(&end, NULL);

		if (!bSilent2)  {
			printf("done\n");
			if (argShowTimes) printf("%cTime:     %.4f ms\n", LINESYMBOL, (double)( (end.tv_sec - start.tv_sec)*1000 + (end.tv_usec - start.tv_usec)/1000.0) );
			printf("%c#Matches: %d\n", LINESYMBOL, iMatchesInTree);
		}

		iTotalMatches += iMatchesInTree;
		iMatches = iMatchesInTree = 0;
	}

	if (!bSilent2) printf("\n%c#Total matches: %d\n", LINESYMBOL, iTotalMatches);

	if(searchParam->uiMinDiffMatches > 0) {
		if (searchParam->bSearchForwardString)
			iTotalMatches -= removeNonsignificantMatches(rbF, uiNumMatchesSeq, searchParam->uiMinDiffMatches, numPatterns);
		if (searchParam->bSearchReverseString)
			iTotalMatches -= removeNonsignificantMatches(rbR, uiNumMatchesRSeq, searchParam->uiMinDiffMatches, numPatterns);
		if (!bSilent2) printf("%c#Total matches after removing unqualified sequences: %d\n", LINESYMBOL, iTotalMatches);
	}

	if (!bSilent1) saveSearchResultsRB(rbF, rbR, searchParam->cTextFile, searchParam->bBED,
			searchParam->bIncludeSeqDesc, patternP, numPatterns, affixArray, chainParam->isactive ? 0 : 1);

	if(chainParam->isactive) {
		if (!bSilent2) {
			printf("\n%cChaining matches... ", LINESYMBOL);
			fflush(stdout);

			if (argShowTimes)
				gettimeofday(&start_chain, NULL);
		}

		if (searchParam->bSearchForwardString)
			chainAll(rbF, rbCR, uiNumMatchesSeq, chainParam, multiPattern, affixArray, false, true);
		if (searchParam->bSearchReverseString)
			chainAll(rbR, rbCR, uiNumMatchesRSeq, chainParam, multiPattern, affixArray, false, false);

		if (!bSilent2) {
			printf("done\n");

			if (argShowTimes) {
				gettimeofday(&end_chain, NULL);
				printf("%cTime:     %.4f ms\n", LINESYMBOL, (double)( (end_chain.tv_sec - start_chain.tv_sec)*1000 + (end_chain.tv_usec - start_chain.tv_usec)/1000.0) );
			}

			printf("\n");
			reportChainingResults(rbCR, chainParam, affixArray);
		}
	}

	if (uiNumMatchesSeq != NULL)
		free(uiNumMatchesSeq);
	if (uiNumMatchesRSeq != NULL)
		free(uiNumMatchesRSeq);
	if (revComplementarityRules != NULL)
		free(revComplementarityRules);

	if (rbCR != NULL)
		rbdestroy(rbCR);
	if (rbF != NULL)
		rbdestroy(rbF);
	if (rbR != NULL)
		rbdestroy(rbR);
}
