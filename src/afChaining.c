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
#include <limits.h>
#include <string.h>
#include <stdio.h>

#include "af.h"
#include "redblack/redblack.h"

#define MAXVALUE 1000
#define MAXLENGTH 25
#define DEFAULTMAXGAP 0

// Global variables for chaining
Chainparam *chainParamGlob;
const void *rbvalCR;
struct rbtree *rbCR; // Tree for storing the output from chain2dim
int matchSeqId, *seqEndPosGlob;
bool bForwardSeqGlob, bShowStructureGlob = false;
struct rbtree *matchesRB;

typedef struct {
	unsigned long chaincounter;
} Counter;

static void chaining2gooutput(void *data, const GtChain2Dimmatchtable *matchtable, const GtChain2Dim *chain) {
	unsigned long idx, chainlength;
	int offset, chainscore;
	Counter *counter = (Counter *) data;
	Chain *chainingResult = NULL;
	GtChain2Dimmatchvalues value;
	Match match;

	const void *rbval;

	chainlength = gt_chain_chainlength(chain);
	chainscore  = gt_chain_chainscore(chain);

	if (counter->chaincounter >= chainParamGlob->topkscoring) {
		return;
	}

	if (chainlength < chainParamGlob->minchainlength || chainscore < chainParamGlob->minchainscore) {
		return;
	}

	chainingResult = (Chain *) malloc(sizeof(Chain));

	chainingResult->seqId        = matchSeqId;
	chainingResult->chainLength  = gt_chain_chainlength(chain);
	chainingResult->chainScore   = chainscore;
	chainingResult->forwardSeq   = bForwardSeqGlob;
	chainingResult->fragmentInfo = NULL;
	chainingResult->bStoredInReverseOrder = gt_chain_storedinreverseorder(chain);

	if ((chainingResult->fragmentInfo = calloc(chainingResult->chainLength, sizeof(ChainFragment))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"fragmentInfo\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	for (idx = 0; idx < chainlength; idx++) {
		gt_chain_extractchainelem(&value, matchtable, chain, idx);
		chainingResult->fragmentInfo[idx].startpos[0] = value.startpos[0];
		chainingResult->fragmentInfo[idx].endpos[0]   = value.endpos[0];
		chainingResult->fragmentInfo[idx].weight      = value.weight;

		if(matchSeqId > 0)
			offset = seqEndPosGlob[matchSeqId - 1] + 1;
		else
			offset = 0;

		if (bForwardSeqGlob) {
			chainingResult->fragmentInfo[idx].startpos[1] = value.startpos[1];
			chainingResult->fragmentInfo[idx].endpos[1]   = value.endpos[1];

			if (bShowStructureGlob) {
				match.pos    = chainingResult->fragmentInfo[idx].startpos[1] + offset;
				match.endPos = chainingResult->fragmentInfo[idx].endpos[1] + offset;
			}
		} else {
			chainingResult->fragmentInfo[idx].startpos[1] = (seqEndPosGlob[chainingResult->seqId] - 1 - value.endpos[1]) - offset;
			chainingResult->fragmentInfo[idx].endpos[1]   = (seqEndPosGlob[chainingResult->seqId] - 1 - value.startpos[1]) - offset;

			if (bShowStructureGlob) {
				match.pos    = chainingResult->fragmentInfo[idx].startpos[1];
				match.endPos = chainingResult->fragmentInfo[idx].endpos[1];
			}
		}

		chainingResult->fragmentInfo[idx].structureString = NULL;

		// Recover/find corresponding match containing the structure string
		if (bShowStructureGlob) {
			match.seqId      = matchSeqId;
			match.patternId  = -1; //info not available

			rbval = rbfind((void *)&match, matchesRB);

			if (rbval != NULL)
				chainingResult->fragmentInfo[idx].structureString = ((Match*)rbval)->structureString;
		}
	}

	// Insert chaining result into rb tree
	rbvalCR = rbsearch((void *) chainingResult, rbCR);
	// if(rbvalCR == NULL) identical chain is already in the tree

	counter->chaincounter++;
}

int** compute2DPositions (MultiPattern *multiPattern) {
	Pattern *patternP = multiPattern->pattern;
	int numPatterns = multiPattern->numPatterns;
	int *patternLength, **patternArtificialInt, i, j, count;
	int lowestInstanceIndex = 999999;
	char *seq, *str;

	// Store length of patterns to avoid recomputation
	patternLength = (int *) calloc(numPatterns, sizeof(int));
	for (i = 0; i < numPatterns; i++) {
		if (patternP[i].seq != NULL) {
			seq = (char*) patternP[i].seq;
			str = (char*) patternP[i].structure;
		} else {
			seq = (char*) patternP[i].seqR;
			str = (char*) patternP[i].structureR;
		}
		patternLength[i] = strlen(seq);

		if (patternP[i].maxstemlength > 0) {
			count = 0;
			for (j = 0; j < patternLength[i]; j++)
				if (str[j] == '(')
					count++;
			patternLength[i] += patternP[i].maxstemlength - count;
		}
		patternLength[i] += patternP[i].maxleftloopextent + patternP[i].maxrightloopextent;
	}

	int iNumInstances = 0;
	int *longestLength, *endPosInstance;
	patternArtificialInt = (int **) calloc(numPatterns, sizeof(int*));

	if (multiPattern->usesStartPosParam) {

		for (i = 0; i < numPatterns; i++) {
			patternArtificialInt[i] = (int *) calloc(2, sizeof(int));

			patternArtificialInt[i][0] = patternP[i].startpos;
			patternArtificialInt[i][1] = patternArtificialInt[i][0] + patternLength[i] - 1;
		}

	} else {

		if (multiPattern->usesInstanceParam) {

			for (i = 0; i < numPatterns; i++) {
				if(iNumInstances < patternP[i].instance) //find last instance
					iNumInstances = patternP[i].instance;
				if (lowestInstanceIndex > patternP[i].instance) { //find lowest instance index
					lowestInstanceIndex = patternP[i].instance;
				}
			}
			if (lowestInstanceIndex == 0) iNumInstances++;

			//find the longest length of the patterns of an instance
			longestLength = (int *) calloc(iNumInstances, sizeof(int));
			for (i = 0; i < iNumInstances; i++) {
				longestLength[i] = 0;
			}
			for (i = 0; i < numPatterns; i++) {
				if(longestLength[patternP[i].instance] < patternLength[i])
					longestLength[patternP[i].instance] = patternLength[i];
			}

			//determine all artificial ending positions of the chained instance
			endPosInstance = (int *) calloc(iNumInstances, sizeof(int*));
			endPosInstance[0] = longestLength[0] - 1;
			for (i = 1; i < iNumInstances; ++i) {
				endPosInstance[i] = endPosInstance[i - 1] + longestLength[i];
			}

			for (i = 0; i < numPatterns; i++) {
				patternArtificialInt[i] = (int *) calloc(2, sizeof(int));

				if(patternP[i].instance == lowestInstanceIndex) {
					patternArtificialInt[i][0] = 0;
					patternArtificialInt[i][1] = endPosInstance[0];
				} else {
					patternArtificialInt[i][0] = endPosInstance[patternP[i].instance - 1] + 1;
					patternArtificialInt[i][1] = patternArtificialInt[i][0] + longestLength[patternP[i].instance] - 1;
				}
			}

		} else {

			for (i = 0; i < numPatterns; i++) {
				patternArtificialInt[i] = (int *) calloc(2, sizeof(int));

				if(i == 0) {
					patternArtificialInt[i][0] = 0;
					patternArtificialInt[i][1] = patternLength[i] - 1;
				} else {
					patternArtificialInt[i][0] = patternArtificialInt[i - 1][1] + 1;
					patternArtificialInt[i][1] = patternArtificialInt[i][0] + patternLength[i] - 1;
				}
			}

		}

	}

	free(patternLength);
	return patternArtificialInt;
}

void chainAll(void *vrb, void *vrbCR, unsigned int *uiNumMatchesSeq, Chainparam *chainParam,
		MultiPattern *multiPattern, AffixArray *affixArray, bool bFreeMatches, bool bForwardSeq) {
	Pattern *patternP = multiPattern->pattern;
	int numPatterns = multiPattern->numPatterns;

	chainParamGlob  = chainParam;
	bForwardSeqGlob = bForwardSeq;
	seqEndPosGlob   = affixArray->multiSeq->seqEndPos;

	RBLIST *rblist;
	struct rbtree *rb = vrb;
	Match *match;
	int **patternArtificialInt, i, offset;

	rbCR = vrbCR;

	matchesRB = vrb;

	bool hasErr;
	GtChain2Dim *chain = NULL;
	GtChain2Dimmode *chainmode = NULL;

	const unsigned int presortdim = 1U;

	if ((rblist=rbopenlist(rb))==NULL) {
		fprintf(stderr, "Insufficient memory from rbopenlist().\n");
		exit(1);
	}

	patternArtificialInt = compute2DPositions(multiPattern);

	GtChain2Dimmatchtable *matchtable = NULL;
	unsigned int uiCount = 0;
	unsigned int uiPrevSeqId = UINT_MAX; // purpose: error checking
	while((match=(Match *)rbreadlist(rblist))) {

		if (uiCount == 0) {
			if (uiPrevSeqId == match->seqId) {
				fprintf(stderr, "The number of matches stored in uiNumMatchesSeq are incorrect.\n");
				exit(1);
			}
			matchtable = gt_chain_matchtable_new(uiNumMatchesSeq[match->seqId]);

			if (match->structureString == NULL)
				bShowStructureGlob = false;
			else
				bShowStructureGlob = true;
		}

		GtChain2Dimmatchvalues fragment;

		fragment.startpos[0] = patternArtificialInt[(match)->patternId][0];
		fragment.endpos[0]   = patternArtificialInt[(match)->patternId][1];
		fragment.weight      = patternP[(match)->patternId].weight == 0 ? 1 : patternP[(match)->patternId].weight;

		if (bForwardSeq) {
			if(match->seqId > 0)
				offset = seqEndPosGlob[match->seqId - 1] + 1;
			else
				offset = 0;
			fragment.startpos[1] = match->pos - offset;
			fragment.endpos[1]   = match->endPos - offset;
		} else {
			fragment.startpos[1] = seqEndPosGlob[match->seqId] - 1 - match->endPos;
			fragment.endpos[1]   = seqEndPosGlob[match->seqId] - 1 - match->pos;
		}

		//printf("%d %d %d %d\n", fragment.startpos[0], fragment.endpos[0], fragment.startpos[1], fragment.endpos[1]);

		gt_chain_matchtable_add(matchtable, &fragment);

		uiCount++;
		if (uiNumMatchesSeq[match->seqId] == uiCount) {
			hasErr = false;
			gt_chain_fillthegapvalues(matchtable);
			gt_chain_applyweight(chainParam->weightfactor, matchtable);
			gt_chain_possiblysortmatches(NULL, matchtable, presortdim);
			chain = gt_chain_chain_new();
			chainmode = gt_chain_chainmode_new(chainParam->maxgap,
							!chainParam->islocal,
							chainParam->globalparm,
							chainParam->islocal,
							chainParam->localparm,
							NULL);

			if (chainmode == NULL) {
				hasErr = true;
				break;
			}
			if (!hasErr) {
				Counter counter;
				counter.chaincounter = 0;

				matchSeqId = match->seqId;

				gt_chain_fastchaining(chainmode,
						chain,
						matchtable,
						true,
						presortdim,
						true,
						chaining2gooutput,
						&counter,
						NULL);
			}

			gt_chain_chain_delete(chain);
			gt_chain_chainmode_delete(chainmode);
			gt_chain_matchtable_delete(matchtable);

			uiCount = 0;
		}
		uiPrevSeqId = match->seqId;
		if (bFreeMatches) {
			rbdelete((void *)&match, rb);
			free(match);
		}
	}

	rbcloselist(rblist);

	for (i = 0; i < numPatterns; i++) {
		free(patternArtificialInt[i]);
	}
	free(patternArtificialInt);
}

