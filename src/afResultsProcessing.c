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
#include <limits.h>
#include <string.h>

#include "af.h"
#include "redblack/redblack.h"

int compareRB(const void *pa, const void *pb, const void *config) {
	Match *ma = (Match *)pa;
	Match *mb = (Match *)pb;

	if (ma->seqId < mb->seqId) return -1;
	if (ma->seqId > mb->seqId) return 1;

	if (ma->pos < mb->pos) return -1;
	if (ma->pos > mb->pos) return 1;

	if (ma->endPos < mb->endPos) return -1;
	if (ma->endPos > mb->endPos) return 1;

	if (ma->patternId == -1)
		return 0;

	if (ma->patternId < mb->patternId) return -1;
	if (ma->patternId > mb->patternId) return 1;

	return 0;
}

int compareRBCR(const void *pa, const void *pb, const void *config) {
	Chain *ma = (Chain *)pa;
	Chain *mb = (Chain *)pb;
	int i=0, longestChain;

	if(ma->chainScore > mb->chainScore) return -1;
	if(ma->chainScore < mb->chainScore) return 1;

	if(ma->seqId < mb->seqId) return -1;
	if(ma->seqId > mb->seqId) return 1;

	if(ma->forwardSeq > mb->forwardSeq) return -1;
	if(ma->forwardSeq < mb->forwardSeq) return 1;

	longestChain = ma->chainLength > mb->chainLength ? ma->chainLength : mb->chainLength;
	while (i < longestChain) {
		if (ma->fragmentInfo[i].startpos[1] > mb->fragmentInfo[i].startpos[1]) return -1;
		if (ma->fragmentInfo[i].startpos[1] < mb->fragmentInfo[i].startpos[1]) return 1;
		i++;
	}

	return 0;
}

int removeNonsignificantMatches(void *vrb, unsigned int *uiNumMatchesSeq, int n, int numPatterns) {
	RBLIST *rblist;
	struct rbtree *rb = vrb;
	Match *match, **matchRef;
	int i, k, numMatches, seqNumPrev, startingPos=0, countDifPatterns, countDeleted=0, temp;
	bool *bPatternOccurs;

	matchRef       = (Match **) calloc(BUFFER1, sizeof(Match*));
	bPatternOccurs = (bool *) calloc(numPatterns, sizeof(bool));

	if((rblist=rbopenlist(rb))==NULL) {
		fprintf(stderr, "Insufficient memory from rbopenlist().\n");
		exit(1);
	}

	numMatches = 0;
	while((match=(Match *)rbreadlist(rblist))) {
		if((numMatches + 1) % BUFFER1 == 0 && (matchRef = realloc(matchRef, (numMatches + BUFFER1 + 1) * sizeof(Match*))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"matchRef\". - %s %d.\n", __FILE__, __LINE__);
			exit(1);
		}
		matchRef[numMatches++] = match;
	}

	seqNumPrev = INT_MAX;
	for(i = 0; i <= numMatches; i++) {
		if(i == numMatches || matchRef[i]->seqId != seqNumPrev) {
			countDifPatterns = 0;
			for(k = 0; k < numPatterns; k++){
				countDifPatterns += bPatternOccurs[k];
			}
			if(countDifPatterns < n) {
				temp = i - startingPos;
				for(k = 0; k < temp; k++) {
					uiNumMatchesSeq[matchRef[startingPos + k]->seqId]--;
					rbdelete((void *)matchRef[startingPos + k], rb);
					free(matchRef[startingPos + k]);
					countDeleted++;
				}
			}

			startingPos = i; // starting position of the matches of a certain sequence in the copied array of matches
			for(k = 0; k < numPatterns; k++){
				bPatternOccurs[k] = 0;
			}
		}

		if(i < numMatches) {
			seqNumPrev = matchRef[i]->seqId;
			bPatternOccurs[matchRef[i]->patternId] = 1;
		}
	}

	free(matchRef);
	free(bPatternOccurs);
	rbcloselist(rblist);

	return countDeleted;
}

void saveSearchResultsRB(void *vrbF, void *vrbR, char *fileName, bool bBED, bool bPrintSeqDesc,
		Pattern *patternP, int numPatterns, AffixArray *affixArray, bool bFreeMatches) {
	RBLIST *rblist;
	struct rbtree *rb;
	Match *match;
	FILE *fp=NULL;
	int seqNumPrev, i, j, k, offset;

	if (vrbF == NULL && vrbR == NULL)
		return;

	if (fileName == NULL) {
		fp = stdout;
		fprintf(fp, "\n");
	} else if((fp = fopen(fileName, "w")) == NULL) {
	   fprintf(stderr,"Error saving file %s.\n", fileName);
	   exit(1);
	}

	if (!bBED) {
		if (bPrintSeqDesc)
			fprintf(fp, "%c[matched substring]	[structure]	[matching pos.]	[pattern id]	[weight]	[strand]\n", LINESYMBOL);
		else
			fprintf(fp, "%c[matched substring]	[structure]	[seq. id]	[matching pos.]	[pattern id]	[weight]	[strand]\n", LINESYMBOL);
	} else {
		fprintf(fp, "#[sequence]	[start]	[end]	[pattern id]	[weight]	[strand]\n");
	}

	for (k = 0; k < 2; k++) {
		if (k == 0 && vrbF != NULL)
			rb = vrbF;
		else if (k == 1 && vrbR != NULL)
			rb = vrbR;
		else
			continue;

		if ((rblist = rbopenlist(rb))==NULL) {
			fprintf(stderr, "Insufficient memory from rbopenlist() - %s %d.\n", __FILE__, __LINE__);
			exit(1);
		}

		seqNumPrev = INT_MAX;
		while((match = (Match *)rbreadlist(rblist))) {
			if (seqNumPrev != match->seqId && bPrintSeqDesc && !bBED)
				fprintf(fp, ">%s\n", affixArray->multiSeq->seqDescription[match->seqId]);
			seqNumPrev = match->seqId;

			if(match->seqId > 0)
				offset = affixArray->multiSeq->seqEndPos[match->seqId - 1] + 1;
			else
				offset = 0;

			if (!bBED) {
				i = match->endPos + 1;
				for (j = match->pos; j < i; j++) { // For printing the matched pattern, the absolute position is required
					fprintf(fp, "%c", affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[j] - 1]);
				}
				fprintf(fp, "	%s	", match->structureString);
			}

			if (!bPrintSeqDesc)
				fprintf(fp, "%d", match->seqId);
			else if (bBED)
				fprintf(fp, "%s", affixArray->multiSeq->seqDescription[match->seqId]);

			fprintf(fp, "	%d", match->pos - offset); // Print position relative to sequence and not to concatenated sequences

			if (bBED) {
				fprintf(fp, "	%d", match->endPos - offset); // Print position relative to sequence and not to concatenated sequences
			}

			fprintf(fp, "	%d", match->patternId);
			fprintf(fp, "	%d", patternP[match->patternId].weight);

			if (bBED) {
				fprintf(fp, "	%c\n", k == 0 ? '+' : '-');
			} else {
				fprintf(fp, "	%c\n", k == 0 ? 'f' : 'r');
			}

			if (bFreeMatches)
				free(match);
		}

		rbcloselist(rblist);
	}


	fprintf(fp, "\r");
	if (fp != stdout)
		fclose(fp);
}

void reportChainingResults(void *vrb, Chainparam *chainparam, AffixArray *affixArray) {
	struct rbtree *rb = vrb;
	RBLIST *rblist;
	Chain *chainingResult;
	int i, j, iDescLegth, absoluteStartMatchPos, absoluteEndMatchPos, offset;
	char *reportFileName = chainparam->reportfile;
	bool bShowChain = chainparam->show;
	bool bChainNotReversed;

	FILE *fp=NULL;
	if (reportFileName == NULL) {
		fp = stdout;
	} else if((fp = fopen(reportFileName, "w")) == NULL) {
		fprintf(stderr,"Error saving file %s.\n", reportFileName);
		exit(1);
	}

	if((rblist=rbopenlist(rb)) == NULL) {
		fprintf(stderr, "Insufficient memory from rbopenlist().\n");
		exit(1);
	}

	fprintf(fp, "%c[sequence]                                	[chain score]	[chain length]	[strand]\n", LINESYMBOL);
	while((chainingResult = (Chain *)rbreadlist(rblist))) {
		bChainNotReversed = !chainingResult->bStoredInReverseOrder;
		iDescLegth = strlen(affixArray->multiSeq->seqDescription[chainingResult->seqId]);
		fprintf(fp, ">");
		for(i = 0; i < 40; i++) {
			if(i < iDescLegth)
				fprintf(fp, "%c", affixArray->multiSeq->seqDescription[chainingResult->seqId][i]);
			else
				fprintf(fp, " ");
		}
		if(i < iDescLegth)
			fprintf(fp, "...");
		else
			fprintf(fp, "   ");
		fprintf(fp, "	");
		fprintf(fp, "%d		", chainingResult->chainScore);
		fprintf(fp, "%d 		", chainingResult->chainLength);
		fprintf(fp, "%c\n", chainingResult->forwardSeq ? 'f' : 'r');

		if (bShowChain && chainingResult->fragmentInfo != NULL) {
			if (bChainNotReversed)
				for (i = 0; i < chainingResult->chainLength; i++) {
					fprintf(fp, "%d %d ", chainingResult->fragmentInfo[i].startpos[0], chainingResult->fragmentInfo[i].endpos[0]);
					fprintf(fp, "%d %d ", chainingResult->fragmentInfo[i].startpos[1], chainingResult->fragmentInfo[i].endpos[1]);
					fprintf(fp, "%d\n", chainingResult->fragmentInfo[i].weight);
				}
			else
				for (i = chainingResult->chainLength; i > 0; i--) {
					fprintf(fp, "%d %d ", chainingResult->fragmentInfo[i - 1].startpos[0], chainingResult->fragmentInfo[i - 1].endpos[0]);
					fprintf(fp, "%d %d ", chainingResult->fragmentInfo[i - 1].startpos[1], chainingResult->fragmentInfo[i - 1].endpos[1]);
					fprintf(fp, "%d\n", chainingResult->fragmentInfo[i - 1].weight);
				}

			if(chainingResult->seqId > 0)
				offset = affixArray->multiSeq->seqEndPos[chainingResult->seqId - 1] + 1;
			else
				offset = 0;

			if (bChainNotReversed)
				for (i = 0; i < chainingResult->chainLength; i++) {
					absoluteStartMatchPos = offset + chainingResult->fragmentInfo[i].startpos[1];
					absoluteEndMatchPos   = offset + chainingResult->fragmentInfo[i].endpos[1];

					for (j = absoluteStartMatchPos; j <= absoluteEndMatchPos; j++) { // For printing the matched pattern, the absolute position is required
						fprintf(fp, "%c", affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[j] - 1]);
					}
					fprintf(fp, " ");
				}
			else
				for (i = chainingResult->chainLength; i > 0; i--) {
					absoluteStartMatchPos = offset + chainingResult->fragmentInfo[i - 1].startpos[1];
					absoluteEndMatchPos   = offset + chainingResult->fragmentInfo[i - 1].endpos[1];

					for (j = absoluteStartMatchPos; j <= absoluteEndMatchPos; j++) { // For printing the matched pattern, the absolute position is required
						fprintf(fp, "%c", affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[j] - 1]);
					}
					fprintf(fp, " ");
				}

			fprintf(fp, "\n");

			if (chainingResult->fragmentInfo[0].structureString != NULL) {
				if (bChainNotReversed)
					for (i = 0; i < chainingResult->chainLength; i++) {
						fprintf(fp, "%s ", chainingResult->fragmentInfo[i].structureString);
					}
				else
					for (i = chainingResult->chainLength; i > 0; i--) {
						fprintf(fp, "%s ", chainingResult->fragmentInfo[i - 1].structureString);
					}
				fprintf(fp, "\n");
			}

			free(chainingResult->fragmentInfo);
		}
		free(chainingResult);
	}

	rbcloselist(rblist);

	if (fp != stdout)
		fclose(fp);
}

