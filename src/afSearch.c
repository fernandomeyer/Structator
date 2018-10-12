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
#include <string.h>
#include <stdio.h>
#include <sys/time.h>

#include "af.h"

#ifndef DEBUG_ST_SEARCH
#define DEBUG_ST_SEARCH 0
#endif

#define rpref(index, d) ((rpref[index] + d) >= (affixArray->length - 1) ? $ : seq[-rpref[index] - (d)])

void printNode(VirtualNode *node, AffixArray *affixArray) {
	int x, y;

	printf("context=%d	i=%d	j=%d	d=%d	type=%d	",node->context, node->i, node->j, node->d, node->type);

	x = node->i;
	if(node->type == 0) {
		for(y = node->context; y < node->d; y++)
			printf("%c",affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[affixArray->esa->xarray[x] + y] - 1]);
	} else if(node->type == 1) {
		for(y = - node->d; y < -node->context; y++)
			printf("%c",affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[(affixArray->length - 1) - affixArray->erpa->xarray[x] + y] - 1]);
	} else {
		for(x = node->i; x <= node->j; x++)
			printf("%c",affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[x] - 1]);
	}
	printf("\n");
}

/* Returns the sequence number to which a suffix belongs. k is the
 * position in the concatenated string of sequences where the suffix
 * begins. The number returned is the order of appearance of the
 * sequence in the fasta file, starting from 0.
 * Function uses binary search. */
int getSeqNumber(int *range, int length, int k) {
	int low = 0;
	int high = length - 1;
	int mid;
	int vMidLeft;
	int vMidRight;

	while (low <= high) {
		mid = (low + high) >> 1;
		vMidLeft = (mid == 0) ? 0 : (range[mid-1] + 1);
		vMidRight = range[mid];

		if (vMidRight < k)
			low = mid + 1;
		else if (vMidLeft > k)
			high = mid - 1;
		else
			return mid;
	}
	return -1;
}

/* Exact pattern search using lcp and skp tables. Returns an array with
 * the absolut position of each match. The first position of the array
 * contains the number of matches.
 * This is an adapted version of the algorithm for matching PSSMs by
 * Beckstette, Homann, Giegerich, and Kurtz, 2006 */
int* findPattern(unsigned char *pattern, int m, AffixArray *affixArray) {
	int depth, i, d, v, *pos, n, numMatches;
	unsigned char *seq, seqChar, patternChar;

	n          = affixArray->length;
	pos        = (int *) calloc(BUFFER1, sizeof(int));
	numMatches = 1;

	depth = i = 0;

	while(i < n) {
		if(n - m < affixArray->esa->xarray[i]) {
			while(n - m < affixArray->esa->xarray[i] && i < n) {
				i++;
				depth =  depth < lcpvalue(i) ? depth : lcpvalue(i);
			}
			if(i >= n) {
				pos[0] = numMatches - 1;
				return pos;
			}
		}

		d   = depth - 1;
		seq = affixArray->multiSeq->convSequences + affixArray->esa->xarray[i];

		do {
			d++;
			seqChar     = *(seq + d);
        	patternChar = *(pattern + d);
        	v = seqChar - patternChar;
		} while((d < m - 1) && (v == 0 || affixArray->alphabet->isWildCard[seqChar - 1] || patternChar == '*'));

		if((d == m - 1) && (v == 0 || affixArray->alphabet->isWildCard[seqChar - 1] || patternChar == '*')) {
			//match at position affixArray->esa->xarray[i]
			pos[numMatches++] = affixArray->esa->xarray[i];
	    	if(numMatches % BUFFER1 == 0)
		    	if((pos = realloc(pos, (numMatches + BUFFER1) * sizeof(int))) == NULL) {
					fprintf(stderr,"Memory allocation failed for \"pos\".\n");
					exit(1);
				}
			while(i < n) {
				i++;
				if(lcpvalue(i) >= m) {
					//match at position affixArray->esa->xarray[i]
					pos[numMatches++] = affixArray->esa->xarray[i];
			    	if(numMatches % BUFFER1 == 0)
				    	if((pos = realloc(pos, (numMatches + BUFFER1) * sizeof(int))) == NULL) {
							fprintf(stderr,"Memory allocation failed for \"pos\".\n");
							exit(1);
						}
				}
				else {
					//since all matches are always in consecutive positions, return
					/* THE STATEMENT ABOVE IS NOT TRUE IF THERE ARE WILDCARDS, SO BREAK ONLY*/
					pos[0] = numMatches - 1;
					return pos;
					break;
				}
			}
		} else {
			/* skipchain */
			if(i < n) {
				i++;
				while(i <= n && lcpvalue(i) > d)
					i = affixArray->esa->xskp[i];
			} else
				i = n;
		}
		depth = lcpvalue(i);
	}

	pos[0] = numMatches - 1;
	return pos;
}

int* findPatternSlow(unsigned char *pattern, int pLength, AffixArray *affixArray) {
	int i, j, v, *pos, k=1;
	unsigned char *seq, seqChar;
	unsigned char *patternCopy, patternChar;
	unsigned char *seq2bSearched;

	seq2bSearched = affixArray->multiSeq->convSequences;

	patternCopy = pattern;

	pos = (int *) calloc(BUFFER1, sizeof(int));

	for(i = 0; i < affixArray->length; i++) {
		seq = seq2bSearched + i;
		pattern = patternCopy;
		j = v = 0;

		do {
			seqChar     = *seq++;
	        patternChar = *pattern++;
	        v = seqChar - patternChar;
	    } while(v == 0 && seqChar != $ && ++j < pLength);

	    if(v == 0) { // pattern found at position i
	    	pos[k++] = i;
	    	if(k % BUFFER1 == 0)
		    	if((pos = realloc(pos, (k + BUFFER1) * sizeof(int))) == NULL) {
					fprintf(stderr,"Memory allocation failed for \"pos\".\n");
					exit(1);
				}
	    }
	}

	pos[0] = k - 1; // #times pattern was found

	return pos;
}

/* Receives an array of absolut sequence positions and returns an
 * array of their relative positions w.r.t. the offset (i.e. #positions
 * shifted due to the symbol $) */
int* getRelativePos(int *pos, int length, AffixArray *affixArray) {
	int i, seqNum, *relPos;

	relPos = (int *) calloc(length, sizeof(int));

	for(i = 0; i < length; i++) {
		seqNum = getSeqNumber(affixArray->multiSeq->seqEndPos, affixArray->multiSeq->numSeqs, pos[i]);
		relPos[i] = pos[i] - seqNum;
	}

	return relPos;
}

int findAffixLink(unsigned char *sequences, unsigned char *pattern, int m, int n, EArray *earray, Alphabet *alphabet) {
	int depth, i, d, v, pos=-1;
	unsigned char *seq, seqChar, patternChar;

	depth = i = 0;

	while(i < n) {
		if(n - m < earray->xarray[i]) {
			while(n - m < earray->xarray[i] && i < n) {
				i++;
				depth =  depth < xlcpvalue(i) ? depth : xlcpvalue(i);
			}
			if(i >= n)
				return pos;
		}

		d   = depth - 1;
		seq = sequences + earray->xarray[i];

		do {
			d++;
			seqChar     = *(seq + d);
			patternChar = *(pattern + d);
			if (alphabet->isWildCard[seqChar - 1] || alphabet->isWildCard[patternChar - 1]) {
				v = -1;
				break;
			}
        	v = seqChar - patternChar;
		} while(d < m - 1 && v == 0);

		if(d == m - 1 && v == 0) {
			//match at position earray->xarray[i]
			pos = i;
			return pos;
		} else {
			/* skipchain */
			if(i < n) {
				i++;
				while(i <= n && xlcpvalue(i) > d)
					i = earray->xskp[i];
			} else
				i = n;
		}
		depth = xlcpvalue(i);
	}

	return pos;
}

int findAffixLink2(unsigned char *sequences, unsigned char *pattern, int m, int n, EArray *earray, Alphabet *alphabet) {
	int depth, i, d, v, pos=-1;
	unsigned char *seq, seqChar, patternChar;

	depth = i = 0;
	pattern = pattern + m - 1;

	while(i < n) {
		if(n - m < earray->xarray[i]) {
			while(n - m < earray->xarray[i] && i < n) {
				i++;
				depth =  depth < xlcpvalue(i) ? depth : xlcpvalue(i);
			}
			if(i >= n)
				return pos;
		}

		d   = depth - 1;

		seq = sequences + (n - 2) - earray->xarray[i];

		do {
			d++;
        	seqChar     = *(seq - d);
        	patternChar = *(pattern - d);
        	if (alphabet->isWildCard[seqChar - 1] || alphabet->isWildCard[patternChar - 1]) {
				v = -1;
				break;
			}
        	v = seqChar - patternChar;
		} while(d < m - 1 && v == 0);

		if(d == m - 1 && v == 0) {
			//match at position earray->xarray[i]
			pos = i;
			return pos;
		} else {
			/* skipchain */
			if(i < n) {
				i++;
				while(i <= n && xlcpvalue(i) > d)
					i = earray->xskp[i];
			} else
				i = n;
		}
		depth = xlcpvalue(i);
	}

	return pos;
}

/*Find pattern using binary search with character comparisons in left-to-right direction. */
int* findPatternBS(unsigned char *pattern, int pLength, int initLeft, int initRight, int d, AffixArray *affixArray) {
	if(affixArray->esa->xarray == NULL) {
		fprintf(stderr,"Suf table was not loaded.\n");
		exit(1);
	}

	unsigned char *convSequences = affixArray->multiSeq->convSequences;
	int left, l, intervalLeft;
	int right, r, r_bak, intervalRight;
	int mid, h=0;
	int dpattern;

	int *result = (int *) calloc(2, sizeof(int));
	int *suf = affixArray->esa->xarray;

	intervalLeft = intervalRight = -1;

	l = r = dpattern = 0;
	while (l < pLength && pattern[l] == convSequences[suf[initLeft] + d + l]) {
		l++;
	}
	while (r < pLength && pattern[r] == convSequences[suf[initRight] + d + r]) {
		r++;
	}
	r_bak = r;

	left = initLeft;
	right = initRight;

	if (l == pLength) {
		intervalLeft = initLeft;
	} else if (pattern[0] >  convSequences[suf[initRight] + d]) {
		intervalLeft = initRight;
	} else {
		while (right > left + 1) {
			mid = ((unsigned int)left + (unsigned int)right) >> 1;

			h = l < r ? l : r;
			while (h < pLength && pattern[h] == convSequences[suf[mid] + d + h]) {
				h++;
			}

			if (h == pLength || pattern[h] <= convSequences[suf[mid] + d + h]) {
				right = mid;
				r = h;
			} else {
				left = mid;
				l = h;
			}

		}

		if (r < pLength) {
			result[0] = -1; //The searched pattern does not exist
			return result;
		}

		intervalLeft = right;
		l = r;
	}

	left = intervalLeft;
	right = initRight;

	r = r_bak;

	if (r == pLength) {
		intervalRight = initRight;
	} else {
		while (right > left + 1) {
			mid = ((unsigned int)left + (unsigned int)right) >> 1;

			h = l < r ? l : r;
			while (h < pLength && pattern[h] == convSequences[suf[mid] + d + h]) {
				h++;
			}

			if (h == pLength || pattern[h] >= convSequences[suf[mid] + d + h]) {
				left = mid;
				l = h;
			} else {
				right = mid;
				r = h;
			}

		}

		if (l < pLength) {
			result[0] = -1; //The searched pattern does not exist
			return result;
		}

		intervalRight = left;
	}

	result[0] = intervalLeft;
	result[1] = intervalRight;

	return result;

}

/* Find pattern using binary search with character comparisons in right-to-left direction.
 * Pointer *pattern points to the last character of the pattern. */
int* findPatternBS2(unsigned char *pattern, int pLength, int initLeft, int initRight, int d, AffixArray *affixArray) {
	unsigned char *seq = affixArray->multiSeq->convSequences + (affixArray->length - 2);
	int left, l, intervalLeft;
	int right, r, r_bak, intervalRight;
	int mid, h=0;
	int dpattern;

	int *result = (int *) calloc(2, sizeof(int));
	int *rpref = affixArray->erpa->xarray;

	intervalLeft = intervalRight = -1;

	l = r = dpattern = 0;
	while (l < pLength && pattern[-l] == rpref(initLeft, d + l)) {
		l++;
	}
	while (r < pLength && pattern[-r] == rpref(initRight, d + r)) {
		r++;
	}
	r_bak = r;

	left = initLeft;
	right = initRight;

	if (l == pLength) {
		intervalLeft = initLeft;
	} else if (pattern[0] >  rpref(initRight, d)) {
		intervalLeft = initRight;
	} else {
		while (right > left + 1) {
			mid = ((unsigned int)left + (unsigned int)right) >> 1;

			h = l < r ? l : r;
			while (h < pLength && pattern[-h] == rpref(mid, d + h)) {
				h++;
			}

			if (h == pLength || pattern[-h] <= rpref(mid, d + h)) {
				right = mid;
				r = h;
			} else {
				left = mid;
				l = h;
			}

		}

		if (r < pLength) {
			result[0] = -1; //The searched pattern does not exist
			return result;
		}

		intervalLeft = right;
		l = r;
	}

	left = intervalLeft;
	right = initRight;

	r = r_bak;

	if (r == pLength) {
		intervalRight = initRight;
	} else {
		while (right > left + 1) {
			mid = ((unsigned int)left + (unsigned int)right) >> 1;

			h = l < r ? l : r;
			while (h < pLength && pattern[-h] == rpref(mid, d + h)) {
				h++;
			}

			if (h == pLength || pattern[-h] >= rpref(mid, d + h)) {
				left = mid;
				l = h;
			} else {
				right = mid;
				r = h;
			}

		}

		if (l < pLength) {
			result[0] = -1; //The searched pattern does not exist
			return result;
		}

		intervalRight = left;
	}

	result[0] = intervalLeft;
	result[1] = intervalRight;

	return result;

}

int getLcpIntervalDepth(int i, int j, int initDepth, AffixArray *affixArray) {
	int depth = initDepth;

	unsigned char *lowerBound = affixArray->multiSeq->convSequences + affixArray->esa->xarray[i] + initDepth;
	unsigned char *upperBound = affixArray->multiSeq->convSequences + affixArray->esa->xarray[j] + initDepth;

	while(*lowerBound == *upperBound && *lowerBound != $ && *upperBound != $
			&& !affixArray->alphabet->isWildCard[*lowerBound -1] && !affixArray->alphabet->isWildCard[*upperBound -1]) {
		lowerBound++;
		upperBound++;
		depth++;
	}

	return depth;
}

int getRLcpIntervalDepth(int i, int j, int initDepth, AffixArray *affixArray) {
	int depth = initDepth;

	unsigned char *lowerBound = affixArray->multiSeq->convSequences + (affixArray->length - 2) - affixArray->erpa->xarray[i] - initDepth;
	unsigned char *upperBound = affixArray->multiSeq->convSequences + (affixArray->length - 2) - affixArray->erpa->xarray[j] - initDepth;

	while((affixArray->erpa->xarray[i] + depth < affixArray->length - 1)
			&& (affixArray->erpa->xarray[j] + depth < affixArray->length - 1) && *lowerBound == *upperBound
			&& *lowerBound != $ && *upperBound != $
			&& !affixArray->alphabet->isWildCard[*lowerBound -1] && !affixArray->alphabet->isWildCard[*upperBound -1]) {
		depth++;
		lowerBound--;
		upperBound--;
	}

	return depth;
}

#define homePos(i, j, data) (lcpvalueX(data, i) < lcpvalueX(data, j + 1) ? j : i)

bool extendRight(unsigned char *charExtension, int plength, VirtualNode *node, AffixArray *affixArray) {
	#if DEBUG_ST_SEARCH
		int x;
		for(x = 0; x < plength; x++)
			printf("'%c'", affixArray->alphabet->classRepresentative[charExtension[x] - 1]);
		if(plength > 0) printf("\n");
		printNode(node, affixArray);
	#endif

	int home, numSubstrings;
	int *result;
	unsigned char *seq=NULL;

	//suf node before the actual extension. NOTE THAT THE ORDER OF THE ASSIGNMENTS MATTERS
	if(node->type == 1) { //If the last node was rpref...
		if(node->context > 0)
			seq = affixArray->multiSeq->convSequences + (affixArray->length - 2) - affixArray->erpa->xarray[node->i] - (node->context - 1);
		while(node->context > 0 && plength > 0) {
			#if DEBUG_ST_SEARCH
				printf("compare: %d %d, %d\n",*charExtension, *seq, affixArray->erpa->xarray[node->i]);
			#endif
			if(*charExtension++ != *seq++)
				return 1; //pattern not found
			plength--;
			node->context--;
		}

		if(node->context == 0  && plength > 0) { //if there is more to be searched...

			/* Modified by Benjamin Albrecht, Jan. 2012 */
			//COMPUTE AFFIX LINK INSTANTLY ---------------------------------
			int * affixLink = NULL;
			unsigned char * prefix;
			int prefixSize;

			if (affixArray->erpa->affixLink == NULL) {
				// compute prefix
				prefixSize = getRLcpIntervalDepth(node->i, node->j, node->d,
						affixArray);
				if ((prefix = (unsigned char *) malloc(prefixSize
						* sizeof(unsigned char))) == NULL)
					return 1;
				int y, k = 0;
				for (y = -prefixSize; y < -node->context; y++) {
					prefix[k]
							= affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[(affixArray->length
									- 1) - affixArray->erpa->xarray[node->i]
									+ y] - 1];
					k++;
				}
				// search prefix in esa...
				if (affixArray->esa->xlcpTree != NULL
						&& affixArray->esa->xlcpException != NULL) {
					//...via improved binary search
					affixLink = cmpAflkViaImprovedBS(0, affixArray->length - 1, -1, -1, affixArray, prefix, prefixSize, 1);
				} else {
					//...via a normal binary seach
					affixLink = cmpAflkViaBS(0, affixArray->length - 1, affixArray, prefix, prefixSize, 1);
				}
				free(prefix);
			}
			//--------------------------------------------------------------

			home           = homePos(node->i, node->j, affixArray->erpa);
			node->context  = getRLcpIntervalDepth(node->i, node->j, node->d, affixArray) - node->d;
			node->d       += node->context;

			numSubstrings  = node->j - node->i;

			/* Modified by Benjamin Albrecht, Jan. 2012 */
			if (affixArray->erpa->affixLink == NULL) {
				node->i    = affixLink[0];
				free(affixLink);
			} else
				node->i    = affixArray->erpa->affixLink[home];

			node->j        = node->i + numSubstrings;
		}
	} else if(node->type > 1) {
		if(node->j < (affixArray->length - 2))
			seq = affixArray->multiSeq->convSequences + node->j + 1;
		else
			return 1; //pattern can no longer be extended to the right
		while(node->j <= (affixArray->length - 2) && plength > 0) {
			if(*charExtension++ != *seq++)
				return 1; //pattern not found
			plength--;
			node->j++;
			node->d++;
		}
		node->type = 2;
		if(plength > 0) return 1; //pattern not found
	}

	if(plength > 0) {
		result = findPatternBS(charExtension,
								plength,
								node->i,
								node->j,
								node->d, //depth for search start
								affixArray);

		if(result[0] == -1) {
			free(result);
			return 1; //pattern not found
		}

		node->i  = result[0];
		node->j  = result[1];
		free(result);
		node->d += plength;
		if(node->i == node->j) { //there is a single occurrence of the pattern
			/* Note1: node->context and node->d will no longer matter for future extensions
			 * Note2: node->i and node->j will no longer be associated with the lexicographic ranking, but
			 * with the direct position in the sequence */
			node->type  = 2;
			node->i     = affixArray->esa->xarray[node->i] + node->context;
			node->j     = node->i - node->context + (node->d - 1);
		} else {
			node->type  = 0;
		}
	}
	#if DEBUG_ST_SEARCH
		printNode(node, affixArray);
	#endif

	return 0;
}

bool extendLeft(unsigned char *charExtension, int plength, VirtualNode *node, AffixArray *affixArray) {
	#if DEBUG_ST_SEARCH
		int x;
		for(x = 0; x < plength; x++) #endif
			printf("'%c'", affixArray->alphabet->classRepresentative[charExtension[x] - 1]);
		if(plength > 0) printf("\n");
		printNode(node, affixArray);
	#endif

	int home, numSubstrings;
	int *result;
	unsigned char *seq=NULL;

	int depth;
	unsigned char *lowerBound, *upperBound;

	//rpref node before the actual extension. NOTE THAT THE ORDER OF THE ASSIGNMENTS MATTERS
	if(node->type == 0) { //If the last node was suf...
		if(node->context > 0) {
			seq = affixArray->multiSeq->convSequences + affixArray->esa->xarray[node->i] + (node->context - 1);
			charExtension += plength - 1;
		}
		while(node->context > 0 && plength > 0) {
			if(*charExtension-- != *seq--)
				return 1; //pattern not found
			plength--;
			node->context--;
		}

		if(node->context == 0  && plength > 0) { //if there is more to be searched...
			home           = homePos(node->i, node->j, affixArray->esa);

			depth = node->d;
			lowerBound = affixArray->multiSeq->convSequences + affixArray->esa->xarray[node->i] + depth;
			upperBound = affixArray->multiSeq->convSequences + affixArray->esa->xarray[node->j] + depth;

			while(*lowerBound == *upperBound && *lowerBound != $ && *upperBound != $
					&& !affixArray->alphabet->isWildCard[*lowerBound -1] && !affixArray->alphabet->isWildCard[*upperBound -1]) {
				lowerBound++;
				upperBound++;
				depth++;
			}

			/* Modified by Benjamin Albrecht, Jan. 2012 */
			//COMPUTE AFFIX LINK INSTANTLY ---------------------------------
			int *affixLink = NULL;
			int prefixSize;
			unsigned char * prefix;

			if (affixArray->esa->affixLink == NULL) {
				// compute prefix
				prefixSize = node->d + (depth - node->d);
				if ((prefix = (unsigned char *) malloc(prefixSize * sizeof(unsigned char))) == NULL)
					return 1;
				int y, k = 0;
				for (y = node->context; y < prefixSize; y++) {
					prefix[k] = affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[affixArray->esa->xarray[node->i]+ y] - 1];
					k++;
				}
				// search prefix in erpa...
				if (affixArray->erpa->xlcpTree != NULL && affixArray->erpa->xlcpException != NULL) {
					//... via improved binary search
					affixLink = cmpAflkViaImprovedBS(0, affixArray->length - 1, -1, -1, affixArray, prefix,prefixSize, 0);
				} else {
					//... via a normal binary search
					affixLink = cmpAflkViaBS(0, affixArray->length - 1, affixArray, prefix, prefixSize, 0);
				}
				free(prefix);
			}
			//--------------------------------------------------------------

			node->context  = depth - node->d;

			node->d       += node->context;
			numSubstrings  = node->j - node->i;

			/* Modified by Benjamin Albrecht, Jan. 2012 */
			if (affixArray->esa->affixLink == NULL) {
				node->i    = affixLink[0];
				free(affixLink);
			} else
				node->i    = affixArray->esa->affixLink[home];

			node->j        = node->i + numSubstrings;
		}
	} else if(node->type > 1) {
		if(node->i > 0) {
			seq = affixArray->multiSeq->convSequences + node->i - 1;
			charExtension += plength - 1;
		} else
			return 1; //pattern can no longer be extended to the left
		while(node->i >= 0 && plength > 0) {
			if(*charExtension-- != *seq--)
				return 1; //pattern not found
			plength--;
			node->i--;
			node->d++;
		}
		node->type = 3;
		if(plength > 0) return 1; //pattern not found
	}

	if(plength > 0) {
		result = findPatternBS2(charExtension,
								plength,
								node->i,
								node->j,
								node->d, //depth for search start
								affixArray);

		if(result[0] == -1) {
			free(result);
			return 1; //pattern not found
		}

		node->i  = result[0];
		node->j  = result[1];
		free(result);
		node->d += plength;
		if(node->i == node->j) { //there is a single occurrence of the pattern
			/* Note1: node->context and node->d will no longer matter for future extensions
			 * Note2: node->i and node->j will no longer be associated with the lexicographic ranking, but
			 * with the direct position in the sequence */
			node->type  = 3;
			node->i     = (affixArray->length - 1) - affixArray->erpa->xarray[node->i] - node->d;
			node->j     = node->i + (node->d - node->context) - 1;
		} else {
			node->type  = 1;
		}
	}

	return 0;
}
