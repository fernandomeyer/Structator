/* Modified by Benjamin Albrecht, Jan. 2012*/

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
#include "divsufsort/config.h"
#include "divsufsort/divsufsort.h"
#include "divsufsort/lfs.h"

/* Modified by Benjamin Albrecht, Jan. 2012 */
void computeAffixArray(AffixArray *affixArray,
					bool argSuf, bool argLcp, bool argLcpTree, bool argSkp, bool argAfsuf,
					bool argRpref, bool argRlcp, bool argRlcpTree, bool argRskp, bool argAfrpref,
					bool argOScreen, char *fileName, bool argShowTimes) {

	/* Modified by Benjamin Albrecht, Jan. 2012 */
	/* Compute Enhanced Suffix Array */
	computeEArray(affixArray->esa,
				affixArray->alphabet,
				affixArray->multiSeq->convSequences,
				affixArray->multiSeq->numSeqs,
				affixArray->length,
				/*isEsa*/ 1,
				argSuf, argLcp, argLcpTree, argSkp,
				argShowTimes);

	if (fileName != NULL)
		if (saveEArray(affixArray->esa, affixArray->length, 1, fileName)){
			fprintf(stderr, "Error occurred while storing the suffix array");
			exit(1);
		}

	/* Modified by Benjamin Albrecht, Jan. 2012 */
	/* Free unnecessary structures */
	if(!argAfsuf && !argAfrpref && !argOScreen) {
		if (affixArray->esa->xarray != NULL) {
			free(affixArray->esa->xarray);
			affixArray->esa->xarray = NULL;
		}
		if(affixArray->esa->xlcp != NULL) {
			free(affixArray->esa->xlcp);
			affixArray->esa->xlcp = NULL;
		}
		if (affixArray->esa->xlcpException != NULL) {
			free(affixArray->esa->xlcpException);
			affixArray->esa->xlcpException = NULL;
		}
		if (affixArray->esa->xlcpTree != NULL) {
			free(affixArray->esa->xlcpTree);
			affixArray->esa->xlcpTree = NULL;
		}
		if (affixArray->esa->xlcpTreeException != NULL) {
			free(affixArray->esa->xlcpTreeException);
			affixArray->esa->xlcpTreeException = NULL;
		}
		if(affixArray->esa->xskp != NULL) {
			free(affixArray->esa->xskp);
			affixArray->esa->xskp = NULL;
		}
	}

	/* Compute Enhanced Reverse Prefix Array */
	computeEArray(affixArray->erpa,
				affixArray->alphabet,
				affixArray->multiSeq->convSequences,
				affixArray->multiSeq->numSeqs,
				affixArray->length,
				/*!isEsa*/ 0,
    			argRpref, argRlcp, argRlcpTree, argRskp,
    			argShowTimes);

	if (fileName != NULL)
		if (saveEArray(affixArray->erpa, affixArray->length, 0, fileName)) {
			fprintf(stderr, "Error occurred while storing the reverse prefix array");
			exit(1);
		}

	/* Modified by Benjamin Albrecht, Jan. 2012 */
	/* Free unnecessary structures */
	if(!argAfsuf && !argAfrpref && !argOScreen) {
		if (affixArray->erpa->xarray != NULL) {
			free(affixArray->erpa->xarray);
			affixArray->erpa->xarray = NULL;
		}
		if(affixArray->erpa->xlcp != NULL) {
			free(affixArray->erpa->xlcp);
			affixArray->erpa->xlcp = NULL;
		}
		if (affixArray->esa->xlcpException != NULL) {
			free(affixArray->esa->xlcpException);
			affixArray->esa->xlcpException = NULL;
		}
		if (affixArray->esa->xlcpTree != NULL) {
			free(affixArray->esa->xlcpTree);
			affixArray->esa->xlcpTree = NULL;
		}
		if (affixArray->esa->xlcpTreeException != NULL) {
			free(affixArray->esa->xlcpTreeException);
			affixArray->esa->xlcpTreeException = NULL;
		}
		if(affixArray->erpa->xskp != NULL) {
			free(affixArray->erpa->xskp);
			affixArray->erpa->xskp = NULL;
		}
	}

	if(argAfsuf) {
		if(affixArray->erpa->xlcp != NULL && affixArray->erpa->xskp != NULL)
			constructAflk(affixArray->esa, affixArray->erpa, true, (int (*)(void *, void *, int, int, void *, void *))(findAffixLink2), affixArray->multiSeq->convSequences, affixArray->length, affixArray->alphabet, argShowTimes);
		else
			constructAflkBS(affixArray->esa, affixArray->erpa, true, (int* (*)(void *, int, int, int, int, void *))(findPatternBS2), affixArray->multiSeq->convSequences, affixArray->length, affixArray, argShowTimes);

		if (fileName != NULL) {
			if (saveAflk(affixArray->esa, affixArray->length, 1, fileName)) {
				fprintf(stderr, "Error occurred while storing the reverse prefix array");
				exit(1);
			}
			if (affixArray->esa->affixLink != NULL && !argOScreen) {
				free(affixArray->esa->affixLink);
				affixArray->esa->affixLink = NULL;
			}
		}
	}

	if(affixArray->erpa->xskp != NULL && !argOScreen) {
		free(affixArray->erpa->xskp);
		affixArray->erpa->xskp = NULL;
	}

	if(argAfrpref) {
		if(affixArray->esa->xlcp != NULL && affixArray->esa->xskp != NULL)
			constructAflk(affixArray->erpa, affixArray->esa, false, (int (*)(void *, void *, int, int, void *, void *))(findAffixLink), affixArray->multiSeq->convSequences, affixArray->length, affixArray->alphabet, argShowTimes);
		else
			constructAflkBS(affixArray->erpa, affixArray->esa, false, (int* (*)(void *, int, int, int, int, void *))(findPatternBS), affixArray->multiSeq->convSequences, affixArray->length, affixArray, argShowTimes);

		if (fileName != NULL) {
			if (saveAflk(affixArray->erpa, affixArray->length, 0, fileName)) {
				fprintf(stderr, "Error occurred while storing the reverse prefix array");
				exit(1);
			}
			if (affixArray->erpa->affixLink != NULL && !argOScreen) {
				free(affixArray->erpa->affixLink);
				affixArray->erpa->affixLink = NULL;
			}
		}
	}

}

void computeEArray(EArray *earray, Alphabet *alphabet, unsigned char *seq, int numSeqs, int length, bool isEsa,
					bool argXarray, bool argXlcp, bool argXlcpTree, bool argXskp, bool argShowTimes) {
	struct timeval start, end;

	if(!isEsa) {
		if(argXarray || argXlcp || argXskp)
			reverseSequences(seq, length);
		else
			return;
	}

	if(argXarray) {
		if(isEsa) {
			printf("\nComputing suf... "); fflush(stdout);
		} else {
			printf("\nComputing sufr... "); fflush(stdout);
		}
		if (argShowTimes) gettimeofday(&start, NULL);
		affXSuf(earray, seq, length);
		if (argShowTimes) gettimeofday(&end, NULL);
		printf("done\n");
		if (argShowTimes) printf("Time: %f ms\n", (double)( (end.tv_sec - start.tv_sec)*1000 + (end.tv_usec - start.tv_usec)/1000.0) );
	}

	if(argXlcp) {
		if(isEsa) {
			printf("\nComputing lcp... "); fflush(stdout);
		} else {
			printf("\nComputing lcpr... "); fflush(stdout);
		}
		if (argShowTimes) gettimeofday(&start, NULL);
		affXLcp(earray, alphabet, seq, length);
		if (argShowTimes) gettimeofday(&end, NULL);
		printf("done\n");
		if (argShowTimes) printf("Time: %f ms\n", (double)( (end.tv_sec - start.tv_sec)*1000 + (end.tv_usec - start.tv_usec)/1000.0) );
	}

	/* Modified by Benjamin Albrecht, Jan. 2012 */
	if (argXlcpTree) {
		if (isEsa) {
			printf("\nComputing lcpTree... ");
			fflush(stdout);
		} else {
			printf("\nComputing lcprTree... ");
			fflush(stdout);
		}
		if (argShowTimes)
			gettimeofday(&start, NULL);
		affXLcpTree(earray, isEsa, length);
		if (argShowTimes)
			gettimeofday(&end, NULL);
		printf("done\n");
		if (argShowTimes)
			printf("Time: %f ms\n", (double) ((end.tv_sec - start.tv_sec)
					* 1000 + (end.tv_usec - start.tv_usec) / 1000.0));
	}

	if(argXskp) {
		if(isEsa) {
			printf("\nComputing skp... "); fflush(stdout);
		} else {
			printf("\nComputing skpr... "); fflush(stdout);
		}
		if (argShowTimes) gettimeofday(&start, NULL);
		affXSkp(earray, length, numSeqs);
		if (argShowTimes) gettimeofday(&end, NULL);
		printf("done\n");
		if (argShowTimes) printf("Time: %f ms\n", (double)( (end.tv_sec - start.tv_sec)*1000 + (end.tv_usec - start.tv_usec)/1000.0) );
	}

	if(!isEsa)
		reverseSequences(seq, length);
}

//Reverse sequences terminating with $
void reverseSequences(unsigned char *seq, int length) {
	int i;
	char temp;

	for(i = 0; i < length / 2; i++) {
		temp = seq[i];
		seq[i] = seq[length - 2 - i];
		seq[length - 2 - i] = temp;
	}
}

//Reverse arbitrary string and return a new pointer
unsigned char* reverseStringNewPointer(unsigned char *seq, unsigned int length) {
	unsigned int i, lengthBy2 = length / 2;
	unsigned char *newSeq = (unsigned char *) malloc((length + 1) * sizeof(unsigned char));

	for(i = 0; i < lengthBy2; i++) {
		newSeq[i] = seq[length - 1 - i];
		newSeq[length - 1 - i] = seq[i];
	}
	if (length % 2 != 0)
		newSeq[lengthBy2] = seq[lengthBy2];
	newSeq[length] = 0;

	return newSeq;
}

#define swap(x, a, b) { temp = x[a]; \
                     x[a] = x[b]; x[b] = temp; }

void sortLcpException(int *index, int *value, int left, int right) {
	int i, last, temp;
	if(left >= right)
		return;
	//i = (left + right)/2;
	i = (left + right) >> 1;
	swap(index, left, i);
	swap(value, left, i);
	last = left;
	for(i = left+1; i <= right; i++){
		if(index[i] < index[left]) {
			++last;
			swap(index, last, i);
			swap(value, last, i);
		}
	}
	swap(index, left, last);
	swap(value, left, last);
	sortLcpException(index, value, left, last-1);
	sortLcpException(index, value, last+1, right);
}

void affXSuf(EArray *earray, unsigned char *seq, int length) {
	if((earray->xarray = (int *) calloc(length, sizeof(int))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"xarray\".\n");

	if(divsufsort(seq, earray->xarray, (saidx_t)length) != 0) {
		fprintf(stderr, "Cannot allocate memory.\n");
		exit(1);
	}
}

int xlcpExceptionValue(LcpException *xlcpException, int k) {
	int low, high, mid;

	/* naive search */
	/*unsigned int i;
	for(i = 0; i < xlcpException->numExceptions; i++)
		if(xlcpException->index[i] == k)
			return xlcpException->value[i];*/

	low = 0;
	high = xlcpException->numExceptions - 1;

	while (low <= high) {
		mid = (low + high) >> 1;
		if (xlcpException->index[mid] < k)
			low = mid + 1;
		else if (xlcpException->index[mid] > k)
			high = mid - 1;
		else
			return xlcpException->value[mid];
	}

	return -1;
}

void affXLcp(EArray *earray, Alphabet *alphabet, unsigned char *seq, int length) {
	if(earray->xarray == NULL) {
		fprintf(stderr,"lcp (lcpr) computation requires table suf (sufr).\n");
		exit(1);
	}

	int i, j, k, h=0;
	int *rank = (int *) calloc(length, sizeof(int));

	if((earray->xlcp = (unsigned char *) calloc(length, sizeof(unsigned char))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"xlcp\".\n");
	if((earray->xlcpException = (LcpException *) malloc(sizeof(LcpException))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"xlcpException\".\n");

	earray->xlcpException->index = (int *) calloc(BUFFER2, sizeof(int));
	earray->xlcpException->value = (int *) calloc(BUFFER2, sizeof(int));
	earray->xlcpException->numExceptions = 0;

	for(i = 0; i < length; i++)
		rank[earray->xarray[i]] = i;

	for(i = 0; i < length; i++){
		k = rank[i];
		if(k == 0)
			earray->xlcp[k] = 0;
		else {
			j = earray->xarray[k - 1];
			while(seq[i + h] != $ && seq[j + h] != $ && seq[i + h] == seq[j + h] && !alphabet->isWildCard[seq[i + h] - 1] && !alphabet->isWildCard[seq[j + h] - 1])
				h++;

			if(h > 254) {
				earray->xlcp[k] = 255;

				if(earray->xlcpException->numExceptions % BUFFER2 == 0) {
					if((earray->xlcpException->index = realloc(earray->xlcpException->index, (earray->xlcpException->numExceptions + BUFFER2)*sizeof(int))) == NULL) {
						fprintf(stderr,"Memory allocation failed for \"xlcpException.index\".\n");
						exit(1);
					}
					if((earray->xlcpException->value = realloc(earray->xlcpException->value, (earray->xlcpException->numExceptions + BUFFER2)*sizeof(int))) == NULL) {
						fprintf(stderr,"Memory allocation failed for \"xlcpException.value\".\n");
						exit(1);
					}
				}
				earray->xlcpException->value[earray->xlcpException->numExceptions] = h;
				earray->xlcpException->index[earray->xlcpException->numExceptions] = k;
				earray->xlcpException->numExceptions++;
			} else
				earray->xlcp[k] = h;
		}
		if(h > 0)
			h--;
	}

	free(rank);

	if(earray->xlcpException->numExceptions > 0)
		sortLcpException(earray->xlcpException->index, earray->xlcpException->value, 0, earray->xlcpException->numExceptions - 1);
}

/* Modified by Benjamin Albrecht, Jan. 2012 */
void setLcpTreeValue(unsigned char *xlcpTree, int index, int value,
		LcpException *xlcpTreeException) {

	if (value > 254) {
		xlcpTree[index] = 255;

		if (xlcpTreeException->numExceptions % BUFFER2 == 0) {
			if ((xlcpTreeException->index = realloc(xlcpTreeException->index,
					(xlcpTreeException->numExceptions + BUFFER2) * sizeof(int)))
					== NULL) {
				fprintf(stderr,
						"Memory allocation failed for \"xlcpException.index\".\n");
				exit(1);
			}
			if ((xlcpTreeException->value = realloc(xlcpTreeException->value,
					(xlcpTreeException->numExceptions + BUFFER2) * sizeof(int)))
					== NULL) {
				fprintf(stderr,
						"Memory allocation failed for \"xlcpException.value\".\n");
				exit(1);
			}
		}

		xlcpTreeException->value[xlcpTreeException->numExceptions] = value;
		xlcpTreeException->index[xlcpTreeException->numExceptions] = index;
		xlcpTreeException->numExceptions++;

		int i = xlcpTreeException->numExceptions - 1;
		while (i > 0 && xlcpTreeException->index[i] < xlcpTreeException->index[i - 1]) {
			int temp = xlcpTreeException->index[i];
			xlcpTreeException->index[i] = xlcpTreeException->index[i - 1];
			xlcpTreeException->index[i - 1] = temp;
			temp = xlcpTreeException->value[i];
			xlcpTreeException->value[i] = xlcpTreeException->value[i - 1];
			xlcpTreeException->value[i - 1] = temp;
			i--;
		}

	} else
		xlcpTree[index] = value;
}

/* Modified by Benjamin Albrecht, Jan. 2012 */
int affXLcpTreeRec(int b, int m, EArray *earray, unsigned char *xlcpTree,
		LcpException *xlcpTreeException) {

	int index = (int) (b + m) / 2;
	int value;

	if (abs(b - m) == 1) {
		int last = b > m ? b : m;
		return lcpvalueX(earray, last);
	} else if (abs(b - m) == 2) {
		int last = b > m ? b : m;
		value = lcpvalueX(earray, last) < lcpvalueX(earray, last-1) ? lcpvalueX(earray, last) : lcpvalueX(earray, last-1);
	} else {
		int newM = (int) ((b + m) / 2);
		if (b < m) {
			int lcpLeft = affXLcpTreeRec(b, newM, earray, xlcpTree,
					xlcpTreeException);
			int lcpRight = affXLcpTreeRec(newM, m, earray, xlcpTree,
					xlcpTreeException);
			value = lcpLeft < lcpRight ? lcpLeft : lcpRight;
		} else {
			int lcpLeft = affXLcpTreeRec(m, newM, earray, xlcpTree,
					xlcpTreeException);
			int lcpRight = affXLcpTreeRec(newM, b, earray, xlcpTree,
					xlcpTreeException);
			value = lcpLeft < lcpRight ? lcpLeft : lcpRight;
		}
	}

	setLcpTreeValue(xlcpTree, index, value, xlcpTreeException);

	unsigned int result = xlcpTree[index];
	return result;

}

/* Modified by Benjamin Albrecht, Jan. 2012 */
void affXLcpTree(EArray *earray, bool isEsa, int length) {

	if (earray->xlcp == NULL) {
		fprintf(stderr,
				"Computation of lcpTree requires table lcp.\n");
		exit(1);
	}

	if ((earray->xlcpTree = (unsigned char *) calloc(length, sizeof(unsigned char))) == NULL)
		fprintf(stderr, "Memory allocation failed for \"xlcpTree\".\n");
	if ((earray->xlcpTreeException = (LcpException *) malloc( sizeof(LcpException))) == NULL)
		fprintf(stderr, "Memory allocation failed for \"xlcpTreeException\".\n");

	earray->xlcpTreeException->numExceptions = 0;
	earray->xlcpTreeException ->index = (int *) calloc(BUFFER2, sizeof(int));
	earray->xlcpTreeException ->value = (int *) calloc(BUFFER2, sizeof(int));

	int m = (int) ((0 + length - 1) / 2);
	int lcpLeft = affXLcpTreeRec(0, m, earray, earray->xlcpTree,
			earray->xlcpTreeException);
	int lcpRight = affXLcpTreeRec(m, length - 1, earray, earray->xlcpTree,
			earray->xlcpTreeException);

	int value = lcpLeft < lcpRight ? lcpLeft : lcpRight;
	setLcpTreeValue(earray->xlcpTree, m, value, earray->xlcpTreeException);

}

void affXSkp(EArray *earray, int length, int numSeqs) {
	if(earray->xarray == NULL || earray->xlcp == NULL) {
		fprintf(stderr,"skp (skpr) computation requires tables suf and lcp (sufr and lcpr).\n");
		exit(1);
	}

	int i, j;

	if((earray->xskp = (int *) calloc(length, sizeof(int))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"xskp\".\n");

	/*Naive computation*/
	/*xskp[0] = length;
	for(i = 1; i < length; i++) {
		xskp[i] = length;
		j = i + 1;
		while(j < length && lcpvalue(i) <= lcpvalue(j)) {
			j++;
		};
		if(xskp[i] > j)
			xskp[i] = j;
	}*/

	/*Computation using a simulated stack*/
	int *stack, si = 0;

	stack = (int *) calloc(BUFFER1, sizeof(int));

	earray->xskp[0] = length;
	for(i = 1; i < length; i++) {
		stack[si++] = i;
		j = i + 1;
		while(j < length && si > 0) {
			while(j < length && xlcpvalue(stack[si-1]) <= xlcpvalue(j)) {
				if(si % BUFFER1 == 0)
					if((stack = realloc(stack, (si + BUFFER1) * sizeof(int))) == NULL) {
						fprintf(stderr,"Memory allocation failed for \"stack\".\n");
						exit(1);
					}
				stack[si++] = j;
				j++;
			}
			//if j == length, lcpvalue(j) is out of bounds
			while(si > 0 && (j == length /*- numSeqs*/ || xlcpvalue(stack[si-1]) > xlcpvalue(j)))
				earray->xskp[stack[--si]] = j;
		}
		i = j-1;
	}
	free(stack);
}
