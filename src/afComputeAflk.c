/*Copyright (C) 2012  Benjamin Albrecht

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

#include "af.h"

/*
 * Benjamin Albrecht, Jan. 2012
 * A normal binary search for finding a pattern in esa or erpa.
 * As result the left border of the corresponding interval is returned.
*/
int* cmpAflkViaBS(int l, int r, AffixArray *affixArray,
		unsigned char * pattern, int pLength, bool isEsa) {

	int mid = ((unsigned int) l + (unsigned int) r) >> 1;
	int searchDir = 0;

	EArray *earray;
	if (isEsa)
		earray = affixArray->esa;
	else
		earray = affixArray->erpa;

	//compute new search-direction
	int i = 0;
	while (i < pLength) {
		char seqChar, patChar;
		if (isEsa) {
			seqChar
					= affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[earray->xarray[mid]
							+ i] - 1];
			patChar = pattern[i];
		} else {
			seqChar
					= affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[(affixArray->length
							- 2) - (earray->xarray[mid] + i)] - 1];
			patChar = pattern[(pLength - 1) - i];
		}

		if (((int) seqChar != (int) patChar) && ((int) seqChar >= 65)
				&& ((int) seqChar <= 90)) {
			searchDir = (int) seqChar < (int) patChar;
			break;
		} else if (((int) seqChar == (int) patChar)) {
			i++;
		} else {
			searchDir = 0;
			break;
		}

	}

	if (mid != r && mid != l) {
		if (i == pLength) {
			//pattern found, now compute left border of the lcp-interval
			int leftLCPborder = mid;
			while (leftLCPborder >= 0) {
				int val = lcpvalueX(earray, leftLCPborder);
				if (val < pLength) {
					char c;
					if (isEsa) {
						while ((c
								= affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[earray->xarray[leftLCPborder] + val] - 1]) == '*')
							val++;
					} else {
						while ((c = affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[(affixArray->length - 1)
										- affixArray->erpa->xarray[leftLCPborder]
										- val] - 1]) == '*')
							val++;
					}
				}
				if (val < pLength)
					break;
				else
					leftLCPborder--;
			}

			int *result = (int *) malloc(sizeof(int));
			result[0] = leftLCPborder;

			return result;

		} else if (searchDir == 1) {
			//continue search in the right part of the interval
			return cmpAflkViaBS(mid, r, affixArray, pattern, pLength, isEsa);
		} else {
			//continue search in the left part of the interval
			return cmpAflkViaBS(l, mid, affixArray, pattern, pLength, isEsa);
		}
	} else {
		//pattern not found, return -1
		int *result = (int *) malloc(sizeof(int));
		result[0] = -1;
		return result;
	}
}

/*
 * Benjamin Albrecht, Jan. 2012
 * An improved binary search for finding a pattern in esa and erpa.
 * The search makes use of the precomputed tables xlcpTree, xlcpTreeException.
 * As result the left border of the corresponding interval is returned.
 * The algorithm follows the publication of Myers and Manber 1993.
*/
int* cmpAflkViaImprovedBS(int l, int r, int lComp, int rComp,
		AffixArray *affixArray, const unsigned char * pattern, int pLength,
		bool isEsa) {

	int mid = (int) (l + r) / 2;

	EArray *earray;
	if (isEsa)
		earray = affixArray->esa;
	else
		earray = affixArray->erpa;

	//only for debugging
	int lMem = l;
	int rMem = r;

	//compute new search-direction
	int i = lComp + 1;
	int j = rComp + 1;

	if (lComp == rComp) {//Case 1
		while (i < pLength) {
			char seqChar, patChar;
			if (isEsa) {
				seqChar = affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[earray->xarray[mid] + i] - 1];
				patChar = pattern[i];
			} else {
				seqChar = affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[(affixArray->length - 2) - (earray->xarray[mid] + i)] - 1];
				patChar = pattern[(pLength - 1) - i];
			}
			if (((int) seqChar != (int) patChar) && ((int) seqChar >= 65)
					&& ((int) seqChar <= 90)) {
				if ((int) seqChar < (int) patChar) {
					lComp = i - 1;
					l = mid;
				} else if ((int) seqChar > (int) patChar) {
					rComp = i - 1;
					r = mid;
				}
				break;
			} else if (((int) seqChar == (int) patChar)) {
				i++;
			} else {
				rComp = i - 1;
				r = mid;
				break;
			}
		}
	} else if (lComp > rComp) {//Case 2
		int lcpLM = mid - l > 1 ? lcpTreevalueX(earray, (int) (l + mid) / 2)
				: lcpvalueX(earray, mid);
		if (lcpLM > lComp + 1) {//Case 2a
			l = mid;
		} else if (lcpLM == lComp + 1) {//Case 2b
			while (i < pLength) {
				char seqChar, patChar;
				if (isEsa) {
					seqChar = affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[earray->xarray[mid] + i] - 1];
					patChar = pattern[i];
				} else {
					seqChar = affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[(affixArray->length - 2) - (earray->xarray[mid] + i)] - 1];
					patChar = pattern[(pLength - 1) - i];
				}

				if (((int) seqChar != (int) patChar) && ((int) seqChar >= 65)
						&& ((int) seqChar <= 90)) {
					if ((int) seqChar < (int) patChar) {
						lComp = i - 1;
						l = mid;
					} else if ((int) seqChar > (int) patChar) {
						rComp = i - 1;
						r = mid;
					}
					break;
				} else if (((int) seqChar == (int) patChar)) {
					i++;
				} else {
					rComp = i - 1;
					r = mid;
					break;
				}
			}
		} else {//Case 2c
			r = mid;
			rComp = lcpLM - 1;
		}
	} else {//Case 3
		int lcpRM = r - mid > 1 ? lcpTreevalueX(earray, (int) (r + mid) / 2)
				: lcpvalueX(earray, r);
		if (lcpRM > rComp + 1) {//Case 3a
			r = mid;
		} else if (lcpRM == rComp + 1) {//Case 3b
			while (j < pLength) {
				char seqChar, patChar;
				if (isEsa) {
					seqChar = affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[earray->xarray[mid] + j] - 1];
					patChar = pattern[j];
				} else {
					seqChar = affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[(affixArray->length - 2) - (earray->xarray[mid] + j)] - 1];
					patChar = pattern[(pLength - 1) - j];
				}

				if (((int) seqChar != (int) patChar) && ((int) seqChar >= 65)
						&& ((int) seqChar <= 90)) {
					if ((int) seqChar < (int) patChar) {
						lComp = j - 1;
						l = mid;
					} else if ((int) seqChar > (int) patChar) {
						rComp = j - 1;
						r = mid;
					}
					break;
				} else if (((int) seqChar == (int) patChar)) {
					j++;
				} else {
					rComp = j - 1;
					r = mid;
					break;
				}
			}
		} else {//Case 3c
			l = mid;
			lComp = lcpRM - 1;
		}
	}

	if (l != lMem || r != rMem || i == pLength || j == pLength) {
		if (i == pLength || j == pLength) {
			//pattern found, now compute left border of the lcp-interval
			int leftLCPborder = mid;
			while (leftLCPborder >= 0) {
				int val = lcpvalueX(earray, leftLCPborder);
				if (val < pLength) {
					char c;
					if (isEsa) {
						while ((c = affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[earray->xarray[leftLCPborder] + val] - 1]) == '*')
							val++;
					} else {
						while ((c = affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[(affixArray->length - 1)
						             - affixArray->erpa->xarray[leftLCPborder] - val] - 1]) == '*')
							val++;
					}
				}
				if (val < pLength)
					break;
				else
					leftLCPborder--;
			}
			int *result = (int *) malloc(sizeof(int));
			result[0] = leftLCPborder;
			return result;
		} else{
			//continue search in the new interval [l,r]
			return cmpAflkViaImprovedBS(l, r, lComp, rComp, affixArray, pattern, pLength, isEsa);
		}
	} else {
		//pattern not found, return -1
		int *result = (int *) malloc(sizeof(int));
		result[0] = -1;
		return result;
	}
}
