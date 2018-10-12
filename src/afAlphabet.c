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
#include <ctype.h>

#include "af.h"

void setPredefinedAlphabet(int type, Alphabet *alphabet) {
	int i, j;
	char *dna     = "AaCcGgTt";
	char *rna     = "AaCcGgUu";
	char *protein = "AaCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtVvWwYy";

	/********* DNA or RNA *********/
	if(type == 0 || type == 1) {
		alphabet->numClasses = 4;
		alphabet->eqClass    = (char **) calloc(alphabet->numClasses, sizeof(char*));
		alphabet->classSize  = (int *) calloc(alphabet->numClasses, sizeof(int));
		alphabet->classRepresentative = (unsigned char *) calloc(alphabet->numClasses, sizeof(unsigned char));
		alphabet->isWildCard = (bool *) calloc(alphabet->numClasses, sizeof(bool));

		for(i = 0, j = 0; i < alphabet->numClasses; i++) {
			alphabet->classSize[i] = 2;
			alphabet->isWildCard[i]= 0;
			alphabet->eqClass[i]   = (char *) calloc(2, sizeof(char));
			if(type == 0) {
				alphabet->eqClass[i][0] = dna[j];
				alphabet->eqClass[i][1] = dna[j + 1];
				alphabet->classRepresentative[i] = dna[j];
			} else {
				alphabet->eqClass[i][0] = rna[j];
				alphabet->eqClass[i][1] = rna[j + 1];
				alphabet->classRepresentative[i] = rna[j];
			}
			j = j + 2;
		}
	}

	/********* Protein *********/
	else if(type == 2) {
		alphabet->numClasses = 20;
		alphabet->eqClass    = (char **) calloc(alphabet->numClasses, sizeof(char*));
		alphabet->classSize  = (int *) calloc(alphabet->numClasses, sizeof(int));
		alphabet->classRepresentative = (unsigned char *) calloc(alphabet->numClasses, sizeof(unsigned char));
		alphabet->isWildCard = (bool *) calloc(alphabet->numClasses, sizeof(bool));

		for(i = 0, j = 0; i < alphabet->numClasses; i++) {
			alphabet->classSize[i] = 2;
			alphabet->isWildCard[i] = 0;
			alphabet->eqClass[i] = (char *) calloc(2, sizeof(char));
			alphabet->eqClass[i][0] = protein[j];
			alphabet->eqClass[i][1] = protein[j + 1];
			alphabet->classRepresentative[i] = protein[j];
			j = j + 2;
		}
	}
}

bool convertToAlphabet(unsigned char *seq, unsigned char *convSeq, int length, bool isPattern, Alphabet* alphabet) {
	int i, j, iS;
	unsigned char *table;

	if (alphabet == NULL) {
		fprintf(stderr, "Error. Alphabet was not initialized. - %s %d.\n", __FILE__, __LINE__);
		return EXIT_FAILURE;
	}

	if((table = (unsigned char *) calloc(256, sizeof(unsigned char))) == NULL) {
		fprintf(stderr, "Memory allocation failed. - %s %d.\n", __FILE__, __LINE__);
		return EXIT_FAILURE;
	}

	for (i = 0; i < 256; i++) {
		table[i] = 0;
	}

	for(i = 0; i < alphabet->numClasses; i++) {
		for(j = 0; j < alphabet->classSize[i]; j++) {
			table[(unsigned int) alphabet->eqClass[i][j]] = i + 1; // Characters in the 1st eqClass receive value 1, in the 2nd eqClass value 2...
		}
	}
	table['$'] = $;
	table[$] = $;

	for(iS = 0; iS < length; iS++) {

		if (isPattern) {
			if(seq[iS] < 127 && alphabet->iupacTable[(int) seq[iS]] == NULL && table[(unsigned char) seq[iS]] == 0) {
				fprintf(stderr,"Character \"%c\" is not defined in the alphabet.\n", seq[iS]);
				return 1;
			}
		} else {
			if(seq[iS] < 127 && table[(unsigned char) seq[iS]] == 0) {
				fprintf(stderr,"Character \"%c\" is not defined in the alphabet.\n", seq[iS]);
				return 1;
			}
		}

		if(isPattern && alphabet->iupacTable[(int) seq[iS]] != NULL)
			convSeq[iS] = seq[iS];
		else {
			convSeq[iS] = table[(unsigned int) seq[iS]];
		}
	}

	return EXIT_SUCCESS;
}

/* Converts each value of an alphabetically transformed array to the
 * value of the representative character of the corresponding class. */
bool convertToRepresentative(char *seq, int length, Alphabet *alphabet) {
	int i;

	for(i = 0; i < length; i++) {
		if(seq[i] == $)
			seq[i] = $;
		else if(seq[i] == '*')
			seq[i] = '*';
		else
			seq[i] = alphabet->classRepresentative[seq[i] - 1];
	}

	return 0;
}

bool** setIupacTable (Alphabet *alphabet) {
	bool **table;
	unsigned int i, j;
	char *characters  = "RYSWKMBDHVN*";
	char *characters2 = "ryswkmbdhvn";

	if((table = (bool **) calloc(127, sizeof(bool*))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"iupacTable\". - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	alphabet->iupacTable = table;

	for (i = 0; i < 127; i++) {
		table[i] = NULL;
	}

	for (i = 0; i < 12; i++) {
		if ((table[(int) characters[i]] = (bool *) calloc(127, sizeof(bool))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"iupacTable\". - %s %d.\n", __FILE__, __LINE__);
			exit(1);
		}
		for (j = 0; j < 127; j++) {
			table[(int) characters[i]][j] = false;
		}
	}

	for (i = 0; i < 127; i++) {
		switch (i) {
			case 'R':
				table[i]['A'] = table[i]['a'] = true;
				table[i]['G'] = table[i]['g'] = true;
				break;
			case 'Y':
				table[i]['C'] = table[i]['c'] = true;
				table[i]['T'] = table[i]['t'] = true;
				table[i]['U'] = table[i]['u'] = true;
				break;
			case 'S':
				table[i]['G'] = table[i]['g'] = true;
				table[i]['C'] = table[i]['c'] = true;
				break;
			case 'W':
				table[i]['A'] = table[i]['a'] = true;
				table[i]['T'] = table[i]['t'] = true;
				table[i]['U'] = table[i]['u'] = true;
				break;
			case 'K':
				table[i]['G'] = table[i]['g'] = true;
				table[i]['T'] = table[i]['t'] = true;
				table[i]['U'] = table[i]['u'] = true;
				break;
			case 'M':
				table[i]['A'] = table[i]['a'] = true;
				table[i]['C'] = table[i]['c'] = true;
				break;
			case 'B':
				table[i]['C'] = table[i]['c'] = true;
				table[i]['G'] = table[i]['g'] = true;
				table[i]['T'] = table[i]['t'] = true;
				table[i]['U'] = table[i]['u'] = true;
				break;
			case 'D':
				table[i]['A'] = table[i]['a'] = true;
				table[i]['G'] = table[i]['g'] = true;
				table[i]['T'] = table[i]['t'] = true;
				table[i]['U'] = table[i]['u'] = true;
				break;
			case 'H':
				table[i]['A'] = table[i]['a'] = true;
				table[i]['C'] = table[i]['c'] = true;
				table[i]['T'] = table[i]['t'] = true;
				table[i]['U'] = table[i]['u'] = true;
				break;
			case 'V':
				table[i]['A'] = table[i]['a'] = true;
				table[i]['C'] = table[i]['c'] = true;
				table[i]['G'] = table[i]['g'] = true;
				break;
			case 'N':
				for (j = 0; j < 127; j++) {
					table[i][j] = true;
				}
				break;
			case '*':
				for (j = 0; j < 127; j++) {
					table[i][j] = true;
				}
		}
	}

	for (i = 0; i < 11; i++) {
		table[(int) characters2[i]] = table[(int) toupper(characters2[i])];
	}

	return table;
}

