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
#include <limits.h>
#include <sys/stat.h>
#include <unistd.h> //includes function getpagesize()
#include <fcntl.h> //includes function open()
#include <sys/mman.h>

#include "af.h"

/* Initializes affixArray. Parameter bConstruct must be set to true if affixArray will be used
 * for affix array construction */
void init(AffixArray *affixArray, Alphabet *alphabet, char *argAlphabetFile, int predefAlphabet,
			bool *argSuf, bool *argLcp, bool *argLcpTree, bool *argSkp, bool *argAfsuf,
			bool *argRpref, bool *argRlcp, bool *argRlcpTree, bool *argRskp, bool *argAfrpref, bool bConstruct) {

	/* Consistency checking */
	/* Modified by Benjamin Albrecht, Jan. 2012  */
	if(*argLcp || *argSkp || *argAfsuf || *argAfrpref || *argLcpTree)
		*argSuf = 1;
	if(*argSkp || *argAfsuf || *argLcpTree)
		*argLcp = 1;
	if(*argRlcp || *argRskp || *argAfrpref || *argAfsuf || *argRlcpTree)
		*argRpref = 1;
	if(*argRskp || *argAfrpref || *argRlcpTree)
		*argRlcp = 1;

	affixArray->alphabet = alphabet;

	if((affixArray->multiSeq = (MultiSeq *) malloc(sizeof(MultiSeq))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"affixArray->multiSeq\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	affixArray->multiSeq->convSequences = NULL;
	affixArray->multiSeq->convSequencesMmapped = false;

	if (*argSuf || bConstruct) {
		if((affixArray->esa = (EArray *) malloc(sizeof(EArray))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"esa\" - %s %d.\n", __FILE__, __LINE__);
			exit(1);
		}
		if((affixArray->erpa = (EArray *) malloc(sizeof(EArray))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"erpa\" - %s %d.\n", __FILE__, __LINE__);
			exit(1);
		}

		affixArray->esa->xarray = NULL;
		affixArray->esa->xlcp = NULL;
		affixArray->esa->xlcpException = NULL;
		affixArray->esa->xskp = NULL;
		affixArray->esa->affixLink = NULL;

		affixArray->erpa->xarray = NULL;
		affixArray->erpa->xlcp = NULL;
		affixArray->erpa->xlcpException = NULL;
		affixArray->erpa->xskp = NULL;
		affixArray->erpa->affixLink = NULL;

		/* Modified by Benjamin Albrecht, Jan. 2012  */
		affixArray->esa->xlcpTree = NULL;
		affixArray->esa->xlcpTreeException = NULL;
		affixArray->erpa->xlcpTree= NULL;
		affixArray->erpa->xlcpTreeException = NULL;

	} else {
		affixArray->esa  = NULL;
		affixArray->erpa = NULL;
    }

	if (alphabet != NULL) {
		if (predefAlphabet > -1) {
			if (argAlphabetFile != NULL) {
				loadAlphabetFile(argAlphabetFile, alphabet);
			} else {
				if (predefAlphabet != 0 && predefAlphabet != 1 && predefAlphabet != 2)
					setPredefinedAlphabet(0, alphabet);
				else
					setPredefinedAlphabet(predefAlphabet, alphabet);
			}
		}
		setIupacTable(alphabet);
	}

}

void freeAll(AffixArray *affixArray, MultiPattern *multiPattern, bool **complementarityRules) {
	unsigned int ui;
	char *iupac_characters  = "RYSWKMBDHVN*";

	if (affixArray->multiSeq != NULL) {
		if (affixArray->multiSeq->convSequences != NULL) {
			if (affixArray->multiSeq->convSequencesMmapped)
				munmap(affixArray->multiSeq->convSequences, affixArray->length * sizeof(unsigned char));
			else
				free(affixArray->multiSeq->convSequences);
		}
		if (affixArray->multiSeq->seqDescLength != NULL)
			free(affixArray->multiSeq->seqDescLength);
		if (affixArray->multiSeq->seqDescription != NULL) {
			for (ui = 0; ui < affixArray->multiSeq->numSeqs; ui++) {
				free(affixArray->multiSeq->seqDescription[ui]);
			}
			free(affixArray->multiSeq->seqDescription);
		}
		if (affixArray->multiSeq->seqEndPos != NULL)
			free(affixArray->multiSeq->seqEndPos);
		if (affixArray->multiSeq->sequences != NULL) {
			if (affixArray->multiSeq->sequencesMmapped)
				munmap(affixArray->multiSeq->sequences, affixArray->length * sizeof(unsigned char));
			else
				free(affixArray->multiSeq->sequences);
		}
		free(affixArray->multiSeq);
	}

	if (affixArray->alphabet != NULL) {
		if (affixArray->alphabet->classRepresentative)
			free(affixArray->alphabet->classRepresentative);
		if (affixArray->alphabet->classSize)
			free(affixArray->alphabet->classSize);
		if (affixArray->alphabet->eqClass) {
			for (ui = 0; ui < affixArray->alphabet->numClasses; ui++) {
				free(affixArray->alphabet->eqClass[ui]);
			}
			free(affixArray->alphabet->eqClass);
		}
		if (affixArray->alphabet->isWildCard)
			free(affixArray->alphabet->isWildCard);

		if (complementarityRules != NULL) {
			for (ui = 0; ui < affixArray->alphabet->numClasses; ui++) {
					free(complementarityRules[ui]);
			}
			free(complementarityRules);
		}

		if (affixArray->alphabet->iupacTable != NULL) {
			for (ui = 0; ui < 12; ui++) {
				free(affixArray->alphabet->iupacTable[(int) iupac_characters[ui]]);
			}
			free(affixArray->alphabet->iupacTable);
		}
	}

	if (affixArray->esa != NULL) {
		if (affixArray->esa->affixLink != NULL)
			munmap(affixArray->esa->affixLink, affixArray->length * sizeof(int));
		if (affixArray->esa->xarray != NULL)
			munmap(affixArray->esa->xarray, affixArray->length * sizeof(int));
		if (affixArray->esa->xlcp != NULL)
			munmap(affixArray->esa->xlcp, affixArray->length * sizeof(unsigned char));
		if (affixArray->esa->xlcpException != NULL) {
			if (affixArray->esa->xlcpException->numExceptions > 0) {
				free(affixArray->esa->xlcpException->index);
				free(affixArray->esa->xlcpException->value);
			}
			free(affixArray->esa->xlcpException);
		}

		/* Modified by Benjamin Albrecht, Jan. 2012  */
		if (affixArray->esa->xlcpTree != NULL)
					munmap(affixArray->esa->xlcpTree, affixArray->length * sizeof(unsigned char));
		if (affixArray->esa->xlcpTreeException != NULL) {
			if (affixArray->esa->xlcpTreeException->numExceptions > 0) {
				free(affixArray->esa->xlcpTreeException->index);
				free(affixArray->esa->xlcpTreeException->value);
			}
			free(affixArray->esa->xlcpTreeException);
		}

		if (affixArray->esa->xskp != NULL)
			munmap(affixArray->esa->xskp, affixArray->length * sizeof(int));
		free(affixArray->esa);
	}

	if (affixArray->erpa != NULL) {
		if (affixArray->erpa->affixLink != NULL)
			munmap(affixArray->erpa->affixLink, affixArray->length * sizeof(int));
		if (affixArray->erpa->xarray != NULL)
			munmap(affixArray->erpa->xarray, affixArray->length * sizeof(int));
		if (affixArray->erpa->xlcp != NULL)
			munmap(affixArray->erpa->xlcp, affixArray->length * sizeof(unsigned char));
		if (affixArray->erpa->xlcpException != NULL) {
			if (affixArray->erpa->xlcpException->numExceptions > 0) {
				free(affixArray->erpa->xlcpException->index);
				free(affixArray->erpa->xlcpException->value);
			}
			free(affixArray->erpa->xlcpException);
		}

		/* Modified by Benjamin Albrecht, Jan. 2012  */
		if (affixArray->erpa->xlcpTree != NULL)
			munmap(affixArray->erpa->xlcpTree, affixArray->length * sizeof(unsigned char));
		if (affixArray->erpa->xlcpTreeException != NULL) {
			if (affixArray->erpa->xlcpTreeException->numExceptions > 0) {
				free(affixArray->erpa->xlcpTreeException->index);
				free(affixArray->erpa->xlcpTreeException->value);
			}
			free(affixArray->erpa->xlcpTreeException);
		}

		if (affixArray->erpa->xskp != NULL)
			munmap(affixArray->erpa->xskp, affixArray->length * sizeof(int));
		free(affixArray->erpa);
	}

	if (multiPattern != NULL) {
		if (multiPattern->pattern != NULL) {
			for (ui = 0; ui < multiPattern->numPatterns; ui++) {
				if (multiPattern->pattern[ui].desc != NULL)
					free(multiPattern->pattern[ui].desc);
				if (multiPattern->pattern[ui].seq != NULL)
					free(multiPattern->pattern[ui].seq);
				if (multiPattern->pattern[ui].structure != NULL)
					free(multiPattern->pattern[ui].structure);
			}
			free(multiPattern->pattern);
		}
		free(multiPattern);
	}
}

bool loadFastaFile(AffixArray *affixArray, char *fileName, bool bAllocConvSeq) {
	FILE *fp;
	char temp, *info;
	int i;

	if((fp=fopen(fileName,"r")) == NULL) {
		fprintf(stderr,"Error opening file \"%s\". If this is an index, please select the desired tables.\n", fileName);
		return 1;
	}

	if (fseek(fp, 0L, SEEK_END)) {
		fprintf(stderr,"Error reading file \"%s\".\n", fileName);
		return 1;
	}

	i=ftell(fp);
	if (i == 0) {
		fprintf(stderr, "File \"%s\" is empty.\n", fileName);
    	return 1;
	}

	rewind(fp);

	if(affixArray->multiSeq == NULL) {
		fprintf(stderr, "Error. \"multiSeq\" was not initialized - %s %d.\n", __FILE__, __LINE__);
		return 1;
	}

	affixArray->multiSeq->convSequences  = NULL;
	if ((affixArray->multiSeq->sequences = (unsigned char *) calloc(BUFFER1, sizeof(unsigned char))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"sequences\" - %s %d.\n", __FILE__, __LINE__);
		return 1;
	}
	if ((affixArray->multiSeq->seqEndPos      = (int *) calloc(BUFFER3, sizeof(int))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"seqEndPos\" - %s %d.\n", __FILE__, __LINE__);
		return 1;
	}
	if ((affixArray->multiSeq->seqDescription = (char **) calloc(BUFFER3, sizeof(char*))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"seqDescription\" - %s %d.\n", __FILE__, __LINE__);
		return 1;
	}
	if ((affixArray->multiSeq->seqDescLength  = (int *) calloc(BUFFER3, sizeof(int))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"seqDescLength\" - %s %d.\n", __FILE__, __LINE__);
		return 1;
	}

	affixArray->multiSeq->sequencesMmapped = false;

	affixArray->length = affixArray->multiSeq->numSeqs = 0;
	temp = getc(fp);
	if(temp != '>') {
		fprintf(stderr,"The file does not have a valid fasta format.\n");
		return 1;
	}
	do {
		if(temp == '>') {
			if(affixArray->length > 0) {
				affixArray->multiSeq->sequences[affixArray->length++] = $;
				if(affixArray->length % BUFFER1 == 0) {
					if((affixArray->multiSeq->sequences = realloc(affixArray->multiSeq->sequences, (affixArray->length+BUFFER1)*sizeof(char))) == NULL) {
						fprintf(stderr,"Memory allocation failed for \"sequences\" - %s %d.\n", __FILE__, __LINE__);
						return 1;
					}
				}
				affixArray->multiSeq->seqEndPos[affixArray->multiSeq->numSeqs-1] = affixArray->length - 1;
			}
			i = 0;
			if((info = (char *) calloc(BUFFER3, sizeof(char))) == NULL){
				fprintf(stderr,"Memory allocation failed for \"info\" - %s %d.\n", __FILE__, __LINE__);
				return 1;
			}
			while((temp = getc(fp)) != EOF && temp != '\r' && temp != '\n') {
				info[i] = temp;
				info[i+1] = 0;
				if(((++i)+1) % BUFFER3 == 0 && (info = realloc(info, (i+BUFFER3+1)*sizeof(char))) == NULL) {
					fprintf(stderr,"Memory allocation failed for \"info\" - %s %d.\n", __FILE__, __LINE__);
					return 1;
				}
			}
			if(affixArray->multiSeq->numSeqs % BUFFER3 == 0) {
				if((affixArray->multiSeq->seqEndPos = realloc(affixArray->multiSeq->seqEndPos, (affixArray->multiSeq->numSeqs+BUFFER3)*sizeof(int))) == NULL) {
					fprintf(stderr,"Memory allocation failed for \"seqEndPos\" - %s %d.\n", __FILE__, __LINE__);
					return 1;
				}
				if((affixArray->multiSeq->seqDescription  = realloc(affixArray->multiSeq->seqDescription, (affixArray->multiSeq->numSeqs+BUFFER3)*sizeof(char*))) == NULL) {
					fprintf(stderr,"Memory allocation failed for \"seqInfo\" - %s %d.\n", __FILE__, __LINE__);
					return 1;
				}
				if((affixArray->multiSeq->seqDescLength = realloc(affixArray->multiSeq->seqDescLength, (affixArray->multiSeq->numSeqs+BUFFER3)*sizeof(int))) == NULL) {
					fprintf(stderr,"Memory allocation failed for \"seqDescLength\" - %s %d.\n", __FILE__, __LINE__);
					return 1;
				}
			}
			affixArray->multiSeq->seqDescription[affixArray->multiSeq->numSeqs] = info;
			affixArray->multiSeq->seqDescLength[affixArray->multiSeq->numSeqs] = i + 1; //1 extra space for 0
			affixArray->multiSeq->numSeqs++;
		}
		do {
			if(temp != '\r' && temp != '\n') {
				affixArray->multiSeq->sequences[affixArray->length++] = temp;
				if(affixArray->length % BUFFER1 == 0){
					if((affixArray->multiSeq->sequences = realloc(affixArray->multiSeq->sequences, (affixArray->length+BUFFER1)*sizeof(unsigned char))) == NULL) {
						fprintf(stderr,"Memory allocation failed for \"sequences\" - %s %d.\n", __FILE__, __LINE__);
						return 1;
					}
				}
			}
		} while((temp = getc(fp)) != EOF && temp != '>');
	} while(temp != EOF);
	affixArray->multiSeq->sequences[affixArray->length++] = $;
	affixArray->multiSeq->seqEndPos[affixArray->multiSeq->numSeqs-1] = affixArray->length - 1;

	fclose(fp);

	if (bAllocConvSeq && (affixArray->multiSeq->convSequences = (unsigned char *) calloc(affixArray->length, sizeof(unsigned char))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"convSequences\" - %s %d.\n", __FILE__, __LINE__);
		return 1;
	}

	return 0;
}

Alphabet* loadAlphabetFile(char *fileName, Alphabet *alphabet) {
	FILE *fp;
	int i, j;
	bool foundRep, isWildCard;
	char temp;

	if((fp=fopen(fileName,"r")) == NULL) {
		fprintf(stderr,"Error opening alphabet file %s.\n", fileName);
		exit(1);
	}

	if (fseek(fp, 0L, SEEK_END)) {
		fprintf(stderr,"Error reading file %s - %s %d.\n", fileName, __FILE__, __LINE__);
		exit(1);
	}

	i=ftell(fp);
	if (i == 0) {
		fprintf(stderr, "File %s is empty - %s %d.\n", fileName, __FILE__, __LINE__);
    	exit(1);
	}

	rewind(fp);

	if ((alphabet->eqClass    = (char **) calloc(127, sizeof(char*))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"eqClass\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}
	if ((alphabet->classSize  = (int *) calloc(127, sizeof(int))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"classSize\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}
	if ((alphabet->isWildCard = (bool *) calloc(127, sizeof(bool))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"isWildCard\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}
	if ((alphabet->classRepresentative = (unsigned char *) calloc(127, sizeof(unsigned char))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"classRepresentative\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	i = 0;
	temp = getc(fp);
	do {
		if(temp == '\r' || temp == '\n') {
			isWildCard = 0;
			temp = getc(fp);
		}
		else if(temp != EOF) {
			j = foundRep = 0;
			alphabet->classRepresentative[i] = ' ';
			alphabet->isWildCard[i] = 0;
			if ((alphabet->eqClass[i] = (char *) calloc(127, sizeof(char))) == NULL) {
				fprintf(stderr,"Memory allocation failed for \"eqClass_i\" - %s %d.\n", __FILE__, __LINE__);
				exit(1);
			}
			do {
				if(temp != ' ') {
					if(foundRep && j > 0)
						alphabet->classRepresentative[i] = temp;
					else if(temp == '*') {
						alphabet->isWildCard[i] = 1;
						alphabet->classRepresentative[i] = '*';
					} else
						alphabet->eqClass[i][j++] = temp;
				} else
					foundRep = 1;
			} while((temp = getc(fp)) != EOF && temp != '\r' && temp != '\n');
			if(j > 0) {
				if(!foundRep && !alphabet->isWildCard[i])
					alphabet->classRepresentative[i] = alphabet->eqClass[i][0];
				alphabet->classSize[i++] = j;
			}
		}
	} while(temp != EOF);
	alphabet->numClasses = i;

	fclose(fp);

	/*printf("numClasses=%d\n", alphabet->numClasses);
	for(i = 0; i < alphabet->numClasses; i++){
		printf("classSize[%d]=%d\n", i, alphabet->classSize[i]);
		printf("classRepresentative[%d]=~%d~\n", i, alphabet->classRepresentative[i]);
		for(j = 0; j < alphabet->classSize[i]; j++){
			printf("eqClass[%d][%d]=%c\n", i, j, alphabet->eqClass[i][j]);
		}
		printf("isWildCard[%d]=%d\n", i, alphabet->isWildCard[i]);
	}*/

	return alphabet;
}

/* Saves aff files
 * .base: numSeqs, totalLength
 * .des: seqDescLength, seqDescription
 * .seq: seqEndPos, sequences
 * .cseq:
 * .suf: suf
 * .lcp: lcp
 * .lcpe:
 * .skp:
 * */
bool saveAffFiles(AffixArray *affixArray, char *fileName, bool bSaveTransfSeq) {
	FILE *fp;
	int i;

	if(affixArray->multiSeq == NULL) {
		fprintf(stderr,"Please load a file.\n");
		return 1;
	}

	int nameLength = strlen(fileName);
	char fileNameCpy[nameLength + 7];
	strcpy(fileNameCpy, fileName);

	/********* base file *********/
	strcat(fileNameCpy, ".base");
	if((fp = fopen(fileNameCpy, "wb")) == NULL) {
	   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
	   return 1;
	}

	if(fwrite(&affixArray->length, sizeof(int), 1, fp) < 1) {
		fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		return 1;
	}
	if(fwrite(&affixArray->multiSeq->numSeqs, sizeof(int), 1, fp) < 1) {
		fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		return 1;
	}
	if(fwrite(affixArray->multiSeq->seqEndPos, sizeof(int), affixArray->multiSeq->numSeqs + 1, fp) < 1) {
		fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		return 1;
	}
	fileNameCpy[nameLength] = 0;
	fclose(fp);

	/********* alphabet file *********/
	if(affixArray->alphabet->numClasses > 0) {
		strcat(fileNameCpy, ".alph");
		if((fp = fopen(fileNameCpy, "wb")) == NULL) {
		   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		   return 1;
		}

		if(fwrite(&affixArray->alphabet->numClasses, sizeof(int), 1, fp) < 1) {
				fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
				return 1;
			}

		for(i=0; i < affixArray->alphabet->numClasses; i++) {
			if(fwrite(&affixArray->alphabet->classSize[i], sizeof(int), 1, fp) < 1) {
				fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
				return 1;
			}
			if(fwrite(&affixArray->alphabet->classRepresentative[i], sizeof(char), 1, fp) < 1) {
				fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
				return 1;
			}
			if(fwrite(affixArray->alphabet->eqClass[i], sizeof(char), affixArray->alphabet->classSize[i], fp) < 1) {
				fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
				return 1;
			}
			if(fwrite(&affixArray->alphabet->isWildCard[i], sizeof(bool), 1, fp) < 1) {
				fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
				return 1;
			}
		}
		fileNameCpy[nameLength] = 0;
		fclose(fp);
	}

	/********* description file *********/
	strcat(fileNameCpy, ".des");
	if((fp = fopen(fileNameCpy, "wb")) == NULL) {
	   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
	   return 1;
	}
	if(fwrite(affixArray->multiSeq->seqDescLength, sizeof(int), affixArray->multiSeq->numSeqs, fp) < 1) {
		fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		return 1;
	}
	for(i=0; i < affixArray->multiSeq->numSeqs; i++)
		if(fwrite(affixArray->multiSeq->seqDescription[i], sizeof(char), affixArray->multiSeq->seqDescLength[i], fp) < 1) {
			fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
			return 1;
		}
	fileNameCpy[nameLength] = 0;
	fclose(fp);

	/********* sequences file *********/
	strcat(fileNameCpy, ".seq");
	if((fp = fopen(fileNameCpy, "wb")) == NULL) {
	   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
	   return 1;
	}
	/*if(fwrite(affixArray->multiSeq->seqEndPos, sizeof(int), affixArray->multiSeq->numSeqs + 1, fp) < 1) {
		fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		return 1;
	}*/
	if(fwrite(affixArray->multiSeq->sequences, sizeof(char), affixArray->length, fp) < 1) {
		fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		return 1;
	}
	fileNameCpy[nameLength] = 0;
	fclose(fp);

	/********* alphabetically transformed sequences file *********/
	if (bSaveTransfSeq) {
		strcat(fileNameCpy, ".tseq");
		if((fp = fopen(fileNameCpy, "wb")) == NULL) {
		   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		   return 1;
		}
		if(fwrite(affixArray->multiSeq->convSequences, sizeof(char), affixArray->length, fp) < 1) {
			fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
			return 1;
		}
		fileNameCpy[nameLength] = 0;
		fclose(fp);
	}

	/*bool saveEArray(EArray* earray, int length, bool isEsa, char *fileName);
	if(saveEArray(affixArray->esa, affixArray->length, 1, fileName)) return 1;
	if(saveEArray(affixArray->erpa, affixArray->length, 0, fileName)) return 1;*/

	return 0;
}

bool saveAflk(EArray* earray, int length, bool isEsa, char *fileName) {
	FILE *fp;
	int nameLength = strlen(fileName);
	char fileNameCpy[nameLength + 7];
	strcpy(fileNameCpy, fileName);

	/********* afsuf | afrpref *********/
	if(earray->affixLink != NULL) {
		if(isEsa)
			strcat(fileNameCpy, ".aflk");
		else
			strcat(fileNameCpy, ".aflkr");
		if((fp = fopen(fileNameCpy, "wb")) == NULL) {
		   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		   return 1;
		}
		if(fwrite(earray->affixLink, sizeof(int), length, fp) < 1) {
			fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
			return 1;
		}
		fileNameCpy[nameLength] = 0;
		fclose(fp);
	}

	return 0;
}

bool saveEArray(EArray* earray, int length, bool isEsa, char *fileName) {
	FILE *fp;
	int i;

	int nameLength = strlen(fileName);
	char fileNameCpy[nameLength + 7];
	strcpy(fileNameCpy, fileName);

	/********* suf | rpref file *********/
	if(earray->xarray != NULL) {
		if(isEsa)
			strcat(fileNameCpy, ".suf");
		else
			strcat(fileNameCpy, ".sufr");
		if((fp = fopen(fileNameCpy, "wb")) == NULL) {
		   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		   return 1;
		}
		if(fwrite(earray->xarray, sizeof(int), length, fp) < 1) {
			fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
			return 1;
		}
		fileNameCpy[nameLength] = 0;
		fclose(fp);
	}

	/********* lcp | rlcp file *********/
	if(earray->xlcp != NULL) {
		if(isEsa)
			strcat(fileNameCpy, ".lcp");
		else
			strcat(fileNameCpy, ".lcpr");
		if((fp = fopen(fileNameCpy, "wb")) == NULL) {
		   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		   return 1;
		}
		if(fwrite(earray->xlcp, sizeof(char), length, fp) < 1) {
			fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
			return 1;
		}
		fileNameCpy[nameLength] = 0;
		fclose(fp);
	}

	/* Modified by Benjamin Albrecht, Jan. 2012 */
	/********* lcpTree | rlcpTree file *********/
	if (earray->xlcpTree != NULL) {
		if (isEsa)
			strcat(fileNameCpy, ".lcpTree");
		else
			strcat(fileNameCpy, ".lcprTree");
		if ((fp = fopen(fileNameCpy, "wb")) == NULL) {
			fprintf(stderr, "Error saving file \"%s\".\n", fileNameCpy);
			return 1;
		}
		if (fwrite(earray->xlcpTree, sizeof(char), length, fp) < 1) {
			fprintf(stderr, "Error saving file \"%s\".\n", fileNameCpy);
			return 1;
		}
		fileNameCpy[nameLength] = 0;
		fclose(fp);
	}

	/********* lcp | rlcp exception file *********/
	if(earray->xlcp != NULL && earray->xlcpException->numExceptions > 0) {
		if(isEsa)
			strcat(fileNameCpy, ".lcpe");
		else
			strcat(fileNameCpy, ".lcper");
		if((fp = fopen(fileNameCpy, "wb")) == NULL) {
		   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		   return 1;
		}

		if(fwrite(&earray->xlcpException->numExceptions, sizeof(int), 1, fp) < 1) {
				fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
				return 1;
			}
		for(i=0; i < earray->xlcpException->numExceptions; i++) {
			if(fwrite(&earray->xlcpException->index[i], sizeof(int), 1, fp) < 1) {
				fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
				return 1;
			}
			if(fwrite(&earray->xlcpException->value[i], sizeof(int), 1, fp) < 1) {
				fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
				return 1;
			}
		}
		fileNameCpy[nameLength] = 0;
		fclose(fp);
	}

	/* Modified by Benjamin Albrecht, Jan. 2012 */
	/********* lcpTree | rlcpTree exception file *********/
	if (earray->xlcpTree != NULL && earray->xlcpTreeException->numExceptions
			> 0) {
		if (isEsa)
			strcat(fileNameCpy, ".lcpeTree");
		else
			strcat(fileNameCpy, ".lcperTree");
		if ((fp = fopen(fileNameCpy, "wb")) == NULL) {
			fprintf(stderr, "Error saving file \"%s\".\n", fileNameCpy);
			return 1;
		}

		if (fwrite(&earray->xlcpTreeException->numExceptions, sizeof(int), 1, fp) < 1) {
			fprintf(stderr, "Error saving file \"%s\".\n", fileNameCpy);
			return 1;
		}
		for (i = 0; i < earray->xlcpTreeException->numExceptions; i++) {
			if (fwrite(&earray->xlcpTreeException->index[i], sizeof(int), 1, fp) < 1) {
				fprintf(stderr, "Error saving file \"%s\".\n", fileNameCpy);
				return 1;
			}
			if (fwrite(&earray->xlcpTreeException->value[i], sizeof(int), 1, fp) < 1) {
				fprintf(stderr, "Error saving file \"%s\".\n", fileNameCpy);
				return 1;
			}
		}
		fileNameCpy[nameLength] = 0;
		fclose(fp);
	}

	/********* skp | rskp file *********/
	if(earray->xskp != NULL) {
		if(isEsa)
			strcat(fileNameCpy, ".skp");
		else
			strcat(fileNameCpy, ".skpr");
		if((fp = fopen(fileNameCpy, "wb")) == NULL) {
		   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		   return 1;
		}
		if(fwrite(earray->xskp, sizeof(int), length, fp) < 1) {
			fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
			return 1;
		}
		fileNameCpy[nameLength] = 0;
		fclose(fp);
	}

	/********* afsuf | afrpref *********/
	if(earray->affixLink != NULL) {
		if(isEsa)
			strcat(fileNameCpy, ".aflk");
		else
			strcat(fileNameCpy, ".aflkr");
		if((fp = fopen(fileNameCpy, "wb")) == NULL) {
		   fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
		   return 1;
		}
		if(fwrite(earray->affixLink, sizeof(int), length, fp) < 1) {
			fprintf(stderr,"Error saving file \"%s\".\n", fileNameCpy);
			return 1;
		}
		fileNameCpy[nameLength] = 0;
		fclose(fp);
	}

	return 0;
}

/* Opens aff files
 * .base: numSeqs, totalLength
 * .des: seqDescLength, seqDescription
 * .seq: seqEndPos, sequences
 * .cseq:
 * .suf: suf
 * .lcp: lcp
 * .lcpe:
 * .skp:
 * */
bool loadAffFiles(AffixArray *affixArray, char *fileName,
					bool argSuf, bool argLcp, bool argLcpTree, bool argSkp, bool argAfsuf,
					bool argRpref, bool argRlcp, bool argRlcpTree, bool argRskp, bool argAfrpref) {
	FILE *fp;
	int i, fd;

	int nameLength = strlen(fileName);
	char  fileNameCpy[100];
	strcpy(fileNameCpy, fileName);

	if(affixArray->multiSeq == NULL) {
		fprintf(stderr, "Error. \"multiSeq\" was not initialized - %s %d.\n", __FILE__, __LINE__);
		return 1;
	}

	/********* base file *********/
	strcat(fileNameCpy, ".base");
	if((fp = fopen(fileNameCpy, "rb")) == NULL) {
	   fprintf(stderr,"Error opening file %s.\n", fileNameCpy);
	   return 1;
	}

	if(!fread(&affixArray->length, sizeof(int), 1, fp) != 0) {
		fprintf(stderr,"Error reading file %s.\n", fileNameCpy);
		return 1;
	}
	if(!fread(&affixArray->multiSeq->numSeqs, sizeof(int), 1, fp) != 0) {
		fprintf(stderr,"Error reading file %s.\n", fileNameCpy);
		return 1;
	}

	if((affixArray->multiSeq->seqEndPos = (int *) calloc(affixArray->multiSeq->numSeqs + 1, sizeof(int))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"seqEndPos\" - %s %d.\n", __FILE__, __LINE__);
			return 1;
		}
	if(!fread(affixArray->multiSeq->seqEndPos, sizeof(int), affixArray->multiSeq->numSeqs + 1, fp) != 0) {
		fprintf(stderr,"Error reading file %s.\n", fileNameCpy);
		return 1;
	}

	fileNameCpy[nameLength] = 0;
	fclose(fp);

	/********* alphabet file *********/
	strcat(fileNameCpy, ".alph");
	if((fp = fopen(fileNameCpy, "rb")) == NULL) {
		fprintf(stderr,"Error opening file %s.\n", fileNameCpy);
		return 1;
	}

	if (affixArray->alphabet == NULL) {
		affixArray->alphabet = (Alphabet *) malloc(sizeof(Alphabet));
	}

	if(!fread(&affixArray->alphabet->numClasses, sizeof(int), 1, fp) != 0) {
		fprintf(stderr,"Error reading file %s - %s %d.\n", fileNameCpy, __FILE__, __LINE__);
		return 1;
	}

	if((affixArray->alphabet->classSize = (int *) calloc(affixArray->alphabet->numClasses, sizeof(int))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"alphabet->classSize\".\n");
	if((affixArray->alphabet->classRepresentative = (unsigned char *) calloc(affixArray->alphabet->numClasses, sizeof(unsigned char))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"alphabet->classRepresentative\".\n");
	if((affixArray->alphabet->eqClass = (char **) calloc(affixArray->alphabet->numClasses, sizeof(char*))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"alphabet->eqClass\".\n");
	if((affixArray->alphabet->isWildCard = (bool *) calloc(affixArray->alphabet->numClasses, sizeof(bool))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"alphabet->isWildCard\".\n");

	for(i = 0; i < affixArray->alphabet->numClasses; i++) {
		if(!fread(&affixArray->alphabet->classSize[i], sizeof(int), 1, fp) != 0) {
			fprintf(stderr,"Error reading file %s - %s %d.\n", fileNameCpy, __FILE__, __LINE__);
			return 1;
		}
		if(!fread(&affixArray->alphabet->classRepresentative[i], sizeof(char), 1, fp) != 0) {
			fprintf(stderr,"Error reading file %s - %s %d.\n", fileNameCpy, __FILE__, __LINE__);
			return 1;
		}

		if((affixArray->alphabet->eqClass[i] = (char *) calloc(affixArray->alphabet->classSize[i], sizeof(char*))) == NULL)
			fprintf(stderr,"Memory allocation failed for \"alphabet->eqClass\".\n");
		if(!fread(affixArray->alphabet->eqClass[i], sizeof(char), affixArray->alphabet->classSize[i], fp) != 0) {
			fprintf(stderr,"Error reading file %s - %s %d.\n", fileNameCpy, __FILE__, __LINE__);
			return 1;
		}
		if(!fread(&affixArray->alphabet->isWildCard[i], sizeof(bool), 1, fp) != 0) {
			fprintf(stderr,"Error reading file %s - %s %d.\n", fileNameCpy, __FILE__, __LINE__);
			return 1;
		}
	}
	fileNameCpy[nameLength] = 0;
	fclose(fp);

	setIupacTable(affixArray->alphabet);

	/********* description file *********/
	strcat(fileNameCpy, ".des");
	if((fp = fopen(fileNameCpy, "rb")) == NULL) {
	   fprintf(stderr,"Error opening file %s.\n", fileNameCpy);
	   return 1;
	}

	if((affixArray->multiSeq->seqDescLength = (int *) calloc(affixArray->multiSeq->numSeqs, sizeof(int))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"seqDescLength\".\n");
	if(!fread(affixArray->multiSeq->seqDescLength, sizeof(int), affixArray->multiSeq->numSeqs, fp) != 0) {
		fprintf(stderr,"Error reading file %s - %s %d.\n", fileNameCpy, __FILE__, __LINE__);
		return 1;
	}

	if((affixArray->multiSeq->seqDescription = (char **) calloc(affixArray->multiSeq->numSeqs, sizeof(char*))) == NULL)
		fprintf(stderr,"Memory allocation failed for \"seqDescription\".\n");
	for(i=0; i < affixArray->multiSeq->numSeqs; i++){
		if((affixArray->multiSeq->seqDescription[i] = (char *) calloc(affixArray->multiSeq->seqDescLength[i], sizeof(char))) == NULL)
			fprintf(stderr,"Memory allocation failed for \"seqDescription\".\n");
		if(!fread(affixArray->multiSeq->seqDescription[i], sizeof(char), affixArray->multiSeq->seqDescLength[i], fp) != 0) {
			fprintf(stderr,"Error reading file %s - %s %d.\n", fileNameCpy, __FILE__, __LINE__);
			return 1;
		}
	}

	fileNameCpy[nameLength] = 0;
	fclose(fp);

	/********* sequences file *********/
	strcat(fileNameCpy, ".seq");
	if((fd = open(fileNameCpy, O_RDONLY)) == -1) {
		fprintf(stderr,"Error reading file %s - %s %d.\n", fileNameCpy, __FILE__, __LINE__);
		return 1;
	}
	affixArray->multiSeq->sequences = mmap(NULL, affixArray->length * sizeof(unsigned char), PROT_READ, MAP_SHARED, fd, 0);
	affixArray->multiSeq->sequencesMmapped = true;

	fileNameCpy[nameLength] = 0;
	//fclose(fd);

	/********* alphabetically transformed sequences file *********/
	strcat(fileNameCpy, ".tseq");
	if((fd = open(fileNameCpy, O_RDONLY)) != -1) {
		affixArray->multiSeq->convSequences = mmap(NULL, affixArray->length * sizeof(unsigned char), PROT_READ, MAP_SHARED, fd, 0);
		affixArray->multiSeq->convSequencesMmapped = true;

		fileNameCpy[nameLength] = 0;
		//fclose(fd);
	} else {
		affixArray->multiSeq->convSequences = NULL;
		affixArray->multiSeq->convSequencesMmapped = false;
	}

	/* Modified by Benjamin Albrecht, Jan. 2012 */
	bool loadEArray(EArray *earray, int length, bool isEsa, char *fileName,
					bool argXarray, bool argXlcp, bool argXlcpTree, bool argXskp, bool argAfx);
	if(loadEArray(affixArray->esa, affixArray->length, 1, fileName,
					argSuf, argLcp, argLcpTree, argSkp, argAfsuf)) return 1;
	if(loadEArray(affixArray->erpa, affixArray->length, 0, fileName,
    				argRpref, argRlcp, argRlcpTree, argRskp, argAfrpref)) return 1;

	return 0;
}

/* Modified by Benjamin Albrecht, Jan. 2012 */
bool loadEArray(EArray *earray, int length, bool isEsa, char *fileName,
				bool argXarray, bool argXlcp, bool argXlcpTree, bool argXskp, bool argAfx) {
	int i, j, fd, *dataInt32;

	int nameLength = strlen(fileName);
	char  fileNameCpy[100];
	strcpy(fileNameCpy, fileName);

	/********* suf | rpref file *********/
	if(argXarray) {
		if(isEsa)
			strcat(fileNameCpy, ".suf");
		else
			strcat(fileNameCpy, ".sufr");
		if((fd = open(fileNameCpy, O_RDONLY)) == -1) {
		   fprintf(stderr,"Error opening file \"%s\".\n", fileNameCpy);
		   return 1;
		}

	    earray->xarray = mmap(NULL, length * sizeof(int), PROT_READ, MAP_SHARED, fd, 0);

		fileNameCpy[nameLength] = 0;
		//close(fd);
	}

	/********* lcp | rlcp file *********/
	if(argXlcp) {
		if(isEsa)
			strcat(fileNameCpy, ".lcp");
		else
			strcat(fileNameCpy, ".lcpr");
		if((fd = open(fileNameCpy, O_RDONLY)) == -1) {
		   fprintf(stderr,"Error opening file \"%s\".\n", fileNameCpy);
		   return 1;
		}

	    earray->xlcp = mmap(NULL, length * sizeof(unsigned char), PROT_READ, MAP_SHARED, fd, 0);

		fileNameCpy[nameLength] = 0;
		//close(fd);

		/********* lcp | rlcpe exception file *********/

		if(isEsa)
			strcat(fileNameCpy, ".lcpe");
		else
			strcat(fileNameCpy, ".lcper");
		if((fd = open(fileNameCpy, O_RDONLY)) == -1) {
			earray->xlcpException = (LcpException *) malloc(sizeof(LcpException));
			earray->xlcpException->numExceptions = 0;
		} else {
			dataInt32 = mmap(NULL, sizeof(int), PROT_READ, MAP_SHARED, fd, 0);
			if((earray->xlcpException = (LcpException *) malloc(sizeof(LcpException))) == NULL)
				fprintf(stderr,"Memory allocation failed for \"lcpException\".\n");
			earray->xlcpException->numExceptions = dataInt32[0];

			dataInt32 = mmap(NULL, earray->xlcpException->numExceptions * sizeof(int) * 2 + sizeof(int), PROT_READ, MAP_SHARED, fd, 0);
			if((earray->xlcpException->index = (int *) calloc(earray->xlcpException->numExceptions, sizeof(int))) == NULL)
				fprintf(stderr,"Memory allocation failed for \"lcpException->index\".\n");
			if((earray->xlcpException->value = (int *) calloc(earray->xlcpException->numExceptions, sizeof(int))) == NULL)
				fprintf(stderr,"Memory allocation failed for \"lcpException->value\".\n");

			j = 1;
			for(i=0; i < earray->xlcpException->numExceptions; i++) {
				earray->xlcpException->index[i] = dataInt32[j++];
				earray->xlcpException->value[i] = dataInt32[j++];
			}
			//close(fd);
		}
		fileNameCpy[nameLength] = 0;

		/*for(i=0; i < earray->xlcpException->numExceptions; i++) {
			printf("%d ==> %d	", i, earray->xlcpException->index[i]);
			printf("%d\n", earray->xlcpException->value[i]);
		}
		printf("\n");*/

	}

	/* Modified by Benjamin Albrecht, Jan. 2012 */
	/********* lcpTree | lcprTree file *********/
	if (argXlcpTree) {
		if (isEsa)
			strcat(fileNameCpy, ".lcpTree");
		else
			strcat(fileNameCpy, ".lcprTree");
		if ((fd = open(fileNameCpy, O_RDONLY)) == -1) {
			fprintf(stderr, "Error opening file \"%s\".\n", fileNameCpy);
			return 1;
		}

		earray->xlcpTree = mmap(NULL, length * sizeof(unsigned char),
				PROT_READ, MAP_SHARED, fd, 0);

		fileNameCpy[nameLength] = 0;
		//close(fd);

		/********* lcperTree | lcpeTree exception file *********/

		if (isEsa)
			strcat(fileNameCpy, ".lcpeTree");
		else
			strcat(fileNameCpy, ".lcperTree");
		if ((fd = open(fileNameCpy, O_RDONLY)) == -1) {
			earray->xlcpTreeException = (LcpException *) malloc(
					sizeof(LcpException));
			earray->xlcpTreeException->numExceptions = 0;
		} else {
			dataInt32 = mmap(NULL, sizeof(int), PROT_READ, MAP_SHARED, fd, 0);
			if ((earray->xlcpTreeException = (LcpException *) malloc(sizeof(LcpException))) == NULL)
				fprintf(stderr, "Memory allocation failed for \"lcpTreeException\".\n");
			earray->xlcpTreeException->numExceptions = dataInt32[0];

			dataInt32 = mmap(NULL, earray->xlcpTreeException->numExceptions* sizeof(int) * 2 + sizeof(int), PROT_READ, MAP_SHARED, fd, 0);
			if ((earray->xlcpTreeException->index = (int *) calloc(earray->xlcpTreeException->numExceptions, sizeof(int))) == NULL)
				fprintf(stderr,
						"Memory allocation failed for \"xlcpTreeException->index\".\n");
			if ((earray->xlcpTreeException->value = (int *) calloc(earray->xlcpTreeException->numExceptions, sizeof(int))) == NULL)
				fprintf(stderr,"Memory allocation failed for \"xlcpTreeException->value\".\n");

			j = 1;
			for (i = 0; i < earray->xlcpTreeException->numExceptions; i++) {
				earray->xlcpTreeException->index[i] = dataInt32[j++];
				earray->xlcpTreeException->value[i] = dataInt32[j++];
			}
			//close(fd);
		}
		fileNameCpy[nameLength] = 0;

		/*for(i=0; i < earray->xlcpException->numExceptions; i++) {
		 printf("%d ==> %d	", i, earray->xlcpException->index[i]);
		 printf("%d\n", earray->xlcpException->value[i]);
		 }
		 printf("\n");*/

	}

	/********* skp | rskp file *********/
	if(argXskp) {
		if(isEsa)
			strcat(fileNameCpy, ".skp");
		else
			strcat(fileNameCpy, ".skpr");
		if((fd = open(fileNameCpy, O_RDONLY)) == -1) {
		   fprintf(stderr,"Error opening file \"%s\".\n", fileNameCpy);
		   return 1;
		}

	    earray->xskp = mmap(NULL, length * sizeof(int), PROT_READ, MAP_SHARED, fd, 0);

		fileNameCpy[nameLength] = 0;
		//close(fd);
	}

	/********* afsuf | afrpref file *********/
	if(argAfx) {
		if(isEsa)
			strcat(fileNameCpy, ".aflk");
		else
			strcat(fileNameCpy, ".aflkr");
		if((fd = open(fileNameCpy, O_RDONLY)) == -1) {
		   fprintf(stderr,"Error opening file \"%s\".\n", fileNameCpy);
		   return 1;
		}

	    earray->affixLink = mmap(NULL, length * sizeof(int), PROT_READ, MAP_SHARED, fd, 0);

		//fileNameCpy[nameLength] = 0;
		//close(fd);
	}

	return 0;
}

void allocConvSequences (AffixArray *affixArray) {
	if ((affixArray->multiSeq->convSequences = (unsigned char *) calloc(affixArray->length, sizeof(unsigned char))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"convSequences\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}
}

bool printEArray(AffixArray *affixArray, bool isEsa) {
	int i;
	unsigned char charValue;
	unsigned char *seqTmp;
	EArray *earray;

	if(affixArray->length == 0) {
		fprintf(stderr,"Please load a file.\n");
		return 1;
	}

	if(isEsa) {
		earray = affixArray->esa;
		if(earray->xarray == NULL) return 0;

		printf("i	suf[i]	");
		if(earray->xlcp != NULL)
			printf("lcp[i]	");
		if(earray->xskp != NULL)
			printf("skp[i]	");
		if(earray->affixLink != NULL)
		 	printf("aflk[i]");
		 printf("Ssuf[i]\n");
	} else {
		earray = affixArray->erpa;
		if(earray->xarray == NULL) return 0;

		printf("i	rpref[i]");
		if(earray->xlcp != NULL)
			printf("rlcp[i]");
		if(earray->xskp != NULL)
			printf("rskp[i]	");
		if(earray->affixLink != NULL)
		 	printf("afrlk[i]");
		 printf("Srpref[i]\n");

		 reverseSequences(affixArray->multiSeq->convSequences, affixArray->length);
	}

	for(i = 0; i < affixArray->length; i++) {
		printf("%d	%d	",i, earray->xarray[i]);
		if(earray->xlcp != NULL)
			printf("%d	", xlcpvalue(i));
		if(earray->xskp != NULL)
			printf("%d	", earray->xskp[i]);
		if(earray->affixLink != NULL)
			earray->affixLink[i] != INT_MAX ? printf("%d	", earray->affixLink[i]) : printf("	");

		seqTmp = affixArray->multiSeq->convSequences + earray->xarray[i];
		do {
			charValue = *seqTmp++;
			if(charValue == 127)
				printf("$");
			else
				printf("%c", affixArray->alphabet->classRepresentative[charValue - 1]);
		} while (charValue != $);
		printf("\n");
	}

	if(!isEsa)
		reverseSequences(affixArray->multiSeq->convSequences, affixArray->length);

	return 0;
}

bool saveTextFile(AffixArray* affixArray, char *fileName, bool isEsa) {
	FILE *fp;
	EArray *earray;
	int i;
	unsigned char charValue;
	unsigned char *seqTmp;

	int nameLength = strlen(fileName);
	char fileNameCpy[nameLength + 7];
	strcpy(fileNameCpy, fileName);

	if(affixArray->length == 0) {
		fprintf(stderr,"Please load a file.\n");
		return 1;
	}

	if(isEsa)
		strcat(fileNameCpy, ".esa");
	else
		strcat(fileNameCpy, ".erpa");

	if((fp = fopen(fileNameCpy, "w")) == NULL) {
	   fprintf(stderr,"Error saving file \"%s\"\n", fileName);
	   return 1;
	}

	if(isEsa) {
		earray = affixArray->esa;
		if(earray->xarray == NULL) return 0;

		fprintf(fp, "i	suf[i]	");
		if(earray->xlcp != NULL)
			fprintf(fp, "lcp[i]	");
		if(earray->xskp != NULL)
			fprintf(fp, "skp[i]	");
		if(earray->affixLink != NULL)
		 	fprintf(fp, "aflk[i]");
		 fprintf(fp, "Ssuf[i]\n");
	} else {
		earray = affixArray->erpa;
		if(earray->xarray == NULL) return 0;

		fprintf(fp, "i	rpref[i]");
		if(earray->xlcp != NULL)
			fprintf(fp, "rlcp[i]");
		if(earray->xskp != NULL)
			fprintf(fp, "rskp[i]	");
		if(earray->affixLink != NULL)
		 	fprintf(fp, "afrlk[i]");
		 fprintf(fp, "Srpref[i]\n");

		 reverseSequences(affixArray->multiSeq->convSequences, affixArray->length);
	}

	for(i = 0; i < affixArray->length; i++) {
		fprintf(fp, "%d	%d	",i, earray->xarray[i]);
		if(earray->xlcp != NULL)
			fprintf(fp, "%d	", xlcpvalue(i));
		if(earray->xskp != NULL)
			fprintf(fp, "%d	", earray->xskp[i]);
		if(earray->affixLink != NULL)
			earray->affixLink[i] != INT_MAX ? fprintf(fp, "%d	", earray->affixLink[i]) : fprintf(fp, "	");

		seqTmp = affixArray->multiSeq->convSequences + earray->xarray[i];
		do {
			charValue = *seqTmp++;
			if(charValue == 127)
				fprintf(fp, "$");
			else
				fprintf(fp, "%c", affixArray->alphabet->classRepresentative[charValue - 1]);
		} while (charValue != $);
		fprintf(fp, "\n");
	}

	if(!isEsa)
		reverseSequences(affixArray->multiSeq->convSequences, affixArray->length);

	fclose(fp);
	return 0;
}

#define swap(x, a, b) { temp = x[a]; \
                     x[a] = x[b]; x[b] = temp; }

void sortSearchResult(int *index, int left, int right) {
	int i, last, temp;
	if(left >= right)
		return;
	//i = (left + right)/2;
	i = (left + right) >> 1;
	swap(index, left, i);
	last = left;
	for(i = left+1; i <= right; i++){
		if(index[i] < index[left]) {
			++last;
			swap(index, last, i);
		}
	}
	swap(index, left, last);
	sortSearchResult(index, left, last-1);
	sortSearchResult(index, last+1, right);
}

int loadPatternFileBidSearchDebug(char *fileName, StringExtension *stringExtension, AffixArray *affixArray) {
	char *p;
	int i, j, pos, fd;
	struct stat sb;

	if((fd = open(fileName, O_RDONLY)) == -1) {
	   fprintf(stderr,"Error opening file \"%s\".\n", fileName);
	   exit(1);
	}
	if (fstat(fd, &sb) == -1) {
		fprintf(stderr,"fstat\n");
		exit(1);
	}
	p = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

	pos = i = 0;
	do {
		if(p[pos] == '>')
			stringExtension[i].toRight = 1;
		else if(p[pos] == '<')
			stringExtension[i].toRight = 0;
		else {
			fprintf(stderr, "File format is invalid.\n");
			exit(1);
		}

		if(++pos >= sb.st_size) return i;

		j = 0;
		do {
			if(p[pos] != '\r' && p[pos] != '\n') {
				stringExtension[i].extension[j++] = p[pos];
			}
		} while(++pos < sb.st_size && p[pos] != '>' && p[pos] != '<');

		if(convertToAlphabet(stringExtension[i].extension,
							stringExtension[i].extension,
							j,
							true,
							affixArray->alphabet))
			exit(1);

		i++;
	} while(pos < sb.st_size);

	return i;
}

MultiPattern* loadPatternFile(char *fileName, AffixArray *affixArray) {
	MultiPattern *multiPattern;
	char *p, cNumber[12], maxLength=12, parameters[100], *pParameters;
	int iWhichPipeSymbol;

	int i, j, k, x, pos, fd;
	int iBrackets = 0;
	Pattern *patternP;
	bool hasBrackets, bParams;
	struct stat sb;

	if((fd = open(fileName, O_RDONLY)) == -1) {
	   fprintf(stderr,"Error opening file \"%s\".\n", fileName);
	   exit(1);
	}
	if (fstat(fd, &sb) == -1) {
		fprintf(stderr,"fstat\n");
		exit(1);
	}
	p = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

	if((multiPattern = (MultiPattern *) malloc(sizeof(MultiPattern))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"multiPattern\".\n");
		exit(1);
	}
	multiPattern->usesInstanceParam = 0;
	multiPattern->usesStartPosParam = 0;
	multiPattern->usesWeightParam   = 0;

	if((patternP = (Pattern *) calloc(BUFFER2, sizeof(Pattern))) == NULL){
		fprintf(stderr,"Memory allocation failed for \"patternP\".\n");
		exit(1);
	}

	pos = i = 0;
	do {
		if(p[pos] != '>') {
			fprintf(stderr, "Format of \"%s\" is invalid. \">\" was expected.\n", fileName);
			exit(1);
		}

		/******************* Read description *******************/
		if((patternP[i].desc = (char *) calloc(BUFFER2, sizeof(char))) == NULL){
			fprintf(stderr,"Memory allocation failed for \"patternP[%d].desc\".\n",i);
			exit(1);
		}

		for (j = 0; j < 100; ++j) {
			parameters[j] = 0;
		}
		for (j = 0; j < maxLength; ++j) {
			cNumber[j] = 0;
		}

		k = j = iWhichPipeSymbol = 0;
		bParams = 0;
		while(++pos < sb.st_size && p[pos] != '\r' && p[pos] != '\n') {
			if(p[pos] == '|') {
				bParams = 1;
			}

			if (!bParams) {
				patternP[i].desc[j++] = p[pos];

				if((j + 1) % BUFFER2 == 0 && (patternP[i].desc = realloc(patternP[i].desc, (j + BUFFER2 + 1) * sizeof(char))) == NULL) {
					fprintf(stderr,"Memory allocation failed for \"patternP[%d].desc\".\n", i);
					exit(1);
				}
			} else {
				parameters[k++] = p[pos];
			}
		}

		if((patternP[i].desc = realloc(patternP[i].desc, (j + 1) * sizeof(char))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"patternP[%d].desc\".\n", i);
			exit(1);
		}

		patternP[i].weight = 1;

		if (k > 0) {
			if ((pParameters = strstr(parameters, "weight")) != NULL) {
				multiPattern->usesWeightParam = 1;
				k = 0;
				while (pParameters[0] != '|' && pParameters[0] != '\0') {
					if (pParameters[0] > 47 && pParameters[0] < 58 && k < maxLength) { //if it is a number
						cNumber[k++] = pParameters[0];
					}
					pParameters++;
				}
				if (k > 0) {
					patternP[i].weight = atoi(cNumber);
					for (x = 0; x <= k; x++) {
						cNumber[x] = 0;
					}
				}
			}
			// repeated code
			if ((pParameters = strstr(parameters, "startpos")) != NULL) {
				multiPattern->usesStartPosParam = 1;
				k = 0;
				while (pParameters[0] != '|' && pParameters[0] != '\0') {
					if (pParameters[0] > 47 && pParameters[0] < 58 && k < maxLength) { //if it is a number
						cNumber[k++] = pParameters[0];
					}
					pParameters++;
				}
				if (k > 0) {
					patternP[i].startpos = atoi(cNumber);
					for (x = 0; x <= k; x++) {
						cNumber[x] = 0;
					}
				}
			}
			if ((pParameters = strstr(parameters, "instance")) != NULL) {
				multiPattern->usesInstanceParam = 1;
				k = 0;
				while (pParameters[0] != '|' && pParameters[0] != '\0') {
					if (pParameters[0] > 47 && pParameters[0] < 58 && k < maxLength) { //if it is a number
						cNumber[k++] = pParameters[0];
					}
					pParameters++;
				}
				if (k > 0) {
					patternP[i].instance = atoi(cNumber);
					for (x = 0; x <= k; x++) {
						cNumber[x] = 0;
					}
				}
			}
			if ((pParameters = strstr(parameters, "maxstemlength")) != NULL || (pParameters = strstr(parameters, "msl")) != NULL) {
				k = 0;
				while (pParameters[0] != '|' && pParameters[0] != '\0') {
					if (pParameters[0] > 47 && pParameters[0] < 58 && k < maxLength) { //if it is a number
						cNumber[k++] = pParameters[0];
					}
					pParameters++;
				}
				if (k > 0) {
					patternP[i].maxstemlength = atoi(cNumber);
					for (x = 0; x <= k; x++) {
						cNumber[x] = 0;
					}
				}
			} else {
				patternP[i].maxstemlength = 0;
			}
			if ((pParameters = strstr(parameters, "maxrightloopextent")) != NULL || (pParameters = strstr(parameters, "mrlex")) != NULL) {
				k = 0;
				while (pParameters[0] != '|' && pParameters[0] != '\0') {
					if (pParameters[0] > 47 && pParameters[0] < 58 && k < maxLength) { //if it is a number
						cNumber[k++] = pParameters[0];
					}
					pParameters++;
				}
				if (k > 0) {
					patternP[i].maxrightloopextent = atoi(cNumber);
					for (x = 0; x <= k; x++) {
						cNumber[x] = 0;
					}
				}
			} else {
				patternP[i].maxrightloopextent = 0;
			}
			if ((pParameters = strstr(parameters, "maxleftloopextent")) != NULL || (pParameters = strstr(parameters, "mllex")) != NULL) {
				k = 0;
				while (pParameters[0] != '|' && pParameters[0] != '\0') {
					if (pParameters[0] > 47 && pParameters[0] < 58 && k < maxLength) { //if it is a number
						cNumber[k++] = pParameters[0];
					}
					pParameters++;
				}
				if (k > 0) {
					patternP[i].maxleftloopextent = atoi(cNumber);
					for (x = 0; x <= k; x++) {
						cNumber[x] = 0;
					}
				}
			} else {
				patternP[i].maxleftloopextent = 0;
			}
			if ((pParameters = strstr(parameters, "maxmispair")) != NULL) {
				k = 0;
				while (pParameters[0] != '|' && pParameters[0] != '\0') {
					if (pParameters[0] > 47 && pParameters[0] < 58 && k < maxLength) { //if it is a number
						cNumber[k++] = pParameters[0];
					}
					pParameters++;
				}
				if (k > 0) {
					patternP[i].maxmispair = atoi(cNumber);
					for (x = 0; x <= k; x++) {
						cNumber[x] = 0;
					}
				}
			} else {
				patternP[i].maxmispair = 0;
			}
		}

		/******************* Read sequence *******************/
		if((patternP[i].seq = (unsigned char *) calloc(BUFFER2, sizeof(unsigned char))) == NULL){
			fprintf(stderr,"Memory allocation failed for \"patternP[%d].seq\".\n", i);
			exit(1);
		}
		j = 0;
		while(++pos < sb.st_size && (!(p[pos] > 39 && p[pos] < 63) || p[pos] == 42)) {
			if(p[pos] != '\r' && p[pos] != '\n') {
				patternP[i].seq[j++] = p[pos];
				if((j + 1) % BUFFER2 == 0 && (patternP[i].seq = realloc(patternP[i].seq, (j + BUFFER2 + 1) * sizeof(char))) == NULL) {
					fprintf(stderr,"Memory allocation failed for \"patternP[%d].seq\".\n", i);
					exit(1);
				}
			}
		}
		if((patternP[i].seq = realloc(patternP[i].seq, (j + 1) * sizeof(char))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"patternP[%d].seq\".\n", i);
			exit(1);
		}

		if(convertToAlphabet(patternP[i].seq,
							patternP[i].seq,
							j,
							true,
							affixArray->alphabet))
			exit(1);

		/******************* Read structure *******************/
		if((patternP[i].structure = (unsigned char *) calloc(BUFFER2, sizeof(unsigned char))) == NULL){
			fprintf(stderr,"Memory allocation failed for \"patternP[%d].structure\".\n", i);
			exit(1);
		}
		j = 0;
		hasBrackets = 0;
		while(pos < sb.st_size && !(p[pos] > 59 && p[pos] < 123)) {
			if(p[pos] != '\r' && p[pos] != '\n') {
				if(p[pos] == '(') {
					iBrackets++;
					hasBrackets = 1;
				} else if(p[pos] == ')')
					iBrackets--;
				if(iBrackets < 0) {
					fprintf(stderr,"Structure of pattern \"%s\" is invalid.\n", patternP[i].desc);
					exit(1);
				}
				patternP[i].structure[j++] = p[pos];
				if((j + 1) % BUFFER2 == 0 && (patternP[i].structure = realloc(patternP[i].structure, (j + BUFFER2 + 1) * sizeof(unsigned char))) == NULL) {
					fprintf(stderr,"Memory allocation failed for \"patternP[%d].structure\".\n", i);
					exit(1);
				}
			}
			pos++;
		}
		if(iBrackets != 0) {
			fprintf(stderr,"Structure of pattern \"%s\" is invalid.\n", patternP[i].desc);
			exit(1);
		}

		if(strlen((char*) patternP[i].structure) != strlen((char*) patternP[i].seq)) {
			fprintf(stderr,"Sequence and structure length of pattern \"%s\" do not match.\n", patternP[i].desc);
			exit(1);
		}
		if(hasBrackets && (patternP[i].structure[0] == '.' || patternP[i].structure[strlen((char*) patternP[i].structure) - 1] == '.')) {
			fprintf(stderr,"This version of afsearch does not support hairpins with dangling ends.\n");
			exit(1);
		}
		if((patternP[i].structure = realloc(patternP[i].structure, (j + 1) * sizeof(char))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"patternP[%d].structure\".\n", i);
			exit(1);
		}

		if(((++i)+1) % BUFFER2 == 0 && (patternP = realloc(patternP, (i + BUFFER2 + 1) * sizeof(Pattern))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"patternP\".\n");
			exit(1);
		}
	} while(pos < sb.st_size);

	multiPattern->numPatterns = i;
	multiPattern->pattern     = patternP;

	return multiPattern;
}

void printSearchResult(AffixArray* affixArray, char *pattern, int *pos, int length) {
	int i, seqNum=0, seqNumPrev, offset=0;
	bool flag=0;

	convertToRepresentative(pattern, strlen(pattern), affixArray->alphabet);
	if(length == 0)
		printf("Pattern \"%s\" was not found in the sequence(s).\n", pattern);
	else {
		printf("Pattern \"%s\" was found %d time(s) in the following sequence(s) in position(s):\n", pattern, length);

		sortSearchResult(pos, 0, length -1);

		seqNumPrev = INT_MAX;
		for(i = 0; i < length; i++) {
			if(!flag || pos[i] >= offset)
				seqNum = getSeqNumber(affixArray->multiSeq->seqEndPos, affixArray->multiSeq->numSeqs, pos[i]);

			if(seqNum > 0)
				offset = affixArray->multiSeq->seqEndPos[seqNum - 1] + 1;
			else
				offset = 0;

			if(seqNumPrev != seqNum) {
				if(flag)
					printf("\n");
				flag = 1;
				printf("%s: ", affixArray->multiSeq->seqDescription[seqNum]);
			}
			printf("%d ", pos[i] - offset);
			seqNumPrev = seqNum;
		}
		printf("\n");

		if(affixArray->multiSeq->numSeqs > 1) {
			printf("Pattern position(s) in the concatenated sequences array: ");
			pos = getRelativePos(pos, length, affixArray);

			for(i = 0; i < length; i++)
				printf("%d ",pos[i]);
			printf("\n");
		}
	}

}

/* This function was used until the introduction of the
 * RedBlack tree */
void saveSearchResults(Pattern *patternP, int *matchPos, int iMatches, char *fileName, bool overwrite, bool needsSorting, AffixArray *affixArray) {
	int i, j, k, seqNum=0, seqNumPrev, offset=0, pLength;
	bool flag=0;
	FILE *fp=NULL;

	if(fileName != NULL && (fp = fopen(fileName, overwrite ? "w" : "a")) == NULL) {
	   fprintf(stderr,"Error saving file %s.\n", fileName);
	   exit(1);
	}

	if(iMatches > 1 && needsSorting)
		sortSearchResult(matchPos, 0, iMatches -1);

	fprintf(fp, "#%s\n", patternP->desc);

	pLength = strlen((char*) patternP->seq);

	seqNumPrev = INT_MAX;
	for(i = 0; i < iMatches; i++) {
		if(!flag || matchPos[i] >= offset)
			seqNum = getSeqNumber(affixArray->multiSeq->seqEndPos, affixArray->multiSeq->numSeqs, matchPos[i]);

		if(seqNum > 0)
			offset = affixArray->multiSeq->seqEndPos[seqNum - 1] + 1;
		else
			offset = 0;

		if(seqNumPrev != seqNum) {
			flag = 1;
			fprintf(fp, ">%s\n", affixArray->multiSeq->seqDescription[seqNum]);
		}

		k = matchPos[i] + pLength;
		fprintf(fp, "%d	", matchPos[i] - offset);
		for (j = matchPos[i]; j < k; j++) {
			fprintf(fp, "%c", affixArray->alphabet->classRepresentative[affixArray->multiSeq->convSequences[j] - 1]);
		}
		fprintf(fp, "\n");

		seqNumPrev = seqNum;
	}
	fprintf(fp, "\r");

	fclose(fp);
}

bool** loadComplementarityFile(char *fileName, AffixArray *affixArray) {
	unsigned char *p, *c;
	int i, j, k, pos, fd;
	struct stat sb;
	bool **bCompCheck;

	if((fd = open(fileName, O_RDONLY)) == -1) {
	   fprintf(stderr,"Error opening file \"%s\".\n", fileName);
	   exit(1);
	}
	if (fstat(fd, &sb) == -1) {
		fprintf(stderr,"fstat\n");
		exit(1);
	}
	p = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

	if ((bCompCheck = (bool **) calloc(affixArray->alphabet->numClasses, sizeof(bool*))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"bCompCheck\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	for (i = 0; i < affixArray->alphabet->numClasses; i++) {
		if ((bCompCheck[i] = (bool *) calloc(affixArray->alphabet->numClasses, sizeof(bool))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"bCompCheck_i\" - %s %d.\n", __FILE__, __LINE__);
			exit(1);
		}
		for (j = 0; j < affixArray->alphabet->numClasses; j++) {
			bCompCheck[i][j] = 0;
		}
	}

	pos = 0;
	do {
		if(p[pos] != '\r' && p[pos] != '\n') {
			c = (unsigned char *) calloc(BUFFER2, sizeof(unsigned char));

			j = 0;
			do {
				c[j++] = p[pos++];
			} while (pos < sb.st_size && p[pos] != '\r' && p[pos] != '\n');

			if(convertToAlphabet(c, c, j, true, affixArray->alphabet))
				exit(1);

			for (k = 1; k < j; k++) {
				bCompCheck[c[0] - 1][c[k] - 1] = 1;
			}

			free(c);
		}

		pos++;
	} while(pos < sb.st_size);

	/* If base X is complement of base Y, make vice-versa true */
	for (i = 0; i < affixArray->alphabet->numClasses; i++) {
		for (j = 0; j < affixArray->alphabet->numClasses; j++) {
			if(bCompCheck[i][j]) {
				bCompCheck[j][i] = 1;
			}
		}
	}

	if (munmap(p, sb.st_size) == -1)
		fprintf(stderr,"Error unmapping file \"%s\".\n", fileName);

	close(fd);

	return bCompCheck;
}

bool** loadDefaultComplementarityRules(bool bSilent2, int predefAlphabet, AffixArray *affixArray) {
	int i, j, k, numComps=/*3*/2, length=2;
	//char c[3][2] = {{'A', 'U'}, {'C', 'G'}, {'G', 'U'}}, cc[length];
	unsigned char c[2][2] = {{'A', 'U'}, {'C', 'G'}}, cc[length];

	bool **bCompCheck = (bool **) calloc(affixArray->alphabet->numClasses, sizeof(bool*));
	if (predefAlphabet > 1)
		return bCompCheck;

	if (!bSilent2) printf("\n%cUsing Watson-Crick complementarity rules A-U, U-A, C-G, G-C\n", LINESYMBOL);

	for (i = 0; i < affixArray->alphabet->numClasses; i++) {
		bCompCheck[i] = (bool *) calloc(affixArray->alphabet->numClasses, sizeof(bool));
		for (j = 0; j < affixArray->alphabet->numClasses; j++) {
			bCompCheck[i][j] = 0;
		}
	}

	for (i = 0; i < numComps; i++) {
		if(convertToAlphabet(c[i], cc, 2, true, affixArray->alphabet))
			exit(1);

		for (k = 1; k < length; k++) {
			bCompCheck[cc[0] - 1][cc[k] - 1] = 1;
		}
	}

	/* If base X is complement of base Y, make vice-versa true */
	for (i = 0; i < affixArray->alphabet->numClasses; i++) {
		for (j = 0; j < affixArray->alphabet->numClasses; j++) {
			if(bCompCheck[i][j]) {
				bCompCheck[j][i] = 1;
			}
		}
	}

	return bCompCheck;
}

bool** loadReverseComplementarityFile(bool **bCompCheck, AffixArray *affixArray) {
	int i, j;
	bool **bRevCompCheck;
	char cA=1, cC=2, cG=3, cT=4, ci, cj;
	Alphabet *alphabet = affixArray->alphabet;

	if ((bRevCompCheck = (bool **) calloc(alphabet->numClasses, sizeof(bool*))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"bRevCompCheck\" - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	for (i = 0; i < alphabet->numClasses; i++) {
		if ((bRevCompCheck[i] = (bool *) calloc(alphabet->numClasses, sizeof(bool))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"bRevCompCheck_i\" - %s %d.\n", __FILE__, __LINE__);
			exit(1);
		}
		for (j = 0; j < alphabet->numClasses; j++) {
			bRevCompCheck[i][j] = 0;
		}
	}

	for (i = 0; i < alphabet->numClasses; i++) {
		switch (toupper(alphabet->classRepresentative[i])) {
			case 'A':
				cA = i;
				break;
			case 'C':
				cC = i;
				break;
			case 'G':
				cG = i;
				break;
			case 'T':
				cT = i;
				break;
			case 'U':
				cT = i;
				break;
		}
	}

	for (i = 0; i < alphabet->numClasses; i++) {
		if (i == cA)
			ci = cT;
		else if (i == cC)
			ci = cG;
		else if (i == cG)
			ci = cC;
		else if (i == cT)
			ci = cA;
		else
			ci = i;

		for (j = 0; j < alphabet->numClasses; j++) {
			if (bCompCheck[i][j]) {

				if (j == cA)
					cj = cT;
				else if (j == cC)
					cj = cG;
				else if (j == cG)
					cj = cC;
				else if (j == cT)
					cj = cA;
				else
					cj = j;

				bRevCompCheck[(int) ci][(int) cj] = 1;
			}
		}
	}

	return bRevCompCheck;
}
