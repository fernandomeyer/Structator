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

#ifndef ESA_H_
#define ESA_H_ 1

#include "stdbool.h"
#include "chaining/chain2dim.h"

#define BUFFER1 1000
#define BUFFER2 100
#define BUFFER3 10
#define $ 127
#define LINESYMBOL '!'

typedef struct {
	bool bSearchForwardString;
	bool bSearchReverseString;
	bool bBED;
	bool bReportAllMatches;
	bool bIncludeSeqDesc;
	char *cTextFile; //file for storing matches
	unsigned int uiMinDiffMatches;
} SearchParam;

typedef struct {
	int  numSeqs;
	unsigned char *sequences;
	unsigned char *convSequences; //alphabetically converted sequences
	bool sequencesMmapped;
	bool convSequencesMmapped;
	int  *seqEndPos; //ending position of each sequence + '$' in the concatenated array
	char **seqDescription;
	int  *seqDescLength;
} MultiSeq;

typedef struct {
	int *index;
	int *value;
	int numExceptions;
} LcpException;

typedef struct {
	char **eqClass;
	int  *classSize;
	unsigned char *classRepresentative;
	int  numClasses;
	bool *isWildCard;
	bool **iupacTable;
} Alphabet;

/* Modified by Benjamin Albrecht, Jan. 2012 */
typedef struct {
	int  *xarray; 					//suf or rpref
	unsigned char *xlcp; 		//lcp or rlcp
	unsigned char *xlcpTree; 	//lcpTree or rlcpTree
	LcpException  *xlcpException; 	//lcpException or rlcpException
	LcpException *xlcpTreeException;//lcpTreeException or rlcpTreeException
	int  *xskp; 					//skp or rskp
	int  *affixLink;
} EArray;

typedef struct {
	MultiSeq *multiSeq;
	Alphabet *alphabet;
	int      length;
	EArray   *esa;
	EArray   *erpa;
} AffixArray;


typedef struct {
	int i;
	//Int32 j; //save space by not storing j
	int lcp;
} LInterval;

typedef struct {
	bool *isPairing;
	bool *isLeftRightExt;
	unsigned int length;
} StructureInfo;

typedef struct {
	int context;
	int i;
	int j;
	int d;
	unsigned char type; //0=suf, 1=rpref, > 1 "single occurrence" node
	StructureInfo structureInfo;
} VirtualNode;

typedef struct {
	char *desc;
	unsigned char *seq;
	unsigned char *structure;
	unsigned char *seqR; // pattern seq for reverse search
	unsigned char *structureR; // pattern structure for reverse search
	int  weight;
	int  startpos;
	int  instance;
	unsigned int maxstemlength;
	unsigned int maxrightloopextent;
	unsigned int maxleftloopextent;
	unsigned int maxmispair;
} Pattern;

typedef struct {
	Pattern *pattern;
	int     numPatterns;
	unsigned char usesWeightParam;
	unsigned char usesStartPosParam;
	unsigned char usesInstanceParam;
} MultiPattern;

typedef struct {
	int seqId;
	unsigned int pos;
	unsigned int endPos;
	int patternId;
	char *structureString;
} Match;

typedef struct {
	unsigned int startpos[2];
	unsigned int endpos[2];
	int weight;
	char *structureString;
} ChainFragment;

typedef struct {
	int  seqId;
	int  chainLength;
	int  chainScore;
	ChainFragment *fragmentInfo;
	bool forwardSeq;
	bool bStoredInReverseOrder;
} Chain;

typedef struct {
	bool isactive;
	char *reportfile;
	bool islocal;
	double weightfactor;
	GtChain2Dimpostype maxgap;
	char *localparm, *globalparm, *origopts;
	bool show;
	unsigned int minchainlength;
	unsigned int minchainscore;
	unsigned int topkscoring;
} Chainparam;

typedef struct {
	unsigned char extension[100];
	unsigned char toRight;
} StringExtension;

/* Modified by Benjamin Albrecht, Jan. 2012 */
#define lcpTreevalueX(data, index) (data->xlcpTree[index] == 255 ? \
		xlcpExceptionValue(data->xlcpTreeException, index) : data->xlcpTree[index])

#define lcpvalueX(data, index) (data->xlcp[index] == 255 ? \
		xlcpExceptionValue(data->xlcpException, index) : data->xlcp[index])

#define lcpvalue(index) (affixArray->esa->xlcp[index] == 255 ? \
		xlcpExceptionValue(affixArray->esa->xlcpException, index) : affixArray->esa->xlcp[index])

#define rlcpvalue(index) (affixArray->erpa->xlcp[index] == 255 ? \
		xlcpExceptionValue(affixArray->erpa->xlcpException, index) : affixArray->erpa->xlcp[index])

#define xlcpvalue(index) (earray->xlcp[index] == 255 ? \
		xlcpExceptionValue(earray->xlcpException, index) : earray->xlcp[index])

#define c(A) affixArray->alphabet->classRepresentative[(A) - 1]

//afCons.c
void reverseSequences(unsigned char*, int);
unsigned char* reverseStringNewPointer(unsigned char*, unsigned int);
void affXSuf(EArray*, unsigned char*, int);
int  xlcpExceptionValue(LcpException*, int);
void affXLcp(EArray*, Alphabet*, unsigned char*, int);
void affXSkp(EArray*, int, int);

/* Modified by Benjamin Albrecht, Jan. 2012 */
void computeAffixArray(AffixArray*, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, char*, bool);
void computeEArray(EArray*, Alphabet*, unsigned char*, int, int, bool, bool, bool, bool, bool, bool);
void affXLcpTree(EArray*, bool, int);

//afConsAfLink.c
void constructAflk(EArray*, EArray*, bool, int (*)(void *, void *, int, int, void *, void *), unsigned char*, int, Alphabet*, bool);
void constructAflkBS(EArray*, EArray*, bool, int* (*)(void *, int, int, int, int, void *), unsigned char*, int, AffixArray*, bool);
void affixLinkSuf(LInterval*, AffixArray*);
void affixLinkRPref(LInterval*, AffixArray*);
void affixLinkSufBS(LInterval*, AffixArray*);
void affixLinkRprefBS(LInterval*, AffixArray*);

//afAlphabet.c
void  setPredefinedAlphabet(int, Alphabet*);
bool convertToAlphabet(unsigned char*, unsigned char*, int, bool, Alphabet*);
bool convertToRepresentative(char*, int, Alphabet*);
bool** setIupacTable (Alphabet*);

//afStreamHandling.c

/* Modified by Benjamin Albrecht, Jan. 2012 */
void init(AffixArray*, Alphabet*, char*, int, bool*, bool*, bool*, bool*, bool*, bool*, bool*, bool*, bool*, bool*, bool);
bool loadAffFiles(AffixArray*, char*, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool);

void freeAll(AffixArray*, MultiPattern*, bool**);
bool loadFastaFile(AffixArray*, char*, bool);
Alphabet*     loadAlphabetFile(char*, Alphabet*);
bool saveAffFiles(AffixArray*, char*, bool);
bool saveEArray(EArray*, int, bool, char*);
bool saveAflk(EArray*, int, bool, char*);
bool saveTextFile(AffixArray*, char*, bool);
void allocConvSequences (AffixArray*);
bool printEArray(AffixArray*, bool);
void printSearchResult(AffixArray*, char*, int*, int);
void saveSearchResults(Pattern*, int*, int, char*, bool, bool, AffixArray*);
int  loadPatternFileBidSearchDebug(char*, StringExtension*, AffixArray*);
MultiPattern* loadPatternFile(char*, AffixArray*);
bool** loadComplementarityFile(char*, AffixArray*);
bool** loadDefaultComplementarityRules(bool, int, AffixArray*);
bool** loadReverseComplementarityFile(bool**, AffixArray*);

/* Modified by Benjamin Albrecht, Jan. 2012 */
//afComputeAflk.c
int* cmpAflkViaBS(int, int, AffixArray *, unsigned char*, int, bool);
int* cmpAflkViaImprovedBS(int, int, int, int, AffixArray*, const unsigned char*, int, bool);

//afSearch.c
void printNode(VirtualNode*, AffixArray*);
int  getSeqNumber(int *range, int length, int k);
int* findPattern(unsigned char*, int, AffixArray*);
int* findPatternSlow(unsigned char*, int, AffixArray*);
int* findPatternBS(unsigned char*, int, int, int, int, AffixArray*);
int* findPatternBS2(unsigned char*, int, int, int, int, AffixArray*);
int getLcpIntervalDepth(int, int, int, AffixArray*);
int getRLcpIntervalDepth(int, int, int, AffixArray*);
int* getRelativePos(int*, int, AffixArray*);
int findAffixLink(unsigned char*, unsigned char*, int, int, EArray*, Alphabet*);
int findAffixLink2(unsigned char*, unsigned char*, int, int, EArray*, Alphabet*);
bool extendLeft(unsigned char*, int, VirtualNode*, AffixArray*);
bool extendRight(unsigned char*, int, VirtualNode*, AffixArray*);

//afSearchStr.c
void setSearchPatterns(MultiPattern*, bool, bool, Alphabet*);
void match(MultiPattern*, SearchParam*, Chainparam*, bool**, bool, bool, bool, AffixArray*);
void matchHairpinLoop(int, int, bool, VirtualNode*, AffixArray*);
void matchStem(int, int, int, int, VirtualNode*, AffixArray*);
void matchSlow(MultiPattern*, SearchParam*, Chainparam*, bool**, bool, bool, bool, AffixArray*);

int BM(char*, int, char*, int, int, int*);

//afResultsProcessing.c
int compareRB(const void *pa, const void *pb, const void *config);
int compareRBCR(const void *pa, const void *pb, const void *config);
int removeNonsignificantMatches(void*, unsigned int*, int, int);
void saveSearchResultsRB(void*, void*, char*, bool, bool, Pattern*, int, AffixArray*, bool);
void reportChainingResults(void*, Chainparam*, AffixArray*);

//afChaining.c
void chainAll(void*, void*, unsigned int*, Chainparam*, MultiPattern *multiPattern, AffixArray *affixArray, bool, bool);

#endif /*ESA_H_*/
