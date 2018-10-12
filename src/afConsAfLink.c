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
#include <limits.h>
#include <sys/time.h>

#include "af.h"

typedef struct {
	LInterval *elem;
	int     pointer;
} Stack;

void push(Stack *stack, LInterval *elem) {
	stack->pointer++;

	if(stack->pointer % BUFFER1 == 0) {
		if((stack->elem = realloc(stack->elem, (stack->pointer + BUFFER1) * sizeof(LInterval))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"stack->elem\".\n");
			exit(1);
		}
	}

	stack->elem[stack->pointer].i       = elem->i;
	stack->elem[stack->pointer].lcp     = elem->lcp;
}

void pop(Stack *stack, LInterval *elem) {
	if(stack->pointer==0) {
		fprintf(stderr,"Stack empty (pop function)\n");
		exit(1);
	} else {
		elem->i       = stack->elem[stack->pointer].i;
		elem->lcp     = stack->elem[stack->pointer].lcp;
		stack->pointer--;
	}

	if(stack->pointer % BUFFER1 == 0) {
		if((stack->elem = realloc(stack->elem, (stack->pointer + BUFFER1) * sizeof(LInterval))) == NULL) {
			fprintf(stderr,"Memory allocation failed for \"stack->elem\".\n");
			exit(1);
		}
	}
}

void top(Stack *stack, LInterval *elem) {
	if(stack->pointer==0) {
		fprintf(stderr,"Stack empty (top function)\n");
		exit(1);
	} else {
		elem->i       = stack->elem[stack->pointer].i;
		elem->lcp     = stack->elem[stack->pointer].lcp;
	}
}

void initStack(Stack *stack) {
	stack->pointer = 0;
	if((stack->elem = (LInterval *) calloc(BUFFER1, sizeof(LInterval))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"stack->elem\".\n");
		exit(1);
	}
}

void constructAflk(EArray *earray, EArray *earrayR, bool isEsa, int (*findAffixLink)(void *, void *, int, int, void *, void *), unsigned char* convSequences, int length, Alphabet *alphabet, bool argShowTimes) {
	int c, j, lb, home;
	LInterval elem;
	Stack stack;
	unsigned char *pattern;
	struct timeval start, end;

	printf("\nComputing %s with %s... ", isEsa ? "aflk" : "aflkr", isEsa ? "skpr" : "skp");
	fflush(stdout);

	if (argShowTimes) gettimeofday(&start, NULL);

	if((earray->affixLink = (int *) calloc(length, sizeof(int))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"affixLink\". - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	initStack(&stack);

	for(c = 0; c < length; c++) {
		earray->affixLink[c] = INT_MAX;
	}

	/* root interval representing the empty string */
	elem.i       = 0;
	elem.lcp     = 0;
	push(&stack, &elem);

	for(c = 1; c < length; c++) {
		lb= c - 1;
		while(xlcpvalue(c) < elem.lcp) {
			j = c - 1;
			pop(&stack, &elem);

			home = j == length - 1 ? length - 1 : (xlcpvalue(elem.i) < xlcpvalue(j + 1) ? j : elem.i);
			if (isEsa)
				pattern = convSequences + earray->xarray[elem.i];
			else
				pattern = convSequences + (length - 1) - earray->xarray[elem.i] - elem.lcp;
			earray->affixLink[home] = (*findAffixLink)(convSequences, pattern, elem.lcp, length, earrayR, alphabet);

			lb = elem.i;
			top(&stack, &elem);
		}
		if(xlcpvalue(c) > elem.lcp) {
			elem.i = lb;
			elem.lcp = xlcpvalue(c);
			push(&stack, &elem);
		}
	}

	free(stack.elem);

	if (argShowTimes) gettimeofday(&end, NULL);
	printf("done\n");
	if (argShowTimes) printf("Time: %f ms\n", (double)( (end.tv_sec - start.tv_sec)*1000 + (end.tv_usec - start.tv_usec)/1000.0) );
}

void constructAflkBS(EArray *earray, EArray *earrayR, bool isEsa, int* (*findPatternBS)(void *, int, int, int, int, void *), unsigned char* convSequences, int length, AffixArray *affixArray, bool argShowTimes) {
	int c, j, lb, home, *result;
	LInterval elem;
	Stack stack;
	unsigned char *pattern;
	struct timeval start, end;

	printf("\nComputing %s with binary search... ", isEsa ? "aflk" : "aflkr");
	fflush(stdout);

	if (argShowTimes) gettimeofday(&start, NULL);

	if((earray->affixLink = (int *) calloc(length, sizeof(int))) == NULL) {
		fprintf(stderr,"Memory allocation failed for \"affixLink\". - %s %d.\n", __FILE__, __LINE__);
		exit(1);
	}

	initStack(&stack);

	for(c = 0; c < length; c++) {
		earray->affixLink[c] = INT_MAX;
	}

	/* root interval representing the empty string */
	elem.i       = 0;
	elem.lcp     = 0;
	push(&stack, &elem);

	for(c = 1; c < length; c++) {
		lb= c - 1;
		while(xlcpvalue(c) < elem.lcp) {
			j = c - 1;
			pop(&stack, &elem);

			home = j == length - 1 ? length - 1 : (xlcpvalue(elem.i) < xlcpvalue(j + 1) ? j : elem.i);
			if (isEsa)
				pattern = convSequences + earray->xarray[elem.i] + elem.lcp - 1;
			else
				pattern = convSequences + (length - 1) - earray->xarray[elem.i] - elem.lcp;
			result = (*findPatternBS)(pattern, elem.lcp, 0, affixArray->length - 1, 0, affixArray);
			earray->affixLink[home] = result[0];
			free(result);

			lb = elem.i;
			top(&stack, &elem);
		}
		if(xlcpvalue(c) > elem.lcp) {
			elem.i = lb;
			elem.lcp = xlcpvalue(c);
			push(&stack, &elem);
		}
	}

	free(stack.elem);

	if (argShowTimes) gettimeofday(&end, NULL);
	printf("done\n");
	if (argShowTimes) printf("Time: %f ms\n", (double)( (end.tv_sec - start.tv_sec)*1000 + (end.tv_usec - start.tv_usec)/1000.0) );
}
