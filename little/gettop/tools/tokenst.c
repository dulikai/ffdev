/*****
This file is part of the Babel Program
Copyright (C) 1992-94 W. Patrick Walters and Matthew T. Stahl
The Babel Program is a product of the Dolata Research Group
Dept. of Chemistry
University of Arizona
Tucson, AZ 85721

For more information please contact :

babel@mercury.aichem.arizona.edu
--------------------------------------------------------------------------------

FILE : tokenst.c
AUTHOR(S) : Pat Walters
DATE : 10-92
PURPOSE : routines to determine the number of tokens in a string
from Turbo Algorithms
by Keith Weiskamp
   Namir Shammas
   Ron Pronk

******/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "tokenst.h"
/*
delimstr means a list of character used as delimiter....
*/
int 
found_token(char *strng, char *delimstr, int zindex)
{
	if (zindex > 0)
		return (((strchr(delimstr, strng[zindex-1])) != NULL) &&
			((strchr(delimstr, strng[zindex])) == NULL));
	else
		return 0;
}


int 
count_tokens(char *tokens, char *delimstr)
{
	int i, count, strlength;
	int all_tokens, no_tokens, is_token;
	char *tempstr;

	if (tokens[0] == '\0')
		return 0;
	if (delimstr[0] == '\0')
		strcpy(delimstr," ");
	tempstr = (char *)malloc(strlen(tokens) + 2);
	tempstr[0] = delimstr[0];
	tempstr[1] = '\0';
	strcat(tempstr,tokens);
	strlength = strlen(tempstr);
	all_tokens = 1;
	no_tokens = 1;
	i = 0;
	while ((i < strlength) && (all_tokens || no_tokens))
	{
		is_token = ((strchr(delimstr, tempstr[i]) != NULL));
		no_tokens = ((!is_token) && no_tokens);
		all_tokens = (is_token && all_tokens);
		i++;
	}
	if (all_tokens || no_tokens)
	{
		if (tempstr) free(tempstr);
		if (no_tokens)
			return 1;
		if (all_tokens)
			return 0;
	}
	count = 0;
	for (i = 0; i <strlength; i ++)
		if (found_token(tempstr, delimstr, i))
			count ++;
	if (tempstr) free(tempstr);
	return count;
}

char *
gettoken(char *tokens, char *delimstr, int tokenindex)
{
	int i,n,strlength,ptr;
	char ch, *tempstr;

	if (tokenindex < 1) tokenindex = 1;
	if (delimstr[0] == '\0') strcpy(delimstr, " ");
	ch = delimstr[0];
	tempstr = (char *)malloc(strlen(tokens) + 3);
	if (!tempstr) 
	{
	  printf("UNABLE TO ALLOCATE MEMORY");
	  exit(1);
	}	  
	tempstr[0] = ch;
	tempstr[1] = '\0';
	strcat(tempstr, tokens);
	strlength = strlen(tempstr);
	tempstr[strlength] = ch;
	tempstr[++strlength] = '\0';
	n = tokenindex;
	i = 1;
	while ((i < strlength) && (n > 0))
	{
		if (found_token(tempstr,delimstr,i)) n --;
		if (n > 0) i++;
	}
	if ( i >= strlength)
	{
		if (tempstr) free(tempstr);
		return NULL;
	}
	else
	{
		ptr = i;
		while ((i < strlength) && (strchr(delimstr, tempstr[i]) == NULL))
			i++;
		tempstr = tempstr + ptr;
		tempstr[i - ptr] = '\0';
		return tempstr;
	}
}




