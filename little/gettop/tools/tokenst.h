
#ifndef _TOKENST_H_
#define _TOKENST_H_

#include<stdio.h>
#include<stdlib.h>

/*
check whether the strng contain one of the character in delimstr,
*/
int found_token(char *strng, char *delimstr, int zindex);


int count_tokens(char *tokens, char *delimstr);

char *gettoken(char *tokens, char *delimstr, int tokenindex);



#endif  /* _TOKENST_H_ */
