
#ifndef _OPSTR_H_
#define _OPSTR_H_

#include<stdio.h>
#include<stdlib.h>
#include<string.h>


/*
 The rtrim() function removes trailing spaces from a string.
 */

char * rtrim(char * str);



/* 

**The ltrim() function removes leading spaces from a string. 
**use rtrim() to realize this function

*/ 

char *ltrim(char * str);


/*

  The rjust() function right-justifies a string. 

*/

char *rjust (char * str);


#endif  /* _OPSTR_H_ */