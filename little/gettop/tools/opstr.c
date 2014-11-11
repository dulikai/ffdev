

#include "opstr.h"

/*
 The rtrim() function removes trailing spaces from a string.
 */
char * rtrim(char * str)
{    
   int n = strlen(str)-1; /* Start at the character BEFORE the null character (\0).  */ 
    while (n>0)  /* Make sure we don't go out of hounds. . . */
    {
        if ( * (str + n) !=' ') /* If we find a nonspace character: */
        {
            * (str+n+1) = '\0' ; /* Put the null character at one character past our current  position. */
			
            break ;  /* Break out of the loop. */
        }
        else  /* Otherwise , keep moving backward in the string. */
		{

            n--;
        }
    }
    return str;    /*Return a pointer to the string*/
}



/* 
**The ltrim() function removes leading spaces from a string. 
**use rtrim() to realize this function
*/ 
char *ltrim(char * str)
{
    strrev(str); /* Call strrevO to reverse the string. */

    rtrim(str); /* Call rtrimO to remvoe the "trailing" spaces. */

    strrev(str); /* Restore the string's original order. */

    return str ; /* Return a pointer to the string. */
}



/*

  The rjust() function right-justifies a string. 

*/

char *rjust (char * str)
{
    int n = strlen(str); /* Save the original length of the string. */

    char* dup_str;

    dup_str = strdup(str); /* Make an exact duplicate of the string. */

    rtrim(dup_str); /* Trim off the trailing spaces. */

    /* Call sprintf () to do a virtual "printf" back into the original

        string. By passing sprintf () the length of the original string,

        we force the output to be the same size as the original, and by

        default the sprintf() right-justifies the output. The sprintf()

        function fills the beginning of the string with spaces to make

        it the same size as the original string.*/ 

    sprintf(str, "%*. * s", n, n, dup_str); 

    free(dup_str) ; /* Free the memory taken by he duplicated string. */

    return str;  /* Return a pointer to the string. */

}


 char *trim(char *str)
 {
      rtrim(str);
      ltrim(str);
      
      return str;
  }

/*void string2array(char *str, char *delimstr, char (*str)[100])
{
      
     }
*/
