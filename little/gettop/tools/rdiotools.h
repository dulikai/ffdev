
#ifndef _RDIOTOOLS_H_
#define _RDIOTOOLS_H_

#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<string.h>
#include<conio.h>
#include<stdarg.h>


/*
  if end return 0;
*/
 extern int check_file_end(FILE *fp);

/*
 ignore all connected space character.
 no matter how many it is.
 the file pointer is move to the first character
 which is not a space character..
*/
 extern void JmpSpace(FILE *fp);
 /*
 jump the all connected lines, 
 if the first non-space character in this line is cht...
 */
 extern void JmpOrNot(FILE *fp,char cht);
 extern void rd_name_lst(FILE *fp,char *name);
 extern int  find_elements_id(char *name,char *Elements[],int element_num);
 extern void rd_dollar_line(FILE *fp,char *line);
 extern int  chk_endnn(FILE *fpin,int n);

 void printfile(char *filename,char *access,char *tips,char *first,...);

/*
 one subroutine get from bable1.1 -- utils.c
 this will jump one line.
*/
 void toss(FILE *file1, int num);
 
 /*
 this will check a string of blank. 
 */
 int is_blank_line(char *str);

//go to a line that contain the character cht, and the decide to go forward/backward 'step' len
 long goto_target_line(FILE *fp, char cht, long step);

#endif  /* _RDIOTOOLS_H_ */
