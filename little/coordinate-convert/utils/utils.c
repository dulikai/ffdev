
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "utils.h"

 int get_element_order(char *element_name)
 {
	 const char *element[] = {"H","He","Li","Be","B","C","N","O","F","Ne",
		 "Na","Mg","Al","Si","P","S","Cl","Ar"};

	 const int len = 18;

	 int i, orderid = 0;

	 for(i = 0; i < len; i++)
	 {
	   if(!stricmp(element[i], element_name))
	   {
	     orderid = i+1;
		 break;
	   }
	 }
     return orderid;
 }


 int to_exist(char *str, int signal, char cmd)
 {
   int flag = 0;
   printf("warning: %s", str);

   ;

   if(cmd != getchar())
   {
     exit(signal);
   }
   else
   {
     flag = 1;
   }
   return flag;

 }