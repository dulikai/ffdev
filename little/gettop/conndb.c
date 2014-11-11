

#include "conndb.h"


 FILE *connect_ff_db(char *file)
 {
   FILE *fp;
   
   fp = fopen(file, "r");

   if(fp == NULL)
   {
     printf("connect to ff db failed!\n, check file: %s\n", file);
	 exit(1);
   }

   return fp;
 }

 void close_ff_db(FILE *fp)
 {
   fclose(fp);
 }

