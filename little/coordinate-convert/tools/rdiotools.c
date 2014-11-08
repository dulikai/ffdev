

#include"rdiotools.h"


/*
  if end return 0;
*/
 int check_file_end(FILE *fp)
 {
  char ch;
  int i = 1;

  ch=fgetc(fp);

  while(isspace(ch))
  {
   ch=fgetc(fp);

   if(feof(fp))
   {
     i = 0;
   }
  }
  fseek(fp,-1L,SEEK_CUR);

  return i;
 }





/*
******************************************
*/

 void JmpSpace(FILE *fp)
 {
  char ch;

  ch=fgetc(fp);

  while(isspace(ch))
  {
   ch=fgetc(fp);
  }
  fseek(fp,-1L,SEEK_CUR);
 }

void JmpOrNot(FILE *fp,char cht)
 {
  char ch,strbox[512];
  
  ch=fgetc(fp);
  while(ch==cht||isspace(ch))
  {
	if(isspace(ch))
	{
     ch=fgetc(fp);
	}
	else
	{
     fgets(strbox,511,fp);
     ch=fgetc(fp);
	}
  }
  fseek(fp,-1L,SEEK_CUR);
 }

 long goto_target_line(FILE *fp, char cht, long step)
 {
   char ch;   

   ch=fgetc(fp);
   while(ch != cht)
   {
      ch=fgetc(fp);
   }

  fseek(fp, -1L, SEEK_CUR);    
   
  fseek(fp, step, SEEK_CUR);    
  
  return step;   
 }

/*
read string and link them.. useless now...
$STO-3G~ STO-6G~ 
*/
 void rd_name_lst(FILE *fp,char *name)
 {
  char ch,tmpstr[20];
  int len;

  strcpy(name,"");

  while((ch=fgetc(fp))!='$')
  {
   ch=fgetc(fp);
  }

  fscanf(fp,"%s",tmpstr);
  len=strlen(tmpstr);

  while(tmpstr[len-1]=='~')
  {
   tmpstr[len-1]=' ';
   strcat(name,tmpstr);

   fscanf(fp,"%s",tmpstr);
   len=strlen(tmpstr);
  }
  strcat(name,tmpstr);
 }


int find_elements_id(char *name,char *Elements[],int element_num)
  {
   int i;
   for(i=0;i<element_num;i++)
   {
    if(!strcmp(Elements[i],name))break;
   }
   if(i>element_num+1){printf("error in Elements Table.\n");getch();return -1;}
   else 
	   return i+1;
  }



 void rd_dollar_line(FILE *fp,char *line)
 {
  char ch;

  strcpy(line,"");

  while((ch=fgetc(fp))!='$')
  {
   ch=fgetc(fp);
  }

  fgets(line,199,fp);
 }



 int chk_endnn(FILE *fpin,int n)
 {
  char ch;
  int flag=1,i=1;

  ch=fgetc(fpin);
  while(isspace(ch))
  {
   ch=fgetc(fpin);
   if(ch=='\n')i++;
   if(i>=n){flag=0;break;}
  }
  if(flag)
    fseek(fpin,-1L,SEEK_CUR);

  return flag;
 }

/*
************************************************************************
*/

void printfile(char *filename,char *access,char *tips,char *first,...)
{
  FILE *fp;
  va_list argp;
  char mode[3];

  strcpy(mode,access);

  if(!strcmp(mode,""))strcpy(mode,"w");

  fp=fopen(filename,mode);

  fprintf(fp,"%s\n",tips);

  va_start(argp,first);

  vfprintf(fp,first,argp);

  va_end(argp);

  fprintf(fp,"\n");

  fclose(fp);

}


/*
jump num lines, and ignore them...
*/
void toss(FILE *file1, int num)
{
  int i;
  char buffer[512];

  for (i = 0; i < num; i++)
    fgets(buffer,sizeof(buffer),file1);
}


/*

*/


int is_blank_line(char *str)
{
    int i = 0;
    int flag = 1;
    while(str[i] != '\0')
    {
       if(!isspace(str[i]))
       {
           flag = 0;
           break;                
        }
       i++;                    
    }   
    return flag; 
 }


