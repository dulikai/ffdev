

#ifndef _FILEINIT_H_
#define _FILEINIT_H_

typedef struct{
	    /* structure contral */
        char structfile[100];
        char structtype[20];
        
		/* db control parameter */
        char dbfile[100];
        char dbtype[10];
		char dbpara[20];

		/* template of output file */
        char templatefile[100];

		/* output files */
        char outitpfile[100];
        char itp_nobfile[100]; 
        
        }FILEControl;
        

 void init_filecontrol(char *file, FILEControl *fc);


#endif
