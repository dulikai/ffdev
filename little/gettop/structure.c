
#include "structure.h"




 void gen_struct(FILEControl *fc, analyzeGraphic *com)
 {
     structGraphic g;
     
     //now only gjf cartesian is aviable..
    read_gjf_cartesian(fc->structfile, &g);
    //print_graphic("gjf.txt", &g);    
   
   
    gen_com(&g, com);
    //print_com("com.txt", com);      
   
   return;
      
      }
