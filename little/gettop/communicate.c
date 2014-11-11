
#include "communicate.h"

 void comm_say_hello(FILEControl *fc, analyzeGraphic *com, aFFP *af, NOBITPFILE *itp)
 {
    ItpControl itpc;  
      
    con_com2ff(com, af);
   // print_com("com.txt", &com);

    //itp...
    init_itp_behavior(fc->templatefile, &itpc);
    //print_itp_behavior("template.txt", &itpc);

    gen_itp(com, af, &itpc, itp);


    return;         
 }
 
