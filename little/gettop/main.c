#include <stdio.h>
#include <stdlib.h>

#include "fileinit.h"
#include "structure.h"
#include "ffpara.h"
#include "communicate.h"

#include "print.h"

//#include "graphic.h"
//#include "ffdb.h"
//#include "handshake.h"


 void process_start(char *file)
 {
    FILEControl fc;
    //structGraphic g;
    //FFP *f;
    aFFP af;
    analyzeGraphic com;
    //ItpControl itpc;
    NOBITPFILE itp;

    init_filecontrol(file, &fc);

    gen_struct(&fc, &com);

    gen_ffpara(&fc, &af);

    comm_say_hello(&fc, &com, &af,  &itp);

    print_top(&fc, &itp);

    //read_gjf_cartesian(fc.structfile, &g);
   //print_graphic("gjf.txt", &g);

    //f = init_ffp();
    //read_ff_db(fc.dbfile, f);
   //print_ff_db("outff.txt", f);

   //convert_ffp2affp(f, &af);
   //print_aff("aff.txt", &af);


    //gen_com(&g, &com);
   //print_com("com.txt", &com);


    //con_com2ff(&com, f);
   // print_com("com.txt", &com);

   //itp...
   //init_itp_behavior(fc.templatefile, &itpc);
   //print_itp_behavior("template.txt", &itpc);



   //gen_itp(&com, &af, &itpc, &itp);

   //print_itp(fc.outitpfile, &itp);

   //print_itp_nob(fc.itp_nobfile, &itp);

   return;


 }

#include "readdb2.h"
#include "ffdb2.h"
int main(int argc, char *argv[])
{

   char *file = "control.txt";

   process_start(file);
/*
   AAFFP aa;
   aFFP af;
   read_aaff_db("amino-acid.ff", &aa);

   convert_aaffp2afp(&aa, &af, "Ala");

   print_aaffp("aaffp.txt", &aa);
   print_aff("aff.txt", &af);
*/
   system("PAUSE");	
   return 0;
}
