
#ifndef _PRINT_H_
#define _PRINT_H_

#include "fileinit.h"

#include "structure.h"
#include "ffpara.h"
#include "communicate.h"

 void print_itp(char *file, NOBITPFILE *itp);

  void print_itp_nob(char *file, NOBITPFILE *itp);


 void print_top(FILEControl *fc, NOBITPFILE *itp);

#endif /* _PRINT_H_ */
