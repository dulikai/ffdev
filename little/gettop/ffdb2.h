
#ifndef _FFDB2_H_
#define _FFDB2_H_

#include "readdb2.h"

#define AA_ATOMS_TABLE 1
#define AA_BONDS_TABLE 2
#define AA_ANGLES_TABLE 3


#define AA_ANGLE_TH 0
#define AA_ANGLE_KTH 1
#define AA_BOND_PARA_R0 0
#define AA_BOND_PARA_KR 1

#define AA_dihedral_PARA_C0 0
#define AA_dihedral_PARA_C1 1



 void convert_aaffp2afp(AAFFP *aa, aFFP *af, char *entity);



#endif  /* _FFDB2_H_ */
