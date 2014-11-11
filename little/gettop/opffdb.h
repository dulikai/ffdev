

#ifndef _OPFFDB_H_
#define _OPFFDB_H_

#include "ffdb.h"

#define ATOMS_TABLE 0
#define TRANSLATIONS_TABLE 1
#define BONDS_TABLE 2
#define ANGLES_TABLE 3
#define DIHEDRALS_TABLE 4
#define IMPROPERS_TABLE 5
#define INDEX_TABLE 6

 void *ff_db_select(aFFP *af, int select_table, int *node);

 void *ff_db_select_atoms(aFFP *af, int inode);

 void *ff_db_select_translations(aFFP *af, int inode);

 void *ff_db_select_bonds(aFFP *af, int *node);

 void *ff_db_select_angles(aFFP *af, int *node);

 void *ff_db_select_dihedrals(aFFP *af, int *node);

 void *ff_db_select_index(aFFP *af, int inode);

 static void get_ff_atom_name(aFFP *af, int inode, char *name);


 void get_ffatom_names(aFFP *af, int *node, int n, char *line, char *link);

 void *ff_db_select_impropers(aFFP *af, int *node);

 void translate_db_inode(int *tnode, aFFP *af, int *node, int n);


 void get_atomtype_name(aFFP *af, int inode, char *name);


#endif  /* _OPFFDB_H_ */
