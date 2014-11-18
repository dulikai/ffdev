      SUBROUTINE i2c
C Converts Internal coordinates to cartesian coordinates.
  
      IMPLICIT NONE
  
      INTEGER MAT, MRS
      PARAMETER (MAT = 500, MRS = 20)
      INTEGER lo
      REAL rot(4,4), x(3,MAT), AM(MAT, 4, 4)
      REAL a, b, r
      INTEGER i, j, k, l, m, i1
  
C                   *** Molecular Topology
C  bl=bond length array <<<<
C  ba=bond angle array <<<<<
C  da= dihedral angle array <<<<<
C  natms= total number of atoms <<<<<
C name=atom names array
C con(4,..)= connectivity table. only con(1,i) is usefule here.<<<<<
C chrg= partial atomic charge array
C type= atomic type array
C ngrup=total number of residues
C grup=group(residue) numbers
C resd=residue names
C X(3,..) = cartesian coorinate array <<<<<<

C ' <<<<<<' marked arrays are used here. rest are for compatibility with
C other programs in MoMoS package and related stuff.
C
C                                Authors 

C  DR.  Mrigank                         Dr. V.Kothekar
C  Bioinformatics Centre                Dept. of Biophysics
C  Institute of Microbial Technology    All India Institute of Medical Sciences
C  Sector 39-A                          Ansari Nagar
C  Chandigarh 160 014                   New Delhi 110 029
C  mrigank@imtech.ernet.in[much preferred]  koth@medinst.ernet.in,
C  mrigank@csimtech.ren.nic.in              vkothek@aiims.ernet.in
C  tele: +91-172-690557                     +91-11-6864851
C 
C This routine is part of MoMoS [MOlecular MOdelling and Simulation] package
C                                --        --            -
      REAL bl(MAT), ba(MAT), da(MAT), chrg(MAT)
      INTEGER con(4,MAT), type(MAT), grup(MAT), natms, ngrp
      CHARACTER title*78, name(MAT)*4, resd(MRS)*4
 
      COMMON
     +      /mtr/ bl, ba, da, chrg
     +      /mti/ con, type, grup, natms, ngrp
     +      /mtc/ title, name, resd
  
      COMMON
     +     /crt/ x
  
      IF(natms.EQ.0) THEN
C
       RETURN
      END IF
  
      DO 10 i = 1, natms
  
       R=bl(i)
       A=0.017453*ba(i)
       B=0.017453*da(i)
       ROT(1,1)=-COS(A)
       ROT(1,2)=-SIN(A)
       ROT(1,3)=0.0
       ROT(1,4) = -R*COS(A)
       ROT(2,1)= SIN(A)*COS(B)
       ROT(2,2)=-COS(A)*COS(B)
       ROT(2,3)=-SIN(B)
       ROT(2,4)=R*SIN(A)*COS(B)
       ROT(3,1)=SIN(A)*SIN(B)
       ROT(3,2)=-COS(A)*SIN(B)
       ROT(3,3)=COS(B)
       ROT(3,4)=R*SIN(A)*SIN(B)
       ROT(4,1)=0.0
       ROT(4,2)=0.0
       ROT(4,3)=0.0
       ROT(4,4)=1.0
       I1=con(1,i)
       DO 70 L=1,4
        DO 70 M=1,4
         AM(i,L,M)=0.0
          IF(i.GT.1)GO TO 30
          AM(i,L,M)=ROT(L,M)
          GO TO 70
 30            DO 60 K=1,4
           AM(i,L,M)= AM(i,L,M)+AM(i1,L,K)*ROT(K,M)
 60             CONTINUE
 70              CONTINUE
      DO 90 J=1,3
       X(j,i)=AM(i,J,4)
  
 90     CONTINUE
 10      CONTINUE
  
 3        CONTINUE
  
c     DO i = 1, natms
c      WRITE(lo, '(A4, 3F10.5, 4I4, 2X, F8.3, 2I2)')
c     +       name(i), (x(j,i), j = 1, 3), (con(j,i), j = 1,4), chrg(i),
c    +       name(i), x(1,i), x(2,i), x(3,i), (con(j,i), j = 1,4),
c    +      chrg(i), type(i), grup(i)
c     END DO
  
 5000       RETURN
 2004         FORMAT(2A2,3F10.5,4I4,F10.3,2I2)
      END
