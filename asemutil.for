cc      MODULE ASEMUTIL
cc      CONTAINS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C         --G E N E R A L  A S S E M B L Y  S U B R O U T I N E S--    C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE NODSTFFASEM                                           C
C                                                                      C
C     This is the actual solver at each increment.                     C
C                                                                      C
C     UG(,)     :Global displacements history at each node             C
C     VG(,)     :Global velocity histories at each node                C
C     AG(,)     :Global acceleration history at each node              C
C     NUMNP     :Number of nodal points                                C
C     NUMEL     :Number of elements                                    C
C     NUMAT     :Number of material profiles                           C
C     NNE()     :Number of nodes per element                           C
C     NMNE      :Maximum number of nodes per element                   C
C     MXDOFDIM  :Maximum degree of freedom dimension                   C
C     MXDOFEL   :Maximum degree of freedom per element                 C
C     IELT()    :Array of element type identifiers                     C
C     IELCON(,) :Elementary connectivity array                         C
C     ILIST()   :Incidence list.Elements connected to the current node C
C     LPLIST()  :Local position of the current node in each one of the C
C                elements of ILIST()                                   C
C     NIEL()    :Number of elements at the current node                C
C     NDOFN()   :Number of degrees of freedom per node                 C
C     NDOFEL()  :Nodal degrees of freedom per element                  C
C     MATP()    :Number of material profile at each element            C
C     NMATP()   :Number of Real material parameters for each profile   C                                     
C     NMPR      :Number of real material profiles                      C
C     NIMTP()   :Number of integer material parameters                 C
C     NIPR      :Number of integer material profiles                   C
C     AMATE(,)  :Real material properties array                        C
C     IMTEI(,)  :Integer material properties array                     C
C     AWAVE()   :Wave field parameters                                 C
C     IWAVE     :Wave type identifier                                  C
c     COORD(,)  :Nodal points coordinates                              C
C     LM(,)     :Global DME operator                                   C
C     ID()      :Equation identifier at each node                      C
C     IDPL      :Point load index                                      C    
C     NEQ       :Number of equations                                   C        
C     DT        :Time step                                             C               
C     NINCR     :Total number of increments                            C
C     KINC      :Current increment                                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE NODSTFFASEM(UTMINT,UT,UTMAST,NUMNP,NUMEL,NUMAT,NNE,
     1                       NMNE,MXDOFDIM,MXDOFEL,IELT,IELCON,ILIST,
     2                       LPLIST,NIEL,NDOFN,NDOFEL,MATP,NMATP,NMPR,
     3                       NIMTP,NIPR,AMATE,IMTEI,AWAVE,IWAVE,COORD,
     4                       LM,ID,IDPL,NPL,NEQ,IOUT,
     5                       DT,NINCR,KINC,NTH)

      USE OMP_LIB
      IMPLICIT REAL*8 (A-H,O-Z)

C     GLOBAL REAL ARRAYS.

      DIMENSION COORD(MXDOFDIM,NUMNP),AMATE(NMPR,NUMAT),AWAVE(5),
     1          MATP(NUMEL),UTMINT(NEQ),UT(NEQ),UTMAST(NEQ),
     2          AG(NEQ),VG(NEQ),VLOADS(NPL)

C     INTEGER ARRAYS

      DIMENSION NDOFN(NUMNP),ID(MXDOFDIM,NUMNP),NMATP(NUMAT),
     1          NIMTP(NUMAT),IMTEI(NIPR,NUMAT),NNE(NUMEL),
     2          NDOFEL(NUMEL),IELCON(NMNE,NUMEL),LM(MXDOFEL,NUMEL),
     3          ILIST(NUMNP,6),LPLIST(NUMNP,6),NIEL(NUMNP),
     4          IELT(NUMEL)

      DIMENSION IDPL(NPL)

      IF (NPL.NE.0) THEN
        DO I=1, NPL
              READ(9999,*) IJKIJK, VLOADS(I)
        END DO
      END IF

C$ CALL OMP_SET_NUM_THREADS(NTH)

C$OMP PARALLEL

C$OMP DO PRIVATE(K)

      DO I=1,NUMNP          !Loops through each node.
        K=NIEL(I)

        CALL ELEASSEM(I,K,ID,UTMINT,UT,UTMAST,NUMNP,NUMEL,NUMAT,
     1                NNE,NMNE,MXDOFDIM,MXDOFEL,IELT,IELCON,ILIST,
     2                LPLIST,NIEL,NDOFN,NDOFEL,MATP,NMATP,NMPR,
     3                NIMTP,NIPR,AMATE,IMTEI,AWAVE,IWAVE,COORD,LM,
     4                IDPL,VLOADS,NPL,NEQ,IOUT,DT,NINCR,KINC)

      END DO

C$OMP END DO

C$OMP END PARALLEL

      UTMINT=UT
      UT=UTMAST

      RETURN

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE ELEASSEM                                              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE ELEASSEM(INOD,NELE,ID,UTMINT,UT,UTMAST,NUMNP,NUMEL,
     1                    NUMAT,NNE,NMNE,MXDOFDIM,MXDOFEL,IELT,IELCON,
     2                    ILIST,LPLIST,NIEL,NDOFN,NDOFEL,MATP,NMATP,
     3                    NMPR,NIMTP,NIPR,AMATE,IMTEI,AWAVE,IWAVE,COORD,
     4                    LM,IDPL,VLOADS,NPL,NEQ,IOUT,DT,NINCR,KINC)

      IMPLICIT REAL*8(A-H,O-Z)

      ALLOCATABLE SM(:),UL(:),VL(:),DUL(:),AL(:),R(:)

C     REAL ARRAYS.

      DIMENSION COORD(MXDOFDIM,NUMNP),AMATE(NMPR,NUMAT),AWAVE(5),
     1          MATP(NUMEL),ELCOOR(MXDOFDIM,NMNE),UTMINT(NEQ),
     2          UT(NEQ),UTMAST(NEQ),AG(NEQ),VG(NEQ),PROPS(NMPR),
     3          AMASS(MXDOFDIM), FORCE(MXDOFDIM),VLOADS(NPL)

C     INTEGER ARRAYS

      DIMENSION NDOFN(NUMNP),ID(MXDOFDIM,NUMNP),NMATP(NUMAT),
     1          NIMTP(NUMAT),IMTEI(NIPR,NUMAT),NNE(NUMEL),
     2          NDOFEL(NUMEL),IELCON(NMNE,NUMEL),LM(MXDOFEL,NUMEL),
     3          ILIST(NUMNP,6),LPLIST(NUMNP,6),NIEL(NUMNP),
     4          LML(MXDOFEL),IPROPS(NIPR),IELT(NUMEL),IDPL(NPL)

      I=INOD
      K=NELE
      AMASS=0.0D0
      FORCE=0.0D0

      DO J=1,K             !Loops through all the elements at the current node I
        IE=ILIST(I,J)
        IM=LPLIST(I,J)
        ND=NDOFEL(IE)
        NN=NNE(IE)
        NPROPS =NMATP(MATP(IE))
        NIPROPS=NIMTP(MATP(IE))
        ALLOCATE(R(ND),SM(ND),UL(ND),DUL(ND),VL(ND),AL(ND))
        CALL CLEARV(SM,ND)
        CALL CLEARV(UL,ND)
        CALL CLEARV(VL,ND)
        CALL CLEARV(R,ND)
        CALL CLEARV(DUL,ND)
        CALL CLEARV(AL,ND)
        CALL CLEARIV(IPROPS,NIPR)

        DO K1=1,NPROPS
          PROPS(K1)=AMATE(K1,MATP(IE))
        END DO

        IF(NIPROPS.NE.0) THEN
          DO K1=1,NIPROPS
            IPROPS(K1)=IMTEI(K1,MATP(IE))
          END DO
        END IF

C       RETRIEVES ELEMENT COORDINATES,
C       LOCAL ASSEMBLY LIST AND ELEMENTAL
C       DISPLACEMENTS AND VELOCITIES.

        K3=0
        DO JI=1,NNE(IE)
          K1=NDOFN(IELCON(JI,IE))
          DO JL=1,MXDOFDIM
            ELCOOR(JL,JI)=COORD(JL,IELCON(JI,IE))
          END DO

          DO K2=1,K1
            LML(K3+K2)=LM(K3+K2,IE)
          END DO
          K3=K3+K1
        END DO

        DO K1=1,ND
          IF(LML(K1).GT.0) THEN
            UL(K1) = UT(LML(K1))
            DUL(K1)= UT(LML(K1))-UTMINT(LML(K1))
            VL(K1) = (UT(LML(K1))-UTMINT(LML(K1)))/DT           ! Predictor
            AL(K1) = (1.0D0/DT/DT)
     1      *(UTMINT(LML(K1))-2.0D0*UT(LML(K1)))                ! Predictor
        ELSE
            UL(K1) =0.0D0
            DUL(K1)=0.0D0
            VL(K1 )=0.0D0
            AL(K1 )=0.0D0
          END IF
        END DO
C
C       CALLS USER ELEMENT SUBROUTINE UEL ACCORDING TO IELT(I)
C
        SELECT CASE (IELT(IE))
C
          CASE(1)
          CALL UELEXP8(R,SM,ND,PROPS,NPROPS,IPROPS,NIPROPS,
     1                 AWAVE,IWAVE,ELCOOR,MXDOFDIM,NN,UL,DUL,VL,
     2                 AL,IE,DT,KINC,NINCR)

          CASE (2)

          CALL UELEXP4(R,SM,ND,PROPS,NPROPS,IPROPS,NIPROPS,
     1                 AWAVE,IWAVE,ELCOOR,MXDOFDIM,NN,UL,DUL,VL,
     2                 AL,IE,DT,KNINC,NINCR)

          CASE (3)
          CALL UELEXP9(R,SM,ND,PROPS,NPROPS,IPROPS,NIPROPS,
     1                 AWAVE,IWAVE,ELCOOR,MXDOFDIM,NN,UL,DUL,VL,
     2                 AL,IE,DT,KINC,NINCR)

          CASE (4)
          CALL UEL8INCOT(R,SM,ND,PROPS,NPROPS,IPROPS,NIPROPS,
     1                   AWAVE,IWAVE,ELCOOR,MXDOFDIM,NN,UL,DUL,VL,
     2                   AL,IE,DT,KINC,NINCR)

          CASE (5)
          CALL UEL9INCOT(R,SM,ND,PROPS,NPROPS,IPROPS,NIPROPS,
     1                 AWAVE,IWAVE,ELCOOR,MXDOFDIM,NN,UL,DUL,VL,
     2                 AL,IE,DT,KINC,NINCR)

          CASE (7)
          CALL UELEDASH3(R,SM,ND,PROPS,NPROPS,IPROPS,NIPROPS,
     1                  AWAVE,IWAVE,ELCOOR,MXDOFDIM,NN,UL,DUL,VL,
     2                  AL,IE,DT,KINC,NINCR)

          CASE (8)
          CALL UELDASHINC(R,SM,2*NN,PROPS,NPROPS,IPROPS,NIPROPS,
     1                    AWAVE,IWAVE,ELCOOR,MXDOFDIM,NN,UL,DUL,VL,
     2                    AL,IE,DT,KINC,NINCR)

        END SELECT

C       RETRIEVES EQN NUMBER AT THE CURRENT NODE OF THE
C       CURRENT ELEMENT

        IJJ=LM(2*IM-1,IE)
        
        DO IK=1, MXDOFDIM
            JK=ID(IK,I)
            
            IF(JK.EQ.IJJ) THEN
              II=2*IM-1
            ELSE
              II=2*IM
            END IF

C       ASSEMBLES MASS AND INERTIAL FORCE

            AMASS(IK)=AMASS(IK)+SM(II)
            FORCE(IK)=FORCE(IK)+R(II )

        END DO
            
        DEALLOCATE(R,SM,UL,DUL,VL, AL)

      END DO                            !Ends loop through the elements.

C     SOLVES FOR THE CURRENT DOF ACCORDING TO THE CURRENT INCREMENT.

      IF(NPL.NE.0)THEN
        DO IK=1, MXDOFDIM
            JK=ID(IK,I)

            DO IJUAN=1, NPL
                IF(IDPL(IJUAN)==JK) THEN
cc                    CALL UPULSE(PP,NINCR,DT,AWAVE,KINC)
                    FORCE(IK)=FORCE(IK)+VLOADS(IJUAN)
                    EXIT
                END IF
            END DO
 
        END DO
      END IF

      DO IK=1, MXDOFDIM
        JK=ID(IK,I)
            
        IF (JK .NE. 0) THEN
            UTMAST(JK)=DT*DT*FORCE(IK)/AMASS(IK)
        END IF
            
      END DO
      
      RETURN

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE ASSEMLIS                                              C
C                                                                      C
C     Creates the DME() operator LM()                                  C
C                                                                      C
C     NUMNP    :Number of nodal points                                 C
C     NUMEL    :Number of elements                                     C
C     MXNE     :Maximum number of nodes per element                    C
C     MXDOFDIM :Maximum degree of freedom dimension                    C
C     NNE()    :Number of nodes at the current element                 C
C     NDOFN()  :Number of degrees of freedom at the current node       C
C     IELCON(,):Nodal connectivity array                               C
C     LM(,)    :DME opeartor                                           C               
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE ASSEMLIS(NUMNP,NUMEL,MXNE,MXDOFDIM,NNE,NDOFN,
     1                    IELCON,LM,ID,IIN,IOUT,NTH)

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION NNE(NUMEL),IELCON(MXNE,NUMEL),LM(MXDOFDIM*MXNE,NUMEL),
     1          ID(MXDOFDIM,NUMNP),NDOFN(NUMNP)

C$ CALL OMP_SET_NUM_THREADS(NTH)
C$OMP PARALLEL
C$OMP DO PRIVATE(K3, J, K2I, K2)
      DO I=1,NUMEL
        K3=0
        DO J=1,NNE(I)
          K2I=NDOFN(IELCON(J,I))
          DO K2=1,K2I
            LM(K3+K2,I)=ID(K2,IELCON(J,I))
          END DO
          K3=K3+K2I
        END DO
      END DO
C$OMP END DO
C$OMP END PARALLEL

      RETURN

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE INCLIST                                               C
C                                                                      C
C     SEARCHING THE DME OPERATOR CREATES INCIDENCE LIST ILIST()        C
C     AND LOCAL POSITION LIST LPLIST()                                 C
C                                                                      C
C     NUMNP     :Number of nodal points                                C
C     NUMEL     :Number of elements                                    C
C     MXNE      :Maximum number of nodes per element                   C
C     IELCON(,) :Nodal connectivity array                              C
C     ILIST()   :Incidence list.Elements connected to the current node C
C     LPLIST()  :Local position of the current node in each one of the C
C                elements of ILIST()                                   C
C     NIEL()    :Number of elements at the current node                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE INCLIST(NUMNP,NUMEL,MXNE,IELCON,ILIST,LPLIST,NIEL,
     1                   IIN,IOUT,NTH)

      USE OMP_LIB
      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION ILIST(NUMNP,6),IELCON(MXNE,NUMEL),LIST(MXNE),
     1          LPLIST(NUMNP,6),NIEL(NUMNP)

      CALL CLEARIM(ILIST, NUMNP,6)
      CALL CLEARIM(LPLIST,NUMNP,6)
      CALL CLEARIV(NIEL,NUMNP)

C$ CALL OMP_SET_NUM_THREADS(NTH)
C$OMP PARALLEL
C$OMP DO PRIVATE(IK, J, K, LIST, IFLAG, IPOS)

      DO I=1,NUMNP
        IK=1
        DO J=1,NUMEL
          DO K=1,MXNE
            LIST(K)=IELCON(K,J)
          END DO
          CALL SEARCHPOS(LIST,MXNE,I,IFLAG,IPOS)
          IF(IFLAG.EQ.1) THEN
            ILIST(I,IK)=J
            LPLIST(I,IK)=IPOS
            IK=IK+1
          END IF
        END DO
        NIEL(I)=IK-1
      END DO

C$OMP END DO
C$OMP END PARALLEL

      RETURN

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SEARCHPOS                                             C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE SEARCHPOS(LIST,NDAT,NVAL,IFLAG,IPOS)

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION LIST(NDAT)

      IFLAG=0
      DO I=1,NDAT
         IF(LIST(I).EQ.NVAL) THEN
            IPOS=I
            IFLAG=1
         END IF
      END DO

      RETURN

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE BELONGS                                               C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE BELONGS(LIST,NDAT,IVAL,IFLAG)

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION LIST(NDAT)

      IFLAG=0
      DO I=1,NDAT
        IF(IVAL.EQ.LIST(I)) THEN
          IFLAG=1
        END IF
      END DO

      RETURN

      END



C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE NODSTFFASEM                                           C
C                                                                      C
C     This is the actual solver at each increment.                     C
C                                                                      C
C     UG(,)     :Global displacements history at each node             C
C     VG(,)     :Global velocity histories at each node                C
C     AG(,)     :Global acceleration history at each node              C
C     NUMNP     :Number of nodal points                                C
C     NUMEL     :Number of elements                                    C
C     NUMAT     :Number of material profiles                           C
C     NNE()     :Number of nodes per element                           C
C     NMNE      :Maximum number of nodes per element                   C
C     MXDOFDIM  :Maximum degree of freedom dimension                   C
C     MXDOFEL   :Maximum degree of freedom per element                 C
C     IELT()    :Array of element type identifiers                     C
C     IELCON(,) :Elementary connectivity array                         C
C     ILIST()   :Incidence list.Elements connected to the current node C
C     LPLIST()  :Local position of the current node in each one of the C
C                elements of ILIST()                                   C
C     NIEL()    :Number of elements at the current node                C
C     NDOFN()   :Number of degrees of freedom per node                 C
C     NDOFEL()  :Nodal degrees of freedom per element                  C
C     MATP()    :Number of material profile at each element            C
C     NMATP()   :Number of Real material parameters for each profile   C                                     
C     NMPR      :Number of real material profiles                      C
C     NIMTP()   :Number of integer material parameters                 C
C     NIPR      :Number of integer material profiles                   C
C     AMATE(,)  :Real material properties array                        C
C     IMTEI(,)  :Integer material properties array                     C
C     AWAVE()   :Wave field parameters                                 C
C     IWAVE     :Wave type identifier                                  C
c     COORD(,)  :Nodal points coordinates                              C
C     LM(,)     :Global DME operator                                   C
C     ID()      :Equation identifier at each node                      C
C     IDPL      :Point load index                                      C    
C     NEQ       :Number of equations                                   C        
C     DT        :Time step                                             C               
C     NINCR     :Total number of increments                            C
C     KINC      :Current increment                                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE NODSTFFASEMINI(UT,NUMNP,NUMEL,NUMAT,NNE,
     1                       NMNE,MXDOFDIM,MXDOFEL,IELT,IELCON,ILIST,
     2                       LPLIST,NIEL,NDOFN,NDOFEL,MATP,NMATP,NMPR,
     3                       NIMTP,NIPR,AMATE,IMTEI,AWAVE,IWAVE,COORD,
     4                       LM,ID,IDPL,NPL,NEQ,IOUT,
     5                       DT,NINCR,KINC,AG,VG,NTH)

      USE OMP_LIB
      IMPLICIT REAL*8 (A-H,O-Z)

C     GLOBAL REAL ARRAYS.

      DIMENSION COORD(MXDOFDIM,NUMNP),AMATE(NMPR,NUMAT),AWAVE(5),
     1          MATP(NUMEL),UT(NEQ),
     2          AG(NEQ),VG(NEQ),VLOADS(NPL)

C     INTEGER ARRAYS

      DIMENSION NDOFN(NUMNP),ID(MXDOFDIM,NUMNP),NMATP(NUMAT),
     1          NIMTP(NUMAT),IMTEI(NIPR,NUMAT),NNE(NUMEL),
     2          NDOFEL(NUMEL),IELCON(NMNE,NUMEL),LM(MXDOFEL,NUMEL),
     3          ILIST(NUMNP,6),LPLIST(NUMNP,6),NIEL(NUMNP),
     4          IELT(NUMEL)

      DIMENSION IDPL(NPL)

      IF (NPL.NE.0) THEN
        DO I=1, NPL
              READ(9999,*) IJKIJK, VLOADS(I)
        END DO
      END IF

C$ CALL OMP_SET_NUM_THREADS(NTH)

C$OMP PARALLEL

C$OMP DO PRIVATE(K)

      DO I=1,NUMNP          !Loops through each node.
        K=NIEL(I)

        CALL ELEASSEMINI(I,K,ID,UT,NUMNP,NUMEL,NUMAT,
     1                NNE,NMNE,MXDOFDIM,MXDOFEL,IELT,IELCON,ILIST,
     2                LPLIST,NIEL,NDOFN,NDOFEL,MATP,NMATP,NMPR,
     3                NIMTP,NIPR,AMATE,IMTEI,AWAVE,IWAVE,COORD,LM,
     4                IDPL,VLOADS,NPL,NEQ,IOUT,DT,NINCR,KINC,AG,VG)

      END DO

C$OMP END DO

C$OMP END PARALLEL
      
      IF (NPL.NE.0) THEN
        CLOSE(9999)
        OPEN(9999,FILE ='Loads.txt')
      END IF
      
      RETURN

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE ELEASSEM                                              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE ELEASSEMINI(INOD,NELE,ID,UT,NUMNP,NUMEL,
     1                    NUMAT,NNE,NMNE,MXDOFDIM,MXDOFEL,IELT,IELCON,
     2                    ILIST,LPLIST,NIEL,NDOFN,NDOFEL,MATP,NMATP,
     3                    NMPR,NIMTP,NIPR,AMATE,IMTEI,AWAVE,IWAVE,COORD,
     4                    LM,IDPL,VLOADS,NPL,NEQ,IOUT,DT,NINCR,KINC,
     5                    AG,VG)

      IMPLICIT REAL*8(A-H,O-Z)

      ALLOCATABLE SM(:),UL(:),VL(:),DUL(:),AL(:),R(:)

C     REAL ARRAYS.

      DIMENSION COORD(MXDOFDIM,NUMNP),AMATE(NMPR,NUMAT),AWAVE(5),
     1          MATP(NUMEL),ELCOOR(MXDOFDIM,NMNE),
     2          UT(NEQ),AG(NEQ),VG(NEQ),PROPS(NMPR),VLOADS(NPL),
     3          AMASS(MXDOFDIM), FORCE(MXDOFDIM)

C     INTEGER ARRAYS

      DIMENSION NDOFN(NUMNP),ID(MXDOFDIM,NUMNP),NMATP(NUMAT),
     1          NIMTP(NUMAT),IMTEI(NIPR,NUMAT),NNE(NUMEL),
     2          NDOFEL(NUMEL),IELCON(NMNE,NUMEL),LM(MXDOFEL,NUMEL),
     3          ILIST(NUMNP,6),LPLIST(NUMNP,6),NIEL(NUMNP),
     4          LML(MXDOFEL),IPROPS(NIPR),IELT(NUMEL),IDPL(NPL)

      I=INOD
      K=NELE
      AMASS=0.0D0
      FORCE=0.0D0

      DO J=1,K             !Loops through all the elements at the current node I
        IE=ILIST(I,J)
        IM=LPLIST(I,J)
        ND=NDOFEL(IE)
        NN=NNE(IE)
        NPROPS =NMATP(MATP(IE))
        NIPROPS=NIMTP(MATP(IE))
        PPV=0.0D0
        ALLOCATE(R(ND),SM(ND),UL(ND),DUL(ND),VL(ND),AL(ND))
        CALL CLEARV(SM,ND)
        CALL CLEARV(UL,ND)
        CALL CLEARV(VL,ND)
        CALL CLEARV(R,ND)
        CALL CLEARV(DUL,ND)
        CALL CLEARV(AL,ND)
        CALL CLEARIV(IPROPS,NIPR)

        DO K1=1,NPROPS
          PROPS(K1)=AMATE(K1,MATP(IE))
        END DO

        IF(NIPROPS.NE.0) THEN
          DO K1=1,NIPROPS
            IPROPS(K1)=IMTEI(K1,MATP(IE))
          END DO
        END IF

C       RETRIEVES ELEMENT COORDINATES,
C       LOCAL ASSEMBLY LIST AND ELEMENTAL
C       DISPLACEMENTS AND VELOCITIES.

        K3=0
        DO JI=1,NNE(IE)
          K1=NDOFN(IELCON(JI,IE))
          DO JL=1,MXDOFDIM
            ELCOOR(JL,JI)=COORD(JL,IELCON(JI,IE))
          END DO

          DO K2=1,K1
            LML(K3+K2)=LM(K3+K2,IE)
          END DO
          K3=K3+K1
        END DO

        DO K1=1,ND
          IF(LML(K1).GT.0) THEN
            AL(K1)=0.0D0
            UL(K1)=UT(LML(K1))
            VL(K1)=VG(LML(K1))
          ELSE
            UL(K1) =0.0D0
            DUL(K1)=0.0D0
            VL(K1 )=0.0D0
            AL(K1 )=0.0D0
          END IF
        END DO
C
C       CALLS USER ELEMENT SUBROUTINE UEL ACCORDING TO IELT(I)
C
        SELECT CASE (IELT(IE))
C
          CASE(1)
          CALL UELEXP8(R,SM,ND,PROPS,NPROPS,IPROPS,NIPROPS,
     1                 AWAVE,IWAVE,ELCOOR,MXDOFDIM,NN,UL,DUL,VL,
     2                 AL,IE,DT,KINC,NINCR)

          CASE (2)

          CALL UELEXP4(R,SM,ND,PROPS,NPROPS,IPROPS,NIPROPS,
     1                 AWAVE,IWAVE,ELCOOR,MXDOFDIM,NN,UL,DUL,VL,
     2                 AL,IE,DT,KNINC,NINCR)

          CASE (3)
          CALL UELEXP9(R,SM,ND,PROPS,NPROPS,IPROPS,NIPROPS,
     1                 AWAVE,IWAVE,ELCOOR,MXDOFDIM,NN,UL,DUL,VL,
     2                 AL,IE,DT,KINC,NINCR)

          CASE (4)
          CALL UEL8INCOT(R,SM,ND,PROPS,NPROPS,IPROPS,NIPROPS,
     1                   AWAVE,IWAVE,ELCOOR,MXDOFDIM,NN,UL,DUL,VL,
     2                   AL,IE,DT,KINC,NINCR)

          CASE (5)
          CALL UEL9INCOT(R,SM,ND,PROPS,NPROPS,IPROPS,NIPROPS,
     1                 AWAVE,IWAVE,ELCOOR,MXDOFDIM,NN,UL,DUL,VL,
     2                 AL,IE,DT,KINC,NINCR)

          CASE (7)
          CALL UELEDASH3(R,SM,ND,PROPS,NPROPS,IPROPS,NIPROPS,
     1                  AWAVE,IWAVE,ELCOOR,MXDOFDIM,NN,UL,DUL,VL,
     2                  AL,IE,DT,KINC,NINCR)

          CASE (8)
          CALL UELDASHINC(R,SM,2*NN,PROPS,NPROPS,IPROPS,NIPROPS,
     1                    AWAVE,IWAVE,ELCOOR,MXDOFDIM,NN,UL,DUL,VL,
     2                    AL,IE,DT,KINC,NINCR)

        END SELECT

C       RETRIEVES EQN NUMBER AT THE CURRENT NODE OF THE
C       CURRENT ELEMENT

        IJJ=LM(2*IM-1,IE)
        
        DO IK=1, MXDOFDIM
            JK=ID(IK,I)
            
            IF(JK.EQ.IJJ) THEN
              II=2*IM-1
            ELSE
              II=2*IM
            END IF

C       ASSEMBLES MASS AND INERTIAL FORCE

            AMASS(IK)=AMASS(IK)+SM(II)
            FORCE(IK)=FORCE(IK)+R(II )

        END DO
            
        DEALLOCATE(R,SM,UL,DUL,VL, AL)

      END DO                            !Ends loop through the elements.

C     SOLVES FOR THE CURRENT DOF ACCORDING TO THE CURRENT INCREMENT.

      IF(NPL.NE.0)THEN
        DO IK=1, MXDOFDIM
            JK=ID(IK,I)

            DO IJUAN=1, NPL
                IF(IDPL(IJUAN)==JK) THEN
cc                    CALL UPULSE(PP,NINCR,DT,AWAVE,KINC)
                    FORCE(IK)=FORCE(IK)+VLOADS(IJUAN)
                    EXIT
                END IF
            END DO

        END DO
      END IF

      DO IK=1, MXDOFDIM
        JK=ID(IK,I)

        IF (JK .NE. 0) THEN
            AG(JK)= FORCE(IK)/AMASS(IK)
        END IF

      END DO

      RETURN

      END
cc      END MODULE ASEMUTIL
