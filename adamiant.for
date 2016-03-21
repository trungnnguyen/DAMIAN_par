CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   --------P R O G R A M   D A M I A N   E X P L I C I T------------  C
C                       Finite Element Method                          C
C                                                                      C
C   EXPLICIT TIME DOMAIN PLANE STRAIN ANALYSIS FOR FINITE DOMAINS.     C
C                                                                      C
C  DYNAMIC MEMORY ALLOCATION                                           C
C  EAFIT UNIVERSITY                                                    C
C  APPLIED MECHANICS LAB                                               C
C  FEBRUARY 25/2015                                                    C
C                                                                      C
C  UNIX VERSION                                                        C
C                                                                      C
C            P R O B L E M   P A R A M E T E R S                       C
C                                                                      C
C  NUMNP     :NUMBER OF NODAL POINTS                                   C
C  NUMEL     :UMBER OF ELEMENTS                                        C
C  NUMAT     :NUMBER OF MATERIAL PROFILES                              C
C  TM        :SIZE OF THE TIME WINDOW.                                 C
C  NINCR     :NUMBER OF INCREMENTS                                     C
C  NDOFDIM   :PROBLEM DIMENSIONALITY                                   C
C  NMNE      :MAXIMUM NUMBER OF NODES PER ELEMENT                      C
C  NMDOFEL   :MAXIMUM NUMBER OF DOF PER ELEMENT                        C
C  NMPR      :MAXIMUM NUMBER OF MATERIAL PROPERTIES IN A PROFILE       C
C  NIPR      :MAXIMUM NUMBER OF INTEGER MAT PROPERTIES IN A PROFILE    C
C  NPL       :NUMBER OF POINT LOADS                                    C
C  NUMPARA   :NUMBER OF ELEMENTS FOR PARAVIEW                          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
***********************************************************************C
C                                                                      C
C                      MAIN PROGRAM STARTS                             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      PROGRAM DAMIAN_EXPLICIT

      USE OMP_LIB
      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION AWAVE(5)
C
      ALLOCATABLE ID(:,:),NDOFN(:),COORD(:,:),MATP(:),NMATP(:),NIMTP(:),
     1            AMATE(:,:),IMTEI(:,:),IELT(:),NNE(:),NDOFEL(:),
     2            IELCON(:,:),LM(:,:),ILIST(:,:),LPLIST(:,:),NIEL(:),
     3            AG(:),VG(:),UTMINT(:),UT(:),UTMAST(:)

      ALLOCATABLE IDPL(:)

      CHARACTER*80 HED
      CHARACTER*10 FILENAME
      CHARACTER*15 FILELOADS
      
      
cccccccccccccccccccccccccccccccccccccccccc
      REAL*8 TIME, ELAPSED_TIME
      INTEGER TCLOCK1, TCLOCK2, CLOCK_RATE
cccccccccccccccccccccccccccccccccccccccccc

      NTH=4

C     *****************************************************************C
C     ***          P R O B L E M  F I L E S                          **C
C     *****************************************************************C

      IIN=5
      IOUT=4
      IDIS=7
      WRITE(*,*) 'INPUT THE JOB NAME(max 10 characters):'
      READ(*,*) FILENAME
cc      FILENAME="barra"
      LST=LEN_TRIM(FILENAME)
      OPEN(UNIT=IIN,FILE =FILENAME(1:LST)//".inp",FORM='FORMATTED')
      OPEN(UNIT=IOUT,FILE=FILENAME(1:LST)//".dat",FORM='FORMATTED')
      OPEN(UNIT=IDIS,FILE=FILENAME(1:LST)//".dis",FORM='FORMATTED')

C     *****************************************************************C
C     ***                 I N P U T   P H A S E                      **C
C     *****************************************************************C

C     READS PROBLEM DEFINITION PARAMETERS.

      READ(IIN,*,IOSTAT=INOUTSTATUS) HED
      IF(INOUTSTATUS.LT.0) STOP "***COULD NOT OPEN FILE"

      CALL SYSTEM_CLOCK(TCLOCK1)
      CALL CPU_TIME(TIME)
      WRITE(*,'(1A23, 1F12.2)') 'Begins the analysis: ', TIME
      WRITE(*,*)

      READ(IIN,     *) NUMNP,NUMEL,NUMAT,TM,NINCR,NDOFDIM,NMNE,NMDOFEL,
     1                 NMPR,NIPR,NPL
     
      IF (NPL.NE.0) THEN
        OPEN(9999,FILE ='Loads.txt')
      END IF

      WRITE(IOUT,1900) HED
      WRITE(IOUT,2000) NUMNP,NUMEL,NUMAT,TM,NINCR,NDOFDIM,NMNE,NMDOFEL,
     1                 NMPR,NIPR,NPL

      ALLOCATE(ID(NDOFDIM,NUMNP),NDOFN(NUMNP),COORD(NDOFDIM,NUMNP),
     1         MATP(NUMEL),NMATP(NUMAT),NIMTP(NUMAT),AMATE(NMPR,NUMAT),
     2         IMTEI(NIPR,NUMAT),IELT(NUMEL),NNE(NUMEL),
     3         NDOFEL(NUMEL),IELCON(NMNE,NUMEL),LM(NMDOFEL,NUMEL),
     4         ILIST(NUMNP,6),LPLIST(NUMNP,6),NIEL(NUMNP))

      ALLOCATE(IDPL(NPL))

C$ CALL OMP_SET_NUM_THREADS(NTH)
C$OMP PARALLEL

C$OMP SINGLE
      CALL CLEARIV(IDPL, NPL)

      CALL CLEARIM(ID,NDOFDIM,NUMNP)
      CALL CLEARIV(NDOFN,NUMNP)
      CALL CLEAR(COORD,NDOFDIM,NUMNP)

      CALL NODINP(NUMNP,ID,NDOFDIM,COORD,NDOFN,NEQ,IIN,IOUT)
C$OMP END SINGLE

C$OMP SECTIONS

C     SOLUTION ARRAYS.

C$OMP SECTION
      ALLOCATE(AG(NEQ),VG(NEQ))

C     CALL CLEARV(AG,NEQ)                   !FOR INITIALIZATION
C     CALL CLEARV(VG,NEQ)
      
      AG=0.0D0
      VG=0.0D0

C$OMP SECTION
      ALLOCATE(UTMINT(NEQ),UT(NEQ))
      UT    =0.0D0
      UTMINT=0.0D0
      
C     CLEARS STORAGE.

C$OMP SECTION
      CALL CLEARIV(MATP,NUMEL)
      CALL CLEARIV(NMATP,NUMAT)
      CALL CLEARIV(NIMTP,NUMAT)
      CALL CLEAR(AMATE,NMPR,NUMAT)
      CALL CLEARIM(IMTEI,NIPR,NUMAT)
      CALL CLEARIV(IELT,NUMEL)
      CALL CLEARIV(NNE,NUMEL)
      CALL CLEARIV(NDOFEL,NUMEL)
      CALL CLEARIM(IELCON,NMNE,NUMEL)
      CALL CLEARIM(LM,NMDOFEL,NUMEL)
      CALL CLEARIV(NIEL,NUMNP)

C     READS MODEL

      CALL MATINP(NUMAT,NMPR,NMATP,AMATE,NIPR,NIMTP,IMTEI,AWAVE,IWAVE,
     1            IIN,IOUT)
      CALL ELEINP(NUMNP,NUMEL,NNE,IELT,NDOFEL,NMNE,MATP,IELCON,
     1              NUMPARA,IIN,IOUT)
      CALL POINTLOAD(NUMNP,NDOFDIM,NPL,ID,IDPL,IIN,FILELOADS)

C$OMP END SECTIONS
C$OMP END PARALLEL

C     CREATES ASSEMBLY LIST (DME OPERATOR) AND LIST
C     OF INCIDENT ELEMENTS FOR EACH NODE.

      CALL ASSEMLIS(NUMNP,NUMEL,NMNE,NDOFDIM,NNE,NDOFN,IELCON,
     1              LM,ID,IIN,IOUT,NTH)

      CALL INCLIST(NUMNP,NUMEL,NMNE,IELCON,ILIST,LPLIST,NIEL,IIN,IOUT,
     1             NTH)

C     *****************************************************************C
C     ***            I N C    S L N   P H A S E                      **C
C     *****************************************************************C

      DT=TM/(NINCR-1)

C     INITIALIZE ACCELERATION
C     KINC=2 equivalent to initial conditions

      KINC=2
      IFLAG=0
     
C     CALL NODSTFFASEM(UTMINT,UT,UTMAST,NUMNP,NUMEL,NUMAT,NNE,NMNE,
C    1                 NDOFDIM,NMDOFEL,IELT,IELCON,ILIST,LPLIST,
C    2                 NIEL,NDOFN,NDOFEL,MATP,NMATP,NMPR,NIMTP,
C    3                 NIPR,AMATE,IMTEI,AWAVE,IWAVE,COORD,LM,ID,
C    4                 IDPL,NPL,NEQ,IOUT,DT,NINCR,KINC,IFLAG,AG,VG)

      CALL NODSTFFASEMINI(UT,NUMNP,NUMEL,NUMAT,NNE,NMNE,
     1                 NDOFDIM,NMDOFEL,IELT,IELCON,ILIST,LPLIST,
     2                 NIEL,NDOFN,NDOFEL,MATP,NMATP,NMPR,NIMTP,
     3                 NIPR,AMATE,IMTEI,AWAVE,IWAVE,COORD,LM,ID,
     4                 IDPL,NPL,NEQ,IOUT,DT,NINCR,KINC,AG,VG,NTH)

C$ CALL OMP_SET_NUM_THREADS(NTH)
C$OMP PARALLEL
C$OMP DO
      DO I=1, NEQ
        UTMINT(I)=UT(I)-DT*VG(I)+(DT*DT/2.0D0)*AG(I)
      END DO
C$OMP END DO

C$OMP SECTIONS

C$OMP SECTION
      ALLOCATE(UTMAST(NEQ))
      UTMAST=0.0D0

C$OMP SECTION
      DEALLOCATE(AG)

C$OMP SECTION
      DEALLOCATE(VG)

C$OMP END SECTIONS

C$OMP END PARALLEL

C     INCREMENTATION BEGINS

      III=0
      JJJ=1

      IFLAG=1

      DO KINC=2,NINCR-1
        
        CALL CPU_TIME(TIME)
        
        WRITE(*,'(1A12, 1I10, 1A8, 1F82.2)') 'INCREMENT', KINC, 
     1                                       'Time: ', TIME

        CALL NODSTFFASEM(UTMINT,UT,UTMAST,NUMNP,NUMEL,NUMAT,NNE,NMNE,
     1                   NDOFDIM,NMDOFEL,IELT,IELCON,ILIST,LPLIST,
     2                   NIEL,NDOFN,NDOFEL,MATP,NMATP,NMPR,NIMTP,
     3                   NIPR,AMATE,IMTEI,AWAVE,IWAVE,COORD,LM,ID,
     4                   IDPL,NPL,NEQ,IOUT,DT,NINCR,KINC,NTH)

        III=III+1
        IF (III==10) THEN

            IF (NMNE==4) THEN
                CALL VTK_4NODES(NUMNP, NUMPARA, NEQ, NDOFDIM, JJJ,
     1                      COORD,IELCON(:,1:NUMPARA),ID,UTMAST)
            END IF

            IF (NMNE==8) THEN
                CALL VTK_8NODES(NUMNP, NUMPARA, NEQ, NDOFDIM, JJJ,
     1                      COORD,IELCON(:,1:NUMPARA),ID,UTMAST)
            END IF

            IF (NMNE==9) THEN
                CALL VTK_9NODES(NUMNP, NUMPARA, NEQ, NDOFDIM, JJJ,
     1                      COORD,IELCON(:,1:NUMPARA),ID,UTMAST)
            END IF
     
            JJJ=JJJ+1
            III=0

        END IF

      END DO

      WRITE(*,3050), FILENAME
      
      IF (NPL.NE.0) THEN
        CLOSE(9999)
      END IF

      WRITE(*,*)
      CALL SYSTEM_CLOCK(TCLOCK2, CLOCK_RATE)
      ELAPSED_TIME = FLOAT(TCLOCK2 - TCLOCK1) / FLOAT(CLOCK_RATE)
      PRINT 11, ELAPSED_TIME
   11 FORMAT("Elapsed time = ",f12.2, " seconds")

      STOP

C     *****************************************************************C
C     ***         S L N  O U T P U T  P H A S E                      **C
C     *****************************************************************C

 1900 FORMAT(//,10X,'P R O B L E M   N A M E',10X,A80,//)
 2000 FORMAT (///,
     1    ' C O N T R O L  I N F O R M A T I O N',//,
     2    '   NUMBER OF NODAL POINTS',10(' .'),' (NUMNP)=',I5,//,
     3    '   NUMBER OF ELEMENTS',12(' .'),'(NUMEL)=',I5,//,
     4    '   NUMBER OF MATERIALS       ',9(' .'),'(NUMAT)=',I5,//,
     5    '   SIZE OF THE TIME WINDOW',9(' .'),'(TM)=',F10.5,//,
     6    '   NUMBER OF LOAD INCREMENTS',9(' .'),'(NINCR)=',I5,//,
     7    '   DEGREE OF FREEDOM DIMENSION',6(' .'),'(NDOFDIM)=',I1,//,
     8    '   MAX.NUMBER OF NODES IN AN ELE.',6(' .'),'(NMNE)=',I2,//,
     9    '   MAX.NUMBER OF DOF IN AN ELE.',6(' .'),'(NMDOFEL)=',I2,//,
     1    '   MAX.NUMBER OF MAT PROPERTIES.',6(' .'),'(NMPR)=',I2,//,
     2    '   MAX.NUMBER OF INT MATP ROPERTIES.',6(' .'),'(NIPR)=',I2,//,
     3    '   MAX.NUMBER OF POINT LOADS.',6(' .'),'(NPL)=',I2)

 3000 FORMAT(1X,'INCREMENT=',I5,5X,'ITERATION=',I5,3X,'ERROR E1=',
     1       3X,F12.7)
 3005 FORMAT(1X,'INCREMENT=',I5,5X,'ITERATION=',I5,3X,'ERROR E3=',
     1       3X,F12.7)
 3016 FORMAT(1X,'TIME STEP COMPLETED=',2X,F5.2,2X,'TOT.TIME COMPLETED
     1       =',2X,F5.2,/)
 3020 FORMAT(5X,I5,10X,I5)
 4000 FORMAT(2X,'INC',6X,'TOT-TIME',6X,'INC-TIME')
 4010 FORMAT(2X,I3,3X,I3,7x,F6.3,9X,F6.3)
 3050 FORMAT(/,'DAMIAN JOB=',A8,1X,'COMPLETED')

      END PROGRAM
