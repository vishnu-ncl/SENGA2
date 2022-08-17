      PROGRAM DECOMP

      IMPLICIT NONE

      INTEGER NX,NY,NZ,NXPROC,NYPROC,NZPROC,NSPEC
      INTEGER NXN,NYN,NZN,NXPROC2,NYPROC2,NZPROC2
C     ===========================================#
C     ORIGINAL DOMAIN CONFIGURATION
      PARAMETER(NX=512,NY=1,NZ=1,
     +          NXPROC=1,NYPROC=1,NZPROC=1,
     +          NSPEC=9) 

C     DESIGN DOMAIN CONFIGURATION
      PARAMETER(NXN=512,NYN=256,NZN=256,
     +          NXPROC2=8,NYPROC2=4,NZPROC2=4)

C     CHOOSE THE POINT TO CENTER THE FLAME
      INTEGER :: PT=0
C     OUTPUT SWITCHES
      INTEGER :: CHECK=1
      INTEGER :: DECOMP_OUT=1

C     ===========================================#

      INTEGER NX1,NY1,NZ1,IL,JL,KL,ISPEC
      PARAMETER(NX1=NX/NXPROC,NY1=NY/NYPROC,NZ1=NZ/NZPROC)

C OLD DOMAIN - INITIALISED AVAIABLE ARRAY
      DOUBLE PRECISION 
     + TRUN(NX,NY,NZ),DRUN(NX,NY,NZ),URUN(NX,NY,NZ),VRUN(NX,NY,NZ),
     + WRUN(NX,NY,NZ),CRUN(NX,NY,NZ,NSPEC)

C TEMP ARRAYS
C      DOUBLE PRECISION 
C     + TEMP1(NX,NY,NZ),TEMP2(NX,NY,NZ),TEMP3(NX,NY,NZ),
C     + TEMP4(NX,NY,NZ),TEMP5(NX,NY,NZ),TEMP6(NX,NY,NZ)
C NEW DOMAIN - INITIALISED AVAIABLE ARRAY
      DOUBLE PRECISION
     + TRUNL(NXN,NYN,NZN),DRUNL(NXN,NYN,NZN),URUNL(NXN,NYN,NZN),
     + VRUNL(NXN,NYN,NZN),WRUNL(NXN,NYN,NZN),CRUNL(NXN,NYN,NZN,NSPEC)

      INTEGER I,J,K,IP,JP,KP,ISTART,IFINISH,JSTART,JFINISH,KSTART,
     +        KFINISH,IPROC,IC,JC,KC,N

      DOUBLE PRECISION GRIDR

      CHARACTER*1 PNFLAG
      CHARACTER*2 FILEFLAG
      CHARACTER*4 IDDUMP,JUSTB3
      CHARACTER*4 PNXDAT,PPROC,JUSTBB
      CHARACTER*6 PNPROC,JUSTB6
      CHARACTER*7 JUSTAA
      CHARACTER*60 FNPROC,JUSTMAT,JUSTMAT2
      PARAMETER(PNXDAT='.dat',PPROC='DMPI',JUSTAA='JUSTLAM')

C BEGIN
C      =================

      DO K = 1,NZN
        DO J = 1,NYN
         DO I = 1,NXN
            DRUNL(I,J,K) = 0.0D0
            URUNL(I,J,K) = 0.0D0
            VRUNL(I,J,K) = 0.0D0
            WRUNL(I,J,K) = 0.0D0
            TRUNL(I,J,K) = 0.0D0
           DO N= 1,NSPEC
            CRUNL(I,J,K,N) = 0.0D0
           END DO
          ENDDO
         ENDDO
       ENDDO



! READ THE 1-D LAMINAR FLAME DATA
      OPEN(UNIT=90,FILE='lamflame_py.dat',STATUS='OLD')
      READ(90,*)
      READ(90,*)
      DO I=1,NX
        READ(90,*)DRUN(I,1,1)
      END DO
      READ(90,*)
      DO I=1,NX
        READ(90,*)TRUN(I,1,1)
      END DO
      READ(90,*)
      DO I=1,NX
        READ(90,*)URUN(I,1,1)
      END DO
      DO J=1,NSPEC
        READ(90,*)
        DO I=1,NX
          READ(90,*)CRUN(I,1,1,J)
        END DO
      END DO
      CLOSE(90)
!SET the V and W Components of the Velocity to zero
      DO I=1,NX
        VRUN(I,1,1)=0.D0
        WRUN(I,1,1)=0.D0
      ENDDO



!TEMP ARRAYS ONLY USE IF YOU WANT TO INVERT THE DIRECTION OF THE FLAME

!!      TEMP1=DRUN(NX:1:-1,:,:)
!!      TEMP2=URUN(NX:1:-1,:,:)       
!!      TEMP3=VRUN(NX:1:-1,:,:)
!!      TEMP4=WRUN(NX:1:-1,:,:)
!!      TEMP5=ERUN(NX:1:-1,:,:)
!!      TEMP6=CRUN(NX:1:-1,:,:)
!!
!!
!!      DO I=1,NX
!!
!!      DRUN(I,1,1)=TEMP1(I,1,1)  
!!      URUN(I,1,1)=TEMP2(I,1,1)       
!!      VRUN(I,1,1)=TEMP3(I,1,1)
!!      WRUN(I,1,1)=TEMP4(I,1,1)
!!      ERUN(I,1,1)=TEMP5(I,1,1)
!!      CRUN(I,1,1)=TEMP6(I,1,1)
!!
!!      ENDDO
        
! CREAT A 3D DATA
      Do N=1,nspec
       DO K = 1,NZN
        DO J = 1,NYN
          DO I = 1,NXN

C         DRUNL(I,J,K) = DRUN(I,1,1)
C         URUNL(I,J,K) = URUN(I,1,1)
C         VRUNL(I,J,K) = VRUN(I,1,1)
C         WRUNL(I,J,K) = WRUN(I,1,1)
C         ERUNL(I,J,K) = ERUN(I,1,1)
C         CRUNL(I,J,K) = CRUN(I,1,1)

C UMOD
!             IF (I.EQ.NXN) THEN
!              DRUNL(I,J,K) = DRUN(NX,1,1)
!              URUNL(I,J,K) = URUN(NX,1,1)
!              VRUNL(I,J,K) = VRUN(NX,1,1)
!              WRUNL(I,J,K) = WRUN(NX,1,1)
!              ERUNL(I,J,K) = ERUN(NX,1,1)
!              CRUNL(I,J,K) = CRUN(NX,1,1)
!             ELSEIF (I.GE.401) THEN
             IF (I.GE.PT+1.AND.I.LE.PT+NX) THEN
              DRUNL(I,J,K) = DRUN(I-PT,1,1)
              URUNL(I,J,K) = URUN(I-PT,1,1)
              VRUNL(I,J,K) = VRUN(I-PT,1,1)
              WRUNL(I,J,K) = WRUN(I-PT,1,1)
              TRUNL(I,J,K) = TRUN(I-PT,1,1)
            CRUNL(I,J,K,N) = CRUN(I-PT,1,1,N)
!            CRUNL(I,J,K,2) = CRUN(I-99,1,1,2)
             ELSEIF (I.LT.PT) THEN
              DRUNL(I,J,K) = DRUN(1,1,1)
              URUNL(I,J,K) = URUN(1,1,1)
              VRUNL(I,J,K) = VRUN(1,1,1)
              WRUNL(I,J,K) = WRUN(1,1,1)
              TRUNL(I,J,K) = TRUN(1,1,1)
            CRUNL(I,J,K,N) = CRUN(1,1,1,N)
!            CRUNL(I,J,K,2) = CRUN(1,1,1,2)
             ELSEIF(I.GE.PT+NX) THEN
!!             IF(I.GT.NX) THEN
              !PRINT*,'TEST'
              DRUNL(I,J,K) = DRUN(NX,1,1)
              URUNL(I,J,K) = URUN(NX,1,1)
              VRUNL(I,J,K) = VRUN(NX,1,1)
              WRUNL(I,J,K) = WRUN(NX,1,1)
              TRUNL(I,J,K) = TRUN(NX,1,1)
            CRUNL(I,J,K,N) = CRUN(NX,1,1,N)
!            CRUNL(I,J,K,2) = CRUN(NX,1,1,2)
             ENDIF
           ENDDO
         ENDDO
        ENDDO
       ENDDO
         

      IF(DECOMP_OUT.EQ.1)THEN
  
      DO IP=1,NXPROC2
       DO JP=1,NYPROC2
        DO KP=1,NZPROC2

         IPROC=(IP-1)+(JP-1)*NXPROC2+(KP-1)*NXPROC2*NYPROC2

         ISTART=1+(IP-1)*(NXN/NXPROC2)
         IFINISH=(IP)*(NXN/NXPROC2)

         JSTART=1+(JP-1)*(NYN/NYPROC2)
         JFINISH=(JP)*(NYN/NYPROC2)

         KSTART=1+(KP-1)*(NZN/NZPROC2)
         KFINISH=(KP)*(NZN/NZPROC2)

         WRITE(JUSTB3,'(I4.4)') IPROC

         JUSTMAT = 'lamsol'//JUSTB3//'.dat'

         OPEN(UNIT=IPROC,FILE=JUSTMAT,STATUS='UNKNOWN')

         DO KC = KSTART,KFINISH
          DO JC = JSTART,JFINISH
           DO IC = ISTART,IFINISH

            WRITE(IPROC,*)DRUNL(IC,JC,KC),
     +                 URUNL(IC,JC,KC),
     +                 VRUNL(IC,JC,KC),
     +                 WRUNL(IC,JC,KC),
     +                 TRUNL(IC,JC,KC),
C     +                 (CRUN(IC,JC,KC,N),N=1,NSPEC)
     +                 CRUNL(IC,JC,KC,1),
     +                 CRUNL(IC,JC,KC,2),
     +                 CRUNL(IC,JC,KC,3),
     +                 CRUNL(IC,JC,KC,4),
     +                 CRUNL(IC,JC,KC,5),
     +                 CRUNL(IC,JC,KC,6),
     +                 CRUNL(IC,JC,KC,7),
     +                 CRUNL(IC,JC,KC,8),
     +                 CRUNL(IC,JC,KC,9)

            END DO
           END DO
         END DO

       CLOSE(IPROC)

       END DO
       END DO
       END DO
      ENDIF

C----------------------------TEC360 output---------------------

      IF(CHECK.EQ.1)THEN       
       OPEN(UNIT=10, FILE='3d.dat',STATUS='unknown')

       WRITE (10,*) 'VARIABLES="X","Y","Z","D","U","V","W","H","H2",',
     +                       '"H20","H2O2","HO2","N2","O","O2","OH","T"'
       
        WRITE (10,*) 'ZONE I=',NXN,'J=',NYN ,'K=',1,'F=POINT'
       
         k=NZN/2!DO K=1,NZN
          DO J=1,NYN
           DO I=1,NXN
       
             WRITE(10,'(24E17.9)')DBLE(I),DBLE(J),DBLE(K),DRUNL(I,J,K),
     +                           URUNL(I,J,K),VRUNL(I,J,K),WRUNL(I,J,K),
     +                           CRUNL(I,J,K,1),CRUNL(I,J,K,2),
     +                           CRUNL(I,J,K,3),CRUNL(I,J,K,4),
     +                           CRUNL(I,J,K,5),CRUNL(I,J,K,6),
     +                           CRUNL(I,J,K,7),CRUNL(I,J,K,8),
     +                           CRUNL(I,J,K,9),TRUNL(I,J,K)
                         
          ENDDO
         ENDDO
       ! ENDDO
       CLOSE(10)
      ENDIF
       STOP
       END

