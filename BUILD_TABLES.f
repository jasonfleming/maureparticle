C NATHAN DILL
C JUNE 12, 2007
C______________________________________________________________________
C======================================================================
      PROGRAM BUILD_TABLES
C----------------------------------------------------------------------
C THIS PROGRAM BUILDS THE NODE2EL AND EL2EL TABLES THAT ARE REQUIRED
C FOR THE MAUREPARTICLE PARTICLE TRACKING CODE
C INPUT IS AN ADCIRC GRID (FORT.14) 
C OUTPUT IS THE TWO FILES el2el.tbl AND node2el.tbl
C THIS CODE IS BASED ON CONNECT.F WRITTEN BY DANA K. SAVIDGE 
C AND BRIAN O. BLANTON, 1990
C----------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER I,J,K,L,CNT,SCORE,NE,NN
      INTEGER, ALLOCATABLE :: HITCNT(:),NOC(:,:),EL2EL(:,:),
     &                        NODE2EL(:,:)
      CHARACTER*80 FILE14
C----------------------------------------------------------------------
C READ ELEMENT CONNECTIVITY FROM FORT.14 FILE
      WRITE(*,*)'WHAT IS THE NAME OF THE FORT.14 FILE?'
      READ(*,*)FILE14
      OPEN(14,FILE=FILE14)
      
      READ(14,*)
      READ(14,*)NE,NN
      ALLOCATE ( HITCNT(NE),NOC(3,NE),EL2EL(3,NE),NODE2EL(12,NN) )
      
      DO I=1,NN
         READ(14,*)
      END DO
      DO I=1,NE
         READ(14,*)J,K,(NOC(L,I),L=1,3)
      END DO
      CLOSE(14)
          
C----------------------------------------------------------------------
C MAKE THE ELEMENT TO ELEMENT NEIGHBOR TABLE BY FINDING ELEMENTS WHICH 
C SHARE TWO NODES, THEN UPDATING ELE2EL TABLE FOR BOTH ELEMENTS
C UPDATING BOTH ELEMENT ENTRIES EACH TIME YOU FIND A SHARED EDGE AND 
C KEEPING TRACK OF THE HIT COUNT MAKES GENERATING THESE TABLES FASTER, 
C BUT REQUIRES MORE MEMORY THAN BLANTON AND SAVIDGE'S CONNECT2D
C . . INITIALIZE EL2EL AND HITCNT WITH ZEROS
      DO I=1,NE
         HITCNT(I)=0
         DO J=1,3
            EL2EL(J,I)=0
         END DO
      END DO
      
      WRITE(*,*)
      WRITE(*,*)'BUILDING ELEMENT TO ELEMENT NEIGHBOR TABLE'
      DO I=1,NE
cc       WRITE(*,*)I,' of ',NE
       IF (HITCNT(I).LT.3) THEN

         
         DO J=I,NE

          IF ((HITCNT(J).LT.3).AND.(I.NE.J)) THEN
            SCORE=0
            DO K=1,3
               DO L=1,3
                  IF (NOC(K,I).EQ.NOC(L,J)) THEN
                    SCORE=SCORE+1
                    GOTO 10
                  END IF
               END DO
 10          CONTINUE
            END DO
            
            IF (SCORE.EQ.2) THEN
C . . . . . THIS IS A HIT, INCREMENT THE HIT COUNT AND UPDATE EL2EL FOR BOTH ELEMENTS               
               HITCNT(I)=HITCNT(I)+1
               HITCNT(J)=HITCNT(J)+1
               EL2EL(HITCNT(I),I)=J
               EL2EL(HITCNT(J),J)=I
CC               GOTO 20
            END IF
          END IF
CC 20       CONTINUE
         END DO   
       END IF
      END DO    
             
C . . NOW WRITE THE RESULTS
      WRITE(*,*)'WRITING EL2EL.TBL'
      OPEN(11,FILE='EL2EL.TBL')
      DO I=1,NE
         WRITE(11,*) (EL2EL(J,I),J=1,3)
      END DO
      CLOSE(11)
      
C----------------------------------------------------------------------
C NOW MAKE THE NODE TO ELEMENT TABLE
C THIS WAS TAKEN ALMOST EXACTLY FROM CONNECT2D
c      WRITE(*,*)'BUILDING AND WRITING NODE TO ELEMENT TABLE'
c      OPEN(12,FILE='node2el.tbl')
c      DO I=1,12
c         NODE2EL(I)=0
c      END DO
       
c      CNT=0
c      DO I=1,NN
c        DO J=1,CNT
c           NODE2EL(J)=0
c        END DO
c        CNT=0
c        DO J=1,NE
c          DO K=1,3
c            IF (NOC(K,J).EQ.I) THEN
c              CNT=CNT+1
c              NODE2EL(CNT)=J
c              GOTO 30
c            END IF
c          END DO
c 30       CONTINUE
c        END DO
c        WRITE(12,*) (NODE2EL(L),L=1,12)
c      END DO
c      CLOSE(12)

C---------------------------------------------------------------------
C BUT THIS WAY IS FASTER
      WRITE(*,*)'BUILDING NODE2EL TABLE'
C . . INITIALIZE NODE2EL AND HITCNT WITH ZEROS
      DO I=1,NN
         HITCNT(I)=0
         DO J=1,12
            NODE2EL(J,I)=0
         END DO
      END DO

C . . LOOP THROUGH THE NOC TABLE AND UPDATE NODE2EL ACCORDINGLY             
      DO I=1,NE
cc         WRITE(*,*)I      
         DO J=1,3
            HITCNT(NOC(J,I))=HITCNT(NOC(J,I))+1
            NODE2EL(HITCNT(NOC(J,I)),NOC(J,I))=I
         END DO
      END DO

      WRITE(*,*)'WRITING NODE2EL.TBL'
      OPEN(12,FILE='NODE2EL.TBL')
      DO I=1,NN
         WRITE(12,*) (NODE2EL(J,I),J=1,12)
      END DO
      CLOSE(12)


      STOP            
      END PROGRAM
C======================================================================

C||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||      
      
C______________________________________________________________________
C======================================================================
      SUBROUTINE REMOVE_WEIRS()
C----------------------------------------------------------------------
C THIS SUBROUTINE WILL CREATE A NEW GRID FILE WITH WEIRS REMOVED FROM 
C THE ADCIRC MESH. THIS NEW MESH MAY BE NECESSARY IF PARTICLES ARE TO
C BE TRACKED OVER WEIRS, IT WILL TRY NOT TO MAKE THIN ELEMENTS
C----------------------------------------------------------------------      
      IMPLICIT NONE

      RETURN
      END SUBROUTINE
C======================================================================
