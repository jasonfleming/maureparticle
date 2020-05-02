C______________________________________________________________________
C======================================================================

C----------------------------------------------------------------------
C
C nld 2018 Aug 24 - lose_wetdry version:
C     modifications to lose particles when they enter a dry element  
C     lost particles are indicated in output by LOCAT = 0
C     modfication is made to reduce the unrealistic accumulation of 
C     particles near wet-dry boundaries that is caused by unrealistic 
C     low current velocities in dry elements, due to interpolation 
C     of zero current from dry nodes.
C
C A particle tracking program that works with 2DDI ADCIRC output
C and potentially any other 2D unstructured grid hydrodynamic model 
C
C This this code is based on the original SCILAB code described in:
C 
C URS Corp. et al., 2007. "Missippi River Reintroduction into Maurepas
C   Swamp Project PO-29,  Volume VII ov VII Diversion Modeling", Report
C   prepared for the Louisiana Department of Natural Resources and 
C   U.S. Environmental Protection Agency. 
C   http://lacoast.gov/reports/project/
C   Vol_VII_Diversion%20Modeling%20Report-Dec%208-FINAL.pdf
C 
C and the enhanced Scilab code described in:
C
C Dill, Nathan, 2007. "Hydrodyamic Modeling of a Hypothetical River
C   Diversion Near Empire, Louisiana",  Thesis submitted in partial
C   fulfillment of the Master of Science Degree in Civil Engineering,
C   Louisiana State University.
C      http://etd.lsu.edu/docs/available/etd-06142007-084318/
C
C This Fortran version of Maureparticle differs somewhat from the
C original SCILAB code. Important differences are: 
C 
C - there is no longer an option for using RMA2 output;
C - the treatment of particles hitting the boundary has been 
C   simplified;
C - in this version you must input an eddy diffusivity to simulate
C   random walk diffusion (use zero if you don't want any diffusion);
C - you must specify the release time for each particle where the 
C   timing coincides with the times in the ADCIRC fort.64 file 
C   (this provides flexibility for setting up continuous releases and
C    such); 
C - to better approximate surface drifters, an option to add a fraction
C   of the wind velocity to the drift velocity was added;
C - and it is much faster than SCILAB.
C
C There is also a slightly less featured, slightly buggy, embarssingly
C parallel MPI version of this program.  Please contact the author
C if you are intersted. 
C
C it should compile pretty easily with any Fortran90 capable compiler
C e.g.
C       gfortran -o maurpt.exe Maureparticle.f
C---------------------------------------------------------------------
C Input:
C
C FORT.14 - ADCIRC grid file 
C
C FORT.64 - ADCIRC global time series of depth averaged current
C           if you want good results, save your FORT.64 output 
C           frequently.
C
C FORT.74 - (optional) ADCIRC global time series of wind velocity
C           if NWS.NE.0 add WFACTOR*(wind velocity) to the drift 
C           velocity. must be saved at same frequency as FORT.64
C
C PARTICLES.INP - run control parameters, initial particle 
C                 positions and release times. see example.
C
C NODE2EL.TBL - describes relation between nodes and their connected 
C               elements. use BUILD_TABLES.f to create this file
C
C EL2EL.TBL - an element to element neighbor table created by
C               BUILD_TABLES.f
C
C
C
C--------------------------------------------------------------------- 
C Copyright (C) 2007, 2008, 2013-2016 Nathan Dill
C
C This program  is free software; you can redistribute it and/or
C modify it under the terms of the GNU General Public License as
C published by the Free Software Foundation; either version 3 of the 
C License, or (at your option) any later version. 
C
C This program is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this program; if not, write to the Free Software 
C Foundation, Inc., 59 Temple Place - Suite 330,Boston, MA  02111-1307,
C USA.\
C
C you can contact the author at natedill(AT)gmail.com 
C
C----------------------------------------------------------------------
      MODULE MAUREPARAMS
      !------------------------- VARIABLES ----------------------------      
      INTEGER NN,NE,NP,NVTS,RK2,FOUND,NSTEPS,DYN,VELCNT,NWS,ICS
      INTEGER, ALLOCATABLE :: NOC(:,:),NOD2EL(:,:),EL2EL(:,:),LOCAT(:)
      INTEGER, ALLOCATABLE :: LOST(:),PID(:),ITRKING(:)
      
      REAL*8 TS,VTIMINC,RNTIM,OUTPER,TIME,OUTTIME,STDY_TIME,FRACTIME
      REAL*8 RELEASEPER,VELTIME1,VELTIME2,EDDY_DIF,WFACTOR,SLAM0,SFEA0
      REAL*8, ALLOCATABLE :: X(:),Y(:),VX(:),VY(:),RLSTIME(:)
cc      ,VXP(:),VYP(:)
      REAL*8, ALLOCATABLE :: XP(:),YP(:),VX2(:),VY2(:)
      
      CHARACTER*80 DESC1,DESC2

      logical :: metonly  ! .true. if particles should track wind instead of water current
      character(len=1024) :: maureParameterInputFile ! particle tracking control plus initial particle location
      character(len=1024) :: meshFile            ! adcirc fort.14
      character(len=1024) :: maureParticleOutputFile ! locations v time
      character(len=1024) :: velocityFile ! locations v time
      character(len=1024) :: windVelocityFile ! locations v time     
      
      END MODULE MAUREPARAMS
      !----------------------------------------------------------------

C----------------------------------------------------------------------
      PROGRAM MAUREPARTICLE
      USE MAUREPARAMS
      IMPLICIT NONE

      integer :: argcount ! number of command line arguments
      character(len=1024) :: cmdlineopt ! command line option
      character(len=1024) :: cmdlinearg ! command line option
      integer :: errorio  ! to detect i/o errors
      !      
      ! Initialize default values
      maureParameterInputFile = 'PARTICLES.INP'
      meshFile = 'FORT.14'
      velocityFile = 'FORT.64'
      windVelocityFile = 'FORT.74'
      maureParticleOutputFile = 'MAUREPT.OUT'
      metonly = .false.
      !
      ! Get command line options
      argcount = command_argument_count() ! count up command line options
      if (argcount.gt.0) then
         i=0
         do while (i.lt.argcount)
            i = i + 1
            call getarg(i, cmdlineopt)
            downcaseCmdlineopt = call downcase(cmdlineopt)
            select case(trim(downcaseCmdlineopt))
               case("--metonly")
                  metonly = .true.
                  call echoCmdLineOpt(cmdlineopt,cmdlinearg)
               case("--meshfile")
                  i = i + 1
                  call getarg(i, cmdlinearg)
                  call echoCmdLineOpt(cmdlineopt,cmdlinearg)
                  meshFile = trim(cmdlinearg)
               case("--velocityfile")
                  i = i + 1
                  call getarg(i, cmdlinearg)
                  call echoCmdLineOpt(cmdlineopt,cmdlinearg)
                  velocityFile = trim(cmdlinearg)
               case("--windvelocityfile")
                  i = i + 1
                  call getarg(i, cmdlinearg)
                  call echoCmdLineOpt(cmdlineopt,cmdlinearg)
                  windVelocityFile = trim(cmdlinearg)
               case("--maureparticleoutputfile")
                  i = i + 1
                  call getarg(i, cmdlinearg)
                 call echoCmdLineOpt(cmdlineopt,cmdlinearg)
                  maureParticleOutputFile = trim(cmdlinearg)
               case("--maureparameterinputfile")
                  i = i + 1
                  call getarg(i, cmdlinearg)
                  call echoCmdLineOpt(cmdlineopt,cmdlinearg)
                  maureParticleInputFile = trim(cmdlinearg)                  
               case default
                  write(6,'(a,a,a)') "WARNING: Command line option '",
     &             TRIM(cmdlineopt),"' was not recognized."
            end select
         end do
      end if      

C . . GET INFO TO ALLOCATE ARRAYS
      inUnit = availableUnitNumber()
      call openFileForRead(inUnit,maureParameterInputFile,errorio)
      mUnit = availableUnitNumber()
      call openFileForRead(inUnit,meshFile,errorio)
      oUnit = availableUnitNumber()
      call openFileForRead(inUnit,maureParticleOutputFile,errorio)
      OPEN(UNIT=11,FILE='PARTICLES.INP')
      OPEN(UNIT=14,FILE='FORT.14')
      OPEN(UNIT=15,FILE='MAUREPT.OUT')
      
      READ(mUnit,'(A80)')DESC1
      READ(mUnit,*)NE,NN
      READ(inUnit,'(A80)')DESC2
      READ(inUnit,*)NP
      READ(inUnit,*)TS
      READ(inUnit,*)RNTIM
      READ(inUnit,*)OUTPER
      READ(inUnit,*)RK2
      READ(inUnit,*)DYN
      READ(inUnit,*)STDY_TIME
      READ(inUnit,*)EDDY_DIF
      READ(inUnit,*)NWS
      READ(inUnit,*)WFACTOR
      READ(inUnit,*)ICS
      READ(inUnit,*)SLAM0,SFEA0
      WRITE(*,*)'NE,NN,NP',NE,NN,NP
      CLOSE(mUnit)                      
      CLOSE(inUnit)
      
      ALLOCATE( NOC(3,NE),NOD2EL(12,NN),EL2EL(3,NE),LOCAT(NP) )
      ALLOCATE( X(NN),Y(NN),VX(NN),VY(NN),VX2(NN),VY2(NN),
     &      ITRKING(NP),RLSTIME(NP),PID(NP),LOST(NP),XP(NP),YP(NP) )    
     
      DO I=1,NP
        ITRKING(I)=0
      END DO
c----------------------------------------------------------------------
     
C . . NOW READ ALL THE INPUT DATA . . . . . . . . . . . . . . . . . . .
      WRITE(*,*)'READING GRID DATA FROM: ',TRIM(DESC1)
      WRITE(*,*)
      WRITE(*,*)'NN = ',NN,' NE = ',NE
      WRITE(*,*)
      WRITE(*,*)'READING PARTICLE DATA FROM: ',TRIM(DESC2)
      WRITE(*,*)
      WRITE(*,*)'NP = ',NP,' TS = ',TS
      WRITE(*,*)
      IF (RK2.EQ.1) THEN
       WRITE(*,*)'VELOCITY INTEGRATION BY 2ND ORDER RUNGE-KUTTA METHOD' 
         WRITE(*,*)
      ELSE
         WRITE(*,*)"VELOCITY INTEGRATION BY EULER'S METHOD"
         WRITE(*,*)
      END IF
      !
      if (metonly.eqv..true.) then
         WRITE(*,*) 
     &    'Only wind velocity will be used for particle tracking.'
      endif
      ! 
      ! read mesh, mapping tables, initial particle positions,
      ! and first datasets from water current velocity (fort.64)
      ! and/or wind velocity (fort.74)
      CALL READ_DATA()
      WRITE(*,*)'ALL INITIAL INPUT DATA HAS BEEN READ SUCCESSFULLY'
      !
      ! time varying velocity
      IF (DYN.EQ.1) THEN
         WRITE(*,*)'VELOCITY SOLUTION WILL BE READ DYNAMICALLY'
         WRITE(*,*)'STARTING AT ',STDY_TIME
         CALL READ_64_DYN()
      !
      ! steady velocity
      ELSE
         WRITE(*,*)'STEADY-STATE SIMULATION, READING FORT.64 UNTIL ',
     &              STDY_TIME    
         CALL READ_64_STDY()
         CLOSE(64)
         WRITE(*,*)'VELOCITY SOLUTION AT ',STDY_TIME,
     &             ' WILL BE TAKEN AS STEADY-STATE SOLUTION'
      END IF
C----------------------------------------------------------------------
      
C . . FIGURE NUMBER OF TIMESTEPS. . . . . . . . . . . . . . . . . . . .
      NSTEPS=INT(RNTIM/TS)
      
C . . FIGURE TOTAL NUMBER OF PARTICLES FOR CONTINUOUS RELEASE . . . . . 
C                 AND ALLOCATE ARRAYS ACCORDINGLY                       
      

C-------------------INITIAL SEARCH FOR PARTICLES-----------------------      
C . . INITIALIZE LOST VECTOR. . . . . . . . . . . . . . . . . . . . . .   
      DO I=1,NP
         PID(I)=I
         LOST(I)=0
      END DO

C . . SEARCH ALL ELEMENTS FOR INITIAL PARTICLE LOCATIONS. . . . . . . .
C . . THIS WILL BE SKIPPED IF LOCAT IS SPECIFIED IN PARTICLE.INP. . . .
      WRITE(*,*)'SEARCHING FOR PARTICLES . . .'
 20   CONTINUE     
      DO I=1,NP
         IF (LOCAT(I).EQ.0) THEN
c            WRITE(*,*)'SEARCHING FOR PARTICLE, ',
            DO J=1,NE
               CALL LOCAT_CHK(J,XP(I),YP(I),NOC,X,Y,FOUND) 
               IF (FOUND.EQ.1) THEN
                   LOCAT(I)=J
                   GOTO 30
               END IF
            END DO
            IF (LOCAT(I).EQ.0) THEN
              WRITE(*,*)'!!! WARNING - PARTICLE ',I,'CANT BE FOUND !!!!'
               WRITE(*,*)'!!!  IT WILL BE LOST FROM THE BEGINNING  !!!!'
               LOST(I)=1
c               STOP
            END IF
         ELSE
            CALL LOCAT_CHK(LOCAT(I),XP(I),YP(I),NOC,X,Y,FOUND)
            IF (FOUND.EQ.0) THEN
               WRITE(*,*)'!! WARNING - BAD LOCAT GIVEN FOR PARTICLE ',I
               WRITE(*,*)'!! WILL NOW SEARCH FOR PROPER LOCAT'
               LOCAT(I)=0
               GOTO 20
            END IF 
         END IF
 30      CONTINUE            
      END DO      
C----------------------------------------------------------------------      
      
C----------------------THIS IS THE TRACKING LOOP-----------------------
C . . INITIALIZE SOME TIME KEEPING VARIABLES. . . . . . . . . . . . . .
      TIME=STDY_TIME
      OUTTIME=STDY_TIME
      VELTIME1=STDY_TIME
      VELTIME2=VELTIME1+VTIMINC
      FRACTIME=0.D0
      
C . . TRAKCING LOOP . . . . . . . . . . . . . . . . . . . . . . . . . .      


      DO J=1,NSTEPS
c         WRITE(*,*)'TRACKING STEP ',J,' OF ',NSTEPS
C . . . .WRITE OUTPUT ? . . . . . . . . . . . . . . . . . . . . . . . .                           
         IF (TIME.EQ.OUTTIME) THEN
           CALL WRITE_DATA(NP,XP,YP,LOCAT,PID,TIME,ICS,SLAM0,SFEA0,LOST)
           OUTTIME=OUTTIME+OUTPER
         END IF

C . . . .RELEASE MORE PARTICLES ?. . . . . . . . . . . . . . . . . . . .        
         DO I=1,NP      
            IF (TIME.GE.RLSTIME(I)) ITRKING(I)=1;
         END DO
         
C . . . .ALL PARTICLES TAKE A STEP. . . . . . . . . . . . . . . . . . .                 
         IF (RK2.EQ.1) THEN
            CALL RK2_STEP(EL2EL,NOD2EL,NP,TS,NOC,X,Y,VX,VY,LOCAT,XP,YP,
     &                     LOST,VX2,VY2,DYN,FRACTIME,EDDY_DIF,ITRKING)
         ELSE      
            CALL EULER_STEP(DYN,NP,TS,NOC,X,Y,VX,VY,VX2,VY2,LOCAT,
     &                     XP,YP,LOST,FRACTIME,EDDY_DIF,ITRKING)
         END IF   
         
C . . . . . UPDATE LOCAT IF NECESSARY . . . . . . . . . . . . .         
         DO I=1,NP
         
            CALL  LOCAT_CHK(LOCAT(I),XP(I),YP(I),NOC,X,Y,FOUND)
            
            IF (FOUND.EQ.0) THEN
               CALL UPDATE_LOCAT(EL2EL,NOD2EL,X,Y,NOC,XP(I),YP(I),
     &                                             LOCAT(I),LOST(I))
            
               IF (LOST(I).EQ.1) THEN
                  XP(I)=-99999
                  YP(I)=-99999
                  WRITE(*,40)' PARTICLE, ',I,' LOST AT TIME, ',TIME,
     &                          ' FROM ELEMENT ',LOCAT(I)          
               END IF      
            END IF
         END DO
 40      Format (A,I9,A,F13.3,A,I9)
C . . . .UPDATE SIMULATION TIME . . . . . . . . . . . . . . . . . . . .            
         TIME=TIME+TS
         
C . . . .CHECK TO SEE IF WE NEED TO GET NEW VELOCITY DATA         
         IF (DYN.EQ.1) THEN
            IF (TIME.GT.VELTIME2) THEN
               CALL READ_64_DYN(STDY_TIME,NN,VX,VY,VX2,VY2,0,
     &                                            NWS,WFACTOR)
               VELTIME1=VELTIME2
               VELTIME2=VELTIME1+VTIMINC
            END IF
            FRACTIME=(TIME-VELTIME1)/VTIMINC
         END IF
         
      END DO
C------------------------END OF TRACKING LOOP--------------------------         

      CLOSE(15)

      STOP
      END PROGRAM
C======================================================================

C||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

C______________________________________________________________________
C======================================================================
      SUBROUTINE EULER_STEP(DYN,NP,TS,NOC,X,Y,VX,VY,VX2,VY2,LOCAT,XP,
     &                      YP,LOST,FRACTIME,EDDY_DIF,ITRKING)
C----------------------------------------------------------------------      
C PARTICLE POSITIONS (XP,YP) ARE INPUT AND OUTPUT, PARTICLE POSITIONS
C ARE MOVED BY TAKING A SIMPLE STEP USING EULER'S METHOD 
C I.E. VELOCITY*TIME=DISPLACEMENT
C INPUT INCLUDES NUMBER OF PARTICLES NP, TIMESTEP TS, NODE CONNECTIVITY
C TABLE [NOC], NODAL POSITIONS  {X} AND {Y} AND VELOCITIES {V} , AND 
C PARTICLE ELEMENTAL LOCATIONS {LOCAT}
C FOR DYNAMIC RUNS(DYN=1) {VX} AND {VY} ARE THE VELOCITIES AT THE 
C BEGINNING OF TIME INTERVAL AND {VX2} AND {VY2} ARE AT THE ENC OF THE
C INTERVAL AND THE VELOCITY IS INTERPOLATED LINEARLY IN BETWEEN USING
C FRACTIME AS THE FRACTION OF THE TIME INTERVAL FROM T1 TO T2.
C----------------------------------------------------------------------
      USE MAUREPARAMS
      IMPLICIT NONE
      
      INTEGER NP,NOC(3,1),LOCAT(1),I,J,K,LOST(1),DYN,ITRKING(1)
      DOUBLE PRECISION TS,  X(1),Y(1),VX(1),VY(1),VX2(1),VY2(1),
     &       EDDY_DIF,R,ANG,XP(1),YP(1),VXP,VYP,VVX(3),VVY(3),FRACTIME
      LOGICAL L_ISWET
           
      DO J=1,NP
       IF (ITRKING(J).EQ.1) THEN
      
       IF (LOST(J).EQ.0) THEN
          I=LOCAT(J)            
          IF (DYN.EQ.1) THEN
C . . . . . .INTERPOLATE VELOCITIES IN TIME FOR DYNAMIC SIMULATION 
             DO K=1,3
                VVX(K)=VX(NOC(K,I))+(VX2(NOC(K,I))-VX(NOC(K,I)))
     &                                                     *FRACTIME
                VVY(K)=VY(NOC(K,I))+(VY2(NOC(K,I))-VY(NOC(K,I)))
     &                                                     *FRACTIME
             END DO
          ELSE   
C . . . . . .USE THE STEADY STATE VELOCITIES                    
             DO K=1,3
                VVX(K)=VX(NOC(K,I))
                VVY(K)=VY(NOC(K,I))
             END DO
          
          END IF
            

C . . . .GET THE X-VELOCITY AT THE PARTICLE POSITION
         CALL VELINTRP(X(NOC(1,I)),Y(NOC(1,I)),VVX(1),
     &                 X(NOC(2,I)),Y(NOC(2,I)),VVX(2),
     &                 X(NOC(3,I)),Y(NOC(3,I)),VVX(3),     
     &                         XP(J),YP(J),VXP)
 
C . . . .GET THE Y-VELOCITY AT THE PARTICLE POSITION      
         CALL VELINTRP(X(NOC(1,I)),Y(NOC(1,I)),VVY(1),
     &                 X(NOC(2,I)),Y(NOC(2,I)),VVY(2),
     &                 X(NOC(3,I)),Y(NOC(3,I)),VVY(3),     
     &                         XP(J),YP(J),VYP)   

C . .    CALCULATE NEW POSITIONS USING VELOCITY AT MIDPOINT AND ORIGINAL POSITIO
C . .    DON'T DIFFUSE IF IN A DRY ELEMENT
C . .    DO STILL ALLOW MOVEMENT IN A DRY ELEMENT WITHOUT DIFFUSION DEPENDING
C . .    ON THE VELOCOITY OF ANY WET NODES IN THE ELEMENT, WHICH SHOULD
C . .    HELP PARTICLES FROM GETTING STUCK ON A WET/DRY BOUNDARY
C . .    ASSUME WE'RE IN A DRY ELEMENT IF ANY OF THE NODES HAVE VELOCITY LESS 
C . .    THAN THE MACHINE PRECISION (FROM EPSILON FUNCTION)
         L_ISWET=.TRUE.
         DO K=1,3
            IF ((ABS(VVX(K)).LE.EPSILON(VVX(K))).AND.
     &          (ABS(VVY(K)).LE.EPSILON(VVY(K)))) THEN
               L_ISWET=.FALSE.
            END IF
         END DO
         IF (.NOT.L_ISWET) LOST(J)=1

         IF (L_ISWET.AND.(EDDY_DIF .GT. 0.D0)) THEN 
           CALL RANDOM_NUMBER(R)
           CALL RANDOM_NUMBER(ANG)
           ANG=6.28318530717959D0*ANG
           R=R*(EDDY_DIF*TS)**0.5D0
           XP(J)=XP(J)+VXP*TS + R * COS(ANG)
           YP(J)=YP(J)+VYP*TS + R * SIN(ANG)
         ELSE
           XP(J)=XP(J)+VXP*TS 
           YP(J)=YP(J)+VYP*TS 
         END IF
       END IF !LOST
       END IF !TRACKING
      END DO   
      
      RETURN
      END SUBROUTINE    
C======================================================================

C||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

C______________________________________________________________________
C======================================================================
      SUBROUTINE RK2_STEP(EL2EL,NOD2EL,NP,TS,NOC,X,Y,VX,VY,LOCAT,XP,YP,
     &                     LOST,VX2,VY2,DYN,FRACTIME,EDDY_DIF,ITRKING)
C----------------------------------------------------------------------
C THIS SUBROUTINE USES THE 2ND ORDER RUNGE-KUTTA INTEGRATION METHOD 
C TO ADVANCE THE PARTICLE POSTIONS OVER ONE TIME STEP, IN ADDITION
C TO THE INPUT ARGUMENTS NEEDED FOR EULER_STEP() THIS ONE ALSO
C NEEDS THE EL2EL AND NOD2EL TABLES.
C----------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER NP,NOC(3,1),LOCAT(1),I,J,K,FOUND,EL2EL(3,1),NOD2EL(12,1)
      INTEGER CLOSEST,LOST(1),DYN,ITRKING(1)
      LOGICAL L_ISWET
      DOUBLE PRECISION TS,X(1),Y(1),VX(1),VY(1),XXP(NP),YYP(NP),VX2(1)
      DOUBLE PRECISION XP(1),YP(1),VXP,VYP,DS1,DS2,DS3,VY2(1)    
      DOUBLE PRECISION FRACTIME,VVX(3),VVY(3),EDDY_DIF,R,ANG
C . . SAVE THE STARTING POSITIONS 
      DO I=1,NP
         IF (LOST(I).EQ.0) THEN
            XXP(I)=XP(I)
            YYP(I)=YP(I)
         END IF
      END DO
      
CC      WRITE(*,*)'STARTING POSITION SAVED'
      
C . . DO AN EULER STEP AS A FIRST GUESS
C . . FOR THIS CALL DIFFUSIVITY IS SET TO ZERO
      CALL EULER_STEP(DYN,NP,TS,NOC,X,Y,VX,VY,VX2,VY2,LOCAT,
     &                     XP,YP,LOST,FRACTIME,0.D0,ITRKING)

CC      WRITE(*,*)'EULER GUESS COMPLETED'
     
C . . FIND THE MID POINT
      DO I=1,NP
         IF (LOST(I).EQ.0) THEN
            XP(I)=(XXP(I)+XP(I))/2
            YP(I)=(YYP(I)+YP(I))/2
         END IF
      END DO
      
CC      WRITE(*,*)'MIDPOINT CALCULATED'
C----------------------------------------------------------------------      
C . . CHECK LOCAT OF MIDPOINT
      DO I=1,NP 
      IF (ITRKING(I).EQ.1) THEN
       IF (LOST(I).EQ.0) THEN    
         CALL LOCAT_CHK(LOCAT(I),XP(I),YP(I),NOC,X,Y,FOUND)
         
cc         WRITE(*,*)'MIDPOINT CHECKED, PARTICLE ',I
         
         IF (FOUND.EQ.0) THEN
        
         CALL UPDATE_LOCAT(EL2EL,NOD2EL,X,Y,NOC,XP(I),YP(I),LOCAT(I),K)
         
cc         WRITE(*,*)'MIDPOINT UPDATED,  PARTICLE ',I
         
            IF (K.EQ.1) THEN
               XP(I)=XXP(I)
               YP(I)=YYP(I)
            END IF        
         END IF
       END IF
      END IF
      END DO 
C----------------------------------------------------------------------         
      
C . . GET THE VELOCITY AT THE MIDPOINT        

      
      DO J=1,NP
      IF (ITRKING(J).EQ.1) THEN
       IF (LOST(J).EQ.0) THEN   
         I=LOCAT(J)
          IF (DYN.EQ.1) THEN
C . . . . . .INTERPOLATE VELOCITIES IN TIME FOR DYNAMIC SIMULATION 
             DO K=1,3
                VVX(K)=VX(NOC(K,I))+(VX2(NOC(K,I))-VX(NOC(K,I)))
     &                                                     *FRACTIME
                VVY(K)=VY(NOC(K,I))+(VY2(NOC(K,I))-VY(NOC(K,I)))
     &                                                     *FRACTIME
             END DO
          ELSE   
C . . . . . .USE THE STEADY STATE VELOCITIES                    
             DO K=1,3
                VVX(K)=VX(NOC(K,I))
                VVY(K)=VY(NOC(K,I))
             END DO
          
          END IF    
cc         WRITE(*,*)'IM AT LINE 359'
C . . . .GET THE X-VELOCITY AT THE MIDPOINT
         CALL VELINTRP(X(NOC(1,I)),Y(NOC(1,I)),VVX(1),
     &                 X(NOC(2,I)),Y(NOC(2,I)),VVX(2),
     &                 X(NOC(3,I)),Y(NOC(3,I)),VVX(3),     
     &                         XP(J),YP(J),VXP)
 
C . . . .GET THE Y-VELOCITY AT THE MIDPOINT      
         CALL VELINTRP(X(NOC(1,I)),Y(NOC(1,I)),VVY(1),
     &                 X(NOC(2,I)),Y(NOC(2,I)),VVY(2),
     &                 X(NOC(3,I)),Y(NOC(3,I)),VVY(3),     
     &                         XP(J),YP(J),VYP)   

                
C . .    CALCULATE NEW POSITIONS USING VELOCITY AT MIDPOINT AND ORIGINAL POSITIO
C . .    DON'T DIFFUSE IF IN A DRY ELEMENT
C . .    DO STILL ALLOW MOVEMENT IN A DRY ELEMENT WITHOUT DIFFUSION DEPENDING
C . .    ON THE VELOCOITY OF ANY WET NODES IN THE ELEMENT, WHICH SHOULD
C . .    HELP PARTICLES FROM GETTING STUCK ON A WET/DRY BOUNDARY
C . .    ASSUME WE'RE IN A DRY ELEMENT IF ANY OF THE NODES HAVE VELOCITY LESS 
C . .    THAN THE MACHINE PRECISION (FROM EPSILON FUNCTION)
         L_ISWET=.TRUE.
         DO K=1,3
            IF ((ABS(VVX(K)).LE.EPSILON(VVX(K))).AND.
     &          (ABS(VVY(K)).LE.EPSILON(VVY(K)))) THEN
               L_ISWET=.FALSE.
            END IF
         END DO
         IF (.NOT.L_ISWET) LOST(J)=1

         IF (L_ISWET.AND.(EDDY_DIF .GT. 0.D0)) THEN 
           CALL RANDOM_NUMBER(R)
           CALL RANDOM_NUMBER(ANG)
           ANG=6.28318530717959D0*ANG
           R=R*(EDDY_DIF*TS)**0.5D0
           XP(J)=XP(J)+VXP*TS + R * COS(ANG)
           YP(J)=YP(J)+VYP*TS + R * SIN(ANG)
         ELSE
           XP(J)=XP(J)+VXP*TS 
           YP(J)=YP(J)+VYP*TS 
         END IF
       END IF
      END IF
      END DO            
      
      
      RETURN
      END SUBROUTINE
C======================================================================

C||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

C______________________________________________________________________
C======================================================================
      SUBROUTINE LOCAT_CHK(LOCAT,XP,YP,NOC,X,Y,FOUND)
C----------------------------------------------------------------------
C THIS SUBROUTINE CHECKS IF A PARTICLE RESIDES WITHIN THE ELEMENT LOCAT
C IT RETURNS THE VALUE FOUND=1 IF THE PARTICLE IS FOUND OR FOUND=0 IF 
C IT WAS NOT FOUND. LOCAT,XP,YP ARE SCALAR INPUT. 
C----------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER LOCAT,FOUND,NOC(3,1),I
      DOUBLE PRECISION XP,YP,X(1),Y(1),X1,X2,X3,Y1,Y2,Y3,
     &                 DS1(2),DS2(2),DS3(2),CROSS,C1,C2,C3
      
     
      FOUND=0
      IF (LOCAT.EQ.0) RETURN
             
C . . GET DISPLACEMENTS FROM PARTICLE TO NODES    
      DS1(1)=X(NOC(1,LOCAT))-XP
      DS1(2)=Y(NOC(1,LOCAT))-YP
      DS2(1)=X(NOC(2,LOCAT))-XP
      DS2(2)=Y(NOC(2,LOCAT))-YP
      DS3(1)=X(NOC(3,LOCAT))-XP
      DS3(2)=Y(NOC(3,LOCAT))-YP

C . . ALL + CROSS PRODS. MEANS PART. IS FOUND IF NOC IS ANTI-CLOCKWISE      
      C1=CROSS(DS1,DS2)
      C2=CROSS(DS2,DS3)
      C3=CROSS(DS3,DS1)
       
      IF ((C1.GE.0).AND.(C2.GE.0).AND.(C3.GE.0)) FOUND=1
      
      RETURN
      END SUBROUTINE
C======================================================================

C||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

C______________________________________________________________________
C======================================================================
       SUBROUTINE VELINTRP(X0,Y0,V0,X1,Y1,V1,X2,Y2,V2,XP,YP,VP)
C----------------------------------------------------------------------
C THIS SUBROUTINE RETURNS (VP), THE VALUE AT POINT (XP,YP) LINEARLY
C INTERPOLATED ON A PLANE BETWEEN THE THREE POINTS (X0,Y0)(X1,Y1)(X2,Y2)
C----------------------------------------------------------------------     
      IMPLICIT NONE
      
      DOUBLE PRECISION X0,Y0,V0,X1,Y1,V1,X2,Y2,V2,XP,YP,VP,T,U,DET,
     &    XX1,YY1,XX2,YY2,XXP,YYP 
      
      XX1=X1-X0
      XX2=X2-X0
      YY1=Y1-Y0
      YY2=Y2-Y0
      XXP=XP-X0
      YYP=YP-Y0
      
      DET=(YY2*XX1-XX2*YY1)

      T=(XXP*YY2-XX2*YYP)/DET
      U=(XX1*YYP-YY1*XXP)/DET

      VP=V0+T*(V1-V0)+U*(V2-V0)
      
      RETURN
      END SUBROUTINE
C======================================================================

C||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

C______________________________________________________________________
C======================================================================
      FUNCTION CROSS(DS1,DS2)
C----------------------------------------------------------------------
C THIS FUNCTION RETURNS THE CROSSPRODUCT OF TWO 2D VECTORS
C----------------------------------------------------------------------
      IMPLICIT NONE
      
      DOUBLE PRECISION CROSS,DS1(2),DS2(2)
      
      CROSS=(DS1(1)*DS2(2))-(DS1(2)*DS2(1))
      
      RETURN
      END FUNCTION
C======================================================================

C||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

C______________________________________________________________________
C======================================================================
      SUBROUTINE UPDATE_LOCAT(EL2EL,NOD2EL,X,Y,NOC,XP,YP,LOCAT,LST)
C----------------------------------------------------------------------
C THIS SUBROUTINE REQUIRES SCALAR INPUT OF PARTICLE POSITION AND LOCAT  
C IT SEARCHES THE EL2EL TABLE AND NOD2EL TABLE, 
C IT MAY RETURN A NEW VALUE LOCAT AND WILL RETURN LST=1 IF PARTICLE 
C LEAVES BOUNDARY OR SKIPS A NODE/ELEMENT GROUP
C----------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER I,J,K,LOCAT,EL2EL(3,1),NOD2EL(12,1),NOC(3,1),FOUND
      INTEGER BEL1,BEL2,CLOSEST,LST,ICNT,NODS1(3),NODS2(3)
      DOUBLE PRECISION X(1),Y(1),DS(3),DSMIN,XP,YP,X1,Y1,X2,Y2

C . . SKIP FOR LOST PARTICLES . . . . . . . . . . . . . . . . . . . . .
      IF (LST.GT.0) THEN
         LST=2
         GOTO 100
      END IF   
      
      LST=0
            
C . . SEARCH EL2EL. . . . . . . . . . . . . . . . . . . . . . . . . . .
      DO J=1,3
          K=EL2EL(J,LOCAT) 
          IF (K.EQ.0) GOTO 40     
          CALL LOCAT_CHK(K,XP,YP,NOC,X,Y,FOUND)
          IF (FOUND.EQ.1) THEN 
             LOCAT=K
             GOTO 100
          END IF   
      END DO
 40   CONTINUE     
C . . IF IT WAS NOT FOUND IN EL2EL SEARCH NOD2EL. . . . . . . . . . . .
C . . FIND THE CLOSEST NODE . . . . . . . . . . . . . . . . . . . . . .
      DS(1)=(X(NOC(1,LOCAT))-XP)**2 + (Y(NOC(1,LOCAT))-YP)**2
      DS(2)=(X(NOC(2,LOCAT))-XP)**2 + (Y(NOC(2,LOCAT))-YP)**2
      DS(3)=(X(NOC(3,LOCAT))-XP)**2 + (Y(NOC(3,LOCAT))-YP)**2
      
      DSMIN=MIN(DS(1),DS(2),DS(3))
      DO J=1,3
         IF (DSMIN.EQ.DS(J)) THEN 
            CLOSEST=NOC(J,LOCAT)
            GOTO 50
         END IF
      END DO     
 50   CONTINUE

      DO J=1,12
         K=NOD2EL(J,CLOSEST)
         IF (K.EQ.0) GOTO 60
         CALL LOCAT_CHK(K,XP,YP,NOC,X,Y,FOUND) 
         IF (FOUND.EQ.1) THEN
            LOCAT=K
            GOTO 100
         END IF
      END DO
 60   CONTINUE
      
C . . STILL NOT FOUND,. . . . . . . . . . . . . . . . . . . . . . . . .
C . . IT MUST HAVE LEFT BOUNDARY OR SKIPPED AN ELEMENT NODE GROUP . . . 
C . . FIND THE ELEMENTS AROUND CLOSEST THAT HAVE ONLY TWO NEIGHBORS . .


c       LOST=1
c instead just put it at the closest node
       XP=X(CLOSEST)
       YP=Y(CLOSEST)

c      BEL1=0
c      BEL2=0
c      ICNT=0
c      DO I=1,12
c         K=NOD2EL(I,CLOSEST)
c         IF (K.EQ.0) GOTO 80
c         IF ( (NOC(3,K).EQ.0) .AND. (BEL1.EQ.0) ) BEL1=K
c         IF ( (NOC(3,K).EQ.0) .AND. (BEL1.NE.0) ) BEL2=K
c      END DO
c 80   CONTINUE
C . . IF THE CLOSEST NODE IS NOT ON THE BOUNDARY...  
C      IF (BEL2.EQ.0) THEN
C         IF (BEL1.EQ.0) THEN
C             LOST=1
C             GOTO 100
C         END IF
C . . . .FIND THE POSITIONS OF THE OTHER TWO NODES IN BEL1          
         
C      END IF

C . . FIND THE NODES WHICH ARE ON EITHER SIDE OF CLOSEST ON THE BOUNDARY            
C      DO I=1,3
C         NODS1(I)=NOC(I,BEL1)
C         NODS2(I)=NOC(I,BEL2)
C      END DO
      
         
      
 
 
 100  CONTINUE
      RETURN
      END SUBROUTINE
C======================================================================   

C||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

C______________________________________________________________________
C======================================================================
      SUBROUTINE READ_DATA(XP,YP,LOCAT,NOC,X,Y,EL2EL,NOD2EL,VX,VY,
     &    NE,NN,NP,NVTS,VTIMINC,RLSTIME,NWS,WFACTOR,ICS,SLAM0,SFEA0,
     &    metonly)
C----------------------------------------------------------------------
C     Reads initial particle positions, mapping tables for tracking
C     particles from one element to the next, mesh file, first datasets
C     from water current velocity and wind velocity files. 
C----------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER I,J,K,NWS,ICS
      INTEGER LOCAT(1),NOC(3,1),EL2EL(3,1),NOD2EL(12,1),NE,NN,NP,NVTS
      REAL*8 XP(1),YP(1),X(1),Y(1),VX(1),VY(1),Z,VTIMINC,RLSTIME(1)
      REAL*8 WFACTOR,WX,WY,SLAM0,SFEA0,SLAM,SFEA
      logical :: metonly ! .true. if only wind data should be read
      
      character*80 chari
      
      OPEN(UNIT=11,FILE='PARTICLES.INP')
      OPEN(UNIT=12,FILE='EL2EL.TBL')
      OPEN(UNIT=13,FILE='NODE2EL.TBL')
      OPEN(UNIT=14,FILE='FORT.14')
      OPEN(UNIT=64,FILE='FORT.64')
      IF (NWS.NE.0) OPEN(UNIT=74,FILE='FORT.74')
      Write(*,*)'reading data opened files'
C . . READ PARTICLE STARTING POSITIONS      
      DO I=1,13
         READ(11,*) 
      END DO
      
      DO I=1,NP
c         READ(11,*)XP(I),YP(I),RLSTIME(I),LOCAT(I)
         READ(11,*)SLAM,SFEA,RLSTIME(I),LOCAT(I)
         IF (ICS.EQ.2) THEN
            CALL CPPD(XP(I),YP(I),SLAM,SFEA,SLAM0,SFEA0)
         ELSE
            XP(I)=SLAM
            YP(I)=SFEA
         END IF
c         WRITE(*,100)'PARTICLE ',I,' WILL BEGIN AT ',XP(I),',',YP(I)
      END DO
 100  FORMAT(1X,A,I8,A,F16.6,A,F16.6)     
C . . READ NEIGHBOR TABLES
      DO I=1,NE
         READ(12,*) (EL2EL(J,I),J=1,3)
      END DO
      WRITE(*,*)
      WRITE(*,*)'EL2EL READ SUCCESSFULLY'
      WRITE(*,*)
      
      DO I=1,NN
         READ(13,*) (NOD2EL(J,I),J=1,12)
      END DO
      WRITE(*,*)'NOD2EL READ SUCCESSFULLY'
      WRITE(*,*)
      
C . . READ GRID INFORMATION
      READ(14,*)
      READ(14,*)
      DO I=1,NN      
c         READ(14,*)K,X(I),Y(I),Z
         READ(14,*)K,SLAM,SFEA,Z
         IF (ICS.EQ.2) THEN
            CALL CPPD(X(I),Y(I),SLAM,SFEA,SLAM0,SFEA0)
         ELSE
            X(I)=SLAM
            Y(I)=SFEA
         END IF
         
      END DO
      IF (K.NE.NN) THEN
        WRITE(*,*)'!!!! ERROR- NODE NUMBERS ARE NOT CONSECUTIVE !!!!'
        WRITE(*,*)'!!!!           STOPING EXECUTION             !!!!'
        STOP
      END IF
      
      DO I=1,NE
         READ(14,*)K,J,(NOC(J,I),J=1,3)
      END DO
      
C . . READ INITIAL VELOCITY DATA
      READ(64,*)
      READ(64,*)NVTS,K,VTIMINC 
      READ(64,*)
      IF (NWS.NE.0) THEN
         READ(74,*)
         READ(74,*) 
         READ(74,*)
      END IF
      IF (K.NE.NN) THEN
        WRITE(*,*)'!!!! FORT.64 AND FORT.14 HAVE DIFFERENT # NODE !!!!'
        WRITE(*,*)'!!!!           STOPING EXECUTION               !!!!'
        STOP
      END IF
      DO I=1,NN
         READ(64,*)K,VX(I),VY(I)
         IF (NWS.NE.0) THEN
            READ(74,*)K,WX,WY
            VX(i)=VX(i)+WFACTOR*WX
            VY(i)=VY(i)+WFACTOR*WY
         END IF  
      END DO
      WRITE(*,*)'VELOCITY DATA READ SUCCESSFULLY'
      WRITE(*,*)
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
      CLOSE(64)
      IF (NWS.NE.0) CLOSE(74)
      RETURN
      END SUBROUTINE
C======================================================================

C||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


C______________________________________________________________________
C======================================================================
      SUBROUTINE WRITE_DATA(NP,XP,YP,LOCAT,PID,TIME,ICS,SLAM0,SFEA0,
     &                                                          LOST)
C----------------------------------------------------------------------
C THIS SUBROUTINE WRITES OUT THE DATA AT THE OUTPUT TIMES
C ALL ARGUMENTS ARE INPUT
C----------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, intent(in) :: I,NP,LOCAT(1),PID(1),ICS,LOST(1),LCT
      REAL*8, intent(in)  ::  XP(1),YP(1),TIME,SLAM0,SFEA0,SLAM,SFEA
      
      DO I=1,NP
         LCT=LOCAT(I)
         IF (LOST(I).ne.0) LCT=0    
         IF (ICS.EQ.2) THEN
            CALL INVCPD(XP(I),YP(I),SLAM,SFEA,SLAM0,SFEA0)
            WRITE(15,101)PID(I),SLAM,SFEA,TIME,LCT
         ELSE    
            WRITE(15,100),PID(I),XP(I),YP(I),LCT
         END IF
      END DO
 100  format(I12,2f14.2,f14.2,I12)   
 101  format(I12,2f14.9,f14.2,I12)   
      RETURN
      END SUBROUTINE
C======================================================================




C      SUBROUTINE BND_LOCAT(EL2EL,NOD2EL,X,Y,NOC,XP,YP,LOCAT,FOUND)
C----------------------------------------------------------------------
C THIS SUBROUTINE WILL "MOVE" THE LOST PARTICLE PERPENDICULARLY
C BACK TO THE BOUNDARY AND ADJUST THE SCALAR XP, AND YP ACCORDINGLY
C----------------------------------------------------------------------
C      IMPLICIT NONE
      
C      INTEGER, EL2EL(3,1),NOD2EL(12,1),NOC(3,1),LOCAT(1),FOUND,I,J,K
C      DOUBLE PRECISION X(1),Y(1),XP,YP,

C______________________________________________________________________
C======================================================================
      SUBROUTINE READ_64_STDY(STDY_TIME,VX,VY,NWS,WFACTOR)
C----------------------------------------------------------------------
C THIS SUBROUTINE READS THROUGH THE FORT.64 FILE UNTIL IF FINDS A TIME 
C GREATER THAN OR EQUAL TO STDY_TIME, THEN RETURNS THE VELOCITY DATA
C FROM THAT TIME.
C----------------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER I,J,K,NN,NTS,NWS
      REAL*8 STDY_TIME,VX(1),VY(1),TIME
      REAL*8 WFACTOR,WX,WY
      
      OPEN(UNIT=64,FILE='FORT.64')
      IF (NWS.NE.0) OPEN(UNIT=74,FILE='FORT.74')
      
      READ(64,*)
      READ(64,*)NTS,NN
      IF (NWS.NE.0) THEN
          READ(74,*)
          READ(74,*)
      END IF
      DO J=1,NTS
         READ(64,*)TIME
         DO I=1,NN
            READ(64,*)K,VX(I),VY(I)
            IF (NWS.NE.0) THEN
               READ(74,*)K,WX,WY
               VX(i)=VX(i)+WFACTOR*WX
               VY(i)=VY(i)+WFACTOR*WY
            END IF  
         END DO
         IF (TIME.GE.STDY_TIME) GOTO 100
      END DO
C . . IF STDY_TIME WAS AFTER THE LAST TS IN FORT.64, USE THE LAST TS
      STDY_TIME=TIME      
      
      
 100  CONTINUE
C      CLOSE(64)
      RETURN
      END SUBROUTINE
C======================================================================  

C|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| READ_64_DYN(STDY_TIME,NN,VX,VY,VX2,VY2,0)

C======================================================================
      SUBROUTINE READ_64_DYN()
C----------------------------------------------------------------------
C THIS SUBROUTINE READS THE NEXT RECORD IN THE FORT.64 TO BE USED
C FOR VELOCITY INTERPOLATIONS IN TIME.  FIRST IS FLAG SET TO ZERO
C FOR THE FIRST TIME THIS SUBROUTINE IS CALLED SO THAT VX,XY ARE NOT
C OVER WRITTEN
C----------------------------------------------------------------------
      USE MAUREPARAMS
      IMPLICIT NONE           
      INTEGER I,J,K
      REAL*8  VX(1),VY(1),TIME,VX2(1),VY2(1)
      logical, save :: firstCall = .true.
     
C . . EXCEPT FOR FIRST TIME REPLACE OLD VELOCITIES WITH NEW ONES      
      IF (firstCall.eqv..false.) THEN
         DO I=1,NN
            VX(I)=VX2(I)
            VY(I)=VY2(I)      
         END DO
      ELSE
         firstCall = .false.
      END IF
      
C . . READ THE NEXT VELOCITY OUTPUT
      IF (firstCall.eqv..true.) THEN
         OPEN(UNIT=64,FILE='FORT.64')
         IF (NWS.NE.0) OPEN(UNIT=74,FILE='FORT.74')
         READ(64,*)
         READ(64,*)NTS,NN
         IF (NWS.NE.0) THEN 
            READ(74,*)
            READ(74,*)
         END IF
         TIME=-90.D0
         DO WHILE (TIME .LT. STDY_TIME)  
            READ(64,*,END=100)TIME
            write(*,*)'Reading Velocity Field at Time ',TIME            
            if (NWS.NE.0) READ(74,*)
            DO I=1,NN
               READ(64,*)K,VX(I),VY(I)
               IF (NWS.NE.0) THEN
                  READ(74,*)K,WX,WY
                  VX(i)=VX(i)+WFACTOR*WX
                  VY(i)=VY(i)+WFACTOR*WY
               END IF  
            END DO  
         END DO 
      END IF
      
      READ(64,*,END=100)TIME
      IF (NWS.NE.0) READ(74,*)
      
      write(*,*)'Reading Velocity Field at Time ',TIME
      DO I=1,NN
         READ(64,*)K,VX2(I),VY2(I)
         IF (NWS.NE.0) THEN
             READ(74,*)K,WX,WY
             VX2(I)=VX2(I)+WFACTOR*WX
             VY2(I)=VY2(I)+WFACTOR*WY
         END IF  
      END DO 
      
      RETURN      
 100  CONTINUE
      WRITE(*,*)'!! END OF FORT.64, STOPPING EXECUTION!!'
      CLOSE(64)
      IF (NWS.NE.0) CLOSE(74)

      STOP
      END SUBROUTINE
C======================================================================     
C the CPP routines below are from the ADCIRC, but modified to take degrees
C instead of radians
C******************************************************************************
C                                                                             *
C    Transform from lon,lat (lamda,phi) coordinates into CPP coordinates.     *
C    Lon,Lat must be in DEGREES.                                              *
C                                                                             *
C******************************************************************************

      SUBROUTINE CPPD(X,Y,RLAMBDA,PHI,RLAMBDA0,PHI0)
      IMPLICIT NONE
      REAL*8 X,Y,RLAMBDA,PHI,RLAMBDA0,PHI0,R
      REAL*8 DEG2RAD
      DEG2RAD=3.141592653589793D0/180.D0
      
      R=6378206.4d0
      X=R*(RLAMBDA-RLAMBDA0)*DEG2RAD*COS(PHI0*DEG2RAD)
      Y=(PHI-PHI0)*DEG2RAD*R
      RETURN
      END SUBROUTINE


C******************************************************************************
C                                                                             *
C    Transform from CPP coordinates to lon,lat (lamda,phi) coordinates        *
C    Lon,Lat is in DEGREES.                                                   *
C                                                                             *
C******************************************************************************

      SUBROUTINE INVCPD(XXCP,YYCP,RLAMBDA,PHI,RLAMBDA0,PHI0)
      IMPLICIT NONE
      REAL*8 XXCP,YYCP,RLAMBDA,PHI,RLAMBDA0,PHI0,R
      REAL*8 DEG2RAD
      DEG2RAD=3.141592653589793D0/180.D0
      
      R=6378206.4d0
      RLAMBDA=RLAMBDA0+XXCP/(R*COS(PHI0*DEG2RAD)) /DEG2RAD
      PHI=YYCP/R /DEG2RAD + PHI0
      RETURN
      END SUBROUTINE


c CPP(X(JKI),Y(JKI),SLAM(JKI),SFEA(JKI),SLAM0,SFEA0)
c INVCP(X(JKI),Y(JKI),SLAM(JKI),SFEA(JKI),SLAM0,SFEA0)

!-----------------------------------------------------------------------
!     S U B R O U T I N E   C H E C K   F I L E   E X I S T E N C E
!-----------------------------------------------------------------------
!     jgf: Just check for the existence of a file. 
!-----------------------------------------------------------------------
      subroutine checkFileExistence(filename, errorIO)
      implicit none
      character(*), intent(in) :: filename ! full pathname of file
      integer, intent(out) :: errorIO 
      logical :: fileFound    ! .true. if the file is present      
      !
      ! Check to see if file exists
      write(6,'("Searching for file ",a," ...")') trim(filename)
      inquire(file=trim(filename),exist=fileFound,iostat=errorIO)
      if (fileFound.eqv..false.) then
         write(6,'("The file ",A," was not found.")') trim(filename)
      else
         write(6,'("The file ",A," was found.")') trim(filename)
      endif
!-----------------------------------------------------------------------
      END SUBROUTINE checkFileExistence
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!      S U B R O U T I N E   F I N D   U N U S E D   U N I T   N U M B E R
!-----------------------------------------------------------------------
!      Dynamic unit numbers prevent collisions. 
!-----------------------------------------------------------------------
      integer function availableUnitNumber()
      implicit none
      logical :: unusedUnitNumberFound ! .true. if an unused unit number is available
      logical :: isOpen ! true if an i/o unit number is connected to an open file
      integer :: errorIO
      integer, parameter :: minUnitNumber = 10 ! lowest fortran i/o unit number
      integer, parameter :: maxUnitNumber = 999 ! highest fortran i/o unit number
      
      do availableUnitNumber = minUnitNumber,maxUnitNumber
         inquire(unit=availableUnitNumber,opened=isOpen,iostat=errorIO)
         if (errorIO.eq.0) then
            if (isOpen.eqv..true.) then
               cycle
            else
               unusedUnitNumberFound = .true.
               exit
            endif
         else
            ! an error occurred
            write(6,'("ERROR: Could not find i/o unit.")') errorIO
            stop
         endif
      end do 
      if (unusedUnitNumberFound.eqv..false.) then
         write(6,'("ERROR: Could not find an unused unit number.")')
         stop
      endif
!-----------------------------------------------------------------------
      end function availableUnitNumber
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E   O P E N  F I L E  F O R  R E A D
!-----------------------------------------------------------------------
!     jgf: Added general subroutine for opening an existing
!     file for reading. Includes error checking.
!-----------------------------------------------------------------------
      subroutine openFileForRead(lun, filename, errorIO)
      implicit none
      integer, intent(in) :: lun   ! fortran logical unit number
      character(*), intent(in) :: filename ! full pathname of file
      integer, intent(out) :: errorIO
      logical unitConnected !.true. if this lun is already being used
      !
      !  Check to see if file exists
      call checkFileExistence(filename, errorIO)
      if (errorIO.ne.0) then
         return
      endif
      !
      ! Check to see if the unit number is already in use
      inquire(unit=lun,opened=unitConnected)
      if (unitConnected.eqv..true.) then
         write(6,'("The i/o unit ",i0," is already connected.")') lun
         errorIO = 1
         return
      endif
      !
      ! Open existing file
      OPEN(lun,FILE=trim(filename),STATUS='OLD',ACTION='READ',
     & IOSTAT=errorIO)
      if (errorIO.ne.0) then
         write(6,'("Could not open the file ",A,".")') trim(filename)
      else
         write(6,'("The file ",A," was opened successfully.")') 
     &    trim(filename)
      endif
      return
!-----------------------------------------------------------------------
      END SUBROUTINE openFileForRead
!-----------------------------------------------------------------------

!----------------------------------------------------------------------
!           S U B R O U T I N E    D O W N C A S E  
!----------------------------------------------------------------------
!     @jasonfleming: return a downcased version of the input string 
!----------------------------------------------------------------------
      subroutine downcase(string)
      implicit none
      character(*), intent(inout) :: string ! character array to downcase
      integer :: asciiCode ! decimal ascii code for a particular character
      integer :: i ! character counter
      !
      ! go through the character array looking for ascii codes between
      ! 65 (uppercase A) and 90 (uppercase Z); replace these characters 
      ! with lowercase 
      do i=1,len_trim(string)
         asciiCode = ichar(string(i:i))
         ! modify uppercase alphabetic characters only
         if ((asciiCode.ge.65).and.(asciiCode.le.90)) then
            asciiCode = asciiCode + 32
            string(i:i) = char(asciiCode)
         endif
      end do
!----------------------------------------------------------------------
      end subroutine downcase
!----------------------------------------------------------------------

      subroutine echoCmdLineOpt(cmdlineopt,cmdlinearg)
      implicit none
      character(*), intent(in) :: cmdlineopt
      character(*), intent(in) :: cmdlinearg
      write(6,'(a,a,a,a,a)') "INFO: Processing ",
     & trim(cmdlineopt)," ",trim(cmdlinearg),"."
      end subroutine echoCmdLineOpt
