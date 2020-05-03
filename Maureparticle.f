
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
C
C       @jasonfleming:      
C       gfortran -o maureparticle.x Maureparticle.f
C       gfortran -g -O0 -Wall -ffixed-line-length-none -fbacktrace -fbounds-check -ffpe-trap=zero,invalid,underflow,overflow,denormal -o maureparticle.x Maureparticle.f
C        
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
C Copyright (C) 2007, 2008, 2013-2016, 2020 Nathan Dill
C Copyright (C) 2020 Jason Fleming
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
      INTEGER :: NVTS   ! number of datasets in fort.64 (according to fort.64 header, which is not reliable)
      INTEGER :: RK2,NSTEPS
      INTEGER :: DYN  ! 1 for time varying velocity, 0 for steady velocity
      INTEGER :: VELCNT,NWS
      INTEGER :: ICS ! coordinate system, 1=cartesian, 2=geographic
      integer :: ndset ! counter for velocity datasets

      INTEGER :: NN ! number of nodes in the mesh
      REAL*8, ALLOCATABLE :: X(:),Y(:)   ! mesh node coordinates (after conversion to cpp if needed)
      REAL*8, ALLOCATABLE :: VX(:),VY(:) ! velocity driving particles (includes wind contrib)
      REAL*8, ALLOCATABLE :: VX2(:),VY2(:) ! next velocity dataset
      REAL*8, ALLOCATABLE :: WX(:),WY(:) ! wind u and v velocity
      INTEGER, ALLOCATABLE :: NOD2EL(:,:) ! elements around each node

      INTEGER :: NE ! number of elements in the mesh
      INTEGER, ALLOCATABLE :: NOC(:,:)    ! element table from mesh file?
      INTEGER, ALLOCATABLE :: EL2EL(:,:)  ! elements sharing an edge with this element

      INTEGER :: NP  ! number of particles total (released and unreleased)
      LOGICAL, ALLOCATABLE :: LOST(:)   !.true. if particle has exited domain
      LOGICAL, ALLOCATABLE :: FOUND(:)  !.true. if particle was found in a known element
      LOGICAL, ALLOCATABLE :: ITRKING(:)!.true. if particle has been released
      INTEGER, ALLOCATABLE :: LOCAT(:) ! index of element containing this PID
      REAL*8, ALLOCATABLE :: RLSTIME(:)  ! release time for the particles 
      REAL*8, ALLOCATABLE :: XP(:),YP(:)   ! current particle positions/T
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: XXP(:),YYP(:) ! rk2 starting particle position        
      
      REAL*8 TS,RNTIM,OUTPER,TIME,OUTTIME,STDY_TIME,FRACTIME,NOW
      REAL*8 VTIMINC   ! time incr btw datasets in fort.64 file (according to the header in that file)
      REAL*8 RELEASEPER,VELTIME1,VELTIME2,EDDY_DIF,WFACTOR,SLAM0,SFEA0
      
      CHARACTER*80 DESC1,DESC2

      logical :: metonly  ! .true. if particles should track wind instead of water current
      character(len=1024) :: maureParameterInputFile ! particle tracking control plus initial particle location
      character(len=1024) :: meshFile            ! adcirc fort.14
      character(len=1024) :: maureParticleOutputFile ! locations v time
      character(len=1024) :: velocityFile 
      character(len=1024) :: windVelocityFile      
      character(len=1024) :: elementLookupTableFile
      character(len=1024) :: nodeLookupTableFile

      logical :: keepDryParticles
      logical :: diffuseDryParticles 
      logical :: steadyState

      contains     
      !________________________________________________________________
      !================================================================
      subroutine initialize()
      implicit none
      integer :: p ! particle loop

C . . GET INFO TO ALLOCATE ARRAYS     

      OPEN(UNIT=14,FILE=trim(meshFile))     
      READ(14,'(A80)')DESC1
      READ(14,*)NE,NN
      CLOSE(14)

      OPEN(UNIT=11,FILE=trim(maureParameterInputFile))      
      READ(11,'(A80)')DESC2
      READ(11,*)NP
      READ(11,*)TS
      READ(11,*)RNTIM
      READ(11,*)OUTPER
      READ(11,*)RK2
      READ(11,*)DYN
      READ(11,*)STDY_TIME
      READ(11,*)EDDY_DIF
      READ(11,*)NWS
      READ(11,*)WFACTOR
      READ(11,*)ICS
      READ(11,*)SLAM0,SFEA0
      WRITE(*,*)'NE,NN,NP',NE,NN,NP
                      
      CLOSE(11)
      !
      ! check input: emit an error message if the output interval is
      ! cannot be divided evenly by timestep
      if ( modulo(outper,ts).ne.0 ) then
         write(*,*) "ERROR: The output interval ",outper,
     &   " cannot be divided evenly by the timestep ",ts,".",
     &   " As a result, output would never be written.",
     &   " Please adjust either the timestep or output interval."
         stop 1
      endif
      !
      ! check input: if --metonly was specified, then nws must be nonzero
      if ((metonly.eqv..true.).and.(nws.eq.0)) then
         write(*,*) "ERROR: The --metonly option was specified ",
     &   " on the command line but the meteorological data parameter ",
     &   " NWS was set to zero. Please remove the --metonly option ",
     &   " from the command line or make the NWS parameter nonzero."
         stop 1
      endif      
   
      ! check if a steady state run was specified
      if (DYN.EQ.0) then
         steadyState=.true.
      endif
     
C . . FIGURE NUMBER OF TIMESTEPS. . . . . . . . . . . . . . . . . . . .
      NSTEPS=INT(RNTIM/TS)
      
C . . FIGURE TOTAL NUMBER OF PARTICLES FOR CONTINUOUS RELEASE . . . . . 
C                 AND ALLOCATE ARRAYS ACCORDINGLY                       
      ALLOCATE(NOC(3,NE),NOD2EL(12,NN),EL2EL(3,NE),LOCAT(NP),FOUND(NP))
      ALLOCATE( X(NN),Y(NN),VX(NN),VY(NN),VX2(NN),VY2(NN),
     &      ITRKING(NP),RLSTIME(NP),LOST(NP),XP(NP),YP(NP))    
      ALLOCATE(XXP(NP),YYP(NP))
      if (nws.ne.0) then
         allocate(wx(nn),wy(nn))
      endif
 
C . . INITIALIZE VECTORS. . . . . . . . . . . . . . . . . . . . . .   
      DO p=1,NP
         ITRKING(p)=.false.
         LOST(p)=.false.
         LOCAT(p)=0
         FOUND(p)=.false.
         VX(p)=0.d0        !nld assume velocity is zero at time zero 
         VY(p)=0.d0
      END DO

      wfactor = 0.d0
      !----------------------------------------------------------------
      end subroutine initialize
      !================================================================

      !________________________________________________________________
      !================================================================
      SUBROUTINE READ_DATA()
      IMPLICIT NONE    
      real*8 slam, sfea
      real*8 z
      integer :: l  ! line/node number/element number counter 
      integer :: p  ! particle loop
      integer :: e  ! element loop
      integer :: n  ! mesh node loop 
      
      OPEN(UNIT=11,FILE=trim(maureParameterInputFile))
      OPEN(UNIT=12,FILE=trim(elementLookupTableFile))
      OPEN(UNIT=13,FILE=trim(nodeLookupTableFile))
      OPEN(UNIT=14,FILE=trim(meshFile))
      write(*,*)'reading data opened files'
      
      DO l=1,13
         READ(11,*) ! skip over control parameters in PARTICLES.INP
      END DO
      
C . . READ PARTICLE STARTING POSITIONS      
      DO p=1,NP
c         READ(11,*)XP(I),YP(I),RLSTIME(I),LOCAT(I)
         READ(11,*)SLAM,SFEA,RLSTIME(p),LOCAT(p)
         ! convert lon lat to cpp coordinate system (meters)
         IF (ICS.EQ.2) THEN
            CALL CPPD(XP(p),YP(p),SLAM,SFEA,SLAM0,SFEA0)
         ELSE
            XP(p)=SLAM
            YP(p)=SFEA
         END IF
c         WRITE(*,100)'PARTICLE ',I,' WILL BEGIN AT ',XP(I),',',YP(I)
      END DO
! 100  FORMAT(1X,A,I8,A,F16.6,A,F16.6)     
C . . READ NEIGHBOR TABLES
      DO e=1,NE
         READ(12,*) (EL2EL(n,e),n=1,3)
      END DO
      WRITE(*,*)
      WRITE(*,*)'EL2EL READ SUCCESSFULLY'
      WRITE(*,*)
      
      DO n=1,NN
         READ(13,*) (NOD2EL(e,n),e=1,12)
      END DO
      WRITE(*,*)'NOD2EL READ SUCCESSFULLY'
      WRITE(*,*)
      
C . . READ GRID INFORMATION
      READ(14,*) ! skip comment line
      READ(14,*) ! skip ne and nn line
      ! read node table
      DO n=1,NN      
c         READ(14,*)K,X(I),Y(I),Z
         READ(14,*) l,SLAM,SFEA,Z  
         IF (ICS.EQ.2) THEN
            CALL CPPD(X(n),Y(n),SLAM,SFEA,SLAM0,SFEA0)
         ELSE
            X(n)=SLAM
            Y(n)=SFEA
         END IF         
      END DO
      !@jasonfleming TODO: this should not be a requirement? 
      IF (l.NE.NN) THEN
         WRITE(*,*)'!!!! ERROR- NODE NUMBERS ARE NOT CONSECUTIVE !!!!'
         WRITE(*,*)'!!!!           STOPPING EXECUTION             !!!!'
         STOP
      END IF
      !
      ! read element table
      DO e=1,NE
         READ(14,*) l,p,(NOC(n,e),n=1,3)
      END DO
      !----------------------------------------------------------------
      END SUBROUTINE READ_DATA
      !================================================================

      !================================================================
      SUBROUTINE READ_VEL(velEnd)
      !----------------------------------------------------------------------
      ! THIS SUBROUTINE READS THE NEXT RECORD IN THE FORT.64 TO BE USED
      ! FOR VELOCITY INTERPOLATIONS IN TIME.  FIRST IS FLAG SET TO ZERO
      ! FOR THE FIRST TIME THIS SUBROUTINE IS CALLED SO THAT VX,XY ARE NOT
      ! OVER WRITTEN
      !----------------------------------------------------------------------
      ! This subroutine assumes that the current and wind velocity data
      ! have the same start time, end time, and time increment
      !----------------------------------------------------------------------
      IMPLICIT NONE
      logical, intent(out) :: velEnd ! true if all velocity data have been read
      logical, save :: first = .true.
      character(2000) :: line
      integer :: n ! node counter

      velEnd = .false. 

      if (first.eqv..true.) then          
         !write(*,*) 'first' !jgfdebug
         first = .false.
         ndset=1
         TIME=-90.D0  !nld - guessing this dont matter as long as TIME .lt. STDY_IIME
C . .    READ VELOCITY METADATA
         if (metonly.eqv..false.) then
            !write(*,*) 'metonly false' !jgfdebug         
            OPEN(UNIT=64,FILE=trim(velocityFile))
            READ(64,*,END=100) line ! skip comment line
            !write(*,*) trim(line) !jgfdebug
            READ(64,*,END=100) NVTS,n,VTIMINC ! number of datasets, number of nodes, time increment of output 
         endif
         if (nws.ne.0) then
            open(unit=74,file=trim(windvelocityfile))
            read(74,*,end=100) ! skip comment
            read(74,*,end=100) nvts,n,vtiminc ! number of datasets, number of nodes, time increment of output 
         endif
         do while (time.lt.stdy_time)
         write(*,*)'b4 first call of READ_VEL time',time !nlddebug
            call read_vel_dataset(velEnd)
         write(*,*)'after first call of READ_VEL time',time !nlddebug
         end do
         if (velEnd.eqv..false.) then
            WRITE(*,*) 'VELOCITY DATA READ SUCCESSFULLY'
         endif
C . . EXCEPT FOR FIRST TIME REPLACE OLD VELOCITIES WITH NEW ONES      
      else
      ! VX,VY will get advanced within the call to read_vel_dataset
      !nld   do n=1,nn
      !nld      vx(n)=vx2(n)
      !nld      vy(n)=vy2(n)      
      !nld   end do
         write(*,*)'b4 later READ_VEL time',time !nlddebug
         call read_vel_dataset(velEnd)
         write(*,*)'after later READ_VEL time',time !nlddebug
      end if
      return
      !--------------------------------------------------------------
      ! jump to here when attempting to read past the end of the file      
 100  CONTINUE
      velEnd = .true.
      write(*,*) 'ERROR: Velocity file ended unexpectedly.'
      CLOSE(64)
      IF (NWS.NE.0) CLOSE(74)
      RETURN      
      !================================================================
      END SUBROUTINE READ_VEL
      !================================================================      

      !================================================================
      SUBROUTINE READ_VEL_DATASET(velEnd)
      ! nld: revise this subroutine so that it returns with 
      ! both VX,VY and VX2,VY2 that account for the various
      ! options to just do current, just met, or blended
      IMPLICIT NONE
      logical, intent(out) :: velEnd ! true if the file has ended
      integer :: n                   ! mesh node counter
      integer :: l                   ! line counter
      !character(2000) :: line
      !
      ! advance VX,VY
      VX(:)=VX2(:)
      VY(:)=VY2(:)
      if (metonly.eqv..false.) then 
            !READ(64,*,END=100) line ! skip comment line
            !write(*,*) trim(line) !jgfdebug
         read(64,*,end=100) time    
      endif
      if (nws.ne.0) then
         read(74,*,end=100) time
      endif         
      do n=1,nn
         if (metonly.eqv..false.) then
            read(64,*) l,vx2(n),vy2(n)
         endif
         if (nws.ne.0) then         
            read(74,*) l,wx(n),wy(n)
         endif
      end do
      ! set the particle-driving velocity according to the current,
      ! wind, or a combination, based on input parameter specifications
      ! ... the new current velocities are set to vx2 and vy2 in read_vel_dataset 
      if (metonly.eqv..true.) then
         vx2(:) = wx(:)
         vy2(:) = wy(:)
      else
         ! blend current velocity with wind velocity if this has been specified 
         if (nws.ne.0) then
            vx2(:)=vx2(:)+wfactor*wx(:)  ! apply wind percent (default 0%)
            vy2(:)=vy2(:)+wfactor*wy(:)
         endif
      endif
      write(6,fmt='(a,f9.1,a)',advance='no')'time[',TIME,']'            
      write(6,fmt='(a,i0,a)',advance='no') '[',ndset,'] '
      ndset=ndset+1  ! jgf: Increment the dataset counter
      RETURN
      !--------------------------------------------------------------
      ! jump to here when attempting to read past the end of the file      
 100  CONTINUE
      velEnd = .true.
      CLOSE(64)
      IF (NWS.NE.0) CLOSE(74)
      RETURN
      !================================================================
      END SUBROUTINE READ_VEL_DATASET
      !================================================================

      !================================================================
      SUBROUTINE EULER_STEP()
      !----------------------------------------------------------------------      
      ! PARTICLE POSITIONS (XP,YP) ARE INPUT AND OUTPUT, PARTICLE POSITIONS
      ! ARE MOVED BY TAKING A SIMPLE STEP USING EULER'S METHOD 
      ! I.E. VELOCITY*TIME=DISPLACEMENT
      ! INPUT INCLUDES NUMBER OF PARTICLES NP, TIMESTEP TS, NODE CONNECTIVITY
      ! TABLE [NOC], NODAL POSITIONS  {X} AND {Y} AND VELOCITIES {V} , AND 
      ! PARTICLE ELEMENTAL LOCATIONS {LOCAT}
      ! FOR DYNAMIC RUNS(DYN=1) {VX} AND {VY} ARE THE VELOCITIES AT THE 
      ! BEGINNING OF TIME INTERVAL AND {VX2} AND {VY2} ARE AT THE ENC OF THE
      ! INTERVAL AND THE VELOCITY IS INTERPOLATED LINEARLY IN BETWEEN USING
      ! FRACTIME AS THE FRACTION OF THE TIME INTERVAL FROM T1 TO T2.
      !----------------------------------------------------------------------
      IMPLICIT NONE
      integer :: p ! particle loop
      integer :: e ! element loop
      integer :: n ! node loop around an element
      DOUBLE PRECISION R,ANG
      DOUBLE PRECISION VXP,VYP !time+space interpolated velocity at particle position
      DOUBLE PRECISION VVX(3),VVY(3)
      LOGICAL L_ISWET

      write(*,*)'inEuler np=',np        !nlddebug
           
      DO p=1,NP
         IF (ITRKING(p).EQV..true.) THEN
            IF (LOST(p).EQV..false.) THEN
               e=LOCAT(p)            
!nld               IF (.not.steadyState) THEN
C . . . . . . . . INTERPOLATE VELOCITIES IN TIME FOR DYNAMIC SIMULATION 
                  DO n=1,3
                     VVX(n)=VX(NOC(n,e))+(VX2(NOC(n,e))-VX(NOC(n,e)))
     &                                                     *FRACTIME
                     VVY(n)=VY(NOC(n,e))+(VY2(NOC(n,e))-VY(NOC(n,e)))
     &                                                     *FRACTIME
                  END DO
!nld               ELSE   
!nldC . . . . . . . . USE THE STEADY STATE VELOCITIES                    
!nld                  DO n=1,3
!                     VVX(n)=VX(NOC(n,e))
!                     VVY(n)=VY(NOC(n,e))
!                  END DO
!nld               END IF
C . . . . . . .GET THE X-VELOCITY AT THE PARTICLE POSITION
               CALL VELINTRP(X(NOC(1,e)),Y(NOC(1,e)),VVX(1),
     &                 X(NOC(2,e)),Y(NOC(2,e)),VVX(2),
     &                 X(NOC(3,e)),Y(NOC(3,e)),VVX(3),     
     &                         XP(p),YP(p),VXP)
 
C . . . . . . .GET THE Y-VELOCITY AT THE PARTICLE POSITION      
               CALL VELINTRP(X(NOC(1,e)),Y(NOC(1,e)),VVY(1),
     &                 X(NOC(2,e)),Y(NOC(2,e)),VVY(2),
     &                 X(NOC(3,e)),Y(NOC(3,e)),VVY(3),     
     &                         XP(p),YP(p),VYP)   
         !write(*,*) p,vxp,vyp !jgfdebug
         write(*,*) 'p,xp,yp',p,XP(p),YP(P) !nlddebug
         write(*,*) 'Euler vxp,vyp',VXP,VYP !nlddebug
         write(*,*) 'fractime, vyn1 vy2n1',fractime, VY(NOC(1,e)), 
     &                       VY2(NOC(1,e))
C . .    CALCULATE NEW POSITIONS USING VELOCITY AT MIDPOINT 
C . .    AND ORIGINAL POSITION
C . .    DON'T DIFFUSE IF IN A DRY ELEMENT
C . .    DO STILL ALLOW MOVEMENT IN A DRY ELEMENT WITHOUT DIFFUSION DEPENDING
C . .    ON THE VELOCITY OF ANY WET NODES IN THE ELEMENT, WHICH SHOULD
C . .    HELP PARTICLES FROM GETTING STUCK ON A WET/DRY BOUNDARY
C . .    ASSUME WE'RE IN A DRY ELEMENT IF ANY OF THE NODES HAVE VELOCITY LESS 
C . .    THAN THE MACHINE PRECISION (FROM EPSILON FUNCTION)
               L_ISWET=.TRUE.
               DO n=1,3
                  IF ((ABS(VVX(n)).LE.EPSILON(VVX(n))).AND.
     &                  (ABS(VVY(n)).LE.EPSILON(VVY(n)))) THEN
                     L_ISWET=.FALSE.
                  END IF
               END DO
          write(*,*)'L_ISWET',L_ISWET ! nlddebug
               !
               ! @jasonfleming: don't mark them as lost if the associated 
               ! command line option was set
               if (keepDryParticles.eqv..false.) then
                  IF (.NOT.L_ISWET) LOST(p)=.true.
               endif

               if (L_ISWET.or.diffuseDryParticles) then
                  IF (EDDY_DIF.GT.0.D0) THEN 
                     CALL RANDOM_NUMBER(R)
                     CALL RANDOM_NUMBER(ANG)
                     ANG=6.28318530717959D0*ANG
                     R=R*(EDDY_DIF*TS)**0.5D0
                     XP(p)=XP(p)+VXP*TS + R * COS(ANG)
                     YP(p)=YP(p)+VYP*TS + R * SIN(ANG)
           write(*,*)'Euler Diffusing, p,x,y',p,XP(p),YP(p) !nlddebug
                  ELSE
                     XP(p)=XP(p)+VXP*TS 
                     YP(p)=YP(p)+VYP*TS 
           write(*,*)'Euler NOT Diffusing, p,x,y',p,XP(p),YP(p) !nlddebug
                  END IF
               ENDIF
            END IF !LOST
         END IF !TRACKING
      END DO   
      
      RETURN
      END SUBROUTINE EULER_STEP
      !================================================================

      !================================================================
      SUBROUTINE RK2_STEP()
      !----------------------------------------------------------------------
      ! THIS SUBROUTINE USES THE 2ND ORDER RUNGE-KUTTA INTEGRATION METHOD 
      ! TO ADVANCE THE PARTICLE POSTIONS OVER ONE TIME STEP, IN ADDITION
      ! TO THE INPUT ARGUMENTS NEEDED FOR EULER_STEP() THIS ONE ALSO
      ! NEEDS THE EL2EL AND NOD2EL TABLES.
      !----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL L_ISWET
      DOUBLE PRECISION VVX(3),VVY(3)   ! time interpolated nodal velocities around an element
      DOUBLE PRECISION R,ANG
      DOUBLE PRECISION VXP,VYP ! time/space interpolated velocity at particle positions
      real*8 eddy_dif_temp

      integer :: p ! particle loop
      integer :: e ! element loop
      integer :: n ! node loop around an element
      integer :: eno ! jgfdebug: element where particle was not found
      
C . . SAVE THE STARTING POSITIONS 
      DO p=1,NP
         IF (LOST(p).EQV..false.) THEN
            XXP(p)=XP(p)
            YYP(p)=YP(p)
      write(*,*)'saving starting posion, p,xp,yp',p,XP(p),YP(p) !nlddebug
         END IF
      END DO
      
      WRITE(*,*)'STARTING POSITION SAVED'  !nlddebug
      
C . . DO AN EULER STEP AS A FIRST GUESS
C . . FOR THIS CALL DIFFUSIVITY IS SET TO ZERO
      eddy_dif_temp = eddy_dif
      eddy_dif = 0.d0
      CALL EULER_STEP()
      eddy_dif = eddy_dif_temp

      WRITE(*,*)'EULER GUESS COMPLETED' !nlddebug
     
C . . relocate particle to the midpoint of the Euler path
      DO p=1,NP
         IF (LOST(p).EQV..false.) THEN
            XP(p)=(XXP(p)+XP(p))*0.5d0
            YP(p)=(YYP(p)+YP(p))*0.5d0
         END IF
      END DO
      
CC      WRITE(*,*)'MIDPOINT CALCULATED'
C----------------------------------------------------------------------      
C . . CHECK LOCAT OF MIDPOINT
      DO p=1,NP 
          WRITE(*,'(a,i0,a,l,a,l)') 
     &              'particle ',p,' itrking ',
     &              itrking(p),' lost ',lost(p) !jgfdebug
         IF (ITRKING(p).EQV..true.) THEN
            IF (LOST(p).EQV..false.) THEN    
               CALL LOCAT_CHK(p,locat(p))
               !WRITE(*,*)'MIDPOINT CHECKED, PARTICLE ',p
               IF (FOUND(p).EQV..false.) THEN
                  eno=locat(p)
                  CALL UPDATE_LOCAT(p, locat(p))
                  WRITE(*,'(a,i0,a,i0,a,i0)') 
     &              'particle ',p,' updated element from ',
     &              eno,' to ',locat(p) !jgfdebug
                  if (lost(p).eqv..true.) then
                     XP(p)=-99999.d0
                     YP(p)=-99999.d0
                     locat(p)=0
                  endif
               ELSE
                 ! WRITE(*,*)'ELEMENT UNCHANGED, PARTICLE ',p !jgfdebug            
               END IF
            END IF
         END IF
      END DO 
C----------------------------------------------------------------------         
      
C . . GET THE VELOCITY AT THE MIDPOINT        
      DO p=1,NP
         IF (ITRKING(p).EQV..true.) THEN
            IF (LOST(p).EQV..false.) THEN   
               e=LOCAT(p)
!nld allow steady time to be between outputs just like the starting time IF (.not.steadyState) THEN
C . . . . . .INTERPOLATE VELOCITIES IN TIME FOR DYNAMIC SIMULATION 
                  DO n=1,3
                     VVX(n)=VX(NOC(n,e))+(VX2(NOC(n,e))-VX(NOC(n,e)))
     &                                                     *FRACTIME
                     VVY(n)=VY(NOC(n,e))+(VY2(NOC(n,e))-VY(NOC(n,e)))
     &                                                     *FRACTIME
                  END DO
!nld               ELSE   
C . . . . . .USE THE STEADY STATE VELOCITIES                    
!                  DO n=1,3
!                     VVX(n)=VX(NOC(n,e))
!                     VVY(n)=VY(NOC(n,e))
!                  END DO
!nld               END IF    
cc         WRITE(*,*)'IM AT LINE 359'
C . . . .GET THE X-VELOCITY AT THE MIDPOINT
               CALL VELINTRP(X(NOC(1,e)),Y(NOC(1,e)),VVX(1),
     &                 X(NOC(2,e)),Y(NOC(2,e)),VVX(2),
     &                 X(NOC(3,e)),Y(NOC(3,e)),VVX(3),     
     &                         XP(p),YP(p),VXP)
 
C . . . .GET THE Y-VELOCITY AT THE MIDPOINT      
               CALL VELINTRP(X(NOC(1,e)),Y(NOC(1,e)),VVY(1),
     &                 X(NOC(2,e)),Y(NOC(2,e)),VVY(2),
     &                 X(NOC(3,e)),Y(NOC(3,e)),VVY(3),     
     &                         XP(p),YP(p),VYP)   

         write (*,*)'in RK2 midpoint vxp,vyp',VXP,VYP !nlddebug       
C . .    CALCULATE NEW POSITIONS USING VELOCITY AT MIDPOINT AND ORIGINAL POSITION
C . .    DON'T DIFFUSE IF IN A DRY ELEMENT
               ! 
               ! @jasonfleming: seems like diffusion in a dry
               ! element should be constrained to be orthogonal to 
               ! (and in the direction of) the nearest wet edge 
               ! 
C . .    DO STILL ALLOW MOVEMENT IN A DRY ELEMENT WITHOUT DIFFUSION DEPENDING
C . .    ON THE VELOCOITY OF ANY WET NODES IN THE ELEMENT, WHICH SHOULD
C . .    HELP PARTICLES FROM GETTING STUCK ON A WET/DRY BOUNDARY
C . .    ASSUME WE'RE IN A DRY ELEMENT IF ANY OF THE NODES HAVE VELOCITY LESS 
C . .    THAN THE MACHINE PRECISION (FROM EPSILON FUNCTION)
               ! 
               ! @jasonfleming: alternatively we could load the fort.63 to see if 
               ! any of the 3 nodes are dry rather than look at velocity
               L_ISWET=.TRUE.
               DO n=1,3
                  IF ((ABS(VVX(n)).LE.EPSILON(VVX(n))).AND.
     &          (ABS(VVY(n)).LE.EPSILON(VVY(n)))) THEN
                     L_ISWET=.FALSE.
                  END IF
               END DO
               !
               ! @jasonfleming: don't mark them as lost if the associated 
               ! command line option was set
               if (keepDryParticles.eqv..false.) then
                  IF (.NOT.L_ISWET) LOST(p)=.true.
               endif

               if (L_ISWET.or.diffuseDryParticles) then
                  IF (EDDY_DIF.GT.0.D0) THEN 
                     CALL RANDOM_NUMBER(R)
                     CALL RANDOM_NUMBER(ANG)
                     ANG=6.28318530717959D0*ANG
                     R=R*(EDDY_DIF*TS)**0.5D0
                     XP(p)=XXP(p)+VXP*TS + R * COS(ANG)
                     YP(p)=YYP(p)+VYP*TS + R * SIN(ANG)
                  ELSE
                     XP(p)=XXP(p)+VXP*TS 
                     YP(p)=YYP(p)+VYP*TS 
                  END IF
               endif
            END IF
         END IF
      END DO                 
      
      RETURN
      END SUBROUTINE RK2_STEP
      !================================================================

      !================================================================
      SUBROUTINE LOCAT_CHK(p, e)
      !----------------------------------------------------------------------
      ! THIS SUBROUTINE CHECKS IF A PARTICLE RESIDES WITHIN THE ELEMENT LOCAT
      ! IT RETURNS THE VALUE FOUND=1 IF THE PARTICLE IS FOUND OR FOUND=0 IF 
      ! IT WAS NOT FOUND. LOCAT,XP,YP ARE SCALAR INPUT. 
      !----------------------------------------------------------------------
      IMPLICIT NONE      
      INTEGER :: e ! the element that the particle might be in
      INTEGER :: p ! the particle under consideration
      DOUBLE PRECISION DS1(2),DS2(2),DS3(2),CROSS,C1,C2,C3
      
      FOUND(p) = .false.
      IF (e.EQ.0) THEN
         RETURN 
      ENDIF
             
C . . GET DISPLACEMENTS FROM PARTICLE TO NODES    
      DS1(1)=X(NOC(1,e))-XP(p)
      DS1(2)=Y(NOC(1,e))-YP(p)
      DS2(1)=X(NOC(2,e))-XP(p)
      DS2(2)=Y(NOC(2,e))-YP(p)
      DS3(1)=X(NOC(3,e))-XP(p)
      DS3(2)=Y(NOC(3,e))-YP(p)

C . . ALL + CROSS PRODS. MEANS PART. IS FOUND IF NOC IS ANTI-CLOCKWISE      
      C1=CROSS(DS1,DS2)
      C2=CROSS(DS2,DS3)
      C3=CROSS(DS3,DS1)
       
      IF ((C1.GE.0).AND.(C2.GE.0).AND.(C3.GE.0)) FOUND(p)=.true.
      
      RETURN
      END SUBROUTINE LOCAT_CHK
      !================================================================
      
      !________________________________________________________________
      !================================================================
      SUBROUTINE VELINTRP(X0,Y0,V0,X1,Y1,V1,X2,Y2,V2,XP,YP,VP)
      ! THIS SUBROUTINE RETURNS (VP), THE VALUE AT POINT (XP,YP) LINEARLY
      ! INTERPOLATED ON A PLANE BETWEEN THE THREE POINTS (X0,Y0)(X1,Y1)(X2,Y2)
      IMPLICIT NONE
      DOUBLE PRECISION, intent(in) :: X0,Y0,V0,X1,Y1,V1,X2,Y2,V2,XP,YP
      DOUBLE PRECISION, intent(out) :: vp
      DOUBLE PRECISION :: T,U,DET,XX1,YY1,XX2,YY2,XXP,YYP 
      
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
      !----------------------------------------------------------------
      END SUBROUTINE VELINTRP
      !================================================================           

      !================================================================
      SUBROUTINE UPDATE_LOCAT(p, e)
      !----------------------------------------------------------------------
      ! THIS SUBROUTINE REQUIRES SCALAR INPUT OF PARTICLE POSITION AND LOCAT  
      ! IT SEARCHES THE EL2EL TABLE AND NOD2EL TABLE, 
      ! IT MAY RETURN A NEW VALUE LOCAT AND WILL RETURN LOST=.true. IF PARTICLE 
      ! LEAVES BOUNDARY OR SKIPS A NODE/ELEMENT GROUP
      !----------------------------------------------------------------------
      implicit none     
      integer, intent(in) :: p    ! particle under consideration
      integer, intent(inout) :: e ! element where particle is/was
      integer :: neighborElementIndex ! index of element around a particular node
      integer :: newElement ! element where particle might be
      integer :: closest    ! closest node to a particular point
      double precision ds(3)! distances of three nodes in an element to particular point

      !INTEGER BEL1,BEL2
      !INTEGER ICNT,NODS1(3),NODS2(3)
      !DOUBLE PRECISION DS(3),DSMIN,X1,Y1,X2,Y2
           
      ! has the particle moved to an element that shares an edge with 
      ! the input element e? 
C . . SEARCH EL2EL. . . . . . . . . . . . . . . . . . . . . . . . . . .
      DO neighborElementIndex=1,3
          newElement=EL2EL(neighborElementIndex,e) 
          IF (newElement.EQ.0) EXIT ! run out of candidates; assume zeroes to the right always?     
          CALL LOCAT_CHK(p, newElement)
          IF (FOUND(p).EQV..true.) THEN 
             locat(p)=newElement
             e=newElement
             RETURN
          END IF   
      END DO

C . . IF IT WAS NOT FOUND IN EL2EL SEARCH NOD2EL. . . . . . . . . . . .

C . . FIND THE CLOSEST NODE . . . . . . . . . . . . . . . . . . . . . .
      DS(1)=(X(NOC(1,e))-XP(p))**2 + (Y(NOC(1,e))-YP(p))**2
      DS(2)=(X(NOC(2,e))-XP(p))**2 + (Y(NOC(2,e))-YP(p))**2
      DS(3)=(X(NOC(3,e))-XP(p))**2 + (Y(NOC(3,e))-YP(p))**2
      closest=noc(minloc(ds,1),e)     
      !
      ! check all the elements around the closest node
      ! @jasonfleming: FIXME? Could this be improved by searching the 
      ! nod2el around each of the three nodes of the locat(p) element instead
      ! of just the closest node? 
      DO neighborElementIndex=1,12
         newElement=NOD2EL(neighborElementIndex,closest)
         IF (newElement.EQ.0) EXIT ! we've run out of candidates
         CALL LOCAT_CHK(p, newElement) 
         IF (FOUND(p).EQV..true.) THEN
            e=newElement
            locat(p)=newElement
            RETURN
         END IF
      END DO
      
C . . STILL NOT FOUND,. . . . . . . . . . . . . . . . . . . . . . . . .

      ! @jasonfleming: try just searching every element like we did 
      ! in the initial particle search ... if we don't find it that
      ! way, then mark it as permanently lost
      write(*,*) 'searching for particle ',p !jgfdebug
      do e=1,ne
         call locat_chk(p,e) ! particle, element 
         if (found(p).eqv..true.) then
            locat(p)=e
            write(*,*) 'particle ',p,' was found in element ',e !jgfdebug
            exit
         end if
      end do
      if (found(p).eqv..false.) then
         write(*,*) 'particle ',p,' cannot be found' !jgfdebug
         locat(p)=0
         xp(p)=-99999.d0
         yp(p)=-99999.d0
         lost(p)=.true.         
      endif

C . . IT MUST HAVE LEFT BOUNDARY OR SKIPPED AN ELEMENT NODE GROUP . . . 
C . . FIND THE ELEMENTS AROUND CLOSEST THAT HAVE ONLY TWO NEIGHBORS . .

c       LOST=1
c instead just put it at the closest node
       !XP=X(CLOSEST)
       !YP=Y(CLOSEST)

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

      RETURN
      END SUBROUTINE UPDATE_LOCAT
      !================================================================

      !________________________________________________________________
      !================================================================
      SUBROUTINE WRITE_DATA()
      !----------------------------------------------------------------
      ! THIS SUBROUTINE WRITES OUT THE DATA AT THE OUTPUT TIMES
      !----------------------------------------------------------------
      IMPLICIT NONE    
      logical, save :: first = .true. 
      integer :: p   ! particle loop
      integer :: lct ! element where the particle is found (if any)
      real*8 slam, sfea 
      !
      ! upon first call, create a new particle output file
      if (first) then
         OPEN(UNIT=15,FILE=trim(maureParticleOutputFile),action='write',
     &           status='replace')
         first = .false.
      else
         ! after the first call, just append to the end of the existing 
         ! output file
         OPEN(UNIT=15,FILE=trim(maureParticleOutputFile),action='write',
     &           status='old',position='append')      
      endif
      
      DO p=1,NP
         LCT=LOCAT(p)
         IF (LOST(p).eqv..true.) LCT=0    
         IF (ICS.EQ.2) THEN
            CALL INVCPD(XP(p),YP(p),SLAM,SFEA,SLAM0,SFEA0)
            WRITE(15,101) p,SLAM,SFEA,TIME,LCT
         ELSE    
            WRITE(15,*) p,XP(p),YP(p),TIME,LCT
         END IF
      END DO
      close(15)
! 100  format(I12,2f14.2,f14.2,I12)   
 101  format(I12,2f14.9,f14.2,I12)   
      RETURN
      !----------------------------------------------------------------
      END SUBROUTINE WRITE_DATA
      !----------------------------------------------------------------
      
      !================================================================
      END MODULE MAUREPARAMS
      !================================================================


C______________________________________________________________________
C======================================================================
      PROGRAM MAUREPARTICLE
      USE MAUREPARAMS
      IMPLICIT NONE
C---------------------------- VARIABLES -------------------------------           
      integer :: argcount ! number of command line arguments
      character(len=1024) :: cmdlineopt ! command line option
      character(len=1024) :: cmdlinearg ! command line option
      integer :: i  ! cmd line option counter
      integer :: p  ! particle loop
      integer :: e  ! element loop
      integer :: s  ! timestepping loop
      logical :: velEnd ! true when all velocity data has been read
C----------------------------------------------------------------------      
      ! set default file names and anything else that can be overridden by cmd 
      ! line arguments
      meshFile = 'FORT.14'
      maureParameterInputFile = 'PARTICLES.INP'
      maureParticleOutputFile = 'MAUREPT.OUT'      
      elementLookupTableFile = 'EL2EL.TBL'
      nodeLookupTableFile = 'NODE2EL.TBL'
      velocityFile = 'FORT.64'
      windVelocityFile = 'FORT.74'
      keepDryParticles = .false.
      diffuseDryParticles = .false.
      metonly = .false.
      steadyState = .false.
      !
      ! Get command line options
      argcount = command_argument_count() ! count up command line options
      if (argcount.gt.0) then
         i=0
         do while (i.lt.argcount)
            i = i + 1
            call getarg(i, cmdlineopt)
            call downcase(cmdlineopt)
            select case(trim(Cmdlineopt))
               case("--keep-dry-particles")
                  keepDryParticles = .true.
                  cmdlinearg = "true"
                  call echoCmdLineOpt(cmdlineopt,cmdlinearg)
               case("--diffuse-dry-particles")
                  diffuseDryParticles = .true.
                  cmdlinearg = "true"
                  call echoCmdLineOpt(cmdlineopt,cmdlinearg)                  
               case("--metonly")
                  metonly = .true.
                  cmdlinearg = "true"
                  call echoCmdLineOpt(cmdlineopt,cmdlinearg)                  
               case("--meshfile")
                  i = i + 1
                  call getarg(i, cmdlinearg)
                  meshFile = trim(cmdlinearg)
                  call echoCmdLineOpt(cmdlineopt,cmdlinearg)                  
               case("--velocityfile")
                  i = i + 1
                  call getarg(i, cmdlinearg)
                  call echoCmdLineOpt(cmdlineopt,cmdlinearg)
                  velocityFile = trim(cmdlinearg)
               case("--elementlookuptablefile")
                  i = i + 1
                  call getarg(i, cmdlinearg)
                  call echoCmdLineOpt(cmdlineopt,cmdlinearg)
                  elementLookupTableFile = trim(cmdlinearg)
               case("--nodelookuptablefile")
                  i = i + 1
                  call getarg(i, cmdlinearg)
                  call echoCmdLineOpt(cmdlineopt,cmdlinearg)
                  nodeLookupTableFile = trim(cmdlinearg)
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
                  maureParameterInputFile = trim(cmdlinearg)           
               case default
                  write(6,'(a,a,a)') "WARNING: Command line option '",
     &             TRIM(cmdlineopt),"' was not recognized."
            end select
         end do
      end if      
      
      ! allocate memory, set default variable values, read parameters
      ! input file, etc 
      call initialize()
      
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

      CALL READ_DATA()
      WRITE(*,*)'ALL INITIAL INPUT DATA HAS BEEN READ SUCCESSFULLY'

      ! Initial Read of the velocity data
      velEnd = .false.      
      CALL READ_VEL(velEnd)
      IF (.not.steadyState) THEN
         WRITE(*,*)'VELOCITY SOLUTION WILL BE READ DYNAMICALLY'
         WRITE(*,*)'STARTING AT ',STDY_TIME
         !nld CALL READ_VEL(velEnd) - moved above IF
         !CALL READ_VEL(velEnd)  !nld thinks this doesn't need to be called twice anymore. 
         if (velEnd) then
            write(*,*) 'ERROR: Could not read velocity data.'
            stop
         else
            WRITE(*,*)'INITIAL VELOCITY DATA HAS BEEN READ SUCCESSFULLY'
         endif
c        CALL READ_64_STDY(STDY_TIME,VX,VY)
C        VELCNT=1
      ELSE
         WRITE(*,*)'STEADY-STATE SIMULATION, READING FORT.64 UNTIL ',
     &              STDY_TIME    
         !nldCALL READ_VEL(velEnd) -moved above IF
         if (velEnd) then
            write(*,*) 'ERROR: Could not read velocity data.'
            stop
         else
            WRITE(*,*)'INITIAL VELOCITY DATA HAS BEEN READ SUCCESSFULLY'
         endif
         CLOSE(64)
         WRITE(*,*)'VELOCITY SOLUTION AT ',STDY_TIME,
     &             ' WILL BE TAKEN AS STEADY-STATE SOLUTION'
      END IF
      
C-------------------INITIAL SEARCH FOR PARTICLES-----------------------      

C . . SEARCH ALL ELEMENTS FOR INITIAL PARTICLE LOCATIONS. . . . . . . .
C . . THIS WILL BE SKIPPED IF LOCAT IS SPECIFIED IN PARTICLE.INP. . . .
      WRITE(*,*) 'SEARCHING FOR PARTICLES . . .'
      found(:)=.false.
      DO p=1,NP
         if (locat(p).ne.0) then
            ! check to make sure the element specified as the initial 
            ! location is actually correct
            call locat_chk(p, locat(p))
            if (found(p).eqv..false.) then
               write(*,*) 'WARNING: Particle ',p,' was not found in the'
     &           //' element specified as its initial location.'
            else
               ! it was found, go to the next one
               write(6,fmt='(a,i0,a)',advance='no') 'p',p,' '
               cycle
            endif
         else
            ! the particles input file for initial locations will normally
            ! contain zeroes for initial element locations 
            !WRITE(*,*) 'SEARCHING FOR PARTICLE ',p !jgfdebug
            do e=1,ne
               call locat_chk(p,e) ! particle, element 
               if (found(p).eqv..true.) then
                  locat(p)=e
                  !write(*,*) 'particle ',p,' was found in element ',e !jgfdebug
                  !write(*,*) 'locat(p)=',locat(p),' found(p)=',found(p) !jgfdebug
                  exit
               end if
            end do
         endif
         write(6,fmt='(a,i0,a)',advance='no') 'p',p,' '
      end do
      ! report any particles that were never found
      do p=1,np
         if ((locat(p).eq.0).or.(found(p).eqv..false.)) then
            !write(*,*) 'locat(p)=',locat(p),' found(p)=',found(p) !jgfdebug
            write(*,*) 'WARNING: Particle ',p,' cannot be found ',
     &         'so it will be marked as lost from the beginning.'
            lost(p)=.true.
            xp(p) = -99999.d0
            yp(p) = -99999.d0
         end if
      end do
C----------------------------------------------------------------------      
      
C----------------------THIS IS THE TRACKING LOOP-----------------------
C . . INITIALIZE SOME TIME KEEPING VARIABLES. . . . . . . . . . . . . .
      OUTTIME=STDY_TIME
      VELTIME2=TIME
      VELTIME1=VELTIME2-VTIMINC
      FRACTIME=(STDY_TIME-VELTIME1)/VTIMINC
      ! time was set to the time of the last read by read_vel_dataset
      ! initial calls, now set it back to the start time incase
      ! the start time was specified in between two current/wind output 
      ! times
      TIME=STDY_TIME  

      write(*,*) 'INFO: Starting particle tracking.'
C . . TRACKING LOOP . . . . . . . . . . . . . . . . . . . . . . . . . .      
      DO s=1,NSTEPS
c         WRITE(*,*)'TRACKING STEP ',J,' OF ',NSTEPS
C . . . .WRITE OUTPUT ? . . . . . . . . . . . . . . . . . . . . . . . .                           
         write(*,*) "time=",time,"    outtime=",outtime !jgfdebug
         IF (TIME.EQ.OUTTIME) THEN
           CALL WRITE_DATA()
           OUTTIME=OUTTIME+OUTPER
         END IF

C . . . .RELEASE MORE PARTICLES ?. . . . . . . . . . . . . . . . . . . .        
         DO p=1,NP      
            IF (TIME.GE.RLSTIME(p)) ITRKING(p)=.true.
         END DO
         
C . . . .ALL PARTICLES TAKE A STEP. . . . . . . . . . . . . . . . . . .                 
         IF (RK2.EQ.1) THEN
            CALL RK2_STEP()
         ELSE      
            CALL EULER_STEP()
         END IF   
         
C . . . . . UPDATE LOCAT IF NECESSARY . . . . . . . . . . . . .         
         do p=1,np
            CALL LOCAT_CHK(p,locat(p))
            ! if it is not in the element where we thought it was
            if (found(p).eqv..false.) then
               ! and it has not been marked permanently lost
               if (lost(p).eqv..false.) THEN
                  ! see if we can find it in a nearby element
                  call update_locat(p, locat(p))
               else
                  ! mark the coordinates as undefined if it is
                  ! permanently lost
                  XP(p)=-99999
                  YP(p)=-99999
               endif
            endif
         end do 
! 40      Format (A,I9,A,F13.3,A,I9)
C . . . .UPDATE SIMULATION TIME . . . . . . . . . . . . . . . . . . . .            
         TIME=TIME+TS
         
C . . . .CHECK TO SEE IF WE NEED TO GET NEW VELOCITY DATA         
         IF (.not.steadyState) THEN
            IF (TIME.GT.VELTIME2) THEN
               ! READ_VEL will call read_vel_dataset, which will set TIME 
               ! to the last read of the current/met file.
               ! Remember what the time is now so we can reset later
               NOW=TIME   
               CALL READ_VEL(velEnd)
               if (velEnd.eqv..true.) then
                  WRITE(*,*) 'INFO: Reached end of velocity file.'
                  stop
               else
                  ! these times correspond to fort.64 timing
                  ! VELTIME2 is the time that was just read
                  ! we can get it from the TIME variable just after read_vel() was called
                  ! VELTIME1 is the previous time.
                  VELTIME2=TIME      
                  VELTIME1=VELTIME2-VTIMINC 
               endif
               ! Since TIME got set to the last fort.64 time
               ! we need to reset it to now, which should be one timestep after VELTIME1 
               write(*,*)'TIMENOW: ', TIME, NOW !nlddebug
               TIME=NOW
            END IF
            FRACTIME=(TIME-VELTIME1)/VTIMINC
         END IF        
      END DO
C------------------------END OF TRACKING LOOP--------------------------         
      STOP
      END PROGRAM MAUREPARTICLE
C----------------------------------------------------------------------
C======================================================================

C||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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


C||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

C      SUBROUTINE BND_LOCAT(EL2EL,NOD2EL,X,Y,NOC,XP,YP,LOCAT,FOUND)
C----------------------------------------------------------------------
C THIS SUBROUTINE WILL "MOVE" THE LOST PARTICLE PERPENDICULARLY
C BACK TO THE BOUNDARY AND ADJUST THE SCALAR XP, AND YP ACCORDINGLY
C----------------------------------------------------------------------
C      IMPLICIT NONE
      
C      INTEGER, EL2EL(3,1),NOD2EL(12,1),NOC(3,1),LOCAT(1),FOUND,I,J,K
C      DOUBLE PRECISION X(1),Y(1),XP,YP,

C______________________________________________________________________


C|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| 

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
