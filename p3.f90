!!
!! Luke McCulloch
!! p3.f90
!! 6160 CFDpart3
!! 10-25-13
!! Potential Flow
!! Nonlinear Free Surface
!! Rankine Source
!! Adaptive Panel Method
!! 

!! Transom stern code:  Raven page 23 - A special form of the free surface condition is used on the transom fs panels.
!! Raven page 108 - hydrostatic pressure doesn't integrate to zero for transom flow
!!

PROGRAM p3


  use precise, only : defaultp           ! Double Precision
  use io                                 ! Input Output (read in fifi)
  use fifthpanel                         ! Get the Initial Panel Geometry vn, vt1, vt2
  use DynamicPanel                       ! Allow Dynamic Panel Geometry
  use computeA                           ! subroutines to create the A matrix   : uses A.f90
  use RHS                                ! compute the RHS vector b
  use computeNonLinearFlowProperties     ! module for Flow properties
  use vtkxmlmod                          ! Create Paraview files
  use slae                               ! invert the Dij "A" matrix
  

  IMPLICIT NONE    

  integer, parameter :: WP=defaultp
  
  INTEGER :: narg   
  INTEGER :: i, j, k, error, Row, Col
  INTEGER :: npoints, npanels, nfspanels, nfsl, nfst
  INTEGER :: ntransompanels, ntransoml, ntransomt
  INTEGER :: nhpl
  INTEGER :: iounit
  
  INTEGER                              :: clock_start, clock_end, clock_rate
  INTEGER                              :: nmax, iter                          ! simquit stuff
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: panels                              ! allocated in the io module

  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: PtsToPanelsConnectivity             ! point-panel connectivity
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: HullPtsUpdate                       ! store hull-fs panel connection
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: TransomPtsUpdate                    ! store hull-fs panel connection
  INTEGER                              :: DyDummy1                            ! Dynamic Panel Dummy Variable
  INTEGER                              :: locK, lock2, locD                   ! locations of max fs error   
  
  INTEGER :: nlIter

  CHARACTER(len=24) :: inputfile   ! file name can be max. 24 char. long
  CHARACTER(len=24) :: outputfile
  CHARACTER(len=24) :: datafile1
  CHARACTER(len=24) :: datafile2
  CHARACTER(len=24) :: datafile3
  CHARACTER(len=24) :: datafile4
  CHARACTER(len=24) :: datafile5
  CHARACTER(len=30) :: title
  CHARACTER(len=24) :: str_Fr      ! Froude number
  CHARACTER(len=24) :: str_Iter    ! number of iterations to perform (testing)

  LOGICAL :: flexists  ! a logical variable. I .true. or .false.


  REAL(WP), parameter :: gravity = 9.81  
  REAL(WP) :: length, Fr
  REAL(WP) :: elapsed_time

  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: points   !allocated in the io module

  ! Arrays for panel geometry
  REAL(WP), DIMENSION(3,4)                                             :: corners
  REAL(wp), ALLOCATABLE, DIMENSION(:,:)                                :: center
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)                              :: coordsys
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)                              :: cornerslocal
  REAL(wp), ALLOCATABLE, DIMENSION(:)                                  :: area
  REAL(wp), ALLOCATABLE, DIMENSION(:)                                  :: ldx

  ! Arrays for flow properties
  REAL(WP), DIMENSION(3)                                               :: vp
  REAL(WP), DIMENSION(3)                                               :: vpp
  REAL(WP), DIMENSION(3)                                               :: vppp
  REAL(WP), DIMENSION(3)                                               :: fieldpoint, fieldpointp, fieldpointpp
  REAL(WP), DIMENSION(3)                                               :: fieldpointppp
  REAL(WP), DIMENSION(3)                                               :: vinf

  REAL(WP),DIMENSION(3)                                                :: force


  ! Arrays for test point flow properties
  REAL(WP), DIMENSION(3)                                               :: vtottest
  REAL(WP), DIMENSION(3)                                               :: testpoint
  REAL(WP), DIMENSION(3)                                               :: vtp
  REAL(WP), DIMENSION(3)                                               :: vtpp
  REAL(WP), DIMENSION(3)                                               :: vtppp
  REAL(WP), DIMENSION(3)                                               :: testpointp, testpointpp
  REAL(WP), DIMENSION(3)                                               :: testpointppp  
  REAL(WP), DIMENSION(1)                                               :: cpt
  !REAL(WP) :: testi, testj, testk

  ! Arrays for matrix solution:
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)                                :: am
  REAL(WP), ALLOCATABLE, DIMENSION(:)                                  :: b
  REAL(WP), ALLOCATABLE, DIMENSION(:)                                  :: sigma


  ! Allocatable Arrays for flow properties
  REAL(WP), ALLOCATABLE, DIMENSION(:,:,:,:)                            :: hessglobal
  REAL(WP), ALLOCATABLE, DIMENSION(:,:,:)                              :: vel

  ! Allocatable Storage for old step flow properties (NLP)
  !REAL(WP), ALLOCATABLE, DIMENSION(:,:,:,:)                            :: hessglobalOld
  !REAL(WP), ALLOCATABLE, DIMENSION(:,:,:)                              :: velOld
  REAL(WP), ALLOCATABLE, DIMENSION(:,:,:)                              :: ddphi
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)                                :: dphi
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)                                :: dphi_temp
  

  ! More flow properties
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)                                :: vtotal
  REAL(WP), ALLOCATABLE, DIMENSION(:)                                  :: vn
  REAL(WP), ALLOCATABLE, DIMENSION(:)                                  :: vt1
  REAL(WP), ALLOCATABLE, DIMENSION(:)                                  :: vt2
  REAL(WP), ALLOCATABLE, DIMENSION(:)                                  :: cp




  REAL(WP), ALLOCATABLE, DIMENSION(:,:)                                :: velt
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)                                :: vtest  

  real(wp)                                                             :: deltax
  !!simquit stuff
  REAL(wp)                                                             :: Amatmax, sing, av, dx
  REAL(WP), DIMENSION(3,3)                                             :: hp ! mirror local hessian
  ! Wave stuff
  REAL(wp)                                                             :: cw, S, sourcesink, waveSOR
  REAL(WP), ALLOCATABLE,DIMENSION(:)                                   :: zeta
  REAL(WP), ALLOCATABLE,DIMENSION(:)                                   :: dzeta_dx
  REAL(WP), ALLOCATABLE,DIMENSION(:)                                   :: dzeta_dy
  REAL(WP), ALLOCATABLE,DIMENSION(:,:)                                 :: profilewave

  ! Nonlinear convergence:
  real(wp)                                                             :: errorD, errorK, errork2

  !  REAL(wp), ALLOCATABLE, DIMENSION(:) :: celldata1
  !  REAL(wp), ALLOCATABLE, DIMENSION(:) :: celldata2
  !  REAL(wp), DIMENSION(3) :: paneloffsets

  CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate
  CALL SYSTEM_CLOCK(COUNT=clock_start) ! Start timing
  


  WRITE(6,'(A)') 'Linear Free Surface Panel Method, v @.0, 4-17-10, Main Program, Luke McCulloch '

  ! Count the number command line arguments (actually Fortran 2003)
  narg = command_argument_count()

  WRITE(6, FMT='(AI3A)') 'we have ', narg, ' command line arguments '
  
  IF ( narg > 3 ) THEN
     CALL get_command_argument(1, inputfile)
     CALL get_command_argument(2, outputfile)
     CALL get_command_argument(3, str_Fr)
     CALL get_command_argument(4, str_Iter)

     READ(str_Fr,*) Fr
     READ(str_Iter,*) nlIter
  ELSE
     WRITE(6,'(A)') ' Output file name and/or Froude number and/or iteration count missing!'
     WRITE(6,'(A)') ' Usage:> ./p3  <inputfile> <outputfile> <FroudeNumber> <ItNumber> '
     STOP
  ENDIF

  INQUIRE(file=inputfile, exist=flexists)

  IF (flexists) THEN
  
    
     WRITE(6,'(AA)') ' input file,  ', inputfile
     WRITE(6,*) 'output file, ', outputfile
  
     !! from the module io
     !! ALLOCATE 2 ARRAYS : "points" and "panels"
     CALL input( inputfile, outputfile, title, &
                 npanels, nfspanels, ntransompanels, npoints, &
                 panels, points,i,j,deltax, nfsl, nfst, ldx, &
                 ntransoml, ntransomt,nhpl )




  Else
     ! If not, the program stops
     WRITE(6,'(A)') ' Input file does not exist!'
     STOP ! abort program
  ENDIF

  

  iounit=1
  CALL VtkXmlPolyDataCellScalar( iounit, 'TLM_boat.vtp', npoints, &
       0, 0, 0, npanels+nfspanels+ntransompanels, &
       RESHAPE(points, (/3*npoints/)), &
       RESHAPE((panels-1), (/4*(npanels+nfspanels+ntransompanels)/)), &
       scalar1=real(panels(1,:),kind=wp), &
       namescalar1='x')







  
  write(6,*) title
  write(6,*) i,j
  write(6,*) npanels, nfspanels, ntransompanels, npoints, deltax

  !! Panel Geometry:  
  ALLOCATE( center(3,npanels+nfspanels+ntransompanels) )
  ALLOCATE( coordsys(3,3,npanels+nfspanels+ntransompanels) )
  ALLOCATE( cornerslocal(2,4,npanels+nfspanels+ntransompanels) )
  ALLOCATE( area(npanels+nfspanels+ntransompanels) )
  


  DO i = 1, npanels+nfspanels+ntransompanels

     ! we store the panel corner points in a 3 by 4 array
     DO j=1, 4
        corners(:,j) = points(:,panels(j,i))
     END DO
 

     ! We call the panel geometry function IN THE  MODULE "fifthpanel"
     CALL panelgeometry( corners, coordsys(:,:,i), cornerslocal(:,:,i), &
                      area(i), center(:,i) ) 

  END DO 



  !! Allocate step 1 gradient and hessian Memory:
  ALLOCATE(vel(3,npanels+nfspanels+ntransompanels,npanels+nfspanels+ntransompanels))
  ALLOCATE(hessglobal(3,3,npanels+nfspanels+ntransompanels,npanels+nfspanels+ntransompanels))

  !! Allocate step 0 gradient and hessian memory:
  ALLOCATE(dphi(3,npanels+nfspanels+ntransompanels))
  ALLOCATE(dphi_temp(3,npanels+nfspanels+ntransompanels))
  ALLOCATE(ddphi(3,3,npanels+nfspanels+ntransompanels))

  !! Allocate Solver Memory:
  ALLOCATE(am(npanels+nfspanels+ntransompanels+2,npanels+nfspanels+ntransompanels+2))
  am = REAL(0.,kind=wp)
  ALLOCATE(b(npanels+nfspanels+ntransompanels))
  ALLOCATE(sigma(npanels+nfspanels+ntransompanels))

  ! Allocate space for flow data 
  !  Note: in some cases extra space is taken for consistent indexing
  ALLOCATE(vtotal(3,npanels+nfspanels+ntransompanels))
  ALLOCATE(vn(npanels+nfspanels+ntransompanels))
  ALLOCATE(vt1(npanels+nfspanels+ntransompanels))
  ALLOCATE(vt2(npanels+nfspanels+ntransompanels))
  ALLOCATE(cp(npanels+nfspanels+ntransompanels))
  ALLOCATE(zeta(npanels+nfspanels+ntransompanels))
  ALLOCATE(dzeta_dx(npanels+nfspanels+ntransompanels))
  ALLOCATE(dzeta_dy(npanels+nfspanels+ntransompanels))
  ALLOCATE(profilewave(2,nfsl))




  !!
  !! Base Flow Settings:
  !!
  waveSOR = 0.7
  length = 1.0
  !Fr = .45 - user input now
  vinf = (/0.,0.,0./)
  vinf(1) = Fr*sqrt(length*gravity)
  print*, 'ub = ',vinf(1)

  ! Set base flow equal to the flow at infinity:
  dphi       = 0.0
  dphi_temp  = 0.0
  dphi(1,:)  = -vinf(1) !start base flow
  ddphi      = 0.0


  vtotal = dphi
  !! normal & tangential components:
  do i=1,npanels+nfspanels+ntransompanels
     vn(i)=Dot_Product(coordsys(:,3,i),vtotal(:,i))
     vt1(i)=dot_product(coordsys(:,1,i),vtotal(:,i))
     vt2(i)=dot_product(coordsys(:,2,i),vtotal(:,i))
  end do

  iounit = 3
  CALL VtkXmlPolyDataCellScalar( iounit, 'TLM_LocalCoords.vtp', npoints, &
        0, 0, 0, npanels+nfspanels+ntransompanels, &
        RESHAPE(points, (/3*npoints/)), &
        RESHAPE((panels-1), (/4*(npanels+nfspanels+ntransompanels)/)), &
        scalar1=vn, &
        scalar2=vt1, &
        scalar3=vt2, &
        namescalar1='vn', &
        namescalar2='vt1', &
        namescalar3='vt2')
        



  ! Set initial wave elevation to mean water level:
  zeta = 0.
  dzeta_dx = 0.
  dzeta_dy = 0.

  !! Additional Control Parameters for the Simqit matrix inverter
  Amatmax = maxval(am)                  !returns the largest value in am
  sing = REAL(1.e-6,kind=wp)*Amatmax    !make the matrix singular if any value exceeds Amatmax
  av = sqrt(1./REAL(npanels+nfspanels))
  dx = 0.000001
  nmax = npanels+nfspanels+2


  call computePtToPanelConnectivity(   npanels    , nfspanels, npoints,              &
                                       nfsl       , nfst     ,                       &           
                                       panels     , points   , PtsToPanelsConnectivity,&
                                       HullPtsUpdate, &
                                       ntransompanels, ntransoml, ntransomt,&
                                       nhpl, TransomPtsUpdate)
  write(*,*)  'point to panel connectivity done'

  !! print connectivity table:
  !do i=1,npoints
  !   DyDummy1 =  PtsToPanelsConnectivity(5,i)-1
  !   write(*,*) 'point',i,'       has',DyDummy1 ,'panels.  They are numbers:',PtsToPanelsConnectivity( 1:DyDummy1 ,i)
  !end do


  write(*,*) '-----------------------------------------------------------------------'
  write(*,*) 'Starting Solver:'

  do k=1,nlIter
     write(*,*) 'starting global solver iteration', k

     ! Set base flow equal to the flow at infinity: pair this with vinf(1) dot n(1) body b.c.
     ! cummlative flows get the updated base flow body b.c.
     dphi       = 0.0
     dphi_temp  = 0.0 
     dphi(1,:)  = -vinf(1) !start base flow
     ddphi      = 0.0

     call  computeBodyBC( npanels, nfspanels, fieldpoint,fieldpointp,center,&
                          vel,coordsys,area,cornerslocal,vp,am, &
                                       ntransompanels, ntransoml, ntransomt)
     write(*,*) 'completed body boundary conditions'

     !call computeLinearizedFreeSurfaceBC( npanels, nfspanels, fieldpoint,fieldpointp,center,deltax,&
     !     vel,coordsys,area,cornerslocal,vp,vinf, hp, hessglobal, am)  
     !write(*,*) 'completed free surface boundary conditions'
     
     call computeNonLinearFreeSurfaceBC( npanels, nfspanels, fieldpoint,fieldpointp,center,deltax,&
          vel,dphi,coordsys,area,cornerslocal,vp,vinf, hp, hessglobal, ddphi, &
          zeta, dzeta_dx, dzeta_dy, am, &
                                       ntransompanels, ntransoml, ntransomt, nhpl)
     write(*,*) 'completed free surface boundary conditions'
     
     
     !call  computeRHS(npanels, nfspanels, coordsys, vinf, b)
     !write(*,*) 'completed RHS'
     
     call computeNonLinearRHS(npanels, nfspanels, coordsys, vinf, dphi, ddphi, b, &
                                       ntransompanels, ntransoml, ntransomt,nhpl,center)
     write(*,*) 'completed RHS'
     
     
     write(*,*) ''
     write(*,*) 'calling matrix inversion routine....'
     call SIMQIT(TRANSPOSE(am), b, npanels+nfspanels, nmax, &
          sing, av, dx, sigma, iter)
     write (*,*) 'Matrix inversion complete'
     write(6,*) 'simquit iter', iter

     print*, 'sigma max = ',maxval(sigma)
     print*, 'sigma min = ',minval(sigma)


     print*, 'sigma max = ',maxval(abs(sigma))
     print*, 'sigma min = ',minval(abs(sigma))
     
     
     CALL computeVelocities(npanels,nfspanels,vtotal,vinf,sigma,vel,vn,vt1,vt2,coordsys, dphi_temp, dphi, &
                                       ntransompanels, ntransoml, ntransomt)
     
     CALL computeWaveHeight(npanels,nfspanels, vinf,zeta,vtotal, waveSOR, dphi_temp, dphi,center, nfsl, nfst,& 
          panels, points, PtsToPanelsConnectivity,HullPtsUpdate, &
                                       ntransompanels, ntransoml, ntransomt, nhpl,&
                                       TransomPtsUpdate)
     
     
     write(*,*) 're-computing local geometry and coordinates...'
     DO i = npanels-nhpl , npanels+nfspanels+ntransompanels
        DO j=1, 4
           corners(:,j) = points(:,panels(j,i))
        END DO
        CALL panelgeometry( corners, coordsys(:,:,i), cornerslocal(:,:,i), area(i), center(:,i) ) 
        !CALL panelgeometry2( corners, coordsys(:,:,i), cornerslocal(:,:,i), area(i) ) 
     END DO
     write(*,*) 'Done'  
     write(*,*) ''
  

     CALL computeDZeta(npanels, nfspanels, nfsl, nfst, zeta, dzeta_dx, dzeta_dy, ldx, center, &
                                       ntransompanels, ntransoml, ntransomt, nhpl) 


     
     
     
     
     CALL updateSurfaceVelocities( npanels, nfspanels, fieldpoint,fieldpointp,center,deltax,&
          vel,dphi, coordsys, area, cornerslocal,vp,vinf, hp, hessglobal, ddphi, &
          zeta, dzeta_dx, dzeta_dy, sigma, vtotal, &
                                       ntransompanels, ntransoml, ntransomt)


     CALL computeFreeSurfaceConvergence(npanels,nfspanels,dphi,zeta, dzeta_dx, dzeta_dy, vinf, nfsl,&
          errorK, errorD , locK, locD, lock2, errork2, &
                                       ntransompanels, ntransoml, ntransomt)
     write(*,*) 'kinematic error = ', errorK
     write(*,*) 'kinematic error 2 = ', errorK2
     write(*,*) 'dynamic error = ', errorD
     
     if (errorK <.002*(vinf(1)) .AND. errorD <.0025*vinf(1)**2) then
        exit
     end if

     locK=locK+npanels
     locD=locD+npanels
     
     print*, 'Print Free surface at max K error: '
     print*, 'i, zeta, dzeta/dx, dzeta/dy, dphi/dz, center'
     print*,  locK-2,zeta(locK-2), dzeta_dx(locK-2), dzeta_dy(locK-2), dphi(3,locK-2), center(1,locK-2) 
     print*,  locK-1,zeta(locK-1), dzeta_dx(locK-1), dzeta_dy(locK-1), dphi(3,locK-1), center(1,locK-1) 
     print*,  locK,zeta(locK), dzeta_dx(locK), dzeta_dy(locK),  dphi(3,locK), center(1,locK) 
     
     print*, 'Print Free surface at max D error: '
     print*, 'i, zeta, dzeta/dx, dzeta/dy, center'
     print*,  locD-2,zeta(locD-2), dzeta_dx(locD-2), dzeta_dy(locD-2), center(1,locD-2) 
     print*,  locD-1,zeta(locD-1), dzeta_dx(locD-1), dzeta_dy(locD-1), center(1,locD-1) 
     print*,  locD,zeta(locD), dzeta_dx(locD), dzeta_dy(locD), center(1,locD) 

  end do

  vtotal = dphi
  !! normal & tangential components:
  do i=1,npanels+nfspanels+ntransompanels
     vn(i)=Dot_Product(coordsys(:,3,i),vtotal(:,i))
     vt1(i)=dot_product(coordsys(:,1,i),vtotal(:,i))
     vt2(i)=dot_product(coordsys(:,2,i),vtotal(:,i))
  end do
  
  !! End of NLP loop until sinkage and trim become available
  
  write(*,*) 'finished loops, computing final flow properties...'

  CALL computePressures(npanels,nfspanels,vtotal,vinf,cp,ntransompanels)
  
  CALL checkSourceSink(npanels, nfspanels, sigma, area, sourcesink,ntransompanels)
    
  CALL computeHullForces(npanels,nfspanels,center,force,area,vinf,cp,coordsys)
  
  CALL computeWaveResistance(npanels,nfspanels,center,S,area,vinf,force,cw)
  
  CALL computeWaveProfile(npanels,nfspanels,profilewave,center,deltax,zeta,nfsl)
  !!
  !!***************************************************************

  call  setPanelLocation(npanels,nfspanels,&
                               vinf,zeta,vtotal, &
                               waveSOR, dphi_temp,dphi, &
                               center, nfsl, nfst,&
                               panels, points,&
                               PtsToPanelsConnectivity,&
                               HullPtsUpdate,&
                               ntransompanels, ntransoml, ntransomt,nhpl,TransomPtsUpdate)

  CALL SYSTEM_CLOCK(COUNT=clock_end) ! Stop timing
  ! Calculate the elapsed time in seconds:
  elapsed_time=REAL((clock_end-clock_start)/clock_rate,WP)
  write(*,*) ''
  write(*,*) 'total time = ',elapsed_time

  !!Write some important data to the output file
  !!----------------------------------------------------------------------
  CALL output( outputfile, iter, sourcesink, cw, profilewave )
  !!----------------------------------------------------------------------
  
  !! Lets try to store the waves on each free surface panel...
  !!***************************************************************

    ! First we store data valid for a panel (cell in VTK jargon)

  ! Create a file name with the appropriate extension
   i = INDEX(outputfile, ".", BACK=.TRUE.)
   IF (i == 0) THEN
      WRITE(datafile1, '(A,A,A)')  TRIM(outputfile), 'wave', '.vtp'
   ELSE
      IF (i > 20) i = 20 ! To make sure the filename does not exceed 24 char
      WRITE(datafile1, '(A,A,A)')  outputfile(1:i-1), 'wave', '.vtp'
   ENDIF 
   PRINT*, datafile1


   iounit = 10
   CALL VtkXmlPolyDataCellScalar( iounit, datafile1, npoints, &
         0, 0, 0, npanels+nfspanels+ntransompanels, &
         RESHAPE(points, (/3*npoints/)), &
         RESHAPE((panels-1), (/4*(npanels+nfspanels+ntransompanels)/)), &
         scalar1=(zeta), &
         scalar2=sigma, &
         namescalar1='zeta', &
         namescalar2='sigma')





  !! Lets try to store the waves on each just the y=0 free surface panels...
  !!***************************************************************

    ! First we store data valid for a panel (cell in VTK jargon)

  ! Create a file name with the appropriate extension
   i = INDEX(outputfile, ".", BACK=.TRUE.)
   IF (i == 0) THEN
      WRITE(datafile2, '(A,A,A)')  TRIM(outputfile), 'waveprofile', '.vtp'
   ELSE
      IF (i > 20) i = 20 ! To make sure the filename does not exceed 24 char
      WRITE(datafile2, '(A,A,A)')  outputfile(1:i-1), 'waveprofile', '.vtp'
   ENDIF 
   PRINT*, datafile2


   print*, nfsl
   iounit = 10
   CALL VtkXmlPolyDataCellScalar( iounit, datafile2, npoints, &
         0, 0, 0, npanels+nfsl, &
         RESHAPE(points, (/3*npoints/)), &
         RESHAPE((panels-1), (/4*(npanels+nfsl)/)), &
         scalar1=zeta, &
         namescalar1='zeta')







  !! Lets try to store the cp on each hull panel...
  !!***************************************************************

    ! First we store data valid for a panel (cell in VTK jargon)

  ! Create a file name with the appropriate extension
   i = INDEX(outputfile, ".", BACK=.TRUE.)
   IF (i == 0) THEN
      WRITE(datafile3, '(A,A,A)')  TRIM(outputfile), 'scalar', '.vtp'
   ELSE
      IF (i > 20) i = 20 ! To make sure the filename does not exceed 24 char
      WRITE(datafile3, '(A,A,A)')  outputfile(1:i-1), 'scalar', '.vtp'
   ENDIF 
   PRINT*, datafile3


   iounit = 10
   CALL VtkXmlPolyDataCellScalar( iounit, datafile3, npoints, &
         0, 0, 0, npanels, &
         RESHAPE(points, (/3*npoints/)), &
         RESHAPE((panels-1), (/4*(npanels)/)), &
         scalar1=cp, &
         namescalar1='cp')       



  !!******************************************************************

  !! Now the total velocity vector on each panel...
  !!***************************************************************
   ! Create a file name with the appropriate extension
  i = INDEX(outputfile, ".", BACK=.TRUE.)
  IF (i == 0) THEN
     WRITE(datafile4, '(A,A,A)')  TRIM(outputfile), 'vtotal', '.vtp'
  ELSE
     IF (i > 20) i = 20 ! To make sure the filename does not exceed 24 char
     WRITE(datafile4, '(A,A,A)')  outputfile(1:i-1), 'vtotal', '.vtp'
  ENDIF 
  PRINT*, datafile4
  

  iounit = 10
  CALL VtkXmlSaveVectorFlowData( iounit, datafile4, npanels+nfspanels+ntransompanels, center, &
        vtotal) 

! This method probably won't work for later versions of the program!

! vt1 and vt2 must be re-converted into vectors to display as vectors
! multiply coordsys()times the scalar vt1 or vt2 to scale the coordsys
! vector.  Then they should work.

  !!******************************************************************


  !! Now the vector data on each panel...
  !!***************************************************************
   ! Create a file name with the appropriate extension
  i = INDEX(outputfile, ".", BACK=.TRUE.)
  IF (i == 0) THEN
     WRITE(datafile5, '(A,A,A)')  TRIM(outputfile), 'components', '.vtp'
  ELSE
     IF (i > 20) i = 20 ! To make sure the filename does not exceed 24 char
     WRITE(datafile5, '(A,A,A)')  outputfile(1:i-1), 'components', '.vtp'
  ENDIF 
  PRINT*, datafile5
  

  iounit = 10
  CALL VtkXmlPolyDataCellScalar( iounit, datafile5, npoints, &
        0, 0, 0, npanels+nfspanels+ntransompanels, &
        RESHAPE(points, (/3*npoints/)), &
        RESHAPE((panels-1), (/4*(npanels+nfspanels+ntransompanels)/)), &
        scalar1=vn, &
        scalar2=vt1, &
        scalar3=vt2, &
        namescalar1='vn', &
        namescalar2='vt1', &
        namescalar3='vt2')
        


! Face specific
! Scalar components of the velocity vector
! correct format
  !!******************************************************************


  ! De-allocate all arrays


!!$  write(*,*) 'cornerslocal memory:' 
!!$  write(*,*) 'cornerslocal, 1st element:'
!!$  !write(*,*) cornerslocal(1,1,1)
!!$  IF (ALLOCATED(cornerslocal)) then
!!$     write(*,*) 'cornerslocal is allocated'
!!$     DEALLOCATE(cornerslocal)
!!$     write(*,*) 'cornerslocal is deallocated'
!!$  else
!!$     write(*,*) 'cornerslocal is not allocarted'
!!$  end if

  write(*,*) 'cornerslocal, first:'
  write(*,*) cornerslocal(1,1,1)
  write(*,*) 'cornerslocal, last:'
  write(*,*) cornerslocal(2,4,npanels+nfspanels)
  !if (ALLOCATED(cornerslocal)) DEALLOCATE(cornerslocal) ! gave error!

  WRITE(*,*) '1'
  IF (ALLOCATED(points)) DEALLOCATE(points)
  IF (ALLOCATED(ldx)) DEALLOCATE(ldx)
  WRITE(*,*) '12'
  IF (ALLOCATED(panels)) DEALLOCATE(panels)
  WRITE(*,*) '123'

  IF (ALLOCATED(center)) DEALLOCATE(center)
  WRITE(*,*) '1234'
  !IF (ALLOCATED(coordsys)) DEALLOCATE(coordsys)  !gave error!
  !WRITE(*,*) '12345'
  IF (ALLOCATED(PtsToPanelsConnectivity)) DEALLOCATE(PtsToPanelsConnectivity)


  !IF (ALLOCATED(area)) DEALLOCATE(area)  !gave error!
  WRITE(*,*) '123456'
  IF (ALLOCATED(vel)) DEALLOCATE(vel)
  !IF (ALLOCATED(velOld)) DEALLOCATE(velOld)
  IF (ALLOCATED(hessglobal)) DEALLOCATE(hessglobal)
  !IF (ALLOCATED(hessglobalOld)) DEALLOCATE(hessglobalOld)
  IF (ALLOCATED(dphi)) DEALLOCATE(dphi)
  IF (ALLOCATED(dphi_temp)) DEALLOCATE(dphi_temp)
  IF (ALLOCATED(ddphi)) DEALLOCATE(ddphi)

  WRITE(*,*) '1234567'
  IF (ALLOCATED(am)) DEALLOCATE(am)
  IF (ALLOCATED(b)) DEALLOCATE(b)
  IF (ALLOCATED(sigma)) DEALLOCATE(sigma)
  WRITE(*,*) '12345678'

  IF (ALLOCATED(vtotal)) DEALLOCATE(vtotal)
  WRITE(*,*) '123456789'
  IF (ALLOCATED(vn)) DEALLOCATE(vn)
  WRITE(*,*) '10'
  IF (ALLOCATED(vt1)) DEALLOCATE(vt1)
  WRITE(*,*) '101'
  IF (ALLOCATED(vt2)) DEALLOCATE(vt2)
  WRITE(*,*) '102'
  IF (ALLOCATED(cp)) DEALLOCATE(cp)
  WRITE(*,*) '103'
  IF (ALLOCATED(profilewave)) DEALLOCATE(profilewave)

  IF (ALLOCATED(vtest)) DEALLOCATE(vtest)
  IF (ALLOCATED(velt)) DEALLOCATE(velt)
  IF (ALLOCATED(zeta)) DEALLOCATE(zeta)
  IF (ALLOCATED(dzeta_dx)) DEALLOCATE(dzeta_dx)
  IF (ALLOCATED(dzeta_dy)) DEALLOCATE(dzeta_dy)

!  IF (ALLOCATED(celldata1)) DEALLOCATE(celldata1)
!  IF (ALLOCATED(celldata2)) DEALLOCATE(celldata2)
  IF (ALLOCATED(HullPtsUpdate)) DEALLOCATE(HullPtsUpdate)

END PROGRAM p3

