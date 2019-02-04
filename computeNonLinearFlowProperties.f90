module computeNonLinearFlowProperties
 

  USE precise, ONLY : defaultp
  USE A
  IMPLICIT NONE    
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp


CONTAINS

  SUBROUTINE computeVelocities(npanels,nfspanels,vtotal,vinf,sigma,vel,vn,vt1,vt2,coordsys, dphi_temp, dphi, &
                                       ntransompanels, ntransoml, ntransomt)
    !! get the total velocity on each panel
    !!
    !! get the normal velocity on each panel (should==0)
    !!
    !! get the tangential components on each panel:
    !!
    
    INTEGER :: i, j,transompt
    INTEGER :: npanels, nfspanels    
    integer :: ntransompanels, ntransoml, ntransomt

    
    REAL(WP), DIMENSION(3)                                    :: vinf
    REAL(WP), ALLOCATABLE, DIMENSION(:,:)                     :: vtotal
    REAL(WP), ALLOCATABLE, DIMENSION(:,:)                     :: dphi_temp
    REAL(WP), ALLOCATABLE, DIMENSION(:,:)                     :: dphi
    REAL(WP), ALLOCATABLE, DIMENSION(:)                       :: sigma
    REAL(WP), ALLOCATABLE, DIMENSION(:,:,:)                   :: vel
    REAL(WP), ALLOCATABLE, DIMENSION(:)                       :: vn
    REAL(WP), ALLOCATABLE, DIMENSION(:)                       :: vt1
    REAL(WP), ALLOCATABLE, DIMENSION(:)                       :: vt2
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)                   :: coordsys
  
    write(*,*) 'updating velocities...'
    
    !vtotal(1,:) = 0.
    !vtotal(2,:) = 0.
    !vtotal(3,:) = 0.
    


    !do j=1,ntransomt
    !   transompt =  npanels+nfspanels+(j-1)*ntransoml+1
    !   sigma(transompt)=0.
    !end do

    Do i=1,npanels+nfspanels+ntransompanels

       !Do j=1, npanels+nfspanels
       !   vtotal(:,i)= vtotal(:,i) + vel(:,i,j)*sigma(j)
       !end do


       !only the perturbation is included:
       dphi_temp(:,i) = 0.


       !total:
       !dphi_temp(:,i) = dphi(:,i)+vtotal(:,i) !old base flow !maybe still vinf... 
       !dphi_temp(:,i) = dphi(:,i)
       !dphi_temp(:,i) = -vinf(:)  !linear case





       Do j=1, npanels+nfspanels+ntransompanels
             dphi_temp(:,i) = dphi_temp(:,i)+sigma(j)*vel(:,i,j)
       End Do




       !vtotal = dphi_temp
       !write(*,*) vtotal(:,i)
       
       !! normal & tangential components:
       !vn(i)=Dot_Product(coordsys(:,3,i),vtotal(:,i))
       !vt1(i)=dot_product(coordsys(:,1,i),vtotal(:,i))
       !vt2(i)=dot_product(coordsys(:,2,i),vtotal(:,i))
       
       ! update the velocity field with the new flow:
       !dphi_temp(:,i) =  dphi(:,i) + vtotal(:,i)
       
    End Do
    
  END SUBROUTINE computeVelocities


  subroutine computeWaveHeight(npanels,nfspanels,&
                               vinf,zeta,vtotal, &
                               waveSOR, dphi_temp,dphi, &
                               center, nfsl, nfst,&
                               panels, points,&
                               PtsToPanelsConnectivity,&
                               HullPtsUpdate, &
                               ntransompanels, ntransoml, ntransomt, nhpl,&
                               TransomPtsUpdate)
    !! subroutine to compute the wave height zeta
    !!  from the dynamic boundary condition

    integer :: i,row, col, hullrow, firstColpt
    integer :: npanels, nfspanels
    integer :: nfsl                                 !! number of free surface panels in the long. dir.
    integer :: nfst                                 !! number of free surface panels in the transvs. dir.
    integer :: ntransompanels, ntransoml, ntransomt

    integer :: transomStart                         !! starting indices for point indexing...
    integer :: transompnl                           !! dummy to hold actual index of the transom panel



    !! Integers specific to Dynamic remesh:
    integer :: pt_index                             !!dummy
    integer :: dydummy1                             !!dummy
    integer :: j,k                                  !!loop dummies    
    integer :: nhpl                                 !!the number of hull panels in the longitudinal direction 
    integer :: fsHpSt                               !! free surface hull panel start index location

    INTEGER, DIMENSION(4)                           :: temp_points
    INTEGER, ALLOCATABLE, DIMENSION(:,:)            :: panels
    INTEGER, ALLOCATABLE, DIMENSION(:,:)            :: PtsToPanelsConnectivity
    INTEGER, ALLOCATABLE, DIMENSION(:,:)            :: HullPtsUpdate  ! store hull-fs panel connection
    INTEGER, ALLOCATABLE, DIMENSION(:,:)            :: TransomPtsUpdate

    integer :: hullStart,  hullpnl, fspnl,localHullPanelIndex, localHullPointIndex,localFsPanelIndex, localFsPointIndex, pi,fi

    INTEGER :: saveDx


    REAL(WP), parameter                             :: gravity = 9.81
    REAL(wp)                                        :: waveSOR
    REAL(WP)                                        :: slope, distance
    REAL(WP), DIMENSION(3)                          :: vinf
    real(wp), ALLOCATABLE, dimension(:)             :: zeta
    REAL(WP), ALLOCATABLE, DIMENSION(:,:)           :: vtotal  
    REAL(WP), DIMENSION(:,:)                        :: dphi_temp
    REAL(WP), DIMENSION(:,:)                        :: dphi
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)           :: center 

    !! Reals specific to Dynamic remesh:
    REAL(WP), ALLOCATABLE, DIMENSION(:,:)           :: points
    REAL(wp), DIMENSION(3)                          :: panelAverage





    write(*,*) 'updating wave height...'

    
    saveDx =  1+npanels+nfspanels+1!use the dx from the first panel instead of a continuous variation which was using "row"
    ! for panel height updates, below:  this variation caused normal vector jumps at fore and aft pts.
    !!
    !! Note:  (11-28-2013)  a large jump in panel length that occurs when perpendicular to hull end points 
    !!                      was making the unit vector normal point more in the direction of the flow
    !!                      and so it vn became large!
    !!  (I think)
    !!

    Do i=1,nfspanels+ntransompanels
       Row=i+npanels
       !zeta(Row) = (vinf(1)/9.81)*(vtotal(1,Row) + vinf(1))

       ! total update using non-linear equation
       !zeta(Row) =  (1.-waveSOR)*zeta(Row) + waveSOR*(.5/gravity)*(vinf(1)**2 &
       !                                                        -dphi_temp(1,Row)**2  &
       !                                                         -dphi_temp(2,Row)**2  &
       !                                                          -dphi_temp(3,Row)**2  )

       ! total update with lienarized eqution:
       !    update with perturation:
       zeta(Row) =  (1.-waveSOR)*zeta(Row) + waveSOR*((.5/gravity)*&
                                             ((vinf(1)**2) &
                                              -dphi(1,Row)**2 -2.*dphi(1,Row)*dphi_temp(1,Row) &
                                               -dphi(2,Row)**2 -2.*dphi(2,Row)*dphi_temp(2,Row) &
                                                -dphi(3,Row)**2 -2.*dphi(3,Row)*dphi_temp(3,Row) ))

       ! modify the centers of every "easy" panel - could improve efficiency by waiting?
       !if ( mod(i-1,nfsl) >0 ) then !we are behind the freeSurf start line so:
       !   center(3,row) = zeta(row)+.65*abs(center(1,row)-center(1,row-1)) !.65 panel length - change to dynamic.
       !end if
    End Do

  
    
    !! fixing the free surface to transom bottom at the transom-wrong
    do j=1,ntransomt 
       row = 1+npanels+nfspanels+(j-1)*ntransoml 
       hullrow = nhpl + (j-1)*nhpl

       slope =  (center(3,hullrow)-center(3,hullrow-1)) /(center(1,hullrow)-center(1,hullrow-1))
       distance = center(1,row)-center(1,hullrow)
       zeta(Row) = center(3,hullrow) + slope*distance

       !center(3,Row)= zeta(Row) + .80*abs(center(1,row)-center(1,hullrow))
       center(3,Row)= zeta(Row) + .80*abs(center(1,saveDX)-center(1,saveDx-1))
    end do



    !! fs panel centers:
    do j=1,nfst
       do i=2,nfsl
          row = i+npanels+(j-1)*nfsl
          !center(3,row) = zeta(row)+.80*abs(center(1,row)-center(1,row-1))
          center(3,row) = zeta(row)+.80*abs(center(1,saveDx)-center(1,saveDx-1))
       end do
    end do



    !! transom fs panel centers:
    do j=1,ntransomt 
       !firstColpt = 1+npanels+nfspanels+(j-1)*ntransoml
       !hullrow = nhpl + (j-1)*nhpl
       !center(3,firstColpt) = zeta(firstColpt)+.80*abs(center(1,firstColpt)-center(1,hullrow)) 
       do i=2,ntransoml
          row = i+npanels+nfspanels+(j-1)*ntransoml
          !center(3,row) = zeta(row)+.80*abs(center(1,row)-center(1,row-1))
          center(3,row) = zeta(row)+.80*abs(center(1,saveDx)-center(1,saveDx-1))
       end do
    end do  
    

    write(*,*) 'wave height computation done'
    write(*,*) ''



    write(*,*) 'Starting Free Surface Panel point Adaptation...'
    do  i=1+npanels,npanels+nfspanels+ntransompanels  !loop over all fs panels
       !print*, i
       !!grab points list from panel index:
       temp_points(:) = panels(:,i) 
       
       !!for each point in the list, points(pt_index) = avg(all panel center points in the list)
       do j=1,4 !each panel has 4 pts, gauranteed.
          pt_index       =    temp_points(j) !current pt index in the list
          DyDummy1       =    PtsToPanelsConnectivity(5,pt_index)-1 !number of panels connected
          panelAverage = 0.0
          do k=1,DyDummy1  !!take the average of those panel centers to be the new value for that point.
             panelAverage(:) = panelAverage(:) + center(:,PtsToPanelsConnectivity(k,pt_index))
          end do
          panelAverage(:) = panelAverage(:)/real(DyDummy1)
          points(3,pt_index) = panelAverage(3) !only update the  Z direction on the panel
       end do

    end do




    !!update corresponding hull-free surface interface points:
    !nhpl      = (nfsl/4)
    !hullStart = npanels - nhpl 


    do i=1,nhpl+1 !1 extra update for the last panel
       !print*, i
       localHullPanelIndex =  HullPtsUpdate(1,i)  ! we are at this hull panel
       localHullPointIndex =  HullPtsUpdate(3,i)  ! we are updating this hull pt
       
       localFsPanelIndex = HullPtsUpdate(2,i)     ! we will use this fs panel
       localFsPointIndex = HullPtsUpdate(4,i)     ! we will use this fs point

       pi = panels(localHullPointIndex,localHullPanelIndex)

       fi = panels(localFsPointIndex,localFsPanelIndex)


       ! We could do several things here:
       ! note we are only moving the panel corner pts of the hull-fs interface
       ! 

       ! 0.) do nothing -> results in bad velocities

       ! 1.) set the "hull" pt == fs pt: best for easy cp calculation... but...
       points(3,pi) = points(3,fi)  !! -> sometimes inverts the panel! 

       ! 2.) set the fs point = hull free surface pt:
       !points(3,fi) = points(3,pi)

       ! 3.) averaging: e.g.
       !points(3,pi) = (points(3,pi)+points(3,fi))/2.
       ! which could include:
       !points(3,fi) = points(3,pi)

       ! 4.) put some flag on points that would be repositioned below the panel bottom
       ! only modify those updates.  otherwise, use method (1.)

       !ie:
       ! (1.)
       !points(3,pi) = points(3,fi)
       !if points(3,pi) 

    end do

    write(*,*) 'hull to free surface panel adaptation done'




    write(*,*) 'transom free surface to free surface panel adaptation done'




  end subroutine computeWaveHeight


  subroutine setPanelLocation(npanels,nfspanels,&
                               vinf,zeta,vtotal, &
                               waveSOR, dphi_temp,dphi, &
                               center, nfsl, nfst,&
                               panels, points,&
                               PtsToPanelsConnectivity,&
                               HullPtsUpdate, &
                                       ntransompanels, ntransoml, ntransomt,nhpl,TransomPtsUpdate)
    !! subroutine to compute the wave height zeta
    !!  from the dynamic boundary condition

    integer :: i,row, col
    integer :: npanels, nfspanels
    integer :: nfsl                             !! number of free surface panels in the long. dir.
    integer :: nfst                             !! number of free surface panels in the transvs. dir.
    integer :: ntransompanels, ntransoml, ntransomt

    !! Integers specific to Dynamic remesh:
    integer :: pt_index                         !! dummy
    integer :: dydummy1                         !! dummy
    integer :: j,k                              !! loop dummies    
    integer :: nhpl                             !! the number of hull panels in the longitudinal direction
    integer :: fsHpSt                           !! free surface hull panel start index location

    INTEGER, DIMENSION(4)                       :: temp_points
    INTEGER, ALLOCATABLE, DIMENSION(:,:)        :: panels
    INTEGER, ALLOCATABLE, DIMENSION(:,:)        :: PtsToPanelsConnectivity
    INTEGER, ALLOCATABLE, DIMENSION(:,:)        :: HullPtsUpdate  ! store hull-fs panel connection
    INTEGER, ALLOCATABLE, DIMENSION(:,:)        :: TransomPtsUpdate

    integer :: hullStart,  hullpnl, fspnl,localHullPanelIndex, localHullPointIndex,localFsPanelIndex, localFsPointIndex, pi,fi

    integer :: transomStart                       !! starting indices for point indexing...
    integer :: transompnl                         !! dummy to hold actual index of the transom panel



    REAL(WP), parameter                         :: gravity = 9.81
    REAL(wp)                                    :: waveSOR
    REAL(WP), DIMENSION(3)                      :: vinf
    real(wp), ALLOCATABLE,dimension(:)          :: zeta
    REAL(WP), ALLOCATABLE,DIMENSION(:,:)        :: vtotal  
    REAL(WP), DIMENSION(:,:)                    :: dphi_temp
    REAL(WP), DIMENSION(:,:)                    :: dphi
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)       :: center 

    !! Reals specific to Dynamic remesh:
    REAL(WP), ALLOCATABLE, DIMENSION(:,:)       :: points
    REAL(wp), DIMENSION(3)                      :: panelAverage


  

       ! modify the centers of every "easy" panel - could improve efficiency by waiting?
       !if ( mod(i-1,nfsl) >0 ) then !we are behind the freeSurf start line so:
       !   center(3,row) = zeta(row)+.65*abs(center(1,row)-center(1,row-1)) !.65 panel length - change to dynamic.
       !end if
  
    !! fs panels:    
    do j=1,nfst
       do i=1,nfsl
          row = i+npanels+(j-1)*nfsl
          center(3,row) = zeta(row)
       end do
    end do

    !! transom fs panels:
    do j=1,ntransomt
       do i=1,ntransoml
          row = i+npanels+nfspanels+(j-1)*ntransoml
          center(3,row) = zeta(row)
       end do
    end do    




    write(*,*) 'Starting Free SurfacePanel Adaptation...'
    do  i=1+npanels,npanels+nfspanels+ntransompanels  !loop over all fs panels
       
       !!grab points list from panel index:
       temp_points(:) = panels(:,i) 
       
       !!for each point in the list, points(pt_index) = avg(all panel center points in the list)
       do j=1,4 !each panel has 4 pts, gauranteed.

          pt_index       =    temp_points(j) !current pt index in the list

          DyDummy1       =    PtsToPanelsConnectivity(5,pt_index)-1 !number of panels connected

          panelAverage = 0.0
          do k=1,DyDummy1  !!take the average of those panel centers to be the new value for that point.
             panelAverage(:) = panelAverage(:) + center(:,PtsToPanelsConnectivity(k,pt_index))
          end do
          panelAverage(:) = panelAverage(:)/real(DyDummy1)

          points(3,pt_index) = panelAverage(3) !only update the  Z direction on the panel
          
       end do

    end do




    !!update corresponding hull-free surface interface points:
    nhpl      = (nfsl/4)
    hullStart = npanels - nhpl 


    do i=1,nhpl+1 !1 extra update for the last panel
       localHullPanelIndex =  HullPtsUpdate(1,i)  ! we are at this hull panel
       localHullPointIndex =  HullPtsUpdate(3,i)  ! we are updating this hull pt
       
       localFsPanelIndex = HullPtsUpdate(2,i)     ! we will use this fs panel
       localFsPointIndex = HullPtsUpdate(4,i)     ! we will use this fs point

       pi = panels(localHullPointIndex,localHullPanelIndex)

       fi = panels(localFsPointIndex,localFsPanelIndex)


       ! We could do several things here:
       ! note we are only moving the panel corner pts of the hull-fs interface
       ! 

       ! 0.) do nothing -> good for end display - let the hull fs panels rise from the water.


    end do







    write(*,*) 'panel adaptation done'




  end subroutine setPanelLocation



  subroutine computeDZeta(npanels,nfspanels, nfsl, nfst, zeta, dzeta_dx, dzeta_dy, ldx, center, &
                                       ntransompanels, ntransoml, ntransomt, nhpl)
    !! upwinding scheme to compute x derivatives
    !! central scheme to compute y derivatives

    !! note proper way to pass previously allocated stuff!
    integer :: i,j, hullpt
    integer :: nhpl                                   !! number of hull pts in the longitudinal dir.
    integer :: npanels, nfspanels                     !! num. panels on hull, num on free surface
    integer :: nfsl                                   !! number of free surface panels in the longi. dir.
    integer :: nfst                                   !! number of free surface panels in the transvs. dir.
    integer :: colSum, currentPos, cPm1, cPm2, cPp1   !! index positioning dummy vbls
    integer :: ntransompanels, ntransoml, ntransomt

    REAL(WP), parameter                         :: gravity = 9.81
    REAL(wp), ALLOCATABLE, dimension(:)         :: zeta
    REAL(wp), ALLOCATABLE, dimension(:)         :: dzeta_dx
    REAL(wp), ALLOCATABLE, dimension(:)         :: dzeta_dy
    REAL(wp), ALLOCATABLE, DIMENSION(:)         :: ldx  !size (nfsl) dx step array
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)       :: center !could have used instead


    write(*,*) 'computing space derivatives of wave height...'

    !! 3 pt stencil, 1st order upwind difference  Order(?)accurate



    !write(*,*) 'computing d zeta / dx'
    dzeta_dx(:)=0.

    !!----------------------------------------------------------------------------------
    !! x - derivatives
    !!----------------------------------------------------------------------------------
    !! 1.) non-transom free surface panels : compute dzeta/dx on the 
    do j=1,nfst !! loop over all trans. panels
       
       
       colSum = (j-1)*nfsl+npanels !sum all elements in columns already traversed.

       !! first two panels in the long. dir.
       dzeta_dx(1+colSum)    = 0.
       dzeta_dx(2+colSum)    = (1.5*zeta(colSum+2)&
                                -2.*zeta(colSum+1)&
                                +0.5*zeta(colSum+1))&
                                / (center(1,colSum+2)-center(1,colSum+1))   !ldx(2)

       !! loop over the remaining long. panels.
       do i=3,nfsl        
          currentPos = (i+colSum)
          cPm1 = (currentPos-1)
          cPm2 = (currentPos-2)

          dzeta_dx(currentPos) = (1.5*zeta(currentPos)&
                                 -2.*zeta(cPm1)&
                                 +0.5*zeta(cPm2)) &
                                 / (center(1,currentPos)-center(1,cPm1) )
       end do
    end do


    !!----------------------------------------------------------------------------------
    !! 2.) transom free surface panels : compute dzeta/dx 
    do j=1,ntransomt
       colSum = (j-1)*ntransoml+npanels+nfspanels !sum all elements in columns already traversed.
       hullpt = nhpl + (j-1)*nhpl

       
       !! first two panels in the long. dir.: Raven, Taylor's Expansion page 138 (141 in pdf):
       !! Very first point is as described 
       dzeta_dx(1+colSum)    =  (2.*(zeta(colSum+1)-center(3,hullpt) ) /(center(1,colSum+2)-center(1,hullpt))) - &
                                                     ((center(3,hullpt)-center(3,hullpt-1)) / (center(1,hullpt)-center(1,hullpt-1)))

       !! second point updated to Raven spec, page 104 at http://www.nap.edu/openbook.php?record_id=9223&page=104
       !! Session 3 - wavy free surface flow: Panel Methods.
       !dzeta_dx(2+colSum)    = (2.*(zeta(colSum+2)-zeta(colSum+1))/(center(1,colSum+2)-center(1,colSum+1))) - &
       !                                              ((center(3,hullpt)-center(3,hullpt-1)) / (center(1,hullpt)-center(1,hullpt-1)))
       dzeta_dx(2+colSum)    = (1.5*zeta(colSum+2)&
                                -2.*zeta(colSum+1)&
                                +0.5*zeta(colSum+1))&
                                / (center(1,colSum+2)-center(1,colSum+1))  

       do i=3,ntransoml

          currentPos = (i+colSum)
          cPm1 = (currentPos-1)
          cPm2 = (currentPos-2)
          dzeta_dx(currentPos) = (1.5*zeta(currentPos)&
                                 -2.*zeta(cPm1)&
                                 +0.5*zeta(cPm2)) &
                                 / (center(1,currentPos)-center(1,cPm1) )
       end do
    end do




    !!----------------------------------------------------------------------------------
    !! y - derivatives
    !!----------------------------------------------------------------------------------
    !! dy central difference
    !! what does this enforce about the centerplane?
    !write(*,*) 'computing d zeta / dy'
    dzeta_dy(:)=0.
    !write(*,*) 'nfsl, nfst ', nfsl, nfst




    !!----------------------------------------------------------------------------------
    !! regular free surface
    j=1 !first column
    colSum = (j-1)*nfsl + npanels

    !  first column, standard
    do i=1,(nfsl-ntransoml) ! move x direction column
       currentPos = colSum+i 
       cPm1 = currentPos
       cPp1 = currentPos + nfsl
       dzeta_dy(currentPos) =  ( zeta(cPp1) - zeta(cPm1) )&
            /(2.*( center(2,cPp1) - center(2,cPm1) ) )
    end do

    ! first column, standard free surface using transom free surface
    do i=nfsl-ntransoml+1,nfsl
       currentPos = colSum+i 
       cPm1 = npanels+nfspanels+ntransompanels - ntransoml + i
       cPp1 = currentPos + nfsl
       dzeta_dy(currentPos) =  ( zeta(cPp1) - zeta(cPm1) )&
            /(2.*( center(2,cPp1) - center(2,cPm1) ) )
    end do    



    !!middle columns, standard free surface
    do j=2,nfst-1 ! moving transverse

       colSum = (j-1)*nfsl + npanels

       do i=1,nfsl ! move x direction column
          currentPos = colSum+i 
          cPm1 = currentPos - nfsl
          cPp1 = currentPos + nfsl
          dzeta_dy(currentPos) =  ( zeta(cPp1) - zeta(cPm1) )&
               /( center(2,cPp1) - center(2,cPm1)  )
       end do
    end do


    j=nfst  !last column, standard free surface
    colSum = (j-1)*nfsl + npanels
    do i=1,nfsl ! move x direction column
       currentPos = colSum+i 
       cPm1 = currentPos - nfsl
       cPp1 = currentPos
       dzeta_dy(currentPos) =  ( zeta(cPp1) - zeta(cPm1) )&
            /(2.*( center(2,cPp1) - center(2,cPm1)  ))
    end do
    


    !!----------------------------------------------------------------------------------
    ! Transom free surface, first column
    j=1
    colSum = (j-1)*ntransoml + npanels+nfspanels
    do i=1,ntransoml
       currentPos = i+colSum
       cPm1 = currentPos
       cPp1 = currentPos + ntransoml
       dzeta_dy(currentPos) =  ( zeta(cPp1) - zeta(cPm1) )&
            /(2.*( center(2,cPp1) - center(2,cPm1) ) )
    end do
    !Transom, middle columns
    do j=2,ntransomt-1
       colSum = (j-1)*ntransoml + npanels+nfspanels
       do i=1,ntransoml
          currentPos=colSum+i
          cPm1 = currentPos - ntransoml
          cPp1 = currentPos + ntransoml
          dzeta_dy(currentPos) =  ( zeta(cPp1) - zeta(cPm1) )&
               /( center(2,cPp1) - center(2,cPm1)  )
       end do
    end do          
    ! Transom, last column
    j=ntransomt
    colSum = (j-1)*ntransoml + npanels+nfspanels
    do i=1,ntransoml
       currentPos=colSum+i
       cPm1 = currentPos - ntransoml
       cPp1 = currentPos
       dzeta_dy(currentPos) =  ( zeta(cPp1) - zeta(cPm1) )&
            /(2.*( center(2,cPp1) - center(2,cPm1)  ))       
    end do

  end subroutine computeDZeta
  




  SUBROUTINE updateSurfaceVelocities( npanels, nfspanels, fieldpoint,fieldpointp,center,deltax,&
       vel,dphi, coordsys, area, cornerslocal,vp,vinf, hp, hessglobal, ddphi, &
       zeta, dzeta_dx, dzeta_dy, sigma, vtotal, &
                                       ntransompanels, ntransoml, ntransomt)
    !! get the total velocity on each panel
    !!
    !! get the normal velocity on each panel (should==0)
    !!
    !! get the tangential components on each panel:
    !!
    INTEGER                                                      :: i,j,k
    INTEGER                                                      :: Row, Col
    INTEGER                                                      :: npanels, nfspanels
    integer :: ntransompanels, ntransoml, ntransomt

    REAL(WP), parameter                 ::    gravity = 9.81
    REAL(WP), DIMENSION(3,3), parameter ::    hpsign  = RESHAPE( (/ &
                                                                    1., -1.,  1., &
                                                                   -1.,  1., -1., &
                                                                    1., -1.,  1.  /), (/3,3/)) 
    REAL(WP), DIMENSION(:,:)                                 :: vtotal
    REAL(wp)                                                 :: deltax
    REAL(WP), DIMENSION(3)                                   :: vinf
    REAL(WP), DIMENSION(3)                                   :: fieldpoint
    REAL(WP), DIMENSION(3)                                   :: fieldpointp
    REAL(WP), DIMENSION(3)                                   :: vp
    REAL(wp), ALLOCATABLE,DIMENSION(:,:)                     :: center
    REAL(WP), ALLOCATABLE,DIMENSION(:,:,:)                   :: cornerslocal 
    REAL(wp), ALLOCATABLE,DIMENSION(:)                       :: area
    REAL(wp), ALLOCATABLE,DIMENSION(:,:,:)                   :: coordsys

    REAL(WP), ALLOCATABLE,DIMENSION(:,:,:)                   :: vel
    REAL(WP), ALLOCATABLE,DIMENSION(:,:)                     :: dphi

    REAL(WP), ALLOCATABLE,DIMENSION(:,:,:,:)                 :: hessglobal
    REAL(WP), ALLOCATABLE,DIMENSION(:,:,:)                   :: ddphi
    REAL(WP), DIMENSION(3,3)                                 :: hp

    REAL(WP), ALLOCATABLE,DIMENSION(:)                       :: zeta
    REAL(WP), ALLOCATABLE,DIMENSION(:)                       :: dzeta_dx
    REAL(WP), ALLOCATABLE,DIMENSION(:)                       :: dzeta_dy

    REAL(WP), ALLOCATABLE, DIMENSION(:)                      :: sigma

    write(*,*) ''
    write(*,*) 'computing velocities at the new free surface......'




    !! Update the free surface collocation PHI
    Do i=1,nfspanels+ntransompanels

       !Get the field points on the first order free surface, include x offset 
       
       Row = npanels + i
       
       fieldpoint(:)=center(:,Row)
       
       !! No mirrors, but need offsets in i and k
       fieldpoint(1)=fieldpoint(1)+deltax  !! shift the collocation (field) points deltax
       fieldpoint(3)=zeta(Row) !0.         !! evaluate the velocities at the surface
       
       fieldpointp(:) = fieldpoint(:)
       fieldpointp(2) = -fieldpoint(2)
       
       !! Calculate the contribution of body panels to the free surface condition
       !! ***********************************************************************
       Do j=1,npanels
          
          Col=j
          
          !  Call what kind of influence function?  Use the nonLinear free surface condition
          vel(:,Row,Col)=hsinfluence(fieldpoint(:),&
                                     center(:,Col),&
                                     area(Col),&
                                     coordsys(:,:,Col),&
                                     cornerslocal(:,:,Col))

          vp(:)         =hsinfluence(fieldpointp(:),&
                                     center(:,Col),&
                                     area(Col),&
                                     coordsys(:,:,Col),&
                                     cornerslocal(:,:,Col))
          
          !Velocities for the half ship (flow is symetric)
          vel(1,Row,Col) = vel(1,Row,Col)+vp(1)
          vel(2,Row,Col) = vel(2,Row,Col)-vp(2)
          vel(3,Row,Col) = vel(3,Row,Col)+vp(3)
          
          
          hessglobal(:,:,Row,Col)=phixxinfluence(fieldpoint(:),&
                                                 center(:,Col),&
                                                 area(Col),&
                                                 coordsys(:,:,Col),&
                                                 cornerslocal(:,:,Col))

          hp(:,:)                =phixxinfluence(fieldpointp(:),&
                                                 center(:,Col),&
                                                 area(Col),&
                                                 coordsys(:,:,Col),&
                                                 cornerslocal(:,:,Col)) !mirrored fieldpoint
          
          hessglobal(:,:,Row,Col) = hessglobal(:,:,Row,Col)+hp(:,:)*hpsign(:,:)
          
           
       END DO
       
       !! Calculate the contribution of free surface panels to the free surface condition
       !!********************************************************************************
       Do j=1,nfspanels+ntransompanels
          
          Col=j+npanels
          
          !   Enforce PHIz 
          vel(:,Row,Col)=hsinfluence(fieldpoint(:),&
                                     center(:,Col),&
                                     area(Col),&
                                     coordsys(:,:,Col),&
                                     cornerslocal(:,:,Col))

          vp(:) =         hsinfluence(fieldpointp(:),&
                                     center(:,Col),&
                                     area(Col),&
                                     coordsys(:,:,Col),&
                                     cornerslocal(:,:,Col))
          
          !Velocities for the half ship (flow is symetric)
          vel(1,Row,Col) = vel(1,Row,Col)+vp(1)
          vel(2,Row,Col) = vel(2,Row,Col)-vp(2)
          vel(3,Row,Col) = vel(3,Row,Col)+vp(3)
          
          ! Enforce PHIxx squared
          hessglobal(:,:,Row,Col)=phixxinfluence(fieldpoint(:),&
                                                 center(:,Col),&
                                                 area(Col),&
                                                 coordsys(:,:,Col),&
                                                 cornerslocal(:,:,Col))

          hp(:,:) =               phixxinfluence(fieldpointp(:),&
                                                 center(:,Col),&
                                                 area(Col),&
                                                 coordsys(:,:,Col),&
                                                 cornerslocal(:,:,Col))
          
          hessglobal(:,:,Row,Col)=hessglobal(:,:,Row,Col)+hp(:,:)*hpsign(:,:)
          

       END DO
       
    END DO














    !! use new influences with known sigma to update dphi and ddphi

    !! edit, or - the new influences are updates to vinf

    do i=1,nfspanels+npanels+ntransompanels
       row = i
       !dphi(:,row)    = 0.
       !dphi(1,row)    = vinf(1)
       !ddphi(:,:,row) = 0.
       do j=1,npanels+nfspanels+ntransompanels !influence of the source at center j on fieldpoint on 
          col=j

          dphi(:,row)    = dphi(:,row)    + vel(:,row,col)*sigma(col)
          ddphi(:,:,row) =  ddphi(:,:,row) + hessglobal(:,:,Row,Col)*sigma(col)
       end do

    end do




    
    
  END SUBROUTINE updateSurfaceVelocities





  SUBROUTINE computeFreeSurfaceConvergence(npanels,nfspanels,dphi,zeta, dzeta_dx, dzeta_dy, vinf, nfsl, &
       errorK, errorD, locK, locD ,lock2,errork2, &
                                       ntransompanels, ntransoml, ntransomt)
    INTEGER :: i, j, row, col
    INTEGER :: npanels, nfspanels, nfsl
    INTEGER :: locK,locD,lock2
    integer :: ntransompanels, ntransoml, ntransomt

    REAL(WP), parameter                                                  :: gravity = 9.81
    REAL(WP), DIMENSION(3)                                               :: vinf
    REAL(WP)                                                             :: errorK, tke, errork2, tke2
    REAL(WP)                                                             :: errorD, tde
    REAL(WP), ALLOCATABLE, DIMENSION(:,:)                                :: dphi
    REAL(WP), ALLOCATABLE,DIMENSION(:)                                   :: zeta
    REAL(WP), ALLOCATABLE,DIMENSION(:)                                   :: dzeta_dx
    REAL(WP), ALLOCATABLE,DIMENSION(:)                                   :: dzeta_dy

    print*, 'Checking convergence criteria'
    errorK=0.
    errork2=0.
    
    errorD=0.

    locD=-1
    locK=-1
    lock2=-1
    print*, 'errorK, errorD'  
    do i=1,nfspanels+ntransompanels
       row = i+npanels

       tke =  abs( -dphi(3,row) - dzeta_dx(row)*dphi(1,row) - dzeta_dy(row)*dphi(2,row) )
       tke2 =  abs( dphi(3,row) - dzeta_dx(row)*dphi(1,row) - dzeta_dy(row)*dphi(2,row) )
       tde =   ( 0.5*(vinf(1)**2 - dphi(1,row)**2 - dphi(2,row)**2 - dphi(3,row)**2)) - gravity*zeta(row) 

       if (tke>errorK) then
          locK=i
       end if

       if (tke2>errorK2) then
          locK2=i
       end if

       if (tde>errorD) then
          locD=i
       end if   

       errorK=max(errorK, tke  )
       errorK2=max(errorK2, tke2  )
       errorD=max(errorD,tde )


    end do
  END SUBROUTINE computeFreeSurfaceConvergence





  SUBROUTINE computePressures(npanels,nfspanels,vtotal,vinf,cp,ntransompanels)
    !! Get the pressure on each panel 
    !! from the steady Bernoulli Eq:

    !! NOT UPDATED FOR NONLINEAR PANEL YET

    INTEGER :: i, j
    INTEGER :: npanels, nfspanels
    integer :: ntransompanels

    REAL(WP), DIMENSION(3)                  :: vinf
    REAL(WP), ALLOCATABLE,DIMENSION(:,:)    :: vtotal    
    REAL(WP), ALLOCATABLE,DIMENSION(:)      :: cp

    print*, 'computing pressures at the still water line... (must change)'
    Do i=1,npanels       
       cp(i)=1.-((sum(vtotal(:,i)**2))/(sum(vinf**2)))
    End do

  END SUBROUTINE computePressures



  SUBROUTINE checkSourceSink(npanels, nfspanels, sigma, area, sourcesink,ntransompanels)
    !! check that the sum of sources and sinks on the panels is equal to 0.
    integer :: i
    integer :: npanels, nfspanels
    integer :: ntransompanels

    REAL(WP)                                                  :: sourcesink
    REAL(WP),  ALLOCATABLE,DIMENSION(:)                       :: sigma
    REAL(WP),  ALLOCATABLE,DIMENSION(:)                       :: area
 

    
    sourcesink=0.
    Do i=1,npanels+nfspanels+ntransompanels
       sourcesink = sourcesink+sigma(i)*area(i)
    end do
    
    write(6,*) 'sourcesink sum (should be near 0.)',sourcesink
  END SUBROUTINE checkSourceSink







  subroutine computeHullForces(npanels,nfspanels,center,force,area,vinf,cp,coordsys)
    !! Coompute the forces on the hull
    integer :: i
    integer :: npanels,nfspanels

    REAL(wp), ALLOCATABLE,DIMENSION(:,:)    :: center
    REAL(WP), DIMENSION(3)                  :: force
    REAL(WP), ALLOCATABLE,DIMENSION(:)      :: area
    REAL(WP), DIMENSION(3)                  :: vinf
    REAL(WP), ALLOCATABLE,DIMENSION(:)      :: cp
    REAL(wp), ALLOCATABLE,DIMENSION(:,:,:)  :: coordsys


    print*, 'Computing LINEAR hull forces... (must change)'
    force(:)=(/0.,0.,0./)
    
    Do i=1,npanels  ! Hull summation only
       
       if(center(3,i)<0.)Then
          
          force(:)=force(:)+ area(i)*( .5*1025.*sum(vinf**2.)*cp(i)-1025.*9.81*center(3,i) )*coordsys(:,3,i)
       Else 
          
       End if
       
    end do
  end subroutine computeHullForces

  

  subroutine computeWaveResistance(npanels,nfspanels,center,S,area,vinf,force,cw)
    !! Subroutine to compute 
    !!   the wave resistance coefficient on the
    !!   panel hull

    integer :: i
    integer :: npanels,nfspanels

    real(wp)                                    :: S
    real(wp)                                    :: cw
    real(wp), DIMENSION(3)                      :: force
    real(wp), ALLOCATABLE,DIMENSION(:,:)        :: center    
    real(wp), ALLOCATABLE,DIMENSION(:)          :: area
    REAL(WP), DIMENSION(3)                      :: vinf

    print*, 'Computing Wave Resistance from Linear Form...(must change)'
    S=0.
    Do i=1,npanels
       if(center(3,i)<0.)Then
          S=S+area(i)
       Else
       End if
    End do
    write(*,*) 'force1',force(1)
    cw = -force(1)/(.5*1000.*(vinf(1)**2.)*S)
    write(6,*) 'cw', cw
    write(6,*) 'force', force
  end subroutine computeWaveResistance


  subroutine computeWaveProfile(npanels,nfspanels,profilewave,center,deltax,zeta, nfsl)
    !! subroutine to compute the wave Profile
    integer :: i,row
    integer :: npanels,nfspanels,  nfsl

    real(wp) :: deltax
    REAL(WP), ALLOCATABLE,DIMENSION(:)      :: zeta
    REAL(WP), DIMENSION(:,:)                :: profilewave
    real(wp), ALLOCATABLE,DIMENSION(:,:)    :: center    
    
    do i=npanels,npanels+nfsl
       row=(i-npanels)
       profilewave(:,row) = (/center(1,i)+deltax, zeta(i)/)
    end do
  end subroutine computeWaveProfile
  
end module computeNonLinearFlowProperties



