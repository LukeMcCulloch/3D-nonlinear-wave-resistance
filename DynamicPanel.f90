module DynamicPanel

  USE precise, ONLY : defaultp
  IMPLICIT NONE    
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp


CONTAINS


  SUBROUTINE computePtToPanelConnectivity(   npanels    , nfspanels, npoints,        &
                                             nfsl       , nfst     ,                 &
                                             panels     , points   ,                 &
                                             PtsToPanelsConnectivity, HullPtsUpdate, &
                                             ntransompanels, ntransoml, ntransomt,   &
                                             nhpl,TransomPtsUpdate)

    !! Poor man's hash 
    !! - array storage 
    !!   of point to panel
    !!   connectivity:

    integer :: i, j, k,l,m, row, col
    integer :: npanels, nfspanels, npoints
    integer :: ntransompanels, ntransoml, ntransomt

    integer :: nfsl                               !! number of free surface panels longitudinal
    integer :: nfst                               !! number of free surface panels transverse

    integer :: breakflag                          !! We found a match - get out of the loops
    integer :: breakflag2                         !! flag (1 of 2 ) to catch the final point - since l-hull loop is over nhl
    integer :: breakflag3                         !! flag (2 of 2) to catch the final point - since l-hull loop is over nhl
    INTEGER, DIMENSION(4)                         :: temp_points, temp_points2
    INTEGER, ALLOCATABLE, DIMENSION(:,:)          :: panels
    INTEGER, ALLOCATABLE, DIMENSION(:,:)          :: PtsToPanelsConnectivity
    INTEGER, ALLOCATABLE, DIMENSION(:,:)          :: HullPtsUpdate  ! store hull-fs panel connection
    INTEGER, ALLOCATABLE, DIMENSION(:,:)          :: TransomPtsUpdate

    integer :: nhpl                               !! the number of hull points in the longitudinal direction -make dynamic!
    integer :: hullStart                          !! starting indices for point indexing...
    integer :: hullpnl, fspnl                     !! dummies to hold the actual indexes of the hull and fs panels

    integer :: transomStart                       !! starting indices for point indexing...
    integer :: transompnl                         !! dummy to hold actual index of the transom panel

    REAL(WP), ALLOCATABLE, DIMENSION(:,:)         :: points
    real(wp), DIMENSION(3)                        :: hullpt, fspt, transompt
    real(wp)                                      :: distance, criteria

    write(*,*)  'computing point to panel connectivity...'

    ALLOCATE(PtsToPanelsConnectivity(5,npoints) ) !the 5th element is the index counter on the panel.
    
    PtsToPanelsConnectivity(1:4,:) = -1  ! flag all array elements as empty unless panels fill them
    PtsToPanelsConnectivity(5,:)   =  1  ! starting index of each panel (max is 4).
    
    do i=1,npanels+nfspanels+ntransompanels !loop all panels

       temp_points(:) = panels(:,i) ! pick off the point indexs of all points on the current panel
       do j=1,4
          PtsToPanelsConnectivity( PtsToPanelsConnectivity(5,temp_points(j))  ,  temp_points(j)  ) = i !add panel i cnctvty for pt j at cnctvty index carried in 5
          PtsToPanelsConnectivity( 5 ,temp_points(j)) = PtsToPanelsConnectivity(5,temp_points(j) ) + 1 !updte index
       end do
    end do
















    !! WARNING,
    !! BRUTE FORCE MERGE AHEAD:
    !! do you really want to do this?1?1?10?
    !   Merge 1-way, panel to point connectivity  - secretly!  
    !   Don't tell the free surface point to panel connectivity
    !   this keeps the update averaging of the free surface easy

    write(*,*) 'Computing Brute Force masked merge of hull-free surface intersection panels...'

    criteria = 0.0000000001
    breakflag=0

    !nhpl      = (nfsl/4)   !number of hull panels on the top of the hull            -static-> make dynamic:done
    hullStart = npanels - nhpl !+1 !upper hull panel start


    !ALLOCATE( HullPtsUpdate(2,2,nhpl)) !store correspondance between hull points and fs points.

    !start the fix HERE!  pick off points and panels - reset flags if needed
    ALLOCATE( HullPtsUpdate(4,nhpl+1)) !change - for each hull panel, 
    !store its local points on the fs, and 
    !for each of those points, 
    !store the panel and point on that panel
    !with which the point is coincident.
    ! currently - size =4
    ! set to get one pt-panelcombo

    !loop along the panels above the hull:
    breakflag2=0
    breakflag3=0
    do i=1,nhpl !npanels !nhpl
       hullpnl = hullstart + i !i  !index of the hull panel
       !print*, hullpnl
       !pause

       if (i==nhpl) then ! custom code for the end point
          breakflag2 = 1
       end if
       !grab each hull-fs panel's points:
       temp_points(:) = panels(:,hullpnl) !4 index values of 4 point's!
       !write(*,*) 'temp_hull_pt=',temp_points
       !loop over the points of each hull panel:
       do j=1,4
          hullpt = points(:,temp_points(j)) !index of hull point

          !compare with panel points on the 1st long. column of the free surface:
          do k=1,nfsl!nfspanels!nfsl !                                               -static-> make dynamic and shrink size
             fspnl = npanels+k  !+1 was done at night !index of the fs panel
             !print*, fspnl

             !grab each fs-hull panel's points index:
             temp_points2(:) = panels(:,fspnl)
             !write(*,*) 'temp_fs_pt=',temp_points2
             do l=1,4 

                !if (breakflag2==1) then
                ! ???????????????????  
                !end if
                fspt = points(:,temp_points2(l)) !x,y,z point from index of fs point
                
                !print*,hullpt
                !print*,fspt
                !write(*,*) 'hullpoint=',hullpt
                !write(*,*) 'fspoint=',fspt
                distance = sqrt(  ( hullpt(1)-fspt(1) )**2 + (hullpt(2)-fspt(2))**2 + (hullpt(3)-fspt(3))**2 )
                !write(*,*) 'distance, criteria:', distance, criteria
                breakflag =0
                if (abs(distance)<criteria) then
                   !! we found a match: 
                
                   write(*,*) 'We found a Match:', i,j,k,l
                   !write(*,*) 'distance, criteria:', distance, criteria
                   !points(:,panels(j,hullpnl))=points(:,panels(k,fspnl))
                   if (breakflag2 == 0 .and. breakflag3==0) then
                      HullPtsUpdate(1,i)=hullpnl !store the index of the hull panel          !maybe overwriting - ok here
                      HullPtsUpdate(2,i)=fspnl   !store the index of the fs panel
                      
                      HullPtsUpdate(3,i)=j       !store hull panel point
                      HullPtsUpdate(4,i)=l       !store fs panel point
                      breakflag = 1
                      
                   end if

                   if (breakflag2 ==1 .and. breakflag3==0) then
                      
                      HullPtsUpdate(1,i)=hullpnl !store the index of the hull panel          !maybe overwriting - ok here
                      HullPtsUpdate(2,i)=fspnl   !store the index of the fs panel
                      
                      HullPtsUpdate(3,i)=j       !store hull panel point
                      HullPtsUpdate(4,i)=l       !store fs panel point
                      breakflag=0
                      breakflag3=1
                   end if
                   

                   if (breakflag2 ==1 .and. breakflag3==1) then
                      
                      HullPtsUpdate(1,i+1)=hullpnl !store the index of the hull panel          !maybe overwriting - ok here
                      HullPtsUpdate(2,i+1)=fspnl   !store the index of the fs panel
                      
                      HullPtsUpdate(3,i+1)=j       !store hull panel point
                      HullPtsUpdate(4,i+1)=l       !store fs panel point
                      breakflag=1
                      breakflag3=1
                   end if                 

                   
                end if
                if (breakflag==1) then !exit fs pt loop
                   exit
                end if
             end do
             if (breakflag==1) then !exit fs pnl loop
                exit
             end if
          end do
          if (breakflag==1 .AND. breakflag2==0) then !exit hull pt loop
             exit
          end if
       end do
    end do !end hull panel loop



    write(*,*) 'hull merge complete'


    


















    write(*,*) 'Computing Brute Force masked merge of transom free surface - standard free surface intersection panels...'



    criteria = 0.0000000001
    breakflag=0
    transomStart = npanels+nfspanels+ntransompanels-ntransoml !upper hull panel start
    ALLOCATE( TransomPtsUpdate(4,ntransoml+1))
    do i=1,ntransoml,2
       transompnl = transomStart + i 
       !grab each transom-fs panel's points:
       temp_points(:) = panels(:,transompnl) !4 index values of 4 point's!
       !loop over the points of each transom fs panel:
       do j=1,4
          transompt = points(:,temp_points(j)) !transom point = Points(index of transom point)


          !compare with panel points on the 1st long. column of the free surface:
          do k=nfsl-ntransoml+1,nfsl
             fspnl = npanels+k !index of the free surface panel
             !grab each fs-hull panel's points index:
             temp_points2(:) = panels(:,fspnl)
             !loop over the points of each fs panel:
             do l=1,4 
                fspt = points(:,temp_points2(l)) !x,y,z point from index of fs point


                distance = sqrt( (transompt(1)-fspt(1))**2 + (transompt(2)-fspt(2))**2 + (transompt(3)-fspt(3))**2 )


                if (abs(distance)<criteria) then
                   !! we found a match: 
                
                   write(*,*) 'We found a Match:', i,j,k,l


                   PtsToPanelsConnectivity( PtsToPanelsConnectivity(5,temp_points(j))  ,  temp_points(j)  ) = fspnl !add panel i cnctvty for pt j at cnctvty index carried in 5
                   PtsToPanelsConnectivity(5,temp_points(j)) = PtsToPanelsConnectivity(5,temp_points(j)) + 1 !updte index   
                   !print*, 'adding panel ', fspnl, 'to point,', temp_points(j)

                   PtsToPanelsConnectivity( PtsToPanelsConnectivity(5,temp_points2(l))  ,  temp_points2(l)  ) = transompnl !add panel i cnctvty for pt j at cnctvty index carried in 5
                   PtsToPanelsConnectivity(5,temp_points2(l)) = PtsToPanelsConnectivity(5,temp_points2(l)) + 1 !updte index 
                   !print*, 'adding panel ', transompnl, 'to point,', temp_points2(j)       
                end if
             end do

          end do
       end do
    end do !end hull panel loop



    write(*,*) 'transom merge complete'




   END SUBROUTINE computePtToPanelConnectivity






END module DynamicPanel
