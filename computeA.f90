MODULE computeA  !! note i think this module will cause cornerslocal to go out of scope if it isn't used again.

  !! Compute the A matrix
  !!---------------------------------------------------------------------
  !!---------------------------------------------------------------------
  !!
  !! procedure to compute the A matrix of v dot a 
  !!
  !! Input:  the locations of panel centers pi and qj 
  !!         the normal vector of each panel center
  !!         the area of each panel
  !!         
  !!
  !!
  !! Output:  the velocity at panel i induced by panel j
  !!          the A matrix
  !!  
  !!
  !!
  !! Note we will use the following  from fifthpanel Output:
  !!     coordsys: a 3x3 matrix containing the two tangent vectors
  !!               and the normal vector
  !!                  t1x   t2x   nx
  !!                  t1y   t2y   ny
  !!                  t1z   t2z   nz
  !!               i.e. coordsys(:,3) is the normal vector
  !!     area: the area of the panel
  !!     center: the 3D coordinates of the panel center
  !!

  USE precise, ONLY : defaultp
  USE A

  IMPLICIT NONE    
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp



!!$  interface
!!$     subroutine BodyBC(a)
!!$       real(WP),allocatable :: cornerslocal(:,:,:)
!!$     end subroutine BodyBC
!!$  end interface
!!$  
!!$  interface
!!$     subroutine freeSurfaceBC(a)
!!$       real(WP),allocatable :: cornerslocal(:,:,:)
!!$     end subroutine freeSurfaceBC
!!$  end interface


CONTAINS


  
  SUBROUTINE computeBodyBC( npanels, nfspanels, fieldpoint,fieldpointp,center,&
       vel,coordsys,area,cornerslocal,vp,am, &
                                       ntransompanels, ntransoml, ntransomt)
    !
    ! Function to compute the A matrix Body Boundary Condition
    !
    INTEGER                                                      :: i,j,k
    INTEGER                                                      :: Row, Col
    INTEGER                                                      :: npanels, nfspanels
    integer :: ntransompanels, ntransoml, ntransomt
    
    REAL(WP), DIMENSION(3)                                       :: fieldpoint
    REAL(WP), DIMENSION(3)                                       :: fieldpointp
    REAL(WP), DIMENSION(3)                                       :: vp
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)                        :: center
    REAL(WP), ALLOCATABLE, DIMENSION(:,:,:)                      :: cornerslocal 
    REAL(wp), ALLOCATABLE, DIMENSION(:)                          :: area
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)                      :: coordsys
    REAL(WP), ALLOCATABLE, DIMENSION(:,:,:)                      :: vel
    REAL(WP), ALLOCATABLE, DIMENSION(:,:)                        :: am

    write(*,*) ''
    write(*,*) 'computing body boundary conditions..........'

    !! First get the image system of mirrored points  (enforce the centerplane!)
    !! in lieu of mirroring panels
    DO i=1,npanels
       
       Row=i
     
       fieldpoint(:)=center(:,Row)
       
       !! First mirror the fieldpoints about the Z-X plane
       fieldpointp(:)=fieldpoint(:)
       fieldpointp(2)=-fieldpoint(2)
       
       
       !! Compute the Influence of the Body Panels on hull fieldpoints
       !!****************************************************************************
       Do j=1,npanels
          
          Col = j
          !
          !          !First get the singular points
          IF(Row==Col)THEN
             vel(:,Row,Col)=-0.5*coordsys(:,3,Row)
             !
             !          !Then get the nonsingular points
          ELSE
             vel(:,Row,Col)=hsinfluence(fieldpoint(:),&
                                        center(:,Col),&
                                        area(Col),&
                                        coordsys(:,:,Col),&
                                        cornerslocal(:,:,Col))
          END IF
          
          !Now get the influence of the panel 
          !on the mirrored points  
          ! - enforce no flow through the centerplane of the ship
          vp(:)=hsinfluence(fieldpointp(:),&
                            center(:,Col),&
                            area(Col),&
                            coordsys(:,:,Col),&
                            cornerslocal(:,:,Col))
          !
          
          !Get Influence of panel sj 
          !and its mirrors @ point pi
          
          !Use no mirrors if modeling the whole geometry.
          
          !Velocities for the half ship (flow is symetric)
          vel(1,Row,Col) = vel(1,Row,Col)+vp(1)
          vel(2,Row,Col) = vel(2,Row,Col)-vp(2)
          vel(3,Row,Col) = vel(3,Row,Col)+vp(3)
          
          
          
          !Finally, make the a(i,j) matrix components satisfying no flow through the body surface  
          !(once the strengths are right)
          am(Row,Col)= DOT_PRODUCT(coordsys(:,3,Row),vel(:,Row,Col))
          
       END DO
       
       !! Compute the Influence of the Free Surface Panels on hull fieldpoints
       !! ***********************************************************************************
       Do j=1,nfspanels+ntransompanels
          
          Col = j+npanels
          !
          ! There are no singular points.
          ! We put the sources above the free surface!
          !  -will this help matrix inversion though?
          !IF(Row==Col)THEN !this never happens
          !   vel(:,Row,Col)=-0.5*coordsys(:,3,Row)
          !   !
          !   !          
          !ELSE ! Get the nonsingular points:
          vel(:,Row,Col)=hsinfluence(fieldpoint(:),&
                                     center(:,Col),&
                                     area(Col),&
                                     coordsys(:,:,Col),&
                                     cornerslocal(:,:,Col))
          !END IF
          
          !Compute influence of the panel on the mirrored points:
          vp(:) =           hsinfluence(fieldpointp(:),&
                                        center(:,Col),&
                                        area(Col),&
                                        coordsys(:,:,Col),&
                                        cornerslocal(:,:,Col))
          
          
          

          !Get Influence of panel sj 
          !and its mirror @ point pi
          
          
          !Velocities for the half ship (flow is symetric)
          vel(1,Row,Col) = vel(1,Row,Col)+vp(1)
          vel(2,Row,Col) = vel(2,Row,Col)-vp(2)
          vel(3,Row,Col) = vel(3,Row,Col)+vp(3)
          
          
          
          !Finally, a(i,j) matrix components satisfying no flow through the body surface (once the strengths are right)
          am(Row,Col)= DOT_PRODUCT(coordsys(:,3,Row),vel(:,Row,Col))
          
       END DO
    END DO
  END SUBROUTINE ComputeBodyBC




  SUBROUTINE computeLinearizedFreeSurfaceBC( npanels, nfspanels, fieldpoint,fieldpointp,center,deltax,&
       vel,coordsys,area,cornerslocal,vp,vinf, hp, hessglobal, am)
    !
    ! Function to compute the A matrix free surface Boundary Condition
    !
    INTEGER                                                      :: i,j,k
    INTEGER                                                      :: Row, Col
    INTEGER                                                      :: npanels, nfspanels
    INTEGER                                                      :: hullpt  !! used to find the free surface at the point of transom departure    

    real(wp)                                                     :: deltax
    REAL(WP), DIMENSION(3)                                       :: vinf
    REAL(WP), DIMENSION(3)                                       :: fieldpoint
    REAL(WP), DIMENSION(3)                                       :: fieldpointp
    REAL(WP), DIMENSION(3)                                       :: vp
    REAL(wp), DIMENSION(3,npanels+nfspanels)                     :: center
    REAL(WP), DIMENSION(2,4,npanels+nfspanels)                   :: cornerslocal 
    !REAL(WP), INTENT(INOUT)                                         :: cornerslocal(2,4,npanels+nfspanels)!(:,:,:)
    REAL(wp), DIMENSION(npanels+nfspanels)                       :: area
    REAL(wp), DIMENSION(3,3,npanels+nfspanels)                   :: coordsys
    REAL(WP), DIMENSION(3,npanels+nfspanels,npanels+nfspanels)   :: vel
    REAL(WP), DIMENSION(npanels+nfspanels+2,npanels+nfspanels+2) :: am
    REAL(WP), DIMENSION(3,3,npanels+nfspanels,npanels+nfspanels) :: hessglobal
    REAL(WP), DIMENSION(3,3)                                     :: hp




    write(*,*) ''
    write(*,*) 'computing linearized free surface boundary conditions..........'

    Do i=1,nfspanels

       !Get the field points on the first order free surface, include x offset 
       
       Row = npanels + i
       
       fieldpoint(:)=center(:,Row)
       
       !      !! No mirrors, but need offsets in i and k
       fieldpoint(1)=fieldpoint(1)+deltax
       fieldpoint(3)=0.
       
       fieldpointp(:) = fieldpoint(:)
       fieldpointp(2) = -fieldpoint(2)
       
       !! Calculate the contribution of body panels to the free surface condition
       !! ***********************************************************************
       Do j=1,npanels
          
          Col=j
          
          !  Call what kind of influence function?  Use the Linear free surface condition
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
          
          hessglobal(1,1,Row,Col)=hessglobal(1,1,Row,Col)+hp(1,1)
          
          !Finally, make the a(i,j) matrix components satisfying the FS BC
          am(Row,Col)= (vinf(1)**2.)*hessglobal(1,1,Row,Col) + 9.81*vel(3,Row,Col)
          
       END DO
       
       !! Calculate the contribution of free surface panels to the free surface condition
       !!********************************************************************************
       Do j=1,nfspanels
          
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
          
          hessglobal(1,1,Row,Col)=hessglobal(1,1,Row,Col)+hp(1,1)
          
          !Finally, make the a(i,j) matrix components satisfying the FS BC
          am(Row,Col)= (vinf(1)**2.)*hessglobal(1,1,Row,Col) + 9.81*vel(3,Row,Col)

       END DO
       
    END DO


  END SUBROUTINE computeLinearizedFreeSurfaceBC





  subroutine  computeNonLinearFreeSurfaceBC( npanels, nfspanels, fieldpoint,fieldpointp,center,deltax,&
       vel,dphi, coordsys, area, cornerslocal,vp,vinf, hp, hessglobal, ddphi, &
       zeta, dzeta_dx, dzeta_dy, am, &
                                       ntransompanels, ntransoml, ntransomt, nhpl)
    !
    ! Function to compute the A matrix non linear free surface Boundary Condition
    !
    INTEGER                                                      :: i,j,k
    INTEGER                                                      :: Row, Col
    INTEGER                                                      :: npanels, nfspanels
    integer :: ntransompanels, ntransoml, ntransomt
    integer :: nhpl                                              !! number of panels longitudinaly on the hull

    REAL(WP), parameter                 ::    gravity = 9.81
    REAL(WP), DIMENSION(3,3), parameter ::    hpsign  = RESHAPE( (/ &
                                                                    1., -1.,  1., &
                                                                   -1.,  1., -1., &
                                                                    1., -1.,  1.  /), (/3,3/)) 
    REAL(wp)                                                     :: deltax
    REAL(WP), DIMENSION(3)                                       :: vinf
    REAL(WP), DIMENSION(3)                                       :: fieldpoint
    REAL(WP), DIMENSION(3)                                       :: fieldpointp
    REAL(WP), DIMENSION(3)                                       :: vp
    REAL(wp), ALLOCATABLE,DIMENSION(:,:)                         :: center
    REAL(WP), ALLOCATABLE,DIMENSION(:,:,:)                      :: cornerslocal 
    REAL(wp), ALLOCATABLE,DIMENSION(:)                           :: area
    REAL(wp), ALLOCATABLE,DIMENSION(:,:,:)                       :: coordsys

    REAL(WP), ALLOCATABLE,DIMENSION(:,:)                         :: am

    REAL(WP), ALLOCATABLE,DIMENSION(:,:,:)                       :: vel
    REAL(WP), ALLOCATABLE,DIMENSION(:,:)                         :: dphi

    REAL(WP), ALLOCATABLE,DIMENSION(:,:,:,:)                     :: hessglobal
    REAL(WP), ALLOCATABLE,DIMENSION(:,:,:)                       :: ddphi
    REAL(WP), DIMENSION(3,3)                                     :: hp

    REAL(WP), ALLOCATABLE,DIMENSION(:)                           :: zeta
    REAL(WP), ALLOCATABLE,DIMENSION(:)                           :: dzeta_dx
    REAL(WP), ALLOCATABLE,DIMENSION(:)                           :: dzeta_dy

    write(*,*) ''
    write(*,*) 'computing full nonlinear free surface boundary conditions......'


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
          
          !Finally, make the a(i,j) matrix components satisfying the Full FS BC:
          !am(Row,Col)= (vinf(1)**2.)*hessglobal(1,1,Row,Col) + 9.81*vel(3,Row,Col)
          am(Row,Col)= &
               dphi(1,Row)       * ddphi(1,1,Row) * vel(1,Row,Col) + &
               (dphi(1,Row)**2)                  * hessglobal(1,1,Row,Col) + &
               dphi(1,Row)       * ddphi(2,1,Row) * vel(2,Row,Col) + &
               dphi(1,Row)       * dphi(2,Row)    * hessglobal(2,1,Row,Col) + &
               dphi(1,Row)       * ddphi(3,1,Row) * vel(3,Row,Col) + &
               dphi(1,Row)       * dphi(3,Row)    * hessglobal(3,1,Row,Col) + &
               dphi(2,Row)       * ddphi(1,2,Row) * vel(1,Row,Col) + &
               dphi(2,Row)       * dphi(1,Row)    * hessglobal(1,2,Row,Col) + &
               dphi(2,Row)       * ddphi(2,2,Row) * vel(2,Row,Col) + &
               (dphi(2,Row)**2)                   * hessglobal(2,2,Row,Col) + &
               dphi(2,Row)       *ddphi(3,2,Row)  * vel(3,Row,Col)+ &
               dphi(2,Row)       *dphi(3,Row)     * hessglobal(3,2,Row,Col)  &
               -gravity * dzeta_dx(Row) * vel(1,Row,Col)  &
               -gravity * dzeta_dy(Row) * vel(2,Row,Col) + &
                gravity * vel(3,Row,Col)

        
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
          
          !Finally, make the a(i,j) matrix components satisfying the Full FS BC
          !am(Row,Col)= (vinf(1)**2.)*hessglobal(1,1,Row,Col) + 9.81*vel(3,Row,Col)
          am(Row,Col)= &
               dphi(1,Row)       * ddphi(1,1,Row) * vel(1,Row,Col) + &
               (dphi(1,Row)**2)                  * hessglobal(1,1,Row,Col) + &
               dphi(1,Row)       * ddphi(2,1,Row) * vel(2,Row,Col) + &
               dphi(1,Row)       * dphi(2,Row)    * hessglobal(2,1,Row,Col) + &
               dphi(1,Row)       * ddphi(3,1,Row) * vel(3,Row,Col) + &
               dphi(1,Row)       * dphi(3,Row)    * hessglobal(3,1,Row,Col) + &
               dphi(2,Row)       * ddphi(1,2,Row) * vel(1,Row,Col) + &
               dphi(2,Row)       * dphi(1,Row)    * hessglobal(1,2,Row,Col) + &
               dphi(2,Row)       * ddphi(2,2,Row) * vel(2,Row,Col) + &
               (dphi(2,Row)**2)                   * hessglobal(2,2,Row,Col) + &
               dphi(2,Row)       *ddphi(3,2,Row)  * vel(3,Row,Col)+ &
               dphi(2,Row)       *dphi(3,Row)     * hessglobal(3,2,Row,Col)  &
               -gravity * dzeta_dx(Row) * vel(1,Row,Col)  &
               -gravity * dzeta_dy(Row) * vel(2,Row,Col) + &
               gravity * vel(3,Row,Col)


       END DO
       
    END DO







    print*, 'Computing transom stern panels'


    !! tack on a modification for the panels just aft of the transom
       
   Do i=1,ntransomt

      
       !Get the field points on the first order free surface, include x offset 
       
       Row = npanels + nfspanels + (i-1)*ntransoml+i

       print*,'Panel number, transom stern= ',row


       fieldpoint(:)=center(:,Row)
       
       !! No mirrors, but need offsets in i and k
       fieldpoint(1)=fieldpoint(1)+deltax  !! shift the collocation (field) points deltax
       fieldpoint(3)=zeta(Row) !0.         !! evaluate the velocities at the surface
       
       fieldpointp(:) = fieldpoint(:)
       fieldpointp(2) = -fieldpoint(2)
       
       !! Calculate the contribution of body panels to the transom condition
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
          
          !Finally, make the a(i,j) matrix components satisfying the Full FS BC:
          !am(Row,Col)= (vinf(1)**2.)*hessglobal(1,1,Row,Col) + 9.81*vel(3,Row,Col)


          !am(Row,Col)= vel(1,Row,col) !Bertram simplification

          !! full Bernoulli, page 244, Bertram:
          am(Row,Col)= vel(1,row,col)**2              + vel(1,row,col)*vel(2,row,col) + vel(1,row,col)*vel(3,row,col) + &
                       vel(2,row,col)*vel(1,row,col)  + vel(2,row,col)**2             + vel(2,row,col)*vel(3,row,col) + &
                       vel(3,row,col)* vel(1,row,col) + vel(3,row,col)*vel(2,row,col) + vel(3,row,col)**2
 
        
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
          
          !Finally, make the a(i,j) matrix components satisfying the Full FS BC
          !am(Row,Col)= (vinf(1)**2.)*hessglobal(1,1,Row,Col) + 9.81*vel(3,Row,Col)

          !am(Row,Col)=vel(1,row,col) !!Bertram Formula linear x only

          !! full Bernoulli, page 244, Bertram:
          am(Row,Col)= vel(1,row,col)**2              + vel(1,row,col)*vel(2,row,col) + vel(1,row,col)*vel(3,row,col) + &
                       vel(2,row,col)*vel(1,row,col)  + vel(2,row,col)**2             + vel(2,row,col)*vel(3,row,col) + &
                       vel(3,row,col)* vel(1,row,col) + vel(3,row,col)*vel(2,row,col) + vel(3,row,col)**2



       END DO
       
    END DO




  end subroutine computeNonLinearFreeSurfaceBC
    
END MODULE computeA
  
