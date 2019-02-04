MODULE A

  USE precise, ONLY : defaultp
  USE io
  USE fourthconst
  USE fifthpanel

!! variable declarations

  IMPLICIT NONE    ! makes sure the compiler complains about 
                   ! undeclared variables
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp

CONTAINS

  FUNCTION hsinfluence( fieldpoint, center, area, coordsys, cornerslocal) RESULT (vel)
    !
    ! Function to compute the del dot phi
    ! It will return a 3D velocity vector
    ! Equivlent to the far approximation
    ! of the velocity induced by a panel
    ! source q on a point p.

    !Note: area is not used!

    INTEGER :: i,j,k
    REAL(wp), DIMENSION(3) :: vel
    REAL(wp), DIMENSION(3) :: dphi           ! Page 3 of the notes
    REAL(wp), DIMENSION(3) :: fieldpoint
    REAL(wp), DIMENSION(3) :: center
    REAL(wp) :: area
    REAL(wp), DIMENSION(3,3) :: coordsys     ! local coordinate system
    REAL(wp), DIMENSION(2,4) :: cornerslocal ! quad corners in local coordinate system
    REAL(wp), DIMENSION(3) :: plocal         ! fieldpoint "p" in local coordinate system
    REAL(wp), DIMENSION(3) :: vellocal

    REAL(wp), DIMENSION(4) :: d,e,h,m,r
    real(wp) :: zero

    integer, dimension(4) :: ip1 = (/ 2, 3, 4, 1 /)

    ! transfer point p into the local coordinate system:
    plocal = MATMUL((TRANSPOSE(coordsys)),fieldpoint-center)
    
    !write(6,*) 'center', center
    !write(6,*) 'plocal', plocal 
    zero=0.

    dphi=zero



    Do, k=1,4     
      
       !write(6,*)'cornerslocal', cornerslocal(:,k)      
       d(k) = sqrt((cornerslocal(1,ip1(k))-cornerslocal(1,k))**2.+(cornerslocal(2,ip1(k))-cornerslocal(2,k))**2.)
       e(k) = (plocal(1)-cornerslocal(1,k))**2. + plocal(3)**2.
       h(k) = ((plocal(1)-cornerslocal(1,k))*(plocal(2)-cornerslocal(2,k)))
       m(k) = (((cornerslocal(2,ip1(k)))-(cornerslocal(2,k))) /(cornerslocal(1,ip1(k))-(cornerslocal(1,k))))
       r(k) = sqrt(((plocal(1))-(cornerslocal(1,k)))**2.+((plocal(2))-(cornerslocal(2,k)))**2.+(plocal(3))**2.)


       !write(6,*) 'd(k)', d(k)
       !write(6,*) 'e(k)', e(k)
       !write(6,*) 'h(k)', h(k)
       !write(6,*) 'm(k)', m(k)
       !write(6,*) 'r(k)', r(k)

     end do

     dphi(:)=0.

     do, k=1,4

       !If (abs(d(k))> 0.000001) then   
       If (abs(d(k))> 0.000000001) then
       
          dphi(1) = dphi(1)+((((cornerslocal(2,ip1(k)))-(cornerslocal(2,k)))/d(k))*log((r(k)+r(ip1(k))-d(k))/(r(k)+r(ip1(k))+d(k))))

          dphi(2) = dphi(2)-(cornerslocal(1,ip1(k))-(cornerslocal(1,k)))/d(k)*log((r(k)+r(ip1(k))-d(k))/(r(k)+r(ip1(k))+d(k)))

          dphi(3) = dphi(3)+(atan((m(k)*e(k)-h(k))/(plocal(3)*r(k)))- atan((m(k)*e(ip1(k))-h(ip1(k)))/(plocal(3)*r(ip1(k)))))

          !write(6,*) 'k', k

       End if

    End do
    dphi(1)=-1./(4.*pi)*dphi(1)
    dphi(2)=-1./(4.*pi)*dphi(2)
    dphi(3)=-1./(4.*pi)*dphi(3)
     

    !write(6,*) 'dphi(1)', dphi(1)
    !write(6,*) 'dphi(2)', dphi(2)
    !write(6,*) 'dphi(3)', dphi(3)





    vellocal=(/ dphi(1),dphi(2),dphi(3) /)
    !write(6,*) 'vellocal'
    !write(6,'(3D16.8)')  vellocal

     vel=MATMUL(coordsys,vellocal)
     !if ( (sum(vel(:)**2)) > 9999.) then
     !   vel = (/0.,0.,0./)
     !end if

     !write(6,*) 'vel'
     !write (6,'(3D16.8)') vel
     !write(6,*) ''

   END FUNCTION hsinfluence











!!
!!*****************************************************************************************************************************************
!!

  FUNCTION phixxinfluenceOLD( fieldpoint, center, area, coordsys, cornerslocal) RESULT (hessglobal)
    !
    ! Function to compute the (del outer del) dot phi
    ! It will return a 3D velocity vector
    ! Equivlent to the far approximation
    ! of the velocity induced by a panel
    ! source q on a point p.
    INTEGER :: i,j,k
    REAL(wp), DIMENSION(3)   :: vel
    REAL(wp), DIMENSION(3,3) :: ddphi           ! Page 4 of the notes
    REAL(wp), DIMENSION(3,3) :: hessglobal
    REAL(wp), DIMENSION(3)   :: fieldpoint
    REAL(wp), DIMENSION(3)   :: center
    REAL(wp)                 :: area, Laplace
    REAL(wp), DIMENSION(3,3) :: coordsys        ! local coordinate system
    REAL(wp), DIMENSION(2,4) :: cornerslocal    ! quad corners in local coordinate system
    REAL(wp), DIMENSION(3)   :: plocal          ! fieldpoint "p" in local coordinate system
    REAL(wp), DIMENSION(3)   :: vellocal
    REAL(wp) :: zero
    REAL(wp), DIMENSION(4)   :: d,e,h,m,r

    integer, dimension(4)    :: ip1 = (/ 2, 3, 4, 1 /)



    plocal = MATMUL((TRANSPOSE(coordsys)),fieldpoint-center)
    
    zero = 0.

    Do, k=1,4     
      
       !write(6,*)'cornerslocal', cornerslocal(:,k)      
       d(k) = sqrt((cornerslocal(1,ip1(k))-cornerslocal(1,k))**2.+(cornerslocal(2,ip1(k))-cornerslocal(2,k))**2.)
       e(k) = (plocal(1)-cornerslocal(1,k))**2. + plocal(3)**2.
       h(k) = ((plocal(1)-cornerslocal(1,k))*(plocal(2)-cornerslocal(2,k)))
       m(k) = (((cornerslocal(2,ip1(k)))-(cornerslocal(2,k))) /(cornerslocal(1,ip1(k))-(cornerslocal(1,k))))
       r(k) = sqrt(((plocal(1))-(cornerslocal(1,k)))**2.+((plocal(2))-(cornerslocal(2,k)))**2.+(plocal(3))**2.)

     end do

     ddphi=zero



     do, k=1,4
   
       !If (abs(d(k))> 0.000001) then
        If (abs(d(k))> 0.000000001) then



          ddphi(1,1) =ddphi(1,1) + ((2.*(cornerslocal(2,ip1(k)) - cornerslocal(2,k))/( (r(k)+r(ip1(K)))**2 - d(k)**2))  *&
               ((plocal(1)-cornerslocal(1,k))/r(k)+(plocal(1)-cornerslocal(1,ip1(k)))/r(ip1(K)) ))

          ddphi(1,2) =ddphi(1,2)+ ((-2.*(cornerslocal(1,ip1(k)) - cornerslocal(1,k))/( (r(k)+r(ip1(K)))**2 - d(k)**2)) *&
               ((plocal(1)-cornerslocal(1,k))/r(k)+(plocal(1)-cornerslocal(1,ip1(k)))/r(ip1(K)) ))

          ddphi(1,3) =ddphi(1,3)+(( plocal(3) * (cornerslocal(2,ip1(k)) - cornerslocal(2,k)) * (r(k) + r(ip1(k))))  /&
               ((r(k)*r(ip1(k)))* ( r(k)*r(ip1(k)) + ((plocal(1)-cornerslocal(1,k))*(plocal(1)-cornerslocal(1,ip1(K)))) &
               +((plocal(2)-cornerslocal(2,k))*(plocal(2)-cornerslocal(2,ip1(k))))  + (plocal(3)**2.)    )))





          ddphi(2,1) = ddphi(2,1)+((2.*(cornerslocal(2,ip1(k)) - cornerslocal(2,k))/( (r(k)+r(ip1(K)))**2. - d(k)**2.)) *&
               ((plocal(2)-cornerslocal(2,k))/r(k)+(plocal(2)-cornerslocal(2,ip1(k)))/r(ip1(K)) ))

          ddphi(2,2) =ddphi(2,2)+((-2.*(cornerslocal(1,ip1(k)) - cornerslocal(1,k))/( (r(k)+r(ip1(K)))**2. - d(k)**2.)) *&
               ((plocal(2)-cornerslocal(2,k))/r(k)+(plocal(2)-cornerslocal(2,ip1(k)))/r(ip1(K)) ))

          ddphi(2,3) =ddphi(2,3)+  ( ((-plocal(3) * (cornerslocal(1,ip1(k)) - cornerslocal(1,k))) * (r(k) + r(ip1(k))))  /&
               ((r(k)*r(ip1(k)))* ( r(k)*r(ip1(k)) + ((plocal(1)-cornerslocal(1,k))*(plocal(1)-cornerslocal(1,ip1(K)))) &
               +((plocal(2)-cornerslocal(2,k))*(plocal(2)-cornerslocal(2,ip1(k))))  + (plocal(3)**2.)    )))




          ddphi(3,1) = ddphi(3,1) + ( 2.*(cornerslocal(2,ip1(k))-cornerslocal(2,k))*( (plocal(3)/r(k))+(plocal(3)/r(ip1(k))) )&
               /( ((r(k)+r(ip1(k)) )**2.) - (d(k)**2.)) )

          ddphi(3,2) =ddphi(3,2) + (-2.*(cornerslocal(1,ip1(k)) - cornerslocal(1,k)) *(( plocal(3)/r(k)) + (plocal(3)/r(ip1(k))) )/&
               ( ((r(k)+r(ip1(K)))**2.) - (d(k)**2.) )   )

          ddphi(3,3) =ddphi(3,3) +   ((  ((plocal(1)-cornerslocal(1,k)) * (plocal(2)-cornerslocal(2,ip1(k))) &
               - ((plocal(1)-cornerslocal(1,ip1(k)))*(plocal(2)-cornerslocal(2,k))))  * (r(k)+r(ip1(k))))  &
               / ((r(k)*r(ip1(k)))* ( r(k)*r(ip1(k)) + ((plocal(1)-cornerslocal(1,k))*(plocal(1)-cornerslocal(1,ip1(K)))) &
               +((plocal(2)-cornerslocal(2,k))*(plocal(2)-cornerslocal(2,ip1(k))))  + (plocal(3)**2.)    )) )




          !ddphi(1,1) =ddphi(1,1) + ((2.*(cornerslocal(2,ip1(k)) - cornerslocal(2,k))/( (r(k)+r(ip1(K)))**2 - d(k)**2))  *&
          !     ((plocal(1)-cornerslocal(1,k))/r(k)+(plocal(1)-cornerslocal(1,ip1(k)))/r(ip1(K)) ))

          !ddphi(1,2) =ddphi(1,2)+ ((-2.*(cornerslocal(1,ip1(k)) - cornerslocal(1,k))/( (r(k)+r(ip1(K)))**2 - d(k)**2)) *&
          !     ((plocal(1)-cornerslocal(1,k))/r(k)+(plocal(1)-cornerslocal(1,ip1(k)))/r(ip1(K)) ))

          !ddphi(1,3) =ddphi(1,3)+(( plocal(3) * (cornerslocal(2,ip1(k)) - cornerslocal(2,k)) * (r(k) + r(ip1(k))))  /&
          !     ((r(k)*r(ip1(k)))* ( r(k)*r(ip1(k)) + ((plocal(1)-cornerslocal(1,k))*(plocal(1)-cornerslocal(1,ip1(K)))) &
          !     +((plocal(2)-cornerslocal(2,k))*(plocal(2)-cornerslocal(2,ip1(k))))  + (plocal(3)**2)    )))






       End if

    End do

    ddphi = -1./(4.*pi) * ddphi

    hessglobal = Matmul(Matmul(coordsys,ddphi),Transpose(coordsys))

    



!!$    write(6,*) ''
!!$    write(6,*) 'The Hessian'
!!$    write(6,*) ddphi(1,1), ddphi(1,2), ddphi(1,3)
!!$    write(6,*) ddphi(2,1), ddphi(2,2), ddphi(2,3)
!!$    write(6,*) ddphi(3,1), ddphi(3,2), ddphi(3,3)

!!$    write(6,*) ''
!!$    Laplace = ddphi(1,1)+ddphi(2,2)+ddphi(3,3)
!!$    write(6,*) 'The Laplacian is...'
!!$    write(6,*) Laplace


!!$    write(6,*) ''
!!$    write(6,*) 'Check the Hessian'
!!$    write(6,*) ddphi(1,2), ddphi(2,1)
!!$    write(6,*) ddphi(1,3), ddphi(3,1)
!!$    write(6,*) ddphi(2,3), ddphi(3,2)



    !vellocal=(/ dphi(1),dphi(2),dphi(3) /)
    !write(6,*) 'vellocal'
    !write(6,'(3D16.8)')  vellocal

     !vel=MATMUL(coordsys,vellocal)
     !if ( (sum(vel(:)**2)) > 9999.) then
     !   vel = (/0.,0.,0./)
     !end if

     !write(6,*) 'vel'
     !write (6,'(3D16.8)') vel
     !write(6,*) ''

   END FUNCTION phixxinfluenceOLD





   FUNCTION phixxinfluence (point, center, area, coordsys, qlocal) result (hessmatrix)

     implicit none
     ! fuildpoint, center of panel, result
       real(8),dimension(3) :: point, center, Plocal
     ! matrix T with voectors i, j, k for panel
       real(8),dimension(3,3) :: coordsys
     ! panel comes in local coordinates
       real(8),dimension(2,4) :: qlocal
     ! the x and y value  
       real(8) :: area
       real(8), dimension (4) :: corner_x ! the x component of each corner
       real(8), dimension (4) :: corner_y ! the y component of each corner 
     ! pi
       real(8), parameter :: pi= 4.d0*datan(1.d0)
     ! name k+1
       integer, dimension(4), parameter :: kp1 = (/ 2, 3, 4, 1 /)
       integer :: k
       real(8), dimension (4) :: d ! length of panel side k between point q_k and q_k+1
       real(8), dimension (4) :: m ! ratios of y and x coordinate difference of panel corners (local coordinates)
       real(8), dimension (4) :: r ! distance from field point to the 4 corners
       real(8), dimension (4) :: e
       real(8), dimension (4) :: h
       real(8), dimension (4) :: rho
       real(8), dimension (4) :: lambda
       real(8) :: phi_xx, phi_yx, phi_zx, phi_xy
       real(8) :: phi_yy, phi_zy, phi_xz, phi_yz, phi_zz
       real(8), dimension (3,3) :: hesslocal ! hessian matrix in local coordinate
       real(8), dimension (3,3) :: hessmatrix !hessian matrix in global coordinate
  
     ! Transfer the point P into the local coordinate system
       Plocal = matmul( transpose( coordsys), point-center)

      ! x value
      corner_x= qlocal(1,:)
      ! y value
      corner_y= qlocal(2,:)
      
      ! the length of panel side k between point q(k) and q(k+1)
      do k = 1,4
         d(k) = dsqrt((corner_x(kp1(k))-corner_x(k))**2 &
                + (corner_y(kp1(k))-corner_y(k))**2)
         
         ! ratio of y and x coordinate differences of panel corners
         if (d(k)==0) then 
              m(k)= 0 
         else
              m(k) = (corner_y(kp1(k))-corner_y(k)) &
                      & / (corner_x(kp1(k))-corner_x(k))
         end if
      end do
      
      do k = 1,4
         ! distance of field points form corner k
         ! x is Plocal(1)
         ! y is Plocal(2)
         ! z is Plocal(3)
         r(k)   = dsqrt( (Plocal(1)-corner_x(k))**2 &
                   & + (Plocal(2)-corner_y(k))**2 &
                   & + Plocal(3)**2)
         e(k)   = (Plocal(1)-corner_x(k))**2 + (Plocal(3))**2
         h(k)   = (Plocal(1)-corner_x(k)) * (Plocal(2)-corner_y(k))
      end do  
      
      do k = 1,4
        rho(k) = r(k)*r(kp1(k)) + &
               & (Plocal(1)-corner_x(k))*(Plocal(1)- corner_x(kp1(k)))&
               & + (Plocal(2)-corner_y(k))*(Plocal(2)-corner_y(kp1(k)))&
               & + Plocal(3)**2
       lambda(k) = ((Plocal(1)-corner_x(k))*(Plocal(2)-corner_y(kp1(k)))&
                & - (Plocal(1) - corner_x(kp1(k)))*(Plocal(2)-corner_y(k)))
      end do
      
      ! the starting value of the second order dericatives
       phi_xx = 0.
       phi_xy = 0.
       phi_xz = 0.
       phi_yx = 0.
       phi_yy = 0.
       phi_yz = 0.
       phi_zx = 0.
       phi_zy = 0.
       phi_zz = 0.       

       do k = 1,4
         phi_xx = phi_xx + 2*(corner_y(kp1(k))-corner_y(k)) / &
                & ((r(k)+r(kp1(k)))**2-d(k)**2)      &
                & *((Plocal(1)-corner_x(k))/r(k) +   &
                & (Plocal(1)-corner_x(kp1(k)))/r(kp1(k)))
       end do
       phi_xx = -1./(4.*pi)*phi_xx
       
       do k= 1,4
         phi_yx = phi_yx + (-2)*(corner_x(kp1(k))-corner_x(k)) &
                & / ((r(k)+r(kp1(k)))**2-d(k)**2)   &
                & *((Plocal(1)-corner_x(k))/r(k)    &
                & + (Plocal(1)-corner_x(kp1(k)))/r(kp1(k)))
       end do
       phi_yx = -1./(4.*pi)*phi_yx
       
       do k = 1,4
         phi_zx = phi_zx + Plocal(3)*(corner_y(kp1(k))-corner_y(k))&
               & *(r(k)+r(kp1(k)))/(r(k)*r(kp1(k))) / rho(k)
       end do
       phi_zx = -1./(4.*pi)*phi_zx    
       
       do k = 1,4
         phi_xy = phi_xy + 2*(corner_y(kp1(k))-corner_y(k)) &
                & /( (r(k)+r(kp1(k)))**2-d(k)**2 )          &
                & *( (Plocal(2)-corner_y(k))/r(k) &
                & +(Plocal(2)-corner_y(kp1(k)))/r(kp1(k)) )
       end do
       phi_xy = -1./(4.*pi)*phi_xy
       
       do k = 1,4
         phi_yy = phi_yy + (-2*(corner_x(kp1(k))-corner_x(k))) &
                & /( (r(k)+r(kp1(k)))**2-d(k)**2 )             &
                & *( (Plocal(2)-corner_y(k))/r(k)              &
                & + (Plocal(2)-corner_y(kp1(k)))/r(kp1(k)) )
       end do
       phi_yy = -1./(4.*pi)*phi_yy 
       
       do k = 1,4
         phi_zy = phi_zy - Plocal(3)*(corner_x(kp1(k))-corner_x(k)) &
                & *(r(k)+r(kp1(k)))/(r(k)*r(kp1(k))) / rho(k)
       end do
       phi_zy = -1./(4.*pi)*phi_zy    
       
       do k = 1,4
         phi_xz = phi_xz + 2*(corner_y(kp1(k))-corner_y(k)) &
                & /( (r(k)+r(kp1(k)))**2-d(k)**2 )          &
                & *( Plocal(3)/r(k) + Plocal(3)/r(kp1(k)) )
       end do
       phi_xz = -1./(4.*pi)*phi_xz
       
       do k = 1,4
         phi_yz = phi_yz - 2*(corner_x(kp1(k))-corner_x(k)) &
                & /( (r(k)+r(kp1(k)))**2-d(k)**2 )          &
                & *( Plocal(3)/r(k)+Plocal(3)/r(kp1(k)) )
       end do
       phi_yz = -1./(4.*pi)*phi_yz
       
       do k = 1,4
         phi_zz = phi_zz+lambda(k)*(r(k)+r(kp1(k))) &
                & /(r(k)*r(kp1(k))*rho(k))
       end do
       phi_zz = -1./(4.*pi)*phi_zz
       
       hesslocal(1,1) = phi_xx
       hesslocal(1,2) = phi_xy
       hesslocal(1,3) = phi_xz
       hesslocal(2,1) = phi_yx
       hesslocal(2,2) = phi_yy
       hesslocal(2,3) = phi_yz
       hesslocal(3,1) = phi_zx
       hesslocal(3,2) = phi_zy
       hesslocal(3,3) = phi_zz
   
       hessmatrix = matmul(coordsys,hesslocal)
       hessmatrix = matmul(hessmatrix,transpose(coordsys))
       
       ! return to caller
       return


  END FUNCTION phixxinfluence

END MODULE A
