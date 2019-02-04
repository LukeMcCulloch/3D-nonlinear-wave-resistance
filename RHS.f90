MODULE RHS
  !Module to compute the right hand side of the Panel Code

  use precise, only : defaultp 

  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp

CONTAINS

  SUBROUTINE computeRHS(npanels, nfspanels, coordsys, vinf, b)
    INTEGER                                     :: i
    INTEGER                                     :: Row, Col
    INTEGER                                     :: npanels, nfspanels
    REAL(WP), DIMENSION(3)                      :: vinf
    REAL(wp), DIMENSION(3,3,npanels+nfspanels)  :: coordsys
    REAL(WP), DIMENSION(npanels+nfspanels)      :: b

    write(*,*) ''
    write(*,*) 'Computing RHS.......'

    Do i=1, npanels
       b(i)= DOT_PRODUCT(coordsys(:,3,i),vinf)   
    End Do
    Do i=1, nfspanels
       row=i+npanels
       b(row)=0.
    End Do
  END SUBROUTINE computeRHS


  SUBROUTINE computeNonLinearRHS(npanels, nfspanels, coordsys, vinf, &
       dphi, ddphi, b, &
                                       ntransompanels, ntransoml, ntransomt,nhpl,center)
    INTEGER                                                      :: i
    INTEGER                                                      :: Row, Col
    INTEGER                                                      :: npanels, nfspanels
    integer :: ntransompanels, ntransoml, ntransomt
    integer :: nhpl,hullpt
    REAL(WP), parameter                                          :: gravity = 9.81
    REAL(WP), DIMENSION(3)                                       :: vinf
    REAL(wp), ALLOCATABLE,DIMENSION(:,:,:)                       :: coordsys
    REAL(WP), ALLOCATABLE,DIMENSION(:)                           :: b
    REAL(WP), ALLOCATABLE,DIMENSION(:,:)                         :: dphi
    REAL(WP), ALLOCATABLE,DIMENSION(:,:,:)                       :: ddphi
    REAL(wp), ALLOCATABLE,DIMENSION(:,:)                         :: center

    write(*,*) ''
    write(*,*) 'Computing full nonlinear RHS.......'

    Do i=1, npanels
       !b(i)= -DOT_PRODUCT(coordsys(:,3,i),dphi(:,i)) 
       b(i)= DOT_PRODUCT(coordsys(:,3,i),vinf)
    End Do
    Do i=1, nfspanels+ntransompanels
       row=i+npanels
       b(row) = -((dphi(1,Row)**2) * ddphi(1,1,Row) + &
                  dphi(1,Row)*dphi(2,Row)*ddphi(2,1,Row) + &
                  dphi(1,Row)*dphi(3,Row)*ddphi(3,1,Row) + &
                  dphi(2,Row)*dphi(1,Row)*ddphi(1,2,Row) + &
                  (dphi(2,Row)**2)*ddphi(2,2,Row) + &
                  dphi(2,Row)*dphi(3,Row)*ddphi(3,2,Row) + &
                  gravity*dphi(3,Row))



    End Do

    Do i=1,ntransomt
       !Get the field points on the first order free surface, include x offset 
       row = npanels + nfspanels + (i-1)*ntransoml+1
       hullpt = nhpl + (i-1)*nhpl
       !b(row) = vinf(1)-sqrt(max(0.,vinf(1)**2.-2.*gravity*center(3,hullpt)))  !!Bertram x-dir linearization
       b(row) = Vinf(1)**2 - 2.*gravity*center(3,hullpt) !!full Bernoulli streamline for the flow
    End Do


  END SUBROUTINE computeNonLinearRHS


END MODULE RHS
