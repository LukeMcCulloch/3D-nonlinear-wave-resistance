module computeFlowProperties
 

  USE precise, ONLY : defaultp

  IMPLICIT NONE    
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp


CONTAINS

  SUBROUTINE computeVelocities(npanels,nfspanels,vtotal,vinf,sigma,vel,vn,vt1,vt2,coordsys)
    !! get the total velocity on each panel
    !!
    !! get the normal velocity on each panel (should==0)
    !!
    !! get the tangential components on each panel:
    !!
    
    INTEGER :: i, j
    INTEGER :: npanels, nfspanels
    
    REAL(WP), DIMENSION(3)                                       :: vinf
    REAL(WP), DIMENSION(3,npanels+nfspanels)                     :: vtotal
    REAL(WP), DIMENSION(npanels+nfspanels)                       :: sigma
    REAL(WP), DIMENSION(3,npanels+nfspanels,npanels+nfspanels)   :: vel
    REAL(WP), DIMENSION(npanels+nfspanels)                       :: vn
    REAL(WP), DIMENSION(npanels+nfspanels)                       :: vt1
    REAL(WP), DIMENSION(npanels+nfspanels)                       :: vt2
    REAL(wp), DIMENSION(3,3,npanels+nfspanels)                   :: coordsys
    
    Do i=1,npanels+nfspanels
       vtotal(:,i) = -vinf(:)
       
       !! total:
       Do j=1, npanels+nfspanels
          vtotal(:,i)=vtotal(:,i)+sigma(j)*vel(:,i,j)
       End Do
       
       !! normal:
       vn(i)=Dot_Product(coordsys(:,3,i),vtotal(:,i))
       
       !! tangents:
       vt1(i)=dot_product(coordsys(:,1,i),vtotal(:,i))
       vt2(i)=dot_product(coordsys(:,2,i),vtotal(:,i))
       
    End Do
    
  END SUBROUTINE computeVelocities



  SUBROUTINE computePressures(npanels,nfspanels,vtotal,vinf,cp)
    !! Get the pressure on each panel 
    !! from the steady Bernoulli Eq:
    INTEGER :: i, j
    INTEGER :: npanels, nfspanels

    REAL(WP), DIMENSION(3)                      :: vinf
    REAL(WP), DIMENSION(3,npanels+nfspanels)    :: vtotal    
    REAL(WP), DIMENSION(npanels+nfspanels)      :: cp

    Do i=1,npanels       
       cp(i)=1.-((sum(vtotal(:,i)**2))/(sum(vinf**2)))
    End do

  END SUBROUTINE computePressures



  SUBROUTINE checkSourceSink(npanels, nfspanels, sigma, area, sourcesink)
    !! check that the sum of sources and sinks on the panels is equal to 0.
    integer :: i
    integer :: npanels, nfspanels

    REAL(WP)                                                     :: sourcesink
    REAL(WP), DIMENSION(npanels+nfspanels)                       :: sigma
    REAL(WP), DIMENSION(npanels+nfspanels)                       :: area
 

    
    sourcesink=0.
    Do i=1,npanels+nfspanels
       sourcesink = sourcesink+sigma(i)*area(i)
    end do
    
    write(6,*) 'sourcesink sum (should be near 0.)',sourcesink
  END SUBROUTINE checkSourceSink


  subroutine computeWaveHeight(npanels,nfspanels,vinf,zeta,vtotal)
    !! subroutine to compute the wave height zeta
    integer :: i,row
    integer :: npanels, nfspanels

    REAL(WP), DIMENSION(3)                      :: vinf
    real(wp), dimension(npanels+nfspanels)      :: zeta
    REAL(WP), DIMENSION(3,npanels+nfspanels)    :: vtotal  


    zeta=0.
    Do i=1,nfspanels
       !zeta=0.
       !     vtotal(:,i) = -vinf(:)
       !Do j=1, nfspanels
       Row=i+npanels
       ! Col=j+npanels
       !zeta(Row) = zeta(Row)+ ((vinf(1)/9.81)*sigma(Col)*vel(1,Row,Col))     ! Free surface wave elevation
       !zeta(Row) = zeta(Row)+ ((vinf(1)/9.81)*vtotal(1,Row))   !4-26-10 old
       zeta(Row) = (vinf(1)/9.81)*(vtotal(1,Row) + vinf(1))
       !End Do
    End Do
  end subroutine computeWaveHeight




  subroutine computeHullForces(npanels,nfspanels,center,force,area,vinf,cp,coordsys)
    !! Coompute the forces on the hull
    integer :: i
    integer :: npanels,nfspanels

    REAL(wp), DIMENSION(3,npanels+nfspanels)    :: center
    REAL(WP), DIMENSION(3)                      :: force
    REAL(WP), DIMENSION(npanels+nfspanels)      :: area
    REAL(WP), DIMENSION(3)                      :: vinf
    REAL(WP), DIMENSION(npanels+nfspanels)      :: cp
    REAL(wp), DIMENSION(3,3,npanels+nfspanels)  :: coordsys

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
    real(wp), DIMENSION(3,npanels+nfspanels)    :: center    
    real(wp), DIMENSION(npanels+nfspanels)      :: area
    REAL(WP), DIMENSION(3)                      :: vinf

    S=0.
    Do i=1,npanels
       if(center(3,i)<0.)Then
          S=S+area(i)
       Else
       End if
    End do
    cw = -force(1)/(.5*1000.*(vinf(1)**2.)*S)
    write(6,*) 'cw', cw
    write(6,*) 'force', force
  end subroutine computeWaveResistance


  subroutine computeWaveProfile(npanels,nfspanels,profilewave,center,deltax,zeta)
    !! subroutine to compute the wave Profile
    integer :: i,row
    integer :: npanels,nfspanels

    real(wp) :: deltax
    REAL(WP), DIMENSION(npanels+nfspanels) :: zeta
    REAL(WP), DIMENSION(2,128) :: profilewave
    real(wp), DIMENSION(3,npanels+nfspanels)    :: center    
    
    do i=npanels,npanels+128
       row=(i-npanels)
       profilewave(:,row) = (/center(1,i)+deltax, zeta(i)/)
    end do
  end subroutine computeWaveProfile
  
end module computeFlowProperties
