module anisotropyroutines

     implicit none
contains

!!!!!               Finds the reynolds stress matrix u_ij at all points     !!!!!
                         !! input :: velocity fluctuations  !!
                         !! output:: instantaneous Reynolds stress value !!

     subroutine reystinv(uu,uv,uw,vv,vw,ww,eta,zeta)
         use channhdvariables
         use prelimcalvar
         use mainparameters
         use mpi
         use mpivariables
         implicit none
         real*8,intent(in),dimension(:,:,:)               :: uu,uv,uw,vv,vw,ww
         real*8,intent(out),dimension(:,:,:)              :: eta,zeta
         integer                                          :: d1,d2,d3,xco,yco,zco,i,j
         real*8,allocatable,dimension(:,:)                :: b,b2,b3
         real*8                                           :: ke,inv,invari2,invari3,fac
         d1 = size(uv,1)
         d2 = size(uv,2)    
         d3 = size(uv,3)
         
         inv = 1.0/3.0
         fac = 1.0/6.0

         do zco = 1,d1
              do yco = 1,d2
                  do xco = 1,d3
                       allocate(b(3,3),b2(3,3),b3(3,3))
                       ke = uu(zco,yco,xco)+vv(zco,yco,xco)+ww(zco,yco,xco)
                       b(1,1) = (uu(zco,yco,xco)/ke)-inv
                       b(1,2) = (uv(zco,yco,xco)/ke)
                       b(1,3) = (uw(zco,yco,xco)/ke)
                       b(2,1) = (uv(zco,yco,xco)/ke)
                       b(2,2) = (vv(zco,yco,xco)/ke)-inv
                       b(2,3) = (vw(zco,yco,xco)/ke)
                       b(3,1) = (uw(zco,yco,xco)/ke)
                       b(3,2) = (vw(zco,yco,xco)/ke)
                       b(3,3) = (ww(zco,yco,xco)/ke)-inv
 
                       b2 = matmul(b,b)
                       b3 = matmul(b2,b)

!                       if(xco.eq.140.and.yco.eq.4.and.zco.eq.7)then
!                              write(*,*)b                              
!                              write(*,*)b2
!                              write(*,*)b3
!                              write(*,*)'the kinetic energy is',ke
!                              write(*,*)'The first invariant is',b(1,1)+b(2,2)+b(3,3)
!                       end if
                       invari2 = b2(1,1)+b2(2,2)+b2(3,3)
                       invari3 = b3(1,1)+b3(2,2)+b3(3,3)
                       eta(zco,yco,xco) = fac*invari2
                       zeta(zco,yco,xco) = fac*invari3
                       deallocate(b,b2,b3)
                   end do
              end do
         end do
      end subroutine reystinv

end module anisotropyroutines
                      
         
               
