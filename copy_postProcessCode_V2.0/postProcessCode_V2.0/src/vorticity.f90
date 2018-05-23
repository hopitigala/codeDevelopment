!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This module consists of all the subroutines required for computing      !!!!
!!! both instantaneous and mean vorticity feilds                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module vorticity 
  implicit none
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This subroutine computes the components of velocity gradient tensor        !!
!! INPUTS:                                                                    !!
!! u   : 4D arrary which has data of all 3 components of velocity             !!
!!        u(:,:,:,1)= w                                                       !!
!!        u(:,:,:,2)= v                                                       !!
!!        u(:,:,:,3)= u                                                       !!
!! xp  : 1D array with x-direction coordinate data                            !! 
!! ypdo: 1D array with y-direction coordinate data for discretized domain     !!
!! zp  : 1D array with x-direction coordinate data                            !!
!!                                                                            !!
!! OUTPUTS:                                                                   !!
!! dudx : 4D array that includes gradients of all 3 components of velocity    !!
!!        in x-direction                                                      !!  
!!        dudx(:,:,:,1)= dw/dx                                                !!
!!        dudx(:,:,:,2)= dv/dx                                                !!
!!        dudx(:,:,:,3)= du/dx                                                !!
!!                                                                            !!     
!! dudy : 4D array that includes gradients of all 3 components of velocity    !!
!!        in y-direction                                                      !!  
!!        dudy(:,:,:,1)= dw/dy                                                !!
!!        dudy(:,:,:,2)= dv/dy                                                !!
!!        dudy(:,:,:,3)= du/dy                                                !!
!!                                                                            !!        
!! dudz : 4D array that includes gradients of all 3 components of velocity    !!
!!        in y-direction                                                      !!  
!!        dudz(:,:,:,1)= dw/dz                                                !!
!!        dudz(:,:,:,2)= dv/dz                                                !!
!!        dudz(:,:,:,3)= du/dz                                                !!
!! central differencing scheme is used to compute gradients at internal nodes !!
!! forward and backward differencing schemes used at the end of the domains   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  subroutine velograd(u,xp,ypdo,zp,dudx,dudy,dudz)
    implicit none
    
    real*8,intent(in) ,dimension(:,:,:,:):: u
    real*8,intent(in) ,dimension(:)      :: xp,zp
    real*8,intent(in) ,dimension(0:)     :: ypdo
    real*8,intent(out),dimension(:,:,:,:):: dudx,dudy,dudz
    real*8,allocatable,dimension(:)      :: inv_x,inv_y,inv_z
    integer                              :: i,j,k,l,d1,d2,d3,d4
    
    ! obtain the sizes of arrays in each dimension
    
    d1=size(u,1)
    d2=size(u,2) ! d2=9 here
    d3=size(u,3)
    d4=size(u,4)
    
    ! find inverse of cordinate differences

    ! allocate arrays
    allocate(inv_x(d3),inv_y(d2),inv_z(d1))

    ! inverse coordinate difference for z cordinate
    
    do i=2,d1-1
       inv_z(i)=1.0/(zp(i+1)-zp(i-1))
    end do
    inv_z(1)=1.0/(zp(2)-zp(1))
    inv_z(d1)=1.0/(zp(d1)-zp(d1-1))
    
    ! inverse coordinate difference for y cordinate
    do j=2,d2-1
       inv_y(j)=1.0/(ypdo(j+1)-ypdo(j-1))
    end do
    inv_y(1)=1.0/(ypdo(2)-ypdo(1))
    inv_y(d2)=1.0/(ypdo(d2)-ypdo(d2-1))

    ! inverse coordinate difference for x cordinate
    do k=2,d3-1
       inv_x(k)=1.0/(xp(k+1)-xp(k-1))
    end do
    inv_x(1)=1.0/(xp(2)-xp(1))
    inv_x(d3)=1.0/(xp(d3)-xp(d3-1))
    
    do i=2,d1-1
       dudz(i,:,:,:)=(u(i+1,:,:,:)-u(i-1,:,:,:))*inv_z(i)
    end do

    dudz(1,:,:,:)=(u(2,:,:,:)-u(1,:,:,:))*inv_z(1)
    dudz(d1,:,:,:)=(u(d1,:,:,:)-u(d1-1,:,:,:))*inv_z(d1)
    
    do j=2,d2-1
       dudy(:,j,:,:)=(u(:,j+1,:,:)-u(:,j-1,:,:))*inv_y(j)
    end do

    dudy(:,1,:,:)=(u(:,2,:,:)-u(:,1,:,:))*inv_y(1)
    dudy(:,d2,:,:)=(u(:,d2,:,:)-u(:,d2-1,:,:))*inv_y(d2)
    
    do k=2,d3-1
       dudx(:,:,k,:)=(u(:,:,k+1,:)-u(:,:,k-1,:))*inv_x(k)
    end do

    dudx(:,:,1,:)=(u(:,:,2,:)-u(:,:,1,:))*inv_x(1)
    dudx(:,:,d3,:)=(u(:,:,d3,:)-u(:,:,d3-1,:))*inv_x(d3)

  end subroutine velograd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This subroutine calculates vorticity vector by using velocity gradients!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Omegax=dW/dy-dV/dz => Omega(I,J,K,3)=dudy(I,J,K,1)-dudz(I,J,K,2) !!!!!!!
!!Omegay=dU/dz-dW/dx => Omega(I,J,K,2)=dudz(I,J,K,3)-dudx(I,J,K,1) !!!!!!!
!!Omegaz=dV/dx-dU/dy => Omega(I,J,K,1)=dudx(I,J,K,2)-dudy(I,J,K,3) !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine instvorticity(dudz,dudy,dudx,omega)
    implicit none
    
    real*8,intent(in) ,dimension(:,:,:,:)::dudz,dudx,dudy
    real*8,intent(out),dimension(:,:,:,:)::omega

    Omega(:,:,:,3)=dudy(:,:,:,1)-dudz(:,:,:,2)!omegax
    Omega(:,:,:,2)=dudz(:,:,:,3)-dudx(:,:,:,1)!omegay
    Omega(:,:,:,1)=dudx(:,:,:,2)-dudy(:,:,:,3)!omegaz

  end subroutine instvorticity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                          !!!
!!! This subroutine computes the components of velocity gradient tensor at   !!!
!!! each location of the domain                                              !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module vorticity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                          !!!
!!! This module contians subroutine that are used to compute vortex          !!!
!!! identification methods: lambda2, Q, and lambda_ci                        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module vortexiden
  implicit none
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                          !!!
!!! This subroutine computes the second largest eigenvalue of the S**2+Omga**2!!
!!! tensor. This eigenvalue is called lambda2. lambda2 is computed at each   !!!
!!! point of the domain.                                                     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine complambda2(dudx,dudy,dudz,lambda2)
  use eigenvalues, only: symmatxeignv 
  implicit none
  
  real*8,intent(in) , dimension(:,:,:,:):: dudx,dudy,dudz
  real*8,intent(out), dimension(:,:,:)  :: lambda2  
  integer                               :: i,j,k,d1,d2,d3
  real*8,allocatable, dimension(:,:)    :: gradu,s,omga,somga2,evec
  real*8,allocatable, dimension(:)      :: eval
  ! s= rate of strain tensor
  ! omga= rate of rotation tensor
  ! somga2=s*s+omga*omga : this is tensor required to find eigenvalues.

  
  
  ! dimensions of the arrays
  d1=size(dudx,1)
  d2=size(dudx,2)
  d3=size(dudx,3)
  ! compute the components of velcotity gradient tensor
  
  do k=1,d3
     do j=1,d2
        do i=1,d1
           allocate(gradu(3,3),s(3,3),omga(3,3),somga2(3,3),evec(3,3),eval(3))
           gradu(1,1)=dudx(i,j,k,3)
           gradu(1,2)=dudx(i,j,k,2)
           gradu(1,3)=dudx(i,j,k,1)
           gradu(2,1)=dudy(i,j,k,3)
           gradu(2,2)=dudy(i,j,k,2)
           gradu(2,3)=dudy(i,j,k,1)
           gradu(3,1)=dudz(i,j,k,3)
           gradu(3,2)=dudz(i,j,k,2)
           gradu(3,3)=dudz(i,j,k,1)
           
           ! then s and omga should be computed at each point
           s = 0.5*(gradu+transpose(gradu))
           omga = 0.5* (gradu-transpose(gradu))
           
           somga2=matmul(s,s)+matmul(omga,omga)
           ! call subroutine to compute eigenvalues. see mathops.f90 file. 
           call symmatxeignv(somga2,'N','U',evec,eval)
           ! lambda2 is the second largest eigenvalue.
           lambda2(i,j,k)=eval(2)
           deallocate(gradu,s,omga,somga2,evec,eval)
        end do
     end do
  end do
end subroutine complambda2

end module vortexiden
