!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! This modules of this file contains all the mathematical operations !!!!!!!
!!!!! required for the code implementations                              !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!****************************************************************************!!
!!****************************************************************************!!
!! This module includes all the derivatives input and out put arrays are 3D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module derivatives
  implicit none
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This subroutine calculates the first derivative of any variable over the !!!!
!! 3D domain. The input vaiable array is 3D and the coordinate array is 1D  !!!!
!! The output 3D array contains the derivative                              !!!!
!! dimen=1 dudz
!! dimen=2 dudy
!! dimen=3 dudx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine fstdervt_3d(u,x,dimen,method,dudx)
    implicit none
    
    real*8          ,intent(in) ,dimension(:,:,:) ::u
    real*8          ,intent(in) ,dimension(:)     ::x
    integer         ,intent(in)                   ::dimen
    character(LEN=*),intent(in)                   ::method
    real*8          ,intent(out),dimension(:,:,:) ::dudx
    integer                                       ::d1,d2,d3,i,j,k,dimx
    real*8, allocatable, dimension(:)             ::inv_x
    
    d1=size(u,1)
    d2=size(u,2)
    d3=size(u,3)

    dimx=size(x,1)
    
    do i=2,dimx-1
       inv_x(i)=1.0/(x(i+1)-x(i-1))
    end do
    
    if (method=='cen') then
       if (dimen==1)then
          dudx(1,:,:)=(u(2,:,:)-u(1,:,:))/(x(2)-x(1))
          do i= 2,d1-1
             dudx(i,:,:)=(u(i+1,:,:)-u(i-1,:,:))*inv_x(i)
          end do
          dudx(d1,:,:)=(u(d1,:,:)-u(d1-1,:,:))/(x(d1)-x(d1-1))
       elseif(dimen==2) then
          dudx(:,1,:)=(u(:,2,:)-u(:,1,:))/(x(2)-x(1))
          do j= 2,d2-1
             dudx(:,j,:)=(u(:,j+1,:)-u(:,j-1,:))*inv_x(j)
          end do
          dudx(:,d2,:)=(u(:,d2,:)-u(:,d2-1,:))/(x(d2)-x(d2-1))
       elseif(dimen==3) then
          dudx(:,:,1)=(u(:,:,2)-u(:,:,1))/(x(2)-x(1))
          do k= 2,d3-1
             dudx(:,:,k)=(u(:,:,k+1)-u(:,:,k-1))*inv_x(k)
          end do
          dudx(:,:,d3)=(u(:,:,d3)-u(:,:,d3-1))/(x(d3)-x(d3-1))
       end if
    end if
  end subroutine fstdervt_3d
  
end module derivatives

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This module contains subroutines that compute eigenvalues of different     !!
!! matrices.                                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module eigenvalues
  implicit none
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! This subroutine computes the eigenvalues of a real symmetric matrix.!!!!!!!
!!!! Routines from intel mkl library is used for the optimized performance !!!!!
!!!! of the code when it is used on intel processors.                       !!!!
!!!! INPUTS                                                                 !!!!
!!!! matx: This should be a real symmetric matrix                           !!!!
!!!! jobz: A character variable that decide whether eigenvectors should be  !!!!
!!!!       computed.                                                        !!!!
!!!!       jobz = 'V' - when eigenvectors need to be computed               !!!!
!!!!       jobz = 'N' - when eigenvectors are not computed. Only eigenvalues !!!
!!!!                    are computed                                        !!!!
!!!! uplo: This character variable tells whether the input is upper         !!!!
!!!!       or lower triangular part of the symmetric matrix.                !!!!
!!!!       uplo = 'U' - upper triangular part is fed.                       !!!!
!!!!       uplo = 'L' - lower triangular part is fed.                       !!!!
!!!!                                                                        !!!!
!!!! OUTPUTS                                                                !!!!
!!!! evec : 2D array that contains eigenvectors                             !!!!
!!!! eval : 1D array that contains eigenvalues. Eignevalues are ordered in  !!!!
!!!!        ascending order.                                                !!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine symmatxeignv(matx,jobz,uplo,evec,eval)
    implicit none

    include 'mkl.fi'
    !  include 'lapack.f90'
    real*8          ,intent(in) , dimension(:,:)  ::matx
    character(len=1),intent(in)                   ::jobz,uplo
    real*8          ,intent(out), dimension(:,:)  ::evec
    real*8          ,intent(out), dimension(:)    ::eval
    ! define variable that are required for mkl routine
    integer                          :: n,lda,lwork,liwork,info,i,j
    real*8,allocatable,dimension(:,:):: a
    real*8,allocatable,dimension(:)  :: work,w
    integer,allocatable,dimension(:) :: iwork
    ! first find the order of the matrix
    n=size(matx,1)

!    write(*,*)'order of the matrix (n) = ',n

    lda=max(1,n)

 !   write(*,*)'leading dimension of the matrix (lda) = ',lda

    if(n.gt.1)then
       if(jobz=='N')then
          lwork=2*n+1
          liwork=1
       elseif(jobz=='V')then
          lwork=2*n*n+6*n+1
          liwork=5*n+3
       end if
    else
       lwork=1
       liwork=1
    end if

    allocate(a(lda,max(1,n)),work(lwork),w(lda),iwork(liwork))
    !define matrix a

    if (uplo=='U')then
       do j=1,n
          do i=1,j
             a(i,j)=matx(i,j)
          end do
       end do
    elseif (uplo=='L') then
       do j=1,n
          do i=j,n
             a(i,j)=matx(i,j)
          end do
       end do
    end if

    call dsyevd(jobz,uplo,n,a,lda,w,work,lwork,iwork,liwork,info)
    do i=1,lda
       evec(:,i)=a(:,lda+1 -i)
    end do
    do i=1,lda
       eval(i)=w(1+lda-i)
    end do
    deallocate(a,w,work,iwork)
  end subroutine symmatxeignv
end module eigenvalues

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! This module contains subroutines that are used to integrate functions!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module integration
  implicit none
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! This subroutine tripple integrate a 3D function over its all the      !!!!
!!!!!! dimensions of the array.                                              !!!!
!!!!!!                                                                       !!!!
!!!!!! INPUTS:                                                               !!!!
!!!!!!       func3d  : 3d array of elements that show the value              !!!!
!!!!!!                 of the function                                       !!!!
!!!!!!       x1,x2,x3: coordinate arrays in three directions                 !!!! 
!!!!!!       dim1    : dimension of the 3d array over which the 1st          !!!!
!!!!!!                 integration is carried out.                           !!!!
!!!!!!       dim2    : dimension of the 2d array resulted from the 1st       !!!!
!!!!!!                 integration over which the 2nd integration is         !!!!
!!!!!!                 carried out.                                          !!!!
!!!!!! Note: x1 and x2 should be consistent with dim1 and dim2               !!!!
!!!!!!                                                                       !!!!
!!!!!! OUTPUTS:                                                              !!!!
!!!!!!        sumd: the result of the 3D integration.                        !!!!
!!!!!! Other variables:                                                      !!!! 
!!!!!!        sum2d: 2D array which is used to store the resultant 2D array  !!!!
!!!!!!               from the 1st integration.                               !!!!
!!!!!!        sum1d: 1D array which is used to store the resultant 1d array  !!!!
!!!!!!               from the 2nd integration.                               !!!!
!!!!!!       d1,d2,d3: dimensions of the 3D array                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  subroutine trapz3d_sum(func3d,x1,x2,x3,dim1,dim2,sumd)
!    implicit none
    
!    real*8,intent(in),dimension(:,:,:)::func3d
!    real*8,intent(in),dimension(:)    ::x1,x2,x3
!    integer,intent(in)                ::dim1,dim2
!    real*8,intent(out)                ::sumd

!    real*8,allocatable,dimension(:,:) ::sum2d
!    integer                           ::d1,d2,d3
    
    ! find the sizes of the input 3D array
!    d1=size(func3d,1)
!    d2=size(func3d,2)
!    d3=size(func3d,3)

    ! allocate sum2d and sum1d arrays according to
    ! dim1 and dim2
    
!    if(dim1==1)then
!       allocate(sum2d(d2,d3))
!    elseif(dim1==2)then
!       allocate(sum2d(d1,d3))
!    else
!       allocate(sum2d(d1,d2))
!    end if
    
    ! 1st integration 3D array to 2D array
!    call trapz3d_2dsum(func3d,x1,dim1,sum2d)
    ! 2nd integration 2D array to 1D array
!    call trapz2d_sum(sum2d,x1,x2,dim2,sumd)
    ! 3rd integration 1D array to a scalar
   ! call trapz1d_sum(sum1d,x3,sumd)
!  end subroutine trapz3d_sum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! This subroutine double integrate a 2D function over its all the       !!!!
!!!!!! dimensions of the array.                                              !!!!
!!!!!!                                                                       !!!!
!!!!!! INPUTS:                                                               !!!!
!!!!!!       func2d  : 2d array of elements that show the value              !!!!
!!!!!!                 of the function                                       !!!!
!!!!!!       x1,x2   : coordinate arrays in two directions.                  !!!! 
!!!!!!       dim1    : dimension of the 2d array over which the 1st          !!!!
!!!!!!                 integration is carried out.                           !!!!
!!!!!! Note: x1should be consistent with dim1.                               !!!!
!!!!!!                                                                       !!!!
!!!!!! OUTPUTS:                                                              !!!!
!!!!!!        sumd: the result of the 3D integration.                        !!!!
!!!!!! Other variables:                                                      !!!! 
!!!!!!                                                                       !!!!
!!!!!!        sum1d: 1D array which is used to store the resultant 1d array  !!!!
!!!!!!               from the 1st integration.                               !!!!
!!!!!!        d1,d2: dimensions of the 2D array.                             !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  subroutine trapz2d_sum(func2d,x1,x2,dim1,sumd)
!    implicit none
    
!    real*8,intent(in),dimension(:,:)::func2d
!    real*8,intent(in),dimension(:)    ::x1,x2
!    integer,intent(in)                ::dim1
!    real*8,intent(out)                ::sumd

!    real*8,allocatable,dimension(:)   ::sum1d
!    integer                           ::d1,d2
    
    ! find the sizes of the input 2D array
!    d1=size(func2d,1)
!    d2=size(func2d,2)

    ! allocate sum1d and sum1d arrays according to dim
    
!    if(dim1==1)then
!       allocate(sum1d(d2))
!    else
!       allocate(sum1d(d1))
!    end if
    
    ! 1st integration 2D array to 1D array
!    call trapz2d_1dsum(func2d,x1,dim1,sum1d)
    ! 2nd integration 1D array to a scalar
!    call trapz1d_sum(sum1d,x2,sumd)
!  end subroutine trapz2d_sum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! This subroutine integrate a 3D function over his 3rd dimension of the!!!!!
!!!!!! array.                                                                !!!!
!!!!!!                                                                       !!!!
!!!!!! INPUTS:                                                               !!!!
!!!!!!       func3d: 3d array of elements that show the value of the function!!!!
!!!!!!       x     : coordinate array                                        !!!! 
!!!!!!       dimen   : dimension of the array over which the integration is    !!!!
!!!!!!               carried out.                                            !!!!
!!!!!!                                                                       !!!!
!!!!!! OUTPUTS:                                                              !!!!
!!!!!!        sum2d: 2d array that shows the summation over the 3rd dimension!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  subroutine trapz3d_2dsum(func3d,x,dimen,sum2d)
!    implicit none
!    real*8,intent(in),dimension(:,:,:)::func3d
!    real*8,intent(in),dimension(:)::x
!    integer,intent(in)::dimen
!    real*8,intent(out),dimension(:,:)::sum2d
!    integer::i,d1,d2,d3
!    d1=size(func3d,1)
!    d2=size(func3d,2)
!    d3=size(func3d,3)
    
!    sum2d=0.0
!    if(dimen==1)then
!       do i=1,d1-1
!          sum2d=sum2d+0.5*(x(i+1)-x(i))*(func3d(i,:,:)+func3d(i+1,:,:))
!       end do
!    elseif(dimen==2)then
!       do i=1,d2-1
!          sum2d=sum2d+0.5*(x(i+1)-x(i))*(func3d(:,i,:)+func3d(:,i+1,:))
!       end do
!    else
!       do i=1,d3-1
!          sum2d=sum2d+0.5*(x(i+1)-x(i))*(func3d(:,:,i)+func3d(:,:,i+1))
!       end do
!    end if
!  end subroutine trapz3d_2dsum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! This subroutine integrate a 2D function over his 3rd dimension of the!!!!!
!!!!!! array.                                                                !!!!
!!!!!!                                                                       !!!!
!!!!!! INPUTS:                                                               !!!!
!!!!!!       func2d: 2d array of elements that show the value of the function!!!!
!!!!!!       x     : coordinate array                                        !!!! 
!!!!!!       dimen   : dimension of the array over which the integration is    !!!!
!!!!!!               carried out.                                            !!!!
!!!!!!                                                                       !!!!
!!!!!! OUTPUTS:                                                              !!!!
!!!!!!        sum1d: 1d array that shows the summation over the 3rd dimension!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  subroutine trapz2d_1dsum(func2d,x,dimen,sum1d)
!    implicit none
!    real*8,intent(in),dimension(:,:)::func2d
!    real*8,intent(in),dimension(:)::x
!    integer,intent(in)::dimen
!    real*8,intent(out),dimension(:)::sum1d
!    integer::i,d1,d2
!    d1=size(func2d,1)
!    d2=size(func2d,2)
    
!    sum1d=0.0
!    if(dimen==1)then
!       do i=1,d1-1
!          sum1d=sum1d+0.5*(x(i+1)-x(i))*(func2d(i,:)+func2d(i+1,:))
!       end do
!    else
!       do i=1,d2-1
!          sum1d=sum1d+0.5*(x(i+1)-x(i))*(func2d(:,i)+func2d(:,i+1))
!       end do
!    end if
!  end subroutine trapz2d_1dsum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! This subroutine integrate a 1D function over his 3rd dimension of the!!!!!
!!!!!! array.                                                                !!!!
!!!!!!                                                                       !!!!
!!!!!! INPUTS:                                                               !!!!
!!!!!!       func1d: 1d array of elements that show the value of the function!!!!
!!!!!!       x     : coordinate array                                        !!!! 
!!!!!!       dimen   : dimension of the array over which the integration is    !!!!
!!!!!!               carried out.                                            !!!!
!!!!!!                                                                       !!!!
!!!!!! OUTPUTS:                                                              !!!!
!!!!!!        sumd: 1d array that shows the summation over the 3rd dimension!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  subroutine trapz1d_sum(func1d,x,sumd)
!    implicit none
!    real*8,intent(in),dimension(:)::func1d
!    real*8,intent(in),dimension(:)::x
!    real*8,intent(out)::sumd
!    integer::i,d1
!    d1=size(func1d,1)
    
!    sumd=0.0
!    do i=1,d1-1
!       sumd=sumd+0.5*(x(i+1)-x(i))*(func1d(i)+func1d(i+1))
!    end do
!  end subroutine trapz1d_sum
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! These routines are copied from integrations.f90!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine trapz3d_eq(func3d,d1,d2,d3,xi,xf,sum2d)
  use channhdvariables
  use prelimcalvar
  implicit none
  integer,intent(in)::d1,d2,d3
  real*8,intent(in)::func3d(d1,d2,d3)
  real*8,intent(out)::sum2d(d1,d2)
  real*8,intent(in)::xi,xf
  !integer::k
  real*8::delx
  integer::k
  ! integrate over the 3rd dimension of the array
 ! write(*,*)'rank',mynode,'is in trapz3d'
  delx=(xf-xi)/real(d3-1)
  sum2d=0.0
  sum2d=0.5*(func3d(:,:,1)+func3d(:,:,d3))
  do k=2,(d3-1)
     sum2d=sum2d+func3d(:,:,k)
  end do
  sum2d=sum2d*delx
end subroutine trapz3d_eq

subroutine trapz2d_eq(func2d,d1,d2,a,b,sum1d)
  implicit none
  integer,intent(in)::d1,d2
  real*8,intent(in)::func2d(d1,d2)
  real*8,intent(out)::sum1d(d1)
  real*8,intent(in)::a,b
  integer::i
  real*8::delz
 ! integrate over the 1st dimension of the array
  delz=(b-a)/real(d2-1)
  sum1d=0.0
  sum1d=sum1d+0.5*(func2d(:,1)+func2d(:,d2))
  do i=2,(d2-1)
     sum1d=sum1d+func2d(:,i)
  end do
  sum1d=sum1d*delz
end subroutine trapz2d_eq

subroutine trapz2d_yz(func2d,d1,d2,a,b,sum1d)
  implicit none
  integer,intent(in)::d1,d2
  real*8,intent(in)::func2d(d1,d2)
  real*8,intent(out)::sum1d(d2)
  real*8,intent(in)::a,b
  integer::i
  real*8::delz
 ! integrate over the 1st dimension of the array
  delz=(b-a)/real(d1-1)
  sum1d=0.0
  sum1d=sum1d+0.5*(func2d(1,:)+func2d(d1,:))
  do i=2,(d1-1)
     sum1d=sum1d+func2d(i,:)
  end do
  sum1d=sum1d*delz
end subroutine trapz2d_yz

subroutine trapz1d_eq(func1d,d1,zi,zf,summ)
  !use MainVariables
  implicit none
  integer,intent(in)::d1
  real*8,intent(in)::func1d(d1)
  real*8,intent(in)::zi,zf
  real*8,intent(out)::summ
  integer::i
  real*8::delz
  !d2=size(func1d)
  !write(*,*)'rank',mynode,'is in trapz1d'
  delz=(zf-zi)/real(d1-1)
  summ=0.0
  summ=summ+0.5*(func1d(1)+func1d(d1))
  do i=2,d1-1
     summ=summ+func1d(i)
  end do
  summ=summ*delz
end subroutine trapz1d_eq

subroutine trapz2d_noneq(func2d,d1,d2,y,sum1d)
  !use MainVariables
  implicit none
  integer,intent(in)::d1,d2
  real*8,intent(in)::func2d(d1,d2)
  real*8,intent(in)::y(0:d2)
  real*8,intent(out)::sum1d(d1)
  integer::j
  !write(*,*)'rank',mynode,'is in trapz2d'
  sum1d=0.0
  sum1d=sum1d+0.5*(func2d(:,1)*(y(2)-y(1))+func2d(:,d2)*(y(d2)-y(d2-1)))
  !write(*,*)'sum1d assigned',mynode
  do j=2,d2-1
     sum1d=sum1d+func2d(:,j)*(y(j)-y(j-1))
  end do
  !write(*,*)'sum1d computed',mynode
end subroutine trapz2d_noneq

subroutine trapz2d_noneq_xy(func2d,d1,d2,y,sum1d)
  !use MainVariables
  implicit none
  integer,intent(in)::d1,d2
  real*8,intent(in)::func2d(d1,d2)
  real*8,intent(in)::y(0:d1)
  real*8,intent(out)::sum1d(d2)
  integer::j
  !write(*,*)'rank',mynode,'is in trapz2d'
  sum1d=0.0
  sum1d=sum1d+0.5*(func2d(1,:)*(y(2)-y(1))+func2d(d1,:)*(y(d1)-y(d1-1)))
  !write(*,*)'sum1d assigned',mynode
  do j=2,d1-1
     sum1d=sum1d+func2d(j,:)*(y(j)-y(j-1))
  end do
  !write(*,*)'sum1d computed',mynode
end subroutine trapz2d_noneq_xy

subroutine trapz1d_eq_xy(func1d,d2,xi,xf,summ)
  !use MainVariables
  implicit none
  integer,intent(in)::d2
  real*8,intent(in)::func1d(d2)
  real*8,intent(in)::xi,xf
  real*8,intent(out)::summ
  integer::i
  real*8::delx
  !d2=size(func1d)
  !write(*,*)'rank',mynode,'is in trapz1d'
  delx=(xf-xi)/real(d2-1)
  summ=0.0
  summ=summ+0.5*(func1d(1)+func1d(d2))
  do i=2,d2-1
     summ=summ+func1d(i)
  end do
  summ=summ*delx
end subroutine trapz1d_eq_xy

subroutine trapezoidal3d(func3d,d1,d2,d3,xi,xf,y,zi,zf,summ)
 ! use MainVariables
  implicit none
  integer,intent(in)::d1,d2,d3
  real*8,intent(in)::func3d(d1,d2,d3)
  real*8,intent(in)::y(0:d2)
  real*8,intent(in)::xi,xf,zi,zf
  real*8,intent(out)::summ
  real*8,allocatable,dimension(:,:)::sum2d
  real*8,allocatable,dimension(:)::sum1d
!  real*8::start_time,end_time
  !write(*,*)'rank',mynode,'is in trapezoidal3d'
 ! call cpu_time(start_time)
  allocate(sum2d(d1,d2))
  call trapz3d_eq(func3d,d1,d2,d3,xi,xf,sum2d)
  allocate(sum1d(d1))
  call trapz2d_noneq(sum2d,d1,d2,y,sum1d)
  deallocate(sum2d)
  call trapz1d_eq(sum1d,d1,zi,zf,summ)
  deallocate(sum1d)
end subroutine trapezoidal3d

SUBROUTINE trapezoidal2d_yz(func,D1,D2,zi,zf,y,prosum)
  !USE mpi
  IMPLICIT NONE
  INTEGER,INTENT(IN)::D1,D2
  REAL*8,INTENT(IN),DIMENSION(D1,D2)::func
  REAL*8,INTENT(IN),DIMENSION(0:D2)::y
  REAL*8,ALLOCATABLE,DIMENSION(:)::sum1d
  REAL*8,INTENT(IN)::zi,zf
  REAL*8,INTENT(OUT)::prosum
  ALLOCATE(sum1d(D1))
  CALL trapz2d_noneq(func,D1,D2,y,sum1d)
  CALL trapz1d_eq(sum1d,D1,zi,zf,prosum)
  DEALLOCATE(sum1d)
END SUBROUTINE trapezoidal2d_yz
SUBROUTINE trapezoidal2d_xz(func,D1,D2,zi,zf,xi,xf,prosum)
  !USE mpi
  IMPLICIT NONE
  INTEGER,INTENT(IN)::D1,D2
  REAL*8,INTENT(IN),DIMENSION(D1,D2)::func
  REAL*8,ALLOCATABLE,DIMENSION(:)::sum1d
  REAL*8,INTENT(IN)::zi,zf,xi,xf
  REAL*8,INTENT(OUT)::prosum
  ALLOCATE(sum1d(D1))
  CALL trapz2d_eq(func,D1,D2,xi,xf,sum1d)
  CALL trapz1d_eq(sum1d,D1,zi,zf,prosum)
  DEALLOCATE(sum1d)
END SUBROUTINE trapezoidal2d_xz
 

end module integration


module compfft
  implicit none
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! SUBROUTINE oneDFFT3Dary: this subroutine computes 1D fourier transformation!!
  !!                          when 3D array is supplied. This uses subroutines  !!
  !!                          of FFTPACK5.1 package. The out put is also a 3D   !!
  !!                          array. Transformation takes place only along the  !!
!!                          first dimension of the array.                     !!
!! INPUTS:                                                                    !!
!! l: INTEGER     :  length of the fourier transform                          !!
!! dimen: integer : dimension of the array that need to do fourier transform    !!
!! ary3d : REAL*8        : 3d arry with dimensions D1,D2, and D3.             !!
!!                                                                            !!
!! OUTPUTS:                                                                   !!
!! ftary3d:COMPLEX(KIND=8): 3d array that contain fourier coefficients        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine onedfft3dary(ary3d,l,dimen,ftary3d)
    implicit none

    interface
       subroutine cfft1i ( n, wsave, lensav, ier )
         integer ( kind = 4 ) lensav
         integer ( kind = 4 ) ier
         integer ( kind = 4 ) iw1
         integer ( kind = 4 ) n
         real ( kind = 8 ) wsave(lensav)
       end subroutine cfft1i

       subroutine cfft1f ( n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )
         integer ( kind = 4 ) lenc
         integer ( kind = 4 ) lensav
         integer ( kind = 4 ) lenwrk
         complex ( kind = 8 ) c(lenc)
         integer ( kind = 4 ) ier
         integer ( kind = 4 ) inc
         integer ( kind = 4 ) iw1
         integer ( kind = 4 ) n
         real ( kind = 8 ) work(lenwrk)
         real ( kind = 8 ) wsave(lensav)
       end subroutine cfft1f

    end interface
    integer,intent(in)::l,dimen
    real*8,intent(in),dimension(:,:,:)::ary3d
    complex(kind=8),intent(out),dimension(:,:,:)::ftary3d
    INTEGER :: LENSAV,LENWRK,IER,INC,LENC,I,J,K,d1,d2,d3
    REAL*8,ALLOCATABLE,DIMENSION(:)::WSAVE,WORK
    COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:)::C
  !  write(*,*)'computing 1D fft'
    if (dimen==1)then
   !    write(*,*)'spectra is computed on the first dimension'
       
       d1=size(ary3d,1)
       d2=size(ary3d,2)
       d3=size(ary3d,3)
       INC=1
       LENSAV=2*L+INT(LOG(REAL(L))/LOG(2.0))+4
       LENWRK=2*L
       LENC=INC*(L-1)+1
       ALLOCATE(WORK(LENWRK))
       ALLOCATE(WSAVE(LENSAV))
       ALLOCATE(C(LENC))
       CALL CFFT1I(L,WSAVE,LENSAV,IER)
       IF(L>D1)THEN
          DO k=1,D3
             DO J=1,D2
                DO I=1,L
                   IF(I<=D1)THEN
                      C(I)=CMPLX(ary3d(I,J,K))
                   ELSE
                      C(I)=CMPLX(0.)
                   END IF
                END DO
                CALL CFFT1F(L,INC,C,LENC,WSAVE,LENSAV,WORK,LENWRK,IER)
                DO I=1,L
                   FTary3d(I,J,K)=C(I)
                END DO
             END DO
          END DO
       ELSE
          DO k=1,d3
             DO j=1,d2
                DO I=1,l
                   C(I)=CMPLX(ary3d(I,J,K))
                END DO
                CALL CFFT1F(L,INC,C,LENC,WSAVE,LENSAV,WORK,LENWRK,IER)
                DO I=1,L
                   FTary3d(I,J,K)=C(I)
                END DO
             END DO
          END DO
       END IF
       DEALLOCATE(WORK)
       DEALLOCATE(WSAVE)
       DEALLOCATE(C)
    elseif(dimen==3)then
    !   write(*,*)'spectra is computed on the third dimension'
       
       d1=size(ary3d,1)
       d2=size(ary3d,2)
       d3=size(ary3d,3)
       INC=1
       LENSAV=2*L+INT(LOG(REAL(L))/LOG(2.0))+4
       LENWRK=2*L
       LENC=INC*(L-1)+1
       ALLOCATE(WORK(LENWRK))
       ALLOCATE(WSAVE(LENSAV))
       ALLOCATE(C(LENC))
       CALL CFFT1I(L,WSAVE,LENSAV,IER)
       IF(L>d3)THEN
          DO j=1,d2
             DO i=1,d1
                DO k=1,L
                   IF(k<=d3)THEN
                      C(k)=CMPLX(ary3d(I,J,K))
                   ELSE
                      C(k)=CMPLX(0.)
                   END IF
                END DO
                CALL CFFT1F(L,INC,C,LENC,WSAVE,LENSAV,WORK,LENWRK,IER)
                DO k=1,l
                   FTary3d(I,J,K)=C(k)
                END DO
             END DO
          END DO
       ELSE
          DO j=1,d2
             DO i=1,d1
                DO k=1,l
                   C(k)=CMPLX(ary3d(I,J,K))
                END DO
                CALL CFFT1F(L,INC,C,LENC,WSAVE,LENSAV,WORK,LENWRK,IER)
                DO k=1,l
                   FTary3d(I,J,K)=C(k)
                END DO
             END DO
          END DO
       END IF
       DEALLOCATE(WORK)
       DEALLOCATE(WSAVE)
       DEALLOCATE(C)
    end if
       
  END SUBROUTINE onedfft3dary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE oneDInvFFT3Dary: this subroutine computes 1D fourier            !!
!!                             invers transformation when 3D array of fourier !!
!!                             coefficinets  is supplied.                     !!
!!                             This uses subroutines                          !!
!!                          of FFTPACK5.1 package. The out put is also a 3D   !!
!!                          array. Transformation takes place only along the  !!
!!                          first dimension of the array.                     !!
!! INPUTS:                                                                    !!
!! D1,D2,D3: INTEGER     : dimensions of the 3d array in z,y, and x directions!!
!!                         respectively                                       !!
!! FTary3d : COMPLEX(KIND=8): 3d arry with dimensions D1,D2, and D3.          !!
!!                                                                            !!
!! OUTPUTS:                                                                   !!
!! invFTary3d: REAL*8: 3d array that contain fourier coefficients             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE oneDInvFFT3Dary(ftary3d,l,dimen,invftary3d)
    IMPLICIT NONE
    interface
       subroutine cfft1i ( n, wsave, lensav, ier )
         integer ( kind = 4 ) lensav
         integer ( kind = 4 ) ier
         integer ( kind = 4 ) iw1
         integer ( kind = 4 ) n
         real ( kind = 8 ) wsave(lensav)
       end subroutine cfft1i

       subroutine cfft1b ( n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )
         integer ( kind = 4 ) lenc
         integer ( kind = 4 ) lensav
         integer ( kind = 4 ) lenwrk
         complex ( kind = 8 ) c(lenc)
         integer ( kind = 4 ) ier
         integer ( kind = 4 ) inc
         integer ( kind = 4 ) iw1
         integer ( kind = 4 ) n
         real ( kind = 8 ) work(lenwrk)
         real ( kind = 8 ) wsave(lensav)
       end subroutine cfft1b

    end interface
    integer,intent(IN)::l,dimen
    REAL*8,INTENT(OUT),dimension(:,:,:)::invftary3d
    COMPLEX(KIND=8),INTENT(IN),dimension(:,:,:)::ftary3d
    INTEGER :: LENSAV,LENWRK,IER,INC,LENC,I,J,K,d1,d2,d3
    REAL*8,ALLOCATABLE,DIMENSION(:)::WSAVE,WORK
    COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:)::C
    if (dimen==1)then
       write(*,*)'spectra is computed on the first dimension'
       
       d1=size(ftary3d,1)
       d2=size(ftary3d,2)
       d3=size(ftary3d,3)
       INC=1
       LENSAV=2*l+INT(LOG(REAL(l))/LOG(2.0))+4
       LENWRK=2*l
       LENC=INC*(l-1)+1
       allocate(WORK(LENWRK))
       allocate(WSAVE(LENSAV))
       allocate(C(LENC))
       call CFFT1I(l,WSAVE,LENSAV,IER)
       if (l>d1)then
          do k=1,d3
             do j=1,d2
                do i=1,l
                   if(i<=d1)then
                      c(i)=ftary3d(I,J,K)
                   else
                      c(i)=cmplx(0.0)
                   end if
                end do
                call CFFT1B(l,INC,C,LENC,WSAVE,LENSAV,WORK,LENWRK,IER)
                do i=1,l
                   invftary3d(i,j,k)=c(i)
                end do
             end do
          end do
       else
          do k=1,d3
             do j=1,d2
                do i=1,l
                   c(i)=ftary3d(i,j,k)
                end do
                call CFFT1B(l,INC,C,LENC,WSAVE,LENSAV,WORK,LENWRK,IER)
                do i=1,l
                   invftary3d(i,j,k)=c(i)
                end do
             end do
          end do
       end if
       deallocate(WORK)
       deallocate(WSAVE)
       deallocate(C)
       
    elseif(dimen==3)then
       write(*,*)'spectra is computed on the third dimension'
       
       d1=size(ftary3d,1)
       d2=size(ftary3d,2)
       d3=size(ftary3d,3)
       INC=1
       LENSAV=2*L+INT(LOG(REAL(L))/LOG(2.0))+4
       LENWRK=2*L
       LENC=INC*(L-1)+1
       allocate(WORK(LENWRK))
       allocate(WSAVE(LENSAV))
       allocate(C(LENC))
       call CFFT1I(L,WSAVE,LENSAV,IER)
       if(L>d3)then
          do j=1,d2
             do i=1,d1
                do k=1,l
                   if (k<=d3)then
                      C(k)=ftary3d(i,j,k)
                   else
                      C(k)=cmplx(0.)
                   endif
                END DO
                CALL CFFT1B(L,INC,C,LENC,WSAVE,LENSAV,WORK,LENWRK,IER)
                DO k=1,l
                   invftary3d(I,J,K)=C(k)
                END DO
             END DO
          END DO
       ELSE
          DO j=1,d2
             DO i=1,d1
                DO k=1,l
                   c(k)=ftary3d(i,j,k)
                END DO
                CALL CFFT1B(l,INC,C,LENC,WSAVE,LENSAV,WORK,LENWRK,IER)
                DO k=1,l
                   invftary3d(i,j,k)=c(k)
                END DO
             END DO
          END DO
       END IF
       DEALLOCATE(WORK)
       DEALLOCATE(WSAVE)
       DEALLOCATE(C)
    end if
END SUBROUTINE oneDInvFFT3Dary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE oneDInvFFTShift3Dary: This subroutine shifts the low wavenumber !!
!!                                  coefficients to the middle of the spectrym!!
!! INPUTS:                                                                    !!
!! D1,D2,D3: INTEGER     : dimensions of the 3d array in z,y, and x directions!!
!!                         respectively                                       !!
!! invFTary3d : COMPLEX(KIND=8): the array contains inverse fourier trasformed!!
!!                               coefficients.                                !!
!! OUTPUTS:                                                                   !!
!! shft3dary: REAL*8: array with shifted spectrum                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE oneDinvFFTShift3Dary(invFTary3d,D1,D2,D3,shft3dary)
  IMPLICIT NONE
  INTEGER,INTENT(IN)::D1,D2,D3
  REAL*8,INTENT(IN)::invFTary3d(D1,D2,D3)
  REAL*8,INTENT(OUT)::shft3dary(D1,D2,D3)
  INTEGER::I,J,K
  DO k=1,d3
     DO J=1,D2
        DO i=1,d1
           IF(I<=D1/2)THEN
              shft3dary(I,J,K)=2.0*invFTary3d(I+D1/2,J,K)
           ELSE
              shft3dary(I,J,K)=2.0*invFTary3d(I-D1/2,J,K)
           END IF
        END DO
     END DO
  END DO
END SUBROUTINE oneDinvFFTShift3Dary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE onedfft2dary: this subroutine computes 1D fourier transformation!!
!!                          when 2D array is supplied. This uses subroutines  !!
!!                          of FFTPACK5.1 package. The out put is also a 2D   !!
!!                          array. Transformation takes any dimension of the  !!
!!                          array which is supplied with dim parameter.       !!
!! INPUTS:                                                                    !!
!! l: INTEGER     : length of the array to be fourier transformed with zero   !!
!!                  padding                                                   !!
!! dimen            : the dimension that needs to to the transformation         !!
!!                                                                            !!
!! ary2d : REAL*8 : 2d arry with                                              !!
!!                                                                            !!
!! OUTPUTS:                                                                   !!
!! ftary3d:COMPLEX(KIND=8): 3d array that contain fourier coefficients        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine onedfft2dary(ary2d,l,dimen,ftary2d)
    implicit none
    
    real*8,intent(in),dimension(:,:)::ary2d
    complex(KIND=8),intent(out),dimension(:,:)::ftary2d
    integer,intent(in)::l,dimen
    integer :: lensav,lenwrk,ier,inc,lenc,i,j,k,d1,d2
    real*8,allocatable,dimension(:)::WSAVE,WORK
    complex(kind=8),allocatable,dimension(:)::C
    write(*,*)'computing 1D fft from 2d ary'
    d1=size(ary2d,1)
    d2=size(ary2d,2)
    
    inc=1
    lensav=2*l+int(log(real(l))/log(2.0))+4
    lenwrk=2*l
    lenc=inc*(l-1)+1
    
    allocate(work(lenwrk))
    allocate(wsave(lensav))
    allocate(c(lenc))
    
    call CFFT1I(l,wsave,lensav,ier)
    if (dimen==1)then
       if(l>d1)then
          do j=1,d2
             do i=1,L
                if(i<=D1)then
                   c(i)=cmplx(ary2d(i,j))
                else
                   c(i)=cmplx(0.)
                end if
             end do
             call CFFT1F(l,inc,c,lenc,wsave,lensav,work,lenwrk,ier)
             do i=1,l
                ftary2d(i,j)=c(i)
             end do
          end do
       else
          do j=1,d2
             do i=1,l
                c(i)=cmplx(ary2d(i,j))
             end do
             call CFFT1F(l,inc,c,lenc,wsave,lensav,work,lenwrk,ier)
             do i=1,l
                ftary2d(i,j)=c(i)
             end do
          end do
       end if
    elseif(dimen==2)then
       if(l>d2)then
          do i=1,d1
             do j=1,l
                if(j<=d2)then
                   c(j)=cmplx(ary2d(i,j))
                else
                   c(j)=cmplx(0.)
                end if
             end do
             call CFFT1F(l,inc,c,lenc,wsave,lensav,work,lenwrk,ier)
             do j=1,l
                ftary2d(i,j)=c(j)
             end do
          end do
       else
          do i=1,d1
             do j=1,l
                c(j)=cmplx(ary2d(i,j))
             end do
             call CFFT1F(l,inc,c,lenc,wsave,lensav,work,lenwrk,ier)
             do j=1,l
                ftary2d(i,j)=c(j)
             end do
          end do
       end if
    end if
    deallocate(work)
    deallocate(wsave)
    deallocate(c)
  end subroutine onedfft2dary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE onedfft2dary: this subroutine computes 1D fourier transformation!!
!!                          when 2D array is supplied. This uses subroutines  !!
!!                          of FFTPACK5.1 package. The out put is also a 2D   !!
!!                          array. Transformation takes any dimension of the  !!
!!                          array which is supplied with dim parameter.       !!
!! INPUTS:                                                                    !!
!! l: INTEGER     : length of the array to be fourier transformed with zero   !!
!!                  padding                                                   !!
!! dimen            : the dimension that needs to to the transformation         !!
!!                                                                            !!
!! ary2d : REAL*8 : 2d arry with                                              !!
!!                                                                            !!
!! OUTPUTS:                                                                   !!
!! ftary3d:COMPLEX(KIND=8): 3d array that contain fourier coefficients        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine onedinvfft2dary(ftary2d,dimen,ary2d)
    implicit none
    
    complex(kind=8),intent(in),dimension(:,:)::ftary2d
    real*8,intent(out),dimension(:,:)        ::ary2d
    integer,intent(in)                       ::dimen
    integer                                  ::lensav,lenwrk,ier,inc,lenc,i,j,k,d1,d2,l
    real*8,allocatable,dimension(:)          ::WSAVE,WORK
    complex(kind=8),allocatable,dimension(:) ::C

    write(*,*)'computing 1D fft from 2d ary'
    d1=size(ary2d,1)
    d2=size(ary2d,2)
    if (dimen==1) then
       l=d1
    else
       l=d2
    end if

    inc=1
    lensav=2*l+int(log(real(l))/log(2.0))+4
    lenwrk=2*l
    lenc=inc*(l-1)+1
    
    allocate(work(lenwrk))
    allocate(wsave(lensav))
    allocate(c(lenc))
    
    call CFFT1I(l,wsave,lensav,ier)
    if (dimen==1)then
       do j=1,d2
          do i=1,l
             c(i)=cmplx(ftary2d(i,j))
          end do
          call CFFT1B(l,inc,c,lenc,wsave,lensav,work,lenwrk,ier)
          do i=1,l
             ary2d(i,j)=c(i)
          end do
       end do
    elseif(dimen==2)then
       do i=1,d1
          do j=1,l
             c(j)=cmplx(ary2d(i,j))
          end do
          call CFFT1B(l,inc,c,lenc,wsave,lensav,work,lenwrk,ier)
          do j=1,l
             ary2d(i,j)=c(j)
          end do
       end do
    end if

    deallocate(work)
    deallocate(wsave)
    deallocate(c)
  end subroutine onedinvfft2dary


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE oneDInvFFTShift2Dary: This subroutine shifts the low wavenumber !!
!!                                  coefficients to the middle of the spectrym!!
!! INPUTS:                                                                    !!
!! dimen     : dimension of the 2d array of which the fft needs be found        !!
!! invFTary2d : real*8 : the array contains inverse fourier trasformed        !!
!!                               coefficients.                                !!
!! OUTPUTS:                                                                   !!
!! shft2dary: REAL*8: array with shifted spectrum                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine onedfftshift2Dary(invftary2d,dimen,shft2dary)
  implicit none
  integer,intent(in)               :: dimen
  real*8,intent(in),dimension(:,:) :: invftary2d
  real*8,intent(out),dimension(:,:):: shft2dary
  integer                          :: i,j,d1,d2
  d1=size(invftary2d,1)
  d2=size(invftary2d,2)

  if(dimen==1)then
     do j=1,d2
        do i=1,d1
           if(i<=d1/2)then
              shft2dary(i,j)=2.0*invftary2d(i+d1/2,j)
           else
              shft2dary(i,j)=2.0*invftary2d(i-d1/2,j)
           end if
        end do
     end do
  elseif(dimen==2)then
     do j=1,d2
        do i=1,d1
           if(j<=d2/2)then
              shft2dary(i,j)=2.0*invftary2d(i,j+d2/2)
           else
              shft2dary(i,j)=2.0*invftary2d(i,j-d2/2)
           end if
        end do
     end do
  end if
end subroutine onedfftshift2Dary



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! twoDFFT: This subroutine finds the 2D FFT of any given 3d array. 2D trans- !!
!! formation is taken in along the dimension 1 and the dimension 3 of the     !!
!! 3darray. Output is a complex 3dFTarray. D1, D2, and D3 are dimensions of   !!
!! the input arrays. L1 and L2 are lengths of transformation. Output array has!!
!! the dimensions of L2,N2,L1. If L1 and L2 are greater than N1 and N3        !!
!! respectively, zeros padding is done for the remaining elelments            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE twodfft(ary3d,l1,l2,ftary3d)
  IMPLICIT NONE
                                                   
  integer,intent(in)::L1,L2                             ! length of the transformation
  real*8,intent(in),dimension(:,:,:)::ary3d           ! input array
  complex(KIND=8),intent(out),dimension(:,:,:)::ftary3d ! out put array of cmplex numbers
  integer::i,j,k,d1,d2,d3                               ! these are array indices

  ! following variables are defined for calling FFT subroutines. Description of each
  ! of these variables is given in fftpack.f90 file where all the FFT subroutines are
  ! located. Additional information of this FFTPACK 5.1 package can be found in
  ! http://people.sc.fsu.edu/~jburkardt/f77_src/fftpack5.1/fftpack5.1_reference.html

  INTEGER :: LENSAV,LENWRK,IER,LDIM
  REAL*8,ALLOCATABLE,DIMENSION(:)::WSAVE,WORK
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:)::C

  d1=size(ary3d,1)
  d2=size(ary3d,2)
  d3=size(ary3d,3)
  
  LENSAV=2*(l1+l2)+INT(LOG(REAL(l1))/LOG(2.0))+INT(LOG(REAL(l1))/LOG(2.0))+8
  LENWRK=2*l1*l2
  LDIM=l1
  
  ALLOCATE(WORK(LENWRK))
  ALLOCATE(WSAVE(LENSAV))
  ALLOCATE(C(LDIM,l2))

  ! Initialize 2DFFT
  call CFFT2I(l1,l2,WSAVE,LENSAV,IER)
  if (l1>d1.and.l2>d3)then
     do j=1,d2
        do k=1,l2
           do i=1,l1
              if(i<=d1.and.k<=d3)then
                 C(i,k)=cmplx(ary3d(i,j,k))
              else
                 C(i,k)=cmplx(0.)
              end if
           end do
        end do
        call CFFT2F(LDIM,L1,L2,C,WSAVE,LENSAV,WORK,LENWRK,IER)
        do K=1,L2
           do I = 1,L1
              ftary3d(i,j,k)=C(i,k)
           end do
        END DO
     END DO
  ELSE
     DO j=1,d2
        DO k=1,l2
           DO i=1,l1
              C(i,k)=CMPLX(ary3d(i,j,k))
           END DO
        END DO
        CALL CFFT2F(LDIM,l1,l2,C,WSAVE,LENSAV,WORK,LENWRK,IER)
        DO k=1,l2
           DO i= 1,l1
              ftary3d(i,j,k)=C(i,k)
           END DO
        END DO
     END DO
  END IF
  DEALLOCATE(WORK)
  DEALLOCATE(WSAVE)
  DEALLOCATE(C)
END SUBROUTINE twoDFFT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE twoDInvFFT: this subroutine computes the 2D inverse fourier     !!
!!                        transformation using FFTPCK5.1 pacakage             !!
!! INPUTS:                                                                    !!
!! D1,D2,D3: INTEGER          : dimensions of the 3d array in z,y, and x dir  !!
!!                              respectively                                  !!
!! ary3dIn : COMPLEX(KIND=8)  : 3d arry with dimensions D1,D2, and D3. This   !!
!!                              contains fourier coefficients                 !!
!!                                                                            !!
!! OUTPUTS:                                                                   !!
!! ary3dinvT:COMPLEX(KIND=8)  : 3d arry with dimensions D1,D2, and D3. inverse!!
!!                              transformed array                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE twoDInvFFT(ary3dIn,ary3dinvT)
  IMPLICIT NONE
  
  COMPLEX(KIND=8),INTENT(IN),dimension(:,:,:)::ary3DIn
  COMPLEX(KIND=8),INTENT(OUT),dimension(:,:,:)::ary3DinvT
  INTEGER::D1,D2,D3
  INTEGER :: LENSAV,LENWRK,IER,LDIM,I,J,K,L1,L2
  REAL*8,ALLOCATABLE,DIMENSION(:)::WSAVE,WORK
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:)::C
  COMPLEX(kind=8),ALLOCATABLE,DIMENSION(:,:,:,:,:)::Rhat_wrk
  
  d1=size(ary3din,1)
  d2=size(ary3din,2)
  d3=size(ary3din,3)

  L1=d1 !x-dir has taken as the dir 1
  L2=d3 !z-dir has take and the dir 2
  LENSAV=2*(L1+L2)+INT(LOG(REAL(L1))/LOG(2.0))+INT(LOG(REAL(L1))/LOG(2.0))+8
  LENWRK=2*L1*L2
  LDIM=L1
  
  ALLOCATE(WORK(LENWRK))
  ALLOCATE(WSAVE(LENSAV))
  ALLOCATE(C(LDIM,L2))

  CALL CFFT2I(L1,L2,WSAVE,LENSAV,IER)
  !write(*,*)'calc inverse fft'
  do J=1,D2
     do k = 1,l2
        do i = 1,l1
           C(i,k)=ary3dIn(i,j,k)
        end do
     end do
     CALL CFFT2B(LDIM,l1,l2,C,WSAVE,LENSAV,WORK,LENWRK,IER)
     do k=1,l2
        do i=1,l1
           ary3dinvT(i,j,k)=4.0*C(i,k)
        end do
     end do
  end do
  DEALLOCATE(work,wsave,c)
END SUBROUTINE twoDInvFFT



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE invFFTshift: This takes zero wavenumber component to the middle !!
!!                         of the array.                                      !!
!! INPUTS:                                                                    !!
!! D1,D2,D3: INTEGER          : dimensions of the 3d array in z,y, and x dir  !!
!!                              respectively                                  !!
!! ary3din : COMPLEX(KIND=8)  : 3d arry with dimensions D1,D2, and D3. This   !!
!!                              contains fourier coefficients                 !!
!!                                                                            !!
!! OUTPUTS:                                                                   !!
!! shftd3dAry:REAL*8 : 3d arry with dimensions D1,D2, and D3. Shifted 3d array!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE invFFTshift(ary3din,shftd3dAry)
  IMPLICIT NONE
  
  COMPLEX(KIND=8),INTENT(IN),dimension(:,:,:)::ary3din
  REAL*8,INTENT(OUT),dimension(:,:,:)::shftd3dAry
  REAL*8,DIMENSION(:,:,:),ALLOCATABLE::work
  INTEGER::I,J,K,d1,d2,d3
  d1=size(ary3din,1)
  d2=size(ary3din,2)
  d3=size(ary3din,3)

  ALLOCATE(work(D1,D2,D3))
  DO K=1,D3
     DO J=1,D2
        DO I=1,D1/2
           work(I,J,K)=REAL(ary3din(I+D1/2,J,K))
           work(I+D1/2,J,K)=REAL(ary3din(I,J,K))
        END DO
     END DO
  END DO
  DO K=1,D3/2
     DO J=1,D2
        DO I=1,D1
           shftd3dAry(I,J,K)=work(I,J,K+D3/2)
           shftd3dAry(I,J,K+D3/2)=work(I,J,K)
        END DO
     END DO
  END DO
  DEALLOCATE(work)
END SUBROUTINE invFFTshift

end module compfft
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module pdf
  implicit none

contains
  
  subroutine count3dfunc(varin,checkval,probcount)
    implicit none
    real*8, intent(in),dimension(:,:,:)::varin
    real*8, intent(in)                 ::checkval
    integer, intent(inout),dimension(:,:,:)::probcount
    integer                            ::d1,d2,d3,i,j,k

    d1=size(varin,1)
    d2=size(varin,2)
    d3=size(varin,3)    
    
    do k=1,d3
       do j=1,d2
          do i=1,d1
             if(varin(i,j,k)<0.and.varin(i,j,k)>=checkval)then
                probcount(i,j,k)=probcount(i,j,k)+1.0
             end if
          end do
       end do
    end do
    
  end subroutine count3dfunc
end module pdf
    




