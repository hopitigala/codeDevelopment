!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This module contains subroutines which are used to compute POD            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module comppod
  implicit none

contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!subroutine fluc_xyplane: extract fluctuating field from 2D plane!!!!!!!
  !!!                                                                       !!!! 
  !!! INPUTS:                                                               !!!!
  !!!       fluc3d: fluctuation field array 3d                              !!!!
  !!!       x0    : the coordinate at which the data is extracted           !!!!
  !!!       plane : this specifies the 2D plane (xy =1, yz=2,xz =3)         !!!!
  !!!       tstp  : the time step at which the data is taken                !!!!
  !!! OUTPUT:                                                               !!!!
  !!!       fluct: array with 2D plane data for each time step.            !!!!
  !!!               the 3rd dimension holds the time step data.             !!!!
  !!!                                                                       !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fluc_2Dplane(fluc3d,x0,plane,tstp,fluct)
    implicit none
    integer,intent(in)::x0,plane,tstp
    real*8,intent(in),dimension(:,:,:)::fluc3d
    real*8,intent(out),dimension(:,:,:)::fluct
    integer :: i,j,k,d1,d2,d3
    
    d1=size(fluc3d,1)
    d2=size(fluc3d,2)
    d3=size(fluc3d,3)
    ! xy=1 ,yz=2, xz=3
    if(plane==1)then
       do k=1,d3
          do j=1,d2
             fluct(j,k,tstp)=fluc3d(x0,j,k)
          end do
       end do
    elseif(plane==2) then
       do j=1,d2
          do i=1,d1             
             fluct(i,j,tstp)=fluc3d(i,j,x0)
          end do
       end do
    else
       do k=1,d3
          do i=1,d1             
             fluct(i,k,tstp)=fluc3d(i,x0,k)
          end do
       end do
    end if
  end subroutine fluc_2Dplane

  subroutine fluc_3d(fluc3d,tstp,fluct)
    implicit none
    integer,intent(in)                   :: tstp
    real*8,intent(in),dimension(:,:,:)   :: fluc3d
    real*8,intent(out),dimension(:,:,:,:):: fluct
    integer                              :: i,j,k,d1,d2,d3
    
    d1=size(fluc3d,1)
    d2=size(fluc3d,2)
    d3=size(fluc3d,3)
   
    do k=1,d3
       do j=1,d2
          do i=1,d1
             fluct(i,j,k,tstp)=fluc3d(i,j,k)
          end do
       end do
    end do
  end subroutine fluc_3d


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!! This subroutine computes the kernal matrix that is used to find !!!!!!
 !!!!!!! eigenvalues and eigen functions.                                !!!!!!
 !!!!!!!                                                                 !!!!!!
 !!!!!!! INPUTS:                                                         !!!!!!
 !!!!!!!      u,v,w = 3d arrays of fluctuating velocity components in    !!!!!!
 !!!!!!!              in 3 directions.                                   !!!!!!
 !!!!!!!      x,y,z = coordinate arrays in 3 directions                  !!!!!!
 !!!!!!!      iplane= integer that defines the plane in which the data   !!!!!!
 !!!!!!!              is obtained. (xy=1, yz=2, xz=3)                    !!!!!!
 !!!!!!!                                                                 !!!!!!
 !!!!!!! OUTPUT:                                                         !!!!!!
 !!!!!!!          c = kernal matrix                                      !!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine compkernal(u,v,w,x,y,z,iplane,c)
!    use integration
    use mpi
    use writedata
    use mpivariables
    implicit none
    real*8,intent(in),dimension(:,:,:):: u,v,w
    real*8,intent(in),dimension(:)    :: x,y,z
    integer,intent(in)                :: iplane
    real*8,intent(out),dimension(:,:) :: c
    
    integer                           :: d1,d2,d3,tm,tn,alocst
    real*8,allocatable,dimension(:,:) :: u2sum
    real*8                            :: prosum,inv

    d1=size(u,1)
    d2=size(u,2)
    d3=size(u,3)
    write(*,*)'d3=',d3
    do tm=1,d3
       do tn=1,d3
          allocate(u2sum(d1,d2),stat=alocst)
          if(alocst/=0)stop "****** NOT ENOUGH MEMORY *******"
          
          u2sum(:,:)=u(:,:,tm)*u(:,:,tn)+v(:,:,tm)*v(:,:,tn)+w(:,:,tm)*w(:,:,tn)

          if(iplane==1)then
             call trapezoidal2d_xy(u2sum,d1,d2,x(1),x(d2),y,prosum)
             
          elseif(iplane==2)then
             call trapezoidal2d_yz(u2sum,d1,d2,z(1),z(d1),y,prosum)
             
          else
             call trapezoidal2d_xz(u2sum,d1,d2,z(1),z(d1),x(mynode*d2+1),x(d2*(mynode+1)),prosum)
          end if

          deallocate(u2sum)

          call MPI_ALLREDUCE(prosum,c(tm,tn),1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
          
       end do
    end do
          
    inv=1.0/real(d3)
    c=inv*c
    ! print kernal matrix
    if (mynode==0)then
       if(iplane==1)then
          call print2DMatx(c,'kernalMatx_xy.dat')
       elseif (iplane==2)then
          call print2DMatx(c,'kernalMatx_yz.dat')
       else
          call print2DMatx(c,'kernalMatx_xz.dat')
       end if
    end if
  end subroutine compkernal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!! This subroutine computes the kernal matrix that is used to find !!!!!!
 !!!!!!! eigenvalues and eigen functions.                                !!!!!!
 !!!!!!!                                                                 !!!!!!
 !!!!!!! INPUTS:                                                         !!!!!!
 !!!!!!!      u,v,w = 4d arrays of fluctuating velocity components in    !!!!!!
 !!!!!!!              in 3 directions.                                   !!!!!!
 !!!!!!!      x,y,z = coordinate arrays in 3 directions                  !!!!!!
 !!!!!!!                                                                 !!!!!!
 !!!!!!! OUTPUT:                                                         !!!!!!
 !!!!!!!          c = kernal matrix                                      !!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine compkernal_3d(u,v,w,x,y,z,c)
!    use integration
    use mpi
    use writedata
    use mpivariables
    implicit none
    real*8,intent(in),dimension(:,:,:,:)  :: u,v,w
    real*8,intent(in),dimension(:)        :: x,y,z
    real*8,intent(out),dimension(:,:)     :: c
    
    integer                               :: d1,d2,d3,d4,tm,tn,alocst
    real*8,allocatable,dimension(:,:,:)   :: u2sum
    real*8                                :: prosum,inv

    d1=size(u,1)
    d2=size(u,2)
    d3=size(u,3)
    d4=size(u,4)
    write(*,*)'d4=',d4
    do tm=1,d4
       do tn=1,d4
          allocate(u2sum(d1,d2,d3),stat=alocst)
          if(alocst/=0)stop "****** NOT ENOUGH MEMORY *******"
          
          u2sum=u(:,:,:,tm)*u(:,:,:,tn)+v(:,:,:,tm)*v(:,:,:,tn)+w(:,:,:,tm)*w(:,:,:,tn)
          
          call trapezoidal3d(u2sum,d1,d2,d3,x(1),x(d3),y,z(1),z(d1),prosum)     

          deallocate(u2sum)

          call MPI_ALLREDUCE(prosum,c(tm,tn),1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
          
       end do
    end do
          
    inv=1.0/real(d4)
    c=inv*c
    ! print kernal matrix
    if (mynode==0)then
       call print2DMatx(c,'kernalMatx.dat')
    end if
  end subroutine compkernal_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine computes POD modes.                                       !
!!!       1. First the 3D velocity array is converted to a 2D array (reshape).!
!!!       2. Then use matmul function to do matrix mutiplication with evec    !
!!!       3. finally 2D array obtained from multiplication is reshaped to a   !
!!!          3D array.                                                        !
!!!                                                                           !
!!! INPUTS:                                                                   !
!!!      u : instantaneous velocity fluctuation array                         ! 
!!!    evec: eigenvectors array                                               !
!!!                                                                           !
!!! OUTPUT:                                                                   !
!!!     phi : 3D array that contain eigenfuctions or POD modes pertaining to  !
!!!           the given velocity field.                                       !
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine podModes_2Dplane(u,evec,phi)
    implicit none
    real*8,intent(in),dimension(:,:,:) :: u
    real*8,intent(in),dimension(:,:)   :: evec
    real*8,intent(out),dimension(:,:,:):: phi
    
    integer                            :: d1,d2,d3
    
    d1=size(u,1)
    d2=size(u,2)
    d3=size(u,3)
    

    phi=reshape(matmul(reshape(u,(/d1*d2,d3/)),evec),(/d1,d2,d3/))
  end subroutine podModes_2Dplane

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine computes POD modes from 3d velocity field.                !
!!!       1. First the 4D velocity array is converted to a 2D array (reshape).!
!!!       2. Then use matmul function to do matrix mutiplication with evec    !
!!!       3. finally 2D array obtained from multiplication is reshaped to a   !
!!!          4D array.                                                        !
!!!                                                                           !
!!! INPUTS:                                                                   !
!!!      u : instantaneous velocity fluctuation array for all time steps      ! 
!!!    evec: eigenvectors array                                               !
!!!                                                                           !
!!! OUTPUT:                                                                   !
!!!     phi : 4D array that contain eigenfuctions or POD modes pertaining to  !
!!!           the given velocity field and for each mode.                      !
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine podmodes_3d(u,evec,phi)
    implicit none
    real*8,intent(in),dimension(:,:,:,:) :: u
    real*8,intent(in),dimension(:,:)     :: evec
    real*8,intent(out),dimension(:,:,:,:):: phi
    
    integer                              :: d1,d2,d3,d4
    
    d1=size(u,1)
    d2=size(u,2)
    d3=size(u,3)
    d4=size(u,4)

    phi=reshape(matmul(reshape(u,(/d1*d2*d3,d4/)),evec),(/d1,d2,d3,d4/))
  end subroutine podmodes_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine normalize POD modes                                        !
!!!                                                                            !
!!! INPUTS:                                                                    !
!!!   x,y,z    : coordinate arrays (1D) for 3 directions (need for integration)!
!!! phi,shi,chi: 3D arrays of POD modes for u,v, and w velocity fields         !
!!!              (not normalized)                                              !
!!! iplane     : (=1,xy; =2, yz; =3,xz planes)                                 !
!!!                                                                            !
!!! OUTPUT:                                                                    !
!!! phi,shi,chi: 3D arrays of POD modes for u,v, and w velocity fields         !
!!!              (normalized)                                                  !
!!!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine norm_modes_2Dplane(x,y,z,iplane,phi,shi,chi)
    use mpi
    use mpivariables
    use writedata
    use integration
    implicit none

    real*8,intent(inout),dimension(:,:,:)::phi,shi,chi
    real*8,intent(in),dimension(:)       ::x,y,z
    integer,intent(in)                   ::iplane
    real*8,allocatable,dimension(:,:)    ::phisum
    real*8,allocatable,dimension(:)      ::norm2

    integer                              ::d1,d2,d3,alocst,m
    real*8                               ::prosum,norm

    d1=size(phi,1)
    d2=size(phi,2)
    d3=size(phi,3)

    allocate(norm2(d3),stat=alocst)
    if(alocst/=0)stop "****** NOT ENOUGH MEMORY *******"

    do m=1,d3
       allocate(phisum(d1,d2),stat=alocst)
       if(alocst/=0)stop "****** NOT ENOUGH MEMORY *******"

       phisum(:,:)=phi(:,:,m)*phi(:,:,m)+shi(:,:,m)*shi(:,:,m)+chi(:,:,m)*chi(:,:,m)

       if(iplane==1)then
          call trapezoidal2d_xy(phisum,d1,d2,x(1),x(d2),y,prosum)

       elseif(iplane==2)then
          call trapezoidal2d_yz(phisum,d1,d2,z(1),z(d1),y,prosum)

       else
          call trapezoidal2d_xz(phisum,d1,d2,z(1),z(d1),x(mynode*d2+1),x(d2*(mynode+1)),prosum)
       end if

       deallocate(phisum)
call MPI_ALLREDUCE(prosum,norm2(m),1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

       norm=1.0/sqrt(norm2(m))
       phi(:,:,m)=phi(:,:,m)*norm
       shi(:,:,m)=shi(:,:,m)*norm
       chi(:,:,m)=chi(:,:,m)*norm

    end do

    !printing norm2

    if (mynode==0)then
       if (iplane==1)then
          call printVector(norm2,'normSquare_xyPlane.dat')
       elseif(iplane==2)then
          call printVector(norm2,'normSquare_yzPlane.dat')
       else
          call printVector(norm2,'normSquare_xzPlane.dat')
       end if
    end if
    deallocate(norm2)

  end subroutine norm_modes_2Dplane

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine normalize POD modes                                        !
!!!                                                                            !
!!! INPUTS:                                                                    !
!!!   x,y,z    : coordinate arrays (1D) for 3 directions (need for integration)!
!!! phi,shi,chi: 4D arrays of POD modes for u,v, and w velocity fields         !
!!!              (not normalized)                                              !
!!! OUTPUT:                                                                    !
!!! phi,shi,chi: 4D arrays of POD modes for u,v, and w velocity fields         !
!!!              (normalized)                                                  !
!!!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine norm_modes_3d(x,y,z,phi,shi,chi)
    use mpi
    use mpivariables
    use writedata
    use integration
    implicit none

    real*8,intent(inout),dimension(:,:,:,:):: phi,shi,chi
    real*8,intent(in),dimension(:)         :: x,y,z
    real*8,allocatable,dimension(:,:,:)    :: phisum
    real*8,allocatable,dimension(:)        :: norm2
   
    integer                                :: d1,d2,d3,d4,alocst,m
    real*8                                 :: prosum,norm
    
    d1=size(phi,1)
    d2=size(phi,2)
    d3=size(phi,3)
    d4=size(phi,4)
    allocate(norm2(d4),stat=alocst)
    if(alocst/=0)stop "****** NOT ENOUGH MEMORY *******"

    do m=1,d4
       allocate(phisum(d1,d2,d3),stat=alocst)
       if(alocst/=0)stop "****** NOT ENOUGH MEMORY *******"
          
       phisum=phi(:,:,:,m)*phi(:,:,:,m)+shi(:,:,:,m)*shi(:,:,:,m)+chi(:,:,:,m)*chi(:,:,:,m)
      
       call trapezoidal3d(phisum,d1,d2,d3,x(1),x(d3),y,z(1),z(d1),prosum)
       deallocate(phisum)

       call MPI_ALLREDUCE(prosum,norm2(m),1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       
       norm=1.0/sqrt(norm2(m))
       phi(:,:,:,m)=phi(:,:,:,m)*norm
       shi(:,:,:,m)=shi(:,:,:,m)*norm
       chi(:,:,:,m)=chi(:,:,:,m)*norm
       
    end do

    !printing norm2

    if (mynode==0)then
       call printVector(norm2,'normSquare.dat')
    end if
    deallocate(norm2)

  end subroutine norm_modes_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine checks whether POD modes are orthonormal                   !
!!! We do not necessarily need this in calculations                            !
!!! INPUTS:                                                                    !
!!!   x,y,z    : coordinate arrays (1D) for 3 directions (need for integration)!
!!! phi,shi,chi: 3D arrays of POD modes for u,v, and w velocity fields         !
!!!              (normalized)                                                  !
!!! iplane     : (1=xy, 2=yz, 3=xz)                                            !
!!!                                                                            !
!!!                                                                            !
!!! OUTPUT:                                                                    !
!!! no outputs creates a file that has magnitudes of POD modes mutliplied by   !
!!! each other                                                                 !
!!!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine checkModeOrthog_2Dpln(x,y,z,iplane,phi,shi,chi)
    use mpivariables
    use integration
    use writedata
    use mpi
  
    implicit none
  
    real*8,intent(in),dimension(:,:,:) :: phi,shi,chi
    real*8,intent(in),dimension(:)     :: x,y,z
    integer,intent(in)                 :: iplane
    
    real*8,dimension(:,:),allocatable::phisum
    real*8,dimension(:,:),allocatable::norm2
    REAL*8::prosum
    integer::d1,d2,d3,m,n,allocStat

    d1=size(phi,1)
    d2=size(phi,2)
    d3=size(phi,3)

    ! Finding the L2 norm of each mode
    allocate(norm2(d3,d3))
    do m=1,d3
       do n=1,d3
          allocate(phisum(d1,d2),STAT=allocStat)
          if(allocStat/=0)STOP "*****NOT ENOUGH MEMORY*****"
          phisum=phi(:,:,m)*phi(:,:,n)+shi(:,:,m)*shi(:,:,n)+chi(:,:,m)*chi(:,:,n)
         
          if(iplane==1)then
             call trapezoidal2d_xy(phisum,d1,d2,x(1),x(d2),y,prosum)
             
          elseif(iplane==2)then
             call trapezoidal2d_yz(phisum,d1,d2,z(1),z(d1),y,prosum)
             
          else
             call trapezoidal2d_xz(phisum,d1,d2,z(1),z(d1),x(mynode*d2+1),x(d2*(mynode+1)),prosum)
          end if
        
          deallocate(phisum)
          
          call MPI_ALLREDUCE(prosum,norm2(m,n),1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       end do
    end do
    if(mynode==0)then
       if(iplane==1)then
          call print2DMatx(norm2,'checkorthog_xy.dat')
       elseif(iplane==2)THEN
          call print2DMatx(norm2,'checkorthog_yz.dat')
       else
          call print2DMatx(norm2,'checkorthog_xz.dat')
       end if
    end if
    deallocate(norm2)
  end subroutine CheckModeOrthog_2Dpln

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine checks whether POD modes are orthonormal                   !
!!! We do not necessarily need this in calculations                            !
!!! INPUTS:                                                                    !
!!!   x,y,z    : coordinate arrays (1D) for 3 directions (need for integration)!
!!! phi,shi,chi: 4D arrays of POD modes for u,v, and w velocity fields         !
!!!              (normalized)                                                  !
!!!                                                                            !
!!!                                                                            !
!!!                                                                            !
!!! OUTPUT:                                                                    !
!!! no outputs creates a file that has magnitudes of POD modes mutliplied by   !
!!! each other                                                                 !
!!!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine checkModeOrthog_3d(x,y,z,phi,shi,chi)
    use mpivariables
    use integration
    use writedata
    use mpi
  
    implicit none
  
    real*8,intent(in),dimension(:,:,:,:) :: phi,shi,chi
    real*8,intent(in),dimension(:)       :: x,y,z
    
    real*8,dimension(:,:,:),allocatable::phisum
    real*8,dimension(:,:),allocatable::norm2
    REAL*8::prosum
    integer::d1,d2,d3,d4,m,n,allocStat

    d1=size(phi,1)
    d2=size(phi,2)
    d3=size(phi,3)
    d4=size(phi,4)

    ! Finding the L2 norm of each mode
    allocate(norm2(d4,d4))
    do m=1,d4
       do n=1,d4
          allocate(phisum(d1,d2,d3),STAT=allocStat)
          if(allocStat/=0)STOP "*****NOT ENOUGH MEMORY*****"
          phisum=phi(:,:,:,m)*phi(:,:,:,n)+shi(:,:,:,m)*shi(:,:,:,n)+chi(:,:,:,m)*chi(:,:,:,n)
          
          call trapezoidal3d(phisum,d1,d2,d3,x(1),x(d3),y,z(1),z(d1),prosum)
        
          deallocate(phisum)
          
          call MPI_ALLREDUCE(prosum,norm2(m,n),1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       end do
    end do
    if(mynode==0)then
       call print2DMatx(norm2,'checkorthog.dat')
    end if
    deallocate(norm2)
  end subroutine CheckModeOrthog_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine computes coefficients for POD modes                        !
!!!                                                                            !
!!! INPUTS:                                                                    !
!!!  u,v,w            : velocit field arrays for all the time steps            !
!!!                     (3rd dimension is for the time)                        !
!!!                                                                            !
!!! modeu,modev,modew : POD modes for u,v, and w fields                        !
!!!                     (3rd dimension is for the number of modes)             !
!!!                                                                            !
!!!   x,y,z           : coordinate arrays (1D) for 3 directions                !
!!!                     (need for integration)                                 !
!!!                                                                            !
!!! iplane            : (1=xy, 2=yz, 3=xz)                                     !
!!!                                                                            !
!!!                                                                            !
!!! OUTPUT:                                                                    !
!!!                                                                            !
!!! coe               : coefficient array 2D                                   !
!!!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine podcoe_2d(x,y,z,u,v,w,modeu,modev,modew,iplane,coe)
    use mpivariables
    use mpi
    use writedata
    use integration
    implicit none
    
    real*8,intent(in),dimension(:,:,:):: u,v,w,modeu,modev,modew
    real*8,intent(in),dimension(:)    :: x,y,z
    integer,intent(in)                :: iplane
    real*8,intent(out),dimension(:,:) :: coe

    real*8,dimension(:,:),allocatable :: modesum
    real*8                            :: prosum
    integer                           :: m,t,alocst,d1,d2,d3

    d1=size(modeu,1)
    d2=size(modeu,2)
    d3=size(modeu,3)

    do t=1,d3
       do m=1,d3
          allocate(modesum(d1,d2),stat=alocst)
          if(alocst/=0)stop"****** NOT ENOUGH MEMORY ********"
          modesum=modeu(:,:,m)*u(:,:,t)+modev(:,:,m)*v(:,:,t)+modew(:,:,m)*w(:,:,t)
         
          if(iplane==1)then
             call trapezoidal2d_xy(modesum,d1,d2,x(1),x(d2),y,prosum)
             
          elseif(iplane==2)then
             call trapezoidal2d_yz(modesum,d1,d2,z(1),z(d1),y,prosum)
             
          else
             call trapezoidal2d_xz(modesum,d1,d2,z(1),z(d1),x(mynode*d2+1),x(d2*(mynode+1)),prosum)
          end if
        
          deallocate(modesum)
          
          call MPI_ALLREDUCE(prosum,coe(t,m),1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       
       end do
    end do
    
    if(mynode==0)then
       if(iplane==1)then
          call print2DMatx(coe,'coeMatx_xy.dat')
       elseif(iplane==2)THEN
          call print2DMatx(coe,'coeMatx_yz.dat')
       else
          call print2DMatx(coe,'coeMatx_xz.dat')
       end if
    end if
  end subroutine podcoe_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine computes coefficients for POD modes                        !
!!!                                                                            !
!!! INPUTS:                                                                    !
!!!  u,v,w            : velocit field arrays for all the time steps            !
!!!                     (4th dimension is for the time)                        !
!!!                                                                            !
!!! modeu,modev,modew : POD modes for u,v, and w fields                        !
!!!                     (4th dimension is for the number of modes)             !
!!!                                                                            !
!!!   x,y,z           : coordinate arrays (1D) for 3 directions                !
!!!                     (need for integration)                                 !
!!!                                                                            !
!!!                                                                            !
!!!                                                                            !
!!! OUTPUT:                                                                    !
!!!                                                                            !
!!! coe               : coefficient array 2D                                   !
!!!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine podcoe_3d(x,y,z,u,v,w,modeu,modev,modew,coe)
    use mpivariables
    use mpi
    use writedata
    use integration
    implicit none
    
    real*8,intent(in),dimension(:,:,:,:):: u,v,w,modeu,modev,modew
    real*8,intent(in),dimension(:)      :: x,y,z
    real*8,intent(out),dimension(:,:)   :: coe

    real*8,dimension(:,:,:),allocatable :: modesum
    real*8                              :: prosum
    integer                             :: m,t,alocst,d1,d2,d3,d4

    d1=size(modeu,1)
    d2=size(modeu,2)
    d3=size(modeu,3)
    d4=size(modeu,4)
    
    do t=1,d4
       do m=1,d4
          allocate(modesum(d1,d2,d3),stat=alocst)
          if(alocst/=0)stop"****** NOT ENOUGH MEMORY ********"
          modesum=modeu(:,:,:,m)*u(:,:,:,t)+modev(:,:,:,m)*v(:,:,:,t)+modew(:,:,:,m)*w(:,:,:,t)
          
          call trapezoidal3d(modesum,d1,d2,d3,x(1),x(d3),y,z(1),z(d1),prosum)
          
          deallocate(modesum)
          
          call MPI_ALLREDUCE(prosum,coe(t,m),1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       
       end do
    end do
    
    if(mynode==0)then
       call print2DMatx(coe,'coeMatx.dat')
    end if
  end subroutine podcoe_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                                    !!!
!!! This subroutine computes the reconstructed velocity field by using given number of !!!
!!! modes.                                                                             !!!
!!! INPUTS:                                                                            !!!
!!!   mode : 3d array contains the POD mode for the given velocity field               !!!
!!!          3rd dimension represents number of modes                                  !!!
!!!   coe  : coefficient matrix (2d)                                                   !!!
!!!                                                                                    !!! 
!!! stmode : integer number that describes the mode number that starts the mode addition!!
!!! endmode: integer number that describes the mode number that ends the mode addition !!!
!!! tstp   :  time step. reconstructed fiels is generated for each time step           !!!
!!!                                                                                    !!!
!!!  OUTPUT:                                                                           !!!
!!!     recnstfld: reconstructed velocity field                                        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                    
  subroutine recnstfield(mode,coe,stmode,endmode,ntstp,recnstfld)
    implicit none
    real*8,dimension(:,:,:),intent(in) ::mode
    real*8,dimension(:,:),intent(in)   ::coe
    integer,intent(in)                 ::stmode,endmode,ntstp
    real*8,dimension(:,:,:),intent(out)::recnstfld
    
    !real*8,dimension(:,:,:),allocatable::temp
    integer::t,m,d1,d2,d3

    d1=size(mode,1)
    d2=size(mode,2)
    d3=size(mode,3)
    do t=1,ntstp
       recnstfld(:,:,t)=0.0
       do m=stmode,endmode
          recnstfld(:,:,t)=recnstfld(:,:,t)+coe(t,m)*mode(:,:,m)
       end do
    end do
  end subroutine recnstfield

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                                    !!!
!!! This subroutine computes the reconstructed velocity field by using given number of !!!
!!! modes.                                                                             !!!
!!! INPUTS:                                                                            !!!
!!!   mode : 3d array contains the POD mode for the given velocity field               !!!
!!!          3rd dimension represents number of modes                                  !!!
!!!   coe  : coefficient matrix (2d)                                                   !!!
!!!                                                                                    !!! 
!!! stmode : integer number that describes the mode number that starts the mode addition!!
!!! endmode: integer number that describes the mode number that ends the mode addition !!!
!!! tstp   :  time step. reconstructed fiels is generated for each time step           !!!
!!!                                                                                    !!!
!!!  OUTPUT:                                                                           !!!
!!!     recnstfld: reconstructed velocity field                                        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                    

  subroutine recnstfield_3d(mode,coe,stmode,endmode,tstp,recnstfld)
    implicit none
    real*8,dimension(:,:,:,:),intent(in)   ::mode
    real*8,dimension(:,:),intent(in)       ::coe
    integer,intent(in)                     ::stmode,endmode,tstp
    real*8,dimension(:,:,:,:),intent(out)  ::recnstfld
    
    !real*8,dimension(:,:,:),allocatable::temp
    integer::t,m,d1,d2,d3,d4

    d1=size(mode,1)
    d2=size(mode,2)
    d3=size(mode,3)
    d4=size(mode,4)
    
    recnstfld(:,:,:,tstp)=0.0
    
    do m=stmode,endmode
       recnstfld(:,:,:,tstp)=recnstfld(:,:,:,tstp)+coe(tstp,m)*mode(:,:,:,m)
    end do
  end subroutine recnstfield_3d
  
  subroutine reyst_ssm_lsm(u,v,ibs,iplane,uv1,uv2)
    use readdata
    implicit none
    integer,intent(in)                :: ibs,iplane
    real*8,intent(in),dimension(:,:,:):: u,v
    real*8,intent(out),dimension(:)   :: uv1
    real*8,intent(out),dimension(:,:) :: uv2
    integer                           :: d1,d2,d3,t,i,k
    real*8                            :: inv,inv1,inv2
    real*8,allocatable,dimension(:,:) :: temp
    
    d1=size(u,1)
    d2=size(u,2)
    d3=size(u,3)
    
    allocate(temp(d1,d2))
    
    temp(:,:)=0.0
    do t=1,d3
       temp(:,:)=temp(:,:)+u(:,:,t)*v(:,:,t)
    end do
    inv=1/real(d3)
    
    temp=inv*temp
    
    if(ibs==0)then
       uv2=0.0
       if(iplane==1)then !xy plane then average over x
          uv1=0.0
          do k=1,d2
             uv1(:)=uv1(:)+temp(:,k)
          end do
          inv1=1.0/real(d2)
          uv1=inv1*uv1
       elseif(iplane==2)then !yz plane then average over z
          uv1=0.0
          uv2=0.0
          do i=1,d1
             uv1(:)=uv1(:)+temp(i,:)
          end do
          inv1=1.0/real(d1)
          uv1=inv1*uv1
       end if
    else !when blowing is activated
       !if(iplane==2)then
       uv1=0.0
       uv2=temp
        !  do i=28,164,34
        !     uv1(:)=uv1(:)+temp(i,:) ! averages taken at the centerline of the blowing holes
        !  end do
        !  inv2=1.0/5.0
        !  uv1=uv1*inv2
        !  do i=45,147,34
        !     uv2(:)=uv2(:)+temp(i,:) ! averages taken at the midline between two blowing holes
        !  end do
        !  uv2=uv2*inv2
       !else
       !   write(*,*)'not an appropriate value for iplane'
       !end if
    end if
    deallocate(temp)
  end subroutine reyst_ssm_lsm

  subroutine reyst_ssm_lsm_3d(u,v,ibs,uv1,uv2)

    implicit none
    integer,intent(in)                   :: ibs
    real*8,intent(in),dimension(:,:,:,:) :: u,v
    real*8,intent(out),dimension(:)      :: uv1
    real*8,intent(out),dimension(:,:,:)  :: uv2
    integer                              :: d1,d2,d3,d4,t,i,k
    real*8                               :: inv,inv1,inv2
    real*8,allocatable,dimension(:,:,:)  :: temp
    real*8,allocatable,dimension(:,:)    :: uv1_x
    
    d1=size(u,1)
    d2=size(u,2)
    d3=size(u,3)
    d4=size(u,4)
    
    allocate(temp(d1,d2,d3))
    
    temp=0.0
    do t=1,d3
       temp=temp+u(:,:,:,t)*v(:,:,:,t)
    end do
    inv=1/real(d4)
    
    temp=inv*temp
    
    if(ibs==0)then
       allocate(uv1_x(d1,d2))
       ! averaging over x direction
       uv1_x=0.0
       do k=1,d3
          uv1_x=uv1_x+temp(:,:,k)
       end do
       inv1=1.0/real(d3)
       uv1_x=inv1*uv1_x
       ! avreaging over z direction
       uv1=0.0
       do i=1,d1
          uv1=uv1+uv1_x(i,:)
       end do
       inv1=1.0/real(d1)
       uv1=inv1*uv1
     
    else 
       uv2=temp
       
    end if
    
  end subroutine reyst_ssm_lsm_3d
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!                                                                                             !!!!!!
  !!!!! This subroutine computes scalar flux profiles along the wall-normal direction               !!!!!!
  !!!!! u - velocit field, t- scalar field. ibs- (=1 blowing, =0 unperturbed),                      !!!!!!
  !!!!! iplane (=1 xy, =2 yz =3 xz planes)                                                          !!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine sclflx_ssm_lsm(u,scl,ibs,iplane,uv1,uv2)
    use readdata
    implicit none
    integer,intent(in)                :: ibs,iplane
    real*8,intent(in),dimension(:,:,:):: u,scl
    real*8,intent(out),dimension(:)   :: uv1
    real*8,intent(out),dimension(:,:) :: uv2
    integer                           :: d1,d2,d3,t,i,k
    real*8                            :: inv,inv1,inv2
    real*8,allocatable,dimension(:,:) :: temp
    
    d1=size(u,1)
    d2=size(u,2)
    d3=size(u,3)
    
    allocate(temp(d1,d2))
    
    temp(:,:)=0.0
    do t=1,d3
       temp(:,:)=temp(:,:)+u(:,:,t)*scl(:,:,t)
    end do
    inv=1/real(d3)
    
    temp=inv*temp
    
    if(ibs==0)then
       uv2=0.0
       if(iplane==1)then !xy plane then average over x
          uv1=0.0
          do k=1,d2
             uv1(:)=uv1(:)+temp(:,k)
          end do
          inv1=1.0/real(d2)
          uv1=inv1*uv1
       elseif(iplane==2)then !yz plane then average over z
          uv1=0.0
          do i=1,d1
             uv1(:)=uv1(:)+temp(i,:)
          end do
          inv1=1.0/real(d1)
          uv1=inv1*uv1
       end if
    else !no spatial averaging for blowing
       uv1=0.0
       uv2=temp
    end if
    deallocate(temp)
  end subroutine sclflx_ssm_lsm

  subroutine sclflx_ssm_lsm_3d(u,scl,ibs,ut1,ut2)
    use readdata
    implicit none
    integer,intent(in)                   :: ibs
    real*8,intent(in),dimension(:,:,:,:) :: u,scl
    real*8,intent(out),dimension(:)      :: ut1
    real*8,intent(out),dimension(:,:,:)  :: ut2
    integer                              :: d1,d2,d3,d4,t,i,k
    real*8                               :: inv,inv1,inv2
    real*8,allocatable,dimension(:,:,:)  :: temp
    real*8,allocatable,dimension(:,:)    :: ut1_x
    
    d1=size(u,1)
    d2=size(u,2)
    d3=size(u,3)
    d4=size(u,4)
    
    allocate(temp(d1,d2,d3))
    
    temp=0.0
    do t=1,d3
       temp=temp+u(:,:,:,t)*scl(:,:,:,t)
    end do
    inv=1/real(d4)
    
    temp=inv*temp
    
    if(ibs==0)then
       allocate(ut1_x(d1,d2))
       ! averaging over x direction
       ut1_x=0.0
       do k=1,d3
          ut1_x=ut1_x+temp(:,:,k)
       end do
       inv1=1.0/real(d3)
       ut1_x=inv1*ut1_x
       ! avreaging over z direction
       ut1=0.0
       do i=1,d1
          ut1=ut1+ut1_x(i,:)
       end do
       inv1=1.0/real(d1)
       ut1=inv1*ut1
       
    else 
       ut2=temp
    
    end if
    
  end subroutine sclflx_ssm_lsm_3d


  subroutine valatcnstx0y0(ary2d,j0,k0,iplane,node,refval)
    use mpi
    use mpivariables
    use channhdvariables
    use prelimcalvar
    implicit none
    real*8,intent(in),dimension(:,:)  :: ary2d
    integer,intent(in)                :: j0,iplane,k0
    real*8, intent(out)               :: refval
    integer,intent(out)               :: node
    integer                           :: yprime

    node=j0/n2do
    if(mynode==node)then
       yprime=j0-node*n2do
       if(iplane==1)then !xy plane k0 is x coordinate
          refval=ary2d(yprime,k0)
       elseif(iplane==2)then !yz plane k0 is z coordinate
          refval=ary2d(k0,yprime)
       end if
    end if
  end subroutine valatcnstx0y0


  subroutine tpcorl_2d_podfield(var1,var2,numtimesteps,i0,j0,k0,iplane)
    use readdata
    use MainParameters
    use channhdvariables
    use mpi
    use mpivariables
    use writedata
    use prelimcalvar
    implicit none

    integer,intent(in)                 ::var1,var2,i0,j0,k0,iplane
    !real*8,intent(in),dimension(0:)    ::ysdo,ypdo
    integer                            ::d1,d2,d3,t,node_l,node_s,node_t
    integer                            ::rmsnode_l,rmsnode_s,rmsnode_t,numtimesteps
    real*8,allocatable,dimension(:,:,:)::v1_lsm,v1_ssm,v2_lsm,v2_ssm
    real*8,allocatable,dimension(:,:)  ::rmsv1_l,rmsv1_s,rmsv2_l,rmsv2_s
    real*8,allocatable,dimension(:,:)  ::rmsv1_t,rmsv2_t,covar_t,covartavg_t
    real*8,allocatable,dimension(:,:)  ::covar_l,covar_s,covartavg_l,covartavg_s
    real*8,allocatable,dimension(:,:)  ::inv_l,inv_s,inv_t 
    real*8                             ::refval_l,refval_s,inv,refval_t
    real*8                             ::rmsrefval_l,rmsrefval_s,rmsrefval_t
    character(len=4)                   ::xprm,zprm
    character(len=3)                   ::yprm
    character(len=2)                   ::prosnum
    character(len=1)                   ::vari1,vari2

     write(prosnum,'(i2.2)')mynode


    if(iplane==1) then!xy plane
       allocate(v1_lsm(n2do+1,n3,numtimesteps+1),v1_ssm(n2do+1,n3,numtimesteps+1))
       if (var1==3)then
          call read3Darray('u_lsm_xy'//trim(prosnum),'unformatted',v1_lsm)
          call read3Darray('u_ssm_xy'//trim(prosnum),'unformatted',v1_ssm)
       elseif (var1==2)then
          call read3Darray('v_lsm_xy'//trim(prosnum),'unformatted',v1_lsm)
          call read3Darray('v_ssm_xy'//trim(prosnum),'unformatted',v1_ssm)
       end if
    elseif(iplane==2)then!yzplane
       allocate(v1_lsm(n1,n2do+1,numtimesteps+1),v1_ssm(n1,n2do+1,numtimesteps+1))
       if (var1==3)then
          call read3Darray('u_lsm_yz'//trim(prosnum),'unformatted',v1_lsm)
          call read3Darray('u_ssm_yz'//trim(prosnum),'unformatted',v1_ssm)
       elseif (var1==2)then
          call read3Darray('v_lsm_yz'//trim(prosnum),'unformatted',v1_lsm)
          call read3Darray('v_ssm_yz'//trim(prosnum),'unformatted',v1_ssm)
       end if
    end if

    if(iplane==1) then!xy plane
       allocate(v2_lsm(n2do+1,n3,numtimesteps+1),v2_ssm(n2do+1,n3,numtimesteps+1))
       if (var2==3)then
          call read3Darray('u_lsm_xy'//trim(prosnum),'unformatted',v2_lsm)
          call read3Darray('u_ssm_xy'//trim(prosnum),'unformatted',v2_ssm)
       elseif (var2==2)then
          call read3Darray('v_lsm_xy'//trim(prosnum),'unformatted',v2_lsm)
          call read3Darray('v_ssm_xy'//trim(prosnum),'unformatted',v2_ssm)
       end if
    elseif(iplane==2)then!yzplane
       allocate(v2_lsm(n1,n2do+1,numtimesteps+1),v2_ssm(n1,n2do+1,numtimesteps+1))
       if (var2==3)then
          call read3Darray('u_lsm_yz'//trim(prosnum),'unformatted',v2_lsm)
          call read3Darray('u_ssm_yz'//trim(prosnum),'unformatted',v2_ssm)
       elseif (var2==2)then
          call read3Darray('v_lsm_yz'//trim(prosnum),'unformatted',v2_lsm)
          call read3Darray('v_ssm_yz'//trim(prosnum),'unformatted',v2_ssm)
       end if
    end if

    d1=size(v1_lsm,1)
    d2=size(v1_lsm,2)
    d3=size(v1_lsm,3)
 
    allocate(covartavg_l(d1,d2),covartavg_s(d1,d2),covartavg_t(d1,d2))
    allocate(rmsv1_l(d1,d2),rmsv1_s(d1,d2),rmsv2_l(d1,d2),rmsv2_s(d1,d2))
    allocate(rmsv1_t(d1,d2),rmsv2_t(d1,d2))

    rmsv1_l=0.
    rmsv2_l=0.
    rmsv1_s=0.
    rmsv2_s=0.
    rmsv1_t=0.
    rmsv2_t=0.

    covartavg_l=0.
    covartavg_s=0.
    covartavg_t=0.

    do t=1,d3

       call valatcnstx0y0(v2_lsm(:,:,t),j0,k0,iplane,node_l,refval_l)
       call valatcnstx0y0(v2_ssm(:,:,t),j0,k0,iplane,node_s,refval_s)
       call valatcnstx0y0(v2_lsm(:,:,t)+v2_ssm(:,:,t),j0,k0,iplane,node_t,refval_t)

       call mpi_bcast(refval_l,1,mpi_real8,node_l,mpi_comm_world,ierr)
       call mpi_bcast(refval_s,2,mpi_real8,node_s,mpi_comm_world,ierr)
       call mpi_bcast(refval_t,3,mpi_real8,node_t,mpi_comm_world,ierr)

       allocate(covar_l(d1,d2),covar_s(d1,d2),covar_t(d1,d2))
       covar_l=v1_lsm(:,:,t)*refval_l
       covar_s=v1_ssm(:,:,t)*refval_s
       covar_t=(v1_lsm(:,:,t)+v1_ssm(:,:,t))*refval_t
          

       covartavg_l=covartavg_l+covar_l
       covartavg_s=covartavg_s+covar_s
       covartavg_t=covartavg_t+covar_t

       deallocate(covar_l,covar_s,covar_t)
       rmsv1_l=rmsv1_l+v1_lsm(:,:,t)*v1_lsm(:,:,t)
       rmsv1_s=rmsv1_s+v1_ssm(:,:,t)*v1_ssm(:,:,t)
       rmsv1_t=rmsv1_t+(v1_lsm(:,:,t)+v1_ssm(:,:,t))*(v1_lsm(:,:,t)+v1_ssm(:,:,t))
       
       rmsv2_l=rmsv2_l+v2_lsm(:,:,t)*v2_lsm(:,:,t)
       rmsv2_s=rmsv2_s+v2_ssm(:,:,t)*v2_ssm(:,:,t)
       rmsv2_t=rmsv2_t+(v2_lsm(:,:,t)+v2_ssm(:,:,t))*(v2_lsm(:,:,t)+v2_ssm(:,:,t))
    end do

    inv=1.0/real(d3)

    
    rmsv1_l=sqrt(rmsv1_l*inv)
    rmsv2_l=sqrt(rmsv2_l*inv)
    rmsv1_s=sqrt(rmsv1_s*inv)
    rmsv2_s=sqrt(rmsv2_s*inv)
    rmsv1_t=sqrt(rmsv1_t*inv)
    rmsv2_t=sqrt(rmsv2_t*inv)

    call valatcnstx0y0(rmsv2_l,j0,k0,iplane,rmsnode_l,rmsrefval_l)
    call valatcnstx0y0(rmsv2_s,j0,k0,iplane,rmsnode_s,rmsrefval_s)
    call valatcnstx0y0(rmsv2_t,j0,k0,iplane,rmsnode_t,rmsrefval_t)

    call mpi_bcast(rmsrefval_l,1,mpi_real8,rmsnode_l,mpi_comm_world,ierr)
    call mpi_bcast(rmsrefval_s,2,mpi_real8,rmsnode_s,mpi_comm_world,ierr)
    call mpi_bcast(rmsrefval_t,3,mpi_real8,rmsnode_s,mpi_comm_world,ierr)

    allocate(inv_l(d1,d2),inv_s(d1,d2),inv_t(d1,d2))
    inv_l=1.0/(rmsv1_l*rmsrefval_l)
    inv_s=1.0/(rmsv1_s*rmsrefval_s)
    inv_t=1.0/(rmsv1_t*rmsrefval_t)


    covartavg_l=covartavg_l*inv*inv_l
    covartavg_s=covartavg_s*inv*inv_s
    covartavg_t=covartavg_t*inv*inv_t

    write(vari1,'(i1.1)')var1
    write(vari2,'(i1.1)')var2


    if(j0<10)then
       write(yprm,'(i1.1)')j0
    elseif(j0<100)then
       write(yprm,'(i2.2)')j0
    else
       write(yprm,'(i3.3)')j0
    end if

    if(k0<10)then
       write(xprm,'(i1.1)')k0
    elseif(k0<100)then
       write(xprm,'(i2.2)')k0
    elseif(k0<1000)then
       write(xprm,'(i3.3)')k0
    else
       write(xprm,'(i4.4)')k0
    end if

    if(i0<10)then
        write(zprm,'(i1.1)')i0
     elseif(i0<100)then
        write(zprm,'(i2.2)')i0
     elseif(i0<1000)then
        write(zprm,'(i3.3)')i0
     else
        write(zprm,'(i4.4)')i0
     end if

    if(iplane==1)then
       call sendrecv2dwrite(covartavg_l,1,iplane,'corcoe_xy_lsm'//trim(vari1)//'_'//trim(vari2)//'_z0_'//trim(zprm)//'_y0_'//trim(yprm)//'_x0_'//trim(xprm)//'.dat')
       call sendrecv2dwrite(covartavg_s,2,iplane,'corcoe_xy_ssm'//trim(vari1)//'_'//trim(vari2)//'_z0_'//trim(zprm)//'_y0_'//trim(yprm)//'_x0_'//trim(xprm)//'.dat')
       call sendrecv2dwrite(covartavg_t,3,iplane,'corcoe_xy_tot'//trim(vari1)//'_'//trim(vari2)//'_z0_'//trim(zprm)//'_y0_'//trim(yprm)//'_x0_'//trim(xprm)//'.dat')
    elseif(iplane==2)then
       call sendrecv2dwrite(covartavg_l,1,iplane,'corcoe_yz_lsm'//trim(vari1)//'_'//trim(vari2)//'_z0_'//trim(xprm)//'_y0_'//trim(yprm)//'_x0_'//trim(zprm)//'.dat')
       call sendrecv2dwrite(covartavg_s,2,iplane,'corcoe_yz_ssm'//trim(vari1)//'_'//trim(vari2)//'_z0_'//trim(xprm)//'_y0_'//trim(yprm)//'_x0_'//trim(zprm)//'.dat')
       call sendrecv2dwrite(covartavg_t,3,iplane,'corcoe_yz_tot'//trim(vari1)//'_'//trim(vari2)//'_z0_'//trim(xprm)//'_y0_'//trim(yprm)//'_x0_'//trim(zprm)//'.dat')
    end if
  end subroutine tpcorl_2d_podfield

  subroutine valatcnstz0y0x0(ary3d,i0,j0,k0,node,refval)
    use mpi
    use mpivariables
    use channhdvariables
    use prelimcalvar
    implicit none
    real*8,intent(in),dimension(:,:,:):: ary3d
    integer,intent(in)                :: j0,i0,k0
    real*8, intent(out)               :: refval
    integer,intent(out)               :: node
    integer                           :: yprime

    node=j0/n2do
    if(mynode==node)then
       yprime=j0-node*n2do
       refval=ary3d(i0,yprime,k0)
    end if
  end subroutine valatcnstz0y0x0

  subroutine tpcorl_3d_podfield(var1,var2,timestcount,i0,j0,k0)
    use readdata
    use MainParameters
    use channhdvariables
    use mpi
    use mpivariables
    use writedata
    use prelimcalvar
    implicit none

    integer,intent(in)                   ::var1,var2,i0,j0,k0,timestcount
    integer                              ::d1,d2,d3,d4,t,node_l,node_s,node_t
    integer                              ::rmsnode_l,rmsnode_s,rmsnode_t
    real*8,allocatable,dimension(:,:,:,:)::v1_lsm,v1_ssm,v2_lsm,v2_ssm
    real*8,allocatable,dimension(:,:,:)  ::rmsv1_l,rmsv1_s,rmsv2_l,rmsv2_s
    real*8,allocatable,dimension(:,:,:)  ::rmsv1_t,rmsv2_t
    real*8,allocatable,dimension(:,:,:)  ::covar_l,covar_s,covartavg_l,covartavg_s
    real*8,allocatable,dimension(:,:,:)  ::covar_t,covartavg_t
    real*8,allocatable,dimension(:,:,:)  ::inv_l,inv_s,inv_t
    real*8                               ::refval_l,refval_s,inv,refval_t
    real*8                               ::rmsrefval_l,rmsrefval_s,rmsrefval_t
    character(len=4)                     ::xprm,zprm
    character(len=3)                     ::yprm
    character(len=2)                     ::prosnum
    character(len=1)                     ::vari1,vari2

    write(prosnum,'(i2.2)')mynode

    allocate(v1_lsm(n1,n2do+1,n3,timestcount),v1_ssm(n1,n2do+1,n3,timestcount))
    if (var1==3)then
       call read4Darray('u_lsm'//trim(prosnum),'unformatted',v1_lsm)
       call read4Darray('u_ssm'//trim(prosnum),'unformatted',v1_ssm)
    elseif (var1==2)then
       call read4Darray('v_lsm'//trim(prosnum),'unformatted',v1_lsm)
       call read4Darray('v_ssm'//trim(prosnum),'unformatted',v1_ssm)
    elseif (var1==1)then
       call read4Darray('w_lsm'//trim(prosnum),'unformatted',v1_lsm)
       call read4Darray('w_ssm'//trim(prosnum),'unformatted',v1_ssm)
    end if
     
    allocate(v2_lsm(n1,n2do+1,n3,timestcount),v2_ssm(n1,n2do+1,n3,timestcount))
    if (var2==3)then
       call read4Darray('u_lsm'//trim(prosnum),'unformatted',v2_lsm)
       call read4Darray('u_ssm'//trim(prosnum),'unformatted',v2_ssm)
    elseif (var2==2)then
       call read4Darray('v_lsm'//trim(prosnum),'unformatted',v2_lsm)
       call read4Darray('v_ssm'//trim(prosnum),'unformatted',v2_ssm)
    elseif (var2==1)then
       call read4Darray('w_lsm'//trim(prosnum),'unformatted',v2_lsm)
       call read4Darray('w_ssm'//trim(prosnum),'unformatted',v2_ssm)
    end if

    d1=size(v1_lsm,1)
    d2=size(v1_lsm,2)
    d3=size(v1_lsm,3)
    d4=size(v1_lsm,4)

    allocate(covartavg_l(d1,d2,d3),covartavg_s(d1,d2,d3),covartavg_t(d1,d2,d3))
    allocate(rmsv1_l(d1,d2,d3),rmsv1_s(d1,d2,d3),rmsv2_l(d1,d2,d3),rmsv2_s(d1,d2,d3))
    allocate(rmsv1_t(d1,d2,d3),rmsv2_t(d1,d2,d3))

    rmsv1_l=0.
    rmsv2_l=0.
    rmsv1_s=0.
    rmsv2_s=0.
    rmsv1_t=0.
    rmsv2_t=0.

    covartavg_l=0.
    covartavg_s=0.
    covartavg_t=0.

    do t=1,d4

       call valatcnstz0y0x0(v2_lsm(:,:,:,t),i0,j0,k0,node_l,refval_l)
       call valatcnstz0y0x0(v2_ssm(:,:,:,t),i0,j0,k0,node_s,refval_s)
       call valatcnstz0y0x0(v2_lsm(:,:,:,t)+v2_ssm(:,:,:,t),i0,j0,k0,node_t,refval_t)

       call mpi_bcast(refval_l,1,mpi_real8,node_l,mpi_comm_world,ierr)
       call mpi_bcast(refval_s,1,mpi_real8,node_s,mpi_comm_world,ierr)
       call mpi_bcast(refval_t,1,mpi_real8,node_s,mpi_comm_world,ierr)

       allocate(covar_l(d1,d2,d3),covar_s(d1,d2,d3),covar_t(d1,d2,d3))
       covar_l=v1_lsm(:,:,:,t)*refval_l
       covar_s=v1_ssm(:,:,:,t)*refval_s
       covar_t=(v1_lsm(:,:,:,t)+v1_ssm(:,:,:,t))*refval_t

       covartavg_l=covartavg_l+covar_l
       covartavg_s=covartavg_s+covar_s
       covartavg_t=covartavg_t+covar_t

       deallocate(covar_l,covar_s,covar_t)

       rmsv1_l=rmsv1_l+v1_lsm(:,:,:,t)*v1_lsm(:,:,:,t)
       rmsv1_s=rmsv1_s+v1_ssm(:,:,:,t)*v1_ssm(:,:,:,t)
       rmsv1_t=rmsv1_t+(v1_lsm(:,:,:,t)+v1_ssm(:,:,:,t))*(v1_lsm(:,:,:,t)+v1_ssm(:,:,:,t))

       rmsv2_l=rmsv2_l+v2_lsm(:,:,:,t)*v2_lsm(:,:,:,t)
       rmsv2_s=rmsv2_s+v2_ssm(:,:,:,t)*v2_ssm(:,:,:,t)
       rmsv2_t=rmsv2_t+(v2_lsm(:,:,:,t)+v2_ssm(:,:,:,t))*(v2_lsm(:,:,:,t)+v2_ssm(:,:,:,t))
    end do

    inv=1.0/real(d4)

    
    rmsv1_l=sqrt(rmsv1_l*inv)
    rmsv2_l=sqrt(rmsv2_l*inv)
    rmsv1_s=sqrt(rmsv1_s*inv)
    rmsv2_s=sqrt(rmsv2_s*inv)
    rmsv1_t=sqrt(rmsv1_t*inv)
    rmsv2_t=sqrt(rmsv2_t*inv)

    call valatcnstz0y0x0(rmsv2_l,i0,j0,k0,rmsnode_l,rmsrefval_l)
    call valatcnstz0y0x0(rmsv2_s,i0,j0,k0,rmsnode_s,rmsrefval_s)
    call valatcnstz0y0x0(rmsv2_t,i0,j0,k0,rmsnode_t,rmsrefval_t)

    call mpi_bcast(rmsrefval_l,1,mpi_real8,rmsnode_l,mpi_comm_world,ierr)
    call mpi_bcast(rmsrefval_s,1,mpi_real8,rmsnode_s,mpi_comm_world,ierr)
    call mpi_bcast(rmsrefval_t,1,mpi_real8,rmsnode_t,mpi_comm_world,ierr)

    allocate(inv_l(d1,d2,d3),inv_s(d1,d2,d3),inv_t(d1,d2,d3))
    inv_l=1.0/(rmsv1_l*rmsrefval_l)
    inv_s=1.0/(rmsv1_s*rmsrefval_s)
    inv_t=1.0/(rmsv1_t*rmsrefval_t)

    !inv=1.0/real(d3)

    covartavg_l=covartavg_l*inv*inv_l
    covartavg_s=covartavg_s*inv*inv_s
    covartavg_t=covartavg_t*inv*inv_t
    
    write(vari1,'(i1.1)')var1
    write(vari2,'(i1.1)')var2


    if(j0<10)then
       write(yprm,'(i1.1)')j0
    elseif(j0<100)then
       write(yprm,'(i2.2)')j0
    else
       write(yprm,'(i3.3)')j0
    end if

    if(k0<10)then
       write(xprm,'(i1.1)')k0
    elseif(k0<100)then
       write(xprm,'(i2.2)')k0
    elseif(k0<1000)then
       write(xprm,'(i3.3)')k0
    else
       write(xprm,'(i4.4)')k0
    end if

    if(i0<10)then
        write(zprm,'(i1.1)')i0
     elseif(i0<100)then
        write(zprm,'(i2.2)')i0
     elseif(i0<1000)then
        write(zprm,'(i3.3)')i0
     end if

     call sendrecv3dwrite(covartavg_l,1,'corcoe_lsm'//trim(vari1)//'_'//trim(vari2)//'_z0_'//trim(zprm)//'_y0_'//trim(yprm)//'_x0_'//trim(xprm)//'.dat')
     call sendrecv3dwrite(covartavg_s,2,'corcoe_ssm'//trim(vari1)//'_'//trim(vari2)//'_z0_'//trim(zprm)//'_y0_'//trim(yprm)//'_x0_'//trim(xprm)//'.dat')
     call sendrecv3dwrite(covartavg_t,3,'corcoe_tot'//trim(vari1)//'_'//trim(vari2)//'_z0_'//trim(zprm)//'_y0_'//trim(yprm)//'_x0_'//trim(xprm)//'.dat')
   end subroutine tpcorl_3d_podfield


end module comppod
