  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file contains routines required for 3d POD computation                                     !!!
module comppod3d
  implicit none
contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!nfils : number of fiels per time step when data is written
!filstopros: number of processors that one file devides to 

  subroutine prelimcal_3dpod(nfils)
    use mainparameters
    use channhdvariables
    use mpivariables
    use prelimcalvar
    implicit none
    integer,intent(in)::nfils
    !integer,intent(out)::filstopros
    n1m=n1-1
    n2m=n2-1
    n3m=n3-1
    n2do=n2m/nfils !the number of files
    filstopros=numprocs/nfils ! number of processor that one file devides to
    n3do=n3m/filstopros
    lx1=lx1d*pi
    lx2=lx2d
    lx3=lx3d*pi
  end subroutine prelimcal_3dpod

subroutine intpolflddata_3dpod(ut,vt,wt,pt,tt,ysdo,ypdo,up,vp,wp,pp,tp)
  implicit none
  real*8,intent(in),dimension(:,0:,:)   :: ut,vt,wt,pt,tt
  real*8,intent(in),dimension(0:)       :: ysdo,ypdo
  real*8,intent(out),dimension(:,:,:)   :: up,vp,wp,pp,tp
  integer                               :: i,j,k,d1,d2,d3
  
  d1=size(ut,1)
  d2=size(ut,2)
  d3=size(ut,3)

  do k=1,d3
     do j=1,d2-1
        do i=1,d1
           wp(i,j,k)=wt(i,j-1,k)+((wt(i,j,k)-wt(i,j-1,k))/(ysdo(j)-ysdo(j-1)))*(ypdo(j)-ysdo(j-1))
           up(i,j,k)=ut(i,j-1,k)+((ut(i,j,k)-ut(i,j-1,k))/(ysdo(j)-ysdo(j-1)))*(ypdo(j)-ysdo(j-1))
           vp(i,j,k)=vt(i,j,k)
           tp(i,j,k)=tt(i,j,k)
           pp(i,j,k)=pt(i,j,k)
        end do
     end do
  end do

end subroutine intpolflddata_3dpod

subroutine intpolcordi_3dpod(yp,xp,ys,ypdo,ysdo,xpdo)
  use prelimcalvar
  use mpivariables
  implicit none
  !integer,intent(in)::mynode,numprocs
  real*8,intent(in),DIMENSION(0:)::yp
  real*8,intent(in),DIMENSION(:)::xp
  real*8, INTENT(OUT),DIMENSION(0:)::ys,ysdo,ypdo
  real*8,intent(out),DIMENSION(:)::xpdo
  integer:: j,jj,d2,d2do,fn,m,k,kk
  character*2                          :: prosnum
  character(*),parameter               :: fmt2='(i2.2)'
    ! number of grid points in y-dir
  d2=size(yp)-1
!  fn=int((mynode-mod(mynode,filstopros))/filstopros)!processor number of the original file
  m=mod(mynode,filstopros)
  !d2do= (d2-1)/numprocs
  ! linear interpolation to find coordinates at the
  do j=1,d2-1
     ys(j)=(yp(j)+yp(j+1))/2.
  end do

  ys(0)=2.*yp(1)-ys(1)
  ys(d2)=2.*yp(d2)-ys(d2-1)
  
  do j=0,n2do+1
     jj=fn*n2do+J
     if(j==0.and.fn==0)then
        ypdo(j)=yp(1)
     end if
     ypdo(j)=yp(jj)
     ysdo(j)=ys(jj)
  end do

  do k=1,n3do
     kk=m*n3do+k
     xpdo(k)=xp(kk)
  end do
!  write(prosnum,fmt2)mynode
!  open(22,file='k_index'//trim(prosnum)//'.dat')
!  do k=1,n3do
!     kk=m*n3do+k
!     write(22,'(2x,i2,2x,i1,2x,i4,2x,i3)')mynode,m,k,kk
!  end do
!  close(22)
end subroutine intpolcordi_3dpod

subroutine readpostdata_3dpod(sttime,numtimesteps,dt,icorlavg,icen,nfils)
  implicit none
  real*8,intent(out)::sttime
  integer,intent(out)::numtimesteps,dt,nfils,icorlavg,icen
  open(11,file='post.d')
  read(11,*)sttime,numtimesteps,dt,nfils
  read(11,*)icorlavg,icen
  close(11)
end subroutine readpostdata_3dpod


subroutine readpod3ddata(icomppod,irecnst,lsmlim,ireyst,iscfx,inst,icorl,var1,var2,i0,j0,k0,i1toM)
  implicit none
  integer,intent(out)::icomppod,irecnst,lsmlim,ireyst,inst,iscfx
  integer,intent(out)::icorl,var1,var2,i0,j0,k0
  integer,intent(out)::i1toM ! use only in anisotripy code for pod 
  open(12,file='pod3ddata.dat')
  read(12,*)icomppod,irecnst,lsmlim,ireyst,inst,iscfx
  read(12,*)icorl,var1,var2,i0,j0,k0
  read(12,*)i1toM
  close(12)
end subroutine readpod3ddata


subroutine readtemp3dpod(mynode,itime,ypdo,ysdo,up,vp,wp,pp,tp)
  use channhdvariables
  use prelimcalvar
  use mpi
  implicit none
  integer,intent(in)                   :: mynode,itime
  real*8,intent(in),dimension(0:)      :: ypdo,ysdo
  real*8,intent(out),dimension(:,:,:)  :: up,vp,wp,pp,tp
  real*8,allocatable,dimension(:,:,:,:):: utemp,ttemp
  real*8,allocatable,dimension(:,:,:)  :: ptemp,utemp_pros,vtemp_pros,wtemp_pros,ttemp_pros,ptemp_pros
  real*8,allocatable,dimension(:,:,:)  :: ut,vt,wt,tt,pt
  integer                              :: i,j,k,kk,l,d1,d2,d3,allocst,m,fn,d3do,ierr
  character*2                          :: prosnum
  character*5                          :: pntime
  character*100                        :: filename
  character(*),parameter               :: fmt1='(i5.5)'
  character(*),parameter               :: fmt2='(i2.2)'
  
  fn=int((mynode-mod(mynode,filstopros))/filstopros)!processor number of the original file
  m=mod(mynode,filstopros)
 
  write(pntime,fmt1)itime
  write(prosnum,fmt2)fn
  !open(23,file='mynode'//trim(prosnum)//'.dat')
 ! write(*,'(2x,i2,2x,i2,2x,i,2x,i3,2x,i4)')mynode,fn,m,n3do,itime
  
  allocate(utemp(n1m,0:n2do+1,n3m,3),ttemp(n1m,0:n2do+1,n3m,1),&
       ptemp(n1m,0:n2do+1,n3m),stat=allocst)
  if(allocst /=0) stop "******Not enough memory*********"
  
  filename='../../../../field/field'//trim(prosnum)//'/field.data'//trim(prosnum)//'_'//trim(pntime)
  
  d1=size(utemp,1)
  d2=size(utemp,2)-1
  d3=size(utemp,3)
  
  open(12,file=filename,form='unformatted')
  read(12)
  read(12)
  if(ipassc==1)then
     read(12)(((utemp(i,j,k,1),i=1,d1),j=0,d2),k=1,d3),&
          (((utemp(i,j,k,2),i=1,d1),j=0,d2),k=1,d3),&
          (((utemp(i,j,k,3),i=1,d1),j=0,d2),k=1,d3),&
          (((ptemp(i,j,k),i=1,d1),j=0,d2),k=1,d3),&
          ((((ttemp(i,j,k,1),i=1,d1),j=0,d2),k=1,d3),l=1,npsc)
  else
     read(12)(((utemp(i,j,k,1),i=1,d1),j=0,d2),k=1,d3),&
          (((utemp(i,j,k,2),i=1,d1),j=0,d2),k=1,d3),&
          (((utemp(i,j,k,3),i=1,d1),j=0,d2),k=1,d3),&
          (((ptemp(i,j,k),i=1,d1),j=0,d2),k=1,d3)
  end if
  close(12)
 ! write(*,*)'reading from original done tstep = ',itime,mynode
 ! call mpi_barrier(MPI_COMM_WORLD,ierr)
  allocate(ut(n1m,0:n2do+1,n3m),vt(n1m,0:n2do+1,n3m),&
       wt(n1m,0:n2do+1,n3m),tt(n1m,0:n2do+1,n3m),&
       pt(n1m,0:n2do+1,n3m),stat=allocst)
  if(allocst /=0) stop "******Not enough memory*********"
  if (ipassc==1)then
     wt=utemp(:,:,:,1)
     vt=utemp(:,:,:,2)
     ut=utemp(:,:,:,3)
     pt=ptemp
     tt=ttemp(:,:,:,1)
  else
     wt=utemp(:,:,:,1)
     vt=utemp(:,:,:,2)
     ut=utemp(:,:,:,3)
     pt=ptemp
  end if
  deallocate(utemp,ttemp,ptemp)
 ! write(*,*)'copying to 3d array is done = ',itime,mynode
 ! call mpi_barrier(MPI_COMM_WORLD,ierr)
  allocate(utemp_pros(n1m,0:n2do+1,n3do),vtemp_pros(n1m,0:n2do+1,n3do),&
       wtemp_pros(n1m,0:n2do+1,n3do),ttemp_pros(n1m,0:n2do+1,n3do),&
       ptemp_pros(n1m,0:n2do+1,n3do),stat=allocst)
  if(allocst /=0) stop "******Not enough memory*********"
  if (ipassc==1)then
     do k=1,n3do
        kk=m*n3do+k
        utemp_pros(:,:,k)=ut(:,:,kk)
        vtemp_pros(:,:,k)=vt(:,:,kk)
        wtemp_pros(:,:,k)=wt(:,:,kk)
        ptemp_pros(:,:,k)=pt(:,:,kk)
        ttemp_pros(:,:,k)=tt(:,:,kk)
     end do
  else
     do k=1,n3do
        kk=m*n3do+k
        utemp_pros(:,:,k)=ut(:,:,kk)
        vtemp_pros(:,:,k)=vt(:,:,kk)
        wtemp_pros(:,:,k)=wt(:,:,kk)
        ptemp_pros(:,:,k)=pt(:,:,kk)
     end do
  end if
  !write(*,*)'dividing data is done = ',itime,mynode
  !call mpi_barrier(MPI_COMM_WORLD,ierr)
  deallocate(ut,vt,wt,tt,pt)
  call intpolflddata_3dpod(utemp_pros,vtemp_pros,wtemp_pros,ptemp_pros,ttemp_pros,ysdo,ypdo,up,vp,wp,pp,tp)
  deallocate(utemp_pros,vtemp_pros,wtemp_pros,ptemp_pros,ttemp_pros)
end subroutine readtemp3dpod

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


subroutine compkernal_3d(u,v,w,x,y,z,c)
!  use integration
  use mpi
  use writedata
  use mpivariables
  use integration
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

subroutine podmodes_3d(u,evec,phi)
  implicit none
  real*8,intent(in),dimension(:,:,:,:) :: u
  real*8,intent(in),dimension(:,:)     :: evec
  real*8,intent(out),dimension(:,:,:,:):: phi
  
  integer                              ::i,j,k,t,m,d1,d2,d3,d4,d5
  
  d1=size(u,1)
  d2=size(u,2)
  d3=size(u,3)
  d4=size(u,4)
  d5=size(evec,2)

  phi=0.0

  do m=1,d5
     phi(:,:,:,m)=0.0
     do t=1,d4
        phi(:,:,:,m)=phi(:,:,:,m)+u(:,:,:,t)*evec(t,m)
     end do
  end do
  
!  phi=reshape(matmul(reshape(u,(/d1*d2*d3,d4/)),evec),(/d1,d2,d3,d4/))
end subroutine podmodes_3d

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
  end subroutine checkModeOrthog_3d
  

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! This subroutine reconstruct the velocity field for the mth mode    !!!!!!!!
!!!! INPUTS:                                                            !!!!!!!!
!!!!        mode      : eigenfunctions                                  !!!!!!!!
!!!!        coe       : coefficients                                    !!!!!!!!
!!!!        m         : mode number                                     !!!!!!!!
!!!!        tstp      : timestep                                        !!!!!!!!
!!!! OUTPUT:                                                            !!!!!!!!
!!!!        recnstfld : reconstructed field                             !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine modalrecnst_3d(mode,coe,m,tstp,i1toM,recnst_old,recnst_new)
    implicit none
    real*8,dimension(:,:,:,:),intent(in)     :: mode
    real*8,dimension(:,:),intent(in)         :: coe
    integer,intent(in)                       :: m,tstp,i1toM
    real*8,dimension(:,:,:),intent(in)       :: recnst_old
    real*8,dimension(:,:,:),intent(out)      :: recnst_new

    integer::d1,d2,d3,d4

    d1=size(mode,1)
    d2=size(mode,2)
    d3=size(mode,3)
    d4=size(mode,4)
    
    recnst_new(:,:,:)=recnst_old(:,:,:)+coe(tstp,m)*mode(:,:,:,m)
  
  end subroutine modalrecnst_3d

  subroutine modalrecnst_3d_1_M(mode,coe,mode_M,tstp,recnstfld)
    implicit none
    real*8,dimension(:,:,:,:),intent(in)   ::mode
    real*8,dimension(:,:),intent(in)       ::coe
    integer,intent(in)                     ::mode_M,tstp
    real*8,dimension(:,:,:),intent(out)  ::recnstfld

    integer::d1,d2,d3,d4,m

    d1=size(mode,1)
    d2=size(mode,2)
    d3=size(mode,3)
    d4=size(mode,4)
    recnstfld(:,:,:)=0.0
    
    do m=1,mode_M
       recnstfld(:,:,:)=recnstfld(:,:,:)+coe(tstp,m)*mode(:,:,:,m)
    end do
  end subroutine modalrecnst_3d_1_M

  
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
    do t=1,d4
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
    do t=1,d4
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
    integer                           :: yprime,xprime,nodex,nodey
    ! code is now 2d decomposed. Therefore, we need to find where the point
    ! is at. 
    !nodey goes from 0 to 23
    !nodex goes from 0 to 3
    nodey=j0/n2do
    nodex=k0/n3do
    node=nodey*filstopros+nodex
    if(mynode==node)then
       yprime=j0-nodey*n2do
       xprime=k0-nodex*n3do
       refval=ary3d(i0,yprime,xprime)
    end if
  end subroutine valatcnstz0y0x0


  subroutine tpcorl_3d_podfield(var1,var2,l_or_s,nfils,tsteps,i0,j0,k0)
    use readdata
    use MainParameters
    use channhdvariables
    use mpi
    use mpivariables
    use writedata
    use prelimcalvar
    implicit none

    integer,intent(in)                   ::var1,var2,i0,j0,k0,tsteps,nfils
    character(len=1)                     ::l_or_s
    integer                              ::d1,d2,d3,d4,t,node,rmsnode,alocst
    real*8,allocatable,dimension(:,:,:,:)::v1,v2
    real*8,allocatable,dimension(:,:,:)  ::rmsv1,rmsv2
    real*8,allocatable,dimension(:,:,:)  ::covar,covartavg
    real*8,allocatable,dimension(:,:,:)  ::inv
    real*8                               ::refval,invt,rmsrefval
    character(len=4)                     ::xprm,zprm
    character(len=3)                     ::yprm
    character(len=2)                     ::prosnum
    character(len=1)                     ::vari1,vari2

    write(prosnum,'(i2.2)')mynode

    allocate(v1(n1m,n2do+1,n3do,tsteps),stat=alocst)
    if(alocst/=0)stop"*****NOT SUFFICIENT MEMORY, V1************"
    if (var1==3)then
       if (l_or_s=='l')then
          call read4Darray('u_lsm'//trim(prosnum),'unformatted',v1)
       elseif(l_or_s=='s')then
          call read4Darray('u_ssm'//trim(prosnum),'unformatted',v1)
       end if
    elseif (var1==2)then
       if (l_or_s=='l')then
          call read4Darray('v_lsm'//trim(prosnum),'unformatted',v1)
       elseif(l_or_s=='s')then
          call read4Darray('v_ssm'//trim(prosnum),'unformatted',v1)
       end if
    elseif (var1==1)then
       if (l_or_s=='l')then
          call read4Darray('w_lsm'//trim(prosnum),'unformatted',v1)
       elseif(l_or_s=='s')then
          call read4Darray('w_ssm'//trim(prosnum),'unformatted',v1)
       end if
    end if

    allocate(v2(n1m,n2do+1,n3do,tsteps),stat=alocst)
    if(alocst/=0)stop"*****NOT SUFFICIENT MEMORY, V2************"
    if (var2==3)then
       if (l_or_s=='l')then
          call read4Darray('u_lsm'//trim(prosnum),'unformatted',v2)
       elseif(l_or_s=='s')then
          call read4Darray('u_ssm'//trim(prosnum),'unformatted',v2)
       end if
    elseif (var2==2)then
       if (l_or_s=='l')then
          call read4Darray('v_lsm'//trim(prosnum),'unformatted',v2)
       elseif(l_or_s=='s')then
          call read4Darray('v_ssm'//trim(prosnum),'unformatted',v2)
       end if
    elseif (var2==1)then
       if (l_or_s=='l')then
          call read4Darray('w_lsm'//trim(prosnum),'unformatted',v2)
       elseif(l_or_s=='s')then
          call read4Darray('w_ssm'//trim(prosnum),'unformatted',v2)
       end if
    end if

    d1=size(v1,1)
    d2=size(v1,2)
    d3=size(v1,3)
    d4=size(v1,4)

    allocate(covartavg(d1,d2,d3),rmsv1(d1,d2,d3),rmsv2(d1,d2,d3),stat=alocst)
    if(alocst/=0)stop"*****NOT SUFFICIENT MEMORY, corltavg************"
    
    rmsv1=0.
    rmsv2=0.
    covartavg=0.

    do t=1,d4

       call valatcnstz0y0x0(v2(:,:,:,t),i0,j0,k0,node,refval)
     !  call valatcnstz0y0x0(v2_ssm(:,:,:,t),i0,j0,k0,node_s,refval_s)
      ! call valatcnstz0y0x0(v2_lsm(:,:,:,t)+v2_ssm(:,:,:,t),i0,j0,k0,node_t,refval_t)

       call mpi_bcast(refval,1,mpi_real8,node,mpi_comm_world,ierr)
       !call mpi_bcast(refval_s,1,mpi_real8,node_s,mpi_comm_world,ierr)
       !call mpi_bcast(refval_t,1,mpi_real8,node_s,mpi_comm_world,ierr)

       !allocate(covar(d1,d2,d3))
       !covar=v1(:,:,:,t)*refval
       !covar_s=v1_ssm(:,:,:,t)*refval_s
       !covar_t=(v1_lsm(:,:,:,t)+v1_ssm(:,:,:,t))*refval_t

       covartavg=covartavg+v1(:,:,:,t)*refval
       !covartavg_s=covartavg_s+covar_s
       !covartavg_t=covartavg_t+covar_t

       !deallocate(covar)

       rmsv1=rmsv1+v1(:,:,:,t)*v1(:,:,:,t)
       !rmsv1_s=rmsv1_s+v1_ssm(:,:,:,t)*v1_ssm(:,:,:,t)
       !rmsv1_t=rmsv1_t+(v1_lsm(:,:,:,t)+v1_ssm(:,:,:,t))*(v1_lsm(:,:,:,t)+v1_ssm(:,:,:,t))

       rmsv2=rmsv2+v2(:,:,:,t)*v2(:,:,:,t)
      ! rmsv2_s=rmsv2_s+v2_ssm(:,:,:,t)*v2_ssm(:,:,:,t)
      ! rmsv2_t=rmsv2_t+(v2_lsm(:,:,:,t)+v2_ssm(:,:,:,t))*(v2_lsm(:,:,:,t)+v2_ssm(:,:,:,t))
    end do

    deallocate(v1,v2)
    
    invt=1.0/real(d4)


    rmsv1=sqrt(rmsv1*invt)
    rmsv2=sqrt(rmsv2*invt)
    !rmsv1_s=sqrt(rmsv1_s*inv)
    !rmsv2_s=sqrt(rmsv2_s*inv)
    !rmsv1_t=sqrt(rmsv1_t*inv)
    !rmsv2_t=sqrt(rmsv2_t*inv)
    call valatcnstz0y0x0(rmsv2,i0,j0,k0,rmsnode,rmsrefval)
    !call valatcnstz0y0x0(rmsv2_s,i0,j0,k0,rmsnode_s,rmsrefval_s)
    !call valatcnstz0y0x0(rmsv2_t,i0,j0,k0,rmsnode_t,rmsrefval_t)

    call mpi_bcast(rmsrefval,1,mpi_real8,rmsnode,mpi_comm_world,ierr)
    !call mpi_bcast(rmsrefval_s,1,mpi_real8,rmsnode_s,mpi_comm_world,ierr)
    !call mpi_bcast(rmsrefval_t,1,mpi_real8,rmsnode_t,mpi_comm_world,ierr)

    allocate(inv(d1,d2,d3))
    inv=1.0/(rmsv1*rmsrefval)
    !inv_s=1.0/(rmsv1_s*rmsrefval_s)
    !inv_t=1.0/(rmsv1_t*rmsrefval_t)


    covartavg=covartavg*invt*inv
    !covartavg_s=covartavg_s*inv*inv_s
    !covartavg_t=covartavg_t*inv*inv_t

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
     !if(l_or_s=='l')then
        !call sendrecv3dwrite_xy_decom(covartavg,1,nfils,'corcoe_lsm'//trim(vari1)//'_'//trim(vari2)//'_z0_'//trim(zprm)//'_y0_'//trim(yprm)//'_x0_'//trim(xprm)//'.dat')
     !elseif(l_or_s=='s')then
        !call sendrecv3dwrite_xy_decom(covartavg,2,nfils,'corcoe_ssm'//trim(vari1)//'_'//trim(vari2)//'_z0_'//trim(zprm)//'_y0_'//trim(yprm)//'_x0_'//trim(xprm)//'.dat')
     !end if
     !call sendrecv3dwrite_xy_decom(covartavg_t,3,nfils,'corcoe_tot'//trim(vari1)//'_'//trim(vari2)//'_z0_'//trim(zprm)//'_y0_'//trim(yprm)//'_x0_'//trim(xprm)//'.dat')
   end subroutine tpcorl_3d_podfield




 subroutine send_recv_elemntadd_3dary_xy_decom(ary3d,tag,nfils,summ)
    use mpi
    use mpivariables
    use prelimcalvar
    use channhdvariables
    implicit none
    real*8,          intent(in),dimension(:,:,:):: ary3d
    integer,         intent(in)                 :: tag,nfils
    real*8,intent(out)                          :: summ
    real*8                                      :: temp
    integer                                     :: i,j,k,d1,d2,d3,p,jj,kk,m,fn
    real*8,allocatable,dimension(:,:,:)         :: globary
    character(len=*),parameter                  :: fmt='(ES13.5E2)'
    d1=size(ary3d,1)
    d2=size(ary3d,2)
    d3=size(ary3d,3)

    

    if (mynode/=0)then
       call MPI_SEND(ary3d,d1*d2*d3,MPI_REAL8,0,tag,MPI_COMM_WORLD,ierr)
    else
       allocate(globary(d1,nfils*(d2-1)+1,filstopros*d3))

       do k=1,d3
          kk=k
          do j=1,d2
             jj=j
             do i=1,d1
                globary(i,jj,kk)=ary3d(i,j,k)
             end do
          end do
       end do
       do p=1,(numprocs-1)
          call MPI_RECV(ary3d,d1*d2*d3,MPI_REAL8,p,tag,MPI_COMM_WORLD,status,ierr)
          fn=int((p-mod(p,filstopros))/filstopros)
          m=mod(p,filstopros)
          do k=1,d3
             kk=k+m*d3
             do j=1,d2
                jj=j+fn*(d2-1)
                do i=1,d1
                   globary(i,jj,kk)=ary3d(i,j,k)
                end do
             end do
          end do
       end do
       temp=0.
       do k=1,size(globary,3)
          do j=1,size(globary,2)
             do i=1,size(globary,1)
                temp=temp+globary(i,j,k)
             end do
          end do
       end do
       summ=temp/real(size(globary,3)*size(globary,2)*size(globary,1))
       deallocate(globary)
    end if
  end subroutine send_recv_elemntadd_3dary_xy_decom
!! This routine computes integral length scales of the field reconstructed from each modes 1 to M
 ! subroutine comp_2dhomo_corel_pod_modes1_M(u,v,M,intlen)
 !   use modules
 !   use compfft
 !   implicit none
 !   real*8,intent(in),dimension(:,:,:)::u,v
 !   integer,intent(in)                ::M
 !   real*8,intent(out)                ::intlen


!    integer                           ::d1,d2,d3,i,j,k
!    complex(kind=8),allocatable,dimension(:,:,:):: uprihat,vprihat
!    complex(kind=8),allocatable,dimension(:,:)::refary
!    d1=size(u,1)
!    d2=size(u,2)
!    d3=size(u,3)
!    
!    allocate (uprihat(2*d1,d2,2*d3),vprihat(2*d1,d2,2*d3))
!    l1=size(uprihat,1)
!    l2=size(vprihat,2)
    
!    call twodfft(u,l1,l2,uprihat)
!    call twodfft(v,l1,l2,vprihat)
    
!    allocate (refary(l1,l2))
!!    
!    call ary2dcnsty0(vprihat,j0,node1,refary)
    
!    call mpi_bcast(refary
       
!  end subroutine comp_2dhomo_corel_pod_modes1_M
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reconstruct velocity field for each mode and each time step and prints streamwise velocity
! at a designated y location (y0) along x-direction
  subroutine velovector_at_y0(timestcount,yloc,filename)
    use mpi
    use mpivariables
    use channhdvariables
    use prelimcalvar
    use readdata
    use writedata
    implicit none
    integer,intent(in)                   ::yloc,timestcount
    character(len=*),intent(in)          ::filename
    real*8,allocatable,dimension(:,:,:,:)::phiu
    real*8,allocatable,dimension(:,:)    ::coe
    real*8,allocatable,dimension(:,:,:)  ::recnst_vel 
    integer                              ::node0,node1,node2,node3,node4,yprime0,d1,d2,d3,d4,m,t,k
    character(len=2)                     ::prosnum          
    !  d1=size(phi,1)
    !  d2=size(phi,2)
    !  d3=size(phi,3)
    !  d4=size(phi,4)
    ! Reconstructig velocity fields
    node0=yloc/n2do
    write(*,*)'node value node0 = ',node0
    node1=4*node0
    write(*,*)'node1 = ',node1
    node2=node1+1
    write(*,*)'node2 = ',node2
    node3=node1+2
    node4=node1+3
    if(mynode==0)then
       write(prosnum,'(i2.2)')node1
    elseif(mynode==1)then
       write(prosnum,'(i2.2)')node2
    elseif(mynode==2)then
       write(prosnum,'(i2.2)')node3
    elseif(mynode==3)then
       write(prosnum,'(i2.2)')node4
    end if
    allocate(phiu(n1m,n2do+1,n3do,timestcount))
    allocate(coe(timestcount,timestcount))
    d1=size(phiu,1)
    d2=size(phiu,2)
    d3=size(phiu,3)
    d4=size(phiu,4)
    call read4Darray('../pod_3d_new/phiu'//(prosnum),'unformatted',phiu)
    call read2DMatx(d4,d4,'../pod_3d_new/coeMatx.dat',coe)
    allocate(recnst_vel(d3,d4,d4))
    yprime0=yloc-node0
    do m=1,d4
       do t=1,d4
          do k=1,d3
             recnst_vel(k,t,m)=coe(t,m)*phiu(97,yprime0,k,m)
          end do
       end do
    end do
    write(prosnum,'(i2.2)')mynode
    open(31,file=filename//trim(prosnum)//'.dat')
    do m=1,d4
       do t=1,d4
          do k=1,d3
             write(31,'(es13.5e2)')recnst_vel(k,t,m)
          end do
       end do
    end do
    close(31)
    deallocate(recnst_vel)
  end subroutine velovector_at_y0
end module comppod3d
