program inst_2d_prog
  
  use channhdvariables
  use mainparameters
  use writedata
  use arrayops
  use readdata
  use mpivariables
  use prelimcalvar
  use interpoldata
  use fluctuation
!  use compfft
  use mpi
!  use specrout
  
  implicit none
  
  real*8,allocatable,dimension(:)::yp,xp,zp,ypdo,ysdo,ys
  real*8,allocatable,dimension(:,:,:)::utavg,vtavg,wtavg,ttavg,ptavg,pp
 ! complex(kind=8),allocatable,dimension(:,:,:)::u_hat,v_hat,w_hat,t_hat
 ! complex(kind=8),allocatable,dimension(:,:,:)::espec_uu,espec_vv,espec_ww,espec_uv
 ! real*8,allocatable,dimension(:,:,:)::espec_uu,espec_vv,espec_ww,espec_uv
 ! real*8,allocatable,dimension(:,:,:)::abs_espct_uu,abs_espct_vv,abs_espct_ww,abs_espct_uv
 ! real*8,allocatable,dimension(:,:,:)::espct_uu,espct_vv,espct_ww,espct_uv
  real*8,allocatable,dimension(:,:,:)::uprime,vprime,wprime,tprime,pprime
  real*8,allocatable,dimension(:,:)::vprime2d,uprime2d,wprime2d,tprime2d,pprime2d,utavg2d 
 ! real*8,allocatable,dimension(:,:,:)::phiu,phiv,phiw
 ! REAL*8,allocatable,dimension(:)::utxzavg,vtxzavg,wtxzavg,ttxzavg
  REAL*8,allocatable,dimension(:,:,:,:)::up,tp
  INTEGER::dt,timestcount,numtimesteps,i,j,k,itime
  integer::icorlavg,icen
  REAL*8::inv,sttime,inv2
  character*2::prosnum
  character*5::time
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mynode,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)


  call readChannelData()
  call readPostData(sttime,numtimesteps,dt,icorlavg,icen)
  call prelimCal()

  allocate(xp(n3),yp(0:n2),zp(n1),ys(0:n2),ypdo(0:n2do+1),ysdo(0:n2do+1))
  call readCoordData(xp,yp,zp)
  call intpolCordi(mynode,numprocs,yp,ys,ypdo,ysdo)


  allocate(wtavg(n1,n2do+1,n3),vtavg(n1,n2do+1,n3),&
       utavg(n1,n2do+1,n3),ttavg(n1,n2do+1,n3),ptavg(n1,n2do+1,n3))
 ! allocate(espct_uu(n1m,n2do+1,n3m),espct_vv(n1m,n2do+1,n3m),&
 !      espct_ww(n1m,n2do+1,n3m),espct_uv(n1m,n2do+1,n3m))
 ! espct_uu=cmplx(0.)
 ! espct_vv=cmplx(0.)
 ! espct_ww=cmplx(0.)
 ! espct_uv=cmplx(0.)
 ! espct_uu=0.
 ! espct_vv=0.
 ! espct_ww=0.
 ! espct_uv=0.
  wtavg=0.
  utavg=0.
  vtavg=0.
  ttavg=0.
  ptavg=0
  timestcount=0
  DO ITime= INT(stTime),(INT(StTime)+(NumTimeSteps*Dt)),Dt
     TimeStCount=TimeStcount+1

     ALLOCATE(up(n1,n2do+1,n3,3),tp(n1,n2do+1,n3,1),pp(n1,n2do+1,n3))

     CALL readTempFieldData(mynode,itime,ypdo,ysdo,up,pp,tp)

     CALL loopAdd3DAry(up(:,:,:,1),wtavg)
     CALL loopAdd3DAry(up(:,:,:,2),vtavg)
     CALL loopAdd3DAry(up(:,:,:,3),utavg)
     CALL loopAdd3DAry(tp(:,:,:,1),ttavg)
     call loopAdd3Dary(pp,ptavg)
     DEALLOCATE(up,tp,pp)
  END DO

  ! Read mean field
!  write(prosnum,'(i2.2)')mynode
!  call read3Darray('../../tmean/1500ts/tmean'//trim(prosnum),'unformatted',ttavg)
!  call read3Darray('../../tmean/1500ts/wmean'//trim(prosnum),'unformatted',wtavg)
!  call read3Darray('../../tmean/1500ts/vmean'//trim(prosnum),'unformatted',vtavg)
!  call read3Darray('../../tmean/1500ts/umean'//trim(prosnum),'unformatted',utavg)
!  call read3Darray('../../tmean/1500ts/pmean'//trim(prosnum),'unformatted',ptavg)
  inv2=1.0/real(timestcount)
  utavg=utavg*inv2
  vtavg=vtavg*inv2
  wtavg=wtavg*inv2
  ttavg=ttavg*inv2
  ptavg=ptavg*inv2
  timestcount=0
 allocate(utavg2d(size(utavg,2),size(utavg,3)))
  do k=2,size(utavg,3)
     do j=1,size(utavg,2)
        utavg2d(j,k)=utavg(96,j,k)
     end do
  end do
  call sendrecv2dwrite(utavg2d,3,1,'u_tavg_xy_96.dat')
 deallocate(utavg2d)

  ! Time loop starts to read data from field.data files
 ! stime=MPI_WTIME()
 ! do ITime= INT(stTime),(INT(StTime)+(NumTimeSteps*Dt)),Dt
 !    TimeStCount=TimeStcount+1
 !    write(time,'(i5.5)')itime
 !    allocate(up(n1,n2do+1,n3,3),tp(n1,n2do+1,n3,1),pp(n1,n2do+1,n3))

 !    call readTempFieldData(mynode,itime,ypdo,ysdo,up,pp,tp)
 !    allocate(tprime(n1,n2do+1,n3),wprime(n1,n2do+1,n3),&
 !         vprime(n1,n2do+1,n3),uprime(n1,n2do+1,n3),pprime(n1,n2do+1,n3))
          

 !    call fluct3Dmean(tp(:,:,:,1),ttavg,tprime)
 !    call fluct3Dmean(up(:,:,:,1),wtavg,wprime)
 !    call fluct3Dmean(up(:,:,:,2),vtavg,vprime)
 !    call fluct3Dmean(up(:,:,:,3),utavg,uprime)
 !    call fluct3Dmean(pp(:,:,:),ptavg,pprime)
     
 !    deallocate(up,tp,pp)
 !    allocate(uprime2d(size(uprime,2),size(uprime,3)),vprime2d(size(vprime,2),size(vprime,3)))
    ! deallocate(uprime,vprime,tprime,wprime,pprime)
  !   do k=1,size(uprime,3)
  !      do j=1,size(uprime,2)
 !          uprime2d(j,k)=uprime(96,j,k)
 !          vprime2d(j,k)=vprime(96,j,k)
 !       end do
 !    end do
! deallocate(uprime,vprime,tprime,wprime,pprime)
!  call sendrecv2dwrite(uprime2d,1,1,'uprime_xy_96_'//trim(time)//'.dat')
!  call sendrecv2dwrite(vprime2d,2,1,'vprime_xy_96_'//trim(time)//'.dat')
! deallocate(uprime2d,vprime2d)
!  end do
!  etime=MPI_WTIME()
!  write(*,*)'time elapsed for data read',etime-stime,'for process',mynode
 ! write(*,*)'no of time steps',timestcount,'read from process',mynode
 ! allocate(utavg2d(size(utavg,2),size(utavg,3)))
 ! do k=2,size(utavg,3)
 !    do j=1,size(utavg,2)
 !       utavg2d(j,k)=utavg(96,j,k)
 !    end do
 ! end do
 ! deallocate(utavg2d)
 ! call sendrecv2dwrite(utavg2d,3,1,'u_tavg_xy_96.dat')
 ! deallocate(wtavg,vtavg,utavg,ttavg)
  
  call mpi_finalize(ierr)

end program inst_2d_prog

