program spectra
  
  use channhdvariables
  use mainparameters
  use writedata
  use arrayops
  use readdata
  use mpivariables
  use prelimcalvar
  use interpoldata
  use fluctuation
  use compfft
  use mpi
  use specrout
  
  implicit none
  
  real*8,allocatable,dimension(:)::yp,xp,zp,ypdo,ysdo,ys
  real*8,allocatable,dimension(:,:,:)::utavg,vtavg,wtavg,ttavg,ptavg,pp
  complex(kind=8),allocatable,dimension(:,:,:)::u_hat,v_hat,w_hat,t_hat
  complex(kind=8),allocatable,dimension(:,:,:)::espec_uu,espec_vv,espec_ww,espec_uv
  real*8,allocatable,dimension(:,:,:)::abs_espct_uu,abs_espct_vv,abs_espct_ww,abs_espct_uv
  complex(kind=8),allocatable,dimension(:,:,:)::espct_uu,espct_vv,espct_ww,espct_uv
  real*8,allocatable,dimension(:,:,:)::uprime,vprime,wprime,tprime,pprime
  real*8,allocatable,dimension(:,:,:)::phiu,phiv,phiw
  REAL*8,ALLOCATABLE,DIMENSION(:)::utxzavg,vtxzavg,wtxzavg,ttxzavg
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:)::up,tp
  INTEGER::timestcount,numtimesteps,i,j,k
  integer::icorlavg,icen
  REAL*8::inv,sttime,inv2,dt,itime
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
  allocate(espct_uu(n1m,n2do+1,n3m),espct_vv(n1m,n2do+1,n3m),&
       espct_ww(n1m,n2do+1,n3m),espct_uv(n1m,n2do+1,n3m))
  espct_uu=cmplx(0.)
  espct_vv=cmplx(0.)
  espct_ww=cmplx(0.)
  espct_uv=cmplx(0.)
  wtavg=0.
  utavg=0.
  vtavg=0.
  ttavg=0.
  ptavg=0
  timestcount=0
  DO ITime= stTime,(StTime)+(NumTimeSteps*dt),dt
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

  ! Time loop starts to read data from field.data files
  stime=MPI_WTIME()
  do ITime= stTime,StTime+(NumTimeSteps*dt),dt
     TimeStCount=TimeStcount+1

     allocate(up(n1,n2do+1,n3,3),tp(n1,n2do+1,n3,1),pp(n1,n2do+1,n3))

     call readTempFieldData(mynode,itime,ypdo,ysdo,up,pp,tp)
     allocate(tprime(n1,n2do+1,n3),wprime(n1,n2do+1,n3),&
          vprime(n1,n2do+1,n3),uprime(n1,n2do+1,n3),pprime(n1,n2do+1,n3))
          

     call fluct3Dmean(tp(:,:,:,1),ttavg,tprime)
     call fluct3Dmean(up(:,:,:,1),wtavg,wprime)
     call fluct3Dmean(up(:,:,:,2),vtavg,vprime)
     call fluct3Dmean(up(:,:,:,3),utavg,uprime)
     call fluct3Dmean(pp(:,:,:),ptavg,pprime)
     
     deallocate(up,tp,pp)
     
     allocate(phiu(n1m,n2do+1,n3m),phiv(n1m,n2do+1,n3m),&
          phiw(n1m,n2do+1,n3m)) 

     call asignphi(uprime,phiu)
     call asignphi(vprime,phiv)
     call asignphi(wprime,phiw)
     
    ! if(mynode==10)then
    !    write(*,*)'phiu(70,4,588) = ',phiu(70,4,588)
    ! end if
     deallocate(uprime,vprime,tprime,wprime,pprime)
     
     allocate(u_hat(n1m,n2do+1,n3m),v_hat(n1m,n2do+1,n3m),&
          w_hat(n1m,n2do+1,n3m),t_hat(n1,n2do+1,n3m))

     call onedfft3dary(phiu,n1m,1,u_hat)
     call onedfft3dary(phiv,n1m,1,v_hat)
     call onedfft3dary(phiw,n1m,1,w_hat)
   ! call onedfft3dary(tprime,n1,1,t_hat)
     !if(mynode==10)then
     !   write(*,*)'u_hat(70,4,588) = ',u_hat(70,4,588)
     !end if
     deallocate(phiu,phiv,phiw)
     
     allocate(espec_uu(n1m,n2do+1,n3m),espec_vv(n1m,n2do+1,n3m),&
          espec_ww(n1m,n2do+1,n3m),espec_uv(n1m,n2do+1,n3m))
     
     call conjmult3dary(u_hat,u_hat,espec_uu)
     call conjmult3dary(v_hat,v_hat,espec_vv)
     call conjmult3dary(w_hat,w_hat,espec_ww)
     call conjmult3dary(u_hat,v_hat,espec_uv)
     deallocate(u_hat,v_hat,w_hat,t_hat)
     
     call loopadd3dcmplxary(espec_uu,espct_uu)
     call loopadd3dcmplxary(espec_vv,espct_vv)
     call loopadd3dcmplxary(espec_ww,espct_ww)
     call loopadd3dcmplxary(espec_uv,espct_uv)
     !if(mynode==10)then
     !   write(*,*)'espct_uu(70,4,588) = ',espct_uu(70,4,588)
     !end if
     deallocate(espec_uu,espec_vv,espec_ww,espec_uv)
  end do
  etime=MPI_WTIME()
!  write(*,*)'time elapsed for data read',etime-stime,'for process',mynode
  write(*,*)'no of time steps',timestcount,'read from process',mynode
  deallocate(wtavg,vtavg,utavg,ttavg)
  inv=1.0/real(timestcount)
  
  allocate(abs_espct_uu(n1m,n2do+1,n3m),abs_espct_vv(n1m,n2do+1,n3m),&
       abs_espct_ww(n1m,n2do+1,n3m),abs_espct_uv(n1m,n2do+1,n3m))
  
  abs_espct_uu=abs(inv*espct_uu)
  abs_espct_vv=abs(inv*espct_vv)
  abs_espct_ww=abs(inv*espct_ww)
  abs_espct_uv=abs(inv*espct_uv)
  
  deallocate(espct_uu,espct_vv,espct_ww)
  
  call spectwrite_2dplane(abs_espct_uu,1,1,96,'espec_uu')
  call spectwrite_2dplane(abs_espct_vv,2,1,96,'espec_vv')
  call spectwrite_2dplane(abs_espct_ww,3,1,96,'espec_ww')
  call spectwrite_2dplane(abs_espct_uv,4,1,96,'espec_uv')

  call mpi_finalize(ierr)

end program spectra

