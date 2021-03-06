program vorticityprogram
  USE channhdvariables
  USE mainparameters
  USE writedata
  USE arrayops
  USE readdata
  USE mpivariables
  USE prelimcalvar
  USE interpoldata
  USE mpi
  use vorticity
  use fluctuation
  use compfft
  use specrout
  IMPLICIT NONE
  
  real*8,allocatable,dimension(:)             ::yp,xp,zp,ypdo,ysdo,ys
  real*8,allocatable,dimension(:,:,:)         ::espec_wx,espec_wy,espec_wz
  real*8,allocatable,dimension(:,:,:)         ::espec_wx_tav,espec_wy_tav,espec_wz_tav
  real*8,allocatable,dimension(:,:,:)         ::wxtavg,wytavg,wztavg
  real*8,allocatable,dimension(:,:,:)         ::pp,dudy3
  real*8,allocatable,dimension(:,:,:)         ::wxprime,wyprime,wzprime
  real*8,allocatable,dimension(:,:,:)         ::phiwx,phiwy,phiwz
  complex(kind=8),allocatable,dimension(:,:,:)::wx_hat,wy_hat,wz_hat
  real*8,allocatable,dimension(:,:)           ::wxrmstx,wyrmstx,wzrmstx 
  real*8,allocatable,dimension(:)             ::wxrmstxz,wyrmstxz,wzrmstxz
  real*8,allocatable,dimension(:,:,:,:)       ::up,tp,omega,dudx,dudy,dudz
  integer                                     ::dt,timestcount,numtimesteps,i,j,k,itime
  integer                                     :: icorlavg,icen
  real*8                                      ::inv,sttime
  character*2                                 ::prosnum

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,mynode,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
  

  CALL readChannelData()
  CALL readPostData(sttime,numtimesteps,dt,icorlavg,icen)
  CALL prelimCal()
  
  ALLOCATE(xp(n3),yp(0:n2),zp(n1),ys(0:n2),ypdo(0:n2do+1),ysdo(0:n2do+1))
  CALL readCoordData(xp,yp,zp)
  CALL intpolCordi(mynode,numprocs,yp,ys,ypdo,ysdo)

  
  ALLOCATE(wxtavg(n1,n2do+1,n3),wytavg(n1,n2do+1,n3),&
       wztavg(n1,n2do+1,n3),dudy3(n1,n2do+1,n3))

  timestcount=0
  wxtavg=0.0
  wytavg=0.0
  wztavg=0.0
  dudy3=0.0
  ! Time loop starts to read data from field.data files
  stime=MPI_WTIME()
  DO ITime= INT(stTime),(INT(StTime)+(NumTimeSteps*Dt)),Dt
     TimeStCount=TimeStcount+1

     allocate(up(n1,n2do+1,n3,3),tp(n1,n2do+1,n3,1),pp(n1,n2do+1,n3))
     CALL readTempFieldData(mynode,itime,ypdo,ysdo,up,pp,tp)
     deallocate(pp,tp)

     allocate(dudx(n1,n2do+1,n3,3),dudy(n1,n2do+1,n3,3),dudz(n1,n2do+1,n3,3))
     call velograd(up,xp,ypdo,zp,dudx,dudy,dudz)
     deallocate(up)
     
     allocate(omega(n1,n2do+1,n3,3))
     call instvorticity(dudz,dudy,dudx,omega)
    
     CALL loopAdd3DAry(omega(:,:,:,1),wztavg)
     CALL loopAdd3DAry(omega(:,:,:,2),wytavg)
     CALL loopAdd3DAry(omega(:,:,:,3),wxtavg)
     CALL loopAdd3DAry(dudy(:,:,:,3),dudy3)
     deallocate(omega,dudx,dudy,dudz)
  END DO
  etime=MPI_WTIME()
  write(*,*)'time elapsed for data read',etime-stime,'for process',mynode
  write(*,*)'no of time steps',timestcount,'read from process',mynode

  inv=1.0/real(timestcount)

  ! Computing time average


  wxtavg=inv*wxtavg
  wytavg=inv*wytavg
  wztavg=inv*wztavg
  dudy3=inv*dudy3
  

  ! time mean without spatial averages.
  !write(prosnum,'(i2.2)')mynode
  !call write3Darray(wxtavg,'wxmean'//trim(prosnum),'unformatted')
  !call write3Darray(wytavg,'wymean'//trim(prosnum),'unformatted')
  !call write3Darray(wztavg,'wzmean'//trim(prosnum),'unformatted')  
  !call write3Darray(dudy3,'dudy3'//trim(prosnum),'unformatted')

  ! if blowing is not activated spatial averages can be taken in both
  ! -x and -y directions
!  call sendrecv3dwrite(wxtavg,1,'wxmean_3d.dat')
!  call sendrecv3dwrite(wytavg,2,'wymean_3d.dat')
!  call sendrecv3dwrite(wztavg,3,'wzmean_3d.dat')
!  call sendrecv3dwrite(dudy3,4,'dudy3mean_3d.dat')

  !deallocate(wxtavg,wytavg,wztavg,dudy3)

  !allocate(wxtavg(n1,n2do+1,n3),wytavg(n1,n2do+1,n3),wztavg(n1,n2do+1,n3))


! This part of the code computes rms of vorticity fluctuations

  !call read3Darray('wxmean'//trim(prosnum),'unformatted',wxtavg)
  !call read3Darray('wymean'//trim(prosnum),'unformatted',wytavg)
  !call read3Darray('wzmean'//trim(prosnum),'unformatted',wztavg)

  allocate(espec_wx_tav(n1m,n2do+1,n3m),espec_wy_tav(n1m,n2do+1,n3m),&
       espec_wz_tav(n1m,n2do+1,n3m))

  timestcount=0
  espec_wx_tav=0.0
  espec_wy_tav=0.0
  espec_wz_tav=0.0
  
  DO ITime= INT(stTime),(INT(StTime)+(NumTimeSteps*Dt)),Dt
     TimeStCount=TimeStcount+1

     allocate(up(n1,n2do+1,n3,3),tp(n1,n2do+1,n3,1),pp(n1,n2do+1,n3))
     CALL readTempFieldData(mynode,itime,ypdo,ysdo,up,pp,tp)
     deallocate(pp,tp)

     allocate(dudx(n1,n2do+1,n3,3),dudy(n1,n2do+1,n3,3),dudz(n1,n2do+1,n3,3))
     call velograd(up,xp,ypdo,zp,dudx,dudy,dudz)
     deallocate(up)
     
     allocate(omega(n1,n2do+1,n3,3))
     call instvorticity(dudz,dudy,dudx,omega)

     allocate(wxprime(n1,n2do+1,n3),wyprime(n1,n2do+1,n3),wzprime(n1,n2do+1,n3))
     
     call fluct3Dmean(omega(:,:,:,1),wztavg,wzprime)
     call fluct3Dmean(omega(:,:,:,2),wytavg,wyprime)
     call fluct3Dmean(omega(:,:,:,3),wxtavg,wxprime)
     
     deallocate(omega,dudx,dudy,dudz)

     allocate(phiwx(n1m,n2do+1,n3m),phiwy(n1m,n2do+1,n3m),phiwz(n1m,n2do+1,n3m))
     call asignphi(wzprime,phiwz)
     call asignphi(wyprime,phiwy)
     call asignphi(wxprime,phiwx)
     deallocate(wzprime,wxprime,wyprime)
     
     allocate(wx_hat(n1m,n2do+1,n3m),wy_hat(n1m,n2do+1,n3m),wz_hat(n1m,n2do+1,n3m))

     call onedfft3dary(phiwz,n1m,1,wz_hat)
     call onedfft3dary(phiwy,n1m,1,wy_hat)
     call onedfft3dary(phiwx,n1m,1,wx_hat)
     
     deallocate(phiwx,phiwy,phiwz)
     allocate(espec_wx(n1m,n2do+1,n3m),espec_wy(n1m,n2do+1,n3m),espec_wz(n1m,n2do+1,n3m))
     
     call conjmult3dary(wz_hat,wz_hat,espec_wz)
     call conjmult3dary(wy_hat,wy_hat,espec_wy)
     call conjmult3dary(wx_hat,wx_hat,espec_wx)
     
     deallocate(wz_hat,wy_hat,wx_hat)
     
     call loopAdd3DAry(espec_wx,espec_wx_tav)
     call loopAdd3DAry(espec_wy,espec_wy_tav)
     call loopAdd3DAry(espec_wz,espec_wz_tav)
     
     deallocate(espec_wx,espec_wy,espec_wz)
  END DO
  
  deallocate(wxtavg,wytavg,wztavg,dudy3)

  inv=1.0/real(timestcount)
  
  espec_wx_tav=(inv*espec_wx_tav)
  espec_wy_tav=(inv*espec_wy_tav)
  espec_wz_tav=(inv*espec_wz_tav)
  
  call spectwrite_2dplane(espec_wx_tav,1,2,158,'espec_wx')
  call spectwrite_2dplane(espec_wy_tav,2,2,158,'espec_wy')
  call spectwrite_2dplane(espec_wz_tav,3,2,158,'espec_wz')
  deallocate(espec_wx_tav,espec_wy_tav,espec_wz_tav)
  CALL MPI_FINALIZE(ierr)
END PROGRAM vorticityprogram
