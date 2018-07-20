program anisotropy
  use channhdvariables
  use mainparameters
  use writedata
  use arrayops
  use readdata
  use mpivariables
  use prelimcalvar
  use interpoldata
  use fluctuation
  use mpi
  use anisotropyroutines

  implicit none

  real *8,allocatable,dimension(:)      ::xp,yp,ypdo,ysdo,ys,zp
  real *8,allocatable,dimension(:,:,:)  ::utavg,vtavg,wtavg,ttavg,ptavg,pp,uu,uv,uw,vv,vw,ww
  real *8,allocatable,dimension(:,:,:)  ::uu2,uv2,uw2,vv2,vw2,ww2,eta,zeta
  real *8,allocatable,dimension(:,:,:)  ::uprime,vprime,wprime,pprime,tprime
  real *8,allocatable,dimension(:,:,:,:)::up,tp
  real *8,allocatable,dimension(:,:)    ::reystins,reystavg
  logical                               ::fuu,fuv,fuw,fvv,fvw,fww

  integer::dt,timestcount,numtimesteps,i,j,k,itime,xcount,ycount,zcount,icorlavg,icen
  real*8::inv,sttime
  character*5::string
  character*2::prosnum
  character*5::time
  
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mynode,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

  call readChannelData()
  call readPostData(sttime,numtimesteps,dt,icorlavg,icen)
  call prelimcal()

  allocate(xp(n3),yp(0:n2),zp(n1),ys(0:n2),ypdo(0:n2do+1),ysdo(0:n2do+1))

  call readCoordData(xp,yp,zp)
  call intpolCordi(mynode,numprocs,yp,ys,ypdo,ysdo)

  allocate(wtavg(n1,n2do+1,n3),vtavg(n1,n2do+1,n3),&
       utavg(n1,n2do+1,n3),ttavg(n1,n2do+1,n3),ptavg(n1,n2do+1,n3))


  write(prosnum,'(i2.2)')mynode

  call read3Darray('tmean'//trim(prosnum),'unformatted',ttavg)
  call read3Darray('wmean'//trim(prosnum),'unformatted',wtavg)
  call read3Darray('vmean'//trim(prosnum),'unformatted',vtavg)
  call read3Darray('umean'//trim(prosnum),'unformatted',utavg)

  timestcount = 0

  stime = MPI_WTIME()
  allocate(uu(n1,n2do+1,n3),uv(n1,n2do+1,n3),uw(n1,n2do+1,n3),vv(n1,n2do+1,n3),vw(n1,n2do+1,n3),ww(n1,n2do+1,n3))

  inquire(file='uu_'//trim(prosnum),exist=fuu)
  inquire(file='uv_'//trim(prosnum),exist=fuv)
  inquire(file='uw_'//trim(prosnum),exist=fuw)
  inquire(file='vv_'//trim(prosnum),exist=fvv)
  inquire(file='vw_'//trim(prosnum),exist=fvw)
  inquire(file='ww_'//trim(prosnum),exist=fww)


  if(fuv.and.fuv.and.fuw.and.fvv.and.fvw.and.fww)then

     write(*,*)'The files exists and loading'
     call read3Darray('uu_'//trim(prosnum),'unformatted',uu)
     call read3Darray('uv_'//trim(prosnum),'unformatted',uv)
     call read3Darray('uw_'//trim(prosnum),'unformatted',uw)
     call read3Darray('vv_'//trim(prosnum),'unformatted',vv)
     call read3Darray('vw_'//trim(prosnum),'unformatted',vw)
     call read3Darray('ww_'//trim(prosnum),'unformatted',ww)

     allocate(eta(n1,n2do+1,n3),zeta(n1,n2do+1,n3))
     call sendrecv3dwrite(uv,3,'uv_3d.dat')
     call reystinv(uu,uv,uw,vv,vw,ww,eta,zeta)
     deallocate(uu,uv,uw,vv,vw,ww)
     call sendrecv3dwrite(eta,1,'eta_3d.dat')
     call sendrecv3dwrite(zeta,2,'zeta_3d.dat')
     deallocate(eta,zeta)

  else 
     write(*,*)'The files do not exist and creating'
     Do itime = int(sttime),(int(sttime)+(numtimesteps*dt)),dt
        allocate(up(n1,n2do+1,n3,3),tp(n1,n2do+1,n3,1),pp(n1,n2do+1,n3))
        timestcount = timestcount+1
        call readTempFieldData(mynode,itime,ypdo,ysdo,up,pp,tp)

        allocate(tprime(n1,n2do+1,n3),wprime(n1,n2do+1,n3),vprime(n1,n2do+1,n3),uprime(n1,n2do+1,n3),pprime(n1,n2do+1,n3))

        call fluct3Dmean(tp(:,:,:,1),ttavg,tprime)
        call fluct3Dmean(up(:,:,:,1),wtavg,wprime)
        call fluct3Dmean(up(:,:,:,2),vtavg,vprime)
        call fluct3Dmean(up(:,:,:,3),utavg,uprime)
        deallocate(up,tp,pp)
        allocate(uu2(n1,n2do+1,n3),uv2(n1,n2do+1,n3),uw2(n1,n2do+1,n3),vv2(n1,n2do+1,n3),vw2(n1,n2do+1,n3),ww2(n1,n2do+1,n3))
        uu2 = uprime*uprime
        uv2 = uprime*vprime
        uw2 = uprime*wprime
        vv2 = vprime*vprime
        vw2 = vprime*wprime
        ww2 = wprime*wprime
        deallocate(uprime,vprime,wprime,tprime,pprime)
        call loopAdd3DAry(uu2,uu)
        call loopAdd3DAry(uv2,uv)
        call loopAdd3DAry(uw2,uw)
        call loopAdd3DAry(vv2,vv)
        CALL loopAdd3DAry(vw2,vw)
        call loopAdd3DAry(ww2,ww)
          
        deallocate(uu2,uv2,uw2,vv2,vw2,ww2)
     end do
     write(*,*)'Completed calculating the reynolds stresses at each point and time averaged'
     deallocate(utavg,wtavg,vtavg,ttavg)
     inv = 1.0/real(timestcount)
     uu = inv*uu
     uv = inv*uv
     uw = inv*uw
     vv = inv*vv
     vw = inv*vw
     ww = inv*ww
     
     call write3Darray(uu,'uu_'//trim(prosnum),'unformatted')
     call write3Darray(uv,'uv_'//trim(prosnum),'unformatted')
     call write3Darray(uw,'uw_'//trim(prosnum),'unformatted')
     call write3Darray(vv,'vv_'//trim(prosnum),'unformatted')
     call write3Darray(vw,'vw_'//trim(prosnum),'unformatted')
     call write3Darray(ww,'ww_'//trim(prosnum),'unformatted')
     
     allocate(eta(n1,n2do+1,n3),zeta(n1,n2do+1,n3))
     call sendrecv3dwrite(uv,3,'uv_3d.dat')
     call reystinv(uu,uv,uw,vv,vw,ww,eta,zeta)
     deallocate(uu,uv,uw,vv,vw,ww)
     call sendrecv3dwrite(eta,1,'eta_3d.dat')
     call sendrecv3dwrite(zeta,2,'zeta_3d.dat')
     deallocate(eta,zeta)

  end if

  deallocate(xp,ypdo,ysdo,yp,zp)
  call MPI_FINALIZE(ierr) 


end program anisotropy
