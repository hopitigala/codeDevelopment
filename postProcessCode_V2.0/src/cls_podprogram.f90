program programpod
  use mpi
  use mpivariables
  use mainparameters
  use comppod
  use readdata
  use channhdvariables
  use interpoldata
  use eigenvalues
  use prelimcalvar
  use fluctuation
  use writedata
  
  implicit none

  real*8,allocatable,dimension(:)      :: yp,xp,zp,ypdo,ysdo,ys,eval
  real*8,allocatable,dimension(:)      :: uv1_lsm,uv1_ssm,uv1_tot
  real*8,allocatable,dimension(:,:)    :: uv2_lsm,uv2_ssm,uv2_tot
  real*8,allocatable,dimension(:,:,:,:):: up,tp
  real*8,allocatable,dimension(:,:,:)  :: wtavg,vtavg,utavg,ttavg,tprime
  real*8,allocatable,dimension(:,:,:)  :: uprime,vprime,wprime,pp
  real*8,allocatable,dimension(:,:,:)  :: uprit,vprit,wprit,tprit,phiu,phiv,phiw
  real*8,allocatable,dimension(:,:,:)  :: u_lsm,v_lsm,u_ssm,v_ssm,w_lsm,w_ssm
  real*8,allocatable,dimension(:,:)    :: c,evec,recnstu,recnstv,coe

  
  integer                              :: timestcount,itime,numtimesteps,iplane,dt
  integer                              :: icorl,x0,t,lsmlim,ireyst,inst,tt,iscfx
  integer                              :: icorlavg,icen,var1,var2,j0,k0
  integer                              :: fexist1,fexist2,fexist3,fexist4,fexist5,fexist6
  real*8                               :: sttime
  character(len=2)                     :: prosnum
  character(len=4)                     :: xloc
  character(len=5)                     :: time
  character(len=30)                    :: filename
  character(len=30)                    :: fname_l
  character(len=30)                    :: fname_s
  
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,mynode,ierr)
  call mpi_comm_size(mpi_comm_world,numprocs,ierr)
  
  call readChannelData()
  call readPostData(sttime,numtimesteps,dt,icorlavg,icen)
  call readpoddata(iplane,x0,lsmlim,ireyst,iscfx,icorl,var1,var2,j0,k0,inst)
  call prelimCal()

  allocate(xp(n3),yp(0:n2),zp(n1),ys(0:n2),ypdo(0:n2do+1),ysdo(0:n2do+1))

  call readCoordData(xp,yp,zp)
  call intpolCordi(mynode,numprocs,yp,ys,ypdo,ysdo)

  allocate(wtavg(n1,n2do+1,n3),vtavg(n1,n2do+1,n3),&
       utavg(n1,n2do+1,n3),ttavg(n1,n2do+1,n3))
  
  ! Read mean field
  write(prosnum,'(i2.2)')mynode
  call read3Darray('../../tmean/1500ts/tmean'//trim(prosnum),'unformatted',ttavg)
  call read3Darray('../../tmean/1500ts/wmean'//trim(prosnum),'unformatted',wtavg)
  call read3Darray('../../tmean/1500ts/vmean'//trim(prosnum),'unformatted',vtavg)
  call read3Darray('../../tmean/1500ts/umean'//trim(prosnum),'unformatted',utavg)
  if(iplane==1)then
     allocate(uprit(n2do+1,n3,numtimesteps+1),vprit(n2do+1,n3,numtimesteps+1),&
          wprit(n2do+1,n3,numtimesteps+1),tprit(n2do+1,n3,numtimesteps+1))
  elseif(iplane==2)then
     allocate(uprit(n1,n2do+1,numtimesteps+1),vprit(n1,n2do+1,numtimesteps+1),&
          wprit(n1,n2do+1,numtimesteps+1),tprit(n1,n2do+1,numtimesteps+1))
  else
     allocate(uprit(n1,n3,numtimesteps+1),vprit(n1,n3,numtimesteps+1),&
          wprit(n1,n3,numtimesteps+1),tprit(n1,n3,numtimesteps+1))
  end if
  
  timestcount=0
  do itime=int(sttime),(int(sttime)+(numtimesteps*dt)),dt
     timestcount=timestcount+1
     
     allocate(up(n1,n2do+1,n3,3),tp(n1,n2do+1,n3,1),pp(n1,n2do+1,n3))
     
     call readTempFieldData(mynode,itime,ypdo,ysdo,up,pp,tp)
     
     allocate(tprime(n1,n2do+1,n3),wprime(n1,n2do+1,n3),vprime(n1,n2do+1,n3),&
          uprime(n1,n2do+1,n3))

     call fluct3Dmean(tp(:,:,:,1),ttavg,tprime)
     call fluct3Dmean(up(:,:,:,1),wtavg,wprime)
     call fluct3Dmean(up(:,:,:,2),vtavg,vprime)
     call fluct3Dmean(up(:,:,:,3),utavg,uprime)
     
     deallocate(up,tp,pp)
  
     call fluc_2Dplane(uprime,x0,iplane,timestcount,uprit)
     call fluc_2Dplane(vprime,x0,iplane,timestcount,vprit)
     call fluc_2Dplane(wprime,x0,iplane,timestcount,wprit)
     call fluc_2Dplane(tprime,x0,iplane,timestcount,tprit)
     
     write(time,'(i5.5)')itime
     
     deallocate(tprime,uprime,vprime,wprime)
     if (iplane==1)then
        call sendrecv2dwrite(uprit(:,:,timestcount),1,iplane,'uprit_xy'//trim(time)//'.dat')
        call sendrecv2dwrite(vprit(:,:,timestcount),2,iplane,'vprit_xy'//trim(time)//'.dat')
        call sendrecv2dwrite(wprit(:,:,timestcount),3,iplane,'wprit_xy'//trim(time)//'.dat')
     elseif(iplane==2)then
        call sendrecv2dwrite(uprit(:,:,timestcount),1,iplane,'uprit_yz'//trim(time)//'.dat')
        call sendrecv2dwrite(vprit(:,:,timestcount),2,iplane,'vprit_yz'//trim(time)//'.dat')
        call sendrecv2dwrite(wprit(:,:,timestcount),3,iplane,'wprit_yz'//trim(time)//'.dat')
     else
        call sendrecv2dwrite(uprit(:,:,timestcount),1,iplane,'uprit_xz'//trim(time)//'.dat')
        call sendrecv2dwrite(vprit(:,:,timestcount),2,iplane,'vprit_xz'//trim(time)//'.dat')
        call sendrecv2dwrite(wprit(:,:,timestcount),3,iplane,'wprit_xz'//trim(time)//'.dat')
     end if
     write(*,*)'reading time step',itime,'is done by',mynode
  end do

  call mpi_finalize(ierr)
end program programpod
        
        
  
