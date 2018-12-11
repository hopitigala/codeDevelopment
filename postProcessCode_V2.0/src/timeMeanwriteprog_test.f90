program programpod
  use mpi
  use mpivariables
  use mainparameters
  use comppod3d
  use readdata
  use channhdvariables
  !use interpoldata
  use eigenvalues
  use prelimcalvar
!  use fluctuation
  use writedata

  implicit none

  real*8,allocatable,dimension(:)        :: yp,xp,zp,ypdo,ysdo,ys,eval,xpdo
  real*8,allocatable,dimension(:)        :: uv_lsm,uv_ssm,uv_tot,norm2
  real*8,allocatable,dimension(:,:,:)    :: uv_lsm_3d,uv_ssm_3d,uv_tot_3d
  real*8,allocatable,dimension(:,:,:)    :: up,tp,vp,wp
  real*8,allocatable,dimension(:,:,:)    :: wtavg,vtavg,utavg,ttavg,tprime
  real*8,allocatable,dimension(:,:,:)    :: uprime,vprime,wprime,pp
  real*8,allocatable,dimension(:,:,:,:)  :: uprit,vprit,wprit,tprit,phiu,phiv,phiw
  real*8,allocatable,dimension(:,:,:,:)  :: u_lsm,v_lsm,u_ssm,v_ssm,w_ssm,w_lsm
  real*8,allocatable,dimension(:,:)      :: c,evec,coe
  real*8,allocatable,dimension(:,:,:)    :: recnstu,recnstv


  integer                              :: timestcount,itime,numtimesteps,dt,nfils,m,i1toM
  integer                              :: icorl,i0,t,lsmlim,ireyst,inst,tt,iscfx,icomppod
  integer                              :: icorlavg,icen,var1,var2,j0,k0,irecnst
  integer                              :: fexist1,fexist2,fexist3,fexist4,fexist5,fexist6
  real*8                               :: sttime,norminv
  character(len=2)                     :: prosnum
  character(len=4)                     :: xloc,time
  character(len=30)                    :: filename

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,mynode,ierr)
  call mpi_comm_size(mpi_comm_world,numprocs,ierr)

  call readChannelData()
  call readpostdata_3dpod(sttime,numtimesteps,dt,icorlavg,icen,nfils)
  !call readpod3ddata(icomppod,irecnst,lsmlim,ireyst,iscfx,inst,icorl,var1,var2,i0,j0,k0,i1toM)
  call prelimcal_3dpod(nfils)

  !allocate(xp(n3),yp(0:n2),zp(n1),ys(0:n2),ypdo(0:n2do+1),ysdo(0:n2do+1),xpdo(n3do))

  !call readCoordData(xp,yp,zp)
  !call intpolcordi_3dpod(yp,xp,ys,ypdo,ysdo,xpdo)

  allocate(utavg(n1m,n2do+1,n3do))

     ! Read mean field
  write(prosnum,'(i2.2)')mynode
  !call read3Darray('../../tmean/1500ts_3dpod/tmean'//trim(prosnum),'unformatted',ttavg)
  !call read3Darray('../../tmean/1500ts_3dpod/wmean'//trim(prosnum),'unformatted',wtavg)
  !call read3Darray('../../tmean/1500ts_3dpod/vmean'//trim(prosnum),'unformatted',vtavg)
  call read3Darray('../../tmean/1500ts_3dpod/umean'//trim(prosnum),'unformatted',utavg)
  
  call sendrecv3dwrite_xy_decom(utavg,1,nfils,'u2_lsm_3d.dat')
  !call sendrecv3dwrite_xy_decom(vtavg,2,nfils,'u2_ssm_3d.dat')
  !call sendrecv3dwrite_xy_decom(wtavg,3,nfils,'u2_tot_3d.dat')
  
  call mpi_finalize(ierr)

end program programpod
