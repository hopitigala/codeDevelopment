program pod3dintlen
  use mpi
  use mpivariables
  implicit none
  integer                           ::d1,d2,d3,i,j,k
    complex(kind=8),allocatable,dimension(:,:,:):: uprihat,vprihat
    complex(kind=8),allocatable,dimension(:,:)::refary
    real*8, allocatable, dimension(:,:,:)  :: recnstu_new,recnstv_new,recnstw_new,eta,zeta,anifac
    real*8, allocatable, dimension(:,:,:,:):: phiu,phiv,phiw
    real*8, allocatable, dimension(:,:,:,:):: recnstu_old,recnstv_old,recnstw_old
    real*8,allocatable, dimension(:,:):: coe
    real*8, allocatable, dimension(:)      :: yp,xp,zp,ypdo,ysdo,ys,xpdo,anifac_int

    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world,mynode,ierr)
    call mpi_comm_size(mpi_comm_world,numprocs,ierr)

    call readChannelData()
    call readpostdata_3dpod(sttime,numtimesteps,dt,icorlavg,icen,nfils)
    call readpod3ddata(icomppod,irecnst,lsmlim,ireyst,iscfx,inst,icorl,var1,var2,i0,j0,k0,i1toM)
    call prelimcal_3dpod(nfils)

    allocate(xp(n3),yp(0:n2),zp(n1),ys(0:n2),ypdo(0:n2do+1),ysdo(0:n2do+1),xpdo(n3do))

    call readCoordData(xp,yp,zp)
    call intpolcordi_3dpod(yp,xp,ys,ypdo,ysdo,xpdo)
    ! loop for mode count starts here                                                                                                

    nt=numtimesteps+1
    inv_nt=1.0/real(nt)

    allocate(phiu(n1m,n2do+1,n3do,nt),phiv(n1m,n2do+1,n3do,nt),phiw(n1m,n2do+1,n3do,nt))
    write(prosnum,'(i2.2)')mynode

    call read4Darray('/scratch/03390/sdharmar/unperturbed_4pi/postProcess/MPI/pod_new/pod_3d_new/phiu'//trim(prosnum),'unformatted',\
phiu)
    call read4Darray('/scratch/03390/sdharmar/unperturbed_4pi/postProcess/MPI/pod_new/pod_3d_new/phiv'//trim(prosnum),'unformatted',\
phiv)
    call read4Darray('/scratch/03390/sdharmar/unperturbed_4pi/postProcess/MPI/pod_new/pod_3d_new/phiw'//trim(prosnum),'unformatted',\
phiw)

! reading coefficients                                                                                                              \
                                                                                                                                     
    allocate(coe(nt,nt))
    call read2DMatx(nt,nt,'coeMatx.dat',coe)

  allocate(recnstu_old(n1m,n2do+1,n3do,nt),recnstv_old(n1m,n2do+1,n3do,nt),recnstw_old(n1m,n2do+1,n3do,nt))
    do m=1,1000
     if (i1toM==0.or.m==1)then
        recnstu_old=0.
        recnstv_old=0.
        recnstw_old=0.

     end if
     do t=1,nt
        allocate(recnstu_new(n1m,n2do+1,n3do),recnstv_new(n1m,n2do+1,n3do))

        call modalrecnst_3d(phiu,coe,m,t,i1toM,recnstu_old(:,:,:,t),recnstu_new)
        call modalrecnst_3d(phiv,coe,m,t,i1toM,recnstv_old(:,:,:,t),recnstv_new)

        recnstu_old(:,:,:,t)=recnstu_new
        recnstv_old(:,:,:,t)=recnstv_new
        recnstw_old(:,:,:,t)=recnstw_new

        deallocate(recnstu_new,recnstv_new,recnstw_new)
     end do

end program pod3dintlen
