program podanisotropyfactor
  use mpi
  use mpivariables
  use mainparameters
  use readdata
  use channhdvariables
  use prelimcalvar
  use writedata
  use anisotropyroutines
  use comppod3d
  use integration

  implicit none
  real*8, allocatable, dimension(:)      :: yp,xp,zp,ypdo,ysdo,ys,xpdo,anifac_int
  real*8, allocatable, dimension(:,:)    :: coe
  real*8, allocatable, dimension(:,:,:)  :: uu_m, vv_m, ww_m, uv_m, uw_m, vw_m
  real*8, allocatable, dimension(:,:,:)  :: recnstu_new,recnstv_new,recnstw_new,eta,zeta,anifac
  real*8, allocatable, dimension(:,:,:,:):: phiu,phiv,phiw
  real*8, allocatable, dimension(:,:,:,:):: recnstu_old,recnstv_old,recnstw_old

  integer  :: nt,numtimesteps,t,m,dt,icorlavg,icen,nfils,i,j,k
  integer  :: icomppod,irecnst,lsmlim,ireyst,iscfx,inst,icorl,var1,var2,i0,j0,k0,i1toM
  real*8   :: inv_nt,prosum,sttime,allsum,summ
  character(len=2)::prosnum
  character(len=4)::noofmodes

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

  nt=numtimesteps+1
  inv_nt=1.0/real(nt)
  
  allocate(phiu(n1m,n2do+1,n3do,nt),phiv(n1m,n2do+1,n3do,nt),phiw(n1m,n2do+1,n3do,nt))
  write(prosnum,'(i2.2)')mynode

  call read4Darray('/scratch/03390/sdharmar/unperturbed_4pi/postProcess/MPI/pod_new/pod_3d_new/phiu'//trim(prosnum),'unformatted',phiu)
  call read4Darray('/scratch/03390/sdharmar/unperturbed_4pi/postProcess/MPI/pod_new/pod_3d_new/phiv'//trim(prosnum),'unformatted',phiv)
  call read4Darray('/scratch/03390/sdharmar/unperturbed_4pi/postProcess/MPI/pod_new/pod_3d_new/phiw'//trim(prosnum),'unformatted',phiw)

! reading coefficients
  allocate(coe(nt,nt))
  call read2DMatx(nt,nt,'coeMatx.dat',coe)

  allocate(recnstu_old(n1m,n2do+1,n3do,nt),recnstv_old(n1m,n2do+1,n3do,nt),recnstw_old(n1m,n2do+1,n3do,nt))
  allocate(anifac_int(nt))
  do m=1,1000

     allocate(uu_m(n1m,n2do+1,n3do),vv_m(n1m,n2do+1,n3do),ww_m(n1m,n2do+1,n3do),&
          uv_m(n1m,n2do+1,n3do),uw_m(n1m,n2do+1,n3do),vw_m(n1m,n2do+1,n3do))

     uu_m=0.0
     vv_m=0.0
     ww_m=0.0
     uw_m=0.0
     vw_m=0.0
     uv_m=0.0
     
     if (i1toM==0.or.m==1)then
        recnstu_old=0.
        recnstv_old=0.
        recnstw_old=0.
        
     end if


     do t=1,nt
        allocate(recnstu_new(n1m,n2do+1,n3do),recnstv_new(n1m,n2do+1,n3do),recnstw_new(n1m,n2do+1,n3do))
        
        call modalrecnst_3d(phiu,coe,m,t,i1toM,recnstu_old(:,:,:,t),recnstu_new)
        call modalrecnst_3d(phiv,coe,m,t,i1toM,recnstv_old(:,:,:,t),recnstv_new)
        call modalrecnst_3d(phiw,coe,m,t,i1toM,recnstw_old(:,:,:,t),recnstw_new)
        
        uu_m=uu_m+recnstu_new*recnstu_new
        vv_m=vv_m+recnstv_new*recnstv_new
        ww_m=ww_m+recnstw_new*recnstw_new
        uw_m=uw_m+recnstu_new*recnstw_new
        uv_m=uv_m+recnstu_new*recnstv_new
        vw_m=vw_m+recnstv_new*recnstw_new
        
        recnstu_old(:,:,:,t)=recnstu_new
        recnstv_old(:,:,:,t)=recnstv_new
        recnstw_old(:,:,:,t)=recnstw_new
       
        deallocate(recnstu_new,recnstv_new,recnstw_new)
     end do
     
     uu_m=uu_m*inv_nt
     vv_m=vv_m*inv_nt
     ww_m=ww_m*inv_nt
     uw_m=uw_m*inv_nt
     uv_m=uv_m*inv_nt
     vw_m=vw_m*inv_nt
     
     allocate(eta(n1m,n2do+1,n3do),zeta(n1m,n2do+1,n3do))
     call reystinv(uu_m,uv_m,uw_m,vv_m,vw_m,ww_m,eta,zeta)
     deallocate(uu_m,uv_m,uw_m,vv_m,vw_m,ww_m)

     allocate(anifac(n1m,n2do+1,n3do))
     anifac=1-27.0*eta+54.0*zeta
     deallocate(eta,zeta)


!    if(m<10)then
!       write(noofmodes,'(i1.1)')m
!    elseif(m<100)then
!       write(noofmodes,'(i2.2)')m
!    elseif(m<1000)then
!       write(noofmodes,'(i3.3)')m
!    else
!       write(noofmodes,'(i4.4)')m
!    end if
!    prosum=0.
    
!    if (mynode<numprocs-4) then
!       do k=1,n3do
!          do j=1,n2do
!             do i=1,n1m
!                prosum=prosum+anifac(i,j,k)
!             end do
!          end do
!       end do
!    else
!       do k=1,n3do
!          do j=1,n2do+1
!             do i=1,n1m
!                prosum=prosum+anifac(i,j,k)
!             end do
!          end do
!       end do
!    end if
 !     if(i1toM/=0)then
 !       call sendrecv3dwrite_xy_decom(anifac,1,nfils,'anifac_for_'//trim(noofmodes)//'modes.dat',summ)
 !    else
 !       call sendrecv3dwrite_xy_decom(anifac,2,nfils,'anifac_for_mode_'//trim(noofmodes)//'.dat',summ)
 !    end if
!     if(i1toM/=0)then
    call sendrecv3dwrite_xy_decom(anifac,1,nfils,summ)
!     else
!        call sendrecv3dwrite_xy_decom(anifac,2,nfils,summ)
!     end if

!     call send_recv_elemntadd_3dary_xy_decom(anifac,1,nfils,summ)
    anifac_int(m)=summ
!     call trapezoidal3d(anifac,n1m,n2do+1,n3do,xpdo(1),xpdo(n3do),ypdo,zp(1),zp(n1m),prosum)
     deallocate(anifac)
!     call mpi_allreduce(prosum,allsum,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
!     anifac_int(m)=allsum/real(n1m*(nfils*(n2do)+1)*filstopros*n3do)
     
     if (mynode==0) then        
        write(*,*)'mode',m,'is done'
     end if
  end do
  
  if(mynode==0)then
     if (i1toM==0)then
        call printVector(anifac_int,'integted_aniso_fac.dat')
     else
        call printVector(anifac_int,'integted_aniso_fac_1toM.dat')
     end if
  end if

  deallocate(anifac_int,recnstu_old,recnstv_old,recnstw_old)
end program podanisotropyfactor
