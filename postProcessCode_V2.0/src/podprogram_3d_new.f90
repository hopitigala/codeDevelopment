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
  use fluctuation
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

  
  integer                              :: timestcount,itime,numtimesteps,dt,nfils,m
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
  call readpod3ddata(icomppod,irecnst,lsmlim,ireyst,iscfx,inst,icorl,var1,var2,i0,j0,k0)
  call prelimcal_3dpod(nfils)

  allocate(xp(n3),yp(0:n2),zp(n1),ys(0:n2),ypdo(0:n2do+1),ysdo(0:n2do+1),xpdo(n3do))

  call readCoordData(xp,yp,zp)
  call intpolcordi_3dpod(yp,xp,ys,ypdo,ysdo,xpdo)

  if(icomppod/=0)then
     allocate(wtavg(n1m,n2do+1,n3do),vtavg(n1m,n2do+1,n3do),&
          utavg(n1m,n2do+1,n3do),ttavg(n1m,n2do+1,n3do))
  
     ! Read mean field
     write(prosnum,'(i2.2)')mynode
     call read3Darray('../../tmean/1500ts_3dpod/tmean'//trim(prosnum),'unformatted',ttavg)
     call read3Darray('../../tmean/1500ts_3dpod/wmean'//trim(prosnum),'unformatted',wtavg)
     call read3Darray('../../tmean/1500ts_3dpod/vmean'//trim(prosnum),'unformatted',vtavg)
     call read3Darray('../../tmean/1500ts_3dpod/umean'//trim(prosnum),'unformatted',utavg)

     allocate(uprit(n1m,n2do+1,n3do,numtimesteps+1),vprit(n1m,n2do+1,n3do,numtimesteps+1),&
          wprit(n1m,n2do+1,n3do,numtimesteps+1),tprit(n1m,n2do+1,n3do,numtimesteps+1))
  
     timestcount=0
     do itime=int(sttime),(int(sttime)+(numtimesteps*dt)),dt
        timestcount=timestcount+1
     
        allocate(up(n1m,n2do+1,n3do),vp(n1m,n2do+1,n3do),wp(n1m,n2do+1,n3do),&
             tp(n1m,n2do+1,n3do),pp(n1m,n2do+1,n3do))
     
        call readtemp3dpod(mynode,itime,ypdo,ysdo,up,vp,wp,pp,tp)
        
        allocate(tprime(n1m,n2do+1,n3do),wprime(n1m,n2do+1,n3do),vprime(n1m,n2do+1,n3do),&
             uprime(n1m,n2do+1,n3do))

        call fluct3Dmean(tp,ttavg,tprime)
        call fluct3Dmean(wp,wtavg,wprime)
        call fluct3Dmean(vp,vtavg,vprime)
        call fluct3Dmean(up,utavg,uprime)
     
        deallocate(up,vp,wp,tp,pp)
  
        call fluc_3d(uprime,timestcount,uprit)
        call fluc_3d(vprime,timestcount,vprit)
        call fluc_3d(wprime,timestcount,wprit)
        call fluc_3d(tprime,timestcount,tprit)

     
        deallocate(tprime,uprime,vprime,wprime)
        write(*,*)'reading time step',itime,'is done by',mynode
     end do
     write(*,*)'number of timesteps = ',timestcount
     
     allocate(c(timestcount,timestcount))
  
     inquire(file='kernalMatx.dat',exist=fexist1)
     
     if(fexist1)then
        write(*,*)'kernal matrix file is available reading file'
        call read2DMatx(timestcount,timestcount,'kernalMatx.dat',c)
     else
        write(*,*)'kernal matrix file is not available. computing'
        call compkernal_3d(uprit,vprit,wprit,xpdo,ypdo,zp,c)
     end if

     allocate(evec(timestcount,timestcount),eval(timestcount))
     
     inquire(file='eigenvalues.dat',exist=fexist1)

     if(fexist1)then
        write(*,*)'eigenvalue file is available reading file'
        call readvector('eigenvalues.dat',eval)
        call read2DMatx(timestcount,timestcount,'eigenvectors.dat',evec)
     else
        write(*,*)'eigenvalue file is not available. Computing'
        call symmatxeignv(c,'V','U',evec,eval)
        if (mynode==0)then
           call printVector(eval,'eigenvalues.dat')
           call print2DMatx(evec,'eigenvectors.dat')
        end if
     end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Computing POD modes and coefficients.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     deallocate(c)
  !if (icomppod/=0)then
     inquire(file='/scratch/03390/sdharmar/unperturbed_4pi/postProcess/MPI/pod_new/pod_3d_new/phiu'//trim(prosnum),exist=fexist1)
     inquire(file='/scratch/03390/sdharmar/unperturbed_4pi/postProcess/MPI/pod_new/pod_3d_new/phiv'//trim(prosnum),exist=fexist2)
     inquire(file='/scratch/03390/sdharmar/unperturbed_4pi/postProcess/MPI/pod_new/pod_3d_new/phiw'//trim(prosnum),exist=fexist3)
     
     if (fexist1.and.fexist2.and.fexist3)then
        write(*,*)'all modes exist no need to compute POD reading files'
        allocate(phiu(n1m,n2do+1,n3do,timestcount),phiv(n1m,n2do+1,n3do,timestcount),&
             phiw(n1m,n2do+1,n3do,timestcount))

        call read4Darray('/scratch/03390/sdharmar/unperturbed_4pi/postProcess/MPI/pod_new/pod_3d_new/phiu'//trim(prosnum),'unformatted',phiu)
        call read4Darray('/scratch/03390/sdharmar/unperturbed_4pi/postProcess/MPI/pod_new/pod_3d_new/phiv'//trim(prosnum),'unformatted',phiv)
        call read4Darray('/scratch/03390/sdharmar/unperturbed_4pi/postProcess/MPI/pod_new/pod_3d_new/phiw'//trim(prosnum),'unformatted',phiw)
     else
        ! allocating arrays to store modes
        allocate(phiu(n1m,n2do+1,n3do,timestcount),phiv(n1m,n2do+1,n3do,timestcount),&
             phiw(n1m,n2do+1,n3do,timestcount))
        ! compute pod modes
        call podmodes_3d(uprit,evec,phiu)
        call podmodes_3d(vprit,evec,phiv)
        call podmodes_3d(wprit,evec,phiw)
        
        write(*,*)'pod modes were computed'
        ! normalizing modes
        inquire(file='normSqure.dat',exist=fexist1)
        if(fexist1)then
           write(*,*)'normSqure.dat is available.Do not want to compute. Reading'
           allocate(norm2(timestcount))
           call readvector('normSqure.dat',norm2)
           do m=1,timestcount
              norminv=1/sqrt(norm2(m))
              phiu(:,:,:,m)=phiu(:,:,:,m)*norminv
              phiv(:,:,:,m)=phiv(:,:,:,m)*norminv
              phiw(:,:,:,m)=phiw(:,:,:,m)*norminv
           end do
           deallocate(norm2)
           write(*,*)'modes were normalized'
        else
           write(*,*)'normSqure.dat is not available. computing'
           call norm_modes_3d(xpdo,ypdo,zp,phiu,phiv,phiw)
           write(*,*)'modes were normalized'
        end if

        call write4Darray(phiu,'phiu'//trim(prosnum),'unformatted')
        call write4Darray(phiv,'phiv'//trim(prosnum),'unformatted')
        call write4Darray(phiw,'phiw'//trim(prosnum),'unformatted')
        

        ! checking orthogonality of pod modes
!        call checkModeOrthog_3d(xp,ypdo,zp,phiu,phiv,phiw)
     end if
     
     inquire(file='coeMatx.dat',exist=fexist1)
     if (fexist1)then
        allocate(coe(timestcount,timestcount))
        write(*,*)'coefficent file exists no need to compute. reading'
        call read2DMatx(timestcount,timestcount,'coeMatx.dat',coe)
     else
        ! compute coefficients
        allocate(coe(timestcount,timestcount))
        call podcoe_3d(xpdo,ypdo,zp,uprit,vprit,wprit,phiu,phiv,phiw,coe)
        if(mynode==0)then
           call print2DMatx(coe,'coeMatx.dat')
           write(*,*)'coefficients are computed'
        end if
     end if
     deallocate(uprit,vprit,wprit,evec,eval)
  end if

  
  !write(*,*)'irecnst before = ',irecnst

  if(irecnst/=0)then
     !velocity reconstruction
     ! u_lsm
     inquire(file='u_lsm'//trim(prosnum),exist=fexist1)
     
     if (fexist1)then
        write(*,*)'u_lsm available no need to compute'

     else
        
        write(*,*)'u_lsm computing'
        allocate(u_lsm(n1m,n2do+1,n3do,timestcount))
        do t=1,timestcount      
           call recnstfield_3d(phiu,coe,1,lsmlim,t,u_lsm)
        end do
        call write4Darray(u_lsm,'u_lsm'//trim(prosnum),'unformatted')
        deallocate(u_lsm)
     end if
     
     !! u_ssm

     inquire(file='u_ssm'//trim(prosnum),exist=fexist2)
     if (fexist2)then
        write(*,*)'u_ssm available no need to compute'

     else
        write(*,*)'u_ssm computing'
        allocate(u_ssm(n1m,n2do+1,n3do,timestcount))
        do t=1,timestcount      
           call recnstfield_3d(phiu,coe,lsmlim+1,timestcount,t,u_ssm)
        end do
        call write4Darray(u_ssm,'u_ssm'//trim(prosnum),'unformatted')
        deallocate(u_ssm,phiu)
     end if

     !! v_lsm
     inquire(file='v_lsm'//trim(prosnum),exist=fexist3)
     if (fexist3)then
        write(*,*)'v_lsm available no need to compute'

     else
        allocate(v_lsm(n1m,n2do+1,n3do,timestcount))
        do t=1,timestcount      
           call recnstfield_3d(phiv,coe,1,lsmlim,t,v_lsm)
        end do
        call write4Darray(v_lsm,'v_lsm'//trim(prosnum),'unformatted')
        deallocate(v_lsm)
     end if

     !! v_ssm
     
     inquire(file='v_ssm'//trim(prosnum),exist=fexist4)
     if (fexist4)then
        write(*,*)'v_ssm available no need to compute'

     else
        allocate(v_ssm(n1m,n2do+1,n3do,timestcount))
        do t=1,timestcount      
           call recnstfield_3d(phiv,coe,lsmlim+1,timestcount,t,v_ssm)
        end do
        call write4Darray(v_ssm,'v_ssm'//trim(prosnum),'unformatted')
        deallocate(v_ssm,phiv)
     end if

     !! w_lsm
     inquire(file='w_lsm'//trim(prosnum),exist=fexist5)
     if (fexist5)then
        write(*,*)'w_lsm available no need to compute'
        
     else
        allocate(w_lsm(n1m,n2do+1,n3do,timestcount))
        do t=1,timestcount      
           call recnstfield_3d(phiw,coe,1,lsmlim,t,w_lsm)
        end do
        call write4Darray(w_lsm,'w_lsm'//trim(prosnum),'unformatted')
        deallocate(w_lsm)
     end if

     !! w_ssm
     inquire(file='w_ssm'//trim(prosnum),exist=fexist6)
     if (fexist6)then
        write(*,*)'w_ssm available no need to compute'

     else   
        allocate(w_ssm(n1m,n2do+1,n3do,timestcount))
        do t=1,timestcount      
           call recnstfield_3d(phiw,coe,lsmlim+1,timestcount,t,w_ssm)
        end do
        call write4Darray(w_ssm,'w_ssm'//trim(prosnum),'unformatted')
        deallocate(w_ssm,phiw)
     end if
  end if


  ! This part of the program compute Reynold stresses
  if(ireyst/=0)then
     allocate(v_lsm(n1m,n2do+1,n3do,timestcount),u_lsm(n1m,n2do+1,n3do,timestcount),& 
          u_ssm(n1m,n2do+1,n3do,timestcount),v_ssm(n1m,n2do+1,n3do,timestcount))
    

     ! This section of the program computes Reynolds stresses.
     ! When ireyst/=0 Reynolds stresses are computed.
     ! ireyst=1 for w2
     ! ireyst=2 for v2
     ! ireyst=3 for u2
     ! ireyst=4 for uv
     
    allocate(uv_lsm(n2do+1),uv_lsm_3d(n1m,n2do+1,n3do),& 
          uv_ssm(n2do+1),uv_ssm_3d(n1m,n2do+1,n3do), & 
          uv_tot(n2do+1),uv_tot_3d(n1m,n2do+1,n3do))
        
     if(ireyst==3)then ! computing u2
        
        call read4Darray('u_lsm'//trim(prosnum),'unformatted',u_lsm)
        call read4Darray('u_ssm'//trim(prosnum),'unformatted',u_ssm)
        
        
        call reyst_ssm_lsm_3d(u_lsm,u_lsm,ibs,uv_lsm,uv_lsm_3d)
        call reyst_ssm_lsm_3d(u_ssm,u_ssm,ibs,uv_ssm,uv_ssm_3d)
        call reyst_ssm_lsm_3d(u_lsm+u_ssm,u_lsm+u_ssm,ibs,uv_tot,uv_tot_3d)
        ! write u2
        if(ibs==0)then  
           call senRev1dwrite(uv_lsm,1,'u2_lsm.dat')
           call senRev1dwrite(uv_ssm,3,'u2_ssm.dat')
           call senRev1dwrite(uv_tot,5,'u2_tot.dat')
        else
           call sendrecv3dwrite_xy_decom(uv_lsm_3d,2,'u2_lsm_3d.dat')
           call sendrecv3dwrite_xy_decom(uv_ssm_3d,4,'u2_ssm_3d.dat')
           call sendrecv3dwrite_xy_decom(uv_tot_3d,6,'u2_tot_3d.dat')
        end if

     elseif(ireyst==2)then ! compute v2
        call read4Darray('v_lsm'//trim(prosnum),'unformatted',u_lsm)
        call read4Darray('v_ssm'//trim(prosnum),'unformatted',u_ssm)
        
        
        call reyst_ssm_lsm_3d(u_lsm,u_lsm,ibs,uv_lsm,uv_lsm_3d)
        call reyst_ssm_lsm_3d(u_ssm,u_ssm,ibs,uv_ssm,uv_ssm_3d)
        call reyst_ssm_lsm_3d(u_lsm+u_ssm,u_lsm+u_ssm,ibs,uv_tot,uv_tot_3d)
        ! write v2
        if(ibs==0)then  
           call senRev1dwrite(uv_lsm,1,'v2_lsm.dat')
           call senRev1dwrite(uv_ssm,3,'v2_ssm.dat')
           call senRev1dwrite(uv_tot,5,'v2_tot.dat')
        else
           call sendrecv3dwrite_xy_decom(uv_lsm_3d,2,'v2_lsm_3d.dat')
           call sendrecv3dwrite_xy_decom(uv_ssm_3d,4,'v2_ssm_3d.dat')
           call sendrecv3dwrite_xy_decom(uv_tot_3d,6,'v2_tot_3d.dat')
        end if
        
     elseif(ireyst==1)then !compute w2
        
        call read4Darray('w_lsm'//trim(prosnum),'unformatted',u_lsm)
        call read4Darray('w_ssm'//trim(prosnum),'unformatted',u_ssm)
        
        
        call reyst_ssm_lsm_3d(u_lsm,u_lsm,ibs,uv_lsm,uv_lsm_3d)
        call reyst_ssm_lsm_3d(u_ssm,u_ssm,ibs,uv_ssm,uv_ssm_3d)
        call reyst_ssm_lsm_3d(u_lsm+u_ssm,u_lsm+u_ssm,ibs,uv_tot,uv_tot_3d)
        ! write w2
        if(ibs==0)then  
           call senRev1dwrite(uv_lsm,1,'w2_lsm.dat')
           call senRev1dwrite(uv_ssm,3,'w2_ssm.dat')
           call senRev1dwrite(uv_tot,5,'w2_tot.dat')
        else
           call sendrecv3dwrite_xy_decom(uv_lsm_3d,2,'w2_lsm_3d.dat')
           call sendrecv3dwrite_xy_decom(uv_ssm_3d,4,'w2_ssm_3d.dat')
           call sendrecv3dwrite_xy_decom(uv_tot_3d,6,'w2_tot_3d.dat')
        end if

     elseif(ireyst==4)then !compute uv
        
        call read4Darray('u_lsm'//trim(prosnum),'unformatted',u_lsm)
        call read4Darray('v_lsm'//trim(prosnum),'unformatted',v_lsm)
        call read4Darray('u_ssm'//trim(prosnum),'unformatted',u_ssm)
        call read4Darray('v_ssm'//trim(prosnum),'unformatted',v_ssm)
        
        call reyst_ssm_lsm_3d(u_lsm,v_lsm,ibs,uv_lsm,uv_lsm_3d)
        call reyst_ssm_lsm_3d(u_ssm,v_ssm,ibs,uv_ssm,uv_ssm_3d)
        call reyst_ssm_lsm_3d(u_lsm+u_ssm,v_lsm+v_ssm,ibs,uv_tot,uv_tot_3d)
        ! write uv
        if(ibs==0)then  
           call senRev1dwrite(uv_lsm,1,'uv_lsm.dat')
           call senRev1dwrite(uv_ssm,3,'uv_ssm.dat')
           call senRev1dwrite(uv_tot,5,'uv_tot.dat')
        else
           call sendrecv3dwrite_xy_decom(uv_lsm_3d,2,'uv_lsm_3d.dat')
           call sendrecv3dwrite_xy_decom(uv_ssm_3d,4,'uv_ssm_3d.dat')
           call sendrecv3dwrite_xy_decom(uv_tot_3d,6,'uv_tot_3d.dat')
        end if
     end if
     deallocate(u_lsm,u_ssm,uv_lsm,uv_ssm,uv_tot,uv_lsm_3d,uv_ssm_3d,uv_tot_3d,v_lsm,v_ssm)
  end if
  
  if(iscfx/=0)then
     allocate(u_lsm(n1m,n2do+1,n3do,timestcount),u_ssm(n1m,n2do+1,n3do,timestcount))
    

     ! This section of the program computes Reynolds stresses.
     ! When ireyst/=0 Reynolds stresses are computed.
     ! ireyst=1 for wt
     ! ireyst=2 for vt
     ! ireyst=3 for ut
     
    allocate(uv_lsm(n2do+1),uv_lsm_3d(n1m,n2do+1,n3do),& 
          uv_ssm(n2do+1),uv_ssm_3d(n1m,n2do+1,n3do), & 
          uv_tot(n2do+1),uv_tot_3d(n1m,n2do+1,n3do))
        
    if(iscfx==3)then ! computing ut
        
       call read4Darray('u_lsm'//trim(prosnum),'unformatted',u_lsm)
       call read4Darray('u_ssm'//trim(prosnum),'unformatted',u_ssm)
       
       
       call reyst_ssm_lsm_3d(u_lsm,tprit,ibs,uv_lsm,uv_lsm_3d)
       call reyst_ssm_lsm_3d(u_ssm,tprit,ibs,uv_ssm,uv_ssm_3d)
       call reyst_ssm_lsm_3d(u_lsm+u_ssm,tprit,ibs,uv_tot,uv_tot_3d)
        ! write u2
       if(ibs==0)then  
          call senRev1dwrite(uv_lsm,1,'ut_lsm.dat')
          call senRev1dwrite(uv_ssm,3,'ut_ssm.dat')
          call senRev1dwrite(uv_tot,5,'ut_tot.dat')
       else
          call sendrecv3dwrite_xy_decom(uv_lsm_3d,2,'ut_lsm_3d.dat')
          call sendrecv3dwrite_xy_decom(uv_ssm_3d,4,'ut_ssm_3d.dat')
          call sendrecv3dwrite_xy_decom(uv_tot_3d,6,'ut_tot_3d.dat')
       end if

    elseif(iscfx==2)then ! compute vt
       call read4Darray('v_lsm'//trim(prosnum),'unformatted',u_lsm)
       call read4Darray('v_ssm'//trim(prosnum),'unformatted',u_ssm)
        
        
       call reyst_ssm_lsm_3d(u_lsm,tprit,ibs,uv_lsm,uv_lsm_3d)
       call reyst_ssm_lsm_3d(u_ssm,tprit,ibs,uv_ssm,uv_ssm_3d)
       call reyst_ssm_lsm_3d(u_lsm+u_ssm,tprit,ibs,uv_tot,uv_tot_3d)
       ! write v2
       if(ibs==0)then  
          call senRev1dwrite(uv_lsm,1,'vt_lsm.dat')
          call senRev1dwrite(uv_ssm,3,'vt_ssm.dat')
          call senRev1dwrite(uv_tot,5,'vt_tot.dat')
       else
          call sendrecv3dwrite_xy_decom(uv_lsm_3d,2,'vt_lsm_3d.dat')
          call sendrecv3dwrite_xy_decom(uv_ssm_3d,4,'vt_ssm_3d.dat')
          call sendrecv3dwrite_xy_decom(uv_tot_3d,6,'vt_tot_3d.dat')
       end if
        
    elseif(iscfx==1)then !compute wt
        
       call read4Darray('w_lsm'//trim(prosnum),'unformatted',u_lsm)
       call read4Darray('w_ssm'//trim(prosnum),'unformatted',u_ssm)
       
        
       call reyst_ssm_lsm_3d(u_lsm,tprit,ibs,uv_lsm,uv_lsm_3d)
       call reyst_ssm_lsm_3d(u_ssm,tprit,ibs,uv_ssm,uv_ssm_3d)
       call reyst_ssm_lsm_3d(u_lsm+u_ssm,tprit,ibs,uv_tot,uv_tot_3d)
       ! write w2
       if(ibs==0)then  
          call senRev1dwrite(uv_lsm,1,'wt_lsm.dat')
          call senRev1dwrite(uv_ssm,3,'wt_ssm.dat')
          call senRev1dwrite(uv_tot,5,'wt_tot.dat')
       else
          call sendrecv3dwrite_xy_decom(uv_lsm_3d,2,'wt_lsm_3d.dat')
          call sendrecv3dwrite_xy_decom(uv_ssm_3d,4,'wt_ssm_3d.dat')
          call sendrecv3dwrite_xy_decom(uv_tot_3d,6,'wt_tot_3d.dat')
       end if
       deallocate(u_lsm,u_ssm,uv_lsm,uv_ssm,uv_tot,uv_lsm_3d,uv_ssm_3d,uv_tot_3d)
    end if
 end if

  ! This part of the program writes the instantaneous field for different planes

 if(icorl/=0)then
    call tpcorl_3d_podfield(var1,var2,'l',nfils,numtimesteps+1,i0,j0,k0)
    call tpcorl_3d_podfield(var1,var2,'s',nfils,numtimesteps+1,i0,j0,k0)
 end if

 !! wirte coodinate files if they are not in the directory
 if(mynode==0)then
    inquire(file='ycoord.dat',exist=fexist1)
    if(fexist1)then
       write(*,*)'ycoord.dat exists and do not need to write'
    else
       call printVector(yp(1:),'ycoord.dat')
    end if
    inquire(file='xcoord.dat',exist=fexist1)
    if(fexist1)then
       write(*,*)'xcoord.dat exists and do not need to write'
    else
       call printVector(xp,'xcoord.dat')
    end if
    inquire(file='zcoord.dat',exist=fexist1)
    if(fexist1)then
       write(*,*)'zcoord.dat exists and do not need to write'
    else
       call printVector(zp,'zcoord.dat')
    end if
 end if
 call mpi_finalize(ierr)
end program programpod
        
        
  
