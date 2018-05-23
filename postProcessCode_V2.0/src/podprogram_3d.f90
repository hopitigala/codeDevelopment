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

  real*8,allocatable,dimension(:)        :: yp,xp,zp,ypdo,ysdo,ys,eval
  real*8,allocatable,dimension(:)        :: uv_lsm,uv_ssm,uv_tot
  real*8,allocatable,dimension(:,:,:)    :: uv_lsm_3d,uv_ssm_3d,uv_tot_3d
  real*8,allocatable,dimension(:,:,:,:)  :: up,tp
  real*8,allocatable,dimension(:,:,:)    :: wtavg,vtavg,utavg,ttavg,tprime
  real*8,allocatable,dimension(:,:,:)    :: uprime,vprime,wprime,pp
  real*8,allocatable,dimension(:,:,:,:)  :: uprit,vprit,wprit,tprit,phiu,phiv,phiw
  real*8,allocatable,dimension(:,:,:,:)  :: u_lsm,v_lsm,u_ssm,v_ssm,w_ssm,w_lsm
  real*8,allocatable,dimension(:,:)      :: c,evec,coe
  real*8,allocatable,dimension(:,:,:)    :: recnstu,recnstv

  
  integer                              :: timestcount,itime,numtimesteps,iplane,dt
  integer                              :: icorl,x0,t,lsmlim,ireyst,inst,tt,iscfx
  integer                              :: icorlavg,icen,var1,var2,j0,k0
  integer                              :: fexist1,fexist2,fexist3,fexist4,fexist5,fexist6
  real*8                               :: sttime
  character(len=2)                     :: prosnum
  character(len=4)                     :: xloc,time
  character(len=30)                    :: filename
  
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

  allocate(uprit(n1,n2do+1,n3,numtimesteps+1),vprit(n1,n2do+1,n3,numtimesteps+1),&
       wprit(n1,n2do+1,n3,numtimesteps+1),tprit(n1,n2do+1,n3,numtimesteps+1))
  
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
  
     call fluc_3d(uprime,timestcount,uprit)
     call fluc_3d(vprime,timestcount,vprit)
     call fluc_3d(wprime,timestcount,wprit)
     call fluc_3d(tprime,timestcount,tprit)

     
     deallocate(tprime,uprime,vprime,wprime)
     write(*,*)'reading time step',itime,'is done by',mynode
  end do
  allocate(c(timestcount,timestcount))
  
  inquire(file='kernalMatx.dat',exist=fexist1)

  if(fexist1)then
     write(*,*)'kernal matrix file is available reading file'
     call read2DMatx(timestcount,timestcount,'kernalMatx.dat',c)
  else     
     call compkernal_3d(uprit,vprit,wprit,xp,ypdo,zp,c)
  end if

  allocate(evec(timestcount,timestcount),eval(timestcount))
  
  inquire(file='eigenvalues.dat',exist=fexist1)

  if(fexist1)then
     call readvector('eigenvalues.dat',eval)
     call read2DMatx(timestcount,timestcount,'eigenvectors.dat',evec)
  else
     call symmatxeignv(c,'V','U',evec,eval)
  
     call printVector(eval,'eigenvalues.dat')
     call print2DMatx(evec,'eigenvectors.dat')
  end if
  
  deallocate(c)
  inquire(file='u_lsm'//trim(prosnum),exist=fexist1)
  inquire(file='u_ssm'//trim(prosnum),exist=fexist2)
  inquire(file='v_lsm'//trim(prosnum),exist=fexist3)
  inquire(file='v_ssm'//trim(prosnum),exist=fexist4)
  inquire(file='w_lsm'//trim(prosnum),exist=fexist5)
  inquire(file='w_ssm'//trim(prosnum),exist=fexist6)
  
  if (fexist1.and.fexist2.and.fexist3.and.fexist4.and.fexist5.and.fexist6)then
     write(*,*)'all reconstructed velocity fields exist no need to compute POD'
  else
     
  ! allocating arrays to store modes
     allocate(phiu(n1,n2do+1,n3,timestcount),phiv(n1,n2do+1,n3,timestcount),&
          phiw(n1,n2do+1,n3,timestcount))
  ! compute pod modes
     call podmodes_3d(uprit,evec,phiu)
     call podmodes_3d(vprit,evec,phiv)
     call podmodes_3d(wprit,evec,phiw)
  ! normalizing modes
     call norm_modes_3d(xp,ypdo,zp,phiu,phiv,phiw)
  ! checking orthogonality of pod modes
     call checkModeOrthog_3d(xp,ypdo,zp,phiu,phiv,phiw)

     ! compute coefficients
     allocate(coe(timestcount,timestcount))
     call podcoe_3d(xp,ypdo,zp,uprit,vprit,wprit,phiu,phiv,phiw,coe)

     deallocate(uprit,vprit,wprit)
  
     !velocity reconstruction
     ! u_lsm
     allocate(u_lsm(n1,n2do+1,n3,timestcount))
     do t=1,timestcount      
        call recnstfield_3d(phiu,coe,1,lsmlim,t,u_lsm)
     end do
     call write4Darray(u_lsm,'u_lsm'//trim(prosnum),'unformatted')
     deallocate(u_lsm)
     
     !! u_ssm
     allocate(u_ssm(n1,n2do+1,n3,timestcount))
     do t=1,timestcount      
        call recnstfield_3d(phiu,coe,lsmlim+1,timestcount,t,u_ssm)
     end do
     call write4Darray(u_ssm,'u_ssm'//trim(prosnum),'unformatted')
     deallocate(u_ssm,phiu)
     
     !! v_lsm
     allocate(v_lsm(n1,n2do+1,n3,timestcount))
     do t=1,timestcount      
        call recnstfield_3d(phiv,coe,1,lsmlim,t,v_lsm)
     end do
     call write4Darray(v_lsm,'v_lsm'//trim(prosnum),'unformatted')
     deallocate(v_lsm)
     
     !! v_ssm
     allocate(v_ssm(n1,n2do+1,n3,timestcount))
     do t=1,timestcount      
        call recnstfield_3d(phiv,coe,lsmlim+1,timestcount,t,v_ssm)
     end do
     call write4Darray(v_ssm,'v_ssm'//trim(prosnum),'unformatted')
     deallocate(v_ssm,phiv)
     
     !! w_lsm
     allocate(w_lsm(n1,n2do+1,n3,timestcount))
     do t=1,timestcount      
        call recnstfield_3d(phiw,coe,1,lsmlim,t,w_lsm)
     end do
     call write4Darray(w_lsm,'w_lsm'//trim(prosnum),'unformatted')
     deallocate(w_lsm)
     
     !! w_ssm
     allocate(w_ssm(n1,n2do+1,n3,timestcount))
     do t=1,timestcount      
        call recnstfield_3d(phiw,coe,lsmlim+1,timestcount,t,w_ssm)
     end do
     call write4Darray(w_ssm,'w_ssm'//trim(prosnum),'unformatted')
     deallocate(w_ssm,phiw)

  end if


  ! This part of the program compute Reynold stresses
  if(ireyst/=0)then
     allocate(v_lsm(n1,n2do+1,n3,timestcount),u_lsm(n1,n2do+1,n3,timestcount),& 
          u_ssm(n1,n2do+1,n3,timestcount),v_ssm(n1,n2do+1,n3,timestcount))
    

     ! This section of the program computes Reynolds stresses.
     ! When ireyst/=0 Reynolds stresses are computed.
     ! ireyst=1 for w2
     ! ireyst=2 for v2
     ! ireyst=3 for u2
     ! ireyst=4 for uv
     
    allocate(uv_lsm(n2do+1),uv_lsm_3d(n1,n2do+1,n3),& 
          uv_ssm(n2do+1),uv_ssm_3d(n1,n2do+1,n3), & 
          uv_tot(n2do+1),uv_tot_3d(n1,n2do+1,n3))
        
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
           call sendrecv3dwrite(uv_lsm_3d,2,'u2_lsm_3d.dat')
           call sendrecv3dwrite(uv_ssm_3d,4,'u2_ssm_3d.dat')
           call sendrecv3dwrite(uv_tot_3d,6,'u2_tot_3d.dat')
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
           call sendrecv3dwrite(uv_lsm_3d,2,'v2_lsm_3d.dat')
           call sendrecv3dwrite(uv_ssm_3d,4,'v2_ssm_3d.dat')
           call sendrecv3dwrite(uv_tot_3d,6,'v2_tot_3d.dat')
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
           call sendrecv3dwrite(uv_lsm_3d,2,'w2_lsm_3d.dat')
           call sendrecv3dwrite(uv_ssm_3d,4,'w2_ssm_3d.dat')
           call sendrecv3dwrite(uv_tot_3d,6,'w2_tot_3d.dat')
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
           call sendrecv3dwrite(uv_lsm_3d,2,'uv_lsm_3d.dat')
           call sendrecv3dwrite(uv_ssm_3d,4,'uv_ssm_3d.dat')
           call sendrecv3dwrite(uv_tot_3d,6,'uv_tot_3d.dat')
        end if
     end if
     deallocate(u_lsm,u_ssm,uv_lsm,uv_ssm,uv_tot,uv_lsm_3d,uv_ssm_3d,uv_tot_3d,v_lsm,v_ssm)
  end if
  
  if(iscfx/=0)then
     allocate(u_lsm(n1,n2do+1,n3,timestcount),u_ssm(n1,n2do+1,n3,timestcount))
    

     ! This section of the program computes Reynolds stresses.
     ! When ireyst/=0 Reynolds stresses are computed.
     ! ireyst=1 for wt
     ! ireyst=2 for vt
     ! ireyst=3 for ut
     
    allocate(uv_lsm(n2do+1),uv_lsm_3d(n1,n2do+1,n3),& 
          uv_ssm(n2do+1),uv_ssm_3d(n1,n2do+1,n3), & 
          uv_tot(n2do+1),uv_tot_3d(n1,n2do+1,n3))
        
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
          call sendrecv3dwrite(uv_lsm_3d,2,'ut_lsm_3d.dat')
          call sendrecv3dwrite(uv_ssm_3d,4,'ut_ssm_3d.dat')
          call sendrecv3dwrite(uv_tot_3d,6,'ut_tot_3d.dat')
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
          call sendrecv3dwrite(uv_lsm_3d,2,'vt_lsm_3d.dat')
          call sendrecv3dwrite(uv_ssm_3d,4,'vt_ssm_3d.dat')
          call sendrecv3dwrite(uv_tot_3d,6,'vt_tot_3d.dat')
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
          call sendrecv3dwrite(uv_lsm_3d,2,'wt_lsm_3d.dat')
          call sendrecv3dwrite(uv_ssm_3d,4,'wt_ssm_3d.dat')
          call sendrecv3dwrite(uv_tot_3d,6,'wt_tot_3d.dat')
       end if
       deallocate(u_lsm,u_ssm,uv_lsm,uv_ssm,uv_tot,uv_lsm_3d,uv_ssm_3d,uv_tot_3d)
    end if
 end if

  ! This part of the program writes the instantaneous field for different planes

 if(icorl/=0)then
    call tpcorl_3d_podfield(var1,var2,timestcount,x0,j0,k0)
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
        
        
  
