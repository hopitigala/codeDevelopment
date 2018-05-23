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
  real*8,allocatable,dimension(:)      :: uv_lsm_cen,uv_lsm_mid,uv_ssm_cen
  real*8,allocatable,dimension(:)      :: uv_ssm_mid,uv_tot_cen,uv_tot_mid
  real*8,allocatable,dimension(:,:,:,:):: up,tp
  real*8,allocatable,dimension(:,:,:)  :: wtavg,vtavg,utavg,ttavg,tprime
  real*8,allocatable,dimension(:,:,:)  :: uprime,vprime,wprime,pp
  real*8,allocatable,dimension(:,:,:)  :: uprit,vprit,wprit,phiu,phiv,phiw
  real*8,allocatable,dimension(:,:,:)  :: u_lsm,v_lsm,u_ssm,v_ssm
  real*8,allocatable,dimension(:,:)    :: c,evec,recnstu,recnstv,coe

  
  integer                              :: timestcount,itime,numtimesteps,iplane,dt
  integer                              :: icorl,x0,t,lsmlim,ireyst,inst,tt
  integer                              :: icorlavg,icen,var1,var2,j0,k0
  integer                              :: fexist1,fexist2, fexist3, fexist4
  real*8                               :: sttime
  character(len=2)                     :: prosnum
  character(len=4)                     :: xloc,time
  character(len=30)                    :: filename
  
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,mynode,ierr)
  call mpi_comm_size(mpi_comm_world,numprocs,ierr)
  
  call readChannelData()
  call readPostData(sttime,numtimesteps,dt,icorlavg,icen)
  call readpoddata(iplane,x0,lsmlim,ireyst,icorl,var1,var2,j0,k0,inst)
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
          wprit(n2do+1,n3,numtimesteps+1))
  elseif(iplane==2)then
     allocate(uprit(n1,n2do+1,numtimesteps+1),vprit(n1,n2do+1,numtimesteps+1),&
          wprit(n1,n2do+1,numtimesteps+1))
  else
     allocate(uprit(n1,n3,numtimesteps+1),vprit(n1,n3,numtimesteps+1),&
          wprit(n1,n3,numtimesteps+1))
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
     
     deallocate(tprime,uprime,vprime,wprime)
     write(*,*)'reading time step',itime,'is done by',mynode
  end do
  !if(iplane==1)then
  !   inquire(file='kernalMatx_yz.dat',exist=fexist1)
  !elseif(iplane==2)then
  !   inquire(file='kernalMatx_xy.dat',exist=fexist1)
  !elseif(iplane==3)
  !   inquire(file='kernalMatx_xz.dat',exist=fexist1)
  !end if
  !if(fexist1)then
  !   write(*,*)'kernal matrix file is available no need to compute'
  !else
     allocate(c(timestcount,timestcount))
     call compkernal(uprit,vprit,wprit,xp,ypdo,zp,iplane,c)
  

  allocate(evec(timestcount,timestcount),eval(timestcount))
  
  call symmatxeignv(c,'V','U',evec,eval)
  
  if(iplane==1)then
     call printVector(eval,'eigenvalues_xy.dat')
     call print2DMatx(evec,'eigenvectors_xy.dat')
  elseif(iplane==2)then
     call printVector(eval,'eigenvalues_yz.dat')
     call print2DMatx(evec,'eigenvectors_yz.dat')
  else
     call printVector(eval,'eigenvalues_xz.dat')
     call print2DMatx(evec,'eigenvectors_xz.dat')
  end if
  
  deallocate(c)
  ! allocating arrays to store modes
   if(iplane==1)then
     allocate(phiu(n2do+1,n3,timestcount),phiv(n2do+1,n3,timestcount),&
          phiw(n2do+1,n3,timestcount))
  elseif(iplane==2)then
     allocate(phiu(n1,n2do+1,timestcount),phiv(n1,n2do+1,timestcount),&
          phiw(n1,n2do+1,timestcount))
  else
     allocate(phiu(n1,n3,timestcount),phiv(n1,n3,timestcount),&
          phiw(n1,n3,timestcount))
  end if
  ! compute pod modes
  call podModes_2Dplane(uprit,evec,phiu)
  call podModes_2Dplane(vprit,evec,phiv)
  call podModes_2Dplane(wprit,evec,phiw)
  ! normalizing modes
  call norm_modes_2Dplane(xp,ypdo,zp,iplane,phiu,phiv,phiw)
  ! checking orthogonality of pod modes
!  call checkModeOrthog_2Dpln(xp,ypdo,zp,iplane,phiu,phiv,phiw)

  ! compute coefficients
  allocate(coe(timestcount,timestcount))
  call podcoe_2d(xp,ypdo,zp,uprit,vprit,wprit,phiu,phiv,phiw,iplane,coe)

  deallocate(uprit,vprit,wprit)
  
  !velocity reconstruction
  if (iplane==1)then
     inquire(file='u_lsm_xy'//trim(prosnum),exist=fexist1)
  elseif(iplane==2)then
     inquire(file='u_lsm_yz'//trim(prosnum),exist=fexist1)
  else
     inquire(file='u_lsm_xz'//trim(prosnum),exist=fexist1)
  end if
  if(fexist1)then
     write(*,*)'file exist no need to reconstruct'
  else
     if (iplane==1)then
        allocate(u_lsm(n2do+1,n3,timestcount))
     elseif(iplane==2)then
        allocate(u_lsm(n1,n2do+1,timestcount))
     else
        allocate(u_lsm(n1,n3,timestcount))
     end if
     
     do t=1,timestcount      
        call recnstfield(phiu,coe,1,lsmlim,t,u_lsm)
     end do
  !write u_lsm
     if(iplane==1)then
        call write3Darray(u_lsm,'u_lsm_xy'//trim(prosnum),'unformatted')
     elseif(iplane==2)then
        call write3Darray(u_lsm,'u_lsm_yz'//trim(prosnum),'unformatted')
     else
        call write3Darray(u_lsm,'u_lsm_xz'//trim(prosnum),'unformatted')
     end if
  
     deallocate(u_lsm)
  end if

  if (iplane==1)then
     inquire(file='u_ssm_xy'//trim(prosnum),exist=fexist2)
  elseif(iplane==2)then
     inquire(file='u_ssm_yz'//trim(prosnum),exist=fexist2)
  else
     inquire(file='u_ssm_xz'//trim(prosnum),exist=fexist2)
  end if
  if(fexist2)then
     write(*,*)'file exist no need to reconstruct'
  else        
     if (iplane==1)then
        allocate(u_ssm(n2do+1,n3,timestcount))
     elseif(iplane==2)then
        allocate(u_ssm(n1,n2do+1,timestcount))
     else
        allocate(u_ssm(n1,n3,timestcount))
     end if
     
     do t=1,timestcount      
        call recnstfield(phiu,coe,lsmlim+1,numtimesteps,t,u_ssm)
     end do
  
     !write u_ssm
     if(iplane==1)then
        call write3Darray(u_ssm,'u_ssm_xy'//trim(prosnum),'unformatted')
     elseif(iplane==2)then
        call write3Darray(u_ssm,'u_ssm_yz'//trim(prosnum),'unformatted')
     else
        call write3Darray(u_ssm,'u_ssm_xz'//trim(prosnum),'unformatted')
     end if
 
     deallocate(u_ssm)
  end if
  deallocate(phiu)
  
  if (iplane==1)then
     inquire(file='v_lsm_xy'//trim(prosnum),exist=fexist3)
  elseif(iplane==2)then
     inquire(file='v_lsm_yz'//trim(prosnum),exist=fexist3)
  else
     inquire(file='v_lsm_xz'//trim(prosnum),exist=fexist3)
  end if
  if(fexist3)then
     write(*,*)'file exist no need to reconstruct'
  else
     ! compute v_lsm
     if (iplane==1)then
        allocate(v_lsm(n2do+1,n3,timestcount))
     elseif(iplane==2)then
        allocate(v_lsm(n1,n2do+1,timestcount))
     else
        allocate(v_lsm(n1,n3,timestcount))
     end if
     
     do t=1,timestcount      
        call recnstfield(phiv,coe,1,lsmlim,t,v_lsm)
     end do
  
     !write v_lsm
     if(iplane==1)then
        call write3Darray(v_lsm,'v_lsm_xy'//trim(prosnum),'unformatted')
     elseif(iplane==2)then
        call write3Darray(v_lsm,'v_lsm_yz'//trim(prosnum),'unformatted')
     else
        call write3Darray(v_lsm,'v_lsm_xz'//trim(prosnum),'unformatted')
     end if
  
     deallocate(v_lsm)
  end if
  
  !compute v_ssm
  if (iplane==1)then
     inquire(file='v_ssm_xy'//trim(prosnum),exist=fexist4)
  elseif(iplane==2)then
     inquire(file='v_ssm_yz'//trim(prosnum),exist=fexist4)
  else
     inquire(file='v_ssm_xz'//trim(prosnum),exist=fexist4)
  end if
  if(fexist4)then
     write(*,*)'file exist no need to reconstruct'
  else
     if (iplane==1)then
        allocate(v_ssm(n2do+1,n3,timestcount))
     elseif(iplane==2)then
        allocate(v_ssm(n1,n2do+1,timestcount))
     else
        allocate(v_ssm(n1,n3,timestcount))
     end if
     
     do t=1,timestcount      
        call recnstfield(phiv,coe,lsmlim+1,numtimesteps,t,v_ssm)
     end do
  
     !write v_lsm
     if(iplane==1)then
        call write3Darray(v_ssm,'v_ssm_xy'//trim(prosnum),'unformatted')
     elseif(iplane==2)then
        call write3Darray(v_ssm,'v_ssm_yz'//trim(prosnum),'unformatted')
     else
        call write3Darray(v_ssm,'v_ssm_xz'//trim(prosnum),'unformatted')
     end if
     deallocate(v_ssm)
  end if

  deallocate(phiv)

  ! This part of the program compute Reynold stresses
  if(ireyst/=0)then
     if (iplane==1)then
        allocate(v_lsm(n2do+1,n3,timestcount),u_lsm(n2do+1,n3,timestcount),& 
             u_ssm(n2do+1,n3,timestcount),v_ssm(n2do+1,n3,timestcount))
     elseif(iplane==2)then
        allocate(v_lsm(n1,n2do+1,timestcount),u_lsm(n1,n2do+1,timestcount),&
             v_ssm(n1,n2do+1,timestcount),u_ssm(n1,n2do+1,timestcount))
     else
        allocate(v_lsm(n1,n3,timestcount),u_lsm(n1,n3,timestcount),& 
             v_ssm(n1,n3,timestcount),u_ssm(n1,n3,timestcount))
     end if

     ! This section of the program computes Reynolds stresses.
     ! When ireyst/=0 Reynolds stresses are computed.
     ! reyst=1 for u2
     ! reyst=2 for v2
     ! reyst=3 for uv
     
     if(x0<10)then
        write(xloc,'(i1.1)')x0
     elseif(x0<100)then
        write(xloc,'(i2.2)')x0
     elseif(x0<1000)then
        write(xloc,'(i3.3)')x0
     else
        write(xloc,'(i4.4)')x0
     end if
     
     if (iplane==2)then ! code works when iplane==2 only
        call read3Darray('u_lsm_yz'//trim(prosnum),'unformatted',u_lsm)
        call read3Darray('v_lsm_yz'//trim(prosnum),'unformatted',v_lsm)
        call read3Darray('u_ssm_yz'//trim(prosnum),'unformatted',u_ssm)
        call read3Darray('v_ssm_yz'//trim(prosnum),'unformatted',v_ssm)
        
        allocate(uv_lsm_cen(n2do+1),uv_lsm_mid(n2do+1),& 
             uv_ssm_cen(n2do+1),uv_ssm_mid(n2do+1), & 
             uv_tot_cen(n2do+1),uv_tot_mid(n2do+1))
        
        if(ireyst==3)then ! computing uv
           
           call reyst_ssm_lsm(u_lsm,v_lsm,ibs,iplane,uv_lsm_cen,uv_lsm_mid)
           call reyst_ssm_lsm(u_ssm,v_ssm,ibs,iplane,uv_ssm_cen,uv_ssm_mid)
           call reyst_ssm_lsm(u_lsm+u_ssm,v_lsm+v_ssm,ibs,iplane,uv_tot_cen,uv_tot_mid)
           ! write uv 
           
           call senRev1dwrite(uv_lsm_cen,1,'uv_lsm_cen'//trim(xloc)//'.dat')
           call senRev1dwrite(uv_lsm_mid,2,'uv_lsm_mid'//trim(xloc)//'.dat')
           call senRev1dwrite(uv_ssm_cen,3,'uv_ssm_cen'//trim(xloc)//'.dat')
           call senRev1dwrite(uv_ssm_mid,4,'uv_ssm_mid'//trim(xloc)//'.dat')
           call senRev1dwrite(uv_tot_cen,5,'uv_tot_cen'//trim(xloc)//'.dat')
           call senRev1dwrite(uv_tot_mid,6,'uv_tot_mid'//trim(xloc)//'.dat')
        
        elseif(ireyst==2)then ! compute v2
        
           call reyst_ssm_lsm(v_lsm,v_lsm,ibs,iplane,uv_lsm_cen,uv_lsm_mid)
           call reyst_ssm_lsm(v_ssm,v_ssm,ibs,iplane,uv_ssm_cen,uv_ssm_mid)
           call reyst_ssm_lsm(v_lsm+v_ssm,v_lsm+v_ssm,ibs,iplane,uv_tot_cen,uv_tot_mid)
           ! write v2
           
           call senRev1dwrite(uv_lsm_cen,1,'v2_lsm_cen'//trim(xloc)//'.dat')
           call senRev1dwrite(uv_lsm_mid,2,'v2_lsm_mid'//trim(xloc)//'.dat')
           call senRev1dwrite(uv_ssm_cen,3,'v2_ssm_cen'//trim(xloc)//'.dat')
           call senRev1dwrite(uv_ssm_mid,4,'v2_ssm_mid'//trim(xloc)//'.dat')
           call senRev1dwrite(uv_tot_cen,5,'v2_tot_cen'//trim(xloc)//'.dat')
           call senRev1dwrite(uv_tot_mid,6,'v2_tot_mid'//trim(xloc)//'.dat')
        
        elseif(ireyst==1)then !compute u2
        
           call reyst_ssm_lsm(u_lsm,u_lsm,ibs,iplane,uv_lsm_cen,uv_lsm_mid)
           call reyst_ssm_lsm(u_ssm,u_ssm,ibs,iplane,uv_ssm_cen,uv_ssm_mid)
           call reyst_ssm_lsm(u_lsm+u_ssm,u_lsm+u_ssm,ibs,iplane,uv_tot_cen,uv_tot_mid)
           ! write u2
           
           call senRev1dwrite(uv_lsm_cen,1,'u2_lsm_cen'//trim(xloc)//'.dat')
           call senRev1dwrite(uv_lsm_mid,2,'u2_lsm_mid'//trim(xloc)//'.dat')
           call senRev1dwrite(uv_ssm_cen,3,'u2_ssm_cen'//trim(xloc)//'.dat')
           call senRev1dwrite(uv_ssm_mid,4,'u2_ssm_mid'//trim(xloc)//'.dat')
           call senRev1dwrite(uv_tot_cen,5,'u2_tot_cen'//trim(xloc)//'.dat')
           call senRev1dwrite(uv_tot_mid,6,'u2_tot_mid'//trim(xloc)//'.dat')
        end if
        if(mynode==0)then
           call printVector(yp(1:),'ycoord.dat')
        end if
     end if
  end if

  ! This part of the program writes the instantaneous field for different planes
  if(inst==1)then  
     if (iplane==1)then
        allocate(v_lsm(n2do+1,n3,numtimesteps+1),u_lsm(n2do+1,n3,numtimesteps+1),&
             u_ssm(n2do+1,n3,numtimesteps+1),v_ssm(n2do+1,n3,numtimesteps+1))
     elseif(iplane==2)then
        allocate(v_lsm(n1,n2do+1,numtimesteps+1),u_lsm(n1,n2do+1,numtimesteps+1),&
             v_ssm(n1,n2do+1,numtimesteps+1),u_ssm(n1,n2do+1,numtimesteps+1))
     else
        allocate(v_lsm(n1,n3,numtimesteps+1),u_lsm(n1,n3,numtimesteps+1),&
             v_ssm(n1,n3,numtimesteps+1),u_ssm(n1,n3,numtimesteps+1))
     end if
    ! write(prosnum,'(i2.2)')mynode
     
     if (iplane==2)then
        call read3Darray('u_lsm_yz'//trim(prosnum),'unformatted',u_lsm)
        call read3Darray('v_lsm_yz'//trim(prosnum),'unformatted',v_lsm)
        call read3Darray('u_ssm_yz'//trim(prosnum),'unformatted',u_ssm)
        call read3Darray('v_ssm_yz'//trim(prosnum),'unformatted',v_ssm)
        do tt=1,numtimesteps+1
         
           if(tt<10)then
              write(time,'(i1.1)')tt
           elseif(tt<100)then
              write(time,'(i2.2)')tt
           elseif(tt<1000)then
              write(time,'(i3.3)')tt
           else
              write(time,'(i4.4)')tt
           end if
           filename='u_lsm_inst_yz_'//trim(time)//'.dat'
           call sendrecv2dwrite(u_lsm(:,:,tt),1,iplane,filename)
           filename='u_ssm_inst_yz_'//trim(time)//'.dat'
           call sendrecv2dwrite(u_ssm(:,:,tt),2,iplane,filename)
        end do
     elseif(iplane==1)then
        call read3Darray('u_lsm_xy'//trim(prosnum),'unformatted',u_lsm)
        call read3Darray('v_lsm_xy'//trim(prosnum),'unformatted',v_lsm)
        call read3Darray('u_ssm_xy'//trim(prosnum),'unformatted',u_ssm)
        call read3Darray('v_ssm_xy'//trim(prosnum),'unformatted',v_ssm)
        do tt=1,numtimesteps+1

           if(tt<10)then
              write(time,'(i1.1)')tt
           elseif(tt<100)then
              write(time,'(i2.2)')tt
           elseif(tt<1000)then
              write(time,'(i3.3)')tt
           else
              write(time,'(i4.4)')tt
           end if
           filename='u_lsm_inst_xy_'//trim(time)//'.dat'
           call sendrecv2dwrite(u_lsm(:,:,tt),1,iplane,filename)
           filename='u_ssm_inst_xy_'//trim(time)//'.dat'
           call sendrecv2dwrite(u_ssm(:,:,tt),2,iplane,filename)
        end do
     end if
  end if

  if(icorl/=0)then
     call tpcorl_2d_podfield(var1,var2,numtimesteps,x0,j0,k0,iplane)
  end if
  if(ireyst==0)then
     if(mynode==0)then
        call printVector(yp(1:),'ycoord.dat')
        call printVector(xp,'xcoord.dat')
        call printVector(zp,'zcoord.dat')
     end if
  end if
  call mpi_finalize(ierr)
end program programpod
        
        
  
