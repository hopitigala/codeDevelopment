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
  character(len=4)                     :: xloc,time
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

     deallocate(tprime,uprime,vprime,wprime)
     write(*,*)'reading time step',itime,'is done by',mynode
  end do

!! computing Kernel Matrix for eigenvalue computation.
  allocate(c(timestcount,timestcount))

  if(iplane==1)then
     inquire(file='kernalMatx_xy.dat',exist=fexist1)
  elseif(iplane==2)then
     inquire(file='kernalMatx_yz.dat',exist=fexist1)
  elseif(iplane==3)then
     inquire(file='kernalMatx_xz.dat',exist=fexist1)
  end if
  if(fexist1)then
     write(*,*)'kernal matrix file is available reading file'
     if(iplane==1)then
        call read2DMatx(timestcount,timestcount,'kernalMatx_xy.dat',c)
     elseif (iplane==2)then
        call read2DMatx(timestcount,timestcount,'kernalMatx_yz.dat',c)
     elseif (iplane==3)then
        call read2DMatx(timestcount,timestcount,'kernalMatx_xz.dat',c)
     end if
  else     
     call compkernal(uprit,vprit,wprit,xp,ypdo,zp,iplane,c)
  end if

  allocate(evec(timestcount,timestcount),eval(timestcount))
  if(iplane==1)then
     inquire(file='eigenvectors_xy.dat',exist=fexist1)
     inquire(file='eigenvalues_xy.dat',exist=fexist2)
  elseif(iplane==2)then
     inquire(file='eigenvectors_yz.dat',exist=fexist1)
     inquire(file='eigenvalues_yz.dat',exist=fexist2)
  elseif(iplane==3)then
     inquire(file='eigenvectors_xz.dat',exist=fexist1)
     inquire(file='eigenvalues_xz.dat',exist=fexist2)
  end if
  if(fexist1.and.fexist2)then
     write(*,*)'eigenvalues and eigenvectors are available and reading from files'
     if(iplane==1)then
        call read2DMatx(timestcount,timestcount,'eigenvectors_xy.dat',evec)
        call readvector('eigenvalues_xy.dat',eval)
     elseif (iplane==2)then
        call read2DMatx(timestcount,timestcount,'eigenvectors_yz.dat',evec)
        call readvector('eigenvalues_yz.dat',eval)
     elseif (iplane==3)then
        call read2DMatx(timestcount,timestcount,'eigenvectors_xz.dat',evec)
        call readvector('eigenvalues_xz.dat',eval)
     end if
  else     
     call symmatxeignv(c,'V','U',evec,eval)
  
     if(iplane==1)then
        call printVector(eval,'eigenvalues_xy.dat')
        call print2DMatx(evec,'eigenvectors_xy.dat')
     elseif(iplane==2)then
        call printVector(eval,'eigenvalues_yz.dat')
        call print2DMatx(evec,'eigenvectors_yz.dat')
     elseif(iplane==3)then
        call printVector(eval,'eigenvalues_xz.dat')
        call print2DMatx(evec,'eigenvectors_xz.dat')
     end if
  end if
  
  deallocate(c)
  ! allocating arrays to store modes
  if (iplane==1)then
     inquire(file='u_lsm_xy'//trim(prosnum),exist=fexist1)
     inquire(file='u_ssm_xy'//trim(prosnum),exist=fexist2)
     inquire(file='v_lsm_xy'//trim(prosnum),exist=fexist3)
     inquire(file='v_ssm_xy'//trim(prosnum),exist=fexist4)
     inquire(file='w_lsm_xy'//trim(prosnum),exist=fexist5)
     inquire(file='w_ssm_xy'//trim(prosnum),exist=fexist6)
  elseif(iplane==2)then
     inquire(file='u_lsm_yz'//trim(prosnum),exist=fexist1)
     inquire(file='u_ssm_yz'//trim(prosnum),exist=fexist2)
     inquire(file='v_lsm_yz'//trim(prosnum),exist=fexist3)
     inquire(file='v_ssm_yz'//trim(prosnum),exist=fexist4)
     inquire(file='w_lsm_yz'//trim(prosnum),exist=fexist5)
     inquire(file='w_ssm_yz'//trim(prosnum),exist=fexist6)
  elseif(iplane==3)then
     inquire(file='u_lsm_yz'//trim(prosnum),exist=fexist1)
     inquire(file='u_ssm_yz'//trim(prosnum),exist=fexist2)
     inquire(file='v_lsm_yz'//trim(prosnum),exist=fexist3)
     inquire(file='v_ssm_yz'//trim(prosnum),exist=fexist4)
     inquire(file='w_lsm_yz'//trim(prosnum),exist=fexist5)
     inquire(file='w_ssm_yz'//trim(prosnum),exist=fexist6)
  end if
  if (fexist1.and.fexist2.and.fexist3.and.fexist4.and.fexist5.and.fexist6)then
     write(*,*)'all reconstructed velocity fields exist no need to compute POD'
  else
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
  !if (iplane==1)then
  !   inquire(file='u_lsm_xy'//trim(prosnum),exist=fexist1)
  !elseif(iplane==2)then
  !   inquire(file='u_lsm_yz'//trim(prosnum),exist=fexist1)
  !else
  !   inquire(file='u_lsm_xz'//trim(prosnum),exist=fexist1)
  !end if
  !if(fexist1)then
  !   write(*,*)'file exist no need to reconstruct'
  !else
     if (iplane==1)then
        allocate(u_lsm(n2do+1,n3,timestcount))
     elseif(iplane==2)then
        allocate(u_lsm(n1,n2do+1,timestcount))
     else
        allocate(u_lsm(n1,n3,timestcount))
     end if
     
     !do t=1,timestcount      
        call recnstfield(phiu,coe,1,lsmlim,timestcount,u_lsm)
     !end do
  !write u_lsm
     if(iplane==1)then
        call write3Darray(u_lsm,'u_lsm_xy'//trim(prosnum),'unformatted')
     elseif(iplane==2)then
        call write3Darray(u_lsm,'u_lsm_yz'//trim(prosnum),'unformatted')
     else
        call write3Darray(u_lsm,'u_lsm_xz'//trim(prosnum),'unformatted')
     end if
  
     deallocate(u_lsm)
  !end if

  !if (iplane==1)then
  !   inquire(file='u_ssm_xy'//trim(prosnum),exist=fexist2)
  !elseif(iplane==2)then
  !   inquire(file='u_ssm_yz'//trim(prosnum),exist=fexist2)
  !else
  !   inquire(file='u_ssm_xz'//trim(prosnum),exist=fexist2)
  !end if
  !if(fexist2)then
  !   write(*,*)'file exist no need to reconstruct'
  !else        
     if (iplane==1)then
        allocate(u_ssm(n2do+1,n3,timestcount))
     elseif(iplane==2)then
        allocate(u_ssm(n1,n2do+1,timestcount))
     else
        allocate(u_ssm(n1,n3,timestcount))
     end if
     
     !do t=1,timestcount      
        call recnstfield(phiu,coe,lsmlim+1,timestcount,timestcount,u_ssm)
     !end do
  
     !write u_ssm
     if(iplane==1)then
        call write3Darray(u_ssm,'u_ssm_xy'//trim(prosnum),'unformatted')
     elseif(iplane==2)then
        call write3Darray(u_ssm,'u_ssm_yz'//trim(prosnum),'unformatted')
     else
        call write3Darray(u_ssm,'u_ssm_xz'//trim(prosnum),'unformatted')
     end if
 
     deallocate(u_ssm)
  !end if
     deallocate(phiu)
  
     !if (iplane==1)then
     !   inquire(file='v_lsm_xy'//trim(prosnum),exist=fexist3)
     !elseif(iplane==2)then
     !  inquire(file='v_lsm_yz'//trim(prosnum),exist=fexist3)
     !else
     !  inquire(file='v_lsm_xz'//trim(prosnum),exist=fexist3)
     !end if
     !if(fexist3)then
     !  write(*,*)'file exist no need to reconstruct'
     !else
     ! compute v_lsm
     if (iplane==1)then
        allocate(v_lsm(n2do+1,n3,timestcount))
     elseif(iplane==2)then
        allocate(v_lsm(n1,n2do+1,timestcount))
     else
        allocate(v_lsm(n1,n3,timestcount))
     end if
     
     !do t=1,timestcount      
        call recnstfield(phiv,coe,1,lsmlim,timestcount,v_lsm)
     !end do
     
     !write v_lsm
     if(iplane==1)then
        call write3Darray(v_lsm,'v_lsm_xy'//trim(prosnum),'unformatted')
     elseif(iplane==2)then
        call write3Darray(v_lsm,'v_lsm_yz'//trim(prosnum),'unformatted')
     else
        call write3Darray(v_lsm,'v_lsm_xz'//trim(prosnum),'unformatted')
     end if
  
     deallocate(v_lsm)
     ! end if
     
  !compute v_ssm
     !if (iplane==1)then
     !   inquire(file='v_ssm_xy'//trim(prosnum),exist=fexist4)
     !elseif(iplane==2)then
     !   inquire(file='v_ssm_yz'//trim(prosnum),exist=fexist4)
     !else
     !   inquire(file='v_ssm_xz'//trim(prosnum),exist=fexist4)
     !end if
     !if(fexist4)then
     !   write(*,*)'file exist no need to reconstruct'
     !else
     if (iplane==1)then
        allocate(v_ssm(n2do+1,n3,timestcount))
     elseif(iplane==2)then
        allocate(v_ssm(n1,n2do+1,timestcount))
     else
        allocate(v_ssm(n1,n3,timestcount))
     end if
     
!     do t=1,timestcount      
        call recnstfield(phiv,coe,lsmlim+1,timestcount,timestcount,v_ssm)
!     end do
     
     !write v_ssm
     if(iplane==1)then
        call write3Darray(v_ssm,'v_ssm_xy'//trim(prosnum),'unformatted')
     elseif(iplane==2)then
        call write3Darray(v_ssm,'v_ssm_yz'//trim(prosnum),'unformatted')
     else
        call write3Darray(v_ssm,'v_ssm_xz'//trim(prosnum),'unformatted')
     end if
     deallocate(v_ssm)
     !end if
     
     deallocate(phiv)
     
     if (iplane==1)then
        allocate(w_lsm(n2do+1,n3,timestcount))
     elseif(iplane==2)then
        allocate(w_lsm(n1,n2do+1,timestcount))
     else
        allocate(w_lsm(n1,n3,timestcount))
     end if
     
    ! do t=1,timestcount      
        call recnstfield(phiw,coe,1,lsmlim,timestcount,w_lsm)
     !end do
  
     !write w_lsm
     if(iplane==1)then
        call write3Darray(w_lsm,'w_lsm_xy'//trim(prosnum),'unformatted')
     elseif(iplane==2)then
        call write3Darray(w_lsm,'w_lsm_yz'//trim(prosnum),'unformatted')
     else
        call write3Darray(w_lsm,'w_lsm_xz'//trim(prosnum),'unformatted')
     end if
     
     deallocate(w_lsm)

     if (iplane==1)then
        allocate(w_ssm(n2do+1,n3,timestcount))
     elseif(iplane==2)then
        allocate(w_ssm(n1,n2do+1,timestcount))
     else
        allocate(w_ssm(n1,n3,timestcount))
     end if
     
    ! do t=1,timestcount      
        call recnstfield(phiw,coe,lsmlim+1,timestcount,timestcount,w_ssm)
     !end do
     
     !write w_ssm
     if(iplane==1)then
        call write3Darray(w_ssm,'w_ssm_xy'//trim(prosnum),'unformatted')
     elseif(iplane==2)then
        call write3Darray(w_ssm,'w_ssm_yz'//trim(prosnum),'unformatted')
     else
        call write3Darray(w_ssm,'w_ssm_xz'//trim(prosnum),'unformatted')
     end if
     deallocate(w_ssm)
     deallocate(phiw)
  end if
  
  ! This part of the program compute Reynold stresses
  if(ireyst/=0)then
     if (iplane==1)then
        allocate(v_lsm(n2do+1,n3,timestcount),u_lsm(n2do+1,n3,timestcount),w_lsm(n2do+1,n3,timestcount),& 
             u_ssm(n2do+1,n3,timestcount),v_ssm(n2do+1,n3,timestcount),w_ssm(n2do+1,n3,timestcount))
     elseif(iplane==2)then
        allocate(v_lsm(n1,n2do+1,timestcount),u_lsm(n1,n2do+1,timestcount),w_lsm(n1,n2do+1,timestcount),&
             v_ssm(n1,n2do+1,timestcount),u_ssm(n1,n2do+1,timestcount),w_ssm(n1,n2do+1,timestcount))
     else
        allocate(v_lsm(n1,n3,timestcount),u_lsm(n1,n3,timestcount),w_lsm(n1,n3,timestcount),& 
             v_ssm(n1,n3,timestcount),u_ssm(n1,n3,timestcount),w_ssm(n1,n3,timestcount))
     end if

     ! This section of the program computes Reynolds stresses.
     ! When ireyst/=0 Reynolds stresses are computed.
     ! ireyst=3 for u2
     ! ireyst=2 for v2
     ! ireyst=4 for uv
     
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
        call read3Darray('w_lsm_yz'//trim(prosnum),'unformatted',w_lsm)
        call read3Darray('u_ssm_yz'//trim(prosnum),'unformatted',u_ssm)
        call read3Darray('v_ssm_yz'//trim(prosnum),'unformatted',v_ssm)
        call read3Darray('w_ssm_yz'//trim(prosnum),'unformatted',w_ssm)
        
        allocate(uv1_lsm(n2do+1),uv1_ssm(n2do+1),uv1_tot(n2do+1),uv2_lsm(n1,n2do+1),uv2_ssm(n1,n2do+1),uv2_tot(n1,n2do+1))
        
        if(ireyst==4)then ! computing uv
           
           call reyst_ssm_lsm(u_lsm,v_lsm,ibs,iplane,uv1_lsm,uv2_lsm)
           call reyst_ssm_lsm(u_ssm,v_ssm,ibs,iplane,uv1_ssm,uv2_ssm)
           call reyst_ssm_lsm(u_lsm+u_ssm,v_lsm+v_ssm,ibs,iplane,uv1_tot,uv2_tot)
           ! write uv 
           if(ibs==0) then! no blowing
              call senRev1dwrite(uv1_lsm,1,'uv_lsm_yz.dat')
              call senRev1dwrite(uv1_ssm,2,'uv_ssm_yz.dat')
              call senRev1dwrite(uv1_tot,3,'uv_tot_yz.dat')
           else
              call sendrecv2dwrite(uv2_lsm,1,iplane,'uv_lsm'//trim(xloc)//'_yz.dat')
              call sendrecv2dwrite(uv2_ssm,2,iplane,'uv_ssm'//trim(xloc)//'_yz.dat')
              call sendrecv2dwrite(uv2_tot,3,iplane,'uv_tot'//trim(xloc)//'_yz.dat')
           end if
        elseif(ireyst==2)then ! compute v2
        
           call reyst_ssm_lsm(v_lsm,v_lsm,ibs,iplane,uv1_lsm,uv2_lsm)
           call reyst_ssm_lsm(v_ssm,v_ssm,ibs,iplane,uv1_ssm,uv2_ssm)
           call reyst_ssm_lsm(v_lsm+v_ssm,v_lsm+v_ssm,ibs,iplane,uv1_tot,uv2_tot)
           ! write v2
           if(ibs==0) then! no blowing
              call senRev1dwrite(uv1_lsm,1,'v2_lsm_yz.dat')
              call senRev1dwrite(uv1_ssm,2,'v2_ssm_yz.dat')
              call senRev1dwrite(uv1_tot,3,'v2_tot_yz.dat')
           else
              call sendrecv2dwrite(uv2_lsm,1,iplane,'v2_lsm'//trim(xloc)//'_yz.dat')
              call sendrecv2dwrite(uv2_ssm,2,iplane,'v2_ssm'//trim(xloc)//'_yz.dat')
              call sendrecv2dwrite(uv2_tot,3,iplane,'v2_tot'//trim(xloc)//'_yz.dat')
           end if        

        elseif(ireyst==3)then !compute u2
        
           call reyst_ssm_lsm(u_lsm,u_lsm,ibs,iplane,uv1_lsm,uv2_lsm)
           call reyst_ssm_lsm(u_ssm,u_ssm,ibs,iplane,uv1_ssm,uv2_ssm)
           call reyst_ssm_lsm(u_lsm+u_ssm,u_lsm+u_ssm,ibs,iplane,uv1_tot,uv2_tot)
           ! write u2
           if(ibs==0) then! no blowing
              call senRev1dwrite(uv1_lsm,1,'u2_lsm_yz.dat')
              call senRev1dwrite(uv1_ssm,2,'u2_ssm_yz.dat')
              call senRev1dwrite(uv1_tot,3,'u2_tot_yz.dat')
           else
              call sendrecv2dwrite(uv2_lsm,1,iplane,'u2_lsm'//trim(xloc)//'_yz.dat')
              call sendrecv2dwrite(uv2_ssm,2,iplane,'u2_ssm'//trim(xloc)//'_yz.dat')
              call sendrecv2dwrite(uv2_tot,3,iplane,'u2_tot'//trim(xloc)//'_yz.dat')
           end if
        elseif(ireyst==1)then !compute w2
        
           call reyst_ssm_lsm(w_lsm,w_lsm,ibs,iplane,uv1_lsm,uv2_lsm)
           call reyst_ssm_lsm(w_ssm,w_ssm,ibs,iplane,uv1_ssm,uv2_ssm)
           call reyst_ssm_lsm(w_lsm+w_ssm,w_lsm+w_ssm,ibs,iplane,uv1_tot,uv2_tot)
           ! write u2
           if(ibs==0) then! no blowing
              call senRev1dwrite(uv1_lsm,1,'w2_lsm_yz.dat')
              call senRev1dwrite(uv1_ssm,2,'w2_ssm_yz.dat')
              call senRev1dwrite(uv1_tot,3,'w2_tot_yz.dat')
           else
              call sendrecv2dwrite(uv2_lsm,1,iplane,'w2_lsm'//trim(xloc)//'_yz.dat')
              call sendrecv2dwrite(uv2_ssm,2,iplane,'w2_ssm'//trim(xloc)//'_yz.dat')
              call sendrecv2dwrite(uv2_tot,3,iplane,'w2_tot'//trim(xloc)//'_yz.dat')
           end if
        end if
        
     elseif(iplane==1)then
        call read3Darray('u_lsm_xy'//trim(prosnum),'unformatted',u_lsm)
        call read3Darray('v_lsm_xy'//trim(prosnum),'unformatted',v_lsm)
        call read3Darray('w_lsm_xy'//trim(prosnum),'unformatted',w_lsm)
        call read3Darray('u_ssm_xy'//trim(prosnum),'unformatted',u_ssm)
        call read3Darray('v_ssm_xy'//trim(prosnum),'unformatted',v_ssm)
        call read3Darray('w_ssm_xy'//trim(prosnum),'unformatted',w_ssm)
        
        allocate(uv1_lsm(n2do+1),uv1_ssm(n2do+1),uv1_tot(n2do+1),uv2_lsm(n2do+1,n3),uv2_ssm(n2do+1,n3),uv2_tot(n2do+1,n3))
        
        if(ireyst==4)then ! computing uv
           
           call reyst_ssm_lsm(u_lsm,v_lsm,ibs,iplane,uv1_lsm,uv2_lsm)
           call reyst_ssm_lsm(u_ssm,v_ssm,ibs,iplane,uv1_ssm,uv2_ssm)
           call reyst_ssm_lsm(u_lsm+u_ssm,v_lsm+v_ssm,ibs,iplane,uv1_tot,uv2_tot)
           ! write uv 
           if(ibs==0) then! no blowing
              call senRev1dwrite(uv1_lsm,1,'uv_lsm_xy.dat')
              call senRev1dwrite(uv1_ssm,2,'uv_ssm_xy.dat')
              call senRev1dwrite(uv1_tot,3,'uv_tot_xy.dat')
           else
              call sendrecv2dwrite(uv2_lsm,1,iplane,'uv_lsm'//trim(xloc)//'_xy.dat')
              call sendrecv2dwrite(uv2_ssm,2,iplane,'uv_ssm'//trim(xloc)//'_xy.dat')
              call sendrecv2dwrite(uv2_tot,3,iplane,'uv_tot'//trim(xloc)//'_xy.dat')
           end if
        elseif(ireyst==2)then ! compute v2
        
           call reyst_ssm_lsm(v_lsm,v_lsm,ibs,iplane,uv1_lsm,uv2_lsm)
           call reyst_ssm_lsm(v_ssm,v_ssm,ibs,iplane,uv1_ssm,uv2_ssm)
           call reyst_ssm_lsm(v_lsm+v_ssm,v_lsm+v_ssm,ibs,iplane,uv1_tot,uv2_tot)
           ! write v2
           if(ibs==0) then! no blowing
              call senRev1dwrite(uv1_lsm,1,'v2_lsm_xy.dat')
              call senRev1dwrite(uv1_ssm,2,'v2_ssm_xy.dat')
              call senRev1dwrite(uv1_tot,3,'v2_tot_xy.dat')
           else
              call sendrecv2dwrite(uv2_lsm,1,iplane,'v2_lsm'//trim(xloc)//'_xy.dat')
              call sendrecv2dwrite(uv2_ssm,2,iplane,'v2_ssm'//trim(xloc)//'_xy.dat')
              call sendrecv2dwrite(uv2_tot,3,iplane,'v2_tot'//trim(xloc)//'_xy.dat')
           end if        

        elseif(ireyst==3)then !compute u2
        
           call reyst_ssm_lsm(u_lsm,u_lsm,ibs,iplane,uv1_lsm,uv2_lsm)
           call reyst_ssm_lsm(u_ssm,u_ssm,ibs,iplane,uv1_ssm,uv2_ssm)
           call reyst_ssm_lsm(u_lsm+u_ssm,u_lsm+u_ssm,ibs,iplane,uv1_tot,uv2_tot)
           ! write u2
           if(ibs==0) then! no blowing
              call senRev1dwrite(uv1_lsm,1,'u2_lsm_xy.dat')
              call senRev1dwrite(uv1_ssm,2,'u2_ssm_xy.dat')
              call senRev1dwrite(uv1_tot,3,'u2_tot_xy.dat')
           else
              call sendrecv2dwrite(uv2_lsm,1,iplane,'u2_lsm'//trim(xloc)//'_xy.dat')
              call sendrecv2dwrite(uv2_ssm,2,iplane,'u2_ssm'//trim(xloc)//'_xy.dat')
              call sendrecv2dwrite(uv2_tot,3,iplane,'u2_tot'//trim(xloc)//'_xy.dat')
           end if
        elseif(ireyst==1)then !compute w2
        
           call reyst_ssm_lsm(w_lsm,w_lsm,ibs,iplane,uv1_lsm,uv2_lsm)
           call reyst_ssm_lsm(w_ssm,w_ssm,ibs,iplane,uv1_ssm,uv2_ssm)
           call reyst_ssm_lsm(w_lsm+w_ssm,w_lsm+w_ssm,ibs,iplane,uv1_tot,uv2_tot)
           ! write u2
           if(ibs==0) then! no blowing
              call senRev1dwrite(uv1_lsm,1,'w2_lsm_yz.dat')
              call senRev1dwrite(uv1_ssm,2,'w2_ssm_yz.dat')
              call senRev1dwrite(uv1_tot,3,'w2_tot_yz.dat')
           else
              call sendrecv2dwrite(uv2_lsm,1,iplane,'w2_lsm'//trim(xloc)//'_yz.dat')
              call sendrecv2dwrite(uv2_ssm,2,iplane,'w2_ssm'//trim(xloc)//'_yz.dat')
              call sendrecv2dwrite(uv2_tot,3,iplane,'w2_tot'//trim(xloc)//'_yz.dat')
           end if
        end if
     end if
     if(mynode==0)then
        call printVector(yp(1:),'ycoord.dat')
        call printVector(xp,'xcoord.dat')
        call printVector(zp,'zcoord.dat')
     end if
  end if

  if(iscfx/=0)then
     if (iplane==1)then
        allocate(v_lsm(n2do+1,n3,timestcount),u_lsm(n2do+1,n3,timestcount),w_lsm(n2do+1,n3,timestcount),& 
             u_ssm(n2do+1,n3,timestcount),v_ssm(n2do+1,n3,timestcount),w_ssm(n2do+1,n3,timestcount))
     elseif(iplane==2)then
        allocate(v_lsm(n1,n2do+1,timestcount),u_lsm(n1,n2do+1,timestcount),w_lsm(n1,n2do+1,timestcount),&
             v_ssm(n1,n2do+1,timestcount),u_ssm(n1,n2do+1,timestcount),w_ssm(n1,n2do+1,timestcount))
     else
        allocate(v_lsm(n1,n3,timestcount),u_lsm(n1,n3,timestcount),w_lsm(n1,n3,timestcount),& 
             v_ssm(n1,n3,timestcount),u_ssm(n1,n3,timestcount),w_ssm(n1,n3,timestcount))
     end if

     ! This section of the program computes Reynolds stresses.
     ! When ireyst/=0 Reynolds stresses are computed.
     ! ireyst=1 for u2
     ! ireyst=2 for v2
     ! ireyst=3 for uv
     
     if(x0<10)then
        write(xloc,'(i1.1)')x0
     elseif(x0<100)then
        write(xloc,'(i2.2)')x0
     elseif(x0<1000)then
        write(xloc,'(i3.3)')x0
     else
        write(xloc,'(i4.4)')x0
     end if
     
     if (iplane==1)then!xy plane
        call read3Darray('u_lsm_xy'//trim(prosnum),'unformatted',u_lsm)
        call read3Darray('v_lsm_xy'//trim(prosnum),'unformatted',v_lsm)
        call read3Darray('w_lsm_xy'//trim(prosnum),'unformatted',w_lsm)
        call read3Darray('u_ssm_xy'//trim(prosnum),'unformatted',u_ssm)
        call read3Darray('v_ssm_xy'//trim(prosnum),'unformatted',v_ssm)
        call read3Darray('w_ssm_xy'//trim(prosnum),'unformatted',w_ssm)
     
        allocate(uv1_lsm(n2do+1),uv1_ssm(n2do+1),uv1_tot(n2do+1),uv2_lsm(n1,n2do+1),uv2_ssm(n1,n2do+1),uv2_tot(n1,n2do+1))
        
        if(iscfx==3)then ! computing ut
           
           call sclflx_ssm_lsm(u_lsm,tprit,ibs,iplane,uv1_lsm,uv2_lsm)
           call sclflx_ssm_lsm(u_ssm,tprit,ibs,iplane,uv1_ssm,uv2_ssm)
           call sclflx_ssm_lsm(u_lsm+u_ssm,tprit,ibs,iplane,uv1_tot,uv2_tot)
           ! write ut
           if(ibs==0) then! no blowing
              call senRev1dwrite(uv1_lsm,1,'ut_lsm_xy.dat')
              call senRev1dwrite(uv1_ssm,2,'ut_ssm_xy.dat')
              call senRev1dwrite(uv1_tot,3,'ut_tot_xy.dat')
           else
              call sendrecv2dwrite(uv2_lsm,1,iplane,'ut_lsm'//trim(xloc)//'_xy.dat')
              call sendrecv2dwrite(uv2_ssm,2,iplane,'ut_ssm'//trim(xloc)//'_xy.dat')
              call sendrecv2dwrite(uv2_tot,3,iplane,'ut_tot'//trim(xloc)//'_xy.dat')
           end if
           
        elseif(iscfx==2)then ! compute vt        
           
           call sclflx_ssm_lsm(v_lsm,tprit,ibs,iplane,uv1_lsm,uv2_lsm)
           call sclflx_ssm_lsm(v_ssm,tprit,ibs,iplane,uv1_ssm,uv2_ssm)
           call sclflx_ssm_lsm(v_lsm+v_ssm,tprit,ibs,iplane,uv1_tot,uv2_tot)
           ! write vt
           if(ibs==0) then! no blowing
              call senRev1dwrite(uv1_lsm,1,'vt_lsm_xy.dat')
              call senRev1dwrite(uv1_ssm,2,'vt_ssm_xy.dat')
              call senRev1dwrite(uv1_tot,3,'vt_tot_xy.dat')
           else
              call sendrecv2dwrite(uv2_lsm,1,iplane,'vt_lsm'//trim(xloc)//'_xy.dat')
              call sendrecv2dwrite(uv2_ssm,2,iplane,'vt_ssm'//trim(xloc)//'_xy.dat')
              call sendrecv2dwrite(uv2_tot,3,iplane,'vt_tot'//trim(xloc)//'_xy.dat')
           end if
           
        elseif(iscfx==1)then !compute wt
           
           call sclflx_ssm_lsm(w_lsm,tprit,ibs,iplane,uv1_lsm,uv2_lsm)
           call sclflx_ssm_lsm(w_ssm,tprit,ibs,iplane,uv1_ssm,uv2_ssm)
           call sclflx_ssm_lsm(w_lsm+w_ssm,tprit,ibs,iplane,uv1_tot,uv2_tot)
           ! write wt
           if(ibs==0) then! no blowing
              call senRev1dwrite(uv1_lsm,1,'wt_lsm_xy.dat')
              call senRev1dwrite(uv1_ssm,2,'wt_ssm_xy.dat')
              call senRev1dwrite(uv1_tot,3,'wt_tot_xy.dat')
           else
              call sendrecv2dwrite(uv2_lsm,1,iplane,'wt_lsm'//trim(xloc)//'_xy.dat')
              call sendrecv2dwrite(uv2_ssm,2,iplane,'wt_ssm'//trim(xloc)//'_xy.dat')
              call sendrecv2dwrite(uv2_tot,3,iplane,'wt_tot'//trim(xloc)//'_xy.dat')
           end if
           
        end if
           
     elseif (iplane==2)then !yz plane
        
        call read3Darray('u_lsm_yz'//trim(prosnum),'unformatted',u_lsm)
        call read3Darray('v_lsm_yz'//trim(prosnum),'unformatted',v_lsm)
        call read3Darray('w_lsm_yz'//trim(prosnum),'unformatted',w_lsm)
        call read3Darray('u_ssm_yz'//trim(prosnum),'unformatted',u_ssm)
        call read3Darray('v_ssm_yz'//trim(prosnum),'unformatted',v_ssm)
        call read3Darray('w_ssm_yz'//trim(prosnum),'unformatted',w_ssm)
     
        allocate(uv1_lsm(n2do+1),uv1_ssm(n2do+1),uv1_tot(n2do+1),uv2_lsm(n2do+1,n3),uv2_ssm(n2do+1,n3),uv2_tot(n2do+1,n3))
     
        if(iscfx==3)then ! computing ut
           
           call sclflx_ssm_lsm(u_lsm,tprit,ibs,iplane,uv1_lsm,uv2_lsm)
           call sclflx_ssm_lsm(u_ssm,tprit,ibs,iplane,uv1_ssm,uv2_ssm)
           call sclflx_ssm_lsm(u_lsm+u_ssm,tprit,ibs,iplane,uv1_tot,uv2_tot)
           ! write ut
           if(ibs==0) then! no blowing
              call senRev1dwrite(uv1_lsm,1,'ut_lsm_yz.dat')
              call senRev1dwrite(uv1_ssm,2,'ut_ssm_yz.dat')
              call senRev1dwrite(uv1_tot,3,'ut_tot_yz.dat')
           else
              call sendrecv2dwrite(uv2_lsm,1,iplane,'ut_lsm'//trim(xloc)//'_yz.dat')
              call sendrecv2dwrite(uv2_ssm,2,iplane,'ut_ssm'//trim(xloc)//'_yz.dat')
              call sendrecv2dwrite(uv2_tot,3,iplane,'ut_tot'//trim(xloc)//'_yz.dat')
           end if
           
        elseif(iscfx==2)then ! compute vt        
           
           call sclflx_ssm_lsm(v_lsm,tprit,ibs,iplane,uv1_lsm,uv2_lsm)
           call sclflx_ssm_lsm(v_ssm,tprit,ibs,iplane,uv1_ssm,uv2_ssm)
           call sclflx_ssm_lsm(v_lsm+v_ssm,tprit,ibs,iplane,uv1_tot,uv2_tot)
           ! write vt
           if(ibs==0) then! no blowing
              call senRev1dwrite(uv1_lsm,1,'vt_lsm_yz.dat')
              call senRev1dwrite(uv1_ssm,2,'vt_ssm_yz.dat')
              call senRev1dwrite(uv1_tot,3,'vt_tot_yz.dat')
           else
              call sendrecv2dwrite(uv2_lsm,1,iplane,'vt_lsm'//trim(xloc)//'_yz.dat')
              call sendrecv2dwrite(uv2_ssm,2,iplane,'vt_ssm'//trim(xloc)//'_yz.dat')
              call sendrecv2dwrite(uv2_tot,3,iplane,'vt_tot'//trim(xloc)//'_yz.dat')
           end if
           
        elseif(iscfx==1)then !compute wt
           
           call sclflx_ssm_lsm(w_lsm,tprit,ibs,iplane,uv1_lsm,uv2_lsm)
           call sclflx_ssm_lsm(w_ssm,tprit,ibs,iplane,uv1_ssm,uv2_ssm)
           call sclflx_ssm_lsm(w_lsm+w_ssm,tprit,ibs,iplane,uv1_tot,uv2_tot)
           ! write wt
           if(ibs==0) then! no blowing
              call senRev1dwrite(uv1_lsm,1,'wt_lsm_yz.dat')
              call senRev1dwrite(uv1_ssm,2,'wt_ssm_yz.dat')
              call senRev1dwrite(uv1_tot,3,'wt_tot_yz.dat')
           else
              call sendrecv2dwrite(uv2_lsm,1,iplane,'wt_lsm'//trim(xloc)//'_yz.dat')
              call sendrecv2dwrite(uv2_ssm,2,iplane,'wt_ssm'//trim(xloc)//'_yz.dat')
              call sendrecv2dwrite(uv2_tot,3,iplane,'wt_tot'//trim(xloc)//'_yz.dat')
           end if
           
        end if
     end if
     if(mynode==0)then
        call printVector(yp(1:),'ycoord.dat')
        call printVector(xp,'xcoord.dat')
        call printVector(zp,'zcoord.dat')
     end if
  end if
  deallocate(tprit)


  ! This part of the program writes the instantaneous field for different planes
  if(inst==1)then  
     if (iplane==1)then
        allocate(u_lsm(n2do+1,n3,numtimesteps+1),u_ssm(n2do+1,n3,numtimesteps+1))
        call read3Darray('w_lsm_xy'//trim(prosnum),'unformatted',u_lsm)
        call read3Darray('w_ssm_xy'//trim(prosnum),'unformatted',u_ssm)
     elseif(iplane==2)then
        allocate(u_lsm(n1,n2do+1,numtimesteps+1),u_ssm(n1,n2do+1,numtimesteps+1))
        call read3Darray('w_lsm_yz'//trim(prosnum),'unformatted',u_lsm)
        call read3Darray('w_ssm_yz'//trim(prosnum),'unformatted',u_ssm)
     else
        allocate(u_lsm(n1,n3,numtimesteps+1),u_ssm(n1,n3,numtimesteps+1))
        call read3Darray('w_lsm_xz'//trim(prosnum),'unformatted',u_lsm)
        call read3Darray('w_ssm_xz'//trim(prosnum),'unformatted',u_ssm)
     end if
        
     do tt=1,numtimesteps+1,500
           
        if(tt<10)then
           write(time,'(i1.1)')tt
        elseif(tt<100)then
           write(time,'(i2.2)')tt
        elseif(tt<1000)then
           write(time,'(i3.3)')tt
        else
           write(time,'(i4.4)')tt
        end if

        if (iplane==1)then
           fname_l='w_lsm_inst_xy_'//trim(time)//'.dat'
           fname_s='w_ssm_inst_xy_'//trim(time)//'.dat'
        elseif(iplane==2)then
           fname_l='w_lsm_inst_yz_'//trim(time)//'.dat'
           fname_s='w_ssm_inst_yz_'//trim(time)//'.dat'
        else
           fname_l='w_lsm_inst_xz_'//trim(time)//'.dat'
           fname_s='w_ssm_inst_xz_'//trim(time)//'.dat'
        end if
        
        call sendrecv2dwrite(u_lsm(:,:,tt),1,iplane,fname_l)
        call sendrecv2dwrite(u_ssm(:,:,tt),2,iplane,fname_s)
        
     end do
        
  elseif(inst==2)then
     
     if (iplane==1)then
        allocate(u_lsm(n2do+1,n3,numtimesteps+1),u_ssm(n2do+1,n3,numtimesteps+1))
        call read3Darray('v_lsm_xy'//trim(prosnum),'unformatted',u_lsm)
        call read3Darray('v_ssm_xy'//trim(prosnum),'unformatted',u_ssm)

     elseif(iplane==2)then
        allocate(u_lsm(n1,n2do+1,numtimesteps+1),u_ssm(n1,n2do+1,numtimesteps+1))
        call read3Darray('v_lsm_yz'//trim(prosnum),'unformatted',u_lsm)
        call read3Darray('v_ssm_yz'//trim(prosnum),'unformatted',u_ssm)
        
     else
        allocate(u_lsm(n1,n3,numtimesteps+1),u_ssm(n1,n3,numtimesteps+1))
        call read3Darray('v_lsm_xz'//trim(prosnum),'unformatted',u_lsm)
        call read3Darray('v_ssm_xz'//trim(prosnum),'unformatted',u_ssm)
     end if
        
     do tt=1,numtimesteps+1,500
           
        if(tt<10)then
           write(time,'(i1.1)')tt
        elseif(tt<100)then
           write(time,'(i2.2)')tt
        elseif(tt<1000)then
           write(time,'(i3.3)')tt
        else
           write(time,'(i4.4)')tt
        end if
        
        if (iplane==1)then
           fname_l='v_lsm_inst_xy_'//trim(time)//'.dat'
           fname_s='v_ssm_inst_xy_'//trim(time)//'.dat'
        elseif(iplane==2)then
           fname_l='v_lsm_inst_yz_'//trim(time)//'.dat'
           fname_s='v_ssm_inst_yz_'//trim(time)//'.dat'
        else
           fname_l='v_lsm_inst_xz_'//trim(time)//'.dat'
           fname_s='v_ssm_inst_xz_'//trim(time)//'.dat'
        end if
        
        call sendrecv2dwrite(u_lsm(:,:,tt),1,iplane,fname_l)
        call sendrecv2dwrite(u_ssm(:,:,tt),2,iplane,fname_s)
        
     end do

  elseif(inst==3)then

     if (iplane==1)then
        allocate(u_lsm(n2do+1,n3,numtimesteps+1),u_ssm(n2do+1,n3,numtimesteps+1))
        call read3Darray('u_lsm_xy'//trim(prosnum),'unformatted',u_lsm)
        call read3Darray('u_ssm_xy'//trim(prosnum),'unformatted',u_ssm)
     elseif(iplane==2)then
        allocate(u_lsm(n1,n2do+1,numtimesteps+1),u_ssm(n1,n2do+1,numtimesteps+1))
        call read3Darray('u_lsm_yz'//trim(prosnum),'unformatted',u_lsm)
        call read3Darray('u_ssm_yz'//trim(prosnum),'unformatted',u_ssm)
     else
        allocate(u_lsm(n1,n3,numtimesteps+1),u_ssm(n1,n3,numtimesteps+1))
        call read3Darray('u_lsm_xz'//trim(prosnum),'unformatted',u_lsm)
        call read3Darray('u_ssm_xz'//trim(prosnum),'unformatted',u_ssm)
     end if
        
     do tt=1,numtimesteps+1,500
           
        if(tt<10)then
           write(time,'(i1.1)')tt
        elseif(tt<100)then
           write(time,'(i2.2)')tt
        elseif(tt<1000)then
           write(time,'(i3.3)')tt
        else
           write(time,'(i4.4)')tt
        end if
        
        if (iplane==1)then
           fname_l='u_lsm_inst_xy_'//trim(time)//'.dat'
           fname_s='u_ssm_inst_xy_'//trim(time)//'.dat'
        elseif(iplane==2)then
           fname_l='u_lsm_inst_yz_'//trim(time)//'.dat'
           fname_s='u_ssm_inst_yz_'//trim(time)//'.dat'
        else
           fname_l='u_lsm_inst_xz_'//trim(time)//'.dat'
           fname_s='u_ssm_inst_xz_'//trim(time)//'.dat'
        end if
        
        call sendrecv2dwrite(u_lsm(:,:,tt),1,iplane,fname_l)
        call sendrecv2dwrite(u_ssm(:,:,tt),2,iplane,fname_s)
        
     end do
        
  end if

  if(icorl/=0)then
     call tpcorl_2d_podfield(var1,var2,numtimesteps,x0,j0,k0,iplane)
  end if

  if(ireyst==0.and.iscfx==0)then
     if(mynode==0)then
        call printVector(yp(1:),'ycoord.dat')
        call printVector(xp,'xcoord.dat')
        call printVector(zp,'zcoord.dat')
     end if
  end if
  call mpi_finalize(ierr)
end program programpod
        
        
  
