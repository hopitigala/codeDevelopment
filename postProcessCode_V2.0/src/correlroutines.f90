module compcorrl
  implicit none
contains
  subroutine valatcnstz0y0x0(ary3d,i0,j0,k0,node,refval)
    use mpi
    use mpivariables
    use channhdvariables
    use prelimcalvar
    implicit none
    real*8,intent(in),dimension(:,:,:)::ary3d
    integer,intent(in)                :: i0,j0,k0
    real*8, intent(out)               :: refval
    integer,intent(out)               :: node
    integer                           :: yprime
    
    node=j0/n2do
    if(mynode==node)then
       yprime=j0-node*n2do
       refval=ary3d(i0,yprime,k0)
    end if
    ! broadcast value to all other processes
    !call mpi_bcast(refval,1,mpi_real8,node,mpi_comm_world,ierr)
  end subroutine valatcnstz0y0x0

  subroutine corlcoe_inhomo(covar,rmsx0y0z0,rmsary,corlcoe)
    implicit none
    real*8,intent(in),dimension(:,:,:) ::covar,rmsary
    real*8,intent(in)                  ::rmsx0y0z0
    real*8,intent(out),dimension(:,:,:)::corlcoe
    real*8,allocatable,dimension(:,:,:)::inv
    integer                            ::d1,d2,d3
    d1=size(rmsary,1)
    d2=size(rmsary,2)
    d3=size(rmsary,3)
    allocate(inv(d1,d2,d3))
    
    inv=1.0/(rmsx0y0z0*rmsary)
    corlcoe=covar*inv
    deallocate(inv)
  end subroutine corlcoe_inhomo

  subroutine tpcorl_inhomo(var1,var2,i0,j0,k0,numtstps,dt,sttime,ysdo,ypdo,corlcoe)
    
    use channhdvariables
    use mpivariables
    use readdata
    use prelimcalvar
    use interpoldata
    use fluctuation
    use mpi
    use arrayops
    
    implicit none
    
    integer,intent(in)                   ::var1,var2,j0,k0,numtstps,dt,i0
    real*8,intent(in)                    ::sttime
    real*8,intent(in),dimension(0:)      ::ysdo,ypdo
    real*8,intent(out),dimension(:,:,:)  ::corlcoe

    real*8,allocatable,dimension(:,:,:)  ::umean,vmean,covar,covartavg
    real*8,allocatable,dimension(:,:,:)  ::urms,vrms,pp,uprime,vprime
    real*8,allocatable,dimension(:,:,:,:)::up,tp
    real*8                               ::inv,refval,rmsrefval
    character*2                          ::prosnum
    integer                              ::d1,d2,d3,timestcount,itime,node1,node2

    d1=size(corlcoe,1)
    d2=size(corlcoe,2)
    d3=size(corlcoe,3)
      
  ! read 3D mean flow data from file
    write(prosnum,'(i2.2)')mynode
    
    allocate(umean(d1,d2,d3),vmean(d1,d2,d3))
    if (var1==1)then
       call read3Darray('../../tmean/1500ts/wmean'//trim(prosnum),'unformatted',umean)
    elseif (var1==2)then
       call read3Darray('../../tmean/1500ts/vmean'//trim(prosnum),'unformatted',umean)
    elseif (var1==3) then
       call read3Darray('../../tmean/1500ts/umean'//trim(prosnum),'unformatted',umean)
    elseif (var1==4)then
       call read3Darray('../../tmean/1500ts/tmean'//trim(prosnum),'unformatted',umean)!umean==tmean
    end if

    if (var2==1)then
       call read3Darray('../../tmean/1500ts/wmean'//trim(prosnum),'unformatted',vmean)
    elseif (var2==2)then
       call read3Darray('../../tmean/1500ts/vmean'//trim(prosnum),'unformatted',vmean)
    elseif (var2==3) then
       call read3Darray('../../tmean/1500ts/umean'//trim(prosnum),'unformatted',vmean)
    else if (var2==4)then
       call read3darray('../../tmean/1500ts/tmean'//trim(prosnum),'unformatted',vmean)!vmean==tmean
    end if
    
    allocate(covartavg(d1,d2,d3))
    
    covartavg=0.
    ! time loop begins
    timestcount=0
    do itime=int(sttime),(int(sttime)+(numtstps*dt)),dt
       timestcount=timestcount+1
       allocate(up(d1,d2,d3,3),tp(d1,d2,d3,1),pp(d1,d2,d3))
       call readTempFieldData(mynode,itime,ypdo,ysdo,up,pp,tp)
       deallocate(pp)
       
       ! find fluctuating field
     
       allocate(uprime(d1,d2,d3),vprime(d1,d2,d3))
       if(var1==4)then !temperatuer fluctuations
          call fluct3Dmean(tp(:,:,:,1),umean,uprime)
       else
          call fluct3Dmean(up(:,:,:,var1),umean,uprime)
       end if
       if (var2==4)then! temperatuer fluctuations
          call fluct3Dmean(tp(:,:,:,1),vmean,vprime)
       else
          call fluct3Dmean(up(:,:,:,var2),vmean,vprime)
       end if
       
       call valatcnstz0y0x0(vprime,i0,j0,k0,node1,refval)
       deallocate(up,tp)
       call mpi_bcast(refval,1,mpi_real8,node1,mpi_comm_world,ierr)
       
       allocate(covar(d1,d2,d3))
       covar=uprime*refval
       deallocate(uprime,vprime)
       call loopAdd3DAry(covar,covartavg)
       deallocate(covar)
    end do
    
    deallocate(umean,vmean)
    
    inv=1.0/real(timestcount)

    covartavg=covartavg*inv

  
    
    allocate(vrms(d1,d2,d3),urms(d1,d2,d3))
    if (var1==1)then
       call read3Darray('../../rms/1500ts/wrms'//trim(prosnum),'unformatted',urms)
    elseif (var1==2)then
       call read3Darray('../../rms/1500ts/vrms'//trim(prosnum),'unformatted',urms)
    elseif (var1==3) then
       call read3Darray('../../rms/1500ts/urms'//trim(prosnum),'unformatted',urms)
    elseif(var1==4)then
       call read3Darray('../../rms/1500ts/trms'//trim(prosnum),'unformatted',urms)
    end if

    if (var2==1)then
       call read3Darray('../../rms/1500ts/wrms'//trim(prosnum),'unformatted',vrms)
    elseif (var2==2)then
       call read3Darray('../../rms/1500ts/vrms'//trim(prosnum),'unformatted',vrms)
    elseif (var2==3) then
       call read3Darray('../../rms/1500ts/urms'//trim(prosnum),'unformatted',vrms)
    elseif(var2==4)then
       call read3Darray('../../rms/1500ts/trms'//trim(prosnum),'unformatted',vrms)
    end if
    
    call valatcnstz0y0x0(vrms,i0,j0,k0,node2,rmsrefval)
  
    call mpi_bcast(rmsrefval,1,mpi_real8,node2,mpi_comm_world,ierr)

    call corlcoe_inhomo(covartavg,rmsrefval,urms,corlcoe)

    deallocate(urms,vrms,covartavg)
  end subroutine tpcorl_inhomo
  
  subroutine corlcoe_inhomo_avg(z0,corlcoe,corlcoe_avg)
    implicit none
    integer,intent(in),dimension(:)::z0
    real*8,intent(in),dimension(:,:,:)::corlcoe
    real*8,intent(out),dimension(:,:,:)::corlcoe_avg
    integer::i0,d1
    d1=size(z0)
    corlcoe_avg=0.0
    do i0=1,d1
       corlcoe_avg=corlcoe_avg+corlcoe
    end do
    corlcoe_avg=(1.0/real(d1))*corlcoe_avg

  end subroutine corlcoe_inhomo_avg

  subroutine corlfilename(var1,var2,i0,j0,k0,icorlavg,icen,filename)
    implicit none
    
    integer,intent(in)::var1,var2,j0,k0,icorlavg,icen
    integer,intent(in),dimension(:)::i0
    character(len=*),intent(out)::filename
    
    character(len=1)::vari1,vari2
    character(len=4)::yprm,xprm,zprm
    
    write(vari1,'(i1.1)')var1
    write(vari2,'(i1.1)')var2
    
    
    if(j0<10)then
       write(yprm,'(i1.1)')j0
    elseif(j0<100)then
       write(yprm,'(i2.2)')j0
    else
       write(yprm,'(i3.3)')j0
    end if
    
    if(k0<10)then
       write(xprm,'(i1.1)')k0
    elseif(k0<100)then
       write(xprm,'(i2.2)')k0
    elseif(k0<1000)then
       write(xprm,'(i3.3)')k0
    else
       write(xprm,'(i4.4)')k0
    end if
    if (icorlavg==1)then
       if (icen==1)then
          filename='corlcoe_'//trim(vari1)//'_'//trim(vari2)//'_cen_avg_y0_'//trim(yprm)//'_x0_'//trim(xprm)//'_.dat'
       else
          filename='corlcoe_'//trim(vari1)//'_'//trim(vari2)//'_mid_avg_y0_'//trim(yprm)//'_x0_'//trim(xprm)//'_.dat'
       end if
    else
       if(i0(1)<10)then
          write(zprm,'(i1.1)')i0(1)
       elseif(i0(1)<100)then
          write(zprm,'(i2.2)')i0(1)
       else
          write(zprm,'(i3.3)')i0(1)
       end if

       filename='corlcoe_'//trim(vari1)//'_'//trim(vari2)//'_z0_'//trim(zprm)//'_y0_'//trim(yprm)//'_x0_'//trim(xprm)//'_.dat'
    end if
  end subroutine corlfilename

  subroutine tpcorl_2dhomo(var1,var2,j0,numtstps,dt,sttime,ysdo,ypdo,corlcoe)
    
    use channhdvariables
    use mpivariables
    use readdata
    use prelimcalvar
    use interpoldata
    use fluctuation
    use mpi
    use arrayops
    use compfft 
    
    implicit none
    
    integer,intent(in)                   ::var1,var2,j0,numtstps,dt
    real*8,intent(in)                    ::sttime
    real*8,intent(in),dimension(0:)      ::ysdo,ypdo
    real*8,intent(out),dimension(:,:,:)  ::corlcoe
    
    complex(kind=8),allocatable,dimension(:,:,:)::uprihat,vprihat
    complex(kind=8),allocatable,dimension(:,:)::refary
    real*8,allocatable,dimension(:,:,:)  ::umean,vmean,covar,covartavg
    real*8,allocatable,dimension(:,:,:)  ::urms,vrms,pp,uprime,vprime
    real*8,allocatable,dimension(:,:,:,:)::up,tp
    real*8                               ::inv,refval,rmsrefval
    character*2                          ::prosnum
    integer                              ::d1,d2,d3,timestcount,itime,node1,node2,l1,l2

    d1=size(corlcoe,1)
    d2=size(corlcoe,2)
    d3=size(corlcoe,3)
      
  ! read 3D mean flow data from file
    write(prosnum,'(i2.2)')mynode
    
    allocate(umean(d1,d2,d3),vmean(d1,d2,d3))
    if (var1==1)then
       call read3Darray('../../tmean/1000tstps/wmean'//trim(prosnum),'unformatted',umean)
    elseif (var1==2)then
       call read3Darray('../../tmean/1000tstps/vmean'//trim(prosnum),'unformatted',umean)
    elseif (var1==3) then
       call read3Darray('../../tmean/1000tstps/umean'//trim(prosnum),'unformatted',umean)
    end if

    if (var2==1)then
       call read3Darray('../../tmean/1000tstps/wmean'//trim(prosnum),'unformatted',vmean)
    elseif (var2==2)then
       call read3Darray('../../tmean/1000tstps/vmean'//trim(prosnum),'unformatted',vmean)
    elseif (var2==3) then
       call read3Darray('../../tmean/1000tstps/umean'//trim(prosnum),'unformatted',vmean)
    end if
    
    allocate(covartavg(2*d1,d2,2*d3))
    
    covartavg=0.
    ! time loop begins
    timestcount=0
    do itime=int(sttime),(int(sttime)+(numtstps*dt)),dt
       timestcount=timestcount+1
       allocate(up(d1,d2,d3,3),tp(d1,d2,d3,1),pp(d1,d2,d3))
       call readTempFieldData(mynode,itime,ypdo,ysdo,up,pp,tp)
       deallocate(tp,pp)
       
       ! find fluctuating field
     
       allocate(uprime(d1,d2,d3),vprime(d1,d2,d3))
       
       call fluct3Dmean(up(:,:,:,var1),umean,uprime)
       call fluct3Dmean(up(:,:,:,var2),vmean,vprime)
       deallocate(up)
       
       allocate(uprihat(2*d1,d2,2*d3),vprihat(2*d1,d2,2*d3))
       l1=size(uprihat,1)
       l2=size(uprihat,3)
       call twodfft(uprime,l1,l2,uprihat)
       call twodfft(vprime,l1,l2,vprihat)
       deallocate(uprime,vprime) 
       allocate(refary(l1,l2))
       
       call ary2dcnsty0(vprihat,j0,node1,refary)
       
       call mpi_bcast(refary,l1*l2,mpi_complex16,node1,mpi_comm_world,ierr)
      
       allocate(covar(d1,d2,d3))
       covar=uprihat()*conjg(refval())
      
       call loopAdd3DAry(covar,covartavg)
       deallocate(covar)
    end do
    
    deallocate(umean,vmean)
    
    inv=1.0/real(timestcount)

    covartavg=covartavg*inv

  
    
    allocate(vrms(d1,d2,d3),urms(d1,d2,d3))
    if (var1==1)then
       call read3Darray('../../rms/1000tstps/wrms'//trim(prosnum),'unformatted',urms)
    elseif (var1==2)then
       call read3Darray('../../rms/1000tstps/vrms'//trim(prosnum),'unformatted',urms)
    elseif (var1==3) then
       call read3Darray('../../rms/1000tstps/urms'//trim(prosnum),'unformatted',urms)
    end if

    if (var2==1)then
       call read3Darray('../../rms/1000tstps/wrms'//trim(prosnum),'unformatted',vrms)
    elseif (var2==2)then
       call read3Darray('../../rms/1000tstps/vrms'//trim(prosnum),'unformatted',vrms)
    elseif (var2==3) then
       call read3Darray('../../rms/1000tstps/urms'//trim(prosnum),'unformatted',vrms)
    end if
    
    !call valatcnstz0y0x0(vrms,i0,j0,k0,node2,rmsrefval)
  
    call mpi_bcast(rmsrefval,1,mpi_real8,node2,mpi_comm_world,ierr)

    call corlcoe_inhomo(covartavg,rmsrefval,urms,corlcoe)

    deallocate(urms,vrms,covartavg)
  end subroutine tpcorl_2dhomo

  subroutine ary2dcnsty0(ary3d,j0,node,refary)
    use mpi
    use mpivariables
    use channhdvariables
    use prelimcalvar
    implicit none
    complex(kind=8),intent(in),dimension(:,:,:)::ary3d
    integer,intent(in)                ::j0
    complex(kind=8),intent(out),dimension(:,:) ::refary
    integer,intent(out)               ::node
    integer                           ::yprime
    
    node=j0/n2do
    if(mynode==node)then
       yprime=j0-node*n2do
       refary(:,:)=ary3d(:,yprime,:)
    end if
    ! broadcast value to all other processes
    !call mpi_bcast(refval,1,mpi_real8,node,mpi_comm_world,ierr)
  end subroutine ary2dcnsty0

end module compcorrl
