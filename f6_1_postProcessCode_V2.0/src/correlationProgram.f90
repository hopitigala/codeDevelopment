program corelcoe
  use mpi
  use mpivariables
  use compcorrl
  use channhdvariables
  use readdata
  use interpoldata
  use prelimcalvar
  use writedata
  use arrayops
  implicit none
  real*8,allocatable,dimension(:,:,:)::corlcoe,corl_tmp
  real*8,allocatable,dimension(:)    ::yp,xp,zp,ypdo,ysdo,ys
  integer,allocatable,dimension(:)    ::z0
  real*8:: sttime
  integer::numtimesteps,dt,var1,var2,i0,j0,k0,x0,iplane,refpoint,lsmlim
  integer::icorlavg,icen,i
  character(len=100)::filename
  
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,mynode,ierr)
  call mpi_comm_size(mpi_comm_world,numprocs,ierr)
  
  call readChannelData()
  call readPostData(sttime,numtimesteps,dt,icorlavg,icen)
  call prelimcal()
  allocate(xp(n3),yp(0:n2),zp(n1),ys(0:n2),ypdo(0:n2do+1),ysdo(0:n2do+1))
  
  call readCoordData(xp,yp,zp)
  call intpolCordi(mynode,numprocs,yp,ys,ypdo,ysdo)
  allocate(corlcoe(n1,n2do+1,n3))
  if (icorlavg==1)then
     if(icen==1)then
        allocate(z0(5))
     else
        allocate(z0(4))
     end if
     call readcorldata(var1,var2,z0,j0,k0)
     corlcoe=0.0
     do i=1,size(z0)
        allocate(corl_tmp(n1,n2do+1,n3))
        call tpcorl_inhomo(var1,var2,z0(i),j0,k0,numtimesteps,dt,sttime,ysdo,ypdo,corl_tmp)
        call loopadd3DAry(corl_tmp,corlcoe)
        deallocate(corl_tmp)
     end do
     corlcoe=(1.0/real(size(z0)))*corlcoe
  else
     allocate(z0(1))
     call readcorldata(var1,var2,z0,j0,k0)
     call tpcorl_inhomo(var1,var2,z0(1),j0,k0,numtimesteps,dt,sttime,ysdo,ypdo,corlcoe)
  end if
  call corlfilename(var1,var2,z0,j0,k0,icorlavg,icen,filename)
  
  call sendrecv3dwrite(corlcoe,1,filename)
  if(mynode==0)then
     call printVector(yp(1:),'ycoord.dat')
     call printVector(xp,'xcoord.dat')
     call printVector(zp,'zcoord.dat')
  end if
  call mpi_finalize(ierr)
end program corelcoe
