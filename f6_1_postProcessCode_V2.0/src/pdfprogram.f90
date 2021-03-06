PROGRAM lambda2program
  USE channhdvariables
  USE mainparameters
  USE writedata
  USE arrayops
  USE readdata
  USE mpivariables
  USE prelimcalvar
  USE interpoldata
  USE mpi
  use vorticity
  use vortexiden
  use fluctuation
  IMPLICIT NONE
  REAL*8,ALLOCATABLE,DIMENSION(:)::yp,xp,zp,ypdo,ysdo,ys
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:)::utavg,vtavg,trms,ttavg,pp,lambda2,tprime,tprimeavg
  INTEGER,ALLOCATABLE,DIMENSION(:,:,:)::tplowcount,tphighcount,lambdanegcount
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:)::pdfhigh,pdflow,probcount
  REAL*8,ALLOCATABLE,DIMENSION(:,:)::utxavg,vtxavg,wtxavg,ttxavg
  REAL*8,ALLOCATABLE,DIMENSION(:)::utxzavg,vtxzavg,wtxzavg,ttxzavg
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:)::up,tp,dudx,dudy,dudz
  INTEGER::itime,dt,timestcount,numtimesteps,i,j,k,icorlavg,icen
  REAL*8::inv,sttime,utau,ttau,retau
  character(len=5)::printim
  character(len=2)::prosnum

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,mynode,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

  CALL readChannelData()
  CALL readPostData(sttime,numtimesteps,dt,icorlavg,icen)
  CALL prelimCal()

  ALLOCATE(xp(n3),yp(0:n2),zp(n1),ys(0:n2),ypdo(0:n2do+1),ysdo(0:n2do+1))
  CALL readCoordData(xp,yp,zp)
  CALL intpolCordi(mynode,numprocs,yp,ys,ypdo,ysdo)
  
  allocate(probcount(n1,n2do+1,n3))
  probcount=0.0

  timestcount=0
  stime=MPI_WTIME()
  do ITime= INT(stTime),(INT(StTime)+(NumTimeSteps*Dt)),Dt
     TimeStCount=TimeStcount+1

     ALLOCATE(up(n1,n2do+1,n3,3),tp(n1,n2do+1,n3,1),pp(n1,n2do+1,n3))
     
     CALL readTempFieldData(mynode,itime,ypdo,ysdo,up,pp,tp)
     
     allocate(dudx(n1,n2do+1,n3,3),dudy(n1,n2do+1,n3,3),dudz(n1,n2do+1,n3,3))
     
     call velograd(up,xp,ypdo,zp,dudx,dudy,dudz)
     deallocate(up,tp,pp)
     allocate(lambda2(n1,n2do+1,n3))
     call complambda2(dudx,dudy,dudz,lambda2)
     
     deallocate(dudx,dudy,dudz)
     do k=1,n3
        do j=1,n2do+1
           do i=1,n1
              if(lambda2(i,j,k)>=-4.0.and.lambda2(i,j,k)<0.0)then
                 probcount(i,j,k)=probcount(i,j,k)+1.0
              end if
           end do
        end do
     end do

     !call sendrecv3dwrite(lambda2,1,'lambda2_'//trim(printim)//'.dat')
     deallocate(lambda2)
     if (mynode==0)then
        write(*,*)'timestep',timestcount,'is done'
     end if
  end do
  probcount=probcount/real(timestcount)
  call sendrecv3dwrite(probcount,1,'probability_4.dat')
  deallocate(probcount)
  if (mynode==0)then
     call printVector(xp,'xcoord.dat')
     call printVector(yp(1:),'ycoord.dat')
     call printVector(zp,'zcoord.dat')
  end if
  call MPI_FINALIZE(ierr)
end PROGRAM lambda2program
