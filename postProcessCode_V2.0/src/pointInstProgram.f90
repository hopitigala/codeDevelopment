PROGRAM pointinst
  USE channhdvariables
  USE mainparameters
  USE writedata
  USE arrayops
  USE readdata
  USE mpivariables
  USE prelimcalvar
  USE interpoldata
  USE mpi
  IMPLICIT NONE
  REAL*8,ALLOCATABLE,DIMENSION(:)::yp,xp,zp,ypdo,ysdo,ys
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:)::pp
  REAL*8,ALLOCATABLE,DIMENSION(:)::u_i_170_1d,u_i_170_3d,u_i_170_5d
  REAL*8,ALLOCATABLE,DIMENSION(:)::u_i_20_1d,u_i_20_3d,u_i_20_5d,u_i_70_1d,u_i_70_3d,u_i_70_5d
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:)::up,tp
  INTEGER::dt,timestcount,numtimesteps,i,j,k,itime
  integer:: icorlavg,icen,node
  REAL*8::inv,sttime,yprime
  character*2::prosnum
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,mynode,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
  

  CALL readChannelData()
  CALL readPostData(sttime,numtimesteps,dt,icorlavg,icen)
  CALL prelimCal()
  
  ALLOCATE(xp(n3),yp(0:n2),zp(n1),ys(0:n2),ypdo(0:n2do+1),ysdo(0:n2do+1))
  CALL readCoordData(xp,yp,zp)
  CALL intpolCordi(mynode,numprocs,yp,ys,ypdo,ysdo)

  
  allocate(u_i_20_1d(numtimesteps+1),u_i_20_3d(numtimesteps+1),u_i_20_5d(numtimesteps+1))
  allocate(u_i_70_1d(numtimesteps+1),u_i_70_3d(numtimesteps+1),u_i_70_5d(numtimesteps+1))
  allocate(u_i_170_1d(numtimesteps+1),u_i_170_3d(numtimesteps+1),u_i_170_5d(numtimesteps+1))
  timestcount=0
  ! Time loop starts to read data from field.data files
  stime=MPI_WTIME()
  DO ITime= INT(stTime),(INT(StTime)+(NumTimeSteps*Dt)),Dt
     TimeStCount=TimeStcount+1

     ALLOCATE(up(n1,n2do+1,n3,3),tp(n1,n2do+1,n3,1),pp(n1,n2do+1,n3))

     CALL readTempFieldData(mynode,itime,ypdo,ysdo,up,pp,tp)
     node=36/n2do
     if(mynode==node)then
        yprime=36-node*n2do
        u_i_20_1d(timestcount)=up(97,yprime,158,3)
        u_i_20_3d(timestcount)=up(97,yprime,184,3)
        u_i_20_5d(timestcount)=up(97,yprime,210,3)
     end if
     
     node=57/n2do
     if(mynode==node)then
        yprime=57-node*n2do
        u_i_70_1d(timestcount)=up(97,yprime,158,3)
        u_i_70_3d(timestcount)=up(97,yprime,184,3)
        u_i_70_5d(timestcount)=up(97,yprime,210,3)
     end if

     node=75/n2do
     if(mynode==node)then
        yprime=75-node*n2do
        u_i_170_1d(timestcount)=up(97,yprime,158,3)
        u_i_170_3d(timestcount)=up(97,yprime,184,3)
        u_i_170_5d(timestcount)=up(97,yprime,210,3)
     end if
     
     DEALLOCATE(up,tp,pp)
  END DO
  node=36/n2do
  if(mynode==node)then
     call printVector(u_i_20_1d,'u_i_20_1d.dat')
     call printVector(u_i_20_3d,'u_i_20_3d.dat')
     call printVector(u_i_20_5d,'u_i_20_5d.dat')
  end if
     
  node=57/n2do
  if(mynode==node)then
     call printVector(u_i_70_1d,'u_i_70_1d.dat')
     call printVector(u_i_70_3d,'u_i_70_3d.dat')
     call printVector(u_i_70_5d,'u_i_70_5d.dat')
  end if

  node=75/n2do
  if(mynode==node)then
     call printVector(u_i_170_1d,'u_i_170_1d.dat')
     call printVector(u_i_170_3d,'u_i_170_3d.dat')
     call printVector(u_i_170_5d,'u_i_170_5d.dat')
  end if
  etime=MPI_WTIME()
  
  CALL MPI_FINALIZE(ierr)
END PROGRAM pointinst
