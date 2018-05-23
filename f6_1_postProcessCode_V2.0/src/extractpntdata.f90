program extpntdata
  USE channhdvariables
  USE mainparameters
  USE writedata
  USE arrayops
  USE readdata
  USE mpivariables
  USE prelimcalvar
  USE interpoldata
  use fluctuation
  use vorticity
  use vortexiden
  USE mpi
  IMPLICIT NONE
  REAL*8,ALLOCATABLE,DIMENSION(:)::yp,xp,zp,ypdo,ysdo,ys
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:)::utavg,vtavg,wtavg,ttavg,ptavg,pp,lambda2
  REAL*8,ALLOCATABLE,DIMENSION(:)::u_20,v_20,w_20,t_20,lam_20
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:)::uprime,vprime,wprime,tprime,pprime
  REAL*8,ALLOCATABLE,DIMENSION(:)::u_70,v_70,w_70,t_70,lam_70
  REAL*8,ALLOCATABLE,DIMENSION(:)::u_170,v_170,w_170,t_170,lam_170 
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:)::up,tp,dudx,dudy,dudz
  INTEGER::dt,timestcount,numtimesteps,i,j,k,itime
  integer::icorlavg,icen,ypri
  REAL*8::inv,sttime,node
  character*2::prosnum
  character*5::time
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,mynode,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
  

  CALL readChannelData()
  CALL readPostData(sttime,numtimesteps,dt,icorlavg,icen)
  CALL prelimCal()
  
  ALLOCATE(xp(n3),yp(0:n2),zp(n1),ys(0:n2),ypdo(0:n2do+1),ysdo(0:n2do+1))
  CALL readCoordData(xp,yp,zp)
  CALL intpolCordi(mynode,numprocs,yp,ys,ypdo,ysdo)

  
  ALLOCATE(wtavg(n1,n2do+1,n3),vtavg(n1,n2do+1,n3),&
       utavg(n1,n2do+1,n3),ttavg(n1,n2do+1,n3),ptavg(n1,n2do+1,n3))
  ! Read mean field
  allocate(u_20(numtimesteps+1),v_20(numtimesteps+1),w_20(numtimesteps+1),&
       t_20(numtimesteps+1),lam_20(numtimesteps+1))
  allocate(u_70(numtimesteps+1),v_70(numtimesteps+1),w_70(numtimesteps+1),&
       t_70(numtimesteps+1),lam_70(numtimesteps+1))
  allocate(u_170(numtimesteps+1),v_170(numtimesteps+1),w_170(numtimesteps+1),&
       t_170(numtimesteps+1),lam_170(numtimesteps+1))

  write(prosnum,'(i2.2)')mynode
  call read3Darray('../../tmean/1500ts/tmean'//trim(prosnum),'unformatted',ttavg)
  call read3Darray('../../tmean/1500ts/wmean'//trim(prosnum),'unformatted',wtavg)
  call read3Darray('../../tmean/1500ts/vmean'//trim(prosnum),'unformatted',vtavg)
  call read3Darray('../../tmean/1500ts/umean'//trim(prosnum),'unformatted',utavg)
  call read3Darray('../../tmean/1500ts/pmean'//trim(prosnum),'unformatted',ptavg)
  timestcount=0

  ! Time loop starts to read data from field.data files
  stime=MPI_WTIME()
  DO ITime= INT(stTime),(INT(StTime)+(NumTimeSteps*Dt)),Dt
     TimeStCount=TimeStcount+1

     ALLOCATE(up(n1,n2do+1,n3,3),tp(n1,n2do+1,n3,1),pp(n1,n2do+1,n3))

     CALL readTempFieldData(mynode,itime,ypdo,ysdo,up,pp,tp)
     allocate(tprime(n1,n2do+1,n3),wprime(n1,n2do+1,n3),&
          vprime(n1,n2do+1,n3),uprime(n1,n2do+1,n3),pprime(n1,n2do+1,n3))
     allocate(dudx(n1,n2do+1,n3,3),dudy(n1,n2do+1,n3,3),dudz(n1,n2do+1,n3,3))
     call velograd(up,xp,ypdo,zp,dudx,dudy,dudz)
     call fluct3Dmean(tp(:,:,:,1),ttavg,tprime)
     call fluct3Dmean(up(:,:,:,1),wtavg,wprime)
     call fluct3Dmean(up(:,:,:,2),vtavg,vprime)
     call fluct3Dmean(up(:,:,:,3),utavg,uprime)
     call fluct3Dmean(pp(:,:,:),ptavg,pprime)
     DEALLOCATE(up,tp,pp)
     
     allocate(lambda2(n1,n2do+1,n3))
     call complambda2(dudx,dudy,dudz,lambda2)
     deallocate(dudx,dudy,dudz)
     node=36/n2do
     if(mynode==node)then
        ypri=36-node*n2do
        u_20(timestcount)=uprime(97,ypri,288)
        v_20(timestcount)=vprime(97,ypri,288)
        w_20(timestcount)=wprime(97,ypri,288)
        t_20(timestcount)=tprime(97,ypri,288)
        lam_20(timestcount)=lambda2(97,ypri,288)
     end if

     node=57/n2do
     if(mynode==node)then
        ypri=57-node*n2do
        u_70(timestcount)=uprime(97,ypri,288)
        v_70(timestcount)=vprime(97,ypri,288)
        w_70(timestcount)=wprime(97,ypri,288)
        t_70(timestcount)=tprime(97,ypri,288)
        lam_70(timestcount)=lambda2(97,ypri,288)
     end if

     node=75/n2do
     if(mynode==node)then
        ypri=75-node*n2do
        u_170(timestcount)=uprime(97,ypri,288)
        v_170(timestcount)=vprime(97,ypri,288)
        w_170(timestcount)=wprime(97,ypri,288)
        t_170(timestcount)=tprime(97,ypri,288)
        lam_170(timestcount)=lambda2(97,ypri,288)
     end if


     deallocate(tprime,wprime,vprime,uprime,pprime,lambda2)

  END DO
  etime=MPI_WTIME()
  write(*,*)'time elapsed for data read',etime-stime,'for process',mynode
  write(*,*)'no of time steps',timestcount,'read from process',mynode
  deallocate(wtavg,vtavg,utavg,ttavg)

   node=36/n2do
  if(mynode==node)then
     call printVector(u_20,'u_20.dat')
     call printVector(v_20,'v_20.dat')
     call printVector(w_20,'w_20.dat')
     call printVector(t_20,'t_20.dat')
     call printVector(lam_20,'lam_20.dat')
  end if

  node=57/n2do
  if(mynode==node)then
     call printVector(u_70,'u_70.dat')
     call printVector(v_70,'v_70.dat')
     call printVector(w_70,'w_70.dat')
     call printVector(t_70,'t_70.dat')
     call printVector(lam_70,'lam_70.dat')
  end if

  node=75/n2do
  if(mynode==node)then
     call printVector(u_170,'u_170.dat')
     call printVector(v_170,'v_170.dat')
     call printVector(w_170,'w_170.dat')
     call printVector(t_170,'t_170.dat')
     call printVector(lam_170,'lam_170.dat')
  end if
  deallocate(u_20,v_20,w_20,t_20,lam_20,u_70,v_70,w_70,t_70,lam_70,u_170,v_170,w_170,&
       t_170,lam_170)
  CALL MPI_FINALIZE(ierr)
END PROGRAM extpntdata
