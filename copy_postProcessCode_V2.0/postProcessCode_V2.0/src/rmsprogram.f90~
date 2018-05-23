program rms
  USE channhdvariables
  USE mainparameters
  USE writedata
  USE arrayops
  USE readdata
  USE mpivariables
  USE prelimcalvar
  USE interpoldata
  use fluctuation
  USE mpi
  IMPLICIT NONE
  REAL*8,ALLOCATABLE,DIMENSION(:)::yp,xp,zp,ypdo,ysdo,ys
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:)::utavg,vtavg,wtavg,ttavg,ptavg,pp
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:)::urms,vrms,wrms,trms,prms
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:)::uprime,vprime,wprime,tprime,pprime
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:)::uprime2,vprime2,wprime2,tprime2,pprime2
  REAL*8,ALLOCATABLE,DIMENSION(:,:)::utxavg,vtxavg,wtxavg,ttxavg 
  REAL*8,ALLOCATABLE,DIMENSION(:)::utxzavg,vtxzavg,wtxzavg,ttxzavg
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:)::up,tp
  INTEGER::dt,timestcount,numtimesteps,i,j,k,itime
  integer::icorlavg,icen
  REAL*8::inv,sttime
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

  
  ALLOCATE(wtavg(n1,n2do+1,n3),vtavg(n1,n2do+1,n3),&
       utavg(n1,n2do+1,n3),ttavg(n1,n2do+1,n3),ptavg(n1,n2do+1,n3))
  ALLOCATE(wrms(n1,n2do+1,n3),vrms(n1,n2do+1,n3),&
       urms(n1,n2do+1,n3),trms(n1,n2do+1,n3),prms(n1,n2do+1,n3))
  ! Read mean field
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
     allocate(tprime2(n1,n2do+1,n3),wprime2(n1,n2do+1,n3),&
          vprime2(n1,n2do+1,n3),uprime2(n1,n2do+1,n3),pprime2(n1,n2do+1,n3))
     call fluct3Dmean(tp(:,:,:,1),ttavg,tprime)
     call fluct3Dmean(up(:,:,:,1),wtavg,wprime)
     call fluct3Dmean(up(:,:,:,2),vtavg,vprime)
     call fluct3Dmean(up(:,:,:,3),utavg,uprime)
     call fluct3Dmean(up(:,:,:,3),ptavg,pprime)
     DEALLOCATE(up,tp,pp)
     uprime2=uprime*uprime
     vprime2=vprime*vprime
     wprime2=wprime*wprime
     tprime2=tprime*tprime
     pprime2=pprime*pprime
     deallocate(tprime,wprime,vprime,uprime,pprime)
     CALL loopAdd3DAry(uprime2,urms)
     CALL loopAdd3DAry(vprime2,vrms)
     CALL loopAdd3DAry(wprime2,wrms)
     CALL loopAdd3DAry(tprime2,trms)
     CALL loopAdd3DAry(pprime2,prms)
     deallocate(tprime2,wprime2,vprime2,uprime2,pprime2)
  END DO
  etime=MPI_WTIME()
  write(*,*)'time elapsed for data read',etime-stime,'for process',mynode
  write(*,*)'no of time steps',timestcount,'read from process',mynode
  deallocate(wtavg,vtavg,utavg,ttavg)
  inv=1.0/real(timestcount)

  ! Computing time average
  ! write Reynolds normal stresses in a 3D domain
  
  call sendrecv3dwrite(inv*wrms,1,'wsq.dat')
  call sendrecv3dwrite(inv*vrms,2,'vsq.dat')
  call sendrecv3dwrite(inv*urms,3,'usq.dat') 
  
  wrms=sqrt(inv*wrms)
  vrms=sqrt(inv*vrms)
  urms=sqrt(inv*urms)
  trms=sqrt(inv*trms)
  prms=sqrt(inv*prms)
  
  
  ! time mean without spatial averages.
  write(prosnum,'(i2.2)')mynode
  call write3Darray(wrms,'wrms'//trim(prosnum),'unformatted')
  call write3Darray(vrms,'vrms'//trim(prosnum),'unformatted')
  call write3Darray(urms,'urms'//trim(prosnum),'unformatted')  
  call write3Darray(trms,'trms'//trim(prosnum),'unformatted')
  call write3Darray(prms,'prms'//trim(prosnum),'unformatted')
  ! if blowing and suction perturbations are not imposed then spatial averaging
  ! can be performed
  call sendrecv3dwrite(urms,4,'urms.dat')
  call sendrecv3dwrite(vrms,5,'vrms.dat')
  call sendrecv3dwrite(wrms,6,'wrms.dat')
  call sendrecv3dwrite(trms,7,'trms.dat')
  call sendrecv3dwrite(prms,8,'prms.dat')
  if (ibs==0)then
  ! computing average over x-direction
     ALLOCATE(wtxavg(n1,n2do+1),vtxavg(n1,n2do+1),utxavg(n1,n2do+1),ttxavg(n1,n2do+1))
  
     CALL spatialAvg3D(wrms,3,wtxavg)
     CALL spatialAvg3D(vrms,3,vtxavg)
     CALL spatialAvg3D(urms,3,utxavg)
     CALL spatialAvg3D(trms,3,ttxavg)

     DEALLOCATE(wrms,vrms,urms,trms)

     ! computing average over z-direction
     ALLOCATE(wtxzavg(n2do+1),vtxzavg(n2do+1),utxzavg(n2do+1),ttxzavg(n2do+1))

     CALL spatialAvg2D(wtxavg,1,wtxzavg)
     CALL spatialAvg2D(vtxavg,1,vtxzavg)
     CALL spatialAvg2D(utxavg,1,utxzavg)
     CALL spatialAvg2D(ttxavg,1,ttxzavg)

     DEALLOCATE(wtxavg,vtxavg,utxavg,ttxavg)

     !CALL writeflowdata(utxzavg,ttxzavg,ypdo)
     ! print y-coordinate
     IF(mynode==0)THEN
        CALL printVector(yp(1:),'ycoord.dat')
     END IF
     ! print mean values
     stime=MPI_WTIME()
     CALL senRev1dWrite(wtxzavg,1,'wrms_txz.dat')
     CALL senRev1dWrite(vtxzavg,2,'vrms_txz.dat')
     CALL senRev1dWrite(utxzavg,3,'urms_txz.dat')
     CALL senRev1dWrite(ttxzavg,4,'trms_txz.dat')
     etime=MPI_WTIME()
     if (mynode == 0) then
        write(*,*)'time elapsed for data write',etime-stime
     end if
     DEALLOCATE(wtxzavg,vtxzavg,utxzavg,ttxzavg)
  end if
  CALL MPI_FINALIZE(ierr)
END PROGRAM rms
