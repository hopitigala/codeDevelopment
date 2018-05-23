PROGRAM timemean
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
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:)::utavg,vtavg,wtavg,ttavg,ptavg,pp
  REAL*8,ALLOCATABLE,DIMENSION(:,:)::utxavg,vtxavg,wtxavg,ttxavg,ptxavg 
  REAL*8,ALLOCATABLE,DIMENSION(:)::utxzavg,vtxzavg,wtxzavg,ttxzavg,ptxzavg
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:)::up,tp
  INTEGER::dt,timestcount,numtimesteps,i,j,k,itime
  integer:: icorlavg,icen
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

  timestcount=0
  wtavg=0.0
  vtavg=0.0
  utavg=0.0
  ttavg=0.0
  ptavg=0.0
  ! Time loop starts to read data from field.data files
  stime=MPI_WTIME()
  DO ITime= INT(stTime),(INT(StTime)+(NumTimeSteps*Dt)),Dt
     TimeStCount=TimeStcount+1

     ALLOCATE(up(n1,n2do+1,n3,3),tp(n1,n2do+1,n3,1),pp(n1,n2do+1,n3))

     CALL readTempFieldData(mynode,itime,ypdo,ysdo,up,pp,tp)
    
     CALL loopAdd3DAry(up(:,:,:,1),wtavg)
     CALL loopAdd3DAry(up(:,:,:,2),vtavg)
     CALL loopAdd3DAry(up(:,:,:,3),utavg)
     CALL loopAdd3DAry(tp(:,:,:,1),ttavg)
     call loopAdd3Dary(pp,ptavg)
     DEALLOCATE(up,tp,pp)
  END DO
  etime=MPI_WTIME()
  write(*,*)'time elapsed for data read',etime-stime,'for process',mynode
  write(*,*)'no of time steps',timestcount,'read from process',mynode

  inv=1.0/real(timestcount)

  ! Computing time average


  wtavg=inv*wtavg
  vtavg=inv*vtavg
  utavg=inv*utavg
  ttavg=inv*ttavg
  ptavg=inv*ptavg

  ! time mean without spatial averages.
  write(prosnum,'(i2.2)')mynode
  call write3Darray(wtavg,'wmean'//trim(prosnum),'unformatted')
  call write3Darray(vtavg,'vmean'//trim(prosnum),'unformatted')
  call write3Darray(utavg,'umean'//trim(prosnum),'unformatted')  
  call write3Darray(ttavg,'tmean'//trim(prosnum),'unformatted')
  call write3Darray(ptavg,'pmean'//trim(prosnum),'unformatted')

  ! if blowing is not activated spatial averages can be taken in both
  ! -x and -y directions
  call sendrecv3dwrite(utavg,1,'umean_3d.dat')
  call sendrecv3dwrite(ptavg,2,'pmean_3d.dat')
  if (ibs==0)then
     ! computing average over x-direction
     ALLOCATE(wtxavg(n1,n2do+1),vtxavg(n1,n2do+1),utxavg(n1,n2do+1),ttxavg(n1,n2do+1),ptxavg(n1,n2do+1))
  
     CALL spatialAvg3D(wtavg,3,wtxavg)
     CALL spatialAvg3D(vtavg,3,vtxavg)
     CALL spatialAvg3D(utavg,3,utxavg)
     CALL spatialAvg3D(ttavg,3,ttxavg)
     CALL spatialAvg3D(ptavg,3,ptxavg)
     
     DEALLOCATE(wtavg,vtavg,utavg,ttavg,ptavg)

  ! computing average over z-direction
     ALLOCATE(wtxzavg(n2do+1),vtxzavg(n2do+1),utxzavg(n2do+1),ttxzavg(n2do+1),ptxzavg(n2do+1))

     CALL spatialAvg2D(wtxavg,1,wtxzavg)
     CALL spatialAvg2D(vtxavg,1,vtxzavg)
     CALL spatialAvg2D(utxavg,1,utxzavg)
     CALL spatialAvg2D(ttxavg,1,ttxzavg)
     CALL spatialAvg2D(ptxavg,1,ptxzavg)

     DEALLOCATE(wtxavg,vtxavg,utxavg,ttxavg,ptxavg)

     CALL writeflowdata(utxzavg,ttxzavg,ypdo)
     ! print y-coordinate
     IF(mynode==0)THEN
        CALL printVector(yp(1:),'ycoord.dat')
     END IF
     ! print mean values
     stime=MPI_WTIME()
     CALL senRev1dWrite(wtxzavg,1,'w_txz_mean.dat')
     CALL senRev1dWrite(vtxzavg,2,'v_txz_mean.dat')
     CALL senRev1dWrite(utxzavg,3,'u_txz_mean.dat')
     CALL senRev1dWrite(ttxzavg,4,'t_txz_mean.dat')
     CALL senRev1dWrite(ptxzavg,4,'p_txz_mean.dat')
     etime=MPI_WTIME()
     if (mynode == 0) then
        write(*,*)'time elapsed for data write',etime-stime
     end if
     DEALLOCATE(wtxzavg,vtxzavg,utxzavg,ttxzavg,ptxzavg)
  end if
  
  CALL MPI_FINALIZE(ierr)
END PROGRAM timemean
