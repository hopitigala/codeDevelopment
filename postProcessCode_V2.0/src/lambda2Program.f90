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
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:)::pdfhigh,pdflow
  REAL*8,ALLOCATABLE,DIMENSION(:,:)::utxavg,vtxavg,wtxavg,ttxavg
  REAL*8,ALLOCATABLE,DIMENSION(:)::utxzavg,vtxzavg,wtxzavg,ttxzavg
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:)::up,tp,dudx,dudy,dudz
  INTEGER::itime,dt,timestcount,numtimesteps,i,j,k
  REAL*8::inv,sttime,utau,ttau,retau
  character(len=5)::printim
  character(len=2)::prosnum

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,mynode,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

  CALL readChannelData()
  CALL readPostData(sttime,numtimesteps,dt)
  CALL prelimCal()

  ALLOCATE(xp(n3),yp(0:n2),zp(n1),ys(0:n2),ypdo(0:n2do+1),ysdo(0:n2do+1))
  CALL readCoordData(xp,yp,zp)
  CALL intpolCordi(mynode,numprocs,yp,ys,ypdo,ysdo)


  !ALLOCATE(ttavg(n1,n2do+1,n3),trms(n1,n2do+1,n3),tprimeavg(n1,n2do+1,n3))
  !tprimeavg=0.0
  !write(prosnum,'(i2.2)')mynode
  !call read3Darray('/scratch/03390/sdharmar/unperturbed_4pi/postProcess/MPI/tMean/3dmean_1000ts/tmean'//trim(prosnum),'unformatted',ttavg)
  !call read3Darray('/scratch/03390/sdharmar/unperturbed_4pi/postProcess/MPI/rms/2000timesteps/trms'//trim(prosnum),'unformatted',trms)
  !call readflowdata(utau,retau,ttau,'/scratch/03390/sdharmar/unperturbed_4pi/postProcess/MPI/tMean/3dmean_1000ts/ChannelFlowStat.dat')
  !if (mynode==0) then
     
   !  write(*,*)utau,ttau,retau
  !end if
!  if (mynode==0) then
!     open(10,file='trms.dat',status='replace')
!     write(10,'(5f10.6)')(trms(97,5,k),k=200,204)
!     close(10)
!  end if
  timestcount=0
!  if (mynode==0)then
!     open(11,file='lambda2.dat',status='replace')
!     close(11)
!     open(12,file='tprime.dat',status='replace')
!     close(12)
!     open(13,file='tprime_norm.dat',status='replace')
!     close(13)
!  end if
  ! Time loop starts to read data from field.data files
  !allocate(tphighcount(n1,n2do+1,n3),tplowcount(n1,n2do+1,n3),lambdanegcount(n1,n2do+1,n3))
  !tphighcount=0.0
  !tplowcount=0.0
  !lambdanegcount=0.0
  stime=MPI_WTIME()
  do ITime= INT(stTime),(INT(StTime)+(NumTimeSteps*Dt)),Dt
     TimeStCount=TimeStcount+1
     
     write(printim,'(I5.5)')itime

     ALLOCATE(up(n1,n2do+1,n3,3),tp(n1,n2do+1,n3,1),pp(n1,n2do+1,n3))
     
     CALL readTempFieldData(mynode,itime,ypdo,ysdo,up,pp,tp)
     
    ! allocate(tprime(n1,n2do+1,n3))
    ! call fluct3Dmean(tp(:,:,:,1),ttavg,tprime)
     
    ! tprimeavg=tprimeavg+tprime
     allocate(dudx(n1,n2do+1,n3,3),dudy(n1,n2do+1,n3,3),dudz(n1,n2do+1,n3,3))
     
     call velograd(up,xp,ypdo,zp,dudx,dudy,dudz)
     deallocate(up,tp,pp)
     allocate(lambda2(n1,n2do+1,n3))
     call complambda2(dudx,dudy,dudz,lambda2)
     
     deallocate(dudx,dudy,dudz)
 !    if (mynode==0)then
 !       open(11,file='lambda2.dat',status='old', position='append',action='write')
 !       write(11,'(5f10.6)')(lambda2(97,5,k),k=200,204)
 !       close(11)
 !       open(12,file='tprime.dat',status='old', position='append',action='write')
 !       write(12,'(5f10.6)')(tprime(97,5,k),k=200,204)
 !       close(12)
 !    end if
     ! normalizing tprime with local rms value
  !!   tprime=tprime/ttau
  !   if (mynode==0)then
  !      open(13,file='tprime_norm.dat',status='old', position='append',action='write')
  !      write(13,'(5f10.6)')(tprime(97,5,k),k=200,204)
  !      close(13)
  !   end if
     ! conditional pdf
   !  do k=1,n3
   !     do j=1,n2do+1
   !        do i=1,n1
   !           if (lambda2(i,j,k)<0)then
   !              lambdanegcount(i,j,k)=lambdanegcount(i,j,k)+1
   !              if(abs(tprime(i,j,k))>=0.5)then
   !                 tphighcount(i,j,k)=tphighcount(i,j,k)+1
   !              else
   !                 tplowcount(i,j,k)=tplowcount(i,j,k)+1
   !              end if
   !           end if
   !        end do
   !     end do
   !  end do

        
     call sendrecv3dwrite(lambda2,1,'lambda2_'//trim(printim)//'.dat')
     deallocate(lambda2)
  end do
 ! allocate(pdfhigh(n1,n2do+1,n3),pdflow(n1,n2do+1,n3))
 ! pdfhigh=real(tphighcount)/real(lambdanegcount)
 ! pdflow=real(tplowcount)/real(lambdanegcount)
 ! deallocate(tphighcount,tplowcount,lambdanegcount)
  !inv=1.0/timestcount
  !tprimeavg=tprimeavg*inv
 ! deallocate(ttavg,trms)
  !call sendrecv3dwrite(pdfhigh,1,'pdfhigh.txt')
  !deallocate(pdfhigh)
  !call sendrecv3dwrite(pdflow,2,'pdflow.txt')
  !deallocate(pdflow)
  if (mynode==0)then
   !  call printVector(tprimeavg(97,5,:),'tprimeavg.dat')
     call printVector(xp,'xcoord.dat')
     call printVector(yp(1:),'ycoord.dat')
     call printVector(zp,'zcoord.dat')
  end if
  call MPI_FINALIZE(ierr)
end PROGRAM lambda2program
