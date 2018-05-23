

MODULE writedata
  IMPLICIT NONE
  
CONTAINS

  SUBROUTINE printVector(Vector,filename)
    IMPLICIT NONE
    REAL*8,INTENT(IN),DIMENSION(:)::Vector
    CHARACTER(LEN=*),INTENT(IN)::filename
    INTEGER::I
    OPEN(11,file=filename)
    DO I=1,size(Vector,1)
       WRITE(11,'(F12.6)')Vector(I)
    END DO
    CLOSE(11)
  END SUBROUTINE printVector
  
  
  SUBROUTINE print2DMatx(Matx,filename)
    IMPLICIT NONE
    
    REAL*8,INTENT(IN),dimension(:,:)::Matx
    CHARACTER(LEN=*),INTENT(IN)::filename
    INTEGER::rows,colms
    CHARACTER*4::num_col
    INTEGER::i,j
    
    colms=size(Matx,2)
    rows=size(Matx,1)

    IF(colms<10)THEN
       WRITE(num_col,'(I1.1)')colms
    ELSEIF(colms<100)THEN
       WRITE(num_col,'(I2.2)')colms
    ELSEIF(colms<1000)THEN
       WRITE(num_col,'(I3.3)')colms
    ELSE
       WRITE(num_col,'(I4.4)')colms
    END IF

    OPEN(10,file=filename)
    DO i=1,rows
       write(10,'('//TRIM(num_col)//'E14.5e2)')(Matx(i,j),j=1,colms)
    END DO
    CLOSE(10)
  END SUBROUTINE print2DMatx

  SUBROUTINE senRev1dwrite(ary1d,tag,filname)
    USE mpi
    USE mpivariables
    IMPLICIT NONE
    REAL*8,INTENT(IN),DIMENSION(:)::ary1d
    REAL*8,ALLOCATABLE,DIMENSION(:)::globary
    CHARACTER(LEN=*),INTENT(IN)::filname
    INTEGER,INTENT(IN)::tag
    INTEGER::j,d1,p,jj
    
    d1=size(ary1d,1)
    
!    IF(mynode/=0)THEN
!       CALL MPI_SEND(ary1d,d1,MPI_REAL8,0,tag,MPI_COMM_WORLD,ierr)
!    ELSE
    if (mynode==0) then
       allocate (globary(numprocs*(d1-1)+1))
    end if
       call MPI_GATHER(ary1d,d1-1,MPI_REAL8,globary,d1-1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!       DO j=1,d1
!          jj=j
!          globary(jj)=ary1d(j)
!       END DO
!       DO p=1,numprocs-1
!          CALL MPI_RECV(ary1d,d1,MPI_REAL8,p,tag,MPI_COMM_WORLD,status,ierr)
!          DO j=1,d1
!             jj=j+p*(d1-1)
!             globary(jj)=ary1d(j)
!          END DO
!       END DO
       if (mynode==0)then
          call printVector(globary,filname)
          deallocate (globary)
       end if
    
  END SUBROUTINE senRev1dwrite
  
  subroutine sendrecv3dwrite(ary3d,tag,filename)
    use mpi
    use mpivariables
    implicit none
    real*8,          intent(in),dimension(:,:,:):: ary3d
    integer,         intent(in)                 :: tag
    character(len=*),intent(in)                 :: filename
    integer                                     :: i,j,k,d1,d2,d3,p,jj
    real*8,allocatable,dimension(:,:,:)         :: globary
    character(len=*),parameter                  :: fmt='(ES13.5E2)'
    d1=size(ary3d,1)
    d2=size(ary3d,2)
    d3=size(ary3d,3)
    
    if (mynode/=0)then
       call MPI_SEND(ary3d,d1*d2*d3,MPI_REAL8,0,tag,MPI_COMM_WORLD,ierr)
    else
       allocate(globary(d1,numprocs*(d2-1)+1,d3))
       
       do k=1,d3
          do j=1,d2
             jj=j
             do i=1,d1
                globary(i,jj,k)=ary3d(i,j,k)
             end do
          end do
       end do
       do p=1,(numprocs-1)
          call MPI_RECV(ary3d,d1*d2*d3,MPI_REAL8,p,tag,MPI_COMM_WORLD,status,ierr)
          do k=1,d3
             do j=1,d2
                jj=j+p*(d2-1)
                do i=1,d1
                   globary(i,jj,k)=ary3d(i,j,k)
                end do
             end do
          end do
       end do
       open(10,file=filename)
       do k=1,d3
          do j=1,(numprocs*(d2-1)+1)
             do i=1,d1
                write(10,fmt)globary(i,j,k)
             end do
          end do
       end do
       deallocate(globary)
    end if
  end subroutine sendrecv3dwrite

  subroutine sendrecv2dwrite(ary2d,tag,iplane,filename)
    use mpi
    use mpivariables
    implicit none
    real*8,          intent(in),dimension(:,:)  :: ary2d
    integer,         intent(in)                 :: tag,iplane
    character(len=*),intent(in)                 :: filename
    integer                                     :: i,j,k,d1,d2,p,jj,ii
    real*8,allocatable,dimension(:,:)           :: globary
    character(len=*),parameter                  :: fmt='(ES13.5E2)'
    d1=size(ary2d,1)
    d2=size(ary2d,2)
    
    if (mynode/=0)then
       call MPI_SEND(ary2d,d1*d2,MPI_REAL8,0,tag,MPI_COMM_WORLD,ierr)
    else
       if(iplane==1) then!xy plane
          allocate(globary(numprocs*(d1-1)+1,d2))
       
          do j=1,d2
             do i=1,d1
                ii=i
                globary(ii,j)=ary2d(i,j)
             end do
          end do
          do p=1,(numprocs-1)
             call MPI_RECV(ary2d,d1*d2,MPI_REAL8,p,tag,MPI_COMM_WORLD,status,ierr)
             do j=1,d2                
                do i=1,d1
                   ii=i+p*(d1-1)
                   globary(ii,j)=ary2d(i,j)
                end do
             end do
          end do
          open(10,file=filename)
          do j=1,d2
             do i=1,(numprocs*(d1-1)+1)
                write(10,fmt)globary(i,j)
             end do
          end do
          deallocate(globary)
       elseif(iplane==2) then!yz-plane
          allocate(globary(d1,numprocs*(d2-1)+1))
       
          do j=1,d2
             jj=j
             do i=1,d1
                globary(i,jj)=ary2d(i,j)
             end do
          end do
          do p=1,(numprocs-1)
             call MPI_RECV(ary2d,d1*d2,MPI_REAL8,p,tag,MPI_COMM_WORLD,status,ierr)
             do j=1,d2
                jj=j+p*(d2-1)
                do i=1,d1                   
                   globary(i,jj)=ary2d(i,j)
                end do
             end do
          end do
          open(10,file=filename)
          do j=1,(numprocs*(d2-1)+1)
             do i=1,d1
                write(10,fmt)globary(i,j)
             end do
          end do
          deallocate(globary)
       end if         
    end if
  end subroutine sendrecv2dwrite
  
  SUBROUTINE writeflowdata(umean,tmean,ypdo)
    USE channhdvariables
    USE mpivariables
    IMPLICIT NONE
    REAL*8,INTENT(IN),DIMENSION(:)::umean,tmean
    REAL*8,INTENT(IN),DIMENSION(0:)::ypdo
    REAL*8::nu,utau,ttau,pranum
    IF (MyNode == 0)THEN
       nu=1.0/reynum
       pranum=shm(1)
       utau=SQRT(abs(nu*(umean(2)-umean(1))/(ypdo(2)-ypdo(1))))
       Write(*,*)'ypdo(1)= ', ypdo(1)
       write(*,*)'ypdo(2)= ', ypdo(2)
       ! Friction temperature
       
       ttau=(nu*((tmean(2)-tmean(1))/(ypdo(2)-ypdo(1))))/(pranum*utau)

       OPEN(10,File='ChannelFlowStat.dat')
       write(10,*)'Friction velocity = ',utau
       write(10,*)'Friction Re. No = ',utau/nu
       write(10,*)'Prandtl Number = ',pranum
       write(10,*)'Friction temperature = ',ttau
       CLOSE(10)
    END IF  
  END SUBROUTINE writeflowdata
  
  subroutine write3Darray(array3d,filename,fmt)
    implicit none
    real*8,intent(in),dimension(:,:,:)::array3d
    character(len=*),intent(in)       ::filename,fmt
    integer                           ::d1,d2,d3,i,j,k
    
    d1=size(array3d,1)
    d2=size(array3d,2)
    d3=size(array3d,3)
    open(11,file=filename,form=fmt)
    do k=1,d3
       do j=1,d2
          do i=1,d1
             write(11)array3d(i,j,k)
          end do
       end do
    end do
    close(11)
  end subroutine write3Darray
    
END MODULE writedata
  
  
  
  
