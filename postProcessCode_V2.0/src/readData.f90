MODULE readdata
  IMPLICIT NONE
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  This file includes subroutines to read data from files                    !!
!!  The code is written by Dr. Suranga Dharmarathne of Texas Tech University  !!
!!  This code is for post processing DNS channel flow data taken from         !!
!!  the DNS code developed by Dr. Stefano Leonardi at UT Dallas               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! the SUBROUTINE reads data from the channht.d file, which contains all the !!
!!! the data that are required to calculate the solution by using the code    !!
!!! See modules.f90 file for the description of variables.                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE readChannelData()
    USE channhdvariables
    IMPLICIT NONE
    INTEGER::l
    OPEN(10,file='../channht.d')
    READ(10,*) N1,N2,N3             ! number of nodes n1,n2,n3=z,y,x-directions
    READ(10,*) Lx1d,Lx2d,Lx3d    ! size of the channel in each direction (this number*Pi)
    READ(10,*) Str2,Istr2,Mpuny,Bhc ! distorsion with in channel, channel type (pg 79 L's thesis)
    READ(10,*) ReyNum,Dtt,totNumTimeSteps  !Reynolds number, time step, number of timesteps.
    READ(10,*) Nstop,Multim,Nread
    READ(10,*) Icfl,Cflc,Tpin,Tprin,Tfin,Timav
    READ(10,*) Ros2,Vper,Omtres
    READ(10,*) Islip,Iflu1,Iflun2
    READ(10,*) (Ibou(l),l=1,3)
    READ(10,*) Ipassc,Npsl,Npsc,(Shm(nps),nps=1,Npsc)
    READ(10,*) Ibody,Bh,Ngr
    READ(10,*) Ibs,Ampl,Angle,Freq
    Close(10)
  END SUBROUTINE readChannelData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This subroutine read data in post.d file which contains parameters for    !!!
!! post processing data                                                      !!!
!! OUTPUTS:                                                                  !!!
!!     StTime      : Start time of the post processing                       !!!
!!     NumTimeSteps: Number of time steps used for the post processing       !!!
!!     Dt          : time step size                                          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE readPostData(sttime,numtimesteps,dt,icorlavg,icen)
    IMPLICIT NONE
    REAL*8,INTENT(OUT):: sttime
    INTEGER,INTENT(OUT):: numtimesteps,dt
    integer,intent(out)::icorlavg,icen
    OPEN(11,file='post.d')
    READ(11,*) sttime,numtimesteps,dt
    READ(11,*) icorlavg,icen
    Close(11)
  END SUBROUTINE readPostData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     iplane      : xy=1,yz=2,xz=3 selected 2d plane to compute pod         !!!
!!     x0          : the point at which the plane was chosen                 !!!
!!     lsmlim      : number of modes that defines LSMs                       !!!
!!     ireyst      : uv=3, v2=2, u2=1                                        !!!
!!     iscfx       : (=1 compute scalar flux; =0 no scalar flux )            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readpoddata(iplane,x0,lsmlim,ireyst,iscfx,icorl,var1,var2,j0,k0,inst)
  implicit none
  integer,intent(out)::iplane,x0,lsmlim,ireyst,inst,iscfx
  integer,intent(out)::icorl,var1,var2,j0,k0
  open(12,file='poddata.dat')
  read(12,*)iplane,x0,lsmlim,ireyst,inst,iscfx
  read(12,*)icorl,var1,var2,j0,k0
  close(12)
end subroutine readpoddata

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This subroutine computes some of the variables required for other          !!
!! subroutines                                                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE prelimcal()
  USE mainparameters
  USE channhdvariables
  USE prelimcalvar
  USE mpivariables
  IMPLICIT NONE
  n1m=n1-1
  n2m=n2-1
  n3m=n3-1
  n2do=n2m/numprocs
  lx1=lx1d*Pi
  lx2=lx2d
  lx3=lx3d*Pi
END SUBROUTINE prelimcal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This subroutine reads data from the field.data file.                       !!
!! Change the Filename according to the location of the file                  !!
!! INPUTS::                                                                   !!
!!       mynode      - processor number                                       !!
!!       itime       - time step number                                       !!
!!                                                                            !!
!!                                                                            !!
!! OUTPUTS::                                                                  !! 
!!       up       - temporary velocity field array cell center                !! 
!!       pp       - temporary pressure field array cell center                !!
!!       tp       - tempprary scalar   field array cell center                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


SUBROUTINE readTempFieldData(mynode,itime,ypdo,ysdo,up,pp,tp)
  USE channhdvariables
  USE prelimcalvar
  USE interpoldata
  
  IMPLICIT NONE
  INTEGER,INTENT(IN)::mynode,itime
  REAL*8,INTENT(IN),DIMENSION(0:)::ypdo,ysdo
  REAL*8,INTENT(OUT)::up(:,:,:,:),pp(:,:,:),tp(:,:,:,:)
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:)::utemp,ttemp
  REAL*8,ALLOCATABLE,DIMENSION(:,:,:)::ptemp
  INTEGER::i,j,k,l,d1,d2,d3,allocst
  CHARACTER*2::prosnum
  CHARACTER*5::pntime
  CHARACTER*100::filename
  CHARACTER(*),PARAMETER ::FMT1='(I5.5)'
  CHARACTER(*),PARAMETER ::FMT2='(I2.2)'
  WRITE(PnTime,FMT1)itime
  WRITE(ProsNum,FMT2)mynode

  ! Allocating arrays to store field data. n1m,n3m, and n2do is declared in module 
  ! channhtdvariables
  
  ALLOCATE(utemp(n1m,0:n2do+1,n3m,3),ttemp(n1m,0:n2do+1,n3m,1),&
       ptemp(n1m,0:n2do+1,n3m),STAT=allocst)
  IF(allocst /= 0) STOP "***Not enough memory***"
  
  d1=size(utemp,1)
  d2=size(utemp,2)-1
  d3=size(utemp,3)
  !write(*,*)'Read data from Procs : ',MyNode
  filename = '../../../../field/field'//TRIM(ProsNum)//'/field.data'//TRIM(ProsNum)//'_'//TRIM(PnTime)
  !Write(*,*)'Read file :',FileName
  OPEN(12,FILE=filename,FORM='UnFormatted')
  READ(12)!N1m,N2do,N3m,Npsc
  READ(12)!writeTime,ReyNum
  IF(ipassc.EQ.1)THEN
     READ(12)(((utemp(i,j,k,1),i=1,d1), j=0,d2),k=1,d3), &
          (((utemp(i,j,k,2),i=1,d1), j=0,d2),k=1,d3), &
          (((utemp(i,j,k,3),i=1,d1), j=0,d2),k=1,d3), &
          (((ptemp(i,j,k),i=1,d1), j=0,d2),k=1,d3),  &
          ((((ttemp(i,j,k,l),i=1,d1), j=0,d2),k=1,d3),l=1,npsc)
  ELSE
     READ(12)(((utemp(i,j,k,1),i=1,d1), j=0,d2),k=1,d3), &
          (((utemp(i,j,k,2),i=1,d1), j=0,d2),k=1,d3), &
          (((utemp(i,j,k,3),i=1,d1), j=0,d2),k=1,d3), &
          (((ptemp(i,j,k),i=1,d1), j=0,d2),k=1,d3)
  ENDIF
  CLOSE(12)
!  CALL printtempData(pntime,prosnum,utemp,ttemp,ptemp)
  CALL intpolFieldData(utemp,ptemp,ttemp,ysdo,ypdo,up,pp,tp)
  DEALLOCATE(utemp,ttemp,ptemp)

END SUBROUTINE readTempFieldData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This subroutine reads coord32.dat and coord13.dat files, which contain     !!
!! cordinate data of x,y and z,x respectively                                 !!
!!                                                                            !!
!! OUTPUTS::                                                                  !! 
!!       xp       - x -cord data                                              !!
!!       yp       - y- cord data                                              !!
!!       zp       - z -cord data                                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE readCoordData(xp,yp,zp)
  IMPLICIT NONE

  REAL*8,INTENT(OUT),DIMENSION(:)::zp,xp
  REAL*8,INTENT(OUT),DIMENSION(0:)::yp
  INTEGER::i,j,k,d1,d2,d3
  OPEN(13,FILE='../coord32.dat')
  READ(13,*)d3,d2
  READ(13,*)((xp(k),k=1,d3),j=1,d2),&
       ((yp(j),k=1,d3),j=1,d2)
  CLOSE(13)

  OPEN(14,FILE='../coord13.dat')
  READ(14,*)d1,d3
  READ(14,*)((zp(i),i=1,d1),k=1,d3),&
       ((xp(k),i=1,d1),k=1,d3)
  CLOSE(14)

END SUBROUTINE readCoordData

SUBROUTINE PrinttempData(pntime,prosnum,u,t,p)
!  USE MainVariables
!  USE MainParameters
!  USE PointFieldArrays
  IMPLICIT NONE
  
  character*5,intent(in)::pntime
  character*2,intent(in)::prosnum
  real*8,intent(in),dimension(:,0:,:,:)::u,t
  real*8,intent(in),dimension(:,0:,:)::p
  integer::j,d2
  character*2::node
  
  CHARACTER(*),PARAMETER :: FMT0='(A8,I3)'
  CHARACTER(*),PARAMETER :: FMT1='(6A15)'
  CHARACTER(*),PARAMETER :: FMT2='(I15,5ES15.5E2)'
 
  d2=size(u,2)-1
  OPEN(16, FILE='TempCheck-'//TRIM(PnTime)//'_'//trim(prosnum)//'.dat')
!  WRITE(16,FMT0)'ProsNo',MyNode
  WRITE(16,FMT1)'J','W','V','U','P','T'
  DO J=0,d2
     WRITE(16,FMT2)J,U(97,J,10,1),U(97,J,10,2),U(97,J,10,3),&
          P(97,J,10),T(97,J,10,1)
  END DO
  CLOSE(16)
END SUBROUTINE PrinttempData


SUBROUTINE read2DMatx(rows,colms,filename,Matx)
  IMPLICIT NONE
  INTEGER,INTENT(IN)::rows,colms
  REAL*8,INTENT(OUT)::Matx(rows,colms)
  CHARACTER(LEN=*),INTENT(IN)::filename
  CHARACTER*4::num_col
  INTEGER::i,j

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
     READ(10,'('//TRIM(num_col)//'E14.5e2)')(Matx(i,j),j=1,colms)
  END DO
  CLOSE(10)
END SUBROUTINE read2DMatx

subroutine readvector(filename,vector)
  implicit none
  character(len=*),intent(in)   ::filename
  real*8,dimension(:),intent(out)::vector
  integer                       ::i
  open(11,file=filename)
  do i=1,size(vector,1)
     read(11,'(f12.6)')vector(i)
  end do
  close(11)
end subroutine readvector

subroutine read3Darray(filename,fmt,array3d)
  implicit none
  character(len=*),intent(in)         :: filename,fmt
  real*8,intent(out),dimension(:,:,:) :: array3d  
  integer                             :: d1,d2,d3,i,j,k
  
  d1=size(array3d,1)
  d2=size(array3d,2)
  d3=size(array3d,3)
  
  open(11,file=filename,form=fmt)
  do k=1,d3
     do j=1,d2
        do i=1,d1
           read(11)array3d(i,j,k)
        end do
     end do
  end do
  close(11)
end subroutine read3Darray

subroutine read4Darray(filename,fmt,array4d)
  implicit none
  character(len=*),intent(in)         :: filename,fmt
  real*8,intent(out),dimension(:,:,:,:) :: array4d  
  integer                             :: d1,d2,d3,d4,i,j,k,m
  
  d1=size(array4d,1)
  d2=size(array4d,2)
  d3=size(array4d,3)
  d4=size(array4d,4)
  
  open(12,file=filename,form=fmt)
  do m=1,d4
     do k=1,d3
        do j=1,d2
           do i=1,d1
              read(12)array4d(i,j,k,m)
           end do
        end do
     end do
  end do
  close(12)
end subroutine read4Darray

subroutine readflowdata(utau,retau,ttau,flowfile)
  implicit none
  real*8,intent(out)::utau,retau,ttau
  character(len=*),intent(in)::flowfile
  character(len=19)::fricvel
  character(len=17)::reytau
  character(len=16)::prandtl
  character(len=22)::frictemp
  real*8:: pranum
  open(10,file=flowfile)
  read(10,*)utau
  read(10,*)retau
  read(10,*)pranum
  read(10,*)ttau
  close(10)
end subroutine readflowdata

subroutine readcorldata(var1,var2,z0,j0,k0)
  implicit none
 ! integer,intent(in)::icorlavg,icen
  integer,intent(out)::var1,var2,j0,k0
  integer,intent(out),dimension(:)::z0
  integer:: i,d1
  d1=size(z0)
  open(15,file='coreldata.dat')
  read(15,*) var1,var2,j0,k0
  read(15,*) (z0(i),i=1,d1)
  !if(icen==1)then
  !   if(icorlavg==1)then
  !      read(15,*) (z0(i),i=1,d1)
  !      read(15,*)
  !   else
  !      read(15,*)
  !      read(15,*)(zo(i),i=1,d2)
  !else
  !   read(15,*) 
  !   read(15,*)
  !   read(15,*) (z0(i),i=1,d1)
  !end if
  close(15)
end subroutine readcorldata

END MODULE readdata




