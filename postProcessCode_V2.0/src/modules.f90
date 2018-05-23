!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!NumProcs - Number of processors in the simulation
!iz - Center grid number of holes in z-direction
!J0 - grid point number at reference y locations in correlation calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MainParameters
  IMPLICIT NONE
  REAL*8 ,PARAMETER    :: Pi=3.14159265359
!  INTEGER , PARAMETER  :: NumProcs=24
  INTEGER, PARAMETER,DIMENSION(5) :: iz=(/28,62,96,130,164/)
!  INTEGER,PARAMETER,DIMENSION(6) :: J0=(/36,57,75,119,138,158/)
  CHARACTER(*),PARAMETER,DIMENSION(11)::labl=(/'ww','vv','uu','tt','uv','ut','vt','wt','pp','pt','uw'/)
END MODULE MainParameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This Module includes declarations of Main variables used in the program. Their
! descriptions are given below.
!
! Mynode          - Processor number of the current task
! l               - An integer variale for indexing
! TimeStCount     - Count the no. of timesteps to calculate the mean
! N1              - No. of grid points in z-dir (spanwise)
! N2              - No. of grid points in y-dir (wall normal)
! N3              - No. of grid points in x-dir (streamwise)
! N1m             - N1-1
! N2m             - N2-1
! N3m             - N3-1
! N2do            - N2-1/NumProcs ( No. of grid points per processor)
! Islip           -
! Iflu1           -
! Iflun2          -
! Nstop           -
! Multim          -
! Nread           -
! Icfl            -
! Istr2           -
! Mpuny           - number of grid points below yp2=-1 till yp2=-1+bhc
! Ipassc          - When = 0 :passive scalar OFF; =1 : passive scalar ON
! Nps             -
! Npsc            -
! Ibody           -
! Ibs             -
! Ngr             -
! Npsl            -
! ReyNum          - Reynolds Number
! StTime          - Starting time to read field.data file (specified in post.d file)
! Cflc            -
! Dtt             - Time step to advace the field calculation
! Tfin            - Final time step to calculate field
! Timav           - Time step to start averaging the field
! Ros2            -
! Vper            -
! Omtres          -
! Tpin            -
! Tprin           -
! Bh              -
! Ampl            - Amplitude of blowing perturbation
! Angle           - Angle of blowing perturbation with wall normal dir
! Freq            - Frequency of blowing perturbation
! Str2            - Constant that is required to calculate y-dir grid distribusion
! Bhc             -
! Lx1d            - z-dir dimension specified in channht.d file
! Lx2d            - y-dir dimension specified in channht.d file
! Lx3d            - x-dir dimension specified in channht.d file
! Lx1             - Lx1d*Pi
! Lx2             - Lx2d
! Lx3             - Lx3d*Pi
! I,J,K,JJ        - Indices variables
! ITime           - Index for time iteration loop
! Dt              - Time interval of field.data file recorded
! NumTimeSteps    - Number of timesteps that field.data was collected
! totNumTimeSteps - Total number of time steps specified in channht.d file
! writeTime       - Time instant that field.data file is written
! PnTime          - Time stamp printed in field.data file name
! ProsNum         - Processor number that outputs the field.data file
! FileName        - character variable to take field.data file
! Utau            - Friction Velocity
! Nu              - Kinematic Viscosity
! Ttau            - Friction temperatutre (dtheta/dy|y=0)*Nu/PraNum*Utau
! PraNum          - Prandtl number(=0.71)
! ypJ0            - yplus values at J0 locations
! compModes       - index to notify POD mode parallel write or read (used in podProgram_trunc.f90)
!flxmodeCal       - index to decide whether prog need to compute flux modes or velocity modes (used in podProgram_trunc.f90, specif\ied in post.d file)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE channhdvariables
  IMPLICIT NONE
  INTEGER ::n1,n2,n3
  REAL*8  ::lx1d,lx2d,lx3d
  INTEGER ::istr2,mpuny,totnumtimesteps
  REAL*8  ::bhc,str2,reynum,dtt,cflc
  INTEGER ::nstop,multim,nread,icfl
  REAL*8  ::tfin,timav,ros2,vper,omtres,tpin,tprin
  INTEGER ::islip,iflu1,iflun2
  INTEGER ::ipassc,nps,npsc,ibody,ibs,ngr,npsl
  REAL*8  ::bh,ampl,angle,freq
  INTEGER,DIMENSION(3)::ibou
  REAL*8,DIMENSION(2)::shm
END MODULE channhdvariables


MODULE mpivariables
  USE mpi
  IMPLICIT NONE
  INTEGER  :: mynode,numprocs
  INTEGER*4:: status(MPI_STATUS_SIZE)
  INTEGER  :: ierr
  REAL*8:: stime, etime
END MODULE mpivariables

MODULE prelimcalvar
  IMPLICIT NONE 
  INTEGER:: n1m,n2m,n3m,n2do,n3do,filstopros
  REAL*8:: lx1,lx2, lx3
END MODULE prelimcalvar
