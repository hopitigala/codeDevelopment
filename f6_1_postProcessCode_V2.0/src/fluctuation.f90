!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  This module contains subroutines that are used to compute fluctuations    !!
!!                                                                            !! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE fluctuation
  IMPLICIT NONE

CONTAINS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! This subroutine computes instantaneous fluctuations when a 3D instantaneous!!
  !! field and 3d mean feild is supplied.                                       !!
  !! INPUTS:                                                                    !!
  !! real*8:: instAry (3d array), mean3dAry                                     !!
  !! OUTPUTS:                                                                   !! 
  !! real*8:: flucAry (3d array)                                                !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE fluct3Dmean(instAry,mean3dAry,flucAry)
    IMPLICIT NONE
    REAL*8,INTENT(IN),DIMENSION(:,:,:)::instAry,mean3dAry
    REAL*8,INTENT(OUT),DIMENSION(:,:,:)::flucAry
    
    flucAry=instAry-mean3dAry
    
  END SUBROUTINE fluct3Dmean

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! This subroutine computes instantaneous fluctuations when a 3D instantaneous!!
  !! field and 1d mean feild is supplied.                                       !!
  !! INPUTS:                                                                    !!
  !! real*8:: instAry (3d array), mean1dAry                                     !!
  !! OUTPUTS:                                                                   !! 
  !! real*8:: flucAry (3d array)                                                !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE fluct1Dmean(instAry,mean1dAry,flucAry)
    IMPLICIT NONE
    REAL*8,INTENT(IN),DIMENSION(:,:,:)::instAry
    REAL*8,INTENT(IN),DIMENSION(:)::mean1dAry
    REAL*8,INTENT(OUT),DIMENSION(:,:,:)::flucAry
    INTEGER::i,j,k,d1,d2,d3

    ! find dimensions of the array
    
    d1=SIZE(instAry,1)
    d2=SIZE(instAry,2)
    d3=SIZE(instAry,3)

    DO k=1,d3
       DO j=1,d2
          DO i=1,d1
             flucAry(i,j,k)=instAry(i,j,k)-mean1dAry(j)
          END DO
       END DO
    END DO
    
  END SUBROUTINE fluct1Dmean

END MODULE fluctuation
