!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This module contains subroutines that do array operations        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE arrayops
  IMPLICIT NONE
  
CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! This subroutine adds 3d array in a loop. This can be used in time!!
    !! loop additions.                                                  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
    SUBROUTINE loopAdd3DAry(aryin,sumAry)
      IMPLICIT NONE
      REAL*8,INTENT(IN),DIMENSION(:,:,:)::aryin
      REAL*8,INTENT(INOUT),DIMENSION(:,:,:)::sumAry
      
      sumAry=sumAry+aryin
    END SUBROUTINE loopAdd3DAry

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! This subroutine multiply two 3d arrays.                          !!
    !!                                                                  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
    SUBROUTINE eleMult3Dary(ary1in,ary2in,aryout)
      IMPLICIT NONE
      REAL*8,INTENT(IN),DIMENSION(:,:,:)::ary1in,ary2in
      REAL*8,INTENT(OUT),DIMENSION(:,:,:)::aryout
      
      aryout=ary1in*ary2in
    END SUBROUTINE eleMult3Dary
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! subroutine finds the average over direction defined by (dim) of  !!
    !! 3D array                                                         !!
    !! INPUTS:                                                          !!
    !!     REAL*8:: ary3d                                               !!
    !!     INTEGER:: dim (1,2, or 3)                                    !!
    !! OUTPUT:                                                          !!
    !!     REAL*8:: avg2d                                               !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    SUBROUTINE spatialAvg3D(ary3d,dim,avg2d)
      IMPLICIT NONE
      REAL*8,INTENT(IN),DIMENSION(:,:,:)::ary3d
      REAL*8,INTENT(OUT),DIMENSION(:,:)::avg2d
      INTEGER,INTENT(IN)::dim
      INTEGER::i,j,k,d1,d2,d3
      REAL*8:: summ
      
      d1=size(ary3d,1)
      d2=size(ary3d,2)
      d3=size(ary3d,3)
      
      IF(dim==1)THEN
         DO k=1,d3
            DO j=1,d2
               summ=0.0
               DO i=1,d1
                  summ=summ+ary3d(i,j,k)
               END DO
               avg2d(j,k)=summ
            END DO
         END DO
         avg2d=(1.0/real(d1))*avg2d
      ELSEIF(dim==2)THEN
         DO k=1,d3
            DO i=1,d1
               summ=0.0
               DO j=1,d2
                  summ=summ+ary3d(i,j,k)
               END DO
               avg2d(i,k)=summ
            END DO
         END DO
         avg2d=(1.0/real(d2))*avg2d
      ElSE
         DO j=1,d2
            DO i=1,d1
               summ=0.0
               DO k=1,d3
                  summ=summ+ary3d(i,j,k)
               END DO
               avg2d(i,j)=summ
            END DO
         END DO
         avg2d=(1.0/real(d3))*avg2d
      END IF
    END SUBROUTINE spatialAvg3D

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! subroutine finds the average over direction defined by (dim) of  !!
    !! 2D array                                                         !!
    !! INPUTS:                                                          !!
    !!     REAL*8:: ary2d                                               !!
    !!     INTEGER:: dim (1,or 2)                                       !!
    !! OUTPUT:                                                          !!
    !!     REAL*8:: avg1d                                               !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

    SUBROUTINE spatialAvg2D(ary2d,dim,avg1d)
      IMPLICIT NONE
      REAL*8,INTENT(IN),DIMENSION(:,:)::ary2d
      REAL*8,INTENT(OUT),DIMENSION(:)::avg1d
      INTEGER,INTENT(IN)::dim
      INTEGER::i,j,d1,d2
      REAL*8:: summ
      
      d1=size(ary2d,1)
      d2=size(ary2d,2)
      
      IF(dim==1)THEN
         DO j=1,d2
            summ=0.0
            DO i=1,d1
               summ=summ+ary2d(i,j)
            END DO
            avg1d(j)=summ
         END DO
         avg1d=(1.0/real(d1))*avg1d
      ELSE
         DO i=1,d1
            summ=0.0
            DO j=1,d2
               summ=summ+ary2d(i,j)
            END DO
            avg1d(i)=summ
         END DO
         avg1d=(1.0/real(d2))*avg1d
      END IF
    END SUBROUTINE spatialAvg2D

END MODULE arrayops
  
