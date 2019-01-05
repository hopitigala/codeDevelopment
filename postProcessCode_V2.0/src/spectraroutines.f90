module specrout
  implicit none
contains
  subroutine velflucforspec(flucary,phiary)
    implicit none
    real*8,intent(in),dimension(:,:,:) ::flucary
    real*8,intent(out),dimension(:,:,:)::phiary
    integer                            ::i,j,k,d1,d2,d3
    
    d1=size(flucary,1)
    d2=size(flucary,2)
    d3=size(flucary,3)
    
    do k=1,d3-1
       do j=1,d2
          do i=1,d1-1
             phiary(i,j,k)=flucary(i,j,k)
          end do
       end do
    end do
  end subroutine velflucforspec

  subroutine velfluc2daryforspec(ary3d,iplane,loc,ary2d)

    implicit none
    real*8,intent(in),dimension(:,:,:)::ary3d
    integer,intent(in)                ::iplane,loc
    real*8, intent(out),dimension(:,:):: ary2d
    integer                           ::i,j,k,d1,d2,d3
    d1=size(ary3d,1)
    d2=size(ary3d,2)
    d3=size(ary3d,3)

    if(iplane==1)then !xy plane
       do k=1,d3-1
          do j=1,d2
             ary2d(j,k)=ary3d(loc,j,k)
          end do
       end do
    elseif(iplane==2)then !yz plane
       do j=1,d2
          do i=1,d1-1
             ary2d(i,j)=ary3d(i,j,loc)
          end do
       end do
    end if
  end subroutine velfluc2daryforspec
  
  
  subroutine one_d_spectra(ary2d,iplane,spec2d)
    use compfft
    implicit none
    real*8,intent(in),dimension(:,:)::ary2d
    real*8,intent(out),dimension(:,:)::spec2d
    integer,intent(in)              ::iplane
    integer                          ::i,j,d1,d2
    real*8,allocatable,dimension(:)  ::ary1d
    complex(kind=8),allocatable,dimension(:)::ftary1d,enerspec
    d1=size(ary2d,1)
    d2=size(ary2d,2)
    
    if(iplane==1)then! xy plane x-cord in 2nd dimension
       do i=1,d1
          allocate(ary1d(d2))
          do j=1,d2
             ary1d(j)=ary2d(i,j)
          end do
          allocate(ftary1d(d2))
          call onedfft1darry(ary1d,d2,ftary1d)
          deallocate(ary1d)
          do j=1,d2
             spec2d(i,j)=abs(ftary1d(j)*conjg(ftary1d(j))*d2)
          end do
          deallocate(ftary1d)
       end do
    elseif(iplane==2)then !yz plane z-cord in 1st dimension
       do j=1,d2
          allocate(ary1d(d1))
          do i=1,d1
             ary1d(i)=ary2d(i,j)
          end do
          allocate(ftary1d(d1))
          call onedfft1darry(ary1d,d1,ftary1d)
          deallocate(ary1d)
          do i=1,d1
             spec2d(i,j)=abs(ftary1d(i)*conjg(ftary1d(i))*d1)
          end do
          deallocate(ftary1d)
       end do
    end if
  end subroutine one_d_spectra
end module specrout
