module specrout
  implicit none
contains
  subroutine asignphi(flucary,phiary)
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
  end subroutine asignphi
end module specrout
