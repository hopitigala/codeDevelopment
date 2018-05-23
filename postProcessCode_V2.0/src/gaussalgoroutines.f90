module gaussroutines
     implicit none
contains
    
!!!!!!!!    		Finds kernel matrix 				!!!!!!!!!
             !! input ::  'standard deviation' value 		  !!
             !! output::  'kernel matrix' used for the convolution!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                 
     subroutine findkernelmatr(sd,outarray)
          use channhdvariables
          use prelimcalvar
          use mainparameters
          implicit none
          real*8,intent(in)                   ::sd
          real*8,intent(out),dimension(:,:)   ::outarray
          real*8,allocatable,dimension(:,:)   ::kernmat
          real*8                              ::tot
          integer                             ::dimn,var1,var2,i,j
          dimn=size(outarray,1)
          var1=-(dimn-1)/2
          var2=-(dimn-1)/2
          allocate(kernmat(dimn,dimn))
          do i=1,dimn
               var2=-(dimn-1)/2
               do j=1,dimn
                    kernmat(j,i)=(1/(2*Pi*sd**2))*exp(-(var1**2+var2**2)/(2*sd**2))
                    var2=var2+1
               end do
               var1=var1+1
          end do 
          
          tot=sum(kernmat)
          outarray=kernmat/tot
          deallocate(kernmat)          
     end subroutine findkernelmatr

!!!!!!!!		finds standard deviation value                         !!!!!!!!
             !! input : ycoordinate value, z node and x node            !!
             !! output: standard deviation value at the particular point!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     subroutine findsd(ycoord,xnode,coeff_matx,val)
          use channhdvariables
          use prelimcalvar
          implicit none
          integer,intent(in)                ::xnode
          real*8,intent(out)                ::val
          real*8,intent(in)                 ::ycoord
          real*8,intent(in),dimension(:,:)  ::coeff_matx 

          val = coeff_matx(xnode,1)*(ycoord**5)+coeff_matx(xnode,2)*(ycoord**4)+coeff_matx(xnode,3)*(ycoord**3)&
                +coeff_matx(xnode,4)*(ycoord**2)+coeff_matx(xnode,5)*(ycoord)+coeff_matx(xnode,6)
     end subroutine findsd

!!!!!!!!                        Reads the standard deviation input file                     !!!!!!!!
                !!                    input : standard deviation data files            !!
                !! output: coefficient matrices in between the jets and behind the jets!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     subroutine read_sd(uucoeff_79,uucoeff_96,vvcoeff_79,vvcoeff_96,wwcoeff_79,wwcoeff_96)
          use channhdvariables
          use prelimcalvar
          implicit none
          integer                            ::i,j
          real*8,intent(out),dimension(:,:)  ::uucoeff_79,uucoeff_96,vvcoeff_79,vvcoeff_96,wwcoeff_79,wwcoeff_96 
          open(23,file='uu_sd_79_coeffiecients.dat')
          do i=1,n3
             Read(23,*) (uucoeff_79(i,j),j=1,6)
          end do
          close(23)
          !write(*,*)'I am here and priniting sd_bwn',coeff_79(3,1:6)
          open(24,file='uu_sd_96_coeffiecients.dat')
          do i=1,n3
             Read(24,*) (uucoeff_96(i,j),j=1,6)
          end do
          close(24)
          !write(*,*)'I am here and printing sd_alo',coeff_96
          open(25,file='vv_sd_79_coeffiecients.dat')
          do i=1,n3
             Read(25,*) (vvcoeff_79(i,j),j=1,6)
          end do
          close(25)
          open(26,file='vv_sd_79_coeffiecients.dat')
          do i=1,n3
             Read(26,*) (vvcoeff_96(i,j),j=1,6)
          end do
          close(26)
          open(27,file='ww_sd_79_coeffiecients.dat')
          do i=1,n3
             Read(27,*) (wwcoeff_79(i,j),j=1,6)
          end do
          close(27)
          open(28,file='ww_sd_96_coeffiecients.dat')
          do i=1,n3
             Read(28,*) (wwcoeff_96(i,j),j=1,6)
          end do
          close(28)
          
     end subroutine read_sd

!!!!!!!!                      Finds the whole matrix                                         !!!!!!!!
              !!    Input : array from each core                                   !!
              !!    Output: Entire array involving the entire domain               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     subroutine whole_maker(part_matx,tag,whole_matx)
          use mpi
          use mpivariables
          implicit none
          integer,intent(in)                  :: tag
          real*8,intent(in),dimension(:,:,:)  :: part_matx
          real*8,intent(out),dimension(:,:,:) :: whole_matx
          integer                             :: i,j,k,d1,d2,d3,p,jj
          d1=size(part_matx,1)
          d2=size(part_matx,2)
          d3=size(part_matx,3)

          if (mynode/=0)then
             call MPI_SEND(part_matx,d1*d2*d3,MPI_REAL8,0,tag,MPI_COMM_WORLD,ierr)
          else 
             do k=1,d3
                  do j=1,d2
                       jj=j
                       do i=1,d1
                            whole_matx(i,jj,k)=part_matx(i,j,k)
                       end do
                  end do
             end do
             do p=1,(numprocs-1)
                  call MPI_RECV(part_matx,d1*d2*d3,MPI_REAL8,p,tag,MPI_COMM_WORLD,status,ierr)
                  do k=1,d3
                       do j=1,d2
                            jj=j+p*(d2-1)
                            do i=1,d1
                                 whole_matx(i,jj,k)=part_matx(i,j,k)
                            end do
                       end do
                  end do
              end do
             !write(*,*)'the size of the whole matrix is',shape(whole_matx)
          end if
     end subroutine whole_maker
    
!!!!!!!!              finds convolution for a given point                                    !!!!!!!!
              !!    Input : kernel matrix, 3D uprime matrix at a given timestep    !!
              !!    Output: filtered point                                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     subroutine convolution(xn,yn,zn,kmatx,datamatx,point)
          use prelimcalvar
          implicit none
          Real*8,intent(in),dimension(:,:)   :: kmatx
          Real*8,intent(in),dimension(:,:,:) :: datamatx
          Real*8,intent(out)                 :: point
          Real*8,allocatable,dimension(:,:)  :: inter_2Draw,rawdata_2D,result_matx
          integer,intent(in)                 :: xn,yn,zn
          integer                            :: i,j,ksize
          ksize = (size(kmatx,1)-1)/2
          !write(*,*) 'size of the ksize',ksize
          allocate(inter_2Draw(size(datamatx(:,:,1),1),size(datamatx(:,:,1),2)))
          inter_2Draw=datamatx(:,:,xn)
          !write(*,*) 'size of raw2D_data before padding is',shape(inter_2Draw)
          allocate(rawdata_2D(size(inter_2Draw,1)+2*ksize,size(inter_2Draw,2)+2*ksize))
          !write(*,*) 'size of raw2D_data after padding is',shape(rawdata_2D) 
          rawdata_2D(1:ksize,:)=0
          rawdata_2D(size(rawdata_2D,1)-ksize+1:size(rawdata_2D,1),:)=0
          rawdata_2D(ksize+1:size(rawdata_2D,1)-ksize,1:ksize)=0
          rawdata_2D(ksize+1:size(rawdata_2D,1)-ksize,size(rawdata_2D,2)-ksize+1:size(rawdata_2D,2))=0
          rawdata_2D(ksize+1:size(rawdata_2D,1)-ksize,ksize+1:size(rawdata_2D,2)-ksize)=inter_2Draw
          !write(*,*) 'size after padding',shape(rawdata_2D)
          deallocate(inter_2Draw)
          allocate(result_matx(size(kmatx,1),size(kmatx,2)))

          do j=1,size(kmatx,2)
               do i=1,size(kmatx,1)
                    result_matx(i,j) = kmatx(i,j)*rawdata_2D(zn+i-1,yn+j-1)
               end do
          end do
          point = sum(result_matx)
          deallocate(rawdata_2D,result_matx)    
     end subroutine convolution
   
!!!!!!!!     makes all the instantaneous velocities in one 4D array                  !!!!!!!!
       !! Input : umean, vmean, wmean, filtered uprime and unfiltered vprime, wprime !!     
       !! Output: instantaneous velocities                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine velmaker(umean,vmean,wmean,fltd_uprm,fltd_vprm,fltd_wprm,vel_tot)
          use channhdvariables
          use prelimcalvar
          implicit none
          Real*8,intent(in),dimension(:,:,:)    :: umean,vmean,wmean,fltd_uprm,fltd_vprm,fltd_wprm
          Real*8,intent(out),dimension(:,:,:,:) :: vel_tot
          Real*8,dimension(n1,n2do+1,n3)        :: u,v,w
          u = umean+fltd_uprm
          v = vmean+fltd_vprm
          w = wmean+fltd_wprm
          vel_tot(:,:,:,1) = w
          vel_tot(:,:,:,2) = v
          vel_tot(:,:,:,3) = u
          write(*,*)'The shape of the vel_tot is', shape(vel_tot)
     end subroutine velmaker
          

end module gaussroutines

                    
