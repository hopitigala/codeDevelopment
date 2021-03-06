program gaussalgo
   USE channhdvariables
   USE mainparameters
   USE writedata
   USE arrayops
   USE readdata
   USE mpivariables
   USE prelimcalvar
   USE interpoldata
   USE fluctuation
   USE mpi
   USE gaussroutines
   USE compfft
   USE vorticity
   USE vortexiden

IMPLICIT NONE
REAL *8,ALLOCATABLE,DIMENSION(:)::xp,yp,ypdo,ysdo,ys,zp
REAL *8,ALLOCATABLE,DIMENSION(:,:,:)::utavg,vtavg,wtavg,ttavg,ptavg,pp,lambda2
REAL *8,ALLOCATABLE,DIMENSION(:,:,:)::uprime,vprime,wprime,pprime,tprime,uprm_whole,vprm_whole,wprm_whole
REAL *8,ALLOCATABLE,DIMENSION(:,:,:)::gauss_filtduprm,low_fltd_ifuprm,gauss_filtdvprm,low_fltd_ifvprm,gauss_filtdwprm,low_fltd_ifwprm
REAL *8,ALLOCATABLE,DIMENSION(:,:,:,:)::up,tp,dudx,dudy,dudz
REAL *8,ALLOCATABLE,DIMENSION(:,:)::kernelmatx
REAL *8,ALLOCATABLE,DIMENSION(:,:)::sdfun_uucoeff_bwn,sdfun_uucoeff_alo,sdfun_vvcoeff_bwn,sdfun_vvcoeff_alo,sdfun_wwcoeff_bwn,sdfun_wwcoeff_alo
COMPLEX(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)::fft_uprm,low_fltd_fuprm,fft_vprm,low_fltd_fvprm,fft_wprm,low_fltd_fwprm


INTEGER::dt,timestcount,numtimesteps,i,j,k,itime,xcount,ycount,zcount
INTEGER::treq,k_ucrtcl,k_vcrtcl,k_wcrtcl,ksize
INTEGER::icorlavg,icen
REAL*8::inv,sttime
REAL*8::sd,filtd_point
CHARACTER*5::string
CHARACTER*2::prosnum
CHARACTER*5::time

CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,mynode,ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)


CALL readChannelData()
CALL readPostData(sttime,numtimesteps,dt,icorlavg,icen,k_ucrtcl,k_vcrtcl,k_wcrtcl)
CALL prelimCal()
!CALL wavenumcritc(k_crtcl)
allocate(sdfun_uucoeff_bwn(n3,6),sdfun_uucoeff_alo(n3,6),sdfun_vvcoeff_bwn(n3,6),sdfun_vvcoeff_alo(n3,6),sdfun_wwcoeff_bwn(n3,6),sdfun_wwcoeff_alo(n3,6))
CALL read_sd(sdfun_uucoeff_bwn,sdfun_uucoeff_alo,sdfun_vvcoeff_bwn,sdfun_vvcoeff_alo,sdfun_wwcoeff_bwn,sdfun_wwcoeff_alo)
ALLOCATE(xp(n3),yp(0:n2),zp(n1),ys(0:n2),ypdo(0:n2do+1),ysdo(0:n2do+1))
CALL readCoordData(xp,yp,zp)
CALL intpolCordi(mynode,numprocs,yp,ys,ypdo,ysdo)

ALLOCATE(wtavg(n1,n2do+1,n3),vtavg(n1,n2do+1,n3),&
     utavg(n1,n2do+1,n3),ttavg(n1,n2do+1,n3),ptavg(n1,n2do+1,n3))

!!!!! Read mean field

write(prosnum,'(i2.2)')mynode

call read3Darray('../tmean'//trim(prosnum),'unformatted',ttavg)
call read3Darray('../wmean'//trim(prosnum),'unformatted',wtavg)
call read3Darray('../vmean'//trim(prosnum),'unformatted',vtavg)
call read3Darray('../umean'//trim(prosnum),'unformatted',utavg)
call read3Darray('../pmean'//trim(prosnum),'unformatted',ptavg)

timestcount = 0
!!!!! time loop starts to read data from the DNS field data files

stime=MPI_WTIME()

     DO itime=int(sttime),(int(sttime)+(numtimesteps*Dt)),Dt
          timestcount=timestcount+1
          Allocate(up(n1,n2do+1,n3,3),tp(n1,n2do+1,n3,1),pp(n1,n2do+1,n3))
          
          call readTempFieldData(mynode,itime,ypdo,ysdo,up,pp,tp)
          allocate(tprime(n1,n2do+1,n3),wprime(n1,n2do+1,n3),&
                  vprime(n1,n2do+1,n3),uprime(n1,n2do+1,n3),pprime(n1,n2do+1,n3))


          call fluct3Dmean(tp(:,:,:,1),ttavg,tprime)
          call fluct3Dmean(up(:,:,:,1),wtavg,wprime)
          call fluct3Dmean(up(:,:,:,2),vtavg,vprime)
          call fluct3Dmean(up(:,:,:,3),utavg,uprime)

          DEALLOCATE(up,tp,pp)
     allocate(uprm_whole(n1,n2,n3),vprm_whole(n1,n2,n3),wprm_whole(n1,n2,n3))
     call whole_maker(uprime,1,uprm_whole)
     call whole_maker(vprime,1,vprm_whole)
     call whole_maker(wprime,1,wprm_whole)
     !write(*,*),'The uprime whole size is',shape(uprm_whole)

!! u prime 

     allocate(gauss_filtduprm(n1,n2do+1,n3))
     Do xcount=1,n3
         Do ycount=1,n2do+1
              Do zcount=1,n1
                   
                   if(zcount.ge.1.and.zcount.le.13 .or. zcount.ge.37.and.zcount.le.53 .or. zcount.ge.71.and.zcount.le.87 .or.&
                      zcount.ge.105.and.zcount.le.121.or.zcount.ge.139.and.zcount.le.155.or.zcount.ge.173.and.zcount.le.193)then
                        call findsd(yp(ycount),xcount,sdfun_uucoeff_bwn,sd)
                        !write(*,*) 'I am in if loop for sd and the value is and z value is',sd,zcount
                   else
                        call findsd(yp(ycount),xcount,sdfun_uucoeff_alo,sd)
                   end if
                   
                   ksize = ceiling(abs(3*sd))
                   allocate(kernelmatx(2*ksize+1,2*ksize+1))
                   call findkernelmatr(sd,kernelmatx)
                   call convolution(xcount,ycount,zcount,kernelmatx,uprm_whole,filtd_point)
                   gauss_filtduprm(zcount,ycount,xcount)=filtd_point
                   deallocate(kernelmatx)
              end Do
         end Do
     end Do
     deallocate(uprm_whole)
     deallocate(uprime,tprime,pprime)
     allocate(fft_uprm(n1,n2do+1,n3))
     call onedfft3dary(gauss_filtduprm,n3,3,fft_uprm)
     deallocate(gauss_filtduprm)
     allocate(low_fltd_fuprm(n1,n2do+1,n3))
     
     fft_uprm(:,:,(k_ucrtcl+1):(n3-k_ucrtcl+1))=0

     low_fltd_fuprm = fft_uprm
  
     deallocate(fft_uprm)
     allocate(low_fltd_ifuprm(n1,n2do+1,n3))
     call oneDInvFFT3Dary(low_fltd_fuprm,n3,3,low_fltd_ifuprm)
     deallocate(low_fltd_fuprm)
     write(*,*) 'The uprm is done by',prosnum,'for timestep',itime 

!! vprime
     call whole_maker(vprime,1,vprm_whole)
     allocate(gauss_filtdvprm(n1,n2do+1,n3))
     Do xcount=1,n3
         Do ycount=1,n2do+1
              Do zcount=1,n1

                   if(zcount.ge.1.and.zcount.le.13 .or. zcount.ge.37.and.zcount.le.53 .or. zcount.ge.71.and.zcount.le.87 .or.&
                      zcount.ge.105.and.zcount.le.121.or.zcount.ge.139.and.zcount.le.155.or.zcount.ge.173.and.zcount.le.193)then
                        call findsd(yp(ycount),xcount,sdfun_vvcoeff_bwn,sd)
                        !write(*,*) 'I am in if loop for sd and the value is and z value is',sd,zcount
                   else
                        call findsd(yp(ycount),xcount,sdfun_vvcoeff_alo,sd)
                   end if

                   ksize = ceiling(abs(3*sd))
                   allocate(kernelmatx(2*ksize+1,2*ksize+1))
                   call findkernelmatr(sd,kernelmatx)
                   call convolution(xcount,ycount,zcount,kernelmatx,vprm_whole,filtd_point)
                   gauss_filtdvprm(zcount,ycount,xcount)=filtd_point
                   deallocate(kernelmatx)
              end Do
         end Do
     end Do
     deallocate(vprm_whole)
     deallocate(vprime)
     allocate(fft_vprm(n1,n2do+1,n3))
     call onedfft3dary(gauss_filtdvprm,n3,3,fft_vprm)
     deallocate(gauss_filtdvprm)
     allocate(low_fltd_fvprm(n1,n2do+1,n3))

     fft_vprm(:,:,(k_vcrtcl+1):(n3-k_vcrtcl+1))=0
     low_fltd_fvprm = fft_vprm
     deallocate(fft_vprm)
     allocate(low_fltd_ifvprm(n1,n2do+1,n3))
     call oneDInvFFT3Dary(low_fltd_fvprm,n3,3,low_fltd_ifvprm)
     deallocate(low_fltd_fvprm)
    write(*,*) 'The vprm is done by',prosnum,'for timestep',itime

!! wprime
     call whole_maker(wprime,1,wprm_whole)
     allocate(gauss_filtdwprm(n1,n2do+1,n3))
     Do xcount=1,n3
         Do ycount=1,n2do+1
              Do zcount=1,n1

                   if(zcount.ge.1.and.zcount.le.13 .or.zcount.ge.37.and.zcount.le.53 .or. zcount.ge.71.and.zcount.le.87 .or.&
                      zcount.ge.105.and.zcount.le.121.or.zcount.ge.139.and.zcount.le.155.or.zcount.ge.173.and.zcount.le.193)then
                        call findsd(yp(ycount),xcount,sdfun_wwcoeff_bwn,sd)
                        !write(*,*) 'I am in if loop for sd and the value is and
                        !z value is',sd,zcount
                   else
                        call findsd(yp(ycount),xcount,sdfun_wwcoeff_alo,sd)
                   end if

                   ksize = ceiling(abs(3*sd))
                   allocate(kernelmatx(2*ksize+1,2*ksize+1))
                   call findkernelmatr(sd,kernelmatx)
                   call convolution(xcount,ycount,zcount,kernelmatx,wprm_whole,filtd_point)
                   gauss_filtdwprm(zcount,ycount,xcount)=filtd_point
                   deallocate(kernelmatx)
              end Do
         end Do
     end Do
     deallocate(wprm_whole)
     deallocate(wprime)
     allocate(fft_wprm(n1,n2do+1,n3))
     call onedfft3dary(gauss_filtdwprm,n3,3,fft_wprm)
     deallocate(gauss_filtdwprm)
     allocate(low_fltd_fwprm(n1,n2do+1,n3))

     fft_wprm(:,:,(k_wcrtcl+1):(n3-k_wcrtcl+1))=0
     low_fltd_fwprm = fft_wprm
     deallocate(fft_wprm)
     allocate(low_fltd_ifwprm(n1,n2do+1,n3))
     call oneDInvFFT3Dary(low_fltd_fwprm,n3,3,low_fltd_ifwprm)
     deallocate(low_fltd_fwprm)
    write(*,*) 'The wprm is done by',prosnum,'for timestep',itime

 
     allocate(up(n1,n2do+1,n3,3))
     call velmaker(utavg,vtavg,wtavg,low_fltd_ifuprm,low_fltd_ifvprm,low_fltd_ifwprm,up)
     deallocate(low_fltd_ifuprm,low_fltd_ifvprm,low_fltd_ifwprm)
     allocate(dudx(n1,n2do+1,n3,3),dudy(n1,n2do+1,n3,3),dudz(n1,n2do+1,n3,3))
     call velograd(up,xp,ypdo,zp,dudx,dudy,dudz)
     deallocate(up)
     allocate(lambda2(n1,n2do+1,n3))
     call complambda2(dudx,dudy,dudz,lambda2)
     deallocate(dudx,dudy,dudz)
       
     ! prints filtered field at 40 multiples of time
     treq = stTime
     write(time,'(i5.5)')itime
     write(string,'(i5.5)')treq
     
     if(time.EQ.string) then
          !call sendrecv3dwrite(low_fltd_ifuprm,1,'filtd_uprm_'//trim(time)//'.dat')
          call sendrecv3dwrite(lambda2,1,'lambda2_'//trim(time)//'.dat')
          treq=treq+20
     end if
     !deallocate(low_fltd_ifuprm)
     deallocate(lambda2)
write(*,*) 'lambda2 found for',itime    
          
     END DO
deallocate(sdfun_uucoeff_bwn,sdfun_uucoeff_alo,sdfun_vvcoeff_bwn,sdfun_vvcoeff_alo,sdfun_wwcoeff_bwn,sdfun_wwcoeff_alo)
deallocate(xp,yp,zp,ys,ysdo,ypdo)
deallocate(wtavg,vtavg,utavg,ttavg,ptavg)
etime=MPI_WTIME()
write(*,*)'The time taken by the gauss filtered process is',etime-stime

         

end program gaussalgo

       
