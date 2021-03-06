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
   USE derivatives


IMPLICIT NONE
REAL *8,ALLOCATABLE,DIMENSION(:)::xp,yp,ypdo,ysdo,ys
REAL *8,ALLOCATABLE,DIMENSION(:,:,:)::utavg,vtavg,wtavg,ttavg,ptavg,pp
REAL *8,ALLOCATABLE,DIMENSION(:,:,:)::uprime,vprime,wprime,pprime,tprime
REAL *8,ALLOCATABLE,DIMENSION(:,:,:)::gauss_filtduprm
REAL *8,ALLOCATABLE,DIMENSION(:,:,:,:)::up,tp
REAL *8,ALLOCATABLE,DIMENSION(:,:)::kernelmatx
REAL *8,ALLOCATABLE,DIMENSION(:,:)::sdfun_coeff_bwn,sdfun_coeff_alo
REAL *8,ALLOCATABLE,DIMENSION(:,:)::fft_uprm,low_fltd_fuprm,low_fltd_ifuprm


INTEGER::dt,timestcount,numtimesteps,i,j,k,itime,xcount,ycount,zcount
INTEGER::treq,kcrtcl
INTEGER::icorlavg,icen
REAL*8::inv,sttime
REAL*8::sd,filtd_point
REAL*8::uprm_whole
CHARACTER*5::string
CHARACTER*2::prosnum
CHARACTER*5::time

CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,mynode,ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numproces,ierr)


CALL readChannelData()
CALL readPostData(sttime,numtimesteps,dt,icorlavg,icen)
CALL prelimCal()
CALL wavenumcritc(k_crtcl)
allocate(sdfun_coeff_bwn(n3,4),sdfun_coeff_alo(n3,4))
CALL read_sd(sdfun_coeff_bwn,sdfun_coeff_alo)

ALLOCATE(xp(n3),yp(0:n2),zp(n1),ys(0:n2),ypdo(0:n2do+1),usdo(0:n2do+1))
CALL readCoordData(xp,yp,zp)
CALL intpolCordi(mynode,numprocs,yp,ys,ypdo,ysdo)

ALLOCATE(wtavg(n1,n2do+1,n3),vtavg(n1,n2do+1,n3),&
     utavg(n1,n2do+1,n3),ttavg(n1,n2do+1,n3),ptavg(n1,n2do+1,n3))

!!!!! Read mean field

write(prosnum,'(i2.2)')mynode

call read3Darray('tmean'//trim(prosnum),'unformatted',ttavg)
call read3Darray('wmean'//trim(prosnum),'unformatted',wtavg)
call read3Darray('vmean'//trim(prosnum),'unformatted',vtavg)
call read3Darray('umean'//trim(prosnum),'unformatted',utavg)
call read3Darray('pmean'//trim(prosnum),'unformatted',ptavg)

timestcount = 0
!!!!! time loop starts to read data from the DNS field data files

stime=MPI_WTIME()

     DO itime=int(sttime),(int(sttime)+(numtimesteps*Dt)),Dt
          timestepcount=timestepcount+1
          Allocate(up(n1,n2do+1,n3,3),tp(n1,n2do+1,n3,1),pp(n1,n2do+1,n3))
          
          call readTempFieldData(mynode,itimes,ypdo,ysdo,up,pp,tp)
          allocate(tprime(n1,n2do+1,n3),wprime(n1,n2do+1,n3),&
                  vprime(n1,n2do+1,n3),uprime(n1,n2do+1,n3),pprime(n1,n2do+1,n3))


          call fluct3Dmean(tp(:,:,:,1),ttavg,tprime)
          call fluct3Dmean(up(:,:,:,1),wtavg,wprime)
          call fluct3Dmean(up(:,:,:,2),vtavg,vprime)
          call fluct3Dmean(up(:,:,:,3),utavg,uprime)

          DEALLOCATE(up,tp,pp)
     allocate(uprm_whole(n1,n2,n3))
     call whole_maker(uprime,1,uprm_whole)
     allocate(gauss_filtduprm(n1,n2do+1,n3))
     Do xcount=1,n3
         Do ycount=1,n2do+1
              Do zcount=1,n1
                   
                   if(zcount.ge.1.and.zcount.le.13 .or. zcount.ge.37.and.zcount.le.53 .or. zcount.ge.71.and.zcount.le.87 .or.&
                      zcount.ge.105.and.zcount.le.121.or.zcount.ge.139.and.zcount.le.155.or.zcount.ge.173.and.zcount.le.193)
                        call findsd(yp(ycount),xcount,sdfun_coeff_bwn,sd)
                   else
                        call findsd(yp(ycount),xcount,sdfun_coeff_alo,sd)
                   end if
                   
                   ksize = ceiling(3*sd)
                   allocate(kernelmatx(2*ksize+1,2*ksize+1))
                   call findkernelmatr(sd,kernelmatx)
                   call convolution(xcount,ycount,zcount,kernelmatx,uprm_whole,filtd_point)
                   gauss_filtduprm(zcount,ycount,xcount)=filtd_point
              end Do
         end Do
     end Do
     deallocate(uprm_whole)
     allocate(fft_uprm)
     call onedfft3dary(gauss_filtuprm,n1,n2do+1,n3,3,fft_uprm)
     deallocate(gauss_filtduprm)
     allocate(low_filtd_uprm(n1,n2do+1,n3))
     
     fft_uprm(:,:,kcrtcl+1:n3-kcrtcl+1)=0
     low_fltd_fuprm = fft_uprm
     deallocate(fft_uprm)
     allocate(low_fltd_ifuprm(n1,n2do+1,n3))
     call onedinvfft3dary(low_fltd_fuprm,low_fltd_ifuprm)
     deallocate(low_fltd_fuprm)  
     ! prints filtered field at 40 multiples of time
     treq = stTime
     write(time,'(i5.5)')itime
     write(string,'(i5.5)')treq
     
     if(time.eq.string) then
          call sendrecv3dwrite(low_fltd_ifuprm,1,'filtd_uprm'//trim(time)//'.dat')
          treq=treq+40
     end if

     deallocate(low_fltd_ifuprm)    
          
     END DO
etime=MPI_WTIME()
write(*,*)'The time taken by the gauss filtered process is',etime-stime

         

end program gaussalgo

       
