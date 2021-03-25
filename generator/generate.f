      program generate

!     gfortran -funderscoring -O3 -Wall -fmessage-length=0 -ffree-line-length-2000 -ffree-form -finit-local-zero -falign-commons generate.f -o generate
      
      implicit none
      save
!     simulation data
      double precision :: maxEnergy, minEnergy, Z, MT_FLAG, C_FLAG
      integer :: v, ios
      character (len=1000) :: line

      write(*,*) '-----------------------------------------------------------------------------------------------------------'
      write(*,*) '-                                                  PolEM                                                  -'
      write(*,*) '-                           Polarized Electron Mott Scattering Model for Geant4                           -'
      write(*,*) '-                                         cross-section generator                                         -'
      write(*,*) '-----------------------------------------------------------------------------------------------------------'
      write(*,*) ''
      
!     read simulation data
      open (unit=55, file = 'generate.inp', form = 'formatted', status = 'old')
      v=0
      ios = 0
      do while ( ios .eq. 0 )
        read (55, fmt = '(a500)', iostat = ios) line
        if ((ios .eq. 0) .and. (line(1:1) .ne. '!')) then
          if (v .eq. 0) then
            read (line, *) maxEnergy
          elseif (v .eq. 1) then
            read (line, *) minEnergy
          elseif (v .eq. 2) then
            read (line, *) Z
          elseif (v .eq. 3) then
            read (line, *) MT_FLAG
          elseif (v .eq. 4) then
            read (line, *) C_FLAG
          end if
          v=v+1
        end if
      end do
      close(55)
      
!     calculate cross sections
      call calc_elastic(maxEnergy,minEnergy,nint(Z),nint(MT_FLAG),nint(C_FLAG))

      end program generate



!     constants
      MODULE constants
      implicit none
      save
      double precision, parameter :: PI=4.d0*datan(1.d0)
      double precision, parameter :: PIx2=8.d0*datan(1.d0)
      double precision, parameter :: d2r=PI/180.0d0
      double precision, parameter :: nm_to_cm=10.0d0**(-7)
      double precision, parameter :: barn=10.0d0**(-24) ! in cm2
      double precision, parameter :: a0_square=2.8002852d0*10.d0**(-17) ! in cm2
      double precision, parameter :: h_bar_c=197.327d0*1000.d0*10.d0**(-13) ! in cm2
      double precision, parameter :: Ry=13.6057d0 ! in eV
      double precision, parameter :: alpha=7.2973525698D-3
      double precision, parameter :: m_e=510.998910d0 ! in keV
      double precision, parameter :: c=29979245800.00019d0
      double precision, parameter :: au=1.660538d0*10.d0**(-24) ! [g]
      double precision, parameter :: r_e=2.8179403267d0*10.d0**(-13.) ! in cm
      double precision, parameter :: val_small=tiny(1.d0)
      double precision, parameter :: val_huge=huge(1.d0)
      integer, parameter :: RegInt_K = selected_int_kind(8)
      integer, parameter :: RegReal_K = selected_real_kind(12)
      integer, parameter :: QR_K = selected_real_kind (32)
      real (kind=QR_K), parameter :: QR_one=real(1.d0,QR_K), QR_zero=real(0.d0,QR_K)
      end module constants

      

!     subroutine calc_elastic
!     Michal Dragowski, Gunter Weber and Marta Wlodarczyk
!     M. Dragowski et al., Nucl. Instr. Meth. B 488, 37 (2021)
      SUBROUTINE calc_elastic(maxEnergy,minEnergy,Z,MT_FLAG,C_FLAG)
      use constants
      implicit none
      double precision, intent(in) :: maxEnergy, minEnergy
      integer, intent(in) :: Z, MT_FLAG, C_FLAG
!     cross sections
      integer, parameter :: nAngularSteps=606
      integer, parameter :: nEnergySteps=500
      double precision :: angle(nAngularSteps), cs(nAngularSteps,nEnergySteps,4), total(nEnergySteps,2), differential(nAngularSteps,nEnergySteps)
      double precision :: en_step
!     correct errors
      double precision :: dTmp, sTmp, tTmp, uTmp, dNew, sNew, tNew, uNew
!     others
      integer :: v, k, I
      double precision :: energy, temp, fact, x
      double precision :: find_value
      character(len=4) :: filename
      
!     set all cross sections to zero
      do I=1, nEnergySteps
        do v=1, nAngularSteps
          do k=1, 4
            cs(v,I,k)=0.d0
          end do
        end do
      end do
      
!     calculate elastic scattering cross sections and polarization correlations for various electron energies
      en_step=dlog(maxEnergy+1.d0-minEnergy)/dble(nEnergySteps-1)
      do I=1, nEnergySteps
        if (I .eq. 1) then
          energy=minEnergy
        end if
        if (I .eq. nEnergySteps) then
          energy=maxEnergy
        end if
        if (I .gt. 1 .and. I .lt. nEnergySteps) then
          energy=minEnergy+dexp(en_step*dble(I-1))-1.d0
        end if
!       fill in energy table
        total(I,1)=energy
!       call ELSEPA
        call ELSCATA((energy*1000.d0),Z,angle,cs(:,I,1),cs(:,I,2),cs(:,I,3),cs(:,I,4),MT_FLAG,0.d0)
!       check S T U functions normalization
        do v=1, nAngularSteps
          do k=2, 4
            cs(v,I,k)=max(min(cs(v,I,k),1.d0),-1.d0)
          end do
          temp=dsqrt(cs(v,I,2)**2+cs(v,I,3)**2+cs(v,I,4)**2)
          do k=2, 4
            if (temp .lt. val_small) then
              cs(v,I,k)=0.d0
            else
              cs(v,I,k)=cs(v,I,k)/temp
            end if
          end do
        end do
!       periodic printing
        if (I .eq. nEnergySteps) then
          write(*,'(a1)') '*'
        else
          write(*,'(a1)',advance='no') '*'
        end if
      end do

!     integration in spherical coordinates along theta axis (sin(theta) factor)
      do I=1, nEnergySteps
        differential(1,I)=0.d0
        do v=2, nAngularSteps
          x=(angle(v)+angle(v-1))/2.0d0
          fact=dsin(x*d2r)
          differential(v,I)=differential(v-1,I)+dabs((angle(v)-angle(v-1)))*find_value(angle,cs(:,I,1),x,nAngularSteps,0,2)*fact*PIx2*d2r
        end do
!       elastic scattering total cross section
        total(I,2)=differential(nAngularSteps,I)
      end do

!     normalize differential cross sections
      do I=1, nEnergySteps
        do v=2, nAngularSteps
          differential(v,I)=differential(v,I)/differential(nAngularSteps,I)
        end do
      end do

!     correct errors from ELSEPA
      if (C_FLAG .eq. 1) then
        do v=1, nAngularSteps
          do I=3, nEnergySteps
            if (total(I,1)>10.d0 .and. dabs(cs(v,I,2)-cs(v,I-1,2)) .gt. 2*dabs(cs(v,I-1,2)-cs(v,I-2,2))) then
              sTmp = cs(v,I-1,2) + (total(I,1)-total(I-1,1)) * (cs(v,I+1,2)-cs(v,I-1,2)) / (total(I+1,1)-total(I-1,1))
              tTmp = cs(v,I-1,3) + (total(I,1)-total(I-1,1)) * (cs(v,I+1,3)-cs(v,I-1,3)) / (total(I+1,1)-total(I-1,1))
              uTmp = cs(v,I-1,4) + (total(I,1)-total(I-1,1)) * (cs(v,I+1,4)-cs(v,I-1,4)) / (total(I+1,1)-total(I-1,1))
              dTmp = differential(v,I-1) + (total(I,1)-total(I-1,1)) * (differential(v,I+1)-differential(v,I-1)) / (total(I+1,1)-total(I-1,1))
              do k=2, 5
                if (dabs(sTmp-cs(v,I-1,2)) .lt. dabs(cs(v,I,2)-cs(v,I-1,2)) .or. I .gt. nEnergySteps-1-k) then
                  goto 1212
                end if
                sNew = cs(v,I-1,2) + (total(I,1)-total(I-1,1)) * (cs(v,I+1+k,2)-cs(v,I-1,2)) / (total(I+1+k,1)-total(I-1,1))
                tNew = cs(v,I-1,3) + (total(I,1)-total(I-1,1)) * (cs(v,I+1+k,3)-cs(v,I-1,3)) / (total(I+1+k,1)-total(I-1,1))
                uNew = cs(v,I-1,4) + (total(I,1)-total(I-1,1)) * (cs(v,I+1+k,4)-cs(v,I-1,4)) / (total(I+1+k,1)-total(I-1,1))
                dNew = differential(v,I-1) + (total(I,1)-total(I-1,1)) * (differential(v,I+1+k)-differential(v,I-1)) / (total(I+1+k,1)-total(I-1,1))
                if (dabs(sNew-cs(v,I-1,2)) .lt. dabs(sTmp-cs(v,I-1,2))) then
                  sTmp=sNew
                  tTmp=tNew
                  uTmp=uNew
                  dTmp=dNew
                end if
              end do
1212          if (dabs(sTmp-cs(v,I-1,2)) .lt. dabs(cs(v,I,2)-cs(v,I-1,2))) then
                cs(v,I,2)=sTmp
                cs(v,I,3)=tTmp
                cs(v,I,4)=uTmp
                differential(v,I)=dTmp
              else
                sTmp = cs(v,I-1,2) + (total(I,1)-total(I-1,1)) * (cs(v,I-1,2)-cs(v,I-2,2)) / (total(I-1,1)-total(I-2,1))
                tTmp = cs(v,I-1,3) + (total(I,1)-total(I-1,1)) * (cs(v,I-1,3)-cs(v,I-2,3)) / (total(I-1,1)-total(I-2,1))
                uTmp = cs(v,I-1,4) + (total(I,1)-total(I-1,1)) * (cs(v,I-1,4)-cs(v,I-2,4)) / (total(I-1,1)-total(I-2,1))
                dTmp = differential(v,I-1) + (total(I,1)-total(I-1,1)) * (differential(v,I-1)-differential(v,I-2)) / (total(I-1,1)-total(I-2,1))
                if (dabs(sTmp-cs(v,I-1,2)) .lt. dabs(cs(v,I,2)-cs(v,I-1,2))) then
                  cs(v,I,2)=sTmp
                  cs(v,I,3)=tTmp
                  cs(v,I,4)=uTmp
                  differential(v,I)=dTmp
                end if
              end if
            end if
          end do
        end do
        do I=1, nEnergySteps
          do v=2, nAngularSteps-1
            if (differential(v,I) .lt. differential(v-1,I)) then
              differential(v,I) = (differential(v-1,I)+differential(v+1,I))/2.d0
              if (differential(v,I) .lt. differential(v-1,I)) then
                write(*,*) 'WARNING: Differential cross section table inconsistent!'
              end if
              cs(v,I,2) = (cs(v-1,I,2)+cs(v+1,I,2))/2.d0
              cs(v,I,3) = (cs(v-1,I,3)+cs(v+1,I,3))/2.d0
              cs(v,I,4) = (cs(v-1,I,4)+cs(v+1,I,4))/2.d0
            end if
          end do
        end do
      end if

!     write total electron-atom cross section as a function of electron kinetic energy
      write(filename, '(I0)') Z
      open (unit=54, file='data/t'//trim(filename)//'.dat', form='unformatted', access='direct', status = 'replace', recl=nEnergySteps*8)
      write(54,rec=1) total(:,1)
      write(54,rec=2) total(:,2)
      close(54)
!     write angular differential electron-atom cross section and S, T, U values for each electron kinetic energy
      open (unit=54, file='data/d'//trim(filename)//'.dat', form='unformatted', access='direct', status = 'replace', recl=nAngularSteps*8)
      write(54,rec=1) angle
      do I=1, nEnergySteps
        write(54,rec=1+I) differential(:,I)
      end do
      do I=1, nEnergySteps
        write(54,rec=1+nEnergySteps+I) cs(:,I,2)
      end do
      do I=1, nEnergySteps
        write(54,rec=1+2*nEnergySteps+I) cs(:,I,3)
      end do
      do I=1, nEnergySteps
        write(54,rec=1+3*nEnergySteps+I) cs(:,I,4)
      end do
      close(54)
      
      RETURN
      END SUBROUTINE calc_elastic



      double precision function find_value(x_array,y_array,x,n_max,n_start,flag)
      use constants, only: val_small
      implicit none       
      double precision x_array(*), y_array(*), x
      integer n_max, n_start, flag
      double precision interpolate, x_int(1), y_int(1)
      integer i, case, start, n
      ! check the order of the x-array
      if (flag .eq. 3) then
        goto 112
      end if
      if (n_start .eq. 0) then
        start=1
      else
        start=n_start
      end if
      if (x_array(start) .lt. x_array(n_max)) then
        case=+1
        if (n_start .gt. 2) then
          i=n_start
        else
          i=2
        end if
        if (x .lt. x_array(start) .or. x .gt. x_array(n_max)) then
          write(*,*) 'ERROR in find_value (+):', x, x_array(start), x_array(n_max), n_max, start, flag
          open (unit=64, file = 'function_error.dat', form = 'formatted', status = 'replace')
          do i=1, n_max
            write(64,*) x_array(i), y_array(i)
          end do
          close(64)
          stop
        end if
      else
        case=-1
        i=n_max-1
        if (x .gt. x_array(1) .or. x .lt. x_array(n_max)) then
          write(*,*) 'ERROR in find_value (-):', x, x_array(start), x_array(n_max), n_max, start, flag
          open (unit=64, file = 'function_error.dat', form = 'formatted', status = 'replace')
          do i=1, n_max
            write(64,*) x_array(i), y_array(i)
          end do
          close(64)
          stop
        end if
      end if
      ! look for the right position within the x-array
      do while (x .gt. x_array(i) .and. i .lt. n_max)
        i=i+case
        if (x_array(i) .lt. x_array(i-case)) then
          write(*,*) 'ERROR in find_value:', i, (i-case), x_array(i), x_array(i-case),x
          open (unit=64, file = 'function_error.dat', form = 'formatted', status = 'replace')
          do i=1, n_max
            write(64,*) x_array(i), y_array(i)
          end do
          close(64)
          stop
        end if
      end do
      if (flag .eq. 0) then ! standard interpolation
        find_value=interpolate(x_array(i),x_array(i-case),y_array(i), y_array(i-case),x)
      end if
      if (flag .eq. 1) then ! for threshold-like y-data
        if (abs(y_array(i-case)) .lt. val_small .and. x .ne. x_array(i)) then
          find_value=0.0d0
        else
          find_value=interpolate(x_array(i),x_array(i-case),y_array(i),y_array(i-case),x)
        end if
      end if
      if (flag .eq. 2) then ! spline interpolation by Akima (1970)
        if (y_array(i) .ne. y_array(i-case) .and. x .ne. x_array(i) .and. x .ne. x_array(i-case)) then
112       n=1
          x_int(1)=x
          call akima_interpol(3,n_max,x_array,y_array,n,x_int,y_int)
          find_value=y_int(1)
          if (isnan(find_value)) then
            write(*,*) 'warning isnan in akima'
            open (unit=64, file = 'function_error.dat', form = 'formatted', status = 'replace')
            do i=1, n_max
              write(64,*) x_array(i), y_array(i)
            end do
            close(64)
            find_value=interpolate(x_array(i),x_array(i-case),y_array(i), y_array(i-case),x)
          end if
        else
          find_value=y_array(i)
        end if
      end if
      if (isnan(find_value)) then
      write(*,*) 'ERROR in find_value:', i, (i-case), x_array(i), x_array(i-case),x, n_max
      open (unit=64, file = 'function_error.dat', form = 'formatted', status = 'replace')
      do i=1, n_max
        write(64,*) x_array(i), y_array(i)
      end do
      close(64)
      stop
      end if
      end



      double precision function interpolate(x1, x2, y1, y2, x)
      implicit none       
      double precision x1, x2, y1, y2, x
      if ((abs(x1-x) .le. abs(x1-x2)) .and. (abs(x2-x) .le. abs(x1-x2))) then
        if (x1 .ne. x2) then
          interpolate=y1+(y2-y1)/(x2-x1)*(x-x1)
        else
          interpolate=(y1+y2)/2.d0
        end if
      else
        write(*,*) 'ERROR in Interpolate:', x1, x2, y1, y2, x
        stop
      end if
      end



!     ELSEPA
!     Francesc Salvat, Aleksander Jablonski, Cedric J. Powell
!     F. Salvat et al., Comput. Phys. Commun. 165, 157 (2005)
      SUBROUTINE ELSCATA(EN,Z,ANGLE,CS,S,T,U,MT_FLAG,MT_RAD)
      use constants, only: RegInt_K, RegReal_K, PI
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX (KIND=RegReal_K) (C), INTEGER (KIND=RegInt_K) (I-N)
!
      PARAMETER (A0B=5.291772083D-9)  ! Bohr radius (cm)
      PARAMETER (A0B2=A0B*A0B)
!  ****  Results from the partial wave calculation.
      PARAMETER (NGT=650)
      double precision, intent(in) :: EN, MT_RAD
      integer, intent(in) :: MT_FLAG, Z
      double precision, intent(out) :: ANGLE(606), CS(606), S(606), T(606), U(606)

      COMMON/DCSTAB/ECS,TCS1,TCS2,TH(NGT),XT(NGT),DCST(NGT),SPOL(NGT),TPOL(NGT),UPOL(NGT),ERROR(NGT),NTAB
      COMMON/CTOTCS/TOTCS,ABCS
!
!  ****  Atomic polarizabilities of free atoms (in cm**3), from
!    Thomas M. Miller, 'Atomic and molecular polarizabilities' in
!    CRC Handbook of Chemistry and Physics, Editor-in-chief David
!    R. Linde. 79th ed., 1998-1999, pp. 10-160 to 10-174.
      DIMENSION ATPOL(103)
      DATA ATPOL/0.666D-24, 0.205D-24, 24.30D-24, 5.600D-24,3.030D-24, 1.760D-24, 1.100D-24, 0.802D-24, 0.557D-24,3.956D-25, 24.08D-24, 10.06D-24, 6.800D-24, 5.380D-24,3.630D-24, 2.900D-24, 2.180D-24, 1.641D-24, 43.40D-24,22.80D-24, 17.80D-24, 14.60D-24, 12.40D-24, 11.60D-24,9.400D-24, 8.400D-24, 7.500D-24, 6.800D-24, 6.100D-24,7.100D-24, 8.120D-24, 6.070D-24, 4.310D-24, 3.770D-24,3.050D-24, 2.484D-24, 47.30D-24, 27.60D-24, 22.70D-24,17.90D-24, 15.70D-24, 12.80D-24, 11.40D-24, 9.600D-24,8.600D-24, 4.800D-24, 7.200D-24, 7.200D-24, 10.20D-24,7.700D-24, 6.600D-24, 5.500D-24, 5.350D-24, 4.044D-24,59.60D-24, 39.70D-24, 31.10D-24, 29.60D-24, 28.20D-24,31.40D-24, 30.10D-24, 28.80D-24, 27.70D-24, 23.50D-24,25.50D-24, 24.50D-24, 23.60D-24, 22.70D-24, 21.80D-24,21.00D-24, 21.90D-24, 16.20D-24, 13.10D-24, 11.10D-24,9.700D-24, 8.500D-24, 7.600D-24, 6.500D-24, 5.800D-24,5.100D-24, 7.600D-24, 6.800D-24, 7.400D-24, 6.800D-24,6.000D-24, 5.300D-24, 48.70D-24, 38.30D-24, 32.10D-24,32.10D-24, 25.40D-24, 24.90D-24, 24.80D-24, 24.50D-24,23.30D-24, 23.00D-24, 22.70D-24, 20.50D-24, 19.70D-24,23.80D-24, 18.20D-24, 17.50D-24, 20.00D-24/
!
!  ****  Ionization energies of neutral atoms (in eV).
!        NIST Physical Reference Data.
!        http://sed.nist.gov/PhysRefData/IonEnergy/tblNew.html
!  For astatine (Z=85), the value given below was calculated with
!  the DHFXA code (Salvat and Fernandez-Varea, UBIR-2003).
      DIMENSION EIONZ(103)
      DATA EIONZ/13.5984D0, 24.5874D0, 5.39170D0, 9.32270D0,8.29800D0, 11.2603D0, 14.5341D0, 13.6181D0, 17.4228D0,21.5646D0, 5.13910D0, 7.64620D0, 5.98580D0, 8.15170D0,10.4867D0, 10.3600D0, 12.9676D0, 15.7596D0, 4.34070D0,6.11320D0, 6.56150D0, 6.82810D0, 6.74620D0, 6.76650D0,7.43400D0, 7.90240D0, 7.88100D0, 7.63980D0, 7.72640D0,9.39420D0, 5.99930D0, 7.89940D0, 9.78860D0, 9.75240D0,11.8138D0, 13.9996D0, 4.17710D0, 5.69490D0, 6.21710D0,6.63390D0, 6.75890D0, 7.09240D0, 7.28000D0, 7.36050D0,7.45890D0, 8.33690D0, 7.57620D0, 8.99380D0, 5.78640D0,7.34390D0, 8.60840D0, 9.00960D0, 10.4513D0, 12.1298D0,3.89390D0, 5.21170D0, 5.57690D0, 5.53870D0, 5.47300D0,5.52500D0, 5.58200D0, 5.64360D0, 5.67040D0, 6.15010D0,5.86380D0, 5.93890D0, 6.02150D0, 6.10770D0, 6.18430D0,6.25420D0, 5.42590D0, 6.82510D0, 7.54960D0, 7.86400D0,7.83350D0, 8.43820D0, 8.96700D0, 8.95870D0, 9.22550D0,10.4375D0, 6.10820D0, 7.41670D0, 7.28560D0, 8.41700D0,9.50000D0, 10.7485D0, 4.07270D0, 5.27840D0, 5.17000D0,6.30670D0, 5.89000D0, 6.19410D0, 6.26570D0, 6.02620D0,5.97380D0, 5.99150D0, 6.19790D0, 6.28170D0, 6.42000D0,6.50000D0, 6.58000D0, 6.65000D0, 4.90000D0/
!
!  ****  First excitation energies of neutral atoms (in eV).
!        NIST Physical Reference Data.
!  The value 0.0D0 indicates that the experimental value for
!  the atom was not available.
      DIMENSION EEX1Z(103)
      DATA EEX1Z/10.20D0, 19.82D0,  1.85D0,  2.73D0,3.58D0,  1.26D0,  2.38D0,  1.97D0, 12.70D0,16.62D0,  2.10D0,  2.71D0,  3.14D0,  0.78D0,1.41D0,  1.15D0,  8.92D0, 11.55D0,  1.61D0,1.88D0,  1.43D0,  0.81D0,  0.26D0,  0.94D0,2.11D0,  0.86D0,  0.43D0,  0.01D0,  1.38D0,4.00D0,  3.07D0,  0.88D0,  1.31D0,  1.19D0,7.87D0,  9.91D0,  0.00D0,  0.00D0,  0.00D0,0.00D0,  0.00D0,  1.34D0,  0.00D0,  0.00D0,0.00D0,  0.00D0,  0.00D0,  3.73D0,  0.00D0,0.00D0,  0.00D0,  0.00D0,  0.00D0,  8.31D0,0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,4.67D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,0.00D0,  0.00D0,  0.00D0,  0.00D0/
!
!  ****  Nearest-neighbour distances (in cm) of the elements,
!    from Ch. Kittel, 'Introduction to Solid State Physics'. 5th
!    ed. (John Wiley and Sons, New York, 1976).
!  The value -1.0D-8 indicates that the experimental value for
!  the element was not available.
      DIMENSION DNNEL(103)
      DATA DNNEL/-1.000D-8,-1.000D-8, 3.124D-8, 2.220D-8,-1.000D-8, 1.540D-8,-1.000D-8,-1.000D-8, 1.440D-8,-1.000D-8, 3.822D-8, 3.200D-8, 2.860D-8, 2.350D-8,-1.000D-8,-1.000D-8,-1.000D-8,-1.000D-8, 4.752D-8,3.950D-8, 3.250D-8, 2.890D-8, 2.620D-8, 2.500D-8,2.240D-8, 2.480D-8, 2.500D-8, 2.490D-8, 2.560D-8,2.660D-8, 2.440D-8, 2.450D-8, 3.160D-8, 2.320D-8,-1.000D-8,-1.000D-8, 4.837D-8, 4.300D-8, 3.550D-8,3.170D-8, 2.860D-8, 2.720D-8, 2.710D-8, 2.650D-8,2.690D-8, 2.750D-8, 2.890D-8, 2.980D-8, 3.250D-8,2.810D-8, 2.910D-8, 2.860D-8, 3.540D-8,-1.000D-8,5.235D-8, 4.350D-8, 3.730D-8, 3.650D-8, 3.630D-8,3.660D-8,-1.000D-8, 3.590D-8, 3.960D-8, 3.580D-8,3.520D-8, 3.510D-8, 3.490D-8, 3.470D-8, 3.540D-8,3.880D-8, 3.430D-8, 3.130D-8, 2.860D-8, 2.740D-8,2.740D-8, 2.680D-8, 2.710D-8, 2.770D-8, 2.880D-8,-1.000D-8, 3.460D-8, 3.500D-8, 3.070D-8, 3.340D-8,-1.000D-8,-1.000D-8,-1.000D-8,-1.000D-8, 3.760D-8,3.600D-8, 3.210D-8, 2.750D-8, 2.620D-8, 3.100D-8,3.610D-8,-1.000D-8,-1.000D-8,-1.000D-8,-1.000D-8,-1.000D-8,-1.000D-8,-1.000D-8,-1.000D-8/
!  ****  Phase shifts.
      PARAMETER (NDM=25000)
      COMMON/PHASES/DP(NDM),DM(NDM),NPH,ISUMP
      COMMON/PHASEI/DPJ(NDM),DMJ(NDM)
!
!  ************  Input data.
!
!  ****  Default model.
! ----+----1----+----2----+----3----+----4----+----5----+----6----+----7
! IZ      80         atomic number                               [none]
! MNUCL   3          rho_n (1=P, 2=U, 3=F, 4=Uu)                  [  3]
! NELEC   80         number of bound electrons                    [ IZ]
! MELEC   4          rho_e (1=TFM, 2=TFD, 3=DHFS, 4=DF, 5=file)   [  4]
! MUFFIN  0          0=free atom, 1=muffin-tin model              [  0]
! RMUF    0          muffin-tin radius (cm)                  [measured]
! IELEC  -1          -1=electron, +1=positron                     [ -1]
! MEXCH   1          V_ex (0=none, 1=FM, 2=TF, 3=RT)              [  1]
! MCPOL   2          V_cp (0=none, 1=B, 2=LDA)                    [  0]
! VPOLA  -1          atomic polarizability (cm**3)           [measured]
! VPOLB  -1          b_pol parameter                          [default]
! MABS    1          W_abs (0=none, 1=LDA)                        [  0]
! VABSA   2.00       absorption-potential strength, Aabs          [2.0]
! VABSD  -1.0        energy gap DELTA (eV)                    [default]
! IHEF    2          high-E factorization (0=no, 1=yes, 2=Born)   [  1]
! EV      1.000E2    kinetic energy (eV)                         [none]
! ----+----1----+----2----+----3----+----4----+----5----+----6----+----7
      IELEC =-1        ! electron
      IZ    = 0        ! no default
      NELEC = 1000     ! =Z (the present value is a flag)
      MNUCL = 3        ! Fermi nuclear charge distribution
      MELEC = 4        ! DF electron density
      MUFFIN= 0        ! free atom
      RMUF  =-200.0D0  ! free atom
      MEXCH = 1        ! FM exchange potential
      MCPOL = 0        ! no correlation-polarization
      VPOLA =-1.0D0    ! atomic polarizability
      VPOLB =-1.0D0    ! polariz. cutoff parameter
      MABS  = 1        ! no absorption
      VABSA = 2.0D0    ! absorption potential strength
      VABSD =-1.0D0    ! energy gap
      IHEF  = 0        ! high-energy factorization on
!
        MUFFIN=MT_FLAG
        RMUF=MT_RAD
        EV=EN
        IZ=Z
        NELEC=Z
        IF(IZ.LT.1.OR.IZ.GT.103) THEN
          WRITE(6,*) 'IZ =',IZ
          STOP 'Wrong atomic number'
        ENDIF
        IF(MNUCL.LT.1.OR.MNUCL.GT.4) MNUCL=3
        IF(NELEC.LT.0) NELEC=IZ
        IF(NELEC.GT.IZ) THEN
          WRITE(6,*) 'NELEC =',NELEC
          STOP 'Negative ion'
        ENDIF
        IF(MELEC.LT.1.OR.MELEC.GT.5) MELEC=4
        IF(MUFFIN.NE.1) MUFFIN=0
        IF(IELEC.NE.+1) IELEC=-1
        IF(MEXCH.LT.0.OR.MEXCH.GT.3) MEXCH=1
        IF(MCPOL.LT.0.OR.MCPOL.GT.2) MCPOL=0
        IF(VPOLA.LT.-1.0D-35) VPOLA=ATPOL(IZ)
        IF(VPOLB.LT.1.0D-10) VPOLB=-10.0D0
        IF(MABS.NE.1) MABS=0
        IF(VABSAI.GE.0.0D0) VABSA=VABSAI
        IF(VABSD.LT.-1.0D-35) VABSD=-1.0D0
        IF(IHEF.NE.0.AND.IHEF.NE.2) IHEF=1
!
!  ****  Potential model parameters.
!
  200 CONTINUE
      IF(NELEC.EQ.1000) NELEC=IZ
!
      IF(NELEC.NE.IZ) MUFFIN=0
      IF(MUFFIN.EQ.1) THEN
        IF(RMUF.LT.1.0D-9) RMUF=0.5D0*DNNEL(IZ)
        IF(RMUF.LT.1.0D-9) STOP 'RMUF is too small'
      ELSE
        MUFFIN=0
        RMUF=200.0D-8
      ENDIF
!
      IF(IELEC.EQ.+1) MEXCH=0
!
      IF(MCPOL.EQ.1) THEN
        IF(VPOLA.LT.-1.0D-35) VPOLA=ATPOL(IZ)  ! Default value.
        IF(VPOLA.LT.-1.0D-35) STOP 'VPOLA must be positive'
        IF(VPOLB.LT.1.0D-10) VPOLB=-10.0D0
      ELSE IF(MCPOL.EQ.2) THEN
        IF(VPOLA.LT.-1.0D-35) VPOLA=ATPOL(IZ)
        IF(VPOLA.LT.-1.0D-35) STOP 'VPOLA must be positive'
        IF(VPOLB.LT.1.0D-10) VPOLB=-10.0D0
      ELSE
        MCPOL=0
        VPOLA=0.0D0
        VPOLB=0.0D0
      ENDIF
!
      IF(MABS.EQ.1) THEN
        IF(VABSA.LT.-1.0D-35) STOP 'VABSA must be positive'
        IF(VABSD.LT.-1.0D-35) THEN
          IF(IELEC.EQ.-1) THEN  ! Default (experimental) value.
            VABSD=EEX1Z(IZ)
            IF(VABSD.LT.-1.0D-35) VABSD=0.0D0
          ELSE
            VABSD=MAX(0.0D0,EIONZ(IZ)-6.8D0)
          ENDIF
        ENDIF
      ELSE
        MABS=0
        VABSA=0.0D0
        VABSD=1.0D0
      ENDIF
      ECUT=MIN(20.0D3*dble(IZ),2.0D6)

!
!  ************  Partial-wave analysis.
!
      IF(EV.LT.10.0D0) STOP 'The kinetic energy is too small'
!
      IF(MCPOL.NE.0) THEN
        IF(EV.GT.1.0D4) THEN
!          WRITE(6,*) 'WARNING: For E>10 keV the correlation-polarization correction is'
!          WRITE(6,*) '         switched off.'
          MCPOLC=0
          VPOLBC=0.0D0
        ELSE
          MCPOLC=MCPOL
          IF(VPOLB.LT.0.0D0) THEN
            VPOLBC=0.25D0*DSQRT(MAX(EV-50.0D0,16.0D0))
          ELSE
            VPOLBC=VPOLB
          ENDIF
        ENDIF
      ELSE
        MCPOLC=0
        VPOLBC=0.0D0
      ENDIF
      IF(MABS.NE.0) THEN
        IF(EV.GT.1.0D6) THEN
!          WRITE(6,*) 'WARNING: For E>1 MeV, the absorption correction is switched off.'
          MABSC=0
        ELSE
          MABSC=MABS
        ENDIF
      ELSE
        MABSC=0
      ENDIF
!
      CALL ELSEPA(IELEC,EV,IZ,NELEC,MNUCL,MELEC,MUFFIN,RMUF,MEXCH,MCPOLC,VPOLA,VPOLBC,MABSC,VABSA,VABSD,IHEF,8)

      DO I=1, 606
        ANGLE(I)=TH(I)
        CS(I)=DCST(I)
        S(I)=SPOL(I)
        T(I)=TPOL(I)
        U(I)=UPOL(I)
      END DO
!
      END SUBROUTINE ELSCATA






!  NOTE: The present subroutine package uses I/O units 97, 98 and 99.
!        Do not use these unit numbers in your main program.
!
!  *********************************************************************
!                       SUBROUTINE ELSEPA
!  *********************************************************************
      SUBROUTINE ELSEPA(IELEC,EV,IZ,NELEC,MNUCL,MELEC,MUFIN,RMUF,MEXCH,MCPOL,VPOLA,VPOLB,MABS,VABSA,VABSD,IHEF,IW)
!
!
!                              F. Salvat, A. Jablonski and !.J. Powell
!                              September 27, 2004
!
!
!     This subroutine computes scattering amplitudes, differential cross
!  sections and total (integrated) cross sections for ELastic Scattering
!  of Electrons and Positrons by neutral Atoms and positive ions.
!
!     The interaction is described through a static (central) field,
!  which consists of the electrostatic potential and, for projectile
!  electrons, an approximate local exchange potential. For slow
!  projectiles, a correlation-polarization potential and an absorptive
!  imaginary potential can optionally be included. The differential
!  cross section is evaluated by means of relativistic (Dirac) partial-
!  wave analysis, or from approximate high-energy factorizations.
!
!  Input arguments:
!    IELEC ..... electron-positron flag;
!                =-1 for electrons,
!                =+1 for positrons.
!    EV ........ projectile's kinetic energy (in eV).
!    IZ ........ atomic number of the target atom or ion.
!    NELEC ..... number of bound atomic electrons.
!    MNUCL ..... nuclear charge density model.
!                  1 --> point nucleus (P),
!                  2 --> uniform distribution (U),
!                  3 --> Fermi distribution (F),
!                  4 --> Helm's uniform-uniform distribution (Uu).
!    MELEC ..... electron density model.
!                  1 --> TFM analytical density,
!                  2 --> TFD analytical density,
!                  3 --> DHFS analytical density,
!                  4 --> DF numerical density, read from 'Z_zzz.DEN',
!                  5 --> density read from file 'density.usr'.
!    MUFIN ..... Aggregation effects...
!                  0 --> free atom,
!                  1 --> muffin-tin model.
!      RMUF .... Muffin-tin radius (in cm).
!    MEXCH ..... exchange correction for electrons.
!                  0 --> no exchange correction,
!                  1 --> Furness-McCarthy (FM),
!                  2 --> Thomas-Fermi (TF),
!                  3 --> Riley-Truhlar (RT).
!    MCPOL ..... correlation-polarization correction.
!                  0 --> no correlation-polarization correction,
!                  1 --> Buckingham potential (B),
!                  2 --> Local density approximation (LDA).
!      VPOLA ... atomic polarizability (in cm**3).
!      VPOLB ... cutoff radius parameter b_pol
!                    (used only when MCPOL>0).
!    MABS ...... absorption correction (imaginary potential).
!                  0 --> no absorption correction,
!                  1 --> LDA (electron-hole excitations only).
!      VABSA ... strength of the absorption potential.
!      VABSD ... energy gap, DELTA (eV).
!                    (used only when MABS is different from 0).
!    IHEF ...... =0: phase shifts are computed for the electrostatic
!                    field of the whole atom (nucleus+electron cloud)
!                    with optional exchange, polarization and absorption
!                    corrections.
!                =1: the differential cross section is obtained from a
!                    high-energy factorization. The phase shifts are
!                    evaluated for the bare nucleus. The screening of
!                    the nuclear charge by the atomic electrons is
!                    accounted for by means of a pre-evaluated high-
!                    energy correction factor, which is read from file
!                    'Z_zzz.DFS'.
!                =2: when the energy is larger than 100 MeV, the DCS is
!                    obtained as the product of the Mott DCS for a point
!                    nucleus, the Helm Uu nuclear form factor (with an
!                    empirical Coulomb correction) and the electron
!                    screening factor.
!    IW ........ output unit (to be defined in the main program).
!
!  The electrostatic potential and electron density of the target atom
!  or ion are calculated by subroutine EFIELD and delivered through the
!  the named common block /CFIELD/.
!
!  Output (through the common block /DCSTAB/):
!     ECS ........ total cross section (cm**2).
!     TCS1 ....... 1st transport cross section (cm**2).
!     TCS2 ....... 2nd transport cross section (cm**2).
!     TH(I) ...... scattering angles (in deg)
!     XT(I) ...... values of (1-COS(TH(I)))/2.
!     DCST(I) .... differential cross section per unit solid
!                  angle at TH(I) (in cm**2/sr).
!     SPOL(I) .... Sherman spin-polarization function at TH(I).
!     ERROR(I) ... relative uncertainty of the computed DCS
!                  values. Estimated from the convergence of the
!                  series.
!     NTAB ....... number of angles in the table.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX (KIND=RegReal_K) (C), INTEGER (KIND=RegInt_K) (I-N)
!
!  ****  The parameter IWR defines the amount of information printed on
!  output files:
!  IWR>0 => the scattering potential is printed on file 'scfield.dat'.
!  IWR>1 => the scattering amplitudes are printed on file 'scatamp.dat'.
      PARAMETER (IWR=2)
!
      PARAMETER (A0B=5.291772083D-9)  ! Bohr radius (cm)
      PARAMETER (HREV=27.2113834D0)  ! Hartree energy (eV)
      PARAMETER (REV=5.10998902D5)  ! Electron rest energy (eV)
      PARAMETER (SL=137.03599976D0)  ! Speed of light (1/alpha)
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (TREV=REV+REV)
      PARAMETER (PI=3.1415926535897932D0,FOURPI=4.0D0*PI)
!
      character (len=34) :: SCFILE,NULL
      character (len=1) :: LIT10(10),LIT1,LIT2,LIT3
      DATA LIT10/'0','1','2','3','4','5','6','7','8','9'/
!
      PARAMETER (NGT=650)
      COMMON/DCSTAB/ECS,TCS1,TCS2,TH(NGT),XT(NGT),DCST(NGT),SPOL(NGT),TPOL(NGT),UPOL(NGT),ERROR(NGT),NTAB
      COMMON/CTOTCS/TOTCS,ABCS
      DIMENSION Q2T(NGT),FQ(NGT)
!
      PARAMETER (NDIM=1000)
      COMMON/CFIELD/R(NDIM),RVN(NDIM),DEN(NDIM),RVST(NDIM),NPOT
      COMMON/FIELD/RAD(NDIM),RV(NDIM),NP
      COMMON/FIELDI/RADI(NDIM),RVI(NDIM),RW(NDIM),IAB,NPI
      DIMENSION RVEX(NDIM),RVPOL(NDIM)
!
      PARAMETER (NPC=1500,NDM=25000)
      COMMON/WORK/XL(NDM),SA(NDM),SB(NDM),SC(NDM),SD(NDM),PL(NDM)
      COMMON/CRMORU/CFM(NPC),CGM(NPC),DPC(NPC),DMC(NPC),CF,CG,RUTHC,WATSC,RK2,ERRFC,ERRGC,NPC1
      DIMENSION DENA(NPC),DENB(NPC),DENC(NPC),DEND(NPC)
!
!  ----  Mott DCS and spin polarization (point unscreened nucleus)
      IF(NELEC.EQ.0.AND.MNUCL.EQ.1) THEN
        CALL MOTTSC(IELEC,IZ,EV,IW)
        RETURN
      ENDIF
!  ----  High-energy Mott-Born approximation.
      IF(EV.GT.100.0D6.AND.IHEF.EQ.2.AND.IZ.EQ.NELEC) THEN
        CALL HEBORN(IELEC,IZ,MNUCL,EV,IW)

        RETURN
      ENDIF
!
      WRITE(IW,1000)
 1000 FORMAT(1X,'#',/1X,'# Subroutine ELSEPA. Elastic scattering of ','electrons and positrons',/1X,'#',20X,'by neutral atoms and positive ions')
      IF(IELEC.EQ.-1) THEN
        WRITE(IW,1100)
 1100   FORMAT(1X,'#',/1X,'# Projectile: electron')
      ELSE
        WRITE(IW,1200)
 1200   FORMAT(1X,'#',/1X,'# Projectile: positron')
      ENDIF
      E=EV/HREV
      WRITE(IW,1300) EV,E
 1300 FORMAT(1X,'# Kinetic energy =',1P,E12.5,' eV =',E12.5,' a.u.')
!
!  ****  You may wish to comment off the next condition to run the
!  program for kinetic energies less that 10 eV. However, the results
!  for these energies may be highly inaccurate.
!
      IF(EV.LT.10.0D0) THEN
        WRITE(*,'(''  STOP. The kinetic energy is too small.'')')
        STOP 'The kinetic energy is too small.'
      ENDIF
      IF(EV.LT.100.0D0) THEN
        WRITE(*,'(1X,''#'',/1X,''#  ***  WARNING: Energy is '',''too low.'',/1X,''#'',16X,''The reliability of the '',''results is questionable.'')')
      ENDIF
!
      IHEF0=0
      IF(EV.GT.20.1D3*IZ.AND.IHEF.GT.0.AND.IZ.EQ.NELEC) THEN
        IF(EV.GT.1.0D6.OR.MABS.EQ.0) THEN
          IHEF0=1
          GO TO 100
        ENDIF
      ENDIF
!
!  ****  Electrostatic field.
!
      IF(MNUCL.LT.1.OR.MNUCL.GT.4) THEN
        WRITE(*,*) 'ELSEPA: incorrect MNUCL value.'
        STOP 'ELSEPA: incorrect MNUCL value.'
      ENDIF
      IF(MELEC.LT.1.OR.MELEC.GT.5) THEN
        WRITE(*,*) 'ELSEPA: incorrect MELEC value.'
        STOP 'ELSEPA: incorrect MELEC value.'
      ENDIF
      CALL EFIELD(IZ,NELEC,NELEC,MNUCL,MELEC,IW)
      ZINF=IELEC*DBLE(IZ-NELEC)
      DO I=1,NPOT
        RV(I)=DBLE(IELEC)*RVST(I)
      ENDDO
      RV(NPOT)=ZINF
!
!  ************  Muffin-tin model for scattering in solids.
!
      NMT=0
      IF(MUFIN.EQ.1) THEN
!        WRITE(*,1600) RMUF
! 1600   FORMAT(1X,'#',/1X,'# Muffin-tin model: Rmt =',1P,E12.5,' cm')
        CALL SPLINE(R,RV,SA,SB,SC,SD,0.0D0,0.0D0,NPOT)
        CALL SPLINE(R,DEN,DENA,DENB,DENC,DEND,0.0D0,0.0D0,NPOT)
        RMT=RMUF/A0B
        IF(RMT.LT.R(NPOT)) THEN
          CALL FINDI(R,RMT,NPOT,J)
          IF(J.LT.5) STOP 'The muffin-tin radius is too small.'
          DENRMT=DENA(J)+RMT*(DENB(J)+RMT*(DENC(J)+RMT*DEND(J)))
        ELSE
          RMT=R(NPOT)
          DENRMT=DEN(NPOT)
        ENDIF
        RHORMT=2.0D0*DENRMT/(FOURPI*RMT**2)
        DO I=1,NPOT
          IF(R(I).GT.RMT) THEN
            IF(RAD(I-1).LT.RMT*0.9999999D0) THEN
              NP=I
              RAD(NP)=RMT
            ELSE
              NP=I-1
              RAD(NP)=RMT
            ENDIF
            RC1=RMT
            CALL FINDI(R,RC1,NPOT,J)
            V1=SA(J)+RC1*(SB(J)+RC1*(SC(J)+RC1*SD(J)))
            DEN1=DENA(J)+RC1*(DENB(J)+RC1*(DENC(J)+RC1*DEND(J)))
            RV(NP)=2.0D0*V1
            DEN(NP)=2.0D0*DEN1
            RVST(NP)=RV(NP)/DBLE(IELEC)
            GO TO 1
          ELSE
            RAD(I)=R(I)
!
            RC1=R(I)
            FD1=FOURPI*RC1**2
            V1=RV(I)
            DEN1=DEN(I)
!
            RC2=2.0D0*RMT-R(I)
            FD2=FOURPI*RC2**2
            CALL FINDI(R,RC2,NPOT,J)
            V2=SA(J)+RC2*(SB(J)+RC2*(SC(J)+RC2*SD(J)))
            DEN2=DENA(J)+RC2*(DENB(J)+RC2*(DENC(J)+RC2*DEND(J)))
!
            IF(I.GT.1) THEN
              RV(I)=V1+RC1*(V2/RC2)
              DEN(I)=DEN1+FD1*(DEN2/FD2)
            ELSE
              RV(I)=V1
              DEN(I)=DEN1
            ENDIF
          ENDIF
          RVST(I)=RV(I)/DBLE(IELEC)
        ENDDO
        NP=NPOT
    1   CONTINUE
!
!  ****  Ensure proper normalization of the muffin-tin electron density.
!
        SUM=RMOM(RAD,DEN,NP,0)
        RHOU=(NELEC-SUM)/(FOURPI*RMT**3/3.0D0)
        DO I=1,NP
          DEN(I)=DEN(I)+RHOU*FOURPI*RAD(I)**2
        ENDDO
        RHORMT=RHORMT+RHOU
        SUM=RMOM(RAD,DEN,NP,0)
!        WRITE(6,*) 'Electron density normalization =',SUM
!
        NMT=NP
        E=E-RV(NMT)/RAD(NMT)
        DO I=2,NMT
          RV(I)=RV(I)-RV(NMT)*RAD(I)/RAD(NMT)
          RVST(I)=RVST(I)-RVST(NMT)*RAD(I)/RAD(NMT)
        ENDDO
        IF(NP.LT.NPOT) THEN
          NP=NP+1
          DO I=NP,NPOT
            IF(I.EQ.NP) RAD(I)=RAD(I-1)
            RV(I)=0.0D0
            RVST(I)=0.0D0
            DEN(I)=0.0D0
          ENDDO
        ENDIF
      ELSE
        NP=NPOT
        DO I=1,NPOT
          RAD(I)=R(I)
        ENDDO
      ENDIF
!
!  ************  Exchange correction for electrons.
!
      IF(IELEC.EQ.-1.AND.MEXCH.NE.0) THEN
        IF(MEXCH.EQ.1) THEN
!  ****  Furness-McCarthy exchange potential.
          WRITE(IW,1500)
 1500     FORMAT(1X,'#',/1X,'# Furness-McCarthy exchange',' potential')
          DO I=2,NP
            AUX=RAD(I)*E*(1.0D0+EV/TREV)+RVST(I)
            AUX2=AUX*AUX
            IF(DEN(I).GT.1.0D-5*AUX2) THEN
              RVEX(I)=0.5D0*(AUX-DSQRT(AUX2+DEN(I)))
            ELSE
              T=DEN(I)/AUX2
              RVEX(I)=-0.5D0*AUX*T*(0.5D0-T*(0.125D0-T*0.065D0))
            ENDIF
            RV(I)=RV(I)+RVEX(I)
          ENDDO
        ELSE IF(MEXCH.EQ.2) THEN
!  ****  Thomas-Fermi exchange potential.
          WRITE(IW,1400)
 1400     FORMAT(1X,'#',/1X,'# Thomas-Fermi exchange potential')
          DO I=1,NP
            RHO=DEN(MAX(2,I))/(FOURPI*RAD(MAX(2,I))**2)
            SKF=(3.0D0*PI*PI*RHO)**3.333333333333333D-1
            EF=0.5D0*SKF*SKF
            SKL=SQRT(2.0D0*(E*(1.0D0+EV/TREV)+EF))
            X=SKF/SKL
            IF(X.LT.0.001D0) THEN
              FX=(2.0D0/3.0D0)*X**3
            ELSE
              FX=X-0.5D0*(1.0D0-X*X)*LOG(ABS((1.0D0+X)/(1.0D0-X)))
            ENDIF
            RVEX(I)=-(SKL/PI)*FX*RAD(I)
            RV(I)=RV(I)+RVEX(I)
          ENDDO
        ELSE IF(MEXCH.EQ.3) THEN
!  ****  Riley-Truhlar exchange potential.
          WRITE(IW,1501)
 1501     FORMAT(1X,'#',/1X,'# Riley-Truhlar exchange potential')
          DO I=1,NP
            AUX=4.0D0*(RAD(I)*E*(1.0D0+EV/TREV)+RVST(I))
            IF(AUX.GT.1.0D-16*DEN(I)) THEN
              RVEX(I)=-DEN(I)/AUX
              RV(I)=RV(I)+RVEX(I)
            ENDIF
          ENDDO
        ELSE
          WRITE(*,*) 'ELSEPA: incorrect MEXCH value.'
          STOP 'ELSEPA: incorrect MEXCH value.'
        ENDIF
        IF(NMT.GT.1) THEN
          E=E-RV(NMT)/RAD(NMT)
          DO I=2,NMT
            RV(I)=RV(I)-RAD(I)*(RV(NMT)/RAD(NMT))
            RVEX(I)=RVEX(I)-RAD(I)*(RVEX(NMT)/RAD(NMT))
          ENDDO
        ENDIF
      ELSE
        DO I=1,NP
          RVEX(I)=0.0D0
        ENDDO
      ENDIF
!
!  ********  Absorption potential (local density approximation).
!
      IF(MABS.EQ.1.AND.VABSA.GT.1.0D-12.AND.EV.LE.1.0D6) THEN
        IAB=1
        WRITE(IW,1502) VABSD,VABSA
 1502   FORMAT(1X,'#',/1X,'# LDA absorption potential (only electr','on-hole excitations):',/1X,'#',27X,'Delta =',1P,E12.5,' eV',/1X,'#',28X,'Aabs =',E12.5)
        DELTA=VABSD/HREV
        AABS=VABSA

!
        RADI(1)=RAD(1)
        RVI(1)=RV(1)
        RW(1)=0.0D0
        RVPOL(1)=0.0D0
        DO I=2,NP
          RADI(I)=RAD(I)
          RVI(I)=RV(I)
          RHO=DEN(I)/(FOURPI*RAD(I)**2)
          EKIN=E-RV(I)/RAD(I)
          IF(RHO.GT.1.0D-16.AND.EKIN.GT.DELTA) THEN
            VEL=DSQRT(2.0D0*EKIN)
            EKEV=EKIN*HREV
            FREL=DSQRT(2.0D0*(EKEV+REV)**2/(REV*(EKEV+2.0D0*REV)))
            CALL XSFEG(RHO,DELTA,IELEC,EKIN,0,XSEC,2)
            RW(I)=-0.5D0*VEL*RHO*XSEC*RAD(I)*AABS*FREL
          ELSE
            RW(I)=0.0D0
          ENDIF
          RVPOL(I)=0.0D0
          WRITE(6,1503) I,RADI(I),RW(I)
 1503     FORMAT(1X,'i=',I4,',   r=',1P,E12.5,',   r*Wabs= ',E12.5)
        ENDDO
        NPI=NP
      ELSE
        IAB=0
        DO I=1,NP
          RADI(I)=RAD(I)
          RVI(I)=RV(I)
          RW(I)=0.0D0
          RVPOL(I)=0.0D0
        ENDDO
        NPI=NP
      ENDIF
!
!  ****  We add a 'constant' tail to the potential to extend the grid
!  up to a point where irregular Coulomb functions can be be calculated.
!
      IF((MCPOL.NE.1.AND.MCPOL.NE.2).OR.EV.GT.1.0D4) THEN
        IF(NP.LT.NDIM-4) THEN
          IF(RAD(NP)-RAD(NP-1).LT.1.0D-16) THEN
            I=NP
            RAD(I)=RAD(I-1)
            RADI(I)=RADI(I-1)
          ELSE

            I=NP+1
            RAD(I)=RAD(I-1)
            RADI(I)=RADI(I-1)
            RV(I)=RV(I-1)
          ENDIF
          RV(I)=ZINF
          RVI(I)=RV(I)
          DEN(I)=0.0D0
          RVEX(I)=0.0D0
          RVPOL(I)=0.0D0
          RW(I)=0.0D0
          IST=I+1
          NADD=1
          DO I=IST,NDIM
            NADD=NADD+1
            RAD(I)=2.0D0*RAD(I-1)
            RADI(I)=2.0D0*RADI(I-1)
            RV(I)=ZINF
            RVI(I)=RV(I)
            DEN(I)=0.0D0
            RVEX(I)=0.0D0
            RVPOL(I)=0.0D0
            RW(I)=0.0D0
            NP=I
            NPI=I
            IF(RAD(I).GT.1.0D4.AND.NADD.GT.4) GO TO 2
          ENDDO
        ELSE
          STOP 'Not enough memory space 1.'
        ENDIF
    2   CONTINUE
      ELSE IF(MCPOL.EQ.1) THEN
!
!  ************  Atomic polarizability correction.
!
!  ****  Buckingham empirical potential.
        WRITE(IW,1700) VPOLA,VPOLB
 1700   FORMAT(1X,'#',/1X,'# Correlation-polarization potential (Buckingham):',/1X,'#',27X,'Alpha =',1P,E12.5,' cm**3',/1X,'#',28X,'Bpol =',E12.5)
        IF(VPOLB.LT.0.01D0) THEN
          WRITE(*,*) 'ELSEPA: VPOLB cannot be less than 0.01.'
          STOP 'ELSEPA: VPOLB cannot be less than 0.01.'
        ENDIF
        ALPHA=VPOLA/A0B**3
        D2=SQRT(0.5D0*ALPHA*VPOLB**2/DBLE(IZ)**3.333333333333333D-1)
        NPOL=NP
        DO I=1,NPOL
          VPOL=-0.5D0*ALPHA/(RAD(I)**2+D2)**2
          RVPOL(I)=VPOL*RAD(I)
          RV(I)=RV(I)+RVPOL(I)
          RVI(I)=RV(I)
        ENDDO
        IF(NPOL.LT.NDIM-5) THEN
          DO I=NPOL+1,NDIM-5
            RAD(I)=1.25D0*RAD(I-1)
            RADI(I)=RAD(I)
            VPOL=-0.5D0*ALPHA/(RAD(I)**2+D2)**2
            RVPOL(I)=VPOL*RAD(I)
            RV(I)=ZINF+RVPOL(I)
            RVI(I)=RV(I)
            DEN(I)=0.0D0
            RVEX(I)=0.0D0
            RW(I)=0.0D0
            NP=I
            NPI=I
            IF(ABS(VPOL).LT.1.0D-6*MAX(E,1.0D1*ABS(ZINF)/RAD(I)).AND.RAD(I).GT.50.0D0) GO TO 3
          ENDDO
        ENDIF
        STOP 'Not enough memory space 2.'
    3   CONTINUE
        IF(NP.LT.NDIM-4) THEN
          I=NP+1
          RAD(I)=RAD(I-1)
          RADI(I)=RADI(I-1)
          RV(I)=ZINF
          RVI(I)=RV(I)
          DEN(I)=0.0D0
          RVEX(I)=0.0D0
          RW(I)=0.0D0
          RVPOL(I)=0.0D0
          DO I=NP+2,NDIM
            RAD(I)=2.0D0*RAD(I-1)
            RADI(I)=2.0D0*RADI(I-1)
            RV(I)=RV(I-1)
            RVI(I)=RV(I)
            DEN(I)=0.0D0
            RVEX(I)=0.0D0
            RW(I)=0.0D0
            RVPOL(I)=0.0D0
            NP=I
            NPI=I
            IF(RAD(I).GT.1.0D4) GO TO 33
          ENDDO
        ELSE
          STOP 'Not enough memory space 3.'
        ENDIF
   33   CONTINUE
!
      ELSE IF(MCPOL.EQ.2) THEN
!  ****  LDA correlation-polarization potential.
        WRITE(IW,1701) VPOLA,VPOLB
 1701   FORMAT(1X,'#',/1X,'# Correlation-polarization potential (LDA','):',/1X,'#',27X,'Alpha =',1P,E12.5,' cm**3',/1X,'#',28X,'Bpol =',E12.5)
        IF(VPOLB.LT.0.01D0) THEN
          WRITE(IW,*) 'ELSEPA: VPOLB cannot be less than 0.01.'
          STOP 'ELSEPA: VPOLB cannot be less than 0.01.'
        ENDIF
        ALPHA=VPOLA/A0B**3
        D2=SQRT(0.5D0*ALPHA*VPOLB**2/DBLE(IZ)**3.333333333333333D-1)
        NPOL=NP
        IMODE=0
        IF(MUFIN.NE.1) THEN
          DO I=NPOL,1,-1
            RIP=RAD(MAX(2,I))
            RHO=DEN(MAX(2,I))/(FOURPI*RIP**2)
            VCO=VCPOL(IELEC,RHO)
            VPAS=-0.5D0*ALPHA/(RIP**2+D2)**2
            IF(IMODE.EQ.0) THEN
              VPOL=VPAS
              IF(VCO.LT.VPAS) IMODE=1
            ELSE
              VPOL=MAX(VCO,VPAS)
            ENDIF
            RVPOL(I)=VPOL*RAD(I)
            RV(I)=RV(I)+RVPOL(I)
            RVI(I)=RV(I)
          ENDDO
          IF(IMODE.EQ.0) THEN
            WRITE(*,1702)
 1702       FORMAT(1X,'#',/1X,'# ERROR: The correlation and polarization potentials do not cross.')
            STOP
          ENDIF
          VCOUT=-1.0D16
        ELSE
          DO I=1,NPOL
            RIP=RAD(MAX(2,I))
            RHO=DEN(MAX(2,I))/(FOURPI*RIP**2)
            IF(RHO.LT.1.0D-35) RHO=RHORMT
            VCO=VCPOL(IELEC,RHO)
            VPAS=-0.5D0*ALPHA/(RIP**2+D2)**2
            VPOL=MAX(VCO,VPAS)
            RVPOL(I)=VPOL*RAD(I)
            RV(I)=RV(I)+RVPOL(I)
            RVI(I)=RV(I)
          ENDDO
          VCOUT=MAX(VCPOL(IELEC,RHORMT),-0.5D0*ALPHA/(RMT**2+D2)**2)
        ENDIF
!
        IF(NPOL.LT.NDIM-5) THEN
          DO I=NPOL+1,NDIM-5
            IF(MUFIN.EQ.1.AND.IMODE.EQ.0) THEN
              RAD(I)=RAD(I-1)+0.05D0
            ELSE
              RAD(I)=1.25D0*RAD(I-1)
            ENDIF
            RADI(I)=RAD(I)
            VPOL=-0.5D0*ALPHA/(RAD(I)**2+D2)**2
            IF(IMODE.EQ.0) THEN
              IF(VPOL.GT.VCOUT) IMODE=1
              VPOL=MAX(VCOUT,VPOL)
            ENDIF
            RVPOL(I)=VPOL*RAD(I)
            RV(I)=ZINF+RVPOL(I)
            RVI(I)=RV(I)
            DEN(I)=0.0D0
            RVEX(I)=0.0D0
            RW(I)=0.0D0
            NP=I
            NPI=I
            IF(ABS(VPOL).LT.1.0D-6*MAX(E,1.0D1*ABS(ZINF)/RAD(I)).AND.RAD(I).GT.50.0D0) GO TO 34
          ENDDO
        ENDIF
        STOP 'Not enough memory space 4.'
   34   CONTINUE
        IF(NP.LT.NDIM-4) THEN
          I=NP+1
          RAD(I)=RAD(I-1)
          RADI(I)=RADI(I-1)
          RV(I)=ZINF
          DEN(I)=0.0D0
          RVI(I)=RV(I)
          RVEX(I)=0.0D0
          RW(I)=0.0D0
          RVPOL(I)=0.0D0
          DO I=NP+2,NDIM
            RAD(I)=2.0D0*RAD(I-1)
            RADI(I)=2.0D0*RADI(I-1)
            RV(I)=RV(I-1)
            RVI(I)=RV(I)
            DEN(I)=0.0D0
            RVEX(I)=0.0D0
            RW(I)=0.0D0
            RVPOL(I)=0.0D0
            NP=I
            NPI=I
            IF(RAD(I).GT.1.0D4) GO TO 35
          ENDDO
        ELSE
          STOP 'Not enough memory space 5.'
        ENDIF
   35   CONTINUE
      ENDIF
!
!  ****  At high energies, we compute the DCS for scattering by the bare
!  (finite) nucleus and multiply it by a pre-evaluated screening factor.
!
  100 CONTINUE
      IF(IHEF0.EQ.1) THEN
        IAB=0
        WRITE(*,1800)
 1800   FORMAT(1X,'#',/1X,'# WARNING: High-energy factorization with free-atom DF screening.',/1X,'#',10X,'Absorption, polarization and exchange corrections are',/1X,'#',10X,'switched off.',/1X,'#',10X,'Phase shifts are calculated for the bare nucleus.'/1X,'#',10X,'Scattering amplitudes are not evaluated.')
!  ****  Read screening function from data files.
        JT=IZ
        J1=JT-10*(JT/10)
        JT=(JT-J1)/10
        J2=JT-10*(JT/10)
        JT=(JT-J2)/10
        J3=JT-10*(JT/10)
        LIT1=LIT10(J1+1)
        LIT2=LIT10(J2+1)
        LIT3=LIT10(J3+1)
        SCFILE='elsepa_data/'//'z_'//LIT3//LIT2//LIT1//'.dfs'
        OPEN(UNIT=99,FILE=SCFILE,STATUS='OLD',ERR=4)
        READ(99,'(1X,A1)') NULL
        READ(99,'(1X,A1)') NULL
        READ(99,'(1X,A1)') NULL
        NQS=0
        DO I=1,NGT
          READ(99,*,END=4) Q2T(I),FQ(I)
          NQS=I
        ENDDO
    4   CONTINUE
        CLOSE(UNIT=99)
        IF(NQS.EQ.0) THEN
          WRITE(*,*) 'ELSEPA: I/O error. SCFILE does not exist.'
          GO TO 5
        ENDIF
!        (The number of atomic electrons is set to zero)
        CALL EFIELD(IZ,NELEC,0,MNUCL,MELEC,IW)
        DO I=1,NPOT
          RAD(I)=R(I)
          RV(I)=DFLOAT(IELEC)*RVN(I)
        ENDDO
        RV(NPOT)=DFLOAT(IELEC)*DFLOAT(IZ)
        NP=NPOT
!  ****  A constant tail is appended to the potential table...
        IF(NP.LT.NDIM-4) THEN
          I=NP+1
          RAD(I)=RAD(I-1)
          RV(I)=RV(I-1)
          DO I=NP+2,NDIM
            RAD(I)=2.0D0*RAD(I-1)
            RV(I)=RV(I-1)
            NP=I
            IF(RAD(I).GT.1.0D4) GO TO 5
          ENDDO
        ELSE
          STOP 'Not enough memory space 6'
        ENDIF
      ENDIF
    5 CONTINUE
!
!  ************  Partial-wave analysis.
!
      NDELTA=25000
      IF(IAB.EQ.0) THEN
        IF(IHEF0.EQ.1) THEN
          CALL DPWA0(EV,NDELTA,1)
        ELSE
          IF(EV.LT.1.001D3) THEN
            ISCH=1
          ELSE
            ISCH=2
          ENDIF
          IF(MCPOL.GT.0) ISCH=1
          CALL DPWA0(EV,NDELTA,ISCH)
        ENDIF
        TOTCS=ECS
        ABCS=0.0D0
      ELSE
        IF(EV.LT.1.001D3) THEN
          ISCH=1
        ELSE
          ISCH=2
        ENDIF
        IF(MCPOL.GT.0.OR.MUFIN.EQ.1) ISCH=1
        CALL DPWAI0(EV,TOTCS,ABCS,NDELTA,ISCH)
      ENDIF
!
!  ************  DCS table.
!
      IF(IHEF0.EQ.1) THEN
        NTABT=NTAB
        DO I=1,NTAB
!  ****  Screening correction and check for numerical artifacts.
          Q2=4.0D0*RK2*XT(I)
          IF(Q2.LE.Q2T(NQS)) THEN
            CALL FINDI(Q2T,Q2,NQS,J)
            F=FQ(J)+(FQ(J+1)-FQ(J))*(Q2-Q2T(J))/(Q2T(J+1)-Q2T(J))
          ELSE
            F=1.0D0
          ENDIF
          IF(TH(I).GT.1.0D0.AND.ERROR(I).LT.1.0D-2) THEN
            DCST(I)=DCST(I)*(Q2*F/(1.0D0+Q2))**2
          ELSE IF(TH(I).LE.10.0D0) THEN
            THRAD=TH(I)*PI/180.0D0
            RMR=DPWAC(THRAD)
            DCST(I)=RUTHC*F*F*RMR/(1.0D0+Q2)**2
            ERROR(I)=1.0D-5
          ELSE
            IF(ERROR(I).GT.1.0D-2) THEN
              NTABT=I-1
              GO TO 6
            ENDIF
          ENDIF
          SPOL(I)=0.0D0
          TPOL(I)=0.0d0
          UPOL(I)=0.0d0
        ENDDO
    6   CONTINUE
        IF(NTABT.LT.NTAB) THEN
          DO I=NTABT,NTAB
            DCST(I)=1.0D-45
            ERROR(I)=1.0D0
          ENDDO
        ENDIF
      ENDIF
!
!  ****  Small-angle DCS for ions.
!
      IMATCH=1
      IF(IZ.NE.NELEC) THEN
        DO I=126,NTAB
          IF(ERROR(I).LT.5.0D-4) THEN
            IMATCH=I
            GO TO 7
          ENDIF
        ENDDO
    7   CONTINUE
        WRITE(*,2112) TH(IMATCH)
 2112   FORMAT(1X,'#',/1X,'# WARNING: Partial-wave scattering amplitudes and DCSs are calculated',/1X,'#',10X,'only for THETA .gt.',1P,E10.3,' deg')
      ENDIF
!
!  ****  Integrated cross sections.
!
      IF(DCST(1).LT.1.0D-35) WRITE(*,2012) DMAX1(0.5D0,TH(IMATCH))
 2012 FORMAT(1X,'#',/1X,'# WARNING: DCSs are integrated only for',' THETA .gt.',1P,E10.3,' deg')
      ECS0=FOURPI*RMOM(XT,DCST,NTAB,0)
      ECS1=FOURPI*RMOM(XT,DCST,NTAB,1)
      ECS2=FOURPI*RMOM(XT,DCST,NTAB,2)
      ECS=ECS0
      TCS1=2.0D0*ECS1
      TCS2=6.0D0*(ECS1-ECS2)
      ECUT=MAX(20.0D3*IZ,2.0D6)
!
!      DO I=1,NTAB
!        WRITE(IW,2018) TH(I),XT(I),DCST(I),DCST(I)/A0B2,SPOL(I),TPOL(I),UPOL(I),ERROR(I)
!      ENDDO
 2018 FORMAT(1X,1P,E10.3,E13.5,5E13.5,E9.1)
!
!  ****  Scattering amplitudes.
!
      RETURN
      END
!  *********************************************************************
!                       SUBROUTINE EFIELD
!  *********************************************************************
      SUBROUTINE EFIELD(IZ,NELEC,N,MNUCL,MELEC,IW)
!
!     Electrostatic field of atoms and ions.
!
!  Input parameters:
!     IZ ....... atomic number (INTEGER).
!     NELEC .... number of electrons (INTEGER, .LE.IZ).
!     N ........ number of electrons used to determine the scattering
!                amplitudes (when the high energy factorization is used,
!                N is set equal to zero).
!     MNUCL .... nuclear density model (INTEGER).
!                 1 --> point nucleus,
!                 2 --> uniform distribution,
!                 3 --> Fermi distribution,
!                 4 --> Helm's uniform-uniform distribution.
!     MELEC .... electron density model (INTEGER).
!                 1 --> TFM analytical density,
!                 2 --> TFD analytical density,
!                 3 --> DHFS analytical density,
!                 4 --> DF density from pre-evaluated files,
!                 5 --> density read from file 'density.usr'.
!    IW ........ output unit (to be defined in the main program).
!
!  Output (through common block /CFIELD/):
!     R(I) ..... radial grid points. R(1)=0.0D0.
!     RVN(I) ... nuclear potential times R.
!     DEN(I) ... radial electron density, i.e. the electron density
!                multiplied by 4*PI*R**2.
!     RVST(I) ... atomic electrostatic potential (nuclear+electronic)
!                times R.
!     NPOT ..... number of grid points where the potential function is
!                tabulated. For I.GT.NPOT, RVST(I) is set equal to
!                RVST(NPOT).
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      character (len=1) :: LIT10(10),LIT1,LIT2,LIT3
      character (len=34) :: ELFILE,NULL
      character (len=2) :: LSYMBL(103)
!
!  ****  Set IWR=1 to print the electrostatic potential on a file.
      PARAMETER (IWR=0)
!
      PARAMETER (A0B=5.291772083D-9)  ! Bohr radius (cm)
      PARAMETER (HREV=27.2113834D0)  ! Hartree energy (eV)
      PARAMETER (F2BOHR=1.0D-13/A0B)
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (PI=3.1415926535897932D0, FOURPI=4.0D0*PI)
!
      PARAMETER (NDIM=1000,NDIN=NDIM-4)
      DIMENSION DIFR(NDIM),DENN(NDIM),RVE(NDIM),ELAW(103)
      COMMON/CFIELD/R(NDIM),RVN(NDIM),DEN(NDIM),RVST(NDIM),NPOT
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/STORE/AUX(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)
!
      DATA LIT10/'0','1','2','3','4','5','6','7','8','9'/
!
      DATA LSYMBL       /' H','He','Li','Be',' B',' C',' N',' O',' F','Ne','Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'/
      DATA ELAW     /1.007900D0,4.002600D0,6.941000D0,9.012200D0,1.081100D1,1.201070D1,1.400670D1,1.599940D1,1.899840D1,2.017970D1,2.298980D1,2.430500D1,2.698150D1,2.808550D1,3.097380D1,3.206600D1,3.545270D1,3.994800D1,3.909830D1,4.007800D1,4.495590D1,4.786700D1,5.094150D1,5.199610D1,5.493800D1,5.584500D1,5.893320D1,5.869340D1,6.354600D1,6.539000D1,6.972300D1,7.261000D1,7.492160D1,7.896000D1,7.990400D1,8.380000D1,8.546780D1,8.762000D1,8.890590D1,9.122400D1,9.290640D1,9.594000D1,9.890630D1,1.010700D2,1.029055D2,1.064200D2,1.078682D2,1.124110D2,1.148180D2,1.187100D2,1.217600D2,1.276000D2,1.269045D2,1.312900D2,1.329054D2,1.373270D2,1.389055D2,1.401160D2,1.409076D2,1.442400D2,1.449127D2,1.503600D2,1.519640D2,1.572500D2,1.589253D2,1.625000D2,1.649303D2,1.672600D2,1.689342D2,1.730400D2,1.749670D2,1.784900D2,1.809479D2,1.838400D2,1.862070D2,1.902300D2,1.922170D2,1.950780D2,1.969666D2,2.005900D2,2.043833D2,2.072000D2,2.089804D2,2.089824D2,2.099871D2,2.220176D2,2.230197D2,2.260254D2,2.270277D2,2.320381D2,2.310359D2,2.380289D2,2.370482D2,2.440642D2,2.430614D2,2.470000D2,2.470000D2,2.510000D2,2.520000D2,2.570000D2,2.580000D2,2.590000D2,2.620000D2/
!
      IF(IZ.LE.0) STOP 'Negative atomic number'
      IF(IZ.GT.103) STOP 'Atomic number larger than 103'
      Z=DBLE(IZ)
      AW=ELAW(IZ)
      IF(IZ.GT.0) THEN
        WRITE(IW,1001) LSYMBL(IZ),IZ,AW
 1001   FORMAT(1X,'#',/1X,'# Element: ',A2,',  Z = ',I3,',  atomic weight =',1P,E12.5,' g/mol')
      ELSE
        WRITE(IW,'(1X,''#'',/1X,''# Only polarization potential'')')
      ENDIF
!
!  ************  Nuclear electrostatic potential (times R).
!
      IF(MNUCL.EQ.1.OR.IZ.EQ.0) THEN
        WRITE(IW,1002)
 1002   FORMAT(1X,'#',/1X,'# Nuclear model: point charge')
        CALL SGRID(R,DIFR,1.0D-6,100.0D0,0.5D0*DBLE(NDIN),NDIN)
        DO I=1,NDIN
          RVN(I)=Z
          DENN(I)=0.0D0
        ENDDO
        GO TO 10
      ELSE IF (MNUCL.EQ.2) THEN
!  ****  The factor F2BOHR=1.889726D-5 transforms from fm to Bohr.
        R1=1.07D0*F2BOHR*AW**0.3333333333333333D0
        R2=2.00D0*F2BOHR
        R1=R1*DSQRT((1.0D0+2.5D0*(R2/R1)**2)/(1.0D0+0.75D0*(R2/R1)**2))
        R0=R1/10.0D0
        CALL SGRID(R,DIFR,R0,100.0D0,0.5D0*DBLE(NDIN),NDIN)
        WRITE(IW,1004) R1*A0B
 1004   FORMAT(1X,'#',/1X,'# Nuclear model: uniform spherical distribution',/1X,'#',16X,'Nuclear radius =',1P,E12.5,' cm')
        DO I=1,NDIN
          X=R(I)
          IF(X.LT.R1) THEN
            RVN(I)=Z*(1.5D0-0.5D0*(X/R1)**2)*(X/R1)
            DENN(I)=Z/(FOURPI*R1**3/3.0D0)
          ELSE
            RVN(I)=Z
            DENN(I)=0.0D0
          ENDIF
        ENDDO
        GO TO 10
      ELSE IF (MNUCL.EQ.3) THEN
        R1=1.07D0*F2BOHR*AW**0.3333333333333333D0
        R2=0.546D0*F2BOHR
        R0=R1/10.0D0
        CALL SGRID(R,DIFR,R0,100.0D0,0.5D0*DBLE(NDIN),NDIN)
        WRITE(IW,1005) R1*A0B,R2*A0B
 1005   FORMAT(1X,'#',/1X,'# Nuclear model: Fermi distribution',/1X,'#',16X,'Average radius =',1P,E12.5, ' cm',/1X,'#',16X,'Skin thickness =',E12.5,' cm')
!  ****  The array DENN contains the nuclear charge density.
!        (unnormalized).
        DO I=1,NDIN
          X=R(I)
          XX=EXP((R1-X)/R2)
          DENN(I)=XX/(1.0D0+XX)
          RVN(I)=DENN(I)*X**2*DIFR(I)
        ENDDO
      ELSE
        RNUC=1.070D0*F2BOHR*AW**0.3333333333333333D0
        R1=0.96219D0*RNUC+0.435D0*F2BOHR
        R2=2.0D0*F2BOHR
        WRITE(IW,1105) R1*A0B,R2*A0B
 1105   FORMAT(1X,'#',/1X,'# Nuclear model: Helm''s Uu distribution',/1X,'#',16X,'  Inner radius =',1P,E12.5, ' cm',/1X,'#',16X,'Skin thickness =',E12.5,' cm')
        IF(R2.GT.R1) THEN
          STORED=R1
          R1=R2
          R2=STORED
        ENDIF
        R0=R1/10.0D0
        CALL SGRID(R,DIFR,R0,100.0D0,0.5D0*DBLE(NDIN),NDIN)
!
        DO I=1,NDIN
          RR=R(I)
          IF(RR.LE.R1-R2) THEN
            V=1.0D0
          ELSE IF(RR.GT.R1+R2) THEN
            V=0.0D0
          ELSE
            T=RR*RR+R1*R1-R2*R2
            V1=(T+4.0D0*RR*R1)*(T-2.0D0*RR*R1)**2
            T=RR*RR+R2*R2-R1*R1
            V2=(T+4.0D0*RR*R2)*(T-2.0D0*RR*R2)**2
            V=(V1+V2)/(32.0D0*(R2*RR)**3)
          ENDIF
          DENN(I)=V
          RVN(I)=DENN(I)*RR**2*DIFR(I)
        ENDDO
      ENDIF
!
      CALL SLAG6(1.0D0,RVN,RVN,NDIN)
      NDIN1=NDIN+1
      DO I=1,NDIN
        K=NDIN1-I
        AUX(I)=DENN(K)*R(K)*DIFR(K)
      ENDDO
      CALL SLAG6(1.0D0,AUX,AUX,NDIN)
      FNORM=Z/RVN(NDIN)
      DO I=1,NDIN
        RVN(I)=FNORM*(RVN(I)+AUX(NDIN1-I)*R(I))
        DENN(I)=FNORM*DENN(I)/FOURPI
        IF(DENN(I).LT.1.0D-35) DENN(I)=0.0D0
      ENDDO
   10 CONTINUE
!
!  ************  Electronic potential (times R).
!
      WRITE(IW,1006) NELEC
 1006 FORMAT(1X,'#',/1X,'# Number of electrons =',I3)
      IF(N.LT.0) STOP 'Negative number of electrons'
      IF(N.GT.IZ) STOP 'Negative ion'
      IF(N.EQ.0) THEN
        DO I=1,NDIN
          DEN(I)=0.0D0
          RVE(I)=0.0D0
        ENDDO
        GO TO 2
      ENDIF
!
      IF(MELEC.LT.4.OR.MELEC.GT.5) THEN
!  ****  Analytical electron density models.
        IF(MELEC.EQ.1) THEN
          CALL TFM(IZ,A1,A2,A3,AL1,AL2,AL3)
          WRITE(IW,1007)
 1007     FORMAT(1X,'#',/1X,'# Electron density: analytical TFM model')
        ELSE IF(MELEC.EQ.2) THEN
          CALL TFD(IZ,A1,A2,A3,AL1,AL2,AL3)
          WRITE(IW,1008)
 1008     FORMAT(1X,'#',/1X,'# Electron density: analytical TFD model')
        ELSE
          CALL DHFS(IZ,A1,A2,A3,AL1,AL2,AL3)
          WRITE(IW,1009)
 1009     FORMAT(1X,'#',/1X,'# Electron density: analytical DHFS model')
        ENDIF
        WRITE(IW,1010) A1,AL1,A2,AL2,A3,AL3
 1010   FORMAT(1X,'#',19X,'A1 = ',1P,D12.5,' ,   ALPHA1 =',D12.5,/1X,'#',19X,'A2 = ',D12.5,' ,   ALPHA2 =',D12.5,/1X,'#',19X,'A3 = ',D12.5,' ,   ALPHA3 =',D12.5)
        XN=DBLE(N)
        DO I=1,NDIN
          DEN(I)=(A1*AL1*AL1*DEXP(-AL1*R(I))+A2*AL2*AL2*DEXP(-AL2*R(I))+A3*AL3*AL3*DEXP(-AL3*R(I)))*XN
        ENDDO
      ELSE
!  ****  Electron density read from a file.
        NE=0
        IF(MELEC.EQ.4) THEN
          JT=IZ
          J1=JT-10*(JT/10)
          JT=(JT-J1)/10
          J2=JT-10*(JT/10)
          JT=(JT-J2)/10
          J3=JT-10*(JT/10)
          LIT1=LIT10(J1+1)
          LIT2=LIT10(J2+1)
          LIT3=LIT10(J3+1)
          ELFILE='elsepa_data/'//'z_'//LIT3//LIT2//LIT1//'.den'
        ELSE
          ELFILE='density.usr '
        ENDIF
!
        WRITE(IW,1011) ELFILE
 1011   FORMAT(1X,'#',/1X,'# Electron density: Read from file ',A12)
        OPEN(99,FILE=ELFILE,STATUS='OLD',ERR=1)
        READ(99,'(A12)') NULL
        READ(99,'(A12)') NULL
        READ(99,'(A12)') NULL
        DO I=1,NDIN
          READ(99,*,END=1) AUX(I),RVE(I)
          NE=I
          RVE(I)=DLOG(RVE(I))
        ENDDO
        STOP 'EFIELD: File is too large.'
    1   CONTINUE
        IF(NE.EQ.0) STOP 'I/O error in EFIELD'
        CLOSE(99)
        WRITE(IW,1012) NE
 1012   FORMAT(1X,'#',19X,'Number of data points = ',I4)
        IF(NE.LT.4) STOP 'SPLINE needs more than 4 points'
!  ****  ... and interpolated (lin-log cubic spline).
        CALL SPLINE(AUX,RVE,A,B,C,D,0.0D0,0.0D0,NE)
        B(NE)=(RVE(NE)-RVE(NE-1))/(AUX(NE)-AUX(NE-1))
        A(NE)=RVE(NE-1)-B(NE)*AUX(NE-1)
        C(NE)=0.0D0
        D(NE)=0.0D0
        DO I=1,NDIN
          X=R(I)
          IF(X.GT.AUX(NE)) THEN
            DEN(I)=0.0D0
          ELSE
            CALL FINDI(AUX,X,NE,J)
            DEN(I)=DEXP(A(J)+X*(B(J)+X*(C(J)+X*D(J))))*X*FOURPI
          ENDIF
        ENDDO
      ENDIF
!  ****  Calculation of the electrostatic potential.
      DO I=1,NDIN
        RVE(I)=DEN(I)*R(I)*DIFR(I)
      ENDDO
      CALL SLAG6(1.0D0,RVE,RVE,NDIN)
      NDIN1=NDIN+1
      DO I=1,NDIN
        K=NDIN1-I
        AUX(I)=DEN(K)*DIFR(K)
      ENDDO
      CALL SLAG6(1.0D0,AUX,AUX,NDIN)
      WRITE(IW,1013) RVE(NDIN)
 1013 FORMAT(1X,'#',19X,'Volume integral =',1P,E12.5)
      FNORM=DBLE(N)/RVE(NDIN)
      DO I=1,NDIN
        RVE(I)=FNORM*(RVE(I)+AUX(NDIN1-I)*R(I))
        DEN(I)=FNORM*DEN(I)*R(I)
      ENDDO
!
    2 CONTINUE
      ZINF=DBLE(IZ-N)
      DO I=1,NDIN
        RVST(I)=RVN(I)-RVE(I)
      ENDDO
      NPOT=NDIN
      DO I=NDIN,6,-1
        IF((ABS(RVST(I)-ZINF).LT.5.0D-12).AND.(ABS(DEN(I)).LT.5.0D-12)) THEN
          RVST(I)=ZINF
          NPOT=I
        ELSE
          GO TO 3
        ENDIF
      ENDDO
    3 CONTINUE
!
      RETURN
      END
!  *********************************************************************
!                       SUBROUTINE TFM
!  *********************************************************************
      SUBROUTINE TFM(IZ,A1,A2,A3,AL1,AL2,AL3)
!
!     Parameters in Moliere's analytical approximation (three Yukawa
!  terms) to the Thomas-Fermi field.
!     Ref.: G. Moliere, Z. Naturforsch. 2a (1947) 133.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT Double Precision (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      Z=DBLE(IZ)
      RTF=0.88534D0/Z**0.33333333333D0
      AL1=6.0D0/RTF
      AL2=1.2D0/RTF
      AL3=0.3D0/RTF
      A1=0.10D0
      A2=0.55D0
      A3=0.35D0
      RETURN
      END
!  *********************************************************************
!                       SUBROUTINE TFD
!  *********************************************************************
      SUBROUTINE TFD(IZ,A1,A2,A3,AL1,AL2,AL3)
!
!     Parameters in the analytical approximation (three Yukawa terms)
!  for the Thomas-Fermi-Dirac field.
!     Ref.: R.A. Bonham and T.G. Strand, J. Chem. Phys. 39 (1963) 2200.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT Double Precision (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      DIMENSION AA1(5),AA2(5),AA3(5),AAL1(5),AAL2(5),AAL3(5)
      DATA AA1/1.26671D-2,-2.61047D-2,2.14184D-2,-2.35686D-3,2.10672D-5/
      DATA AA2/5.80612D-2,2.93077D-2,8.57135D-2,-2.23342D-2,1.64675D-3/
      DATA AA3/9.27968D-1,-1.64643D-3,-1.07685D-1,2.47998D-2,-1.67822D-3/
      DATA AAL1/1.64564D2,-1.52192D2,6.23879D1,-1.15005D1,8.08424D-1/
      DATA AAL2/1.13060D1,-6.31902D0,2.26025D0,-3.70738D-1,2.61151D-2/
      DATA AAL3/1.48219D0,-5.57601D-2,1.64387D-2,-4.39703D-3,9.97225D-4/
!
      IF(IZ.LE.0) THEN
        WRITE(6,100)
  100   FORMAT(5X,'*** TFD: Negative atomic number. STOP.')
        STOP 'TFD: Negative atomic number.'
      ENDIF
!
      X=DLOG(DBLE(IZ))
      A1=AA1(1)+X*(AA1(2)+X*(AA1(3)+X*(AA1(4)+X*AA1(5))))
      A2=AA2(1)+X*(AA2(2)+X*(AA2(3)+X*(AA2(4)+X*AA2(5))))
      A3=AA3(1)+X*(AA3(2)+X*(AA3(3)+X*(AA3(4)+X*AA3(5))))
      AL1=AAL1(1)+X*(AAL1(2)+X*(AAL1(3)+X*(AAL1(4)+X*AAL1(5))))
      AL2=AAL2(1)+X*(AAL2(2)+X*(AAL2(3)+X*(AAL2(4)+X*AAL2(5))))
      AL3=AAL3(1)+X*(AAL3(2)+X*(AAL3(3)+X*(AAL3(4)+X*AAL3(5))))
      A3=1.0D0-A1-A2
!
      RETURN
      END
!  *********************************************************************
!                       SUBROUTINE DHFS
!  *********************************************************************
      SUBROUTINE DHFS(IZ,A1,A2,A3,AL1,AL2,AL3)
!
!     DHFS analytical screening function parameters for free neutral
!  atoms. The input argument is the atomic number.
!
!     Ref.: F. Salvat et al., Phys. Rev. A36 (1987) 467-474.
!     Elements from Z=93 to 103 added in march 1992.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      DIMENSION B1(103),B2(103),BL1(103),BL2(103),BL3(103)
      DATA B1/-7.05665D-6,-2.25920D-1,6.04537D-1,3.27766D-1,2.32684D-1,1.53676D-1,9.95750D-2,6.25130D-2,3.68040D-2,1.88410D-2,7.44440D-1,6.42349D-1,6.00152D-1,5.15971D-1,4.38675D-1,5.45871D-1,7.24889D-1,2.19124D+0,4.85607D-2,5.80017D-1,5.54340D-1,1.11950D-2,3.18350D-2,1.07503D-1,4.97556D-2,5.11841D-2,5.00039D-2,4.73509D-2,7.70967D-2,4.00041D-2,1.08344D-1,6.09767D-2,2.11561D-2,4.83575D-1,4.50364D-1,4.19036D-1,1.73438D-1,3.35694D-2,6.88939D-2,1.17552D-1,2.55689D-1,2.69313D-1,2.20138D-1,2.75057D-1,2.71053D-1,2.78363D-1,2.56210D-1,2.27100D-1,2.49215D-1,2.15313D-1,1.80560D-1,1.30772D-1,5.88293D-2,4.45145D-1,2.70796D-1,1.72814D-1,1.94726D-1,1.91338D-1,1.86776D-1,1.66461D-1,1.62350D-1,1.58016D-1,1.53759D-1,1.58729D-1,1.45327D-1,1.41260D-1,1.37360D-1,1.33614D-1,1.29853D-1,1.26659D-1,1.28806D-1,1.30256D-1,1.38420D-1,1.50030D-1,1.60803D-1,1.72164D-1,1.83411D-1,2.23043D-1,2.28909D-1,2.09753D-1,2.70821D-1,2.37958D-1,2.28771D-1,1.94059D-1,1.49995D-1,9.55262D-2,3.19155D-1,2.40406D-1,2.26579D-1,2.17619D-1,2.41294D-1,2.44758D-1,2.46231D-1,2.55572D-1,2.53567D-1,2.43832D-1,2.41898D-1,2.44050D-1,2.40237D-1,2.34997D-1,2.32114D-1,2.27937D-1,2.29571D-1/
      DATA B2/-1.84386D+2,1.22592D+0,3.95463D-1,6.72234D-1,7.67316D-1,8.46324D-1,9.00425D-1,9.37487D-1,9.63196D-1,9.81159D-1,2.55560D-1,3.57651D-1,3.99848D-1,4.84029D-1,5.61325D-1,-5.33329D-1,-7.54809D-1,-2.2852D0,7.75935D-1,4.19983D-1,4.45660D-1,6.83176D-1,6.75303D-1,7.16172D-1,6.86632D-1,6.99533D-1,7.14201D-1,7.29404D-1,7.95083D-1,7.59034D-1,7.48941D-1,7.15671D-1,6.70932D-1,5.16425D-1,5.49636D-1,5.80964D-1,7.25336D-1,7.81581D-1,7.20203D-1,6.58088D-1,5.82051D-1,5.75262D-1,5.61797D-1,5.94338D-1,6.11921D-1,6.06653D-1,6.50520D-1,6.15496D-1,6.43990D-1,6.11497D-1,5.76688D-1,5.50366D-1,5.48174D-1,5.54855D-1,6.52415D-1,6.84485D-1,6.38429D-1,6.46684D-1,6.55810D-1,7.05677D-1,7.13311D-1,7.20978D-1,7.28385D-1,7.02414D-1,7.42619D-1,7.49352D-1,7.55797D-1,7.61947D-1,7.68005D-1,7.73365D-1,7.52781D-1,7.32428D-1,7.09596D-1,6.87141D-1,6.65932D-1,6.46849D-1,6.30598D-1,6.17575D-1,6.11402D-1,6.00426D-1,6.42829D-1,6.30789D-1,6.21959D-1,6.10455D-1,6.03147D-1,6.05994D-1,6.23324D-1,6.56665D-1,6.42246D-1,6.24013D-1,6.30394D-1,6.29816D-1,6.31596D-1,6.49005D-1,6.53604D-1,6.43738D-1,6.48850D-1,6.70318D-1,6.76319D-1,6.65571D-1,6.88406D-1,6.94394D-1,6.82014D-1/
      DATA BL1/ 4.92969D+0,5.52725D+0,2.81741D+0,4.54302D+0,5.99006D+0,8.04043D+0,1.08122D+1,1.48233D+1,2.14001D+1,3.49994D+1,4.12050D+0,4.72663D+0,5.14051D+0,5.84918D+0,6.67070D+0,6.37029D+0,6.21183D+0,5.54701D+0,3.02597D+1,6.32184D+0,6.63280D+0,9.97569D+1,4.25330D+1,1.89587D+1,3.18642D+1,3.18251D+1,3.29153D+1,3.47580D+1,2.53264D+1,4.03429D+1,2.01922D+1,2.91996D+1,6.24873D+1,8.78242D+0,9.33480D+0,9.91420D+0,1.71659D+1,5.52077D+1,3.13659D+1,2.20537D+1,1.42403D+1,1.40442D+1,1.59176D+1,1.43137D+1,1.46537D+1,1.46455D+1,1.55878D+1,1.69141D+1,1.61552D+1,1.77931D+1,1.98751D+1,2.41540D+1,3.99955D+1,1.18053D+1,1.65915D+1,2.23966D+1,2.07637D+1,2.12350D+1,2.18033D+1,2.39492D+1,2.45984D+1,2.52966D+1,2.60169D+1,2.54973D+1,2.75466D+1,2.83460D+1,2.91604D+1,2.99904D+1,3.08345D+1,3.16806D+1,3.13526D+1,3.12166D+1,3.00767D+1,2.86302D+1,2.75684D+1,2.65861D+1,2.57339D+1,2.29939D+1,2.28644D+1,2.44080D+1,2.09409D+1,2.29872D+1,2.37917D+1,2.66951D+1,3.18397D+1,4.34890D+1,2.00150D+1,2.45012D+1,2.56843D+1,2.65542D+1,2.51930D+1,2.52522D+1,2.54271D+1,2.51526D+1,2.55959D+1,2.65567D+1,2.70360D+1,2.72673D+1,2.79152D+1,2.86446D+1,2.93353D+1,3.01040D+1,3.02650D+1/


      DATA BL2/ 2.00272D+0,2.39924D+0,6.62463D-1,9.85154D-1,1.21347D+0,1.49129D+0,1.76868D+0,2.04035D+0,2.30601D+0,2.56621D+0,8.71798D-1,1.00247D+0,1.01529D+0,1.17314D+0,1.34102D+0,2.55169D+0,3.38827D+0,4.56873D+0,3.12426D+0,1.00935D+0,1.10227D+0,4.12865D+0,3.94043D+0,3.06375D+0,3.78110D+0,3.77161D+0,3.79085D+0,3.82989D+0,3.39276D+0,3.94645D+0,3.47325D+0,4.12525D+0,4.95015D+0,1.69671D+0,1.79002D+0,1.88354D+0,3.11025D+0,4.28418D+0,4.24121D+0,4.03254D+0,2.97020D+0,2.86107D+0,3.36719D+0,2.73701D+0,2.71828D+0,2.61549D+0,2.74124D+0,3.08408D+0,2.88189D+0,3.29372D+0,3.80921D+0,4.61191D+0,5.91318D+0,1.79673D+0,2.69645D+0,3.45951D+0,3.46574D+0,3.48193D+0,3.50982D+0,3.51987D+0,3.55603D+0,3.59628D+0,3.63834D+0,3.73639D+0,3.72882D+0,3.77625D+0,3.82444D+0,3.87344D+0,3.92327D+0,3.97271D+0,4.09040D+0,4.20492D+0,4.24918D+0,4.24261D+0,4.23412D+0,4.19992D+0,4.14615D+0,3.73461D+0,3.69138D+0,3.96429D+0,3.24563D+0,3.62172D+0,3.77959D+0,4.25824D+0,4.92848D+0,5.85205D+0,2.90906D+0,3.55241D+0,3.79223D+0,4.00437D+0,3.67795D+0,3.63966D+0,3.61328D+0,3.43021D+0,3.43474D+0,3.59089D+0,3.59411D+0,3.48061D+0,3.50331D+0,3.61870D+0,3.55697D+0,3.58685D+0,3.64085D+0/
      DATA BL3/ 1.99732D+0,1.00000D+0,1.00000D+0,1.00000D+0,1.00000D+0,1.00000D+0,1.00000D+0,1.00000D+0,1.00000D+0,1.00000D+0,1.00000D+0,1.00000D+0,1.00000D+0,1.00000D+0,1.00000D+0,1.67534D+0,1.85964D+0,2.04455D+0,7.32637D-1,1.00000D+0,1.00000D+0,1.00896D+0,1.05333D+0,1.00137D+0,1.12787D+0,1.16064D+0,1.19152D+0,1.22089D+0,1.14261D+0,1.27594D+0,1.00643D+0,1.18447D+0,1.35819D+0,1.00000D+0,1.00000D+0,1.00000D+0,7.17673D-1,8.57842D-1,9.47152D-1,1.01806D+0,1.01699D+0,1.05906D+0,1.15477D+0,1.10923D+0,1.12336D+0,1.43183D+0,1.14079D+0,1.26189D+0,9.94156D-1,1.14781D+0,1.28288D+0,1.41954D+0,1.54707D+0,1.00000D+0,6.81361D-1,8.07311D-1,8.91057D-1,9.01112D-1,9.10636D-1,8.48620D-1,8.56929D-1,8.65025D-1,8.73083D-1,9.54998D-1,8.88981D-1,8.96917D-1,9.04803D-1,9.12768D-1,9.20306D-1,9.28838D-1,1.00717D+0,1.09456D+0,1.16966D+0,1.23403D+0,1.29699D+0,1.35350D+0,1.40374D+0,1.44284D+0,1.48856D+0,1.53432D+0,1.11214D+0,1.23735D+0,1.25338D+0,1.35772D+0,1.46828D+0,1.57359D+0,7.20714D-1,8.37599D-1,9.33468D-1,1.02385D+0,9.69895D-1,9.82474D-1,9.92527D-1,9.32751D-1,9.41671D-1,1.01827D+0,1.02554D+0,9.66447D-1,9.74347D-1,1.04137D+0,9.90568D-1,9.98878D-1,1.04473D+0/
!
      IIZ=IABS(IZ)
      IF(IIZ.GT.103) IIZ=103
      IF(IIZ.EQ.0) IIZ=1
      A1=B1(IIZ)
      A2=B2(IIZ)
      A3=1.0D0-(A1+A2)
      IF(ABS(A3).LT.1.0D-15) A3=0.0D0
      AL1=BL1(IIZ)
      AL2=BL2(IIZ)
      AL3=BL3(IIZ)
      RETURN
      END
!  *********************************************************************
!                       SUBROUTINE MOTTCS
!  *********************************************************************
      SUBROUTINE MOTTSC(IELEC,IZ,EV,IW)
!
!     Mott cross section for elastic scattering of high-energy electrons
!  and positrons by unscreened point nuclei.
!
!  Input parameters:
!    IELEC ..... electron-positron flag;
!                =-1 for electrons,
!                =+1 for positrons.
!    IZ ........ atomic number of the target atom.
!    EV ........ projectile's kinetic energy (in eV).
!    IW ........ output unit (to be defined in the main program).
!
!  The Mott DCS and spin polarization function are printed on unit IW.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX (KIND=RegReal_K) (C),INTEGER (KIND=RegInt_K) (I-N)
!
      PARAMETER (SL=137.03599976D0)  ! Speed of light (1/alpha)
      PARAMETER (A0B=5.291772083D-9)  ! Bohr radius (cm)
      PARAMETER (HREV=27.2113834D0)  ! Hartree energy (eV)
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (F2BOHR=1.0D-13/A0B)
      PARAMETER (PI=3.1415926535897932D0,FOURPI=4.0D0*PI)
!
      PARAMETER (NGT=650)
      COMMON/DCSTAB/ECS,TCS1,TCS2,TH(NGT),XT(NGT),DCST(NGT),SPOL(NGT),TPOL(NGT),UPOL(NGT),ERROR(NGT),NTAB
!
      PARAMETER (NPC=1500)
      COMMON/CRMORU/CFM(NPC),CGM(NPC),DPC(NPC),DMC(NPC),CF,CG,RUTHC,WATSC,RK2,ERRF,ERRG,NPC1
!
      WRITE(IW,1000)
 1000 FORMAT(1X,'#',/1X,'# Subroutine MOTTSC. Elastic scattering of ','electrons and positrons',/1X,'#',20X,'by unscreened Coulomb fields')
      IF(IELEC.EQ.-1) THEN
        WRITE(IW,1100)
 1100   FORMAT(1X,'#',/1X,'# Projectile: electron')
      ELSE
        WRITE(IW,1200)
 1200   FORMAT(1X,'#',/1X,'# Projectile: positron')
      ENDIF
      E=EV/HREV
      WRITE(IW,1300) EV,E
 1300 FORMAT(1X,'# Kinetic energy =',1P,E12.5,' eV =',E12.5,' a.u.')
!
      WRITE(*,'(1X,''#'',/1X,''#  ***  WARNING: Mott '',''scattering (point unscreend nucleus).'')')
!
      IF(IZ.LE.0) STOP 'Negative atomic number.'
      WRITE(IW,1001) IZ
 1001 FORMAT(1X,'#',/1X,'# Z = ',I3)
!
      Z=DFLOAT(IZ*IELEC)
      CALL DPWAC0(Z,EV)
!
      TH(1)=0.0D0
      TH(2)=1.0D-4
      I=2
   10 CONTINUE
      I=I+1
      IF(TH(I-1).LT.0.9999D-3) THEN
        TH(I)=TH(I-1)+2.5D-5
      ELSE IF(TH(I-1).LT.0.9999D-2) THEN
        TH(I)=TH(I-1)+2.5D-4
      ELSE IF(TH(I-1).LT.0.9999D-1) THEN
        TH(I)=TH(I-1)+2.5D-3
      ELSE IF(TH(I-1).LT.0.9999D+0) THEN
        TH(I)=TH(I-1)+2.5D-2
      ELSE IF(TH(I-1).LT.0.9999D+1) THEN
        TH(I)=TH(I-1)+1.0D-1
      ELSE IF(TH(I-1).LT.2.4999D+1) THEN
        TH(I)=TH(I-1)+2.5D-1
      ELSE
        TH(I)=TH(I-1)+5.0D-1
      ENDIF
      IF(TH(I).LT.180.0D0) GO TO 10
      NTAB=I
!
      NTABT=NTAB
      DO I=1,NTAB
        THR=TH(I)*PI/180.0D0
        XT(I)=(1.0D0-COS(THR))/2.0D0
        IF(TH(I).GT.1.0D-5) THEN
          Q2=4.0D0*RK2*XT(I)
          RMR=DPWAC(THR)
          DCST(I)=RUTHC*RMR/Q2**2
!  ****  Spin polarization (Sherman) function.
          CF=CF*A0B
          CG=CG*A0B
          ACF=CDABS(CF)**2
          ACG=CDABS(CG)**2
          DCS=ACF+ACG
          IF(DCS.GT.1.0D-45) THEN
            ERR=2.0D0*(ACF*ERRF+ACG*ERRG)/DCS
          ELSE
            ERR=1.0D0
          ENDIF
          ERROR(I)=ERR
!
          CSPL1=DCMPLX(0.0D0,1.0D0)*CF*DCONJG(CG)
          CSPL2=DCMPLX(0.0D0,1.0D0)*CG*DCONJG(CF)
          TST=CDABS(CSPL1-CSPL2)/(CDABS(CSPL1)+1.0D-30)
          IF(TST.GT.1.0D-6.AND.ERROR(I).LT.1.0D-5) THEN
            SPOL(I)=dble(CSPL1-CSPL2)/DCS
            TPOL(I)=dble(ACF-ACG)/DCS
            UPOL(I)=dble(CF*DCONJG(CG)+CG*DCONJG(CF))/DCS
            IF (TH(I) .EQ. 0.d0 .OR. TH(I) .EQ. 180.d0) THEN
              SPOL(I)=0.d0
              TPOL(I)=1.d0
              UPOL(I)=0.d0
            ENDIF
            IF (SPOL(I) .NE. 0.d0 .OR. TPOL(I) .NE. 0.d0 .OR. UPOL(I) .NE. 0.d0) THEN
              FACT=1.d0/(SPOL(I)**2+TPOL(I)**2+UPOL(I)**2)
              SPOL(I)=SPOL(I)*FACT
              TPOL(I)=TPOL(I)*FACT
              UPOL(I)=UPOL(I)*FACT
            END IF
          ELSE
            SPOL(I)=0.0D0
            TPOL(I)=1.d0
            UPOL(I)=0.d0
          ENDIF
        ELSE
          DCST(I)=1.0D-45
          SPOL(I)=0.0D0
          TPOL(I)=0.0d0
          UPOL(I)=0.0d0
          ERROR(I)=1.0D0
        ENDIF
      ENDDO
!
      WRITE(IW,'(1X,''#'')')
      WRITE(IW,'(1X,''# Differential cross section'',6X,''MU=(1-COS(THETA))/2'')')
      WRITE(IW,'(1X,''#'',/1X,''#  THETA'',8X,''MU'',10X,''DCS'',10X,''DCS'',22X,''Sherman'',18X,''error''/1X,''#  (deg)'',17X,''(cm**2/sr)'',3X,''(a0**2/sr)'',9X,''S'',11X,''T'',12X,''U'',/1X,''#'',70(''-''))')
!      DO I=1,NTAB
!        WRITE(IW,2018) TH(I),XT(I),DCST(I),DCST(I)/A0B2,SPOL(I),TPOL(I),UPOL(I),ERROR(I)
!      ENDDO
 2018 FORMAT(1X,1P,E10.3,E13.5,5E13.5,2X,E8.1)
!
      ECS=1.0D35
      TCS1=1.0D35
      TCS2=1.0D35
!
      RETURN
      END
!  *********************************************************************
!                       SUBROUTINE HEBORN
!  *********************************************************************
      SUBROUTINE HEBORN(IELEC,IZ,MNUCL,EV,IW)
!
!     Mott-Born cross section for elastic scattering of high-energy
!  electrons and positrons by neutral atoms.
!
!    The DCS is obtained as the product of the Mott DCS for a point
!  nucleus, the Helm uniform-uniform nuclear form factor (with an
!  empirical Coulomb correction) and the high-energy DF screening
!  factor.
!
!  Input parameters:
!    IELEC ..... electron-positron flag;
!                =-1 for electrons,
!                =+1 for positrons.
!    IZ ........ atomic number of the target atom.
!    MNUCL ..... nuclear charge density model.
!                  1 --> point nucleus (P),
!                  2 --> uniform distribution (U),
!                3,4 --> Helm's uniform-uniform distribution (Uu).
!    EV ........ projectile's kinetic energy (in eV).
!    IW ........ output unit (to be defined in the main program).
!
!  Output (through the common block /DCSTAB/):
!     ECS ........ total cross section (cm**2).
!     TCS1 ....... 1st transport cross section (cm**2).
!     TCS2 ....... 2nd transport cross section (cm**2).
!     TH(I) ...... scattering angles (in deg)
!     XT(I) ...... values of (1-COS(TH(I)))/2.
!     DCST(I) .... differential cross section per unit solid angle at
!                  TH(I) (in cm**2/sr).
!     ERROR(I) ... relative uncertainty of the computed DCS values.
!                  Estimated from the convergence of the series.
!     NTAB ....... number of angles in the table.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX (KIND=RegReal_K) (C),INTEGER (KIND=RegInt_K) (I-N)
!
      PARAMETER (SL=137.03599976D0)  ! Speed of light (1/alpha)
      PARAMETER (A0B=5.291772083D-9)  ! Bohr radius (cm)
      PARAMETER (HREV=27.2113834D0)  ! Hartree energy (eV)
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (F2BOHR=1.0D-13/A0B)
      PARAMETER (PI=3.1415926535897932D0,FOURPI=4.0D0*PI)
!
      PARAMETER (NGT=650)
      COMMON/DCSTAB/ECS,TCS1,TCS2,TH(NGT),XT(NGT),DCST(NGT),SPOL(NGT),TPOL(NGT),UPOL(NGT),ERROR(NGT),NTAB
      COMMON/CTOTCS/TOTCS,ABCS
      COMMON/CDCSHE/Q2T(NGT),FQ(NGT),U1,U2,NQS,MOM
!
      DIMENSION ELAW(103)
!
      PARAMETER (NPC=1500)
      COMMON/CRMORU/CFM(NPC),CGM(NPC),DPC(NPC),DMC(NPC),CF,CG,RUTHC,WATSC,RK2,ERRF,ERRG,NPC1
!
      character (len=34) :: SCFILE,NULL
      character (len=1) :: LIT10(10),LIT1,LIT2,LIT3
      character (len=2) :: LSYMBL(103)
      DATA LIT10/'0','1','2','3','4','5','6','7','8','9'/
!
      DATA LSYMBL/' H','He','Li','Be',' B',' C',' N',' O',' F','Ne','Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'/
!
      DATA ELAW/1.007900D0,4.002600D0,6.941000D0,9.012200D0,1.081100D1,1.201070D1,1.400670D1,1.599940D1,1.899840D1,2.017970D1,2.298980D1,2.430500D1,2.698150D1,2.808550D1,3.097380D1,3.206600D1,3.545270D1,3.994800D1,3.909830D1,4.007800D1,4.495590D1,4.786700D1,5.094150D1,5.199610D1,5.493800D1,5.584500D1,5.893320D1,5.869340D1,6.354600D1,6.539000D1,6.972300D1,7.261000D1,7.492160D1,7.896000D1,7.990400D1,8.380000D1,8.546780D1,8.762000D1,8.890590D1,9.122400D1,9.290640D1,9.594000D1,9.890630D1,1.010700D2,1.029055D2,1.064200D2,1.078682D2,1.124110D2,1.148180D2,1.187100D2,1.217600D2,1.276000D2,1.269045D2,1.312900D2,1.329054D2,1.373270D2,1.389055D2,1.401160D2,1.409076D2,1.442400D2,1.449127D2,1.503600D2,1.519640D2,1.572500D2,1.589253D2,1.625000D2,1.649303D2,1.672600D2,1.689342D2,1.730400D2,1.749670D2,1.784900D2,1.809479D2,1.838400D2,1.862070D2,1.902300D2,1.922170D2,1.950780D2,1.969666D2,2.005900D2,2.043833D2,2.072000D2,2.089804D2,2.089824D2,2.099871D2,2.220176D2,2.230197D2,2.260254D2,2.270277D2,2.320381D2,2.310359D2,2.380289D2,2.370482D2,2.440642D2,2.430614D2,2.470000D2,2.470000D2,2.510000D2,2.520000D2,2.570000D2,2.580000D2,2.590000D2,2.620000D2/
!
!      EXTERNAL DCSHB
!
      WRITE(IW,1000)
 1000 FORMAT(1X,'#',/1X,'# Subroutine HEBORN. Elastic scattering of ','electrons and positrons',/1X,'#',20X,'by neutral atoms')
      IF(IELEC.EQ.-1) THEN
        WRITE(IW,1100)
 1100   FORMAT(1X,'#',/1X,'# Projectile: electron')
      ELSE
        WRITE(IW,1200)
 1200   FORMAT(1X,'#',/1X,'# Projectile: positron')
      ENDIF
      E=EV/HREV
      WRITE(IW,1300) EV,E
 1300 FORMAT(1X,'# Kinetic energy =',1P,E12.5,' eV =',E12.5,' a.u.')
!
      WRITE(*,'(1X,''#'',/1X,''#  ***  WARNING: High-energy Mott-Born approximation. Neutral atom.'')')
!
      IF(IZ.LE.0) STOP 'Negative atomic number.'
      IF(IZ.GT.103) STOP 'Atomic number larger than 103.'
      AW=ELAW(IZ)
      WRITE(IW,1001) LSYMBL(IZ),IZ,AW
 1001 FORMAT(1X,'#',/1X,'# Element: ',A2,',  Z = ',I3,',  atomic weight =',1P,E12.5,' g/mol')
!
      Z=DFLOAT(IZ*IELEC)
      CALL DPWAC0(Z,EV)
!
!  ****  Read screening function from data files.
      JT=IZ
      J1=JT-10*(JT/10)
      JT=(JT-J1)/10
      J2=JT-10*(JT/10)
      JT=(JT-J2)/10
      J3=JT-10*(JT/10)
      LIT1=LIT10(J1+1)
      LIT2=LIT10(J2+1)
      LIT3=LIT10(J3+1)
      SCFILE='elsepa_data/'//'z_'//LIT3//LIT2//LIT1//'.dfs'
      OPEN(UNIT=99,FILE=SCFILE,STATUS='OLD',ERR=4)
      READ(99,'(1X,A1)') NULL
      READ(99,'(1X,A1)') NULL
      READ(99,'(1X,A1)') NULL
      NQS=0
      DO I=1,NGT
        READ(99,*,END=4) Q2T(I),FQ(I)
        NQS=I
      ENDDO
    4 CONTINUE
      CLOSE(UNIT=99)
      IF(NQS.EQ.0) THEN
        WRITE(*,*) 'ELSEPA: I/O error. SCFILE does not exist.'
        RETURN
      ENDIF
!
!  ****  Nuclear charge density parameters.
!
      IF(MNUCL.EQ.1) THEN
!  ****  Point nucleus.
        WRITE(IW,1002)
 1002   FORMAT(1X,'#',/1X,'# Nuclear model: point charge')
        U1=0.0D0
        U2=0.0D0
      ELSE IF (MNUCL.EQ.2) THEN
!  ****  Uniform distribution..
        R1=1.07D0*F2BOHR*AW**0.3333333333333333D0
        R2=2.00D0*F2BOHR
        R1=R1*DSQRT((1.0D0+2.5D0*(R2/R1)**2)/(1.0D0+0.75D0*(R2/R1)**2))
        WRITE(IW,1004) R1*A0B
 1004   FORMAT(1X,'#',/1X,'# Nuclear model: uniform spherical distribution',/1X,'#',16X,'Nuclear radius =',1P,E12.5,' cm')
        U1=R1**2
        U2=0.0D0
      ELSE IF(MNUCL.EQ.4.OR.MNUCL.EQ.3) THEN
!  ****  Helm's Uu distribution.
        RNUC=1.070D0*F2BOHR*AW**0.3333333333333333D0
        R1=0.962D0*RNUC+0.435D0*F2BOHR
        R2=2.0D0*F2BOHR
        IF(R2.GT.R1) THEN
          STORE=R1
          R1=R2
          R2=STORE
        ENDIF
        WRITE(IW,1105) R1*A0B,R2*A0B
 1105   FORMAT(1X,'#',/1X,'# Nuclear model: Helm''s Uu distribution',/1X,'#',16X,'  Inner radius =',1P,E12.5, ' cm',/1X,'#',16X,'Skin thickness =',E12.5,' cm')
        U1=R1**2
        U2=R2**2
      ELSE
        WRITE(IW,1003)
 1003   FORMAT(1X,'#',/1X,'# Undefined nuclear charge density model.',/1X,'# The calculation was aborted by subroutine HEBRON.')
      ENDIF
!
      TH(1)=0.0D0
      TH(2)=1.0D-4
      I=2
   10 CONTINUE
      I=I+1
      IF(TH(I-1).LT.0.9999D-3) THEN
        TH(I)=TH(I-1)+2.5D-5
      ELSE IF(TH(I-1).LT.0.9999D-2) THEN
        TH(I)=TH(I-1)+2.5D-4
      ELSE IF(TH(I-1).LT.0.9999D-1) THEN
        TH(I)=TH(I-1)+2.5D-3
      ELSE IF(TH(I-1).LT.0.9999D+0) THEN
        TH(I)=TH(I-1)+2.5D-2
      ELSE IF(TH(I-1).LT.0.9999D+1) THEN
        TH(I)=TH(I-1)+1.0D-1
      ELSE IF(TH(I-1).LT.2.4999D+1) THEN
        TH(I)=TH(I-1)+2.5D-1
      ELSE
        TH(I)=TH(I-1)+5.0D-1
      ENDIF
      IF(TH(I).LT.180.0D0) GO TO 10
      NTAB=I
!
      NTABT=NTAB
      DO I=1,NTAB
        THR=TH(I)*PI/180.0D0
        XT(I)=(1.0D0-COS(THR))/2.0D0
!  ****  Screening correction.
        Q2=4.0D0*RK2*XT(I)
        IF(Q2.LE.Q2T(NQS)) THEN
          CALL FINDI(Q2T,Q2,NQS,J)
          F=FQ(J)+(FQ(J+1)-FQ(J))*(Q2-Q2T(J))/(Q2T(J+1)-Q2T(J))
        ELSE
          F=1.0D0
        ENDIF
!  ****  Nuclear form factor.
        QR2=Q2*U1
        QR=DSQRT(QR2)
        IF(QR2.LT.1.0D-8) THEN
          FR=1.0D0+QR2*(-0.1D0+QR2*3.5714285714285714D-3)
        ELSE
          FR=3.0D0*(DSIN(QR)-QR*DCOS(QR))/(QR*QR2)
        ENDIF
        QU2=Q2*U2
        QU=DSQRT(QU2)
        IF(QU2.LT.1.0D-8) THEN
          FU=1.0D0+QU2*(-0.1D0+QU2*3.5714285714285714D-3)
        ELSE
          FU=3.0D0*(DSIN(QU)-QU*DCOS(QU))/(QU*QU2)
        ENDIF
        FN=FR*FU
!
        RMR=DPWAC(THR)
        DCST(I)=RUTHC*(F*FN)**2*RMR/(1.0D0+Q2)**2
      ENDDO
!
!  ****  Integrated cross sections.
!
      SUM0=0.0D0
      SUM1=0.0D0
      SUM2=0.0D0
      RMUL=0.0D0
      RMUU=1.0D-16
   20 CONTINUE
      MOM=0
      CALL GABQ(DCSHB,RMUL,RMUU,SUMP0,1.0D-9,IER)
      MOM=1
      CALL GABQ(DCSHB,RMUL,RMUU,SUMP1,1.0D-9,IER)
      MOM=2
      CALL GABQ(DCSHB,RMUL,RMUU,SUMP2,1.0D-9,IER)
      SUM0=SUM0+SUMP0
      SUM1=SUM1+SUMP1
      SUM2=SUM2+SUMP2
      RMUL=RMUU
      RMUU=MIN(2.0D0*RMUL,1.0D0)
      IF(RMUL.LT.0.9999999D0) GO TO 20
      ECS0=FOURPI*SUM0
      ECS1=FOURPI*SUM1
      ECS2=FOURPI*SUM2
!
      ECS=ECS0
      TCS1=2.0D0*ECS1
      TCS2=6.0D0*(ECS1-ECS2)
!      WRITE(IW,2013) ECS,ECS/A0B2
 2013 FORMAT(1X,'#',/1X,'# Total elastic cross section =',1P,E12.5,' cm**2 =',E12.5,' a0**2')
!      WRITE(IW,2014) TCS1,TCS1/A0B2
 2014 FORMAT(1X,'# 1st transport cross section =',1P,E12.5,' cm**2 =',E12.5,' a0**2')
!      WRITE(IW,2015) TCS2,TCS2/A0B2
 2015 FORMAT(1X,'# 2nd transport cross section =',1P,E12.5,' cm**2 =',E12.5,' a0**2',/1X,'#')
!
      WRITE(IW,'(1X,''#'')')
      WRITE(IW,'(1X,''# Differential cross section'',6X,''MU=(1-COS(THETA))/2'')')
      WRITE(IW,'(1X,''#'',/1X,''#  THETA'',8X,''MU'',10X,''DCS'',10X,''DCS'',22X,''Sherman'',18X,''error''/1X,''#  (deg)'',17X,''(cm**2/sr)'',3X,''(a0**2/sr)'',9X,''S'',11X,''T'',12X,''U'',/1X,''#'',70(''-''))')
      DO I=1,NTAB
        SPOL(I)=0.0D0
        TPOL(I)=0.0D0
        UPOL(I)=0.0D0
        ERROR(I)=1.0D-5
        WRITE(IW,2018) TH(I),XT(I),DCST(I),DCST(I)/A0B2,SPOL(I),TPOL(I),UPOL(I),ERROR(I)
      ENDDO
 2018 FORMAT(1X,1P,E10.3,E13.5,5E13.5,2X,E8.1)
!
      RETURN
      END
!  *********************************************************************
      FUNCTION DCSHB(RMU)
!     Mott-Born DCS for elastic scattering of high-energy electrons and
!  positrons by neutral atoms.
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX (KIND=RegReal_K) (C),INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER (NPC=1500)
      COMMON/CRMORU/CFM(NPC),CGM(NPC),DPC(NPC),DMC(NPC),CF,CG,RUTHC,WATSC,RK2,ERRF,ERRG,NPC1
      PARAMETER (NGT=650)
      COMMON/CDCSHE/Q2T(NGT),FQ(NGT),U1,U2,NQS,MOM
!  ****  Screening correction.
      Q2=4.0D0*RK2*RMU
      IF(Q2.LE.Q2T(NQS)) THEN
        CALL FINDI(Q2T,Q2,NQS,J)
        F=FQ(J)+(FQ(J+1)-FQ(J))*(Q2-Q2T(J))/(Q2T(J+1)-Q2T(J))
      ELSE
        F=1.0D0
      ENDIF
!  ****  (nuclear form factor)**2.
      QR2=Q2*U1
      QR=DSQRT(QR2)
      IF(QR2.LT.1.0D-8) THEN
        FR=1.0D0+QR2*(-0.1D0+QR2*3.5714285714285714D-3)
      ELSE
        FR=3.0D0*(DSIN(QR)-QR*DCOS(QR))/(QR*QR2)
      ENDIF
      QU2=Q2*U2
      QU=DSQRT(QU2)
      IF(QU2.LT.1.0D-8) THEN
        FU=1.0D0+QU2*(-0.1D0+QU2*3.5714285714285714D-3)
      ELSE
        FU=3.0D0*(DSIN(QU)-QU*DCOS(QU))/(QU*QU2)
      ENDIF
      FN=FR*FU
!
      RMR=DPWAC(DACOS(1.0D0-2.0D0*RMU))
      DCSHB=(RUTHC*(F*FN)**2*RMR/(1.0D0+Q2)**2)*RMU**MOM
      RETURN
      END
!  *********************************************************************
!                        SUBROUTINE SGRID
!  *********************************************************************
      SUBROUTINE SGRID(R,DIFR,R0,RN,C,N)
!
!  This subroutine sets up a radial grid R(I) (I=1:N) such that
!     1) A*(R(I)+R0)+B*LOG(R(I)+R0)+C=I.
!     2) The grid spacing, R(I+1)-R(I), increases monotonously with I.
!     3) The parameters A and B are determined by requiring that R(1)=0
!        and R(N)=RN. They are required to be positive.
!
!     To describe the wave functions of an electron in the field of an
!  atom, the following parameter values should be adequate: R0 about
!  1.0D-5 or smaller, RN of the order of 10, N=400 or larger, C=N/2 or
!  so. Then, C is approximately equal to the number of grid points
!  between 0 and 1, R(2) is smaller than R0 in a factor of about 10, and
!  R(N)-R(N-1)=RN/(N-C), approximately.
!
!     The output arrays contain the grid points R(1:N) and the relative
!  increments DIFR(1:N), i.e. the values of the derivative of R(I) with
!  respect to I.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
!  ****  Set IWR=1 to print the grid on a file.
      PARAMETER (IWR=0)
      DIMENSION R(N),DIFR(N)
!
 1000 FORMAT(4X,'R0=',1P,E18.11,'  RN=',E18.11,/4X,' C=',E18.11,'   N=',I5,/)
!  ****  R0 must be positive (gt. 1.0d-35).
      IF(R0.LT.1.0D-35) THEN
        WRITE(6,1001)
 1001 FORMAT(1X,'** Error in SGRID: R0 is too small.')
        WRITE(6,1000) R0,RN,C,N
        STOP 'SGRID: R0 is too small.'
      ENDIF
!  ****  RN has to be larger than R0.
      IF(RN.LE.R0) THEN
        WRITE(6,1002)
 1002 FORMAT(1X,'** Error in SGRID: RN is less than or equal to R0.')
        WRITE(6,1000) R0,RN,C,N
        STOP 'SGRID: RN is less than or equal to R0.'
      ENDIF
!  ****  N should be larger than 10.
      IF(N.LT.10) THEN
        WRITE(6,1003)
 1003 FORMAT(1X,'** WARNING in SGRID: N is too small.')
        WRITE(6,1000) R0,RN,C,N
      ENDIF
!  ****  Risk of round off errors if RN/R0 is too large.
      RN0=RN/R0
      IF(RN0.GT.1.0D12) THEN
        WRITE(6,1004)
 1004 FORMAT(1X,'** WARNING in SGRID: RN/R0 is too large.')
        WRITE(6,1000) R0,RN,C,N
      ENDIF
!
      CC=C
      B=(1.0D0-CC)/DLOG(R0)
      RPR0=RN+R0
      A=(DBLE(N)-CC-B*DLOG(RPR0))/RPR0
      IF(B.LT.1.0D-15.OR.A.LT.1.0D-15) THEN
        A=0.0D0
        B=(DBLE(N)-1.0D0)/DLOG(RN0+1.0D0)
        CC=1.0D0-B*DLOG(R0)
      ENDIF
!
      R(1)=0.0D0
      RPR0=R(1)+R0
      DIFR(1)=RPR0/(A*RPR0+B)
      RR=1.0D-35
      DO I=2,N
        RL=RR
        RU=RL
    1   RU=2.0D0*RU
        RPR0=RU+R0
        FU=A*RPR0+B*DLOG(RPR0)+CC-DBLE(I)
        IF(FU.LT.0.0D0) GO TO 1
    2   RR=0.5D0*(RU+RL)
        RPR0=RR+R0
        FR=A*RPR0+B*DLOG(RPR0)+CC-DBLE(I)
        IF(FR.GT.0.0D0) THEN
          RU=RR
        ELSE
          RL=RR
        ENDIF
        IF(RU-RL.GT.1.0D-15*RR) GO TO 2
        R(I)=RR
        RPR0=RR+R0
        DIFR(I)=RPR0/(A*RPR0+B)
      ENDDO
!
!  ****  Print the grid on a file.
!
      RETURN
      END
!  *********************************************************************
!                       SUBROUTINE SLAG6
!  *********************************************************************
      SUBROUTINE SLAG6(H,Y,S,N)
!
!     Piecewise six-point Lagrange integration of a uniformly tabulated
!  function.
!
!     H ...... step length,
!     Y ...... array of function values (ordered abscissas),
!     S ...... array of integral values defined as
!              S(I)=INTEGRAL(Y) from X(1) to X(I)=X(1)+(I-1)*H,
!     N ...... number of data points.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT Double Precision (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      DIMENSION Y(N),S(N)
      IF(N.LT.6) STOP
      HR=H/1440.0D0
      Y1=0.0D0
      Y2=Y(1)
      Y3=Y(2)
      Y4=Y(3)
      S(1)=0.0D0
      S(2)=HR*(475*Y2+1427*Y3-798*Y4+482*Y(4)-173*Y(5)+27*Y(6))
      S(3)=S(2)+HR*(-27*Y2+637*Y3+1022*Y4-258*Y(4)+77*Y(5)-11*Y(6))
      DO I=4,N-2
        Y1=Y2
        Y2=Y3
        Y3=Y4
        Y4=Y(I)
        S(I)=S(I-1)+HR*(11*(Y1+Y(I+2))-93*(Y2+Y(I+1))+802*(Y3+Y4))
      ENDDO
      Y5=Y(N-1)
      Y6=Y(N)
      S(N-1)=S(N-2)+HR*(-27*Y6+637*Y5+1022*Y4-258*Y3+77*Y2-11*Y1)
      S(N)=S(N-1)+HR*(475*Y6+1427*Y5-798*Y4+482*Y3-173*Y2+27*Y1)
      RETURN
      END
!  *********************************************************************
!                       FUNCTION VCPOL
!  *********************************************************************
      FUNCTION VCPOL(IELEC,DEN)
!
!     This function gives the correlation potential of an electron
!  (IELEC=-1) or positron (IELEC=+1) in an homogeneous electron gas of
!  density DEN (electrons per unit volume).
!
!  ****  All quantities are in atomic units.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT Double Precision (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER (PI=3.1415926535897932D0,FOURPI=4.0D0*PI)
!
      IF(DEN.LT.1.0D-12) THEN
        VCPOL=0.0D0
        RETURN
      ENDIF
      RS=(3.0D0/(FOURPI*MAX(DEN,1.0D-12)))**3.333333333333D-1
      RSL=DLOG(RS)
!
      IF(IELEC.EQ.-1) THEN
!  ****  Electron exchange-correlation potential.
!        Ref:  Padial and Norcross, Phys. Rev. A 29(1984)1742.
!              Perdew and Zunger, Phys. Rev. B 23(1981)5048.
        IF(RS.LT.1.0D0) THEN
          VCPOL=0.0311D0*RSL-0.0584D0+0.00133D0*RS*RSL-0.0084D0*RS
        ELSE
          GAM=-0.1423D0
          BET1=1.0529D0
          BET2=0.3334D0
          RSS=DSQRT(RS)
          VCPOL=GAM*(1.0D0+(7.0D0/6.0D0)*BET1*RSS+(4.0D0/3.0D0)*BET2*RS)/(1.0D0+BET1*RSS+BET2*RS)**2
        ENDIF
      ELSE
!  ****  Positron correlation potential.
!        Ref:  Jain, Phys. Rev. A 41(1990)2437.
        IF(RS.LT.0.302D0) THEN
          VCPOL=(-1.82D0/SQRT(RS))+(0.051D0*RSL-0.115D0)*RSL+1.167D0
        ELSE IF(RS.LT.0.56D0) THEN
          VCPOL=-0.92305D0-0.09098D0/RS**2
        ELSE IF(RS.LT.8.0D0) THEN
          RSD=1.0D0/(RS+2.5D0)
          VCPOL=(-8.7674D0*RS*RSD**3)+(-13.151D0+0.9552D0*RS)*RSD**2+2.8655D0*RSD-0.6298D0
        ELSE
          VCPOL=-179856.2768D0*3.0D0*DEN**2+186.4207D0*2.0D0*DEN-0.524D0
        ENDIF
        VCPOL=0.5D0*VCPOL
      ENDIF
      RETURN
      END
!  *********************************************************************
!                        SUBROUTINE XSFEG
!  *********************************************************************
      SUBROUTINE XSFEG(DEN,DELTA,IELEC,EK,MORD,XSEC,IMODE)
!
!     This subroutine computes restricted (W>DELTA) total cross sections
!  for interactions of electrons (IELEC=-1) or positrons (IELEC=+1) with
!  a degenerate free electron gas (per electron in the gas). The DCS
!  is obtained from Lindhard's dielectric function (i.e. within the
!  first Born approximation), with the Ochkur exchange correction for
!  electrons.
!
!  Ref.: F. Salvat, Phys. Rev. A 68 (2003) 012708.
!
!
!  Input arguments:
!     DEN ...... density of the electron gas (electrons per unit
!                volume).
!     DELTA .... energy gap (or minimum energy loss).
!     IELEC .... kind of projectile.
!                =-1, electron; =+1, positron.
!     EK ....... kinetic energy of the projectile.
!     MORD ..... order of the calculated cross section;
!                =0, total cross section,
!                =1, stopping cross section,
!                =2, energy straggling cross section.
!     XSEC ..... total integrated cross section.
!     IMODE .... =1, the complete DCS (for electron-hole and plasmon
!                    excitations) is calculated.
!                =2, the output value XSEC corresponds to electron-hole
!                    excitations (binary collisions) only.
!
!                                  (All quantities in atomic units).
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER(F2O3=2.0D0/3.0D0,F3O16=3.0D0/16.0D0)
      PARAMETER(PI=3.1415926535897932D0, FOURPI=4.0D0*PI)
      COMMON/CXSFEG/EF,EP,CHI2,XP,ZC,XC,XE,SXE,IELPO
      PARAMETER(NPM=50)
      COMMON/CXSPL0/Z
      COMMON/CXSPL1/ZPT(NPM),XPT(NPM),FPL(NPM),ZANAL,XANAL,AX(NPM),BX(NPM),CX(NPM),DX(NPM),AF(NPM),BF(NPM),CF(NPM),DF(NPM),MOM
!  ****  Contributions from electron-hole and plasmon excitations.
      COMMON/XSFEGO/XSEH,XSPL
!
!      EXTERNAL PLSTR
!  ****  Plasmon excitation functions are printed if IWR=1.
      IWR=0
!
      IF(MORD.LT.0.OR.MORD.GT.2) THEN
          STOP 'Wrong MORD value.'
      ENDIF
!
!  ****  Constants and energy-independent parameters.
!
      DENM=MAX(DEN,1.0D-5)
      IELPO=IELEC
      EP2=FOURPI*DENM
      EP=DSQRT(EP2)
      EF=0.5D0*(3.0D0*PI**2*DENM)**F2O3
      XE=EK/EF
      IF(XE.LT.1.0D-3) THEN
        XSEC=0.0D0
        RETURN
      ENDIF
      SXE=DSQRT(XE)
      XP=EP/EF
      CHI2=F3O16*XP*XP
!  ****  Plasmon cutoff momentum.
        ZL=XP-1.0D-3
    1   CONTINUE
        FL=ZL*ZL+CHI2*F1(ZL,4.0D0*ZL*(ZL+1.0D0))
        IF(FL.GT.0.0D0) THEN
          ZL=0.5D0*ZL
          GO TO 1
        ENDIF
        ZU=XP+1.0D-2
    2   CONTINUE
        FU=ZU*ZU+CHI2*F1(ZU,4.0D0*ZU*(ZU+1.0D0))
        IF(FU.LT.0.0D0) THEN
          ZU=ZU+ZU
          GO TO 2
        ENDIF
    3   ZC=0.5D0*(ZL+ZU)
        FT=ZC*ZC+CHI2*F1(ZC,4.0D0*ZC*(ZC+1.0D0))
        IF(FT.LT.0.0D0) THEN
          ZL=ZC
        ELSE
          ZU=ZC
        ENDIF
        IF(ABS(ZL-ZU).GT.1.0D-15*ZC) GO TO 3
        XC=4.0D0*ZC*(ZC+1.0D0)
!
!  ************  Electron-hole contribution.
!
      CALL SEH0(XSEH,DELTA,MORD,IWR)
!
      IF(IMODE.EQ.2) THEN
        XSEC=XSEH
        RETURN
      ENDIF
!
!  ************  Plasmon contribution.
!
      IF(XE.LT.XP) THEN
        XSPL=0.0D0
      ELSE
!  ****  Plasmon line.
        ZPT(1)=0.0D0
        XPT(1)=XP
        FPL(1)=1.0D0
        ZANAL=0.0D0
        XANAL=0.0D0
!  ****  Varying step: 2*DZ for I<NPH.
        NPH=2*NPM/3
        DFZ=0.999999999D0*ZC/DBLE(NPM+NPH-3)
        DO I=2,NPM
          IF(I.LT.NPH) THEN
            Z=ZPT(I-1)+DFZ*2.0D0
          ELSE
            Z=ZPT(I-1)+DFZ
          ENDIF
          IF(Z.GT.0.02D0*ZC) THEN
!  The starting endpoints must be outside the Lindhard continuum.
            XL=MAX(4.0D0*Z*(Z+1.0D0)+1.0D-9,0.9D0*XP)
            XU=1.1D0*XC
    4       X=0.5D0*(XL+XU)
            FT=Z*Z+CHI2*F1(Z,X)
            IF(FT.GT.0.0D0) THEN
              XU=X
            ELSE
              XL=X
            ENDIF
!           WRITE(6,'('' X,FT ='',1P,3E18.11)') X,FT
            IF(FT.GT.1.0D-6) GO TO 4
            IF(ABS(XL-XU).GT.1.0D-13*X) GO TO 4
          ELSE
            X=DSQRT(XP**2+(48.0D0/5.0D0)*Z**2+16.0D0*Z**4)
          ENDIF
          XPT(I)=X
          ZPT(I)=Z
        ENDDO
        DO I=2,NPM-1
          Z=ZPT(I)
          XUP=4.0D0*Z*(Z+1.0D0)-1.0D-9
          CALL GABQ(PLSTR,1.0D-10,XUP,SUM,1.0D-6,IER)
          IF(IER.EQ.1) THEN
            WRITE(6,*) 'GABQ error in XSFEG.'
            STOP
          ENDIF
          FPL(I)=1.0D0-SUM*(6.0D0/(16.0D0*PI))
          XAP=XP+(24.0D0/5.0D0)*Z*Z/XP
          IF(ABS(XAP-XPT(I)).LT.1.0D-3*XPT(I).AND.FPL(I).GT.0.999D0) THEN
            ZANAL=ZPT(I)
            XANAL=XPT(I)
          ENDIF
        ENDDO
        FPL(NPM)=FPL(NPM-1)+(FPL(NPM-1)-FPL(NPM-2))*(XPT(NPM)-XPT(NPM-1))/(XPT(NPM-1)-XPT(NPM-2))
!
        CALL SPLINE(ZPT,XPT,AX,BX,CX,DX,0.0D0,0.0D0,NPM)
        CALL SPLINE(ZPT,FPL,AF,BF,CF,DF,0.0D0,0.0D0,NPM)
        CALL SPL0(XSPL,DELTA,MORD)
      ENDIF
!
      XSEC=XSEH+XSPL
      RETURN
      END
!  *********************************************************************
!                       FUNCTION PLSTR
!  *********************************************************************
      FUNCTION PLSTR(X)
!
!     Integrand of the DDCS for a point (Z,X) within the Lindhard
!  continuum.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER(PI=3.1415926535897932D0)
      COMMON/CXSFEG/EF,EP,CHI2,XP,ZC,XC,XE,SXE,IELPO
      COMMON/CXSPL0/Z
!
      PLSTR=0.0D0
      IF(Z.LT.1.0D-8) RETURN
!  ****  F2 function.
      IF(X.LT.1.0D0) THEN
        IF(X.LT.4.0D0*Z*(1.0D0-Z)) THEN
          F2=PI*X*0.125D0/Z  ! Region a.
        ELSE
          ZIN=1.0D0/Z
          ZM=Z-X*ZIN*0.25D0
          F2=PI*0.125D0*ZIN*(1.0D0-ZM*ZM)  ! Region b.
        ENDIF
      ELSE
        ZIN=1.0D0/Z
        ZM=Z-X*ZIN*0.25D0
        F2=PI*0.125D0*ZIN*(1.0D0-ZM*ZM)  ! Region b.
      ENDIF
!
      PLSTR=X*Z*Z*F2/((Z*Z+CHI2*F1(Z,X))**2+(CHI2*F2)**2)
      RETURN
      END
!  *********************************************************************
!                       SUBROUTINE SPL0
!  *********************************************************************
      SUBROUTINE SPL0(XSPL,DELTA,MORD)
!
!     Restricted total cross sections for plasmon excitation.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER(PI=3.1415926535897932D0)
      PARAMETER(NPM=50)
      COMMON/CXSFEG/EF,EP,CHI2,XP,ZC,XC,XE,SXE,IELPO
      COMMON/CXSPL1/ZPT(NPM),XPT(NPM),FPL(NPM),ZANAL,XANAL,AX(NPM),BX(NPM),CX(NPM),DX(NPM),AF(NPM),BF(NPM),CF(NPM),DF(NPM),MOM
!      EXTERNAL SPL1
!
      XSPL=0.0D0
      IF(XE.LT.XP+1.0D-8) THEN
        WRITE(6,*) 'WARNING: X is less than XP (SPL0).'
        RETURN
      ENDIF
!  ****  Minimum and maximum allowed Z-values.
      I1=0
      IN=NPM
      DO I=2,NPM
        IF(IELPO.EQ.-1) THEN
          XUP=MIN(4.0D0*ZPT(I)*(SXE-ZPT(I)),XE-1.0D0)
        ELSE
          XUP=4.0D0*ZPT(I)*(SXE-ZPT(I))
        ENDIF
        IF(XUP.GT.XPT(I)) THEN
          IF(I1.EQ.0) I1=I
        ENDIF
        IF(XUP.LT.XPT(I).AND.I1.GT.0) THEN
          IN=I
          GO TO 1
        ENDIF
      ENDDO
    1 CONTINUE
      IF(I1.EQ.0) RETURN
!
      I=I1-1
      ZL=ZPT(I)
      ZU=ZPT(I+1)
    2 Z=0.5D0*(ZL+ZU)
      X=AX(I)+Z*(BX(I)+Z*(CX(I)+Z*DX(I)))
      IF(IELPO.EQ.-1) THEN
        XMIN=MIN(4.0D0*Z*(SXE-Z),XE-1.0D0)
      ELSE
        XMIN=4.0D0*Z*(SXE-Z)
      ENDIF
      IF(XMIN.GT.X) THEN
        ZU=Z
      ELSE
        ZL=Z
      ENDIF
!       WRITE(6,'('' Z1,X-XCON ='',1P,3E18.11)') Z,X-XMIN
      IF(ABS(ZU-ZL).GT.1.0D-14*Z) GO TO 2
      ZMIN=Z
!
      IF(IN.LT.NPM) THEN
        I=IN-1
        ZL=ZPT(I)
        ZU=ZPT(I+1)
    3   Z=0.5D0*(ZL+ZU)
        X=AX(I)+Z*(BX(I)+Z*(CX(I)+Z*DX(I)))
        IF(IELPO.EQ.-1) THEN
          XMAX=MIN(4.0D0*Z*(SXE-Z),XE-1.0D0)
        ELSE
          XMAX=4.0D0*Z*(SXE-Z)
        ENDIF
        IF(XMAX.LT.X) THEN
          ZU=Z
        ELSE
          ZL=Z
        ENDIF
!         WRITE(6,'('' Z2,X-XCON ='',1P,3E18.11)') Z,X-XMAX
        IF(ABS(ZU-ZL).GT.1.0D-14*Z) GO TO 3
        ZMAX=Z
      ELSE
        XMAX=XC
        ZMAX=ZC
      ENDIF
!
      XDEL=DELTA/EF
      IF(XDEL.GE.XMAX) RETURN
      IF(XDEL.GT.XMIN) THEN
        CALL FINDI(XPT,XDEL,NPM,I)
        ZL=ZPT(I)
        ZU=ZPT(I+1)
    4   Z=0.5D0*(ZL+ZU)
        X=AX(I)+Z*(BX(I)+Z*(CX(I)+Z*DX(I)))
        IF(XDEL.LT.X) THEN
          ZU=Z
        ELSE
          ZL=Z
        ENDIF
!         WRITE(6,'('' Z1,X-XCON ='',1P,3E18.11)') Z,X-XDEL
        IF(ABS(ZU-ZL).GT.1.0D-14*Z) GO TO 4
        ZMIN=Z
        XMIN=XDEL
      ENDIF
!
      IF(XMIN.GT.XMAX) RETURN
!
!  ****  Soft plasmon excitation.
!
      FACT= 3.0D0/(16.0D0*CHI2)
      IF(XMIN.LT.XANAL.AND.XMAX.GT.XANAL) THEN
        IF(MORD.EQ.0) THEN
          X=XANAL
          S0U=X+(XP/2.0D0)*DLOG((X-XP)/(X+XP))
          X=XMIN
          S0L=X+(XP/2.0D0)*DLOG((X-XP)/(X+XP))

          SUMP=FACT*(S0U-S0L)
          ZMIN=ZANAL
        ELSE IF(MORD.EQ.1) THEN
          X=XANAL
          S1U=(X**2/2.0D0)+(XP**2/2.0D0)*DLOG(X*X-XP*XP)
          X=XMIN
          S1L=(X**2/2.0D0)+(XP**2/2.0D0)*DLOG(X*X-XP*XP)
          SUMP=FACT*(S1U-S1L)
          ZMIN=ZANAL
        ELSE IF(MORD.EQ.2) THEN
          X=XANAL
          S2U=(X**3/3.0D0)+XP**2*X+(XP**3/2.0D0)*DLOG((X-XP)/(X+XP))
          X=XMIN
          S2L=(X**3/3.0D0)+XP**2*X+(XP**3/2.0D0)*DLOG((X-XP)/(X+XP))
          SUMP=FACT*(S2U-S2L)
          ZMIN=ZANAL
        ELSE
          STOP 'Wrong MORD value.'
        ENDIF
      ELSE
        SUMP=0.0D0
      ENDIF
!
      IF(ZMIN.LT.ZMAX) THEN
        MOM=MORD
        CALL GABQ(SPL1,ZMIN,ZMAX,SUM,1.0D-6,IER)
        IF(IER.NE.0) THEN
          OPEN(9,FILE='plasma.dat')
          DO I=1,NPM
            WRITE(9,'(I4,1P,5E14.6)') I,XPT(I),ZPT(I),FPL(I)
          ENDDO
          CLOSE(9)
          WRITE(6,*) 'Accumulated numerical errors...'
          WRITE(6,*) 'GABQ error in SPL0.'
          STOP
        ENDIF
      ELSE
        SUMP=0.0D0
      ENDIF
      XSPL=(2.0D0*PI/(XE*EF*EF))*(SUM+SUMP)*EF**MORD
      RETURN
      END
!  *********************************************************************
!                       FUNCTION SPL1
!  *********************************************************************
      FUNCTION SPL1(Z)
!
!     DCS for plasmon excitations.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER(R96O5=96.0D0/5.0D0,R3O4=3.0D0/4.0D0)
      COMMON/CXSFEG/EF,EP,CHI2,XP,ZC,XC,XE,SXE,IELPO
      PARAMETER(NPM=50)
      COMMON/CXSPL1/ZPT(NPM),XPT(NPM),FPL(NPM),ZANAL,XANAL,AX(NPM),BX(NPM),CX(NPM),DX(NPM),AF(NPM),BF(NPM),CF(NPM),DF(NPM),MOM
!
      CALL FINDI(ZPT,Z,NPM,I)
      X=AX(I)+Z*(BX(I)+Z*(CX(I)+Z*DX(I)))
      FP=AF(I)+Z*(BF(I)+Z*(CF(I)+Z*DF(I)))
      SPL1=FP*X**MOM/(Z*X)
      RETURN
      END
!  *********************************************************************
!                       SUBROUTINE SEH0
!  *********************************************************************
      SUBROUTINE SEH0(XSEH,DELTA,MORD,IWR)
!
!  Restricted total cross sections for electron-hole excitations.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER(PI=3.1415926535897932D0)
      PARAMETER(NHM=150)
      DIMENSION XT(NHM),DW(NHM)
      COMMON/CXSFEG/EF,EP,CHI2,XP,ZC,XC,XE,SXE,IELPO
!
      IF(IELPO.EQ.-1) THEN
        XMAX=0.5D0*(XE-1.0D0)
!       XMAX=XE-1.0D-8  !!!! NO EXCHANGE !!!!
      ELSE
        XMAX=XE-1.0D-8
      ENDIF
      XMIN=MAX(DELTA/EF,1.0D-10)
      IF(XMIN.GE.XMAX) THEN
        XSEH=0.0D0
        RETURN
      ENDIF
!
      FACTL=6.0D0/(16.0D0*PI)
      FACTR=0.5D0
      NP=1
      XT(1)=XMIN
      DW(1)=SEH1(XT(1))
      IF(XMIN.LT.1.2D0*XC) THEN
        NS1=2*NHM/3
        DX=(MIN(1.2D0*XC,XMAX)-XMIN)/DBLE(NS1-1)
        DO I=2,NS1
          NP=I
          XT(I)=XT(I-1)+DX
          DW(I)=FACTL*SEH1(XT(I))
        ENDDO
      ENDIF
      IF(XT(NP).LT.XMAX-1.0D-10) THEN
        DFX=EXP(LOG((XMAX)/XT(NP))/DBLE(NHM-NP))
        NP1=NP+1
        ICALC=0
        DO I=NP1,NHM
          NP=I
          XT(I)=XT(I-1)*DFX
          IF(ICALC.EQ.0) THEN
            DW(I)=FACTL*SEH1(XT(I))
            DWA=FACTR/XT(I)**2
            IF(IELPO.EQ.-1) THEN  ! Exchange correction.
              FEXP=XT(I)/(XE-XT(I))
!             FEXP=1.0D0  !!!! NO EXCHANGE !!!!
              DWA=DWA*(1.0D0-FEXP*(1.0D0-FEXP))
            ENDIF
            IF(ABS(DW(I)-DWA).LT.1.0D-4*DWA) ICALC=1
          ELSE
!  ****  High-Z electron-hole excitations. Moller or Rutherford
!        differential cross section.
            DW(I)=FACTR/XT(I)**2
            IF(IELPO.EQ.-1) THEN  ! Exchange correction.
              FEXP=XT(I)/(XE-XT(I))
!             FEXP=1.0D0  !!!! NO EXCHANGE !!!!
              DW(I)=DW(I)*(1.0D0-FEXP*(1.0D0-FEXP))
            ENDIF
          ENDIF
        ENDDO
      ENDIF
      IF(NP.LT.3) THEN
        XSEH=0.0D0
        WRITE(6,*) 'WARNING: NP is too small (SEH0).'
        RETURN
      ENDIF
      DW(NP)=DEXP(DLOG(DW(NP-1))+DLOG(DW(NP-1)/DW(NP-2))*(XT(NP)-XT(NP-1))/(XT(NP-1)-XT(NP-2)))
!
      FACT=2.0D0*PI/(XE*EF*EF)
      XSEH=FACT*RMOM(XT,DW,NP,MORD)*EF**MORD
      RETURN
      END
!  *********************************************************************
!                       FUNCTION SEH1
!  *********************************************************************
      FUNCTION SEH1(X)
!
!     Integral of the DDCS over Z within the Lindhard continuum.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      COMMON/CXSFEG/EF,EP,CHI2,XP,ZC,XC,XE,SXE,IELPO
      COMMON/CXSEH1/XX
!      EXTERNAL SEH2
!
      SEH1=0.0D0
      SXP1=DSQRT(X+1.0D0)
      SXEX=DSQRT(XE-X)
      ZMIN=MAX(0.5D0*(SXE-SXEX),0.5D0*(SXP1-1.0D0))+1.0D-10
      ZMAX=MIN(0.5D0*(SXE+SXEX),0.5D0*(SXP1+1.0D0))-1.0D-10
      IF(ZMIN.GE.ZMAX) RETURN
!
      XX=X
      IF(ABS(X-XC).LT.2.0D-2*XC) THEN
        DZ=ZMAX-ZMIN
        ZMINM=ZMIN+1.0D-7*(ZMAX-ZMIN)
        CALL GABQ(SEH2,ZMINM,ZMAX,SUM,1.0D-6,IER)
        IF(IER.EQ.1) THEN
          WRITE(6,*) 'GABQ error in SEH1.'
          STOP
        ENDIF
        SEH1=SUM
      ELSE
        CALL GABQ(SEH2,ZMIN,ZMAX,SUM,1.0D-6,IER)
        IF(IER.EQ.1) THEN
          WRITE(6,*) 'GABQ error in SEH1.'
          STOP
        ENDIF
        SEH1=SUM
      ENDIF
      RETURN
      END
!  *********************************************************************
!                       FUNCTION SEH2
!  *********************************************************************
      FUNCTION SEH2(Z)
!
!     Integrand of the DDCS for a point (Z,X) within the Lindhard
!  continuum.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER(PI=3.1415926535897932D0)
      COMMON/CXSFEG/EF,EP,CHI2,XP,ZC,XC,XE,SXE,IELPO
      COMMON/CXSEH1/X
!
      SEH2=0.0D0
      IF(Z.LT.1.0D-8) RETURN
!  ****  F2 function.
      IF(X.LT.1.0D0) THEN
        IF(X.LT.4.0D0*Z*(1.0D0-Z)) THEN
          F2=PI*X*0.125D0/Z  ! Region a.
        ELSE
          ZIN=1.0D0/Z
          ZM=Z-X*ZIN*0.25D0
          F2=PI*0.125D0*ZIN*(1.0D0-ZM*ZM)  ! Region b.
        ENDIF
      ELSE
        ZIN=1.0D0/Z
        ZM=Z-X*ZIN*0.25D0
        F2=PI*0.125D0*ZIN*(1.0D0-ZM*ZM)  ! Region b.
      ENDIF
!
      SEH2=Z*F2/((Z*Z+CHI2*F1(Z,X))**2+(CHI2*F2)**2)
      IF(IELPO.EQ.-1) THEN  ! Exchange correction for electrons.
        FEXP=4.0D0*Z*Z/(XE-X)
!       FEXP=1.0D0  !!!! NO EXCHANGE !!!!
        SEH2=SEH2*(1.0D0-FEXP*(1.0D0-FEXP))
      ENDIF
      RETURN
      END
!  *********************************************************************
!                       FUNCTION F1
!  *********************************************************************
      FUNCTION F1(Z,X)
!
!     Lindhard's f_1(z,x) function.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
!
      IF(Z.LT.1.0D-5*X) THEN
        R=(Z/X)**2
        F1=-((16.0D0/3.0D0)+(256.0D0/5.0D0)*R)*R
        RETURN
      ENDIF
!
      ZIN=1.0D0/Z
      ZM=Z-X*ZIN*0.25D0
      IF(ABS(ZM).LT.1.0D-8) THEN
        AUX1=2.0D0*ZM-(4.0D0/3.0D0)*ZM**3-(4.0D0/15.0D0)*ZM**5
      ELSE
        ARGL=ABS((1.0D0+ZM)/(1.0D0-ZM))
        IF(ARGL.LT.1.0D-25.OR.ARGL.GT.1.0D25) THEN
          AUX1=0.0D0
        ELSE
          AUX1=(1.0D0-ZM**2)*LOG(ARGL)
        ENDIF
      ENDIF
!
      ZP=Z+X*ZIN*0.25D0
      IF(ABS(ZP).LT.1.0D-8) THEN
        AUX2=2.0D0*ZP-(4.0D0/3.0D0)*ZP**3-(4.0D0/15.0D0)*ZP**5
      ELSE
        ARGL=ABS((1.0D0+ZP)/(1.0D0-ZP))
        IF(ARGL.LT.1.0D-25.OR.ARGL.GT.1.0D25) THEN
          AUX2=0.0D0
        ELSE
          AUX2=(1.0D0-ZP**2)*LOG(ARGL)
        ENDIF
      ENDIF
!
      F1=0.5D0+0.125D0*(AUX1+AUX2)*ZIN
      RETURN
      END
!  *********************************************************************
!                       SUBROUTINE GABQ

!  *********************************************************************

      SUBROUTINE GABQ(FCT,XL,XU,SUM,TOL,IER)
!
!     This subroutine calculates the value SUM of the integral of the
!  (external) function FCT over the interval (XL,XU) using the 20-point
!  Gauss quadrature method with an adaptive bipartition scheme.
!
!     TOL is the tolerance, i.e. the maximum allowed relative error; it
!  should not exceed 1.0D-13. IER is an error flag; its output value is
!  0 when the required accuracy has been attained and 1 otherwise.
!
!                              Francesc Salvat. Barcelona, January 2002.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER(NP=10,NST=128,NCALLS=20000)
      DIMENSION X(NP),W(NP),S(NST),SN(NST),XR(NST),XRN(NST)
!  ****  Gauss 20-point integration formula.
!  Abscissas.
      DATA X/7.6526521133497334D-02,2.2778585114164508D-01,3.7370608871541956D-01,5.1086700195082710D-01,6.3605368072651503D-01,7.4633190646015079D-01,8.3911697182221882D-01,9.1223442825132591D-01,9.6397192727791379D-01,9.9312859918509492D-01/
!  Weights.
      DATA W/1.5275338713072585D-01,1.4917298647260375D-01,1.4209610931838205D-01,1.3168863844917663D-01,1.1819453196151842D-01,1.0193011981724044D-01,8.3276741576704749D-02,6.2672048334109064D-02,4.0601429800386941D-02,1.7614007139152118D-02/
!  ****  Error control.
      CTOL=MIN(MAX(TOL,1.0D-13),1.0D-2)
      PTOL=0.01D0*CTOL
      ERR=1.0D35
      IER=0
!  ****  Gauss integration from XL to XU.
      H=XU-XL
      SUM=0.0D0
      A=0.5D0*(XU-XL)
      B=0.5D0*(XL+XU)
      C=A*X(1)
      D=W(1)*(FCT(B+C)+FCT(B-C))
      DO I1=2,NP
        C=A*X(I1)
        D=D+W(I1)*(FCT(B+C)+FCT(B-C))
      ENDDO
      ICALL=NP+NP
      LH=1
      S(1)=D*A
      XR(1)=XL
!  ****  Adaptive bipartition scheme.
    1 CONTINUE
      HO=H
      H=0.5D0*H
      SUMR=0.0D0
      LHN=0
      DO I=1,LH
        SI=S(I)
        XA=XR(I)
        XB=XA+H
        XC=XA+HO
        A=0.5D0*(XB-XA)
        B=0.5D0*(XB+XA)
        C=A*X(1)
        D=W(1)*(FCT(B+C)+FCT(B-C))
        DO I2=2,NP
          C=A*X(I2)
          D=D+W(I2)*(FCT(B+C)+FCT(B-C))
        ENDDO
        S1=D*A
        A=0.5D0*(XC-XB)
        B=0.5D0*(XC+XB)
        C=A*X(1)
        D=W(1)*(FCT(B+C)+FCT(B-C))
        DO I3=2,NP
          C=A*X(I3)
          D=D+W(I3)*(FCT(B+C)+FCT(B-C))
        ENDDO
        S2=D*A
        ICALL=ICALL+4*NP
        S12=S1+S2
        IF(ABS(S12-SI).LE.MAX(PTOL*ABS(S12),1.0D-25)) THEN
          SUM=SUM+S12
        ELSE
          SUMR=SUMR+S12
          LHN=LHN+2
          IF(LHN.GT.NST) GO TO 2
          SN(LHN)=S2
          XRN(LHN)=XB
          SN(LHN-1)=S1
          XRN(LHN-1)=XA
        ENDIF
        IF(ICALL.GT.NCALLS) GO TO 2
      ENDDO
      ERR=ABS(SUMR)/MAX(ABS(SUMR+SUM),1.0D-25)
      IF(ERR.LT.CTOL.OR.LHN.EQ.0) RETURN
      LH=LHN
      DO I=1,LH
        S(I)=SN(I)
        XR(I)=XRN(I)
      ENDDO
      GO TO 1
!  ****  Warning (low accuracy) message.
    2 CONTINUE
      IER=1
      WRITE(6,11)
   11 FORMAT(/2X,'>>> GABQ. Gauss adaptive-bipartition quadrature.')
      WRITE(6,12) XL,XU,TOL
   12 FORMAT(2X,'XL =',1P,E19.12,',  XU =',E19.12,',  TOL =',E8.1)
      WRITE(6,13) ICALL,SUM,ERR,LHN
   13 FORMAT(2X,'NCALLS = ',I5,',  SUM =',1P,E20.13,',  ERR =',E8.1,/2X,'Number of open subintervals =',I3)
      WRITE(6,14)
   14 FORMAT(2X,'WARNING: the required accuracy has not been attained.'/)
      RETURN
      END


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!                  *****************************
!                  *  SUBROUTINE PACKAGE DPWA  *
!                  *****************************
!
!
!                                          Francesc Salvat.
!                                          Universitat de Barcelona.
!                                          July 19, 2003
!
!
!     Dirac Partial Wave Analysis for elastic scattering of electrons
!  and positrons by Coulomb fields with short-range central
!  modifications. The radial Dirac equation is solved by using
!  subroutines taken from the RADIAL FORTRAN package described in the
!  reference
!    F. Salvat, J. M. Fernandez-Varea and W. Williamson, Jr.
!       Comput. Phys. Commun. 90 (1995) 151.
!    with modifications by F. Salvat et al. (internal report,
!       University of Barcelona, 2001).
!
!     The calling sequence from the main program is:
!
!****   CALL DPWA0(EV,NDELTA,ISCH)
!
!    This subroutine determines the phase shifts. It acts as the
!  initialization routine for the evaluation of scattering amplitudes
!  and differential cross sections.
!
!****   CALL DPWA(TH,CF,CG,DCS,SPL,ERRF,ERRG)
!
!    Subroutine DPWA gives elastic scattering functions at the
!  scattering angle TH obtained from the phase shifts calculated
!  previously by subroutine DPWA0.

!
!
!            ****  All I/O energies and lengths in eV and cm, resp.
!
!  *********************************************************************
!                      SUBROUTINE DPWA0
!  *********************************************************************
      SUBROUTINE DPWA0(EV,NDELTA,ISCH)
!
!     This subroutine computes Dirac phase shifts, differential cross
!  sections and scattering amplitudes for elastic scattering of
!  electrons in central fields.
!
!  Input arguments:
!     EV ....... effective kinetic energy of the projectile (eV).
!     NDELTA ... number of required phase shifts (LT.25000).
!     ISCH ..... =1: all phase shifts are computed by solving the radial
!                    equation.
!                =2: only phase shifts of selected orders are computed
!                    from the solution of the radial equation, the
!                    others are obtained by lin-log natural cubic spline
!                    interpolation. For high energies, ISCH=2 leads to a
!                    considerable reduction of the calculation time.
!
!  Input (through the common block /FIELD/):
!     R(I) .... radial grid points (radii in increasing order). The
!               first point in the grid must be the origin, i.e. R(1)=0.
!               Repeated values are interpreted as discontinuities.
!     RV(I).... R(I) times the potential energy at R=R(I). The last
!               component, RV(NP), is assumed to be equal to the
!               asymptotic value.
!     NP ...... number of input grid points.
!
! *** NOTE: The radii and potential values, R(I) and RV(I), are in
!           atomic units.
!
!  Output (through the common block /DCSTAB/):
!     ECS ........ total cross section (cm**2)
!                    (only for finite range fields).
!     TCS1 ....... 1st transport cross section (cm**2)
!                    (only for finite range fields).
!     TCS2 ....... 2nd transport cross section (cm**2)
!                    (only for finite range fields).
!     TH(I) ...... scattering angles (in deg)
!     XT(I) ...... values of (1-COS(TH(I)))/2.0D0.
!     DCST(I) .... differential cross section per unit solid angle at
!                    TH(I) (cm**2/sr).
!     ERROR(I) ... estimated relative uncertainty of the computed DCS
!                    value.
!     NTAB ....... number of angles in the table.
!
!  NOTE: The values of ECS, TCS1 and TCS2 are computed from the DCS
!  table. This introduces a certain error (of the order of 0.01 per
!  cent) but ensures consistency of multiple scattering simulations
!  using the DCS table.
!
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX (KIND=RegReal_K) (C),INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER (A0B=5.291772083D-9)  ! Bohr radius (cm)
      PARAMETER (HREV=27.2113834D0)  ! Hartree energy (eV)
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (SL=137.03599976D0)  ! Speed of light (1/alpha)
      PARAMETER (PI=3.1415926535897932D0,FOURPI=4.0D0*PI)
!  ****  Input-output.
      PARAMETER (NDIM=1000)
      COMMON/FIELD/R(NDIM),RV(NDIM),NP
      PARAMETER (NGT=650)
      COMMON/DCSTAB/ECS,TCS1,TCS2,TH(NGT),XT(NGT),DCST(NGT),SPOL(NGT),TPOL(NGT),UPOL(NGT),ERROR(NGT),NTAB
!  ****  Link with the RADIAL package.
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RGRID/RRR(NPTG),P(NPTG),Q(NPTG),INDD(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),VD(NPPG),NVT
!  ****  Phase shifts and partial wave series coefficients.
      PARAMETER (NPC=1500,NDM=25000)
      DIMENSION DPI(NDM),DMI(NDM)
      COMMON/PHASES/DP(NDM),DM(NDM),NPH,ISUMP
      COMMON/WORK/XL(NDM),SA(NDM),SB(NDM),SC(NDM),SD(NDM),P1(NDM)
      COMMON/CSA/CFL(NDM),CGL(NDM),CFM(NDM),CGM(NDM),NPHM,IZINF
      COMMON/CRMORU/CFMC(NPC),CGMC(NPC),DPC(NPC),DMC(NPC),CFC,CGC,RUTHC,WATSC,RK2,ERRFC,ERRGC,NPC1
!
      CI=DCMPLX(0.0D0,1.0D0)
      ISUMP=0
!
      OPEN(UNIT=98,FILE='dpwa.dat')
      WRITE(98,2000)
 2000 FORMAT(//2X,'**** PARTIAL WAVE ANALYSIS (DPWA0) ',42('*')/)
!
      NDELT=NDELTA
      IF(NDELT.GT.NDM) THEN
        WRITE(98,2001)
 2001   FORMAT(/2X,'WARNING: NDELTA IS TOO LARGE')
        NDELT=NDM
      ENDIF
      IF(NDELT.LT.6) NDELT=6
      EPS=1.0D-14
      EPSCUT=1.0D-9
!
      E=EV/HREV
!
!  ****  Initialization of the RADIAL package.
!
      CALL VINT(R,RV,NP)
      ZINF=RVG(NVT)
      IF(DABS(ZINF).GT.1.0D-10) THEN
        IZINF=1
        CALL DPWAC0(ZINF,EV)
        IF(NDELT.GT.8000) NDELT=8000
      ELSE
        IZINF=0
      ENDIF
      DO I=1,NVT
        RRR(I)=RG(I)
        INDD(I)=I
      ENDDO
      NRT=NVT
!
      WRITE(98,2002) EV
 2002 FORMAT(/2X,'KINETIC ENERGY =',1P,E12.5,' eV')
      IF(E.LT.0.0D0) THEN
        WRITE(98,2003)
        WRITE(6,2003)
 2003   FORMAT(//2X,'NEGATIVE ENERGY. STOP.')
        STOP
      ENDIF
      RK=DSQRT(E*(E+2.0D0*SL*SL))/SL
!
      IF(IZINF.EQ.1) WRITE(98,2004)
 2004 FORMAT(/2X,'ONLY INNER PHASE SHIFTS ARE TABULATED')
      WRITE(98,2005)
 2005 FORMAT(/6X,'L',7X,'PHASE(SPIN UP)',5X,'PHASE(SPIN DOWN)',/2X,47('-'))
      ISCH0=ISCH
      IF(ISCH0.EQ.2.AND.EV.GT.1000.0D0) GO TO 1
!
!  ****  ISCH0=1, all phase shifts are computed by solving the radial
!        equation.
!
      L=0
      CALL DPHASE(E,EPS,PHP,-1,IER)
      IF(IER.NE.0) STOP
      PHM=0.0D0
      WRITE(98,2006) L,PHP,PHM
!      WRITE(6,2006) L,PHP,PHM
 2006 FORMAT(3X,I5,4X,1P,E16.8,4X,E16.8)
      DP(1)=PHP
      DM(1)=0.0D0
!
      IFIRST=2
   33 CONTINUE
      ISUMP=1
      TST=0.0D0
      DO I=IFIRST,NDELT
        L=I-1
        CALL DPHASE(E,EPS,PHP,-L-1,IER)
        IF(IER.NE.0) STOP
        CALL DPHASE(E,EPS,PHM,L,IER)
        IF(IER.NE.0) STOP
        DP(I)=PHP
        DM(I)=PHM
        TST=DMAX1(DABS(PHP),DABS(PHM),DABS(DP(I-1)))
        NPH=I
        WRITE(98,2006) L,PHP,PHM
!        WRITE(6,2006) L,PHP,PHM
        IF(TST.LT.EPSCUT.AND.L.GT.10) GO TO 6
!  ****  When the last phase shift (spin up) differs in more than 20 per
!  cent from the quadratic extrapolation, accumulated roundoff errors
!  may be important and the calculation of phase shifts is disc!ntinued.
        IF(I.GT.500) THEN
          DPEXT=DP(I-3)+3.0D0*(DP(I-1)-DP(I-2))
          DPMAX=MAX(ABS(DP(I-3)),ABS(DP(I-2)),ABS(DP(I-1)),ABS(DP(I)))
          IF(ABS(DP(I)-DPEXT).GT.0.20D0*DPMAX) THEN
            NPH=I-1
            WRITE(*,*) 'WARNING: Possible accumulation of round-off errors. Energy [keV]:', EV/1000.
            GO TO 6
          ENDIF
        ENDIF
      ENDDO
      WRITE(*,*) 'WARNING: TST =',TST,'. CHECK CONVERGENCE. Energy [keV]:', EV/1000.
      GO TO 6
!
!  ****  ISCH0=2, only inner phase shifts of orders L in a given grid
!        are computed from the solution of the radial equation. Phase
!        shifts of orders not included in this grid are obtained by
!        lin-log cubic spline interpolation.
!          The adopted grid is: 0(1)100(5)300(10) ...
!
!        This is a somewhat risky procedure, which is based on the
!        observed variation of the calculated phase shifts with L for
!        atomic scattering fields. When a change of sign is found, all
!        the phases are recalculated.
!
    1 L=0
      CALL DPHASE(E,EPS,PHP,-1,IER)
      IF(IER.NE.0) STOP
      PHM=0.0D0
      WRITE(98,2006) L,PHP,PHM
!      WRITE(6,2006) L,PHP,PHM
      DP(1)=PHP
      DM(1)=0.0D0
!
      LMAX=NDELT-1
      IND=0
      IADD=1
      LPP=1
    2 L=LPP
      CALL DPHASE(E,EPS,PHP,-L-1,IER)
      IF(IER.NE.0) STOP
      CALL DPHASE(E,EPS,PHM,L,IER)
      IF(IER.NE.0) STOP
!      WRITE(6,2006) L,PHP,PHM
!
      DP(L+1)=PHP
      DM(L+1)=PHM
!
      IF(L.LT.95) THEN
        WRITE(98,2006) L,PHP,PHM
      ELSE
        IF(DMAX1(DABS(PHP),DABS(PHM)).LT.EPSCUT) GO TO 3
        IND=IND+1
        XL(IND)=L
        DPI(IND)=PHP
        DMI(IND)=PHM
        IF(IND.GT.1) THEN
          S1=SIGN(PI,DPI(IND))
          S0=SIGN(PI,DPI(IND-1))
          IF(S1*S0.LT.0.0D0.AND.L.LT.500) THEN
            IF(DABS(DPI(IND-1)).LT.1.0D-6.AND.L.GT.300) THEN
              IND=IND-1
              L=XL(IND)+0.5D0
              GO TO 3
            ENDIF
            ISCH0=1
            IFIRST=MIN(L,94)
            GO TO 33
          ENDIF
          S1=SIGN(PI,DMI(IND))
          S0=SIGN(PI,DMI(IND-1))
          IF(S1*S0.LT.0.0D0.AND.L.LT.500) THEN
            IF(DABS(DMI(IND-1)).LT.1.0D-6.AND.L.GT.300) THEN
              IND=IND-1
              L=XL(IND)+0.5D0
              GO TO 3
            ENDIF
            ISCH0=1
            IFIRST=MIN(L,94)
            GO TO 33
          ENDIF
        ENDIF
        IF(IND .GT. 3) THEN
          IF(L.GT.500.AND.DABS(DPI(IND-1)).LT.1.0D-5) THEN
            I=IND
            DPEXT=DPI(I-3)+3.0D0*(DPI(I-1)-DPI(I-2))
            DPMAX=MAX(ABS(DPI(I-3)),ABS(DPI(I-2)),ABS(DPI(I-1)),ABS(DPI(I)))
            IF(ABS(DPI(I)-DPEXT).GT.0.20D0*DPMAX) THEN
              IND=I-1
              WRITE(*,*) 'WARNING: Possible accumulation of round-off errors. Energy [keV]:', EV/1000.
              GO TO 3
            ENDIF
            DMEXT=DMI(I-3)+3.0D0*(DMI(I-1)-DMI(I-2))
            DMMAX=MAX(ABS(DMI(I-3)),ABS(DMI(I-2)),ABS(DMI(I-1)),ABS(DMI(I)))
            IF(ABS(DMI(I)-DMEXT).GT.0.20D0*DMMAX) THEN
              IND=I-1
              WRITE(*,*) 'WARNING: Possible accumulation of round-off errors. Energy [keV]:', EV/1000.
              GO TO 3
            ENDIF
          ENDIF
        END IF
      ENDIF
      TST=DMAX1(DABS(PHP),DABS(PHM))
      IF(TST.LT.EPSCUT.AND.L.GT.3) GO TO 3
      IF(L.GE.LMAX) GO TO 3
!
      IF(L.GT.99) IADD=5
      IF(L.GT.299) IADD=10
      IF(L.GT.599) IADD=20
      IF(L.GT.1199) IADD=50
      IF(L.GT.2999) IADD=100
      IF(L.GT.9999) IADD=250
      LPP=L+IADD
      IF(LPP.GT.LMAX) LPP=LMAX
      GO TO 2
!
!  ****  Check consistency of sparsely tabulated phase shifts.
!        A discontinuity larger than 0.25*PI is considered as
!        a symptom of numerical inconsistencies.
!
    3 CONTINUE
      NPH=L+1
      TST=0.0D0
      DO I=1,IND
        WRITE(98,2008) INT(XL(I)+0.5D0),DPI(I),DMI(I)
 2008   FORMAT(3X,I5,4X,1P,E16.8,4X,E16.8,'  i')
        IF(I.GT.1) THEN
          TST=MAX(TST,DABS(DPI(I)-DPI(I-1)),DABS(DMI(I)-DMI(I-1)))
        ENDIF
      ENDDO
      IF(IND.LT.4) GO TO 6
      IF(TST.GT.0.25D0*PI) THEN
        WRITE(98,2009)
        WRITE(6,2009)
 2009   FORMAT(/2X,'ERROR: DIRECTLY COMPUTED PHASE SHIFTS SHOW LARGE DISCONTINUITIES.')
        STOP
      ENDIF
!
!  ****  Interpolated phase shifts (lin-log cubic spline).
!
      IF(DPI(IND).GT.0.0D0) THEN
        ITRAN=+1.0D0
      ELSE
        ITRAN=-1.0D0
      ENDIF
      DO I=1,IND
        DPI(I)=DLOG(DABS(DPI(I)))
      ENDDO
      CALL SPLINE(XL,DPI,SA,SB,SC,SD,0.0D0,0.0D0,IND)
      DO 4 I=2,NPH
        L=I-1
        IF(L.LT.95) GO TO 4
        RL=L
        CALL FINDI(XL,RL,IND,J)
        DP(I)=SA(J)+RL*(SB(J)+RL*(SC(J)+RL*SD(J)))
        DP(I)=real(ITRAN)*DEXP(DP(I))
    4 CONTINUE
!
      IF(DMI(IND).GT.0.0D0) THEN
        ITRAN=+1.0D0
      ELSE
        ITRAN=-1.0D0
      ENDIF
      DO I=1,IND
        DMI(I)=DLOG(DABS(DMI(I)))
      ENDDO
      CALL SPLINE(XL,DMI,SA,SB,SC,SD,0.0D0,0.0D0,IND)
      DO 5 I=2,NPH
        L=I-1
        IF(L.LT.95) GO TO 5
        RL=L
        CALL FINDI(XL,RL,IND,J)
        DM(I)=SA(J)+RL*(SB(J)+RL*(SC(J)+RL*SD(J)))
        IF(ITRAN.NE.0) DM(I)=real(ITRAN)*DEXP(DM(I))
    5 CONTINUE
      TST=DMAX1(DABS(DP(NPH)),DABS(DM(NPH)))
      IF(TST.GT.10.0D0*EPSCUT) THEN
        WRITE(*,*) 'WARNING: TST =',TST,'. CHECK CONVERGENCE. Energy [keV]:', EV/1000.
      ENDIF
!
!  ************  Coefficients in the partial-wave expansion.
!
    6 CONTINUE
      CFACT=1.0D0/(2.0D0*CI*RK)
      IF(IZINF.EQ.1) THEN
        CXP=CDEXP(2.0D0*CI*DP(1))
        CXPC=CDEXP(2.0D0*CI*DPC(1))
        CFL(1)=CXPC*(CXP-1)*CFACT
        CGL(1)=0.0D0
        DO I=2,NPH
          L=I-1
          CXP=CDEXP(2.0D0*CI*DP(I))
          CXM=CDEXP(2.0D0*CI*DM(I))
          CXPC=CDEXP(2.0D0*CI*DPC(I))
          CXMC=CDEXP(2.0D0*CI*DMC(I))
          CFL(I)=((L+1)*CXPC*(CXP-1)+L*CXMC*(CXM-1))*CFACT
          CGL(I)=(CXMC*(CXM-1)-CXPC*(CXP-1))*CFACT
        ENDDO
      ELSE
        CXP=CDEXP(2*CI*DP(1))
        CFL(1)=(CXP-1.0D0)*CFACT
        CGL(1)=0.0D0
        DO I=2,NPH
          L=I-1
          CXP=CDEXP(2.0D0*CI*DP(I))
          CXM=CDEXP(2.0D0*CI*DM(I))
          CFL(I)=((L+1)*(CXP-1)+L*(CXM-1))*CFACT
          CGL(I)=CXM*(1.0D0-CDEXP(2.0D0*CI*(DP(I)-DM(I))))*CFACT
        ENDDO
      ENDIF
!
!  ****  Reduced series (two iterations).
!
      IF(NPH.GE.250.AND.ISUMP.EQ.0) THEN
        DO I=1,NPH
          CFM(I)=CFL(I)
          CGM(I)=CGL(I)
        ENDDO
!
        NPHM=NPH
        DO 7 NTR=1,2
          NPHM=NPHM-1
          CFC=0.0D0
          CFP=CFM(1)
          CGC=0.0D0
          CGP=CGM(1)
          DO I=1,NPHM
            RL=I-1
            CFA=CFC
            CFC=CFP
            CFP=CFM(I+1)
            CFM(I)=CFC-CFP*(RL+1)/(RL+RL+3)-CFA*RL/(RL+RL-1)
            CGA=CGC
            CGC=CGP
            CGP=CGM(I+1)
            CGM(I)=CGC-CGP*(RL+2)/(RL+RL+3)-CGA*(RL-1)/(RL+RL-1)
          ENDDO
    7   CONTINUE
      ENDIF
!
!  ****  Scattering amplitudes and DCS.
!
      WRITE(98,2010)
 2010 FORMAT(//2X,'*** SCATTERING AMPLITUDES AND DIFFERENTIAL CROSS SECTION ***')
      WRITE(98,2011)
 2011 FORMAT(/4X,'ANGLE',6X,'DCS',7X,'ASYMMETRY',4X,'DIRECT AMPLITUDE',7X,'SPIN-FLIP AMPLITUDE',5X,'ERROR',/4X,'(deg)',3X,'(cm**2/sr)',22X,'(cm)',20X,'(cm)',/2X,91('-'))
!
!  ****  Angular grid (TH in deg).
!
      TH(1)=0.0D0
      TH(2)=1.0D-4
      I=2
   10 CONTINUE
      I=I+1
      IF(TH(I-1).LT.0.9999D-3) THEN
        TH(I)=TH(I-1)+2.5D-5
      ELSE IF(TH(I-1).LT.0.9999D-2) THEN
        TH(I)=TH(I-1)+2.5D-4
      ELSE IF(TH(I-1).LT.0.9999D-1) THEN
        TH(I)=TH(I-1)+2.5D-3
      ELSE IF(TH(I-1).LT.0.9999D+0) THEN
        TH(I)=TH(I-1)+2.5D-2
      ELSE IF(TH(I-1).LT.0.9999D+1) THEN
        TH(I)=TH(I-1)+1.0D-1
      ELSE IF(TH(I-1).LT.2.4999D+1) THEN
        TH(I)=TH(I-1)+2.5D-1
      ELSE
        TH(I)=TH(I-1)+5.0D-1
      ENDIF
      IF(I.GT.NGT) STOP 'DPWA0. The NGT parameter is too small.'
      IF(TH(I).LT.180.0D0) GO TO 10
      NTAB=I
!
      DO I=1,NTAB
        THR=TH(I)*PI/180.0D0
        XT(I)=(1.0D0-DCOS(THR))/2.0D0
        CALL DPWA(THR,CF,CG,DCS,SPL,TPL,UPL,ERRF,ERRG)
        IF(DMAX1(ERRF,ERRG).GT.0.95D0) THEN
          ERR=1.0D0
        ELSE
          ACF=CDABS(CF)**2
          ACG=CDABS(CG)**2
          ERR=2.0D0*(ACF*ERRF+ACG*ERRG)/DMAX1(DCS,1.0D-45)
        ENDIF
        DCST(I)=DCS
        ERROR(I)=DMAX1(ERR,1.0D-7)
        SPOL(I)=SPL
        TPOL(I)=TPL
        UPOL(I)=UPL
         IF (TH(I) .EQ. 0.d0 .OR. TH(I) .EQ. 180.d0) THEN
          SPOL(I)=0.d0
          TPOL(I)=1.d0
          UPOL(I)=0.d0
        ENDIF
        WRITE(98,2012) TH(I),DCST(I),SPOL(I),TPOL(I),UPOL(I),CF,CG,ERROR(I)
 2012   FORMAT(1X,1P,E10.3,E12.5,1X,3E10.3,2(1X,'(',E10.3,',',E10.3,')'),E10.2)
      ENDDO
!
!  ************  Total and momentum transfer cross sections.
!                Convergence test (only for finite range fields).
!
      IF(IZINF.EQ.0) THEN
        INC=5
        IF(ISUMP.EQ.1) INC=1
        TST1=0.0D0
        TST2=0.0D0
        ECS=4.0D0*PI*CFL(1)*DCONJG(CFL(1))
        TCS=0.0D0
        ECSO=ECS
        TCSO=TCS
        DO I=2,NPH
          L=I-1
          RL=L
          DECS=CFL(I)*DCONJG(CFL(I))+RL*(L+1)*CGL(I)*DCONJG(CGL(I))
          DECS=4.0D0*PI*DECS/(L+L+1)
          DTCS=dble(CFL(L)*DCONJG(CFL(I))+DCONJG(CFL(L))*CFL(I)+(L-1)*(RL+1)*(CGL(L)*DCONJG(CGL(I))+DCONJG(CGL(L))*CGL(I)))
          DTCS=4.0D0*PI*DTCS*L/((RL+L-1)*(L+L+1))
          ECS=ECS+DECS
          TCS=TCS+DTCS
!  ****  Convergence test.
          ITW=L-(L/INC)*INC
          IF(ITW.EQ.0) THEN
            TST1=DABS(ECS-ECSO)/(DABS(ECS)+1.0D-35)
            TST2=DABS(TCS-TCSO)/(DABS(TCS)+1.0D-35)
            ECSO=ECS
            TCSO=TCS
          ENDIF
        ENDDO
        TST=DMAX1(TST1,TST2)
        TCS=ECS-TCS
        IF(TST.GT.1.0D-5.AND.NPH.GT.40) THEN
          WRITE(*,*) 'WARNING: TST =',TST,'. CHECK CONVERGENCE. Energy [keV]:', EV/1000.
        ENDIF
        ECS=ECS*A0B2
        TCS=TCS*A0B2
!
!  ****  ECS and TCSs are evaluated from the DCS table.
!
        ECS0=FOURPI*RMOM(XT,DCST,NTAB,0)
        ECS1=FOURPI*RMOM(XT,DCST,NTAB,1)
        ECS2=FOURPI*RMOM(XT,DCST,NTAB,2)
        TST1=DABS(ECS-ECS0)/(DABS(ECS)+1.0D-35)
        WRITE(98,2013) ECS,ECS0,TST1
!        WRITE(6,2013) ECS,ECS0,TST1
 2013   FORMAT(/2X,'TOTAL ELASTIC CROSS SECTION =',1P,E13.6,' cm**2',/2X,'             FROM DCS TABLE =',E13.6,'  (REL. DIF. =',E9.2,')')
        TCS1=2.0D0*ECS1
        TCS2=6.0D0*(ECS1-ECS2)
        TST2=DABS(TCS-TCS1)/(DABS(TCS)+1.0D-35)
        WRITE(98,2014) TCS,TCS1,TST2
!        WRITE(6,2014) TCS,TCS1,TST2
 2014   FORMAT(/2X,'1ST TRANSPORT CROSS SECTION =',1P,E13.6,' cm**2',/2X,'             FROM DCS TABLE =',E13.6,'  (REL. DIF. =',E9.2,')')
        WRITE(98,2015) TCS2
!        WRITE(6,2015) TCS2
 2015   FORMAT(/2X,'2ND TRANSPORT CROSS SECTION =',1P,E13.6,' cm**2')
        TST=DMAX1(TST1,TST2)
        IF(TST.GT.2.0D-3) THEN
          WRITE(*,*) 'WARNING: RELATIVE DIFFERENCES ARE TOO LARGE THE DCS TABLE IS NOT CONSISTENT. Energy [keV]:', EV/1000.
        ENDIF
      ENDIF
!
      WRITE(98,2017)
 2017 FORMAT(/2X,'**** DPWA0 ENDED ',60('*')/)
      CLOSE(UNIT=98)
!
      RETURN
      END
!  *********************************************************************
!                       SUBROUTINE DPWA
!  *********************************************************************
      SUBROUTINE DPWA(TH,CF,CG,DCS,SPL,TPL,UPL,ERRF,ERRG)
!
!    This subroutine gives various elastic scattering functions at the
!  scattering angle TH (in radians) computed from Dirac phase shifts.
!  It should be previously initialized by calling subroutine DPWA0.
!
!  Input argument:
!     TH ....... scattering angle (in rad)
!
!  Output arguments:
!     CF ....... F scattering amplitude (cm).
!     CG ....... G scattering amplitude (cm).
!     DCS ...... differential cross section per unit solid angle for
!                unpolarized beams.
!     SPL ...... asymmetry function.
!     ERRF ..... relative uncertainty of CF.
!     ERRG ..... relative uncertainty of CG.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z),COMPLEX (KIND=RegReal_K) (C),INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER (A0B=5.291772083D-9)  ! Bohr radius (cm)
      PARAMETER (HREV=27.2113834D0)  ! Hartree energy (eV)
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (NDM=25000,NPC=1500,TOL=5.0D-8)
!  ****  Phase shifts and partial wave series coefficients.
      COMMON/PHASES/DP(NDM),DM(NDM),NPH,ISUMP
      COMMON/CSA/CFL(NDM),CGL(NDM),CFM(NDM),CGM(NDM),NPHM,IZINF
      COMMON/CRMORU/CFMC(NPC),CGMC(NPC),DPC(NPC),DMC(NPC),CFC,CGC,RUTHC,WATSC,RK2,ERRFC,ERRGC,NPC1
!
      X=DCOS(TH)
      Y=DSIN(TH)
      TST=1.0D35
!
!  ************  Reduced series method. Only when TH is greater
!                than 1.0 deg and NPH.ge.250.
!
      IF(TH.LT.1.74D-2.OR.NPH.LT.250.OR.ISUMP.EQ.1) THEN
        CFO=0.0D0
        CGO=0.0D0
        ERRFO=1.0D10
        ERRGO=1.0D10
        GO TO 10
      ENDIF
      FACT=1.0D0/(1.0D0-X)**2
!
!  ****  F scattering amplitude.
!
      P2=1.0D0
      P3=X
      CFS=CFM(1)
      CFSO=CFS
      CFS=CFS+CFM(2)*P3
      DO I=3,NPHM
        L=I-1
        P1=P2
        P2=P3
        P3=((L+L-1)*X*P2-(L-1)*P1)/L
        CTERM=CFM(I)*P3
        CFS=CFS+CTERM
!  ****  Convergence test.
        IF(L.LT.149) THEN
          INC=1
        ELSE IF(L.LT.999) THEN
          INC=5
        ELSE
          INC=25
        ENDIF
        ITW=L-(L/INC)*INC
        IF(ITW.EQ.0) THEN
          TST=CDABS(CFS-CFSO)/DMAX1(CDABS(CFS),1.0D-45)

          CFSO=CFS
        ENDIF
      ENDDO
      CF=FACT*CFS
      ERRF=TST
!
!  ****  G scattering amplitude.
!
      IF(Y.LT.1.0D-30) THEN
        CG=0.0D0
        ERRG=0.0D0
      ELSE
        P2=1.0D0
        P3=3*X
        CGS=CGM(2)
        CGSO=CGS
        CGS=CGS+CGM(3)*P3
        DO I=4,NPHM
          L=I-1
          P1=P2
          P2=P3
          P3=((L+L-1)*X*P2-L*P1)/(L-1)
          CTERM=CGM(I)*P3
          CGS=CGS+CTERM
!  ****  Convergence test.
          IF(L.LT.149) THEN
            INC=1
          ELSE IF(L.LT.999) THEN
            INC=5
          ELSE
            INC=25
          ENDIF
          ITW=L-(L/INC)*INC
          IF(ITW.EQ.0) THEN
            TST=CDABS(CGS-CGSO)/DMAX1(CDABS(CGS),1.0D-45)
            CGSO=CGS
          ENDIF
        ENDDO
        CG=FACT*Y*CGS
        ERRG=TST
      ENDIF
!
      IF(ERRF.LT.TOL.AND.ERRG.LT.TOL) GO TO 20
      CFO=CF
      ERRFO=ERRF
      CGO=CG
      ERRGO=ERRG
!
!  ************  TH smaller than 1.0 deg or NPH.LT.250 or ISUMP=1.
!
   10 CONTINUE
!  ****  If IZINF=1, scattering functions are calculated only for
!        TH larger than 0.5 deg.
      IF(IZINF.EQ.1.AND.TH.LT.0.008726D0) THEN
        CF=0.0D0
        CG=0.0D0
        ERRF=1.0D0
        ERRG=1.0D0
        DCS=1.0D-45
        SPL=0.0D0
        TPL=0.0D0
        UPL=0.0D0
        RETURN
      ENDIF
!
!  ****  F scattering amplitude.
!
      P2=1.0D0
      P3=X
      CFS=CFL(1)
      CFSO=CFS
      CFS=CFS+CFL(2)*P3
      DO I=3,NPH
        L=I-1
        P1=P2
        P2=P3
        P3=((L+L-1)*X*P2-(L-1)*P1)/L
        CTERM=CFL(I)*P3
        CFS=CFS+CTERM
!  ****  Convergence test.
        IF(L.LT.149) THEN
          INC=1
        ELSE IF(L.LT.999) THEN
          INC=5
        ELSE
          INC=25
        ENDIF
        ITW=L-(L/INC)*INC
        IF(ITW.EQ.0) THEN
          TST=CDABS(CFS-CFSO)/DMAX1(CDABS(CFS),1.0D-45)
          CFSO=CFS
        ENDIF
      ENDDO
      CF=CFS
      ERRF=TST
!
!  ****  G scattering amplitude.
!
      IF(Y.LT.1.0D-30) THEN
        CG=0.0D0
        ERRG=0.0D0
      ELSE
        P2=1.0D0
        P3=3*X
        CGS=CGL(2)
        CGSO=CGS
        CGS=CGS+CGL(3)*P3
        DO I=4,NPH
          L=I-1
          P1=P2
          P2=P3
          P3=((L+L-1)*X*P2-L*P1)/(L-1)
          CTERM=CGL(I)*P3
          CGS=CGS+CTERM
!  ****  Convergence test.
          IF(L.LT.149) THEN
            INC=1
          ELSE IF(L.LT.999) THEN
            INC=5
          ELSE
            INC=25
          ENDIF
          ITW=L-(L/INC)*INC
          IF(ITW.EQ.0) THEN
            TST=CDABS(CGS-CGSO)/DMAX1(CDABS(CGS),1.0D-45)
            CGSO=CGS
          ENDIF
        ENDDO
        CG=Y*CGS
        ERRG=TST
      ENDIF
!  ****  The following four sentences are introduced to prevent abnormal
!        termination of the calculation when the number of (inner) phase
!        shifts is small. This solves the problem found by M. Berger.
      IF(NPH.LT.20.AND.CDABS(CTERM).LT.TOL) THEN
        ERRF=0.0D0
        ERRG=0.0D0
      ENDIF
!
!  ****  Select the most accurate method.
!
      IF(ERRFO.LT.ERRF) THEN
        CF=CFO
        ERRF=ERRFO
      ENDIF
      IF(ERRGO.LT.ERRG) THEN
        CG=CGO
        ERRG=ERRGO
      ENDIF
!
!  ****  Differential cross section (unpolarized beam).
!
   20 CONTINUE
      CF=CF*A0B
      CG=CG*A0B
      IF(IZINF.EQ.1) THEN
        XAUX=DPWAC(TH)
        CFC=CFC*A0B
        CGC=CGC*A0B
        DCSM=CDABS(CFC)**2+CDABS(CGC)**2
        CF=CF+CFC
        CG=CG+CGC
        ACF=CDABS(CF)**2
        ACG=CDABS(CG)**2
        DCS=ACF+ACG
!  ****  Scattering amplitudes that are much smaller than the Coulomb
!        ones may not be correct due to rounding off.
!        (Modified Coulomb fields only).
        IF(DCS.LT.1.0D-10*DCSM.OR.ERRFC+ERRGC.GT.1.0D0.OR.TH.LT.1.74D-2) THEN
          CF=0.0D0
          CG=0.0D0
          ERRF=1.0D0
          ERRG=1.0D0
          DCS=1.0D-45
          SPL=0.0D0
          TPL=0.0d0
          UPL=0.d0
          RETURN
        ENDIF
        ERRF=ERRF+ERRFC
        ERRG=ERRG+ERRGC
      ELSE
        ACF=CDABS(CF)**2
        ACG=CDABS(CG)**2
        DCS=ACF+ACG
      ENDIF
!
      ERR=2.0D0*(ACF*ERRF+ACG*ERRG)/DMAX1(DCS,1.0D-45)
      IF(ERR.GT.0.10D0) THEN
        CF=0.0D0
        CG=0.0D0
        ERRF=1.0D0
        ERRG=1.0D0
        DCS=1.0D-45
        SPL=0.0D0
        TPL=0.d0
        UPL=0.d0
        RETURN
      ENDIF
!
!  ****  Asymmetry function.
!
      CSPL1=DCMPLX(0.0D0,1.0D0)*CF*DCONJG(CG)
      CSPL2=DCMPLX(0.0D0,1.0D0)*CG*DCONJG(CF)
      TST=CDABS(CSPL1-CSPL2)/DMAX1(CDABS(CSPL1),1.0D-45)
      IF(TST.GT.1.0D-6.AND.ERR.LT.0.01D0) THEN
        SPL=dble(CSPL1-CSPL2)/DCS
        TPL=dble(ACF-ACG)/DCS
        UPL=dble(CF*DCONJG(CG)+CF*DCONJG(CG))/DCS
        IF (SPL .NE. 0.d0 .OR. TPL .NE. 0.d0 .OR.UPL .NE. 0.d0) THEN
          FACT=1.d0/(SPL**2+TPL**2+UPL**2)
          SPL=SPL*FACT
          TPL=TPL*FACT
          UPL=UPL*FACT
        END IF
      ELSE
        SPL=0.0D0
        TPL=1.d0
        UPL=0.d0
      ENDIF
      RETURN
      END
!  *********************************************************************
!                         SUBROUTINE DPHASE
!  *********************************************************************
      SUBROUTINE DPHASE(E,EPS,PHASE,K,IER)
!
!     This subroutine computes Dirac phase shifts.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER (NDIM=1000,NPPG=NDIM+1,NPTG=NDIM+NPPG)
      PARAMETER (SL=137.03599976D0)  ! Speed of light (1/alpha)
      PARAMETER (PI=3.1415926535897932D0,PIH=0.5D0*PI)
      COMMON/RGRID/RT(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT
      COMMON/VGRID/R(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),VD(NPPG),NVT
      COMMON/STORE/PA(NPTG),QA(NPTG),PB(NPTG),QB(NPTG),D(NPTG)
!      EXTERNAL BESJN
!
      IER=0
      IF(K.EQ.0) THEN
        WRITE(6,2001)
 2001   FORMAT(1X,'*** ERROR IN DPHASE: K.EQ.0.')
        STOP
      ENDIF
!
      IF(E.LE.0.0D0) THEN
        IER=7
        WRITE(6,2002)
 2002   FORMAT(1X,'*** ERROR 7  (IN DPHASE): E.LE.0.')
        RETURN
      ENDIF
      EPSCUT=1.0D-9
!  ****  Orbital angular momentum quantum number.
      IF(K.LT.0) THEN
        L=-K-1
        KSIGN=1
      ELSE
        L=K
        KSIGN=-1
      ENDIF
      FL1=0.5D0*L*(L+1)
      RK=DSQRT(E*(E+2.0D0*SL*SL))/SL
!
!  ****  Asymptotic solution.
!
      ZINF=RV(NVT)
      IF(DABS(ZINF).LT.EPS) THEN
!  ****  Finite range fields.
        FACTOR=DSQRT(E/(E+2.0D0*SL*SL))
        ILAST=NRT+1
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          RVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
          T=EPS*RN*DABS(E*RN-FL1/RN)
          X=RK*RN
          IF(DABS(RVN).GT.T) GO TO 1
          BNL=BESJN_MOD(2,L,X)
          IF(DABS(BNL).GT.100.0D0) GO TO 1
          BNL1=BESJN_MOD(2,L+KSIGN,X)
          IF(DABS(BNL1).GT.100.0D0) GO TO 1
          BJL=BESJN_MOD(1,L,X)
          BJL1=BESJN_MOD(1,L+KSIGN,X)
          ILAST=IL
          PA(ILAST)=X*BJL
          PB(ILAST)=-X*BNL
          QA(ILAST)=-FACTOR*KSIGN*X*BJL1
          QB(ILAST)=FACTOR*KSIGN*X*BNL1
        ENDDO
    1   CONTINUE
        IF(ILAST.EQ.NRT+1) THEN
          IER=8
          WRITE(6,2003)
 2003   FORMAT(1X,'*** ERROR 8  (IN DPHASE): RAD(NGP) TOO SMALL.'/5X,'(EXTEND THE GRID TO LARGER RADII).')
          RETURN
        ENDIF
      ELSE
!  ****  Coulomb fields.
        TAS=DMAX1(1.0D-11,EPS)*DABS(ZINF)

        ILAST=NRT+1
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          RVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))
          IF(DABS(RVN-ZINF).GT.TAS) GO TO 2
          CALL DCOUL(ZINF,E,K,RN,P0,Q0,P1,Q1,ERR)
          IF(ERR.GT.EPSCUT.OR.DABS(P1).GT.100.0D0) GO TO 2
          ILAST=IL
          PA(ILAST)=P0
          PB(ILAST)=P1
          QA(ILAST)=Q0
          QB(ILAST)=Q1
        ENDDO
    2   CONTINUE
        IF(ILAST.EQ.NRT+1) THEN
          IER=8
          WRITE(6,2003)
          RETURN
        ENDIF
      ENDIF
!
!  ****  Outward solution of the radial equation.
!
      CALL DOUTW(E,EPS,K,1,NZERO,ILAST)
!
!  ****  Phase shift.
!
      RM=R(ILAST)
      IL=IND(ILAST-1)
      VF=VA(IL)/RM+VB(IL)+RM*(VC(IL)+RM*VD(IL))
      FG=(E-VF+2.0D0*SL*SL)/SL
      PO=P(ILAST)
      POP=-K*PO/RM+FG*Q(ILAST)
      IL=IND(ILAST)
      VF=VA(IL)/RM+VB(IL)+RM*(VC(IL)+RM*VD(IL))
      FG=(E-VF+2.0D0*SL*SL)/SL
      PIA=PA(ILAST)
      PIAP=-K*PIA/RM+FG*QA(ILAST)
      PIB=PB(ILAST)
      PIBP=-K*PIB/RM+FG*QB(ILAST)
!
      IF(DABS(PO).GT.EPS) THEN
        RATIO=POP/PO
        PHASE=DATAN2(RATIO*PIA-PIAP,PIBP-RATIO*PIB)
      ELSE
        PHASE=DATAN2(-PIA,PIB)
      ENDIF
      TT=DABS(PHASE)
      IF(TT.GT.PIH) PHASE=PHASE*(1.0D0-PI/TT)
!
      RETURN
      END
!  *********************************************************************
!                       SUBROUTINE DPWAC0
!  *********************************************************************
      SUBROUTINE DPWAC0(ZZP,EV)
!
!     This subroutine computes Coulomb phase shifts and initializes the
!  calculation of the Mott differential cross section for electron or
!  positron elastic scattering by a bare point nucleus.
!
!  Input:
!     ZZP....... product of nuclear and projectile charges, that is, R
!                times the interaction energy at the distance R.
!                Negative for electrons, positive for positrons.
!     EV ....... kinetic energy of the projectile (eV).
!
!  After calling DPWAC0, the function DPWAC(TH) delivers the ratio
!  (Mott DCS / Rutherford DCS) for the scattering angle TH (rad).
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z),COMPLEX (KIND=RegReal_K) (C),INTEGER (KIND=RegInt_K) (I-N)
      double precision DELTAC
!
!  ****  Set IWR=1 to print Coulomb phase shifts on a file.
      PARAMETER (IWR=0)
!
      PARAMETER (A0B=5.291772083D-9)  ! Bohr radius (cm)
      PARAMETER (HREV=27.2113834D0)  ! Hartree energy (eV)
      PARAMETER (SL=137.03599976D0)  ! Speed of light (1/alpha)
      PARAMETER (SL2=SL*SL)
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (PI=3.1415926535897932D0,PIH=0.5D0*PI)
!  ****  Phase shifts and partial-wave series coefficients.
      PARAMETER (NPC=1500)
      COMMON/CRMORU/CFM(NPC),CGM(NPC),DPC(NPC),DMC(NPC),CF,CG,RUTHC,WATSC,RK2,ERRFC,ERRGC,NPC1
!
      CI=DCMPLX(0.0D0,1.0D0)
!
!  ************  Coulomb phase shifts.
!
      E=EV/HREV
      RUTHC=(A0B*2.0D0*ZZP*(1.0D0+E/SL2))**2
      PC=DSQRT(E*(E+2.0D0*SL*SL))
      RK=PC/SL
      RK2=RK*RK
      ZETA=ZZP/SL
      W=E+SL2
      ETA=ZETA*W/PC
      RNUR=ZETA*(W+SL2)
!  ****  Negative kappa.
      DO I=1,NPC
        L=I-1
        K=-L-1
        RLAMB=DSQRT(K*K-ZETA*ZETA)
        RNUI=-1.0D0*(K+RLAMB)*PC
        RNU=DATAN2(RNUI,RNUR)
        DELTAC=-CI*CLGAM(RLAMB+CI*ETA)
        DPC(I)=RNU-(RLAMB-(L+1))*PIH+DELTAC
      ENDDO
!  ****  Positive kappa.
      DMC(1)=0.0D0
      DO I=2,NPC
        L=I-1
        K=L
        RLAMB=DSQRT(K*K-ZETA*ZETA)
        RNUI=-1.0D0*(K+RLAMB)*PC
        RNU=DATAN2(RNUI,RNUR)
        DELTAC=-CI*CLGAM(RLAMB+CI*ETA)
        DMC(I)=RNU-(RLAMB-(L+1))*PIH+DELTAC
      ENDDO
!
!  ****  Prints Coulomb phase shifts in file CPHASES.DAT
!        if the parameter IWR equals 1.
!
!  ************  Coefficients in the partial wave expansion.
!
      CXP=CDEXP(2*CI*DPC(1))
      CFACT=1.0D0/(2.0D0*CI*RK)
      CFM(1)=(CXP-1.0D0)*CFACT
      CGM(1)=0.0D0
      DO I=2,NPC
        L=I-1
        RL=L
        CXP=CDEXP(2.0D0*CI*DPC(I))
        CXM=CDEXP(2.0D0*CI*DMC(I))
        CFM(I)=((L+1)*(CXP-1)+L*(CXM-1))*CFACT
        CGM(I)=(CXM-CXP)*CFACT
      ENDDO
!
!  ****  Reduced series.
!
      NPC1=NPC
      DO NTR=1,2
        NPC1=NPC1-1
        CFC=0.0D0
        CFP=CFM(1)
        CGC=0.0D0
        CGP=CGM(1)
        DO I=1,NPC1
          RL=I-1
          CFA=CFC
          CFC=CFP
          CFP=CFM(I+1)
          CFM(I)=CFC-CFP*(RL+1)/(RL+RL+3)-CFA*RL/(RL+RL-1)
          CGA=CGC
          CGC=CGP
          CGP=CGM(I+1)
          CGM(I)=CGC-CGP*(RL+2)/(RL+RL+3)-CGA*(RL-1)/(RL+RL-1)
        ENDDO
      ENDDO
!
!  ****  Bartlett and Watson's formula for small angles.
!
      TARG=-2.0D0*CLGAM(DCMPLX(0.5D0,ETA))*CI
      C5=CDEXP(TARG*CI)
      TARG=-2.0D0*CLGAM(DCMPLX(1.0D0,ETA))*CI
      C1=CDEXP(TARG*CI)
      BETA2=E*(E+2.0D0*SL2)/(E+SL2)**2
      WATSC=-PI*BETA2*ETA*(C5/C1)
!
      RETURN
      END
!  *********************************************************************
!                       FUNCTION DPWAC
!  *********************************************************************
      FUNCTION DPWAC(TH)
!
!     Ratio (Mott DCS / Rutherford DCS) for collisions with scattering
!  angle TH (rad). Additional information is provided through the common
!  block /CRMORU/.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z),COMPLEX (KIND=RegReal_K) (C),INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER (A0B=5.291772083D-9)  ! Bohr radius (cm)
      PARAMETER (HREV=27.2113834D0)  ! Hartree energy (eV)
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (NPC=1500)
      COMMON/CRMORU/CFM(NPC),CGM(NPC),DPC(NPC),DMC(NPC),CF,CG,RUTHC,WATSC,RK2,ERRF,ERRG,NPC1
!
      X=DCOS(TH)
      Q2=2.0D0*RK2*(1.0D0-X)
      NTEST=NPC1-5
!
!  ****  TH greater than 0.5 deg.
!
      IF(TH.LT.0.008726D0) GO TO 2
      FACT=1.0D0/(1.0D0-X)**2
!  ****  Direct scattering amplitude.
      P2=1.0D0
      P3=X
      CF=CFM(1)
      CF=CF+CFM(2)*P3
      CFA=0.0D0
      DO I=3,NPC1
        L=I-1
        P1=P2
        P2=P3
        P3=((L+L-1)*X*P2-(L-1)*P1)/L
        CF=CF+CFM(I)*P3
        IF(I.EQ.NTEST) CFA=CF
      ENDDO
      ERRF=CDABS(CFA-CF)/DMAX1(CDABS(CF),1.0D-15)
      CF=FACT*CF
!  ****  Spin-flip scattering amplitude.
      Y=DSIN(TH)
      IF(Y.LT.1.0D-20) THEN
        CG=0.0D0
        ERRG=ERRF
        GO TO 1
      ENDIF
      P2=1.0D0
      P3=3*X
      CG=CGM(2)
      CG=CG+CGM(3)*P3
      CGA=0.0D0
      DO I=4,NPC1
        L=I-1
        P1=P2
        P2=P3
        P3=((L+L-1)*X*P2-L*P1)/(L-1)
        CG=CG+CGM(I)*P3
        IF(I.EQ.NTEST) CGA=CG
      ENDDO
      ERRG=CDABS(CGA-CG)/DMAX1(CDABS(CG),1.0D-15)
    1 CG=FACT*Y*CG
      PAV1=CDABS(CF)**2
      PAV2=CDABS(CG)**2
      ERR=2.0D0*(PAV1*ERRF+PAV2*ERRG)/(PAV1+PAV2)
      DCS=(PAV1+PAV2)*A0B2
      DPWAC=DCS*(Q2*Q2/RUTHC)
      IF(ERR.LT.1.0D-3.OR.TH.GT.0.08726D0) RETURN
!
!  ****  Bartlett and Watson's formula; used only for TH less than
!        5 deg, if needed. The computed DPWAC value may have slight
!        discontinuities, of the order of 0.01 per cent, between
!        0.5 and 5 deg.
!
    2 CONTINUE
      CF=0.0D0
      CG=0.0D0
      ERRF=1.0D10
      ERRG=1.0D10
      DPWAC=1.0D0+WATSC*DSIN(0.5D0*TH)
      RETURN
      END
!  *********************************************************************
!                       FUNCTION RMOM
!  *********************************************************************
      FUNCTION RMOM(X,PDF,NP,N)
!
!     Calculation of momenta of a pdf, PDF(X), obtained from linear
!  log-log interpolation on a given table. The independent variable X
!  is assumed to take only positive values.
!
!     X ..... array of variable values (in increasing order).
!     PDF ... corresponding pdf values.
!     NP .... number of points in the table.
!     N ..... moment order.
!     RMOM = INTEGRAL (X**N)*PDF(X) dX   if N.GT.-100,
!          = INTEGRAL LOG(X)*PDF(X) dX   if N.LT.-100.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER (EPS=1.0D-35)
      DIMENSION X(NP),PDF(NP)
!
      IF(NP.LT.2) STOP 'RMOM. Error code 1.'
      IF(X(1).LT.0.0D0.OR.PDF(1).LT.0.0D0) THEN
        WRITE(6,*) 'X(1),PDF(1) =',X(1),PDF(1)
        STOP 'RMOM. Error code 2'
      ENDIF
      DO I=2,NP
        IF(X(I).LT.0.0D0.OR.PDF(I).LT.0.0D0) THEN
          WRITE(6,*) 'I,X(I),PDF(I) =',I,X(I),PDF(I)
          STOP 'RMOM. Error code 3'
        ENDIF
        IF(X(I).LT.X(I-1)) STOP 'RMOM. Error code 4.'
      ENDDO
!
      IF(N.LT.-100) GO TO 1
      RMOM=0.0D0
      X2=MAX(X(1),EPS)
      Y2=PDF(1)*X(1)**N
      DO I=2,NP
        X1=X2
        Y1=Y2
        X2=X(I)
        Y2=PDF(I)*X(I)**N
        IF(Y1.GT.EPS.AND.Y2.GT.EPS) THEN
          DXL=DLOG(X2)-DLOG(X1)
          DYL=DLOG(Y2)-DLOG(Y1)
          IF(ABS(DXL).GT.1.0D-14*ABS(DYL)) THEN
            AP1=1.0D0+(DYL/DXL)
            IF(ABS(AP1).GT.1.0D-12) THEN
              DS=(Y2*X2-Y1*X1)/AP1
            ELSE
              DS=Y1*X1*DXL
            ENDIF
          ELSE
            DS=0.5D0*(Y1+Y2)*(X2-X1)
          ENDIF
          RMOM=RMOM+DS
        ENDIF
      ENDDO
      RETURN
!
    1 CONTINUE
      RMOM=0.0D0
      X2=MAX(X(1),EPS)
      Y2=PDF(1)
      DO I=2,NP
        X1=X2
        Y1=Y2
        X2=X(I)
        Y2=PDF(I)
        IF(Y1.GT.EPS.AND.Y2.GT.EPS) THEN
          DXL=DLOG(X2)-DLOG(X1)
          DYL=DLOG(Y2)-DLOG(Y1)
          IF(ABS(DXL).GT.1.0D-14*ABS(DYL)) THEN
            AP1=1.0D0+(DYL/DXL)
            IF(ABS(AP1).GT.1.0D-12) THEN
              APREC=1.0D0/AP1
              DS=(Y2*X2*(DLOG(X2)-APREC)-Y1*X1*(DLOG(X1)-APREC))*APREC
            ELSE
              DS=Y1*X1*0.5D0*(DLOG(X2)**2-DLOG(X1)**2)
            ENDIF
          ELSE
            DS=0.5D0*(Y1*DLOG(X1)+Y2*DLOG(X2))*(X2-X1)
          ENDIF
          RMOM=RMOM+DS
        ENDIF
      ENDDO
      RETURN
      END

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!  The following subroutines perform Dirac partial-wave calculations of
!  scattering of electrons and positrons in a complex central field with
!  an imaginary (absorptive) part.
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!  *********************************************************************
!                      SUBROUTINE DPWAI0
!  *********************************************************************
      SUBROUTINE DPWAI0(EV,TOTCS,ABCS,NDELTA,ISCH)
!
!     This subroutine computes Dirac phase shifts, differential cross
!  sections and scattering amplitudes for elastic scattering of
!  electrons in central fields with an imaginary (absorptive) part.
!
!  Input/output arguments:
!     EV ....... effective kinetic energy of the projectile (eV).
!     TOTCS .... total cross section (cm**2).
!     ABCS ..... absorption cross section (cm**2).
!     NDELTA ... number of required phase shifts (LT.25000).
!     ISCH ..... =1: all phase shifts are computed by solving the radial
!                    equation.
!                =2: only phase shifts of selected orders are computed
!                    from the solution of the radial equation, the
!                    others are obtained by lin-log natural cubic spline
!                    interpolation. For high energies, ISCH=2 leads to a
!                    considerable reduction of the calculation time.
!
!  Input (through the common block /FIELDI/):
!     R(I) .... radial grid points (radii in increasing order). The
!               first point in the grid must be the origin, i.e. R(1)=0.
!               Repeated values are interpreted as discontinuities.
!     RV(I).... R(I) times the potential energy at R=R(I). The last
!               component, RV(NP), is assumed to be equal to the
!               asymptotic value.
!     RW(I).... R(I) times the imaginary potential (it must be negative
!               or zero).
!     IAB ..... 0 if the potential is real, 1 if it has an imaginary
!               part.
!     NP ...... number of input grid points.
!
! *** NOTE: The radii and potential values, R(I) and RV(I), are in
!           atomic units.
!
!  Output (through the common block /DCSTAB/):
!     ECS ........ total cross section (cm**2)
!                    (only for finite range fields).
!     TCS1 ....... 1st transport cross section (cm**2)
!                    (only for finite range fields).
!     TCS2 ....... 2nd transport cross section (cm**2)
!                    (only for finite range fields).
!     TH(I) ...... scattering angles (in deg)
!     XT(I) ...... values of (1-COS(TH(I)))/2.0D0.
!     DCST(I) .... differential cross section per unit solid angle at
!                    TH(I) (cm**2/sr).
!     ERROR(I) ... estimated relative uncertainty of the computed DCS
!                    value.
!     NTAB ....... number of angles in the table.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX (KIND=RegReal_K) (C),INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER (A0B=5.291772083D-9)  ! Bohr radius (cm)
      PARAMETER (HREV=27.2113834D0)  ! Hartree energy (eV)
      PARAMETER (SL=137.03599976D0)  ! Speed of light (1/alpha)
      PARAMETER (A0B2=A0B*A0B)
      PARAMETER (PI=3.1415926535897932D0,FOURPI=4.0D0*PI)
!  ****  Input-output.
      PARAMETER (NDIM=1000)
      COMMON/FIELDI/R(NDIM),RV(NDIM),RW(NDIM),IAB,NP
      PARAMETER (NGT=650)
      COMMON/DCSTAB/ECS,TCS1,TCS2,TH(NGT),XT(NGT),DCST(NGT),SPOL(NGT),TPOL(NGT),UPOL(NGT),ERROR(NGT),NTAB
!  ****  Link with the RADIAL package.
      PARAMETER (NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RGRID/RRR(NPTG),P(NPTG),Q(NPTG),INDD(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),VD(NPPG),NVT
      COMMON/VGRIDI/RWG(NPPG),WA(NPPG),WB(NPPG),WC(NPPG),WD(NPPG)
!  ****  Phase shifts and partial wave series coefficients.
      PARAMETER (NPC=1500,NDM=25000)
      DIMENSION CXP(NDM),CXM(NDM)
      COMMON/PHASES/DP(NDM),DM(NDM),NPH,ISUMP
      COMMON/PHASEI/DPJ(NDM),DMJ(NDM)
      DIMENSION DPI(NDM),DMI(NDM),DPJI(NDM),DMJI(NDM)
      COMMON/WORK/XL(NDM),SA(NDM),SB(NDM),SC(NDM),SD(NDM),P1(NDM)
      COMMON/CSA/CFL(NDM),CGL(NDM),CFM(NDM),CGM(NDM),NPHM,IZINF
      COMMON/CRMORU/CFMC(NPC),CGMC(NPC),DPC(NPC),DMC(NPC),CFC,CGC,RUTHC,WATSC,RK2,ERRFC,ERRGC,NPC1
!
      CI=DCMPLX(0.0D0,1.0D0)
      ISUMP=0
!
      OPEN(UNIT=98,FILE='dpwai.dat')
      WRITE(98,2000)
 2000 FORMAT(//2X,'**** PARTIAL WAVE ANALYSIS (DPWAI0) ',42('*')/)
!
      NDELT=NDELTA
      IF(NDELT.GT.NDM) THEN
        WRITE(98,2001)
 2001   FORMAT(/2X,'WARNING: NDELTA IS TOO LARGE')
        NDELT=NDM
      ENDIF
      IF(NDELT.LT.6) NDELT=6
      EPS=1.0D-14
      EPSCUT=1.0D-9
!
      E=EV/HREV
!
!  ****  Initialization of the RADIAL package.
!
      IF(IAB.EQ.0) THEN
        CALL VINT(R,RV,NP)
      ELSE
        CALL VINTI(R,RV,RW,NP)
      ENDIF
!
      ZINF=RVG(NVT)
      IF(DABS(ZINF).GT.1.0D-10) THEN
        IZINF=1
        CALL DPWAC0(ZINF,EV)
        IF(NDELT.GT.8000) NDELT=8000
      ELSE
        IZINF=0
      ENDIF
      DO I=1,NVT
        RRR(I)=RG(I)
        INDD(I)=I
      ENDDO
      NRT=NVT
!
      WRITE(98,2002) EV
 2002 FORMAT(/2X,'KINETIC ENERGY =',1P,E12.5,' eV')
      IF(E.LT.0.0D0) THEN
        WRITE(98,2003)
        WRITE(6,2003)
 2003   FORMAT(//2X,'NEGATIVE ENERGY. STOP.')
        STOP
      ENDIF
      RK=DSQRT(E*(E+2.0D0*SL*SL))/SL
!
      IF(IZINF.EQ.1) WRITE(98,2004)
 2004 FORMAT(/2X,'ONLY INNER PHASE SHIFTS ARE TABULATED')
      WRITE(98,2005)
 2005 FORMAT(/14X,'--------- SPIN UP ---------',6X,'-------- SPIN DOWN --------',/6X,'L',9X,'Re(phase)      Im(phase)',9X,'Re(phase)      Im(phase)',/2X,74('-'))
      ISCH0=ISCH
      IF(ISCH0.EQ.2.AND.EV.GT.1000.0D0) GO TO 1
!
!  ****  ISCH0=1, all phase shifts are computed by solving the radial
!        equation.
!
      L=0
      IF(IAB.EQ.0) THEN
        CALL DPHASE(E,EPS,PHP,-1,IER)
        CXP(1)=CDEXP(2.0D0*CI*PHP)
      ELSE
        CALL ZDPHAS(E,EPS,CXP(1),-1,IER)
      ENDIF
      IF(IER.NE.0) STOP
      CXM(1)=1.0D0
      CDP=CDLOG(CXP(1))/(2.0D0*CI)
      CDM=CDLOG(CXM(1))/(2.0D0*CI)
      WRITE(98,2006) L,CDP,CDM
!      WRITE(6,2006) L,CDP,CDM
 2006 FORMAT(3X,I5,4X,1P,E16.8,1X,E12.5,4X,E16.8,1X,E12.5)
      DP(1)=CDP
      DM(1)=0.0D0
!
      IFIRST=2
   33 CONTINUE
      ISUMP=1
      TST=0.0D0
      DO I=IFIRST,NDELT
        L=I-1
        IF(IAB.EQ.0) THEN
          CALL DPHASE(E,EPS,PHP,-L-1,IER)
          CXP(I)=CDEXP(2.0D0*CI*PHP)
        ELSE
          CALL ZDPHAS(E,EPS,CXP(I),-L-1,IER)
        ENDIF
        IF(IER.NE.0) STOP
        IF(IAB.EQ.0) THEN
          CALL DPHASE(E,EPS,PHM,L,IER)
          CXM(I)=CDEXP(2.0D0*CI*PHM)
        ELSE
          CALL ZDPHAS(E,EPS,CXM(I),L,IER)
        ENDIF
        IF(IER.NE.0) STOP
        CDP=CDLOG(CXP(I))/(2.0D0*CI)
        CDM=CDLOG(CXM(I))/(2.0D0*CI)
        DP(I)=CDP
        DM(I)=CDM
        TST=MAX(CDABS(CDP),CDABS(CDM),DSQRT(DP(I-1)**2+DM(I-1)**2))
        NPH=I
        WRITE(98,2006) L,CDP,CDM
!        WRITE(6,2006) L,CDP,CDM
        IF(TST.LT.EPSCUT.AND.L.GT.30) GO TO 6
!  ****  When the last phase shift (spin up) differs in more than 20 per
!  cent from the quadratic extrapolation, accumulated roundoff errors
!  may be important and the calculation of phase shifts is discontinued.
        IF(I.GT.500) THEN
          DPEXT=DP(I-3)+3.0D0*(DP(I-1)-DP(I-2))
          DPMAX=MAX(ABS(DP(I-3)),ABS(DP(I-2)),ABS(DP(I-1)),ABS(DP(I)))
          IF(ABS(DP(I)-DPEXT).GT.0.20D0*DPMAX) THEN
            NPH=I-1
            WRITE(*,*) 'WARNING: Possible accumulation of round-off errors. Energy [keV]:', EV/1000.
            GO TO 6
          ENDIF
        ENDIF
      ENDDO
      WRITE(*,*) 'WARNING: TST =',TST,'. CHECK CONVERGENCE. Energy [keV]:', EV/1000.
      GO TO 6
!
!  ****  ISCH0=2, only inner phase shifts of orders L in a given grid
!        are computed from the solution of the radial equation. Phase
!        shifts of orders not included in this grid are obtained by
!        lin-log cubic spline interpolation.
!          The adopted grid is: 0(1)100(5)300(10) ...
!
!        This is a somewhat risky procedure, which is based on the
!        observed variation of the calculated phase shifts with L for
!        atomic scattering fields. When a change of sign is found, all
!        the phases are recalculated.
!
    1 L=0
      IF(IAB.EQ.0) THEN
        CALL DPHASE(E,EPS,PHP,-1,IER)
        CXP(1)=CDEXP(2.0D0*CI*PHP)
      ELSE
        CALL ZDPHAS(E,EPS,CXP(1),-1,IER)
      ENDIF
      IF(IER.NE.0) STOP
      CXM(1)=1.0D0
      CDP=CDLOG(CXP(1))/(2.0D0*CI)
      CDM=CDLOG(CXM(1))/(2.0D0*CI)
      WRITE(98,2006) L,CDP,CDM
!      WRITE(6,2006) L,CDP,CDM
      DP(1)=CDP
      DPJ(1)=-CI*CDP
      DM(1)=0.0D0
      DMJ(1)=0.0D0
!
      LMAX=NDELT-1
      IND=0
      IADD=1
      LPP=1
    2 L=LPP
      IF(IAB.EQ.0) THEN
        CALL DPHASE(E,EPS,PHP,-L-1,IER)
        CXP(L+1)=CDEXP(2.0D0*CI*PHP)
      ELSE
        CALL ZDPHAS(E,EPS,CXP(L+1),-L-1,IER)
      ENDIF
      IF(IER.NE.0) STOP
      IF(IAB.EQ.0) THEN
        CALL DPHASE(E,EPS,PHM,L,IER)
        CXM(L+1)=CDEXP(2.0D0*CI*PHM)
      ELSE
        CALL ZDPHAS(E,EPS,CXM(L+1),L,IER)
      ENDIF
      IF(IER.NE.0) STOP
      CDP=CDLOG(CXP(L+1))/(2.0D0*CI)
      CDM=CDLOG(CXM(L+1))/(2.0D0*CI)
!      WRITE(6,2006) L,CDP,CDM
!
      DP(L+1)=CDP
      DPJ(L+1)=-CI*CDP
      DM(L+1)=CDM
      DMJ(L+1)=-CI*CDM
!
      IF(L.LT.95) THEN
        WRITE(98,2006) L,CDP,CDM
      ELSE
        IF(DMAX1(CDABS(CDP),CDABS(CDM)).LT.EPSCUT) GO TO 3
        IND=IND+1
        XL(IND)=L
        DPI(IND)=DP(L+1)
        DPJI(IND)=DPJ(L+1)
        DMI(IND)=DM(L+1)
        DMJI(IND)=DMJ(L+1)
        IF(IND.GT.1) THEN
          S1=SIGN(PI,DPI(IND))
          S0=SIGN(PI,DPI(IND-1))
          IF(S1*S0.LT.0.0D0) THEN
            IF(DABS(DPI(IND-1)).LT.1.0D-6.AND.L.GT.300) THEN
              IND=IND-1
              L=XL(IND)+0.5D0
              write(*,*) L, XL(IND)
              GO TO 3
            ENDIF
            ISCH0=1
            IFIRST=MIN(L,94)
            GO TO 33
          ENDIF
          S1=SIGN(PI,DMI(IND))
          S0=SIGN(PI,DMI(IND-1))
          IF(S1*S0.LT.0.0D0) THEN
            IF(DABS(DMI(IND-1)).LT.1.0D-6.AND.L.GT.300) THEN
              IND=IND-1
              L=XL(IND)+0.5D0
              write(*,*) L, XL(IND)
              GO TO 3
            ENDIF
            ISCH0=1
            IFIRST=MIN(L,94)
            GO TO 33
          ENDIF
        ENDIF
        IF(L.GT.500.AND.DABS(DPI(IND-1)).LT.1.0D-5) THEN
          I=IND
          DPEXT=DPI(I-3)+3.0D0*(DPI(I-1)-DPI(I-2))
          DPMAX=MAX(ABS(DPI(I-3)),ABS(DPI(I-2)),ABS(DPI(I-1)),ABS(DPI(I)))
          IF(ABS(DPI(I)-DPEXT).GT.0.20D0*DPMAX) THEN
            IND=I-1
            WRITE(*,*) 'WARNING: Possible accumulation of round-off errors. Energy [keV]:', EV/1000.
            GO TO 3
          ENDIF
          DMEXT=DMI(I-3)+3.0D0*(DMI(I-1)-DMI(I-2))
          DMMAX=MAX(ABS(DMI(I-3)),ABS(DMI(I-2)),ABS(DMI(I-1)),ABS(DMI(I)))
          IF(ABS(DMI(I)-DMEXT).GT.0.20D0*DMMAX) THEN
            IND=I-1
            WRITE(*,*) 'WARNING: Possible accumulation of round-off errors. Energy [keV]:', EV/1000.
            GO TO 3
          ENDIF
        ENDIF
      ENDIF
      TST=DMAX1(CDABS(CDP),CDABS(CDM))
      IF(TST.LT.EPSCUT.AND.L.GT.3) GO TO 3
      IF(L.GE.LMAX) GO TO 3
!
      IF(L.GT.99) IADD=5
      IF(L.GT.299) IADD=10
      IF(L.GT.599) IADD=20
      IF(L.GT.1199) IADD=50
      IF(L.GT.2999) IADD=100
      IF(L.GT.9999) IADD=250
      LPP=L+IADD
      IF(LPP.GT.LMAX) LPP=LMAX
      GO TO 2
!
!  ****  Check consistency of sparsely tabulated phase shifts.
!        A discontinuity larger than 0.25*PI is considered as
!        a symptom of numerical inconsistencies.
!
    3 CONTINUE
      NPH=L+1
      TST=0.0D0
      DO I=1,IND
        WRITE(98,2008) INT(XL(I)+0.5D0),DPI(I),DPJI(I),DMI(I),DMJI(I)
 2008   FORMAT(3X,I5,4X,1P,E16.8,1X,E12.5,4X,E16.8,1X,E12.5,'  i')
        IF(I.GT.1) THEN
          TST=MAX(TST,DABS(DPI(I)-DPI(I-1)),DABS(DMI(I)-DMI(I-1)))
        ENDIF
      ENDDO
      IF(IND.LT.4) GO TO 6
      IF(TST.GT.0.25D0*PI) THEN
        WRITE(98,2009)
        WRITE(6,2009)
 2009   FORMAT(/2X,'ERROR: DIRECTLY COMPUTED PHASE SHIFTS SHOW LARGE DISCONTINUITIES.')
        STOP
      ENDIF
!
!  ****  Interpolated phase shifts (lin-log cubic spline).
!
      IF(DPI(IND).GT.0.0D0) THEN
        ITRAN=+1.0D0
      ELSE
        ITRAN=-1.0D0
      ENDIF
      DO I=1,IND
        DPI(I)=DLOG(DABS(DPI(I)))
      ENDDO
      CALL SPLINE(XL,DPI,SA,SB,SC,SD,0.0D0,0.0D0,IND)
      DO 4 I=2,NPH
        L=I-1
        IF(L.LT.95) GO TO 4
        RL=L
        CALL FINDI(XL,RL,IND,J)
        DP(I)=SA(J)+RL*(SB(J)+RL*(SC(J)+RL*SD(J)))
        DP(I)=real(ITRAN)*DEXP(DP(I))
    4 CONTINUE
!
      IF(DPJI(IND).GT.0.0D0) THEN
        ITRAN=+1.0D0
      ELSE
        ITRAN=-1.0D0
      ENDIF
      NHIM=0
      DO I=1,IND
        IF(DABS(DPJI(I)).LT.1.0D-14.AND.NHIM.GT.4) GO TO 44
        DPJI(I)=DLOG(DABS(DPJI(I)))
        NHIM=I
      ENDDO
   44 CONTINUE
      CALL SPLINE(XL,DPJI,SA,SB,SC,SD,0.0D0,0.0D0,NHIM)
      DO 444 I=2,NPH
        L=I-1
        IF(L.LT.95) GO TO 444
        RL=L
        IF(RL.LE.XL(NHIM)) THEN

          CALL FINDI(XL,RL,IND,J)
          DPJ(I)=SA(J)+RL*(SB(J)+RL*(SC(J)+RL*SD(J)))
          DPJ(I)=real(ITRAN)*DEXP(DPJ(I))
        ELSE
          DPJ(I)=0.0D0
        ENDIF
  444 CONTINUE
!

      IF(DMI(IND).GT.0.0D0) THEN
        ITRAN=+1.0D0
      ELSE
        ITRAN=-1.0D0
      ENDIF
      DO I=1,IND
        DMI(I)=DLOG(DABS(DMI(I)))
      ENDDO
      CALL SPLINE(XL,DMI,SA,SB,SC,SD,0.0D0,0.0D0,IND)
      DO 5 I=2,NPH
        L=I-1
        IF(L.LT.95) GO TO 5
        RL=L
        CALL FINDI(XL,RL,IND,J)
        DM(I)=SA(J)+RL*(SB(J)+RL*(SC(J)+RL*SD(J)))
        IF(ITRAN.NE.0) DM(I)=real(ITRAN)*DEXP(DM(I))
    5 CONTINUE
!
      IF(DMJI(IND).GT.0.0D0) THEN
        ITRAN=+1.0D0
      ELSE
        ITRAN=-1.0D0
      ENDIF
      NHIM=0
      DO I=1,IND
        IF(DABS(DPJI(I)).LT.1.0D-14.AND.NHIM.GT.4) GO TO 55
        DMJI(I)=DLOG(DABS(DMJI(I)))

        NHIM=I
      ENDDO
   55 CONTINUE
      CALL SPLINE(XL,DMJI,SA,SB,SC,SD,0.0D0,0.0D0,NHIM)
      DO 555 I=2,NPH
        L=I-1
        IF(L.LT.95) GO TO 555
        RL=L
        IF(RL.LE.XL(NHIM)) THEN
          CALL FINDI(XL,RL,IND,J)
          DMJ(I)=SA(J)+RL*(SB(J)+RL*(SC(J)+RL*SD(J)))
          DMJ(I)=real(ITRAN)*DEXP(DMJ(I))
        ELSE
          DMJ(I)=0.0D0
        ENDIF
  555 CONTINUE
!
      TST=DMAX1(DABS(DP(NPH)),DABS(DM(NPH)))
      IF(TST.GT.EPSCUT) THEN
        WRITE(*,*) 'WARNING: TST =',TST,'. CHECK CONVERGENCE. Energy [keV]:', EV/1000.
      ENDIF
      DO I=1,NPH
        CXP(I)=CDEXP(2.0D0*CI*DCMPLX(DP(I),DPJ(I)))
        CXM(I)=CDEXP(2.0D0*CI*DCMPLX(DM(I),DMJ(I)))
      ENDDO
!
!  ************  Coefficients in the partial-wave expansion.
!
    6 CONTINUE
      CFACT=1.0D0/(2.0D0*CI*RK)
      IF(IZINF.EQ.1) THEN
        CXPC=CDEXP(2*CI*DPC(1))
        CFL(1)=CXPC*(CXP(1)-1)*CFACT
        CGL(1)=0.0D0
        DO I=2,NPH
          L=I-1
          CXPC=CDEXP(2.0D0*CI*DPC(I))
          CXMC=CDEXP(2.0D0*CI*DMC(I))
          CFL(I)=((L+1)*CXPC*(CXP(I)-1)+L*CXMC*(CXM(I)-1))*CFACT
          CGL(I)=(CXMC*(CXM(I)-1)-CXPC*(CXP(I)-1))*CFACT
        ENDDO
      ELSE
        CFL(1)=(CXP(1)-1.0D0)*CFACT
        CGL(1)=0.0D0
        DO I=2,NPH
          L=I-1
          CFL(I)=((L+1)*(CXP(I)-1)+L*(CXM(I)-1))*CFACT
          CGL(I)=(CXM(I)-CXP(I))*CFACT
        ENDDO
      ENDIF
!
!  ****  Reduced series (two iterations).
!
      IF(NPH.GE.250.AND.ISUMP.EQ.0) THEN
        DO I=1,NPH
          CFM(I)=CFL(I)
          CGM(I)=CGL(I)
        ENDDO
!
        NPHM=NPH
        DO 7 NTR=1,2
          NPHM=NPHM-1
          CFC=0.0D0
          CFP=CFM(1)
          CGC=0.0D0
          CGP=CGM(1)
          DO I=1,NPHM
            RL=I-1
            CFA=CFC
            CFC=CFP
            CFP=CFM(I+1)
            CFM(I)=CFC-CFP*(RL+1)/(RL+RL+3)-CFA*RL/(RL+RL-1)
            CGA=CGC
            CGC=CGP
            CGP=CGM(I+1)
            CGM(I)=CGC-CGP*(RL+2)/(RL+RL+3)-CGA*(RL-1)/(RL+RL-1)
          ENDDO
    7   CONTINUE
      ENDIF
!
!  ****  Scattering amplitudes and DCS.
!
      WRITE(98,2010)
 2010 FORMAT(//2X,'*** SCATTERING AMPLITUDES AND DIFFERENTIAL CROSS SECTION ***')
      WRITE(98,2011)
 2011 FORMAT(/4X,'ANGLE',6X,'DCS',7X,'ASYMMETRY',4X,'DIRECT AMPLITUDE',7X,'SPIN-FLIP AMPLITUDE',5X,'ERROR',/4X,'(deg)',3X,'(cm**2/sr)',22X,'(cm)',20X,'(cm)',/2X,91('-'))
!
!  ****  Angular grid (TH in deg).
!
      TH(1)=0.0D0
      TH(2)=1.0D-4
      I=2
   10 CONTINUE
      I=I+1
      IF(TH(I-1).LT.0.9999D-3) THEN
        TH(I)=TH(I-1)+2.5D-5
      ELSE IF(TH(I-1).LT.0.9999D-2) THEN
        TH(I)=TH(I-1)+2.5D-4
      ELSE IF(TH(I-1).LT.0.9999D-1) THEN
        TH(I)=TH(I-1)+2.5D-3
      ELSE IF(TH(I-1).LT.0.9999D+0) THEN
        TH(I)=TH(I-1)+2.5D-2
      ELSE IF(TH(I-1).LT.0.9999D+1) THEN
        TH(I)=TH(I-1)+1.0D-1
      ELSE IF(TH(I-1).LT.2.4999D+1) THEN
        TH(I)=TH(I-1)+2.5D-1
      ELSE
        TH(I)=TH(I-1)+5.0D-1
      ENDIF
      IF(I.GT.NGT) STOP 'DPWA0. The NGT parameter is too small.'
      IF(TH(I).LT.180.0D0) GO TO 10
      NTAB=I
!
      DO I=1,NTAB
        THR=TH(I)*PI/180.0D0
        XT(I)=(1.0D0-DCOS(THR))/2.0D0
        CALL DPWA(THR,CF,CG,DCS,SPL,TPL,UPL,ERRF,ERRG)
        IF(DMAX1(ERRF,ERRG).GT.0.95D0) THEN
          ERR=1.0D0
        ELSE
          ACF=CDABS(CF)**2
          ACG=CDABS(CG)**2
          ERR=2.0D0*(ACF*ERRF+ACG*ERRG)/DMAX1(DCS,1.0D-45)
        ENDIF
        DCST(I)=DCS
        ERROR(I)=DMAX1(ERR,1.0D-7)
        SPOL(I)=SPL
        TPOL(I)=TPL
        UPOL(I)=UPL
        WRITE(98,2012) TH(I),DCST(I),SPOL(I),TPOL(I),UPOL(I),CF,CG,ERROR(I)
 2012   FORMAT(1X,1P,E10.3,E12.5,1X,3E10.3,2(1X,'(',E10.3,',',E10.3,')'),E10.2)
      ENDDO
!
!  ************  Total and momentum transfer cross sections.
!                Convergence test (only for finite range fields).
!
      IF(IZINF.EQ.0) THEN
        INC=5
        IF(ISUMP.EQ.1) INC=1
        TST1=0.0D0
        TST2=0.0D0
        ECS=4.0D0*PI*CFL(1)*DCONJG(CFL(1))
        TCS=0.0D0
        ECSO=ECS
        TCSO=TCS
        DO I=2,NPH
          L=I-1
          RL=L
          DECS=CFL(I)*DCONJG(CFL(I))+RL*(L+1)*CGL(I)*DCONJG(CGL(I))
          DECS=4.0D0*PI*DECS/(L+L+1)
          DTCS=CFL(L)*DCONJG(CFL(I))+DCONJG(CFL(L))*CFL(I)+(L-1)*(RL+1)*(CGL(L)*DCONJG(CGL(I))+DCONJG(CGL(L))*CGL(I))
          DTCS=4.0D0*PI*DTCS*L/((RL+L-1)*(L+L+1))
          ECS=ECS+DECS
          TCS=TCS+DTCS
!  ****  Convergence test.
          ITW=L-(L/INC)*INC
          IF(ITW.EQ.0) THEN
            TST1=DABS(ECS-ECSO)/(DABS(ECS)+1.0D-35)
            TST2=DABS(TCS-TCSO)/(DABS(TCS)+1.0D-35)
            ECSO=ECS
            TCSO=TCS
          ENDIF
        ENDDO
        TST=DMAX1(TST1,TST2)
        TCS=ECS-TCS
        IF(TST.GT.1.0D-5.AND.NPH.GT.40) THEN
          WRITE(*,*) 'WARNING: TST =',TST,'. CHECK CONVERGENCE. Energy [keV]:', EV/1000.
        ENDIF
        ECS=ECS*A0B2
        TCS=TCS*A0B2
!
!  ****  ECS and TCSs are evaluated from the DCS table.
!
        ECS0=FOURPI*RMOM(XT,DCST,NTAB,0)
        ECS1=FOURPI*RMOM(XT,DCST,NTAB,1)
        ECS2=FOURPI*RMOM(XT,DCST,NTAB,2)
        TST1=DABS(ECS-ECS0)/(DABS(ECS)+1.0D-35)
        WRITE(98,2013) ECS,ECS0,TST1
        WRITE(6,2013) ECS,ECS0,TST1
 2013   FORMAT(/2X,'TOTAL ELASTIC CROSS SECTION =',1P,E13.6,' cm**2',/2X,'             FROM DCS TABLE =',E13.6,'  (REL. DIF. =',E9.2,')')
        TCS1=2.0D0*ECS1
        TCS2=6.0D0*(ECS1-ECS2)
        TST2=DABS(TCS-TCS1)/(DABS(TCS)+1.0D-35)
        WRITE(98,2014) TCS,TCS1,TST2
        WRITE(6,2014) TCS,TCS1,TST2
 2014   FORMAT(/2X,'1ST TRANSPORT CROSS SECTION =',1P,E13.6,' cm**2',/2X,'             FROM DCS TABLE =',E13.6,'  (REL. DIF. =',E9.2,')')
        WRITE(98,2015) TCS2
        WRITE(6,2015) TCS2
 2015   FORMAT(/2X,'2ND TRANSPORT CROSS SECTION =',1P,E13.6,' cm**2')
        TST=DMAX1(TST1,TST2)
        IF(TST.GT.2.0D-3) THEN
          WRITE(*,*) 'WARNING: RELATIVE DIFFERENCES ARE TOO LARGE THE DCS TABLE IS NOT CONSISTENT. Energy [keV]:', EV/1000.
        ENDIF
!
!  ****  Absorption cross section.
!
        CALL DPWA(0.0D0,CF,CG,DCS,SPL,TPL,UPL,ERRF,ERRG)
        TOTCS=(FOURPI*(A0B/RK))*(-CI*CF)
        ABCS=TOTCS-ECS
        WRITE(98,2018) TOTCS
        WRITE(6,2018) TOTCS
 2018   FORMAT(/2X,'  GRAND TOTAL CROSS SECTION =',1P,E13.6,' cm**2')
        WRITE(98,2019) ABCS
        WRITE(6,2019) ABCS
 2019   FORMAT(2X,'   ABSORPTION CROSS SECTION =',1P,E13.6,' cm**2')
      ENDIF
!
      WRITE(98,2017)
 2017 FORMAT(/2X,'**** DPWAI0 ENDED ',60('*')/)
      CLOSE(UNIT=98)
!
      RETURN
      END
!  **************************************************************
!                        SUBROUTINE VINTI
!  **************************************************************
      SUBROUTINE VINTI(R,RV,RW,NV)
!
!     NATURAL CUBIC SPLINE INTERPOLATION FOR R*V(R) FROM THE
!  INPUT RADII AND POTENTIAL VALUES.
!
!  ****  Complex potential.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER (NDIM=1000,NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),VD(NPPG),NVT
      COMMON/VGRIDI/RWG(NPPG),WA(NPPG),WB(NPPG),WC(NPPG),WD(NPPG)
      COMMON/RGRID/X(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT
      COMMON/STORE/Y(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)
      DIMENSION Z(NPTG),CA(NPTG),CB(NPTG),CC(NPTG),CD(NPTG)
      DIMENSION R(NDIM),RV(NDIM),RW(NDIM)
!
      IF(R(1).LT.0.0D0) THEN
      WRITE(6,2101)
 2101 FORMAT(1X,'*** ERROR IN VINTI: R(1).LT.0.')
      STOP
      ENDIF
      IF(NV.GT.NDIM) THEN
      WRITE(6,2102) NDIM
 2102 FORMAT(1X,'*** ERROR IN VINTI: INPUT POTENTIAL GRID WITH MORE THAN ',I5,' DATA POINTS.')
      STOP
      ENDIF
      R(1)=0.0D0
!
      IO=0
      I=0
      K=0
    1 I=I+1
      K=K+1
      X(K)=R(I)
      Y(K)=RV(I)
      Z(K)=RW(I)
      IF(I.EQ.NV) GO TO 2
      IF(R(I).LT.R(I+1)-1.0D-12) GO TO 1
    2 CONTINUE
!
      CALL SPLINE(X,Y,A,B,C,D,0.0D0,0.0D0,K)
      CALL SPLINE(X,Z,CA,CB,CC,CD,0.0D0,0.0D0,K)
!
      K=K-1
      DO 3 J=1,K
      IO=IO+1
      RG(IO)=X(J)
      RVG(IO)=Y(J)
      RWG(IO)=Z(J)
      VA(IO)=A(J)
      VB(IO)=B(J)
      VC(IO)=C(J)
      VD(IO)=D(J)
      WA(IO)=CA(J)
      WB(IO)=CB(J)
      WC(IO)=CC(J)
      WD(IO)=CD(J)
    3 CONTINUE
      IF(I.LT.NV) THEN
        K=0
        GO TO 1
      ENDIF
!  ****  AN EXTRA POINT IS ADDED TO THE GRID, AND R*V(R) IS SET
!        EQUAL TO RV(NV) FOR R.GE.R(NV)
      IO=IO+1
      RG(IO)=X(K+1)
      RVG(IO)=Y(K+1)
      RWG(IO)=Z(K+1)
      VA(IO)=RVG(IO)
      VB(IO)=0.0D0
      VC(IO)=0.0D0
      VD(IO)=0.0D0
      WA(IO)=RWG(IO)
      WB(IO)=0.0D0
      WC(IO)=0.0D0
      WD(IO)=0.0D0
      NVT=IO+1
      RG(NVT)=2.0D0*RG(IO)
      RVG(NVT)=RVG(IO)
      RWG(NVT)=RWG(IO)
      VA(NVT)=RVG(IO)
      VB(NVT)=0.0D0
      VC(NVT)=0.0D0
      VD(NVT)=0.0D0
      WA(NVT)=RWG(IO)
      WB(NVT)=0.0D0
      WC(NVT)=0.0D0
      WD(NVT)=0.0D0
      RETURN
      END
!  *********************************************************************
!                         SUBROUTINE ZDPHAS
!  *********************************************************************
      SUBROUTINE ZDPHAS(E,EPS,ZPHASE,K,IER)
!
!   This subroutine computes Dirac phase shifts for complex potentials.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER (KIND=RegInt_K) (I-N),COMPLEX (KIND=RegReal_K) (Z)
      PARAMETER (NDIM=1000,NPPG=NDIM+1,NPTG=NDIM+NPPG)
      PARAMETER (SL=137.03599976D0)  ! Speed of light (1/alpha)
      PARAMETER (PI=3.1415926535897932D0,PIH=0.5D0*PI)
      COMMON/RGRID/RRR(NPTG),P(NPTG),Q(NPTG),INDD(NPTG),NRTT
      COMMON/ZRGRID/R(NPTG),ZP(NPTG),ZQ(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),VD(NPPG),NVT
      COMMON/VGRIDI/RWG(NPPG),WA(NPPG),WB(NPPG),WC(NPPG),WD(NPPG)
      DIMENSION ZPA(NPTG),ZQA(NPTG),ZPB(NPTG),ZQB(NPTG)
      COMMON/OCOUL/WAVNUM,ETA,DELTA
!      EXTERNAL BESJN
!
      NRT=NRTT
      DO I=1,NRT
        R(I)=RRR(I)
        IND(I)=INDD(I)
      ENDDO
!
      ZI=DCMPLX(0.0D0,1.0D0)
      IER=0
      IF(K.EQ.0) THEN
        WRITE(6,2001)
 2001   FORMAT(1X,'*** ERROR IN ZDPHAS: K.EQ.0.')
        STOP
      ENDIF
!
      IF(E.LE.0.0D0) THEN
        IER=7
        WRITE(6,2002)
 2002   FORMAT(1X,'*** ERROR 7  (IN ZDPHAS): E.LE.0.')
        RETURN
      ENDIF
      EPSCUT=1.0D-9
!  ****  Orbital angular momentum quantum number.
      IF(K.LT.0) THEN
        L=-K-1
        KSIGN=1
      ELSE
        L=K
        KSIGN=-1
      ENDIF
      FL1=0.5D0*L*(L+1)
      RK=DSQRT(E*(E+2.0D0*SL*SL))/SL
!
!  ****  Asymptotic solution.
!
      RZINF=RV(NVT)
      IF(DABS(RZINF).LT.EPS) THEN
!  ****  Finite range fields.
        FACTOR=DSQRT(E/(E+2.0D0*SL*SL))
        ILAST=NRT+1
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          ZRVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))+(WA(INJ)+RN*(WB(INJ)+RN*(WC(INJ)+RN*WD(INJ))))*ZI
          T=EPS*RN*DABS(E*RN-FL1/RN)
          X=RK*RN
          IF(ABS(ZRVN).GT.T) GO TO 1
          BNL=BESJN_MOD(2,L,X)
          IF(DABS(BNL).GT.100.0D0) GO TO 1
          BNL1=BESJN_MOD(2,L+KSIGN,X)
          IF(DABS(BNL1).GT.100.0D0) GO TO 1
          BJL=BESJN_MOD(1,L,X)
          BJL1=BESJN_MOD(1,L+KSIGN,X)
          ILAST=IL
          ZPA(ILAST)=X*BJL
          ZPB(ILAST)=-X*BNL
          ZQA(ILAST)=-FACTOR*KSIGN*X*BJL1
          ZQB(ILAST)=FACTOR*KSIGN*X*BNL1
        ENDDO
    1   CONTINUE
        IF(ILAST.EQ.NRT+1) THEN
          IER=8
          WRITE(6,2003)
 2003   FORMAT(1X,'*** ERROR 8  (IN ZDPHAS): RAD(NGP) TOO SMALL.'/5X,'(EXTEND THE GRID TO LARGER RADII).')
          RETURN
        ENDIF
      ELSE
!  ****  Coulomb fields.
        TAS=MAX(1.0D-11,EPS)*DABS(RZINF)
        ILAST=NRT+1
        DO I=4,NRT
          IL=ILAST-1
          RN=R(IL)
          INJ=IND(IL)
          ZRVN=VA(INJ)+RN*(VB(INJ)+RN*(VC(INJ)+RN*VD(INJ)))+(WA(INJ)+RN*(WB(INJ)+RN*(WC(INJ)+RN*WD(INJ))))*ZI
          IF(ABS(ZRVN-RZINF).GT.TAS) GO TO 2
          CALL DCOUL(RZINF,E,K,RN,P0,Q0,P1,Q1,ERR)
          IF(ERR.GT.EPSCUT.OR.DABS(P1).GT.100.0D0) GO TO 2
          ILAST=IL
          ZPA(ILAST)=P0
          ZPB(ILAST)=P1
          ZQA(ILAST)=Q0
          ZQB(ILAST)=Q1
        ENDDO
    2   CONTINUE
        IF(ILAST.EQ.NRT+1) THEN
          IER=8
          WRITE(6,2003)
          RETURN
        ENDIF
      ENDIF
!
!  ****  Outward solution of the radial equation.
!
      CALL ZDOUTW(E,EPS,K,1,NZERO,ILAST)
!
!  ****  Phase shift.
!
      RM=R(ILAST)
      IL=IND(ILAST-1)
      ZVF=VA(IL)/RM+VB(IL)+RM*(VC(IL)+RM*VD(IL))+(WA(IL)/RM+WB(IL)+RM*(WC(IL)+RM*WD(IL)))*ZI
      ZFG=(E-ZVF+2.0D0*SL*SL)/SL
      ZPO=ZP(ILAST)
      ZPOP=-K*ZPO/RM+ZFG*ZQ(ILAST)
      IL=IND(ILAST)
      ZVF=VA(IL)/RM+VB(IL)+RM*(VC(IL)+RM*VD(IL))+(WA(IL)/RM+WB(IL)+RM*(WC(IL)+RM*WD(IL)))*ZI
      ZFG=(E-ZVF+2.0D0*SL*SL)/SL
      ZPIA=ZPA(ILAST)
      ZPIAP=-K*ZPIA/RM+ZFG*ZQA(ILAST)
      ZPIB=ZPB(ILAST)
      ZPIBP=-K*ZPIB/RM+ZFG*ZQB(ILAST)
!
      ZRATIO=ZPOP/ZPO
      ZPHASE=(ZRATIO*(ZPIA+ZI*ZPIB)-(ZPIAP+ZI*ZPIBP))/((ZPIAP-ZI*ZPIBP)-ZRATIO*(ZPIA-ZI*ZPIB))
!
      RETURN
      END
!  **************************************************************
!                       SUBROUTINE ZDOUTW
!  **************************************************************
      SUBROUTINE ZDOUTW(E,EPS,K,NR,NZERO,IOTP)
!
!     OUTWARD SOLUTION OF THE DIRAC RADIAL EQUATION FOR A COMPLEX
!  PIECEWISE CUBIC POTENTIAL. POWER SERIES METHOD.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER (KIND=RegInt_K) (I-N),COMPLEX (KIND=RegReal_K) (Z)
      PARAMETER (NDIM=1000,NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/ZRADWF/RAD(NDIM),ZPIO(NDIM),ZQIO(NDIM),NGP,ILAST,IER
      COMMON/ZRGRID/R(NPTG),ZP(NPTG),ZQ(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),VD(NPPG),NVT
      COMMON/VGRIDI/RWG(NPPG),WA(NPPG),WB(NPPG),WC(NPPG),WD(NPPG)
      COMMON/ZPOTEN/RV0,RV1,RV2,RV3,RW0,RW1,RW2,RW3
      COMMON/ZDINOU/ZPI,ZQI,ZPF,ZQF,RA,RB,RLN,NSTEP,NCHS
      COMMON/ZDSAVE/ZP0,ZQ0,ZP1,ZQ1,ZCA(60),ZCB(60),R0,R1,NSUM
      COMMON/NZT/NZMAX
      NZERO=0
      NZMAX=0
      AK=K
      IF(E.LT.0.0D0) THEN

        N1=NRT
      ELSE
        N1=IOTP-1
      ENDIF
!
      ZP(1)=0.0D0
      ZQ(1)=0.0D0
      DO 2 I=1,N1
      RA=R(I)
      RB=R(I+1)
      IN=IND(I)
      RV0=VA(IN)
      RV1=VB(IN)
      RV2=VC(IN)
      RV3=VD(IN)
      RW0=WA(IN)
      RW1=WB(IN)
      RW2=WC(IN)
      RW3=WD(IN)
      ZPI=ZP(I)
      ZQI=ZQ(I)
      CALL ZDIR(E,AK,EPS)
      NZERO=NZERO+NCHS
      IF(NCHS.GT.NZMAX) NZMAX=NCHS
      IF(NZERO.GT.NR.AND.E.LT.0.0D0) RETURN
      ZP(I+1)=ZPF
      ZQ(I+1)=ZQF
      IF(E.LT.0.0D0) THEN
!  ****  TCONV IS THE PRODUCT OF P AND ITS SECOND DERIVATIVE AT
!        THE I-TH GRID POINT (POSITIVE IF P IS CONVEX).
        TCONV=dble(2.0D0*ZCA(3)*ZPI)
        IF(I.GE.IOTP.AND.TCONV.GT.1.0D-15) THEN
          IOTP=I+1
          RETURN
        ENDIF
      ENDIF
      IF(I.EQ.1) GO TO 2
!  ****  RENORMALIZATION.
      IF(RLN.GT.0.0D0) THEN
      FACT=DEXP(-RLN)
      DO 1 J=1,I
      ZP(J)=ZP(J)*FACT
      ZQ(J)=ZQ(J)*FACT
    1 CONTINUE
      ENDIF
    2 CONTINUE
      RETURN
      END
!  **************************************************************
!                       SUBROUTINE ZDIR
!  **************************************************************
      SUBROUTINE ZDIR(E,AK,EPS)
!
!  THIS SUBROUTINE SOLVES THE RADIAL DIRAC EQUATION FOR A COMPLEX
!  CENTRAL FIELD V(R) SUCH THAT
!              R*V(R) = RV0+RV1*R+RV2*R**2+RV3*R**3
!                     +ZI*(RW0+RW1*R+RW2*R**2+RW3*R**3)
!     GIVEN THE BOUNDARY CONDITIONS (I.E. THE VALUE OF THE
!  LARGE AND SMALL RADIAL FUNCTIONS) AT RA, THE SOLUTION IN
!  THE INTERVAL BETWEEN RA AND RB IS GENERATED BY USING A
!  PIECEWISE POWER SERIES EXPANSION FOR A PARTITION OF THE
!  INTERVAL, SUITABLY CHOSEN TO ALLOW FAST CONVERGENCE OF THE
!  SERIES.
!
!   INPUT ARGUMENTS:
!      E ..................... PARTICLE KINETIC ENERGY
!      AK .................... RELATIVISTIC ANGULAR MOMENTUM
!                              QUANTUM NUMBER
!
!   INPUT (COMMON POTEN):
!      RV0, RV1, RV2, RV3 .... REAL POTENTIAL PARAMETERS
!      RW0, RW1, RW2, RW3 .... IMAGINARY POTENTIAL PARAMETERS
!
!   INPUT-OUTPUT (COMMON DINOUT):
!      RA, RB ................ INTERVAL END POINTS (INPUT)
!      ZPI, ZQI ................ VALUES OF THE LARGE AND SMALL
!                              RADIAL FUNCTIONS AT RA (INPUT)
!      ZPF, ZQF ................ VALUES OF THE LARGE AND SMALL
!                              RADIAL FUNCTIONS AT RB (OUTPUT)
!      RLN ................... DLOG OF THE RE-NORMALIZING FACTOR
!      EPS ................... ESTIMATE OF THE GLOBAL ERROR IN
!                              PF AND QF
!      NSTEP ................. NUMBER OF STEPS

!      NCHS .................. NUMBER OF ZEROS OF P(R) IN (RA,RB)
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER (KIND=RegInt_K) (I-N),COMPLEX (KIND=RegReal_K) (Z)
      COMMON/ZPOTEN/RV0,RV1,RV2,RV3,RW0,RW1,RW2,RW3
      COMMON/ZDINOU/ZPI,ZQI,ZPF,ZQF,RA,RB,RLN,NSTEP,NCHS
      COMMON/ZDSAVE/ZP0,ZQ0,ZP1,ZQ1,ZCA(60),ZCB(60),R0,R1,NSUM
      NCHS=0
      RLN=0.0D0
!
      H=RB-RA
      IF(H.LT.0.0D0) THEN
      DIRECT=-1.0D0
      ELSE
      DIRECT=1.0D0
      ENDIF
      K=-2
      NSTEP=0
!
      R1=RA
      ZP1=ZPI
      ZQ1=ZQI
    1 R0=R1
      ZP0=ZP1
      ZQ0=ZQ1
    2 IOUT=0
      R1=R0+H
      IF(DIRECT*(RB-R1).LT.DIRECT*1.0D-1*H) THEN
      R1=RB
      H=RB-R0
      IOUT=1
      ENDIF
      CALL ZDIR0(E,AK,EPS)
!
      K=K+1
      IF(NSUM.GT.15) GO TO 3
      IF(K.LT.0) GO TO 4
      H=H+H
      K=0
      GO TO 4
    3 IF(NSUM.LT.60) GO TO 4
      H=0.5D0*H
      K=-4
      GO TO 2
    4 NSTEP=NSTEP+1
      TST=ABS(ZP1)
      IF(TST.GT.1.0D2) THEN
!  ****  RENORMALIZATION.
      RLN=RLN+DLOG(TST)
      ZP1=ZP1/TST
      ZQ1=ZQ1/TST
      ENDIF
      TSTN=dble(ZP0*ZP1)
      IF(TSTN.LT.0.0D0.AND.R0.GT.0.0D0) NCHS=NCHS+1
      IF(IOUT.EQ.0) GO TO 1
!  ****  OUTPUT.
      ZPF=ZP1
      ZQF=ZQ1
      RETURN
      END
!  **************************************************************
!                       SUBROUTINE ZDIR0
!  **************************************************************
      SUBROUTINE ZDIR0(E,AK,EPS)
!
!  Power series solution of the Dirac eq. for a central potential
!  with an imaginary component (negative for absorptive interac-
!  tions).
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Y), INTEGER (KIND=RegInt_K) (I-N),COMPLEX (KIND=RegReal_K) (Z)
!  ****  SPEED OF LIGHT AND OVERFLOW LEVEL.
      PARAMETER (SL=137.03599976D0)  ! Speed of light (1/alpha)
      PARAMETER (OVER=1.0D15)
      COMMON/ZPOTEN/RV0,RV1,RV2,RV3,RW0,RW1,RW2,RW3
      COMMON/ZDSAVE/ZP0,ZQ0,ZP1,ZQ1,ZCA(60),ZCB(60),R0,R1,NSUM
!
      ZI=DCMPLX(0.0D0,1.0D0)
!
      ISIG=1
      IF(AK.GT.0.0D0) ISIG=-1
      H=R1-R0
      H2=H*H
      ZRVE=RV1+ZI*RW1-E
      ZRV0=RV0+ZI*RW0
      ZRV1=RV1+ZI*RW1
      ZRV2=RV2+ZI*RW2
      ZRV3=RV3+ZI*RW3
!
      IF(R0.GT.1.0D-10) GO TO 7
!
!  **** FIRST INTERVAL.
!
      ZU0=ZRV0/SL
      ZU1=ZRVE*R1/SL
      ZU2=ZRV2*R1**2/SL
      ZU3=ZRV3*R1**3/SL
      ZUT=ZU0+ZU1+ZU2+ZU3
      ZUQ=ZUT-2*SL*R1
      ZUH=ZU1-2*SL*R1
      IF(ABS(ZU0).LT.1.0D-10) GO TO 2
!
!  ****  U0.NE.0.
      ZS=SQRT(AK*AK-ZU0*ZU0)
      ZDS=ZS+ZS
      ZCA(1)=1.0D0
      ZCB(1)=-(ZS+AK)/ZU0
      ZCAI=ZU1*ZCA(1)
      ZCBI=ZUH*ZCB(1)
      ZCA(2)=(-ZU0*ZCAI-(ZS+1-AK)*ZCBI)/(ZDS+1)
      ZCB(2)=((ZS+1+AK)*ZCAI-ZU0*ZCBI)/(ZDS+1)
      ZCAI=ZU1*ZCA(2)+ZU2*ZCA(1)
      ZCBI=ZUH*ZCB(2)+ZU2*ZCB(1)
      ZCA(3)=(-ZU0*ZCAI-(ZS+2-AK)*ZCBI)/(2*(ZDS+2))
      ZCB(3)=((ZS+2+AK)*ZCAI-ZU0*ZCBI)/(2*(ZDS+2))
      ZP1=ZCA(1)+ZCA(2)+ZCA(3)
      ZPP1=ZS*ZCA(1)+(ZS+1)*ZCA(2)+(ZS+2)*ZCA(3)
      ZQ1=ZCB(1)+ZCB(2)+ZCB(3)
      ZQP1=ZS*ZCB(1)+(ZS+1)*ZCB(2)+(ZS+2)*ZCB(3)
!
      DO 1 I=4,60
      K=I-1
      ZCAI=ZU1*ZCA(K)+ZU2*ZCA(I-2)+ZU3*ZCA(I-3)
      ZCBI=ZUH*ZCB(K)+ZU2*ZCB(I-2)+ZU3*ZCB(I-3)
      ZCA(I)=(-ZU0*ZCAI-(ZS+K-AK)*ZCBI)/(K*(ZDS+K))
      ZCB(I)=((ZS+K+AK)*ZCAI-ZU0*ZCBI)/(K*(ZDS+K))
      ZP1=ZP1+ZCA(I)
      ZPP1=ZPP1+(ZS+K)*ZCA(I)
      ZQ1=ZQ1+ZCB(I)
      ZQP1=ZQP1+(ZS+K)*ZCB(I)
!  ****  CHECK OVERFLOW LIMIT.
      TST=DMAX1(ABS(ZP1),ABS(ZQ1),ABS(ZPP1),ABS(ZQP1))
      IF(TST.GT.OVER) THEN
      NSUM=100
      RETURN
      ENDIF
      T1A=ABS(R1*ZPP1+H*(AK*ZP1+ZUQ*ZQ1))
      T1B=ABS(R1*ZQP1-H*(AK*ZQ1+ZUT*ZP1))
      T1=MAX(T1A,T1B)
      T2=MAX(ABS(ZCA(I)),ABS(ZCB(I)))
      TST=EPS*MAX(ABS(ZP1),ABS(ZQ1))
      IF(T1.LT.TST.AND.T2.LT.TST) GO TO 6
    1 CONTINUE
      GO TO 6
!
!  ****  ZU0.EQ.0 AND SIG=1.
    2 IF(ISIG.LT.0) GO TO 4
      ZS=ABS(AK)
      ZDS1=ZS+ZS+1
      ZCA(1)=1.0D0
      ZCB(1)=ZU1*ZCA(1)/ZDS1
      ZCA(2)=0.0D0
      ZCB(2)=ZU2*ZCA(1)/(ZDS1+1)
      ZCA(3)=-ZUH*ZCB(1)/2
      ZCB(3)=(ZU1*ZCA(3)+ZU3*ZCA(1))/(ZDS1+2)
      ZCA(4)=-(ZUH*ZCB(2)+ZU2*ZCB(1))/3
      ZCB(4)=(ZU1*ZCA(4)+ZU2*ZCA(3))/(ZDS1+3)
      ZP1=ZCA(1)+ZCA(2)+ZCA(3)+ZCA(4)
      ZPP1=ZS*ZCA(1)+(ZS+1)*ZCA(2)+(ZS+2)*ZCA(3)+(ZS+3)*ZCA(4)
      ZQ1=ZCB(1)+ZCB(2)+ZCB(3)+ZCB(4)
      ZQP1=(ZS+1)*ZCB(1)+(ZS+2)*ZCB(2)+(ZS+3)*ZCB(3)
!
      DO 3 I=5,60
      K=I-1
      ZCA(I)=-(ZUH*ZCB(I-2)+ZU2*ZCB(I-3)+ZU3*ZCB(I-4))/K
      ZCB(I)=(ZU1*ZCA(I)+ZU2*ZCA(K)+ZU3*ZCA(I-2))/(ZDS1+K)
      ZP1=ZP1+ZCA(I)
      ZPP1=ZPP1+(ZS+K)*ZCA(I)
      ZQ1=ZQ1+ZCB(I)
      ZQP1=ZQP1+(ZS+I)*ZCB(I)
!  ****  CHECK OVERFLOW LIMIT.
      TST=MAX(ABS(ZP1),ABS(ZQ1),ABS(ZPP1),ABS(ZQP1))
      IF(TST.GT.OVER) THEN
      NSUM=100
      RETURN
      ENDIF
      T1A=ABS(R1*ZPP1+H*(AK*ZP1+ZUQ*ZQ1))
      T1B=ABS(R1*ZQP1-H*(AK*ZQ1+ZUT*ZP1))
      T1=MAX(T1A,T1B)
      T2=MAX(ABS(ZCA(I)),ABS(ZCB(I)))
      TST=EPS*MAX(ABS(ZP1),ABS(ZQ1))
      IF(T1.LT.TST.AND.T2.LT.TST) GO TO 6
    3 CONTINUE
      GO TO 6
!
!  ****  ZU0.EQ.0 AND SIG=-1.
    4 S=DABS(AK)+1
      DS1=S+DABS(AK)
      RZUH=dble(ZUH)
      IF(RZUH.GT.0.0D0) THEN
      ZCB(1)=-1.0D0
      ELSE
      ZCB(1)=1.0D0
      ENDIF
      ZCA(1)=-ZUH*ZCB(1)/DS1
      ZCB(2)=0.0D0
      ZCA(2)=-ZU2*ZCB(1)/(DS1+1)
      ZCB(3)=ZU1*ZCA(1)/2
      ZCA(3)=-(ZUH*ZCB(3)+ZU3*ZCB(1))/(DS1+2)
      ZCB(4)=(ZU1*ZCA(2)+ZU2*ZCA(1))/3
      ZCA(4)=-(ZUH*ZCB(4)+ZU2*ZCB(3))/(DS1+3)
      ZP1=ZCA(1)+ZCA(2)+ZCA(3)+ZCA(4)
      ZPP1=S*ZCA(1)+(S+1)*ZCA(2)+(S+2)*ZCA(3)+(S+3)*ZCA(4)
      ZQ1=ZCB(1)+ZCB(2)+ZCB(3)+ZCB(4)
      ZQP1=(S-1)*ZCB(1)+S*ZCB(2)+(S+1)*ZCB(3)
!
      DO 5 I=5,60
      K=I-1
      ZCB(I)=(ZU1*ZCA(I-2)+ZU2*ZCA(I-3)+ZU3*ZCA(I-4))/K
      ZCA(I)=-(ZUH*ZCB(I)+ZU2*ZCB(K)+ZU3*ZCB(I-2))/(DS1+K)
      ZP1=ZP1+ZCA(I)
      ZPP1=ZPP1+(S+K)*ZCA(I)
      ZQ1=ZQ1+ZCB(I)
      ZQP1=ZQP1+(S+K-1)*ZCB(I)
!  ****  CHECK OVERFLOW LIMIT.
      TST=MAX(ABS(ZP1),ABS(ZQ1),ABS(ZPP1),ABS(ZQP1))
      IF(TST.GT.OVER) THEN
      NSUM=100
      RETURN
      ENDIF
      T1A=ABS(R1*ZPP1+H*(AK*ZP1+ZUQ*ZQ1))
      T1B=ABS(R1*ZQP1-H*(AK*ZQ1+ZUT*ZP1))
      T1=MAX(T1A,T1B)
      T2=MAX(ABS(ZCA(I)),ABS(ZCB(I)))
      TST=EPS*MAX(ABS(ZP1),ABS(ZQ1))
      IF(T1.LT.TST.AND.T2.LT.TST) GO TO 6
    5 CONTINUE
!  ****  RENORMALIZATION.
    6 NSUM=K+1
      ZQ1=ZQ1/ABS(ZP1)
      ZP1=ZP1/ABS(ZP1)
      RETURN
!
!  **** MIDDLE REGION.
!
    7 CONTINUE
      RHO=H/R0
      ZU0=(ZRV0+R0*(ZRVE+R0*(ZRV2+R0*ZRV3)))/SL
      ZU1=(ZRVE+R0*(2*ZRV2+R0*3*ZRV3))*H/SL
      ZU2=(ZRV2+R0*3*ZRV3)*H2/SL
      ZU3=ZRV3*H*H2/SL
      ZUB=ZU0-2*SL*R0
      ZUH=ZU1-2*SL*H
      ZUT=ZU0+ZU1+ZU2+ZU3
      ZUQ=ZUT-2*SL*R1
!
      ZCA(1)=ZP0
      ZCB(1)=ZQ0
      ZCA(2)=-RHO*(AK*ZCA(1)+ZUB*ZCB(1))
      ZCB(2)=RHO*(AK*ZCB(1)+ZU0*ZCA(1))
      ZCA(3)=-RHO*((AK+1)*ZCA(2)+ZUB*ZCB(2)+ZUH*ZCB(1))/2
      ZCB(3)=RHO*((AK-1)*ZCB(2)+ZU0*ZCA(2)+ZU1*ZCA(1))/2
      ZCA(4)=-RHO*((AK+2)*ZCA(3)+ZUB*ZCB(3)+ZUH*ZCB(2)+ZU2*ZCB(1))/3
      ZCB(4)=RHO*((AK-2)*ZCB(3)+ZU0*ZCA(3)+ZU1*ZCA(2)+ZU2*ZCA(1))/3
!
      ZP1=ZCA(1)+ZCA(2)+ZCA(3)+ZCA(4)
      ZPP1=ZCA(2)+2*ZCA(3)+3*ZCA(4)
      ZQ1=ZCB(1)+ZCB(2)+ZCB(3)+ZCB(4)
      ZQP1=ZCB(2)+2*ZCB(3)+3*ZCB(4)
!
      DO 9 I=5,60
      K=I-1
      ZCA(I)=-RHO*((AK+K-1)*ZCA(K)+ZUB*ZCB(K)+ZUH*ZCB(I-2)+ZU2*ZCB(I-3)+ZU3*ZCB(I-4))/K
      ZCB(I)=RHO*((AK-K+1)*ZCB(K)+ZU0*ZCA(K)+ZU1*ZCA(I-2)+ZU2*ZCA(I-3)+ZU3*ZCA(I-4))/K
      ZP1=ZP1+ZCA(I)
      ZPP1=ZPP1+K*ZCA(I)
      ZQ1=ZQ1+ZCB(I)
      ZQP1=ZQP1+K*ZCB(I)
!  ****  CHECK OVERFLOW LIMIT.
      TST=MAX(ABS(ZP1),ABS(ZQ1),ABS(ZPP1),ABS(ZQP1))
      IF(TST.GT.OVER) THEN
      NSUM=100
      RETURN
      ENDIF
      T1A=ABS(R1*ZPP1+H*(AK*ZP1+ZUQ*ZQ1))
      T1B=ABS(R1*ZQP1-H*(AK*ZQ1+ZUT*ZP1))
      T1=MAX(T1A,T1B)
      T2=MAX(ABS(ZCA(I)),ABS(ZCB(I)))
      TST=EPS*MAX(ABS(ZP1),ABS(ZQ1))
      IF(T1.LT.TST.AND.T2.LT.TST) GO TO 10
    9 CONTINUE
!
   10 NSUM=K+1
      RETURN
      END
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! Subroutines from the RADIAL package (modified to follow Roses's phase
! convention). Here all quantities are in atomic units.
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!  **************************************************************
!                        SUBROUTINE VINT
!  **************************************************************
      SUBROUTINE VINT(R,RV,NV)
!
!     NATURAL CUBIC SPLINE INTERPOLATION FOR R*V(R) FROM THE
!  INPUT RADII AND POTENTIAL VALUES (128).
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER (NDIM=1000,NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/VGRID/RG(NPPG),RVG(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),VD(NPPG),NVT
      COMMON/RGRID/X(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT
      COMMON/STORE/Y(NPTG),A(NPTG),B(NPTG),C(NPTG),D(NPTG)
      DIMENSION R(NDIM),RV(NDIM)
!
      IF(R(1).LT.0.0D0) THEN
      WRITE(6,2101)
 2101 FORMAT(1X,'*** ERROR IN VINT: R(1).LT.0.')
      STOP
      ENDIF
      IF(NV.GT.NDIM) THEN
      WRITE(6,2102) NDIM
 2102 FORMAT(1X,'*** ERROR IN VINT: INPUT POTENTIAL GRID WITH MORE THAN ',I5,' DATA POINTS.')
      STOP
      ENDIF
      R(1)=0.0D0
!
      IO=0
      I=0
      K=0
    1 I=I+1
      K=K+1
      X(K)=R(I)
      Y(K)=RV(I)
      IF(I.EQ.NV) GO TO 2
      IF(R(I).LT.R(I+1)-1.0D-12) GO TO 1
    2 CONTINUE
!
      CALL SPLINE(X,Y,A,B,C,D,0.0D0,0.0D0,K)
!
      K=K-1
      DO 3 J=1,K
      IO=IO+1
      RG(IO)=X(J)
      RVG(IO)=Y(J)
      VA(IO)=A(J)
      VB(IO)=B(J)
      VC(IO)=C(J)
      VD(IO)=D(J)
    3 CONTINUE
      IF(I.LT.NV) THEN
        K=0
        GO TO 1
      ENDIF
!  ****  AN EXTRA POINT IS ADDED TO THE GRID, AND R*V(R) IS SET
!        EQUAL TO RV(NV) FOR R.GE.R(NV)
      IO=IO+1
      RG(IO)=X(K+1)
      RVG(IO)=Y(K+1)
      VA(IO)=RVG(IO)
      VB(IO)=0.0D0
      VC(IO)=0.0D0
      VD(IO)=0.0D0
      NVT=IO+1
      RG(NVT)=2.0D0*RG(IO)
      RVG(NVT)=RVG(IO)
      VA(NVT)=RVG(IO)
      VB(NVT)=0.0D0
      VC(NVT)=0.0D0
      VD(NVT)=0.0D0
      RETURN
      END
!  **************************************************************
!                       SUBROUTINE DOUTW
!  **************************************************************
      SUBROUTINE DOUTW(E,EPS,K,NR,NZERO,IOTP)
!
!     OUTWARD SOLUTION OF THE DIRAC RADIAL EQUATION FOR A PIECE-
!   WISE CUBIC FIELD. POWER SERIES METHOD.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER (NDIM=1000,NPPG=NDIM+1,NPTG=NDIM+NPPG)
      COMMON/RADWF/RAD(NDIM),PIO(NDIM),QIO(NDIM),NGP,ILAST,IER
      COMMON/RGRID/R(NPTG),P(NPTG),Q(NPTG),IND(NPTG),NRT
      COMMON/VGRID/RG(NPPG),RV(NPPG),VA(NPPG),VB(NPPG),VC(NPPG),VD(NPPG),NVT
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/DINOUT/PI,QI,PF,QF,RA,RB,RLN,NSTEP,NCHS
      COMMON/DSAVE/P0,Q0,P1,Q1,CA(60),CB(60),R0,R1,NSUM
      COMMON/NZT/NZMAX
      NZERO=0
      NZMAX=0
      AK=K
      IF(E.LT.0.0D0) THEN
        N1=NRT
      ELSE
        N1=IOTP-1
      ENDIF
!
      P(1)=0.0D0
      Q(1)=0.0D0
      DO 2 I=1,N1
      RA=R(I)
      RB=R(I+1)
      IN=IND(I)
      RV0=VA(IN)
      RV1=VB(IN)
      RV2=VC(IN)
      RV3=VD(IN)
      PI=P(I)
      QI=Q(I)
      CALL DIR(E,AK,EPS)
      NZERO=NZERO+NCHS
      IF(NCHS.GT.NZMAX) NZMAX=NCHS
      IF(NZERO.GT.NR.AND.E.LT.0.0D0) RETURN
      P(I+1)=PF
      Q(I+1)=QF
      IF(E.LT.0.0D0) THEN
!  ****  TCONV IS THE PRODUCT OF P AND ITS SECOND DERIVATIVE AT
!        THE I-TH GRID POINT (POSITIVE IF P IS CONVEX).
        TCONV=2.0D0*CA(3)*PI
        IF(I.GE.IOTP.AND.TCONV.GT.1.0D-15) THEN
          IOTP=I+1
          RETURN
        ENDIF
      ENDIF
      IF(I.EQ.1) GO TO 2
!  ****  RENORMALIZATION.
      IF(RLN.GT.0.0D0) THEN
      FACT=DEXP(-RLN)
      DO 1 J=1,I
      P(J)=P(J)*FACT
      Q(J)=Q(J)*FACT
    1 CONTINUE
      ENDIF
    2 CONTINUE
      RETURN
      END
!  **************************************************************
!                       SUBROUTINE DIR
!  **************************************************************
      SUBROUTINE DIR(E,AK,EPS)
!
!      THIS SUBROUTINE SOLVES THE RADIAL DIRAC EQUATION FOR A
!   CENTRAL FIELD V(R) SUCH THAT
!              R*V(R) = RV0+RV1*R+RV2*R**2+RV3*R**3
!      GIVEN THE BOUNDARY CONDITIONS (I.E. THE VALUE OF THE
!   LARGE AND SMALL RADIAL FUNCTIONS) AT RA, THE SOLUTION IN
!   THE INTERVAL BETWEEN RA AND RB IS GENERATED BY USING A
!   PIECEWISE POWER SERIES EXPANSION FOR A PARTITION OF THE
!   INTERVAL, SUITABLY CHOSEN TO ALLOW FAST CONVERGENCE OF THE
!   SERIES.
!
!   INPUT ARGUMENTS:
!      E ..................... PARTICLE KINETIC ENERGY
!      AK .................... RELATIVISTIC ANGULAR MOMENTUM
!                              QUANTUM NUMBER
!
!   INPUT (COMMON POTEN):
!      RV0, RV1, RV2, RV3 .... POTENTIAL PARAMETERS
!
!   INPUT-OUTPUT (COMMON DINOUT):
!      RA, RB ................ INTERVAL END POINTS (INPUT)
!      PI, QI ................ VALUES OF THE LARGE AND SMALL
!                              RADIAL FUNCTIONS AT RA (INPUT)
!      PF, QF ................ VALUES OF THE LARGE AND SMALL
!                              RADIAL FUNCTIONS AT RB (OUTPUT)
!      RLN ................... DLOG OF THE RE-NORMALIZING FACTOR
!      EPS ................... ESTIMATE OF THE GLOBAL ERROR IN
!                              PF AND QF
!      NSTEP ................. NUMBER OF STEPS
!      NCHS .................. NUMBER OF ZEROS OF P(R) IN (RA,RB)
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/DINOUT/PI,QI,PF,QF,RA,RB,RLN,NSTEP,NCHS
      COMMON/DSAVE/P0,Q0,P1,Q1,CA(60),CB(60),R0,R1,NSUM
      NCHS=0
      RLN=0.0D0
!
      H=RB-RA
      IF(H.LT.0.0D0) THEN
      DIRECT=-1.0D0
      ELSE
      DIRECT=1.0D0
      ENDIF
      K=-2
      NSTEP=0
!
      R1=RA
      P1=PI
      Q1=QI
    1 R0=R1
      P0=P1
      Q0=Q1
    2 IOUT=0
      R1=R0+H
      IF(DIRECT*(RB-R1).LT.DIRECT*1.0D-1*H) THEN
      R1=RB
      H=RB-R0
      IOUT=1
      ENDIF
      CALL DIR0(E,AK,EPS)
!
      K=K+1
      IF(NSUM.GT.15) GO TO 3
      IF(K.LT.0) GO TO 4
      H=H+H
      K=0
      GO TO 4
    3 IF(NSUM.LT.60) GO TO 4
      H=0.5D0*H
      K=-4
      GO TO 2
    4 NSTEP=NSTEP+1
      TST=DABS(P1)
      IF(TST.GT.1.0D2) THEN
!  ****  RENORMALIZATION.
      RLN=RLN+DLOG(TST)
      P1=P1/TST
      Q1=Q1/TST
      ENDIF
      IF(P0*P1.LT.0.0D0.AND.R0.GT.0.0D0) NCHS=NCHS+1
      IF(IOUT.EQ.0) GO TO 1
!  ****  OUTPUT.
      PF=P1
      QF=Q1
      RETURN
      END
!  **************************************************************
!                       SUBROUTINE DIR0
!  **************************************************************
      SUBROUTINE DIR0(E,AK,EPS)
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
!  ****  SPEED OF LIGHT AND OVERFLOW LEVEL.
      PARAMETER (SL=137.03599976D0)  ! Speed of light (1/alpha)
      PARAMETER (OVER=1.0D15)
      COMMON/POTEN/RV0,RV1,RV2,RV3
      COMMON/DSAVE/P0,Q0,P1,Q1,CA(60),CB(60),R0,R1,NSUM
!
      ISIG=1
      IF(AK.GT.0.0D0) ISIG=-1
      H=R1-R0
      H2=H*H
      RVE=RV1-E
!
      IF(R0.GT.1.0D-10) GO TO 7
!
!  **** FIRST INTERVAL. (147-164)
!
      U0=RV0/SL
      U1=RVE*R1/SL
      U2=RV2*R1**2/SL
      U3=RV3*R1**3/SL
      UT=U0+U1+U2+U3
      UQ=UT-2*SL*R1
      UH=U1-2*SL*R1
      IF(DABS(U0).LT.1.0D-10) GO TO 2
!
!  ****  U0.NE.0. (155-159)
      S=DSQRT(AK*AK-U0*U0)
      DS=S+S
      CA(1)=1.0D0
      CB(1)=-(S+AK)/U0
      CAI=U1*CA(1)
      CBI=UH*CB(1)
      CA(2)=(-U0*CAI-(S+1-AK)*CBI)/(DS+1)
      CB(2)=((S+1+AK)*CAI-U0*CBI)/(DS+1)
      CAI=U1*CA(2)+U2*CA(1)
      CBI=UH*CB(2)+U2*CB(1)
      CA(3)=(-U0*CAI-(S+2-AK)*CBI)/(2*(DS+2))
      CB(3)=((S+2+AK)*CAI-U0*CBI)/(2*(DS+2))
      P1=CA(1)+CA(2)+CA(3)
      PP1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)
      Q1=CB(1)+CB(2)+CB(3)
      QP1=S*CB(1)+(S+1)*CB(2)+(S+2)*CB(3)
!
      DO 1 I=4,60
      K=I-1
      CAI=U1*CA(K)+U2*CA(I-2)+U3*CA(I-3)
      CBI=UH*CB(K)+U2*CB(I-2)+U3*CB(I-3)
      CA(I)=(-U0*CAI-(S+K-AK)*CBI)/(K*(DS+K))
      CB(I)=((S+K+AK)*CAI-U0*CBI)/(K*(DS+K))
      P1=P1+CA(I)
      PP1=PP1+(S+K)*CA(I)
      Q1=Q1+CB(I)
      QP1=QP1+(S+K)*CB(I)
!  ****  CHECK OVERFLOW LIMIT.
      TST=DMAX1(DABS(P1),DABS(Q1),DABS(PP1),DABS(QP1))
      IF(TST.GT.OVER) THEN
      NSUM=100
      RETURN
      ENDIF
      T1A=DABS(R1*PP1+H*(AK*P1+UQ*Q1))
      T1B=DABS(R1*QP1-H*(AK*Q1+UT*P1))
      T1=DMAX1(T1A,T1B)
      T2=DMAX1(DABS(CA(I)),DABS(CB(I)))
      TST=EPS*DMAX1(DABS(P1),DABS(Q1))
      IF(T1.LT.TST.AND.T2.LT.TST) GO TO 6
    1 CONTINUE
      GO TO 6
!
!  ****  U0.EQ.0 AND SIG=1. (160,161)
    2 IF(ISIG.LT.0) GO TO 4
      S=DABS(AK)
      DS1=S+S+1
      CA(1)=1.0D0
      CB(1)=U1*CA(1)/DS1
      CA(2)=0.0D0
      CB(2)=U2*CA(1)/(DS1+1)
      CA(3)=-UH*CB(1)/2
      CB(3)=(U1*CA(3)+U3*CA(1))/(DS1+2)
      CA(4)=-(UH*CB(2)+U2*CB(1))/3
      CB(4)=(U1*CA(4)+U2*CA(3))/(DS1+3)
      P1=CA(1)+CA(2)+CA(3)+CA(4)
      PP1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)+(S+3)*CA(4)
      Q1=CB(1)+CB(2)+CB(3)+CB(4)
      QP1=(S+1)*CB(1)+(S+2)*CB(2)+(S+3)*CB(3)
!
      DO 3 I=5,60
      K=I-1
      CA(I)=-(UH*CB(I-2)+U2*CB(I-3)+U3*CB(I-4))/K
      CB(I)=(U1*CA(I)+U2*CA(K)+U3*CA(I-2))/(DS1+K)
      P1=P1+CA(I)
      PP1=PP1+(S+K)*CA(I)
      Q1=Q1+CB(I)
      QP1=QP1+(S+I)*CB(I)
!  ****  CHECK OVERFLOW LIMIT.
      TST=DMAX1(DABS(P1),DABS(Q1),DABS(PP1),DABS(QP1))
      IF(TST.GT.OVER) THEN
      NSUM=100
      RETURN
      ENDIF
      T1A=DABS(R1*PP1+H*(AK*P1+UQ*Q1))
      T1B=DABS(R1*QP1-H*(AK*Q1+UT*P1))
      T1=DMAX1(T1A,T1B)
      T2=DMAX1(DABS(CA(I)),DABS(CB(I)))
      TST=EPS*DMAX1(DABS(P1),DABS(Q1))
      IF(T1.LT.TST.AND.T2.LT.TST) GO TO 6
    3 CONTINUE
      GO TO 6
!
!  ****  U0.EQ.0 AND SIG=-1. (162,163)
    4 S=DABS(AK)+1
      DS1=S+DABS(AK)
      IF(UH.GT.0.0D0) THEN
      CB(1)=-1.0D0
      ELSE
      CB(1)=1.0D0
      ENDIF
      CA(1)=-UH*CB(1)/DS1
      CB(2)=0.0D0
      CA(2)=-U2*CB(1)/(DS1+1)
      CB(3)=U1*CA(1)/2
      CA(3)=-(UH*CB(3)+U3*CB(1))/(DS1+2)
      CB(4)=(U1*CA(2)+U2*CA(1))/3
      CA(4)=-(UH*CB(4)+U2*CB(3))/(DS1+3)
      P1=CA(1)+CA(2)+CA(3)+CA(4)
      PP1=S*CA(1)+(S+1)*CA(2)+(S+2)*CA(3)+(S+3)*CA(4)
      Q1=CB(1)+CB(2)+CB(3)+CB(4)
      QP1=(S-1)*CB(1)+S*CB(2)+(S+1)*CB(3)
!
      DO 5 I=5,60
      K=I-1
      CB(I)=(U1*CA(I-2)+U2*CA(I-3)+U3*CA(I-4))/K
      CA(I)=-(UH*CB(I)+U2*CB(K)+U3*CB(I-2))/(DS1+K)
      P1=P1+CA(I)
      PP1=PP1+(S+K)*CA(I)
      Q1=Q1+CB(I)
      QP1=QP1+(S+K-1)*CB(I)
!  ****  CHECK OVERFLOW LIMIT.
      TST=DMAX1(DABS(P1),DABS(Q1),DABS(PP1),DABS(QP1))
      IF(TST.GT.OVER) THEN
      NSUM=100
      RETURN
      ENDIF
      T1A=DABS(R1*PP1+H*(AK*P1+UQ*Q1))
      T1B=DABS(R1*QP1-H*(AK*Q1+UT*P1))
      T1=DMAX1(T1A,T1B)
      T2=DMAX1(DABS(CA(I)),DABS(CB(I)))
      TST=EPS*DMAX1(DABS(P1),DABS(Q1))
      IF(T1.LT.TST.AND.T2.LT.TST) GO TO 6
    5 CONTINUE
!  ****  RENORMALIZATION. (164)
    6 NSUM=K+1
      Q1=Q1/DABS(P1)
      P1=P1/DABS(P1)
      RETURN
!
!  **** MIDDLE REGION. (148-152)
!
    7 CONTINUE
      RHO=H/R0
      U0=(RV0+R0*(RVE+R0*(RV2+R0*RV3)))/SL
      U1=(RVE+R0*(2*RV2+R0*3*RV3))*H/SL
      U2=(RV2+R0*3*RV3)*H2/SL
      U3=RV3*H*H2/SL
      UB=U0-2*SL*R0
      UH=U1-2*SL*H
      UT=U0+U1+U2+U3
      UQ=UT-2*SL*R1
!
      CA(1)=P0
      CB(1)=Q0
      CA(2)=-RHO*(AK*CA(1)+UB*CB(1))
      CB(2)=RHO*(AK*CB(1)+U0*CA(1))
      CA(3)=-RHO*((AK+1)*CA(2)+UB*CB(2)+UH*CB(1))/2
      CB(3)=RHO*((AK-1)*CB(2)+U0*CA(2)+U1*CA(1))/2
      CA(4)=-RHO*((AK+2)*CA(3)+UB*CB(3)+UH*CB(2)+U2*CB(1))/3
      CB(4)=RHO*((AK-2)*CB(3)+U0*CA(3)+U1*CA(2)+U2*CA(1))/3
!
      P1=CA(1)+CA(2)+CA(3)+CA(4)
      PP1=CA(2)+2*CA(3)+3*CA(4)
      Q1=CB(1)+CB(2)+CB(3)+CB(4)
      QP1=CB(2)+2*CB(3)+3*CB(4)
!
      DO 9 I=5,60
      K=I-1
      CA(I)=-RHO*((AK+K-1)*CA(K)+UB*CB(K)+UH*CB(I-2)+U2*CB(I-3)+U3*CB(I-4))/K
      CB(I)=RHO*((AK-K+1)*CB(K)+U0*CA(K)+U1*CA(I-2)+U2*CA(I-3)+U3*CA(I-4))/K
      P1=P1+CA(I)
      PP1=PP1+K*CA(I)
      Q1=Q1+CB(I)
      QP1=QP1+K*CB(I)
!  ****  CHECK OVERFLOW LIMIT.
      TST=DMAX1(DABS(P1),DABS(Q1),DABS(PP1),DABS(QP1))
      IF(TST.GT.OVER) THEN
      NSUM=100
      RETURN
      ENDIF
      T1A=DABS(R1*PP1+H*(AK*P1+UQ*Q1))
      T1B=DABS(R1*QP1-H*(AK*Q1+UT*P1))
      T1=DMAX1(T1A,T1B)
      T2=DMAX1(DABS(CA(I)),DABS(CB(I)))
      TST=EPS*DMAX1(DABS(P1),DABS(Q1))
      IF(T1.LT.TST.AND.T2.LT.TST) GO TO 10
    9 CONTINUE
!
   10 NSUM=K+1
      RETURN
      END
!  **************************************************************
!                       SUBROUTINE SPLINE
!  **************************************************************
      SUBROUTINE SPLINE(X,Y,A,B,C,D,S1,SN,N)
!
!      CUBIC SPLINE INTERPOLATION BETWEEN TABULATED DATA.
!   INPUT:
!     X(I) (I=1, ...,N) ........ GRID POINTS.
!                    (THE X VALUES MUST BE IN INCREASING ORDER).
!     Y(I) (I=1, ...,N) ........ CORRESPONDING FUNCTION VALUES.
!     S1,SN ..... SECOND DERIVATIVES AT X(1) AND X(N).
!            (THE NATURAL SPLINE CORRESPONDS TO TAKING S1=SN=0).
!     N ........................ NUMBER OF GRID POINTS.
!      THE INTERPOLATING POLYNOMIAL IN THE I-TH INTERVAL, FROM
!   X(I) TO X(I+1), IS
!            PI(X) = A(I)+X*(B(I)+X*(C(I)+X*D(I)))
!   OUTPUT:
!     A(I),B(I),C(I),D(I) ...... SPLINE COEFFICIENTS.
!
!      REF.: M.J. MARON, 'NUMERICAL ANALYSIS: A PRACTICAL
!            APPROACH', MACMILLAN PUBL. CO., NEW YORK 1982.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      DIMENSION X(N),Y(N),A(N),B(N),C(N),D(N)
      IF(N.LT.4) THEN
      WRITE(6,10) N
   10 FORMAT(5X,'SPLINE INTERPOLATION CANNOT BE PERFORMED WITH',I4,' POINTS. STOP.')
      STOP
      ENDIF
      N1=N-1
      N2=N-2
!  ****  AUXILIARY ARRAYS H(=A) AND DELTA(=D).
      DO 1 I=1,N1
      IF(X(I+1)-X(I).LT.1.0D-13) THEN
      WRITE(6,11)
   11 FORMAT(5X,'SPLINE X VALUES NOT IN INCREASING ORDER. STOP.')
      STOP
      ENDIF
      A(I)=X(I+1)-X(I)
    1 D(I)=(Y(I+1)-Y(I))/A(I)
!  ****  SYMMETRIC COEFFICIENT MATRIX (AUGMENTED).
      DO 2 I=1,N2
      B(I)=2.0D0*(A(I)+A(I+1))
      K=N1-I+1
    2 D(K)=6.0D0*(D(K)-D(K-1))
      D(2)=D(2)-A(1)*S1
      D(N1)=D(N1)-A(N1)*SN
!  ****  GAUSS SOLUTION OF THE TRIDIAGONAL SYSTEM.
      DO 3 I=2,N2
      R=A(I)/B(I-1)
      B(I)=B(I)-R*A(I)
    3 D(I+1)=D(I+1)-R*D(I)
!  ****  THE SIGMA COEFFICIENTS ARE STORED IN ARRAY D.
      D(N1)=D(N1)/B(N2)
      DO 4 I=2,N2
      K=N1-I+1
    4 D(K)=(D(K)-A(K)*D(K+1))/B(K-1)
      D(N)=SN
!  ****  SPLINE COEFFICIENTS.
      SI1=S1
      DO 5 I=1,N1
      SI=SI1
      SI1=D(I+1)
      H=A(I)
      HI=1.0D0/H
      A(I)=(HI/6.0D0)*(SI*X(I+1)**3-SI1*X(I)**3)+HI*(Y(I)*X(I+1)-Y(I+1)*X(I))+(H/6.0D0)*(SI1*X(I)-SI*X(I+1))
      B(I)=(HI/2.0D0)*(SI1*X(I)**2-SI*X(I+1)**2)+HI*(Y(I+1)-Y(I))+(H/6.0D0)*(SI-SI1)
      C(I)=(HI/2.0D0)*(SI*X(I+1)-SI1*X(I))
    5 D(I)=(HI/6.0D0)*(SI1-SI)
      RETURN
      END
!  **************************************************************
!                       SUBROUTINE FINDI
!  **************************************************************
      SUBROUTINE FINDI(X,XC,N,I)
!
!      FINDS THE INTERVAL (X(I),X(I+1)) CONTAINING THE VALUE XC.
!   INPUT:
!     X(I) (I=1, ...,N) ........ GRID POINTS.
!                    (THE X VALUES MUST BE IN INCREASING ORDER).
!     XC ....................... POINT TO BE LOCATED.
!     N ........................ NUMBER OF GRID POINTS.
!   OUTPUT:
!     I ........................ INTERVAL INDEX.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      DIMENSION X(N)
      IF(XC.GT.X(N)) THEN
      I=N-1
      RETURN
      ENDIF
      IF(XC.LT.X(1)) THEN
      I=1
      RETURN
      ENDIF
      I=1
      I1=N
    1 IT=(I+I1)/2
      IF(XC.GT.X(IT)) I=IT
      IF(XC.LE.X(IT)) I1=IT
      IF(I1-I.GT.1) GO TO 1
      RETURN
      END
!  **************************************************************
!                       SUBROUTINE DCOUL
!  **************************************************************
      SUBROUTINE DCOUL(Z,E,K,R,FU,FL,GU,GL,ERR)
!
!     THIS SUBROUTINE COMPUTES RADIAL DIRAC-COULOMB WAVE FUNC-
!  TIONS FOR FREE STATES.
!
!  **** ALL QUANTITIES IN ATOMIC UNITS.
!
!  INPUT ARGUMENTS:
!     Z ........ FIELD STRENGTH, I.E. VALUE OF R*V(R) (ASSUMED
!                CONSTANT).
!     E ........ PARTICLE KINETIC ENERGY (POSITIVE).
!     K ........ ANGULAR MOMENTUM QUANTUM NUMBER KAPPA (.NE.0).
!     R ........ RADIAL DISTANCE (POSITIVE).
!
!  OUTPUT ARGUMENTS:
!     FU, FL ... UPPER AND LOWER COMPONENTS OF THE REGULAR RADIAL
!                COULOMB FUNCTION.
!     GU, GL ... UPPER AND LOWER COMPONENTS OF THE IRREGULAR RA-
!                DIAL COULOMB FUNCTION.
!     ERR ...... ACCURACY OF THE COMPUTED FUNCTIONS (RELATIVE UN-
!                CERTAINTY).
!  OUTPUT THROUGH COMMON/OCOUL/:
!     WAVNUM ... WAVE NUMBER.
!     ETA ...... SOMMERFELD'S PARAMETER.
!     DELTA .... COULOMB PHASE SHIFT (MODULUS 2*PI).
!
!     RADIAL FUNCTIONS ARE NORMALIZED SO THAT, FOR LARGE R, THE
!  UPPER COMPONENT OSCILLATES WITH UNIT AMPLITUDE.
!
!     OTHER SUBPROGRAMS REQUIRED: SUBROUTINES FCOUL AND SUM2F0,
!                                 AND FUNCTION CLGAM.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX (KIND=RegReal_K) (C),INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER (SL=137.03599976D0)  ! Speed of light (1/alpha)
      PARAMETER (SL2=SL*SL,TSL2=SL2+SL2,ALPHA=1.0D0/SL)
      PARAMETER (PI=3.1415926535897933D0,PIH=0.5D0*PI,TPI=PI+PI)
      COMMON/OCOUL/WAVNUM,ETA,DELTA
!
      IF(DABS(Z).GT.0.00001D0) THEN
        ZETA=Z*ALPHA
        ICAL=0
      ELSE
        ZETA=0.0D0
        ICAL=1
      ENDIF
      RLAMBS=K*K-ZETA*ZETA
      RLAMB=DSQRT(RLAMBS)
      PC=DSQRT(E*(E+TSL2))
      WAVNUM=PC/SL
      X=WAVNUM*R
!
      IF(E.LT.0.0001D0.OR.K.EQ.0) THEN
        FU=0.0D0
        FL=0.0D0
        GU=1.0D35
        GL=-1.0D35
        ERR=1.0D0
        DELTA=0.0D0
        IF(E.LT.0.0001D0) WRITE(6,2101)
 2101   FORMAT(1X,'*** ERROR IN DCOUL: E IS TOO SMALL.')
        IF(K.EQ.0) WRITE(6,2102)
 2102   FORMAT(1X,'*** ERROR IN DCOUL: K.EQ.0.')
        RETURN
      ENDIF
      IF(ICAL.EQ.1) GO TO 1
!
!  ****  PARAMETERS.
!
      RLAMB1=RLAMB-1.0D0
      W=E+SL2
      ETA=ZETA*W/PC
      RLA=DSQRT(RLAMBS+ETA*ETA)
      P1=K+RLAMB
      P2=RLAMB*SL2-K*W
      RNUR=ZETA*(W+SL2)
      RNUI=-P1*PC
      RNU=DATAN2(RNUI,RNUR)
      RNORM=1.0D0/(DSQRT(RNUR*RNUR+RNUI*RNUI)*RLAMB)
!
!  ****  COULOMB PHASE SHIFT.
!
      IF(K.GT.0) THEN
        L=K
      ELSE
        L=-K-1
      ENDIF
      DELTA0=DELTAC(ETA,RLAMB1)
      DELTA=RNU-(RLAMB-L-1)*PIH+DELTA0
      IF(Z.LT.0.0D0.AND.K.LT.0) THEN
        RNORM=-RNORM
        DELTA=DELTA-PI
      ENDIF
      IF(DELTA.GE.0.0D0) THEN
        DELTA=DMOD(DELTA,TPI)
      ELSE
        DELTA=-DMOD(-DELTA,TPI)
      ENDIF
!
!  ****  COULOMB FUNCTIONS.
!
      CALL FCOUL(ETA,RLAMB1,X,FM1,FPM1,GM1,GPM1,ERR)
      IF(ERR.GT.1.0D-6) THEN
        FU=0.0D0
        FL=0.0D0
        GU=1.0D35
        GL=-1.0D35
        ERR=1.0D0
        RETURN
      ENDIF
      SLA=(RLAMB/X)+(ETA/RLAMB)
      F=RLAMB*(SLA*FM1-FPM1)/RLA
      G=RLAMB*(SLA*GM1-GPM1)/RLA
!
!  ****  DIRAC-COULOMB WAVE FUNCTIONS.
!
      Q2=P1*P2*RNORM
      Q1=RLA*PC*RNORM
      P1=P1*Q1
      Q1=ZETA*Q1
      P2=ZETA*P2*RNORM
!
      FU=P1*F+P2*FM1
      FL=-Q1*F-Q2*FM1
      GU=P1*G+P2*GM1
      GL=-Q1*G-Q2*GM1
      RETURN
!
!  ****  Z=0. SPHERICAL BESSEL FUNCTIONS.
!
    1 CONTINUE
      RLAMB=IABS(K)
      CALL FCOUL(0.0D0,RLAMB,X,F,FP,G,GP,ERR)
      DELTA=0.0D0
      IF(ERR.GE.1.0D-6) THEN
        FU=0.0D0
        FL=0.0D0
        GU=1.0D35
        GL=-1.0D35
        ERR=1.0D0
        RETURN
      ENDIF
      FM1=(RLAMB*F/X)+FP
      GM1=(RLAMB*G/X)+GP
      FACT=DSQRT(E/(E+TSL2))
      IF(K.LT.0) THEN
        FU=FM1
        FL=-FACT*F
        GU=GM1
        GL=-FACT*G
      ELSE
        FU=F
        FL=FACT*FM1
        GU=G
        GL=FACT*GM1
      ENDIF
      RETURN
      END
!  **************************************************************
!                       SUBROUTINE FCOUL
!  **************************************************************
      SUBROUTINE FCOUL(ETA,RLAMB,X,F,FP,G,GP,ERR)
!
!     CALCULATION OF COULOMB FUNCTIONS FOR REAL ETA, RLAMB.GT.-1
!  AND X LARGER THAN, OR OF THE ORDER OF XTP0 (THE TURNING POINT
!  FOR RLAMB=0). STEED'S CONTINUED FRACTION METHOD IS COMBINED
!  WITH RECURSION RELATIONS AND AN ASYMPTOTIC EXPANSION. THE
!  OUTPUT VALUE ERR=1.0D0 INDICATES THAT THE ADOPTED EVALUATION
!  ALGORITHM IS NOT APPLICABLE (X IS TOO SMALL).
!
!  INPUT ARGUMENTS:
!     ETA ...... SOMMERFELD'S PARAMETER.
!     RLAMB .... ANGULAR MOMENTUM.
!     X ........ VARIABLE (=WAVE NUMBER TIMES RADIAL DISTANCE).
!
!  OUTPUT ARGUMENTS:
!     F, FP .... REGULAR FUNCTION AND ITS DERIVATIVE.
!     G, GP .... IRREGULAR FUNCTION AND ITS DERIVATIVE.
!     ERR ...... RELATIVE NUMERICAL UNCERTAINTY. A VALUE OF THE
!                ORDER OF 10**(-N) MEANS THAT THE CALCULATED
!                FUNCTIONS ARE ACCURATE TO N DECIMAL FIGURES.
!                THE MAXIMUM ACCURACY ATTAINABLE WITH DOUBLE
!                PRECISION ARITHMETIC IS ABOUT 1.0D-15.
!
!     OTHER SUBPROGRAMS REQUIRED: SUBROUTINE SUM2F0 AND
!                                 FUNCTIONS DELTAC AND CLGAM.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX (KIND=RegReal_K) (C),INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER (PI=3.1415926535897932D0,PIH=0.5D0*PI,TPI=PI+PI,EPS=1.0D-16,TOP=1.0D5,NTERM=1000)
!
      IF(RLAMB.LT.-0.999D0) THEN
        WRITE(6,'(1X,''*** ERROR IN FCOUL: RLAMB.LT.-0.999'')')
        STOP
      ENDIF
      IF(X.LT.EPS) GO TO 10
!
!  ****  NUMERICAL CONSTANTS.
!
      CI=DCMPLX(0.0D0,1.0D0)
      CI2=2.0D0*CI
      CIETA=CI*ETA
      X2=X*X
      ETA2=ETA*ETA
!
!  ****  TURNING POINT (XTP). (44)
!
      IF(RLAMB.GE.0.0D0) THEN
        XTP=ETA+DSQRT(ETA2+RLAMB*(RLAMB+1.0D0))
      ELSE
        XTP=EPS
      ENDIF
      ERRS=10.0D0
      IF(X.LT.XTP) GO TO 1
!
!  ************  ASYMPTOTIC EXPANSION. (71-75)
!
!  ****  COULOMB PHASE-SHIFT.
      DELTA=DELTAC(ETA,RLAMB)
!
      CPA=CIETA-RLAMB
      CPB=CIETA+RLAMB+1.0D0
      CPZ=CI2*X
      CALL SUM2F0(CPA,CPB,CPZ,C2F0,ERR1)
      CQA=CPA+1.0D0
      CQB=CPB+1.0D0
      CALL SUM2F0(CQA,CQB,CPZ,C2F0P,ERR2)
      C2F0P=CI*C2F0P*CPA*CPB/(2.0D0*X2)
!  ****  FUNCTIONS.
      THETA=X-ETA*DLOG(2.0D0*X)-RLAMB*PIH+DELTA
      IF(THETA.GT.1.0D4) THETA=DMOD(THETA,TPI)
      CEITH=CDEXP(CI*THETA)
      CGIF=C2F0*CEITH
      G=dble(CGIF)
      F=dble(-CI*CGIF)
!  ****  DERIVATIVES.
      CGIFP=(C2F0P+CI*(1.0D0-ETA/X)*C2F0)*CEITH
      GP=dble(CGIFP)
      FP=dble(-CI*CGIFP)
!  ****  GLOBAL UNCERTAINTY. THE WRONSKIAN MAY DIFFER FROM 1 DUE
!        TO TRUNCATION AND ROUNDOFF ERRORS.
      ERR=DMAX1(ERR1,ERR2,DABS(G*FP-F*GP-1.0D0))
      IF(ERR.LE.EPS) RETURN
      ERRS=ERR
!
!  ************  STEED'S CONTINUED FRACTION METHOD.
!
    1 CONTINUE
      CIETA2=CIETA+CIETA
      ETAX=ETA*X
!
!  ****  CONTINUED FRACTION FOR F. (60-70)
!
      INULL=0
      RLAMBN=RLAMB+1.0D0
      A1=-(RLAMBN+1.0D0)*(RLAMBN**2+ETA2)*X/RLAMBN
      B0=(RLAMBN/X)+(ETA/RLAMBN)
      B1=(2.0D0*RLAMBN+1.0D0)*(RLAMBN*(RLAMBN+1.0D0)+ETAX)
      FA3=B0
      FA2=B0*B1+A1
      FB3=1.0D0
      FB2=B1
      RF=FA3
!
      DO 2 N=2,NTERM
      RFO=RF
      DAF=DABS(RF)
      RLAMBN=RLAMB+N
      AN=-(RLAMBN**2-1.0D0)*(RLAMBN**2+ETA2)*X2
      BN=(2.0D0*RLAMBN+1.0D0)*(RLAMBN*(RLAMBN+1.0D0)+ETAX)
      FA1=FA2*BN+FA3*AN
      FB1=FB2*BN+FB3*AN
      TST=DABS(FB1)
!
      IF(TST.LT.1.0D-25) THEN
        IF(INULL.GT.0) STOP
        INULL=1
        FA3=FA2
        FA2=FA1
        FB3=FB2
        FB2=FB1
        RF=RFO
      ELSE
        FA3=FA2/TST
        FA2=FA1/TST
        FB3=FB2/TST
        FB2=FB1/TST
        RF=FA2/FB2
        IF(DABS(RF-RFO).LT.EPS*DAF) GO TO 3
      ENDIF
    2 CONTINUE
    3 CONTINUE
      IF(DAF.GT.1.0D-25) THEN
        ERRF=DABS(RF-RFO)/DAF
      ELSE
        ERRF=EPS
      ENDIF
      IF(ERRF.GT.ERRS) THEN
        ERR=ERRS
        IF(ERR.GT.1.0D-6) GO TO 10
        RETURN
      ENDIF
!
!  ****  DOWNWARD RECURSION FOR F AND FP. ONLY IF RLAMB.GT.1 AND
!        X.LT.XTP. (48,49)
!
      RLAMB0=RLAMB
      IF(X.GE.XTP.OR.RLAMB0.LT.1.0D0) THEN
        ISHIFT=0
        XTPC=XTP
        RFM=0.0D0
      ELSE
        FT=1.0D0
        FTP=RF
        IS0=nint(RLAMB0+1.0D-6)
        TST=X*(X-2.0D0*ETA)
        RL1T=0.0D0
        DO 4 I=1,IS0
        ETARL0=ETA/RLAMB0
        RL=DSQRT(1.0D0+ETARL0**2)
        SL=(RLAMB0/X)+ETARL0
        RLAMB0=RLAMB0-1.0D0
        FTO=FT
        FT=(SL*FT+FTP)/RL
        FTP=SL*FT-RL*FTO
        IF(FT.GT.1.0D10) THEN
          FTP=FTP/FT
          FT=1.0D0
        ENDIF
        RL1T=RLAMB0*(RLAMB0+1.0D0)
        IF(TST.GT.RL1T) THEN
          ISHIFT=I
          GO TO 5
        ENDIF
    4   CONTINUE
        ISHIFT=IS0
    5   CONTINUE
        XTPC=ETA+DSQRT(ETA2+RL1T)
        RFM=FTP/FT
      ENDIF
!
!  ****  CONTINUED FRACTION FOR P+CI*Q WITH RLAMB0. (76-79)
!
      INULL=0
      CAN=CIETA-ETA2-RLAMB0*(RLAMB0+1.0D0)
      CB0=X-ETA
      CBN=2.0D0*(X-ETA+CI)
      CFA3=CB0
      CFA2=CB0*CBN+CAN
      CFB3=1.0D0
      CFB2=CBN
      CPIQ=CFA3
!
      DO 6 N=2,NTERM
      CPIQO=CPIQ
      DAPIQ=CDABS(CPIQ)
      CAN=CAN+CIETA2+(N+N-2)
      CBN=CBN+CI2
      CFA1=CFA2*CBN+CFA3*CAN
      CFB1=CFB2*CBN+CFB3*CAN
      TST=CDABS(CFB1)
!
      IF(TST.LT.1.0D-25) THEN
        IF(INULL.GT.0) STOP
        INULL=1
        CFA3=CFA2
        CFA2=CFA1
        CFB3=CFB2
        CFB2=CFB1
        CPIQ=CPIQO
      ELSE
        CFA3=CFA2/TST
        CFA2=CFA1/TST
        CFB3=CFB2/TST
        CFB2=CFB1/TST
        CPIQ=CFA2/CFB2
        IF(CDABS(CPIQ-CPIQO).LT.EPS*DAPIQ) GO TO 7
      ENDIF
    6 CONTINUE
    7 CONTINUE
      IF(DAPIQ.GT.1.0D-25) THEN
        ERRPIQ=CDABS(CPIQ-CPIQO)/DAPIQ
      ELSE
        ERRPIQ=EPS
      ENDIF
      IF(ERRPIQ.GT.ERRS) THEN
        ERR=ERRS
        IF(ERR.GT.1.0D-6) GO TO 10
        RETURN
      ENDIF
      CPIQ=CI*CPIQ/X
!
      RP=dble(CPIQ)
      RQ=dble(-CI*CPIQ)
      IF(RQ.LE.1.0D-25) GO TO 10
      ERR=DMAX1(ERRF,ERRPIQ)
!
!  ****  INVERTING STEED'S TRANSFORMATION. (57,58)
!
      IF(ISHIFT.LT.1) THEN
        RFP=RF-RP
        F=DSQRT(RQ/(RFP**2+RQ**2))
        IF(FB2.LT.0.0D0) F=-F
        FP=RF*F
        G=RFP*F/RQ
        GP=(RP*RFP-RQ**2)*F/RQ
        IF(X.LT.XTP.AND.G.GT.TOP*F) GO TO 10
      ELSE
        RFP=RFM-RP
        FM=DSQRT(RQ/(RFP**2+RQ**2))
        G=RFP*FM/RQ
        GP=(RP*RFP-RQ**2)*FM/RQ
        IF(X.LT.XTPC.AND.G.GT.TOP*FM) GO TO 10
!  ****  UPWARD RECURSION FOR G AND GP (IF ISHIFT.GT.0). (50,51)
        DO 8 I=1,ISHIFT
        RLAMB0=RLAMB0+1.0D0
        ETARL0=ETA/RLAMB0
        RL=DSQRT(1.0D0+ETARL0**2)
        SL=(RLAMB0/X)+ETARL0
        GO=G
        G=(SL*GO-GP)/RL
        GP=RL*GO-SL*G
        IF(G.GT.1.0D35) GO TO 10
    8   CONTINUE
    9   W=RF*G-GP
        F=1.0D0/W
        FP=RF/W
      ENDIF
!  ****  THE WRONSKIAN MAY DIFFER FROM 1 DUE TO ROUNDOFF ERRORS.
      ERR=DMAX1(ERR,DABS(FP*G-F*GP-1.0D0))
      IF(ERR.LT.1.0D-6) RETURN
!
   10 F=0.0D0
      FP=0.0D0
      G=1.0D35
      GP=-1.0D35
      ERR=1.0D0
      RETURN
      END
!  **************************************************************
!                       SUBROUTINE SUM2F0
!  **************************************************************
      SUBROUTINE SUM2F0(CA,CB,CZ,CF,ERR)
!
!     SUMMATION OF THE 2F0(CA,CB;CS) HYPERGEOMETRIC ASYMPTOTIC
!  SERIES. THE POSITIVE AND NEGATIVE CONTRIBUTIONS TO THE REAL
!  AND IMAGINARY PARTS ARE ADDED SEPARATELY TO OBTAIN AN ESTIMATE
!  OF ROUNDING ERRORS.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX (KIND=RegReal_K) (C),INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER (EPS=1.0D-16,ACCUR=0.5D-15,NTERM=75)
      RRP=1.0D0
      RRN=0.0D0
      RIP=0.0D0
      RIN=0.0D0
      CDF=1.0D0
      ERR2=0.0D0
      ERR3=1.0D0
      AR=0.0D0
      AF=0.0D0
      DO 1 I=1,NTERM
      J=I-1
      CDF=CDF*(CA+J)*(CB+J)/(I*CZ)
      ERR1=ERR2
      ERR2=ERR3
      ERR3=CDABS(CDF)
      IF(ERR1.GT.ERR2.AND.ERR2.LT.ERR3) GO TO 2
      AR=dble(CDF)
      IF(AR.GT.0.0D0) THEN
        RRP=RRP+AR
      ELSE
        RRN=RRN+AR
      ENDIF
      AI=dble(DCMPLX(0.0D0,-1.0D0)*CDF)
      IF(AI.GT.0.0D0) THEN
        RIP=RIP+AI
      ELSE
        RIN=RIN+AI
      ENDIF
      CF=DCMPLX(RRP+RRN,RIP+RIN)
      AF=CDABS(CF)
      IF(AF.GT.1.0D25) THEN
        CF=0.0D0
        ERR=1.0D0
        RETURN
      ENDIF
      IF(ERR3.LT.1.0D-25*AF.OR.ERR3.LT.EPS) THEN
         ERR=EPS
         RETURN
      ENDIF
    1 CONTINUE
!  ****  ROUNDOFF ERROR.
    2 CONTINUE
      TR=DABS(RRP+RRN)
      IF(TR.GT.1.0D-25) THEN
        ERRR=(RRP-RRN)*ACCUR/TR
      ELSE
        ERRR=1.0D0
      ENDIF
      TI=DABS(RIP+RIN)
      IF(TI.GT.1.0D-25) THEN
        ERRI=(RIP-RIN)*ACCUR/TI
      ELSE
        ERRI=1.0D0
      ENDIF
!  ****  ... AND TRUNCATION ERROR.
      IF(AR.GT.1.0D-25) THEN
      ERR=DMAX1(ERRR,ERRI)+ERR2/AF
      ELSE
      ERR=DMAX1(ERRR,ERRI)
      ENDIF
      RETURN
      END
!  **************************************************************
!                         FUNCTION DELTAC
!  **************************************************************
      FUNCTION DELTAC(ETA,RLAMB)
!
!     CALCULATION OF COULOMB PHASE SHIFT (MODULUS 2*PI). (47)
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX (KIND=RegReal_K) (C), INTEGER (KIND=RegInt_K) (I-N)
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI)
      CI=DCMPLX(0.0D0,1.0D0)
!  ****  COULOMB PHASE-SHIFT.
      DELTAC=dble(-CI*CLGAM(RLAMB+1.0D0+CI*ETA))
      IF(DELTAC.GE.0.0D0) THEN
        DELTAC=DMOD(DELTAC,TPI)
      ELSE
        DELTAC=-DMOD(-DELTAC,TPI)
      ENDIF
      RETURN
      END
!  **************************************************************
!                       FUNCTION CLGAM
!  **************************************************************
      FUNCTION CLGAM(CZ)
!
!     THIS FUNCTION GIVES LOG(GAMMA(CZ)) FOR COMPLEX ARGUMENTS.
!
!   REF.: M. ABRAMOWITZ AND I.A. STEGUN, 'HANDBOOK OF MATHEMATI-
!         CAL FUNCTIONS'. DOVER, NEW YORK (1974). PP 255-257.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX (KIND=RegReal_K) (C), INTEGER (KIND=RegInt_K) (I-N)
      CZA=CZ
      ICONJ=0
      AR=dble(CZA)
      CLGAM=36.84136149D0
      IF(CDABS(CZA).LT.1.0D-16) RETURN
!
      AI=dble(CZA*DCMPLX(0.0D0,-1.0D0))
      IF(AI.GT.0.0D0) THEN
        ICONJ=0
      ELSE
        ICONJ=1
        CZA=DCONJG(CZA)
      ENDIF
!
      CZFAC=1.0D0
      CZFL=0.0D0
    1 CZFAC=CZFAC/CZA
      IF(CDABS(CZFAC).GT.1.0D8) THEN
        CZFL=CZFL+CDLOG(CZFAC)
        CZFAC=1.0D0
      ENDIF
      CZA=CZA+1.0D0
      AR=dble(CZA)
      IF(CDABS(CZA).LT.1.0D-16) RETURN
      IF(CDABS(CZA).GT.15.0D0.AND.AR.GT.0.0D0) GO TO 2
      GO TO 1
!  ****  STIRLING'S EXPANSION OF CDLOG(GAMMA(CZA)).
    2 CZI2=1.0D0/(CZA*CZA)
      CZS=(43867.0D0/244188.0D0)*CZI2
      CZS=(CZS-3617.0D0/122400.0D0)*CZI2
      CZS=(CZS+1.0D0/156.0D0)*CZI2
      CZS=(CZS-691.0D0/360360.0D0)*CZI2
      CZS=(CZS+1.0D0/1188.0D0)*CZI2
      CZS=(CZS-1.0D0/1680.0D0)*CZI2
      CZS=(CZS+1.0D0/1260.0D0)*CZI2
      CZS=(CZS-1.0D0/360.0D0)*CZI2
      CZS=(CZS+1.0D0/12.0D0)/CZA
      CLGAM=(CZA-0.5D0)*CDLOG(CZA)-CZA+9.1893853320467274D-1+CZS+CZFL+CDLOG(CZFAC)
      IF(ICONJ.EQ.1) CLGAM=DCONJG(CLGAM)
      RETURN
      END
!  **************************************************************
!                         FUNCION BESJN_MOD
!  **************************************************************
      FUNCTION BESJN_MOD(JY,N,X)
!
!      THIS FUNCTION COMPUTES THE SPHERICAL BESSEL FUNCTIONS OF
!   THE FIRST KIND AND SPHERICAL BESSEL FUNCTIONS OF THE SECOND
!   KIND (ALSO KNOWN AS SPHERICAL NEUMANN FUNCTIONS) FOR REAL
!   POSITIVE ARGUMENTS.
!
!      INPUT:
!         JY ...... KIND: 1(BESSEL) OR 2(NEUMANN).
!         N ....... ORDER (INTEGER).
!         X ....... ARGUMENT (REAL AND POSITIVE).
!
!   REF.: M. ABRAMOWITZ AND I.A. STEGUN, 'HANDBOOK OF MATHEMATI-
!         CAL FUNCTIONS'. DOVER, NEW YORK (1974). PP 435-478.
!
      use constants, only: RegInt_K, RegReal_K
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (KIND=RegInt_K) (I-N)
      double precision F1
      IF(X.LT.0) THEN
        WRITE(6,1000)
 1000   FORMAT(1X,'*** NEGATIVE ARGUMENT IN FUNCTION BESJN.')
        STOP
      ENDIF
!  ****  ORDER AND PHASE CORRECTION FOR NEUMANN FUNCTIONS.
!        ABRAMOWITZ AND STEGUN, EQ. 10.1.15.
      IF(JY.EQ.2) THEN
        NL=-N-1
        IPH=2*MOD(IABS(N),2)-1
      ELSE
        NL=N
        IPH=1
      ENDIF
!  ****  SELECTION OF CALCULATION MODE.
      IF(NL.LT.0) GO TO 10
      IF(X.GT.1.0D0*NL) GO TO 7
      XI=X*X
      IF(XI.GT.NL+NL+3.0D0) GO TO 4
!  ****  POWER SERIES FOR SMALL ARGUMENTS AND POSITIVE ORDERS.
!        ABRAMOWITZ AND STEGUN, EQ. 10.1.2.
      F1=1.0D0
      IP=1
      IF(NL.NE.0) THEN
        DO 1 I=1,NL
        IP=IP+2
    1   F1=F1*X/IP
      ENDIF
      XI=0.5D0*XI
      BESJN_MOD=1.0D0
      PS=1.0D0
      DO 2 I=1,500
      IP=IP+2
      PS=-PS*XI/(I*IP)
      BESJN_MOD=BESJN_MOD+PS
      IF(DABS(PS).LT.1.0D-18*DABS(BESJN_MOD)) GO TO 3
    2 CONTINUE
    3 BESJN_MOD=IPH*F1*BESJN_MOD
      RETURN
!  ****  MILLER'S METHOD FOR POSITIVE ORDERS AND INTERMEDIATE
!        ARGUMENTS. ABRAMOWITZ AND STEGUN, EQ. 10.1.19.
    4 XI=1.0D0/X
      F2=0.0D0
      F3=1.0D-35
      IP=2*(NL+31)+3
      DO 5 I=1,31
      F1=F2
      F2=F3
      IP=IP-2
      F3=IP*XI*F2-F1
      IF(DABS(F3).GT.1.0D30) THEN
        F2=F2/F3
        F3=1.0D0
      ENDIF
    5 CONTINUE
      BESJN_MOD=1.0D0
      F2=F2/F3
      F3=1.0D0
      DO 6 I=1,NL
      F1=F2
      F2=F3
      IP=IP-2
      F3=IP*XI*F2-F1
      IF(DABS(F3).GT.1.0D30) THEN
        BESJN_MOD=BESJN_MOD/F3
        F2=F2/F3
        F3=1.0D0
      ENDIF
    6 CONTINUE
      BESJN_MOD=IPH*XI*DSIN(X)*BESJN_MOD/F3
      RETURN
!  ****  RECURRENCE RELATION FOR ARGUMENTS GREATER THAN ORDER.
!        ABRAMOWITZ AND STEGUN, EQ. 10.1.19.
    7 XI=1.0D0/X
      F3=XI*DSIN(X)
      IF(NL.EQ.0) GO TO 9
      F2=F3
      F3=XI*(F2-DCOS(X))
      IF(NL.EQ.1) GO TO 9
      IP=1
      DO 8 I=2,NL
      F1=F2
      F2=F3
      IP=IP+2
    8 F3=IP*XI*F2-F1
    9 BESJN_MOD=IPH*F3
      RETURN
!  ****  RECURRENCE RELATION FOR NEGATIVE ORDERS.
!        ABRAMOWITZ AND STEGUN, EQ. 10.1.19.
   10 NL=IABS(NL)
      IF(X.LT.7.36D-1*(NL+1)*1.0D-35**(1.0D0/(NL+1))) THEN
        BESJN_MOD=-1.0D35
        RETURN
      ENDIF
      XI=1.0D0/X
      F3=XI*DSIN(X)
      F2=XI*(F3-DCOS(X))
      IP=3
      DO 11 I=1,NL
      F1=F2
      F2=F3
      IP=IP-2
      F3=IP*XI*F2-F1
      IF(DABS(F3).GT.1.0D35) THEN
        BESJN_MOD=-1.0D35
        RETURN
      ENDIF
   11 CONTINUE
      BESJN_MOD=IPH*F3
      RETURN
      END



      SUBROUTINE  akima_interpol(NP,ND,XD,YD,NI,XI,YI)
! Univariate Interpolation (Improved Akima Method)
!
! Hiroshi Akima
! U.S. Department of Commerce, NTIA/ITS
! Version of 89/07/04
!
! This subroutine performs univariate interpolation.  It is based
! on the improved A method developed by Hiroshi Akima, 'A method
! of univariate interpolation that has the accuracy of a third-
! degree polynomial,' ACM TOMS, vol. xx, pp. xxx-xxx, 19xx.  (The
! equation numbers referred to in the comments below are those in
! the paper.)
!
! In this method, the interpolating function is a piecewise
! function composed of a set of polynomials applicable to
! successive intervals of the given data points.  This method
! uses third-degree polynomials as the default, but the user has
! an option to use higher-degree polynomial to reduce undulations
! in resulting curves.
!
! This method has the accuracy of a third-degree polynomial if
! the degree of the polynomials for the interpolating function is
! set to three.
!
! The input arguments are
!   NP = degree of the polynomials for the interpolating
!        function,
!   ND = number of input data points
!        (must be equal to 2 or greater),
!   XD = array of dimension ND, containing the abscissas of
!        the input data points
!        (must be in a monotonic increasing order),
!   YD = array of dimension ND, containing the ordinates of
!        the input data points,
!   NI = number of points for which interpolation is desired
!        (must be equal to 1 or greater),
!   XI = array of dimension NI, containing the abscissas of

!        the desired points.
!
! The output argument is
!   YI = array of dimension NI, where the ordinates of the
!        desired points are to be stored.
!
! If an integer value smaller than 3 is given to the NP argument,
! this subroutine assumes NP = 3.
!
! The XI array elements need not be monotonic, but this
! subroutine interpolates faster if the XI array elements are
! given in a monotonic order.
!
! If the XI array element is less than XD(1) or greater than
! XD(ND), this subroutine linearly extrapolates the YI value.
!
!
      INTEGER            NP, ND, NI
      DOUBLE PRECISION   XD(ND), YD(ND), XI(NI), YI(NI)
      INTEGER IINTPV, ID1, ID2, ID3
! Initialization added by ThPe
      DOUBLE PRECISION X0, Y0, A1, YP0, YP1, DX, DY, DY1, DY2, DY3
      DOUBLE PRECISION A0, A2, A3, AA0, AA1, A12, A13
      DOUBLE PRECISION SMPEF, SMWTF, SMPEI, SMWTI
      DOUBLE PRECISION DNM, B0, B1, WT, EPSLN, U, UC, V, VOL
      DOUBLE PRECISION XX, XII, X1, X2, X3, Y1, Y2, Y3, YP, T0, T1, DLT, PE, SX, SY, SXX, SXY, DY0
      X0=0.d0
      Y0=0.d0
      YP0=0.d0
      YP1=0.d0
      DX=0.d0
      DY=0.d0
      A0=0.d0
      A1=0.d0
      A2=0.d0
      A3=0.d0
      AA0=0.d0
      AA1=0.d0
      IINTPV=0
      ID1=0
      ID2=0
      ID3=0
! Error check
   10 IF (ND.LE.1)   GO TO 90
      IF (NI.LE.0)   GO TO 91
      DO 11  ID=2,ND
        IF (XD(ID).LE.XD(ID-1))     GO TO 92
   11 CONTINUE
! Branches off special cases.
      IF (ND.LE.4)   GO TO 50
! General case  --  Five data points of more
! Calculates some local variables.
   20 NP0=MAX(3,NP)
      NPM1=NP0-1
      RENPM1=NPM1
      RENNM2=NP0*(NP0-2)
! Main calculation for the general case
! First (outermost) DO-loop with respect to the desired points
   30 DO 39  II=1,NI
        IF (II.EQ.1)      IINTPV=-1
        XII=XI(II)
! Locates the interval that includes the desired point by binary
! search.
        IF (XII.LE.XD(1))  THEN
          IINT=0
        ELSE IF (XII.LT.XD(ND))  THEN
          IDMN=1
          IDMX=ND
          IDMD=(IDMN+IDMX)/2
   31     IF (XII.GE.XD(IDMD))  THEN
            IDMN=IDMD
          ELSE
            IDMX=IDMD
          END IF
          IDMD=(IDMN+IDMX)/2
          IF (IDMD.GT.IDMN)    GO TO 31
          IINT=IDMD
        ELSE
          IINT=ND
        END IF
! End of locating the interval of interest
! Interpolation or extrapolation in one of the three subcases
        IF (IINT.LE.0)  THEN
! Subcase 1  --  Linear extrapolation when the abscissa of the
!                desired point is equal to that of the first data
!                point or less.
! Estimates the first derivative when the interval is not the
! same as the one for the previous desired point.  --
! cf. Equation (8)
          IF (IINT.NE.IINTPV)  THEN
            IINTPV=IINT
            X0=XD(1)
            X1=XD(2)-X0
            X2=XD(3)-X0
            X3=XD(4)-X0
            Y0=YD(1)
            Y1=YD(2)-Y0
            Y2=YD(3)-Y0
            Y3=YD(4)-Y0
            DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
            A1=(((X2*X3)**2)*(X3-X2)*Y1+((X3*X1)**2)*(X1-X3)*Y2+((X1*X2)**2)*(X2-X1)*Y3)/DLT
          END IF
! Evaluates the YI value.
          YI(II)=Y0+A1*(XII-X0)
! End of Subcase 1
        ELSE IF (IINT.GE.ND)  THEN
! Subcase 2  --  Linear extrapolation when the abscissa of the
!                desired point is equal to that of the last data
!                point or greater.
! Estimates the first derivative when the interval is not the
! same as the one for the previous desired point.  --
! cf. Equation (8)
          IF (IINT.NE.IINTPV)  THEN
            IINTPV=IINT
            X0=XD(ND)
            X1=XD(ND-1)-X0
            X2=XD(ND-2)-X0
            X3=XD(ND-3)-X0
            Y0=YD(ND)
            Y1=YD(ND-1)-Y0
            Y2=YD(ND-2)-Y0
            Y3=YD(ND-3)-Y0
            DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
            A1=(((X2*X3)**2)*(X3-X2)*Y1+((X3*X1)**2)*(X1-X3)*Y2+((X1*X2)**2)*(X2-X1)*Y3)/DLT
          END IF
! Evaluates the YI value.
          YI(II)=Y0+A1*(XII-X0)
! End of Subcase 2
        ELSE
! Subcase 3  --  Interpolation when the abscissa of the desired
!                point is  between those of the first and last
!                data points.
! Calculates the coefficients of the third-degree polynomial (for
! NP.LE.3) or the factors for the higher-degree polynomials (for
! NP.GT.3), when the interval is not the same as the one for the
! previous desired point.
          IF (IINT.NE.IINTPV)  THEN
            IINTPV=IINT
! The second DO-loop with respect to the two endpoints of the
! interval
            DO 37  IEPT=1,2
! Calculates the estimate of the first derivative at an endpoint.
! Initial setting for calculation
              ID0=IINT+IEPT-1
              X0=XD(ID0)
              Y0=YD(ID0)
              SMPEF=0.d0
              SMWTF=0.d0
              SMPEI=0.d0
              SMWTI=0.d0
! The third (innermost) DO-loop with respect to the four primary
! estimate of the first derivative
              DO 36  IPE=1,4
! Selects point numbers of four consecutive data points for
! calculating the primary estimate of the first derivative.
                IF (IPE.EQ.1)  THEN
                  ID1=ID0-3
                  ID2=ID0-2
                  ID3=ID0-1
                ELSE IF (IPE.EQ.2)  THEN
                  ID1=ID0+1
                ELSE IF (IPE.EQ.3)  THEN
                  ID2=ID0+2
                ELSE
                  ID3=ID0+3
                END IF
! Checks if any point number falls outside the legitimate range
! (between 1 and ND).  Skips calculation of the primary estimate
! if any does.
                IF (ID1.LT.1.OR.ID2.LT.1.OR.ID3.LT.1.OR.ID1.GT.ND.OR.ID2.GT.ND.OR.ID3.GT.ND) GO TO 36
! Calculates the primary estimate of the first derivative  --
! cf. Equation (8)
                X1=XD(ID1)-X0
                X2=XD(ID2)-X0
                X3=XD(ID3)-X0
                Y1=YD(ID1)-Y0
                Y2=YD(ID2)-Y0
                Y3=YD(ID3)-Y0
                DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
                PE=(((X2*X3)**2)*(X3-X2)*Y1+((X3*X1)**2)*(X1-X3)*Y2+((X1*X2)**2)*(X2-X1)*Y3)/DLT
! Calculates the volatility factor, VOL, and distance factor,
! SXX, for the primary estimate.  --  cf. Equations (9) and (11)
                SX=X1+X2+X3
                SY=Y1+Y2+Y3
                SXX=X1*X1+X2*X2+X3*X3
                SXY=X1*Y1+X2*Y2+X3*Y3
                DNM=4.d0*SXX-SX*SX
                B0=(SXX*SY-SX*SXY)/DNM
                B1=(4.d0*SXY-SX*SY)/DNM
                DY0=-B0
                DY1=Y1-(B0+B1*X1)
                DY2=Y2-(B0+B1*X2)
                DY3=Y3-(B0+B1*X3)
                VOL=DY0*DY0+DY1*DY1+DY2*DY2+DY3*DY3
! Calculates the EPSLN value, which is used to decide whether or
! not the volatility factor, VOL, is essentially zero.
                EPSLN=(YD(ID0)**2+YD(ID1)**2+YD(ID2)**2+YD(ID3)**2)*1.0d-12
! Accumulates the weighted primary estimates.  --
! cf. Equations (13) and (14)
                IF (VOL.GT.EPSLN)  THEN
! - For finite weight.
                  WT=1.d0/(VOL*SXX)
                  SMPEF=SMPEF+PE*WT
                  SMWTF=SMWTF+WT
                ELSE
! - For infinite weight.
                  SMPEI=SMPEI+PE
                  SMWTI=SMWTI+1.0
                END IF
   36         CONTINUE
! End of the third DO-loop
! Calculates the final estimate of the first derivative.  --
! cf. Equation (14)
              IF (SMWTI.LT.0.5)  THEN
! - When no infinite weights exist.
                YP=SMPEF/SMWTF
              ELSE
! - When infinite weights exist.
                YP=SMPEI/SMWTI
              END IF
              IF (IEPT.EQ.1)  THEN
                YP0=YP
              ELSE
                YP1=YP
              END IF
! End of the calculation of the estimate of the first derivative
! at an endpoint
   37       CONTINUE
! End of the second DO-loop
            IF (NP0.LE.3)  THEN
! Calculates the coefficients of the third-degree polynomial
! (when NP.LE.3).  --  cf. Equation (4)
              DX=XD(IINT+1)-XD(IINT)
              DY=YD(IINT+1)-YD(IINT)
              A0=YD(IINT)
              A1=YP0
              YP1=YP1-YP0
              YP0=YP0-DY/DX
              A2=-(3.d0*YP0+YP1)/DX
              A3= (2.d0*YP0+YP1)/(DX*DX)
            ELSE
! Calculates the factors for the higher-degree polynomials
! (when NP.GT.3).  --  cf. Equation (20)
              DX=XD(IINT+1)-XD(IINT)
              DY=YD(IINT+1)-YD(IINT)
              T0=YP0*DX-DY
              T1=YP1*DX-DY
              AA0= (T0+RENPM1*T1)/RENNM2
              AA1=-(RENPM1*T0+T1)/RENNM2
            END IF
          END IF
! End of the calculation of the coefficients of the third-degree
! polynomial (when NP.LE.3) or the factors for the higher-degree
! polynomials (when NP.GT.3), when the interval is not the same
! as the one for the previous desired point.
! Evaluates the YI value.
          IF (NP0.LE.3)  THEN
! - With a third-degree polynomial.  --  cf. Equation (3)
            XX=XII-XD(IINT)
            YI(II)=A0+XX*(A1+XX*(A2+XX*A3))
          ELSE
! - With a higher-degree polynomial.  --  cf. Equation (19)
            U=(XII-XD(IINT))/DX
            UC=1.0-U
            V=AA0*((U**NP0)-U)+AA1*((UC**NP0)-UC)
            YI(II)=YD(IINT)+DY*U+V
          END IF
! End of Subcase 3
        END IF
   39 CONTINUE
! End of the first DO-loop
! End of general case
      RETURN
! Special cases  --  Four data points or less
! Preliminary processing for the special cases
   50 X0=XD(1)
      Y0=YD(1)
      X1=XD(2)-X0
      Y1=YD(2)-Y0
      IF (ND.EQ.2)   GO TO 60
      X2=XD(3)-X0
      Y2=YD(3)-Y0
      IF (ND.EQ.3)   GO TO 70
      X3=XD(4)-X0
      Y3=YD(4)-Y0
      GO TO 80
! Special Case 1  --  Two data points
! (Linear interpolation and extrapolation)
   60 A1=Y1/X1
      DO 61  II=1,NI
        YI(II)=Y0+A1*(XI(II)-X0)
   61 CONTINUE
! End of Special Case 1
      RETURN
! Special Case 2  --  Three data points
! (Quadratic interpolation and linear extrapolation)
   70 DLT=X1*X2*(X2-X1)
      A1=(X2*X2*Y1-X1*X1*Y2)/DLT
      A2=(X1*Y2-X2*Y1)/DLT
      A12=2.d0*A2*X2+A1
      DO 71  II=1,NI
        XX=XI(II)-X0
        IF (XX.LE.0.0)  THEN
          YI(II)=Y0+A1*XX
        ELSE IF (XX.LT.X2) THEN
          YI(II)=Y0+XX*(A1+XX*A2)
        ELSE
          YI(II)=Y0+Y2+A12*(XX-X2)
        END IF
   71 CONTINUE
! End of Special Case 2
      RETURN
! Special Case 3  --  Four data points
! (Cubic interpolation and linear extrapolation)
   80 DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
      A1=(((X2*X3)**2)*(X3-X2)*Y1+((X3*X1)**2)*(X1-X3)*Y2+((X1*X2)**2)*(X2-X1)*Y3)/DLT
      A2=(X2*X3*(X2*X2-X3*X3)*Y1+X3*X1*(X3*X3-X1*X1)*Y2+X1*X2*(X1*X1-X2*X2)*Y3)/DLT
      A3=(X2*X3*(X3-X2)*Y1+X3*X1*(X1-X3)*Y2+X1*X2*(X2-X1)*Y3)/DLT
      A13=(3.d00*A3*X3+2.0*A2)*X3+A1
      DO 81  II=1,NI
        XX=XI(II)-X0
        IF (XX.LE.0.0)  THEN
          YI(II)=Y0+A1*XX
        ELSE IF (XX.LT.X3) THEN
          YI(II)=Y0+XX*(A1+XX*(A2+XX*A3))
        ELSE
          YI(II)=Y0+Y3+A13*(XX-X3)
        END IF
   81 CONTINUE
! End of Special Case 3
      RETURN
! Error exit
   90 WRITE (*,99090) ND
      GO TO 99
   91 WRITE (*,99091) NI
      GO TO 99
   92 WRITE (*,99092) ID,XD(ID-1),XD(ID),XI(1)
   99 WRITE (*,99099)
      write(*,*) XD
      write(*,*) YD
      stop
      RETURN
! Format statements for error messages
99090 FORMAT (1X/ ' ***   Insufficient data points.'7X,'ND =',I3)
99091 FORMAT (1X/ ' ***   No desired points.'7X,'NI =',I3)
99092 FORMAT (1X/ ' ***   Two data points identical or out of sequence.'/7X,'ID, XD(ID-1), XD(ID), XI(1) =',I5,3F15.10)
99099 FORMAT (' Error detected in the UVIP3P subroutine'/)
      END SUBROUTINE akima_interpol


