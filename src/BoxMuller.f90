Module BoxMuller  

Contains 

    subroutine GaussRand(gauss_rand) 
    ! Box-Muller transform
    ! return a random number from gaussian distribution
      Implicit None
      !OUTPUT
      real(4)   gauss_rand
      !WORK
      real(4)   PI
      parameter( PI=3.1415926535897934)

      save     compute,store ! define this variables as static

      real(4)   phi, r  ! temporary parameters

      real(4)   store
      logical  compute
      data     compute  / .true. /
      !BODY
      if (compute) then
        r = dble(randm())
        ! to avoid NaN - 
        ! r must be large enough
        do while (r < 0.01)
            r = dble(randm())
        enddo
        phi = dble(randm())

        ! Box-Muller transform (first variant)
        gauss_rand = cos(2.0*PI*phi)*sqrt(-2.0*log(r))
        store = sin(2.0*PI*phi)*sqrt(-2.0*log(r))
        compute = .false.
      else
         gauss_rand  = store
         compute  = .true.
      end if

      return
      end subroutine GaussRand
!**********************************************************************


!**********************************************************************
    real*4 function randm()
!======================================================================
! randmon number generator (0,1]
!======================================================================
      save
      real*4 temp(2048)
      data   key /2049/, len /2048/

      key = key + 1
      if ( key .gt. len ) then
         if ( key .eq. 2050 ) then
            call ctime00 ( tx, ty )
            ij=1 + int( 31327. * ( ty-aint( ty ) ) )
            do i = 1, 10000
               call ctime00 ( tx, ty )
            end do
            kl = 1 + int(30080.*(10000.*(tx+ty)-aint(10000.*(tx+ty))))
            if ( ij .lt. 0  .or.  ij .gt. 31328 ) ij = 1
            if ( kl .lt. 0  .or.  kl .gt. 30081 ) kl = 1
            call rmarin ( ij, kl )
         end if
         call ranmar ( temp, len )
         key = 1
      end if
      randm = temp ( key )
      return
      end function randm
!=============================================================================

!=============================================================================
!>	@brief Creating initial array for the random number generation
!!	@param	ij	[IN]	от 0 до 31328 
!!	@param	kl	[IN}	от 0 до 30081
!!	@author P.G. Khalatur
    subroutine rmarin (ij, kl)
!=============================================================================
! this is the initialization routine for the randmom number generator ranmar()
! note: the seed variables can have values between:    0 <= ij <= 31328
!                                                      0 <= kl <= 30081
! the randmom number sequences created by these two seeds are of sufficient
! length to complete an entire calculation with. for example, if sveral
! different groups are working on different parts of the same calculation,
! each group could be assigned its own ij seed. this would leave each group
! with 30000 choices for the second seed. that is to say, this randmom
! number generator can create 900 million different subsequences -- with
! each subsequence having a length of approximately 10^30.
!
! use ij = 1802 & kl = 9373 to test the randmom number generator. the
! subroutine ranmar should be used to generate 20000 randmom numbers.
! then display the next six randmom numbers generated multiplied by 4096*4096
! if the randmom number generator is working properly, the randmom numbers
! should be:
!           6533892.0  14220222.0  7275067.0
!           6172232.0  8354498.0   10633180.0
!==============================================================================
      real u(97), c, cd, cm
      integer i97, j97
      logical test
      common /raset1/ u, c, cd, cm, i97, j97, test

      test=.false.


      i = mod(ij/177, 177) + 2
      j = mod(ij    , 177) + 2
      k = mod(kl/169, 178) + 1
      l = mod(kl,     169)

      do ii = 1, 97
         s = 0.0
         t = 0.5
         do jj = 1, 24
            m = mod(mod(i*j, 179)*k, 179)
            i = j
            j = k
            k = m
            l = mod(53*l+1, 169)
            if (mod(l*m, 64) .ge. 32) then
               s = s + t
            endif
            t = 0.5 * t
ENDDO
         u(ii) = s
ENDDO

      c = 362436.0 / 16777216.0
      cd = 7654321.0 / 16777216.0
      cm = 16777213.0 /16777216.0

      i97 = 97
      j97 = 33

      test = .true.
      return
      end
!==========================================================================

!==========================================================================
!>	@brief Randim numbers generator developed by G. Marsaglia
!!	@param	rvec	[OUT]	array of the random numbers
!!	@param	len 	[IN]	size of the output array
!!	@author P.G. Khalatur
      subroutine ranmar ( rvec, len )
!==========================================================================
! this is the randmom number generator proposed by george marsaglia in
! florida state university report: fsu-scri-87-50
! it was slightly modified by f. james to produce an array of pseudorandmom
! numbers.
!==========================================================================
      real rvec(*)
      real u(97), c, cd, cm
      integer i97, j97
      logical test
      common /raset1/ u, c, cd, cm, i97, j97, test
      integer ivec

      do ivec = 1, len
         uni = u(i97) - u(j97)
         if( uni .lt. 0.0 ) uni = uni + 1.0
         u(i97) = uni
         i97 = i97 - 1
         if(i97 .eq. 0) i97 = 97
         j97 = j97 - 1
         if(j97 .eq. 0) j97 = 97
         c = c - cd
         if( c .lt. 0.0 ) c = c + cm
         uni = uni - c
         if( uni .lt. 0.0 ) uni = uni + 1.0
         rvec(ivec) = uni
ENDDO
      return
      end
!============================================================================
!!	@param	tx	[OUT]	number of seconds
!!	@param	ty	[OUT]	number of miliseconds
!!	@author P.G. Khalatur
!!
subroutine ctime00 ( tx, ty )
!******************************************************************
!  returns the number of seconds and hundredths of seconds elapsed
!  since midnight
!******************************************************************
      save
      integer ih, im, is, ihu
      integer hms(3)
      data    key /0/

!     use the Unix-style "itime" and "idate" intrinsic functions,
!     this code works for all compilers except those noted below

!     call gettim (ih, im, is, ihu)
!     tx = ih*3600.0 + im*60 + is + ihu/100.0

!     call System_Clock (ih, im, is)
!     tx = real (ih) / real (im)

      call itime (hms)
      ih = hms(1)
      im = hms(2)
      is = hms(3)

      tx = (ih * 3600.0 + im * 60.0 + is) / 60.

      if (key .gt. 0) go to 10
      key = 1
      tr  = tx
      ihu = 0
   10 tx  = tx - tr
      ty  = tr
      return
      end
      
End Module BoxMuller
