
module pywmap

use wmap_likelihood_7yr
use WMAP_OPTIONS
use WMAP_UTIL

logical :: pywmap_initalized

contains

    subroutine WMAPInit(tt_min,tt_max,te_min,te_max,data_dir)
        integer :: ttmin, ttmax, temin, temax

		character(len=*) :: data_dir
        WMAP_data_dir = trim(data_dir)
        ttmin = tt_min
        ttmax = tt_max
        temin = te_min
        temax = te_max

        call wmap_likelihood_init
        pywmap_initalized = .true.

    end subroutine


    ! To call from Python use WMAPLnLike(cl_tt) where cl_tt[0] is ell=2
    function WMAPLnLike(cltt, clte, clee, clbb, lmax)
      real(8), dimension(10) :: WMAPLnLike
      real(8), dimension(num_WMAP) :: like
      real(8), dimension(2:lmax) :: cltt,clte,clee,clbb

      WMAPLnLike = 1e30
      if (lmax < ttmax) then
         print *, "WMAPLnLike: Please provide c_ell's up to", ttmax
         return
      end if

      if (.not. pywmap_initalized) then
         print *, "WMAPLnLike: WMAP likelihood called before initialization"
         return
      end if
         
      call wmap_likelihood_compute(cltt(2:ttmax),clte(2:ttmax),clee(2:ttmax),clbb(2:ttmax),like)
      
      if (wmap_likelihood_ok) then
         WMAPLnLike = 0
         WMAPLnLike(1:num_WMAP) = like
      else
         print *, "WMAPLnLike: Error computing WMAP likelihood."
         call wmap_likelihood_error_report
      end if

    end function

end module
