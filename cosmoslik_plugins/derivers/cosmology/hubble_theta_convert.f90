! -*- f90 -*-

!! This is a routine to convert between Theta and Hubble and vice-versa.  
!! It is taken from the CAMB and CosmoMC codes, but avoids the need for the 
!! CAMB/CMB types. Everything interesting was there in CAMB/CosmoMC already.
!! Important: This also assumes that there are no massive neutrinos.

module hubble_theta_convert
implicit none
integer, parameter :: r8b = selected_real_kind(12,200)
real(r8b), parameter :: pi = 3.1415926535897932384626433832795_r8b
private

         !Constants in SI units
   real(r8b), parameter :: Mpc = 3.085678e22_r8b, G = 6.6742e-11_r8b,   &
                         kappa = 8.0_r8b*pi*G,    c = 2.99792458e8_r8b, & 
                   sigma_boltz = 5.67051e-8_r8b

   real(r8b), parameter :: tcmb = 2.725_r8b,  &
                        omkflat = 5.0e-7_r8b, & !! When to consider model flat
                      TOLERANCE = 1.0e-5_r8b

   real(r8b) :: omegab, omegac, omegav, omegak, omegan, hubble, wlambda, &
                omegabh2, omegach2, omeganh2, Num_Nu_massless, Num_Nu_massive, &
                theta_hold, grhom, grhog, grhor, grhob, grhoc, grhov,    &
                grhok, grhornomass, grhormass

   logical :: flat, integration_error, nu_initialized = .false.

      !! extra massive neutrino parameters (see CAMB:modules:MassiveNu)
   integer, parameter  :: nrhopn = 2000
   real(r8b), parameter :: am_min = 0.01_r8b,                     &
                           am_max = 600.0_r8b,                    &
                          am_minp = am_min*1.1_r8b,               &
                          am_maxp = am_max*0.9_r8b,               &
                            const = 7.0_r8b/120.0_r8b*pi**4,      &
                           const2 = 5.0_r8b/7.0_r8b/pi**2,        &
                            zeta3 = 1.2020569031595942853997_r8b, &
                            zeta5 = 1.0369277551433699263313_r8b
   real(r8b) :: dlnam, nu_mass
   real(r8b), dimension(1:nrhopn) :: r1, dr1
   logical :: has_massive_nu


   public :: f_hubble2theta, f_theta2hubble

contains

   !!=========================================================================!!

   function f_hubble2theta(hub, omb, omc, omv, omn, w, num_nu, num_nu_mass)

         real(r8b), intent(in) :: hub, omb, omc, omv, omn, w, num_nu, num_nu_mass
         real(r8b) :: f_hubble2theta

         real(r8b) :: ombh2, ommh2, rs, DA, astar, omdmh2, atol


         if (hub < 2.0) write(*,'(A,A,A,E15.5)') 'hubble2theta: warning: ', &
                        'Your hub may be in the wrong units. They should ', &
                        'be km/s/Mpc but you have hub =', hub
         hubble = hub
         omegab = omb
         omegac = omc
         omegav = omv
         omegan = omn
         wlambda = w
         Num_Nu_massless = num_nu
         Num_Nu_massive = num_nu_mass

         omegak = 1.0_r8b - omegab - omegac - omegav - omegan
         flat = abs(omegak) <= omkflat

         ombh2 = omb*(hubble/100.0_r8b)**2
         ommh2 = (omb+omc+omn)*(hubble/100.0_r8b)**2

         astar = scalefactor_at_decoupling(ombh2,ommh2)

         call fixrhos(hub)
         if (has_massive_nu .and. (.not. nu_initialized)) call nu_init

         !atol = 1e-8
         rs = rombint(dsoundda,1d-8,astar,TOLERANCE)
         DA = AngularDiameterDistance(astar)/astar


         f_hubble2theta = rs/DA

   end function

   !!=========================================================================!!

   function f_theta2hubble(theta, ombh2, omch2, omk, omnh2, w, num_nu, &
                         num_nu_mass)

         real(r8b), intent(in) :: theta, ombh2, omch2, omk, omnh2, w, num_nu, &
                                  num_nu_mass
         logical :: error
         real(r8b) :: f_theta2hubble

         real(r8b), parameter :: h_guess = 70.0_r8b,  &
                                h_upper = 181.0_r8b, &
                                h_lower = 19.0_r8b
   !      real(r8b), parameter :: tol = 1.0e-10_r8b
         real(r8b) :: htemp, err


         if (theta > 1.0) write(*,'(A,A,A,E15.5)') 'theta2hubble: warning: ', &
                        'Your theta is stange. You should be passing theta ', &
                        'not 100*theta? You have theta = ', theta


         theta_hold = theta
         omegabh2 = ombh2
         omegach2 = omch2
         omegak = omk
         omeganh2 = omnh2
         wlambda = w
         Num_Nu_massless = num_nu
         Num_Nu_massive = num_nu_mass

         !!=== Assuming prior 20 < H < 180.
         integration_error = .false.
         nu_initialized = .false.
         htemp = minimize(h_minimum,h_lower,h_guess,h_upper,TOLERANCE)
         err = h_minimum(htemp)
         f_theta2hubble = htemp

         nu_initialized = .false.

   end function

   !!=========================================================================!!
      !! Only private routines are below.
   !!=========================================================================!!

   function h_minimum(h)

         real(r8b), intent(in) :: h
         real(r8b) :: h_minimum 

         real(r8b) :: tmp


         hubble = h
         omegab = omegabh2/(h/100.0_r8b)**2
         omegac = omegach2/(h/100.0_r8b)**2
         omegan = omeganh2/(h/100.0_r8b)**2
         omegav = 1.0_r8b - omegab - omegac - omegak - omegan
         if (.not. integration_error) then
            tmp = f_hubble2theta(h,omegab,omegac,omegav,omegan,wlambda,&
                               Num_Nu_massless,Num_Nu_massive)
         else
            tmp = theta_hold
         end if
         h_minimum = (tmp - theta_hold)**2

   end function h_minimum

   !!=========================================================================!!

   function scalefactor_at_decoupling(ombh2_in,ommh2_in)
      ! Given Omega_m h^2 and Omega_b h^2, returns 
      ! the scale factor at decoupling.
      ! (from cmbwarp, physical parameter routine)

         real(r8b), intent(in) :: ombh2_in, ommh2_in
         real(r8b) :: scalefactor_at_decoupling

         real(r8b) :: zstar


         zstar = 1048*(1+0.00124*ombh2_in**(-0.738))*(1+ &
                 (0.0783*ombh2_in**(-0.238)/(1+39.5*ombh2_in**0.763)) * &
                 ommh2_in**(0.560/(1+21.1*ombh2_in**1.81)))
         scalefactor_at_decoupling = 1.0_r8b/(1.0_r8b+zstar)

   end function scalefactor_at_decoupling

   !!=========================================================================!!

   subroutine fixrhos(hub)
     !! See CAMB: modules.f90: CAMBParams_Set for more information.
         real(r8b), intent(in) :: hub


         grhom = 3.0_r8b*hub**2/(c/1000.0_r8b)**2
         grhob = grhom*omegab
         grhoc = grhom*omegac
         grhov = grhom*omegav
         grhok = grhom*omegak

         grhog = kappa/c**2*4*sigma_boltz/c**3*tcmb**4*Mpc**2
         grhor = 7._r8b/8*(4._r8b/11)**(4._r8b/3)*grhog

         if (omegan == 0) then
            has_massive_nu = .false.
            Num_Nu_massless = Num_Nu_massless + Num_Nu_massive
            Num_Nu_massive = 0.0_r8b
         else
            nu_mass = const/(1.5_r8b*zeta3)*grhom/grhor*omegan/Num_Nu_massive
            has_massive_nu = .true.
         end if

         grhornomass = grhor*Num_Nu_massless
         grhormass = grhor*Num_Nu_massive


   end subroutine fixrhos

   !!=========================================================================!!

   function AngularDiameterDistance(a)

         real(r8b), intent(in) :: a
         real(r8b) ::  AngularDiameterDistance

         real(r8b) :: temp, r


         temp = rombint(dtauda,a,1.0_r8b,TOLERANCE)
         if (flat) then
            AngularDiameterDistance = a*temp
         else
            r = -omegak/((c/1000)/hubble)**2
            r = 1._r8b/sqrt(abs(r))
            if (omegak < 0.0_r8b) then
               temp = sin(temp/r)
            else
               temp = sinh(temp/r)
            end if
            AngularDiameterDistance = a*r*temp
         end if

   end function AngularDiameterDistance

   !!=========================================================================!!

   function dsoundda(a)

         real(r8b), intent(in) :: a
         real(r8b) :: dsoundda

         real(r8b) :: R, cs, tmp

         R = 3.0e4_r8b*a*omegab*(hubble/100.0_r8b)**2
         cs = 1.0_r8b/sqrt(3.0_r8b*(1.0_r8b+R))
         tmp =  dtauda(a) 
         dsoundda = tmp*cs

   end function dsoundda

   !!=========================================================================!!

   function dtauda(a)
         !get d tau / d a

         real(r8b), intent(in) :: a
         real(r8b) ::  dtauda

         integer :: i 
         real(r8b) :: grhoa2, a2, rhonu

         a2=a**2        
         grhoa2=grhok*a2+(grhoc+grhob)*a+grhog+grhornomass
         grhoa2=grhoa2+grhov*a**(1.0_r8b-3.0_r8b*wlambda)

         if (has_massive_nu) then
            rhonu = get_rhonu(nu_mass*a)
            grhoa2 = grhoa2 + rhonu*grhormass
         end if

         if (grhoa2 < 0.0_r8b) then
            integration_error = .true.
            dtauda = 0.0_r8b
         else
            dtauda = sqrt(3.0_r8b/grhoa2)
         end if

   end function dtauda

   !!=========================================================================!!

   subroutine nu_init

         integer :: i
         real(r8b) :: rhonu, am
         real(r8b), dimension(1:nrhopn) :: spline_data


         dlnam = -(log(am_min/am_max))/real(nrhopn-1,r8b)
         do i = 1, nrhopn
            am = am_min*exp(real(i-1,r8b)*dlnam)
            call nuRhoPres(am,rhonu)
            r1(i) = log(rhonu)
         end do

         call splini(spline_data,nrhopn)
         call splder(r1,dr1,nrhopn,spline_data)

         nu_initialized = .true.

   contains

      subroutine nuRhoPres(am,rhonu)

            real(r8b), intent(in) :: am
            real(r8b), intent(out) :: rhonu

            integer, parameter :: nq=1000
            real(r8b),  parameter :: qmax=30.0_r8b
      
            integer :: i 
            real(r8b) :: q, aq, v, adq, aqdn
            real(r8b), dimension(1:nq+1) :: dum1


            adq = qmax/real(nq,r8b)
            dum1(1) = 0.0_r8b
            do i = 1, nq
               q = real(i,r8b)*adq
               aq = am/q
               v = 1.0_r8b/sqrt(1.0_r8b+aq*aq)
               aqdn = adq*q*q*q/(exp(q)+1.0_r8b)
               dum1(i+1) = aqdn/v
            end do
            call splint(dum1,rhonu,nq+1)
            rhonu = (rhonu+dum1(nq+1)/adq)/const

        end subroutine nuRhoPres

   end subroutine nu_init

   !!=========================================================================!!

   function get_rhonu(nu_mass) result(rhonu)

         real(r8b), intent(in) :: nu_mass
         real(r8b) :: rhonu

         integer :: i
         real(r8b) :: am, d


         if (nu_mass <= am_minp) then
            rhonu = 1.0_r8b + const2*am**2
            return
         else if (nu_mass >= am_maxp) then
            rhonu = 3.0_r8b / (2.0_r8b*const) * &
                    (zeta3*nu_mass + (15.0_r8b*zeta5)/2.0_r8b/nu_mass)
            return
         end if
        
         d = log(nu_mass/am_min)/real(dlnam,r8b) + 1.0_r8b
         i = int(d)
         d = d - real(i,r8b)
       
         rhonu = r1(i)+d*(dr1(i)+d*(3._r8b*(r1(i+1)-r1(i))-2._r8b*dr1(i) &
               -dr1(i+1)+d*(dr1(i)+dr1(i+1)+2._r8b*(r1(i)-r1(i+1)))))
         rhonu = exp(rhonu)

   end function get_rhonu

   !!=========================================================================!!

   function minimize(func, int_left, int_guess, int_right, tol) result(minimum)

         real(r8b), intent(in) :: int_left, int_right, int_guess, tol
         real(r8b) :: minimum

         integer, parameter :: MAX_ITERS = 10000
!         real(r8b), parameter :: RATIO = 0.5_r8b*(3.0_r8b - sqrt(5.0_r8b)), &
!                                 CONST = 1.0_r8b - RATIO
         integer :: i
         real(r8b) :: x0, x1, x2, x3, f1, f2, RATIO, CONST

         interface
            function func(x)
            implicit none
            integer, parameter :: r8b = selected_real_kind(12,200)
               real(r8b), intent(in) :: x
               real(r8b) :: func
            end function func
         end interface


         RATIO = 0.5_r8b*(3.0_r8b - sqrt(5.0_r8b))
         CONST = 1.0_r8b - RATIO

         x0 = int_left
         x3 = int_right
         if (abs(int_right - int_guess) > abs(int_guess - int_left)) then
            x1 = int_guess
            x2 = int_guess + CONST*(int_right - int_guess)
         else
            x1 = int_guess - CONST*(int_guess - int_left)
            x2 = int_guess
         end if
         f1 = func(x1)
         f2 = func(x2)
         
         do i = 1, MAX_ITERS
            if (x3-x0 <= tol*(abs(x1)+abs(x2))) exit
            if (f2 < f1) then
               x0 = x1
               x1 = x2
               x2 = RATIO*x1 + CONST*x3 !x2=x1-(1-RATIO)*(x3-x1)
               f1 = f2
               f2 = func(x2)
            else
               x3 = x2
               x2 = x1
               x1 = RATIO*x2 + CONST*x0 !x1=x2-(1-RATIO)*(x2-x0)
               f2 = f1
               f1 = func(x1)
            end if
         end do

         if (i >= MAX_ITERS) then
            print '(A)', 'minimize failed. Exceeded maximum iterations.'
            stop
         end if

         if (f1 < f2) then
            minimum = x1
         else
            minimum = x2
         end if
            
   end function minimize

   !!=========================================================================!!

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        function rombint(f,a,b,tol)
!  Rombint returns the integral from a to b of using Romberg integration.
!  The method converges provided that f(x) is continuous in (a,b).
!  f must be real(r8b) and must be declared external in the calling
!  routine.  tol indicates the desired relative accuracy in the integral.
!
        implicit none
        integer, parameter :: MAXITER=20,MAXJ=5
        dimension g(MAXJ+1)
        real(r8b) f
        external f
        real(r8b) :: rombint
        real(r8b), intent(in) :: a,b,tol
        integer :: nint, i, k, jmax, j
        real(r8b) :: h, gmax, error, g, g0, g1, fourj
!
        h=0.5_r8b*(b-a)
        gmax=h*(f(a)+f(b))
        g(1)=gmax
        nint=1
        error=1.0e20_r8b
        i=0
10        i=i+1
          if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
            go to 40
!  Calculate next trapezoidal rule approximation to integral.
          g0=0._r8b
            do 20 k=1,nint
            g0=g0+f(a+(k+k-1)*h)
20        continue
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint=nint+nint
          jmax=min(i,MAXJ)
          fourj=1._r8b
            do 30 j=1,jmax
!  Use Richardson extrapolation.
            fourj=4._r8b*fourj
            g1=g0+(g0-g(j))/(fourj-1._r8b)
            g(j)=g0
            g0=g1
30        continue
          if (abs(g0).gt.tol) then
            error=1._r8b-gmax/g0
          else
            error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
        go to 10
40      rombint=g0
        if (i.gt.MAXITER.and.abs(error).gt.tol)  then
          write(*,*) 'Warning: Rombint failed to converge; '
          write (*,*)'integral, error, tol:', rombint,error, tol
        end if

        end function rombint



       subroutine splder(y,dy,n, g)
!  Splder fits a cubic spline to y and returns the first derivatives at
!  the grid points in dy.  Dy is equivalent to a 4th-order Pade
!  difference formula for dy/di.
        implicit none
        integer, intent(in) :: n
        real(r8b), intent(in) :: y(n),g(n)
        real(r8b), intent(out) :: dy(n)
        integer :: n1, i
        real(r8b) :: f(n)

        n1=n-1
!  Quartic fit to dy/di at boundaries, assuming d3y/di3=0.
        f(1)=(-10._r8b*y(1)+15._r8b*y(2)-6._r8b*y(3)+y(4))/6._r8b
        f(n)=(10._r8b*y(n)-15._r8b*y(n1)+6._r8b*y(n-2)-y(n-3))/6._r8b
!  Solve the tridiagonal system
!  dy(i-1)+4*dy(i)+dy(i+1)=3*(y(i+1)-y(i-1)), i=2,3,...,n1,
!  with dy(1)=f(1), dy(n)=f(n).
        do i=2,n1
          f(i)=g(i)*(3._r8b*(y(i+1)-y(i-1))-f(i-1))
        end do
        dy(n)=f(n)
        do i=n1,1,-1
          dy(i)=f(i)-g(i)*dy(i+1)
        end do

        end subroutine splder
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine splini(g,n)
!  Splini must be called before splder to initialize array g in common.
        implicit none
        integer, intent(in) :: n
        real(r8b), intent(out):: g(n)
        integer :: i

        g(1)=0._r8b
        do i=2,n
          g(i)=1/(4._r8b-g(i-1))
        end do
        end subroutine splini


       subroutine splint(y,z,n)
!  Splint integrates a cubic spline, providing the output value
!  z = integral from 1 to n of s(i)di, where s(i) is the spline fit
!  to y(i).
!
        implicit none
        integer, intent(in) :: n
        real(r8b), intent(in) :: y(n)
        real(r8b), intent(out) :: z

        integer :: n1
        real(r8b) :: dy1, dyn
!
        n1=n-1
!  Cubic fit to dy/di at boundaries.
!       dy1=(-11._r8b*y(1)+18._r8b*y(2)-9._r8b*y(3)+2._r8b*y(4))/6._r8b
        dy1=0._r8b
        dyn=(11._r8b*y(n)-18._r8b*y(n1)+9._r8b*y(n-2)-2._r8b*y(n-3))/6._r8b
!
        z=0.5d0*(y(1)+y(n))+(dy1-dyn)/12._r8b
        z= z + sum(y(2:n1))
        end subroutine splint

end module hubble_theta_convert

