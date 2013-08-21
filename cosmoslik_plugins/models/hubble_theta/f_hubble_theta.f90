module f_hubble_theta

contains

    function hubble2theta(hub, omb, omc, omv, omn, w, num_nu, num_nu_mass)
        use hubble_theta_convert
        real(8) :: hub, omb, omc, omv, omn, w, num_nu, num_nu_mass
        real(8) :: hubble2theta

        hubble2theta = f_hubble2theta(hub, omb, omc, omv, omn, w, num_nu, num_nu_mass)
    end function

    function theta2hubble(theta, ombh2, omch2, omk, omnh2, w, num_nu, num_nu_mass)
        use hubble_theta_convert
        real(8) :: theta, ombh2, omch2, omk, omnh2, w, num_nu, num_nu_mass
        real(8) :: theta2hubble

        theta2hubble = f_theta2hubble(theta, ombh2, omch2, omk, omnh2, w, num_nu, num_nu_mass)
    end function

end module
