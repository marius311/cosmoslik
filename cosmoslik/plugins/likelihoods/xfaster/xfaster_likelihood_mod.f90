
MODULE xfaster_likelihood_mod
  !Code     : Xfaster likelihood mod
  !Copyright: Carlo R. Contaldi
  !Author   : Carlo R. Contladi, Imperial College London
  !Email    : c.contaldi@imperial.ac.uk
  !Date     : 08 April 2008
  !Author1  : Graca Rocha
  !Email    : graca@caltech.edu
  !Date     : 08 April 2008
  !
  !This code reads a "_save" file produced by an xfaster run to calculate the 
  !power spectrum from data. It uses the stored information to compute an
  !approximate likelihood wrt a model passed to it by e.g. cosmomc.
  !There is no need for window functions or the power spectrum itself, the 
  !inputs are the raw pseudo-cl of the observations (or simulated ones) and 
  !all the info required in xfaster to relate the pseudo-cl to the full-sky
  !cl.
  !
  !Use at your own risk, no guarantee is implied etc etc...
  !
  !mods:
  !04/04/08 - Basic implementation
  !08/04/08 - Fixed matrix algebra, TE inclusion works now!

  IMPLICIT NONE

  REAL(8), PARAMETER, PUBLIC :: myPI    = 3.141592653589793238462643383279502884197

  INTEGER                                 :: lmax_xfaster, nbins,ncl_xfaster,ndim
  INTEGER                                 :: nmaps,ncl_tot
  REAL(8)                                 :: kern_norm

  REAL(4), ALLOCATABLE, DIMENSION(:,:,:)  :: kern,xkern,pkern,mkern
  REAL(4), ALLOCATABLE, DIMENSION(:,:)    :: transf,cl_shape,cl_noise,cl_obs
  REAL(4), ALLOCATABLE, DIMENSION(:)      :: gdensity
  REAL(4), ALLOCATABLE, DIMENSION(:,:)    :: gdensity_multi
  REAL(8), ALLOCATABLE, DIMENSION(:,:)    :: fsky_multi

  REAL(4), ALLOCATABLE, DIMENSION(:,:,:)  :: cbl
  REAL(4), ALLOCATABLE, DIMENSION(:,:,:)  :: xcbl
  REAL(4), ALLOCATABLE, DIMENSION(:,:,:)  :: pcblee,mcblee
  REAL(4), ALLOCATABLE, DIMENSION(:,:,:)  :: pcblbb,mcblbb
  REAL(4), ALLOCATABLE, DIMENSION(:,:,:)  :: xtbcbl
  REAL(4), ALLOCATABLE, DIMENSION(:,:,:)  :: ebcbl

  INTEGER, ALLOCATABLE, DIMENSION(:,:)    :: ltab
  INTEGER, ALLOCATABLE, DIMENSION(:)      :: nbins_type
  
  INTEGER :: low_l_cut = 31

  LOGICAL                                 :: xfaster_init_call = .TRUE.

CONTAINS
  !==================================================================
  SUBROUTINE xfaster_init(save_file)
    !==================================================================
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: save_file
    INTEGER :: recl

    WRITE(*,'(a)') 'Initializing xfaster likelihood'

    !INQUIRE (IOLENGTH=recl)  lmax_xfaster, nbins, ncl_xfaster, ndim, nmaps, ncl_tot
    OPEN(unit=21,file=TRIM(ADJUSTL(save_file)),status='unknown',form='unformatted')
    READ(21) lmax_xfaster, nbins, ncl_xfaster, ndim, nmaps, ncl_tot

    ALLOCATE(kern(0:lmax_xfaster,0:lmax_xfaster,1:nmaps*(nmaps+1)/2))
    ALLOCATE(mkern(0:lmax_xfaster,0:lmax_xfaster,1:nmaps*(nmaps+1)/2))
    ALLOCATE(pkern(0:lmax_xfaster,0:lmax_xfaster,1:nmaps*(nmaps+1)/2))
    ALLOCATE(xkern(0:lmax_xfaster,0:lmax_xfaster,1:nmaps*(nmaps+1)/2))

    ALLOCATE(transf(0:lmax_xfaster,1:ncl_tot))
    ALLOCATE(cl_noise(0:lmax_xfaster,1:ncl_tot))      
    ALLOCATE(cl_obs(0:lmax_xfaster,1:ncl_tot))
    ALLOCATE(cl_shape(0:lmax_xfaster,1:ncl_tot))
    ALLOCATE(gdensity(0:lmax_xfaster))
    ALLOCATE(gdensity_multi(0:lmax_xfaster,1:nmaps))
    ALLOCATE(fsky_multi(1:nmaps,1:ndim))
    ALLOCATE(ltab(0:1,0:nbins-1))
    ALLOCATE(nbins_type(1:7))

    read(21) ltab, nbins_type, cl_noise, cl_obs, transf
    read(21) gdensity, gdensity_multi, fsky_multi
    read(21) kern_norm
    read(21) kern,mkern,pkern,xkern
    CLOSE(21)

    ALLOCATE(cbl(0:nbins_type(1)-1, 0:lmax_xfaster,1:nmaps*(nmaps+1)/2))
    ALLOCATE(xcbl(0:nbins_type(4)-1, 0:lmax_xfaster,1:nmaps*(nmaps+1)/2))
    ALLOCATE(pcblee(0:nbins_type(2)-1, 0:lmax_xfaster,1:nmaps*(nmaps+1)/2))
    ALLOCATE(mcblee(0:nbins_type(2)-1, 0:lmax_xfaster,1:nmaps*(nmaps+1)/2))
    ALLOCATE(pcblbb(0:nbins_type(3)-1, 0:lmax_xfaster,1:nmaps*(nmaps+1)/2))
    ALLOCATE(mcblbb(0:nbins_type(3)-1, 0:lmax_xfaster,1:nmaps*(nmaps+1)/2))
    ALLOCATE(xtbcbl(0:nbins_type(5)-1, 0:lmax_xfaster,1:nmaps*(nmaps+1)/2))
    ALLOCATE(ebcbl(0:nbins_type(6)-1, 0:lmax_xfaster,1:nmaps*(nmaps+1)/2))

    xfaster_init_call = .FALSE.
    WRITE(*,'(a)') 'Done xfaster initializing'
  END SUBROUTINE xfaster_init

  !==================================================================
  SUBROUTINE xfaster_lnL(cl_tt,cl_te,cl_ee,cl_bb,like)
    !==================================================================
    IMPLICIT NONE

    REAL(8), INTENT(in)  :: cl_tt(2:*),cl_te(2:*),cl_ee(2:*),cl_bb(2:*)
    REAL(8), INTENT(out) :: like

    REAL(4),     DIMENSION(0:lmax_xfaster,1:ncl_tot) :: cl_theory, cl_mod
    REAL(4),     DIMENSION(1:ndim*nmaps,1:ndim*nmaps) :: a,b,c,d,cinv,cmat
    REAL(4),     DIMENSION(1:ndim*nmaps) :: w
    REAL(4),     DIMENSION(1:10*ndim*nmaps) :: work
    REAL(8) :: trace,trace_l
    INTEGER :: i,j,k,l,ib,ib1,l0,l1,info
    INTEGER :: ix,iy,b1,icl,i_block,j_block
    !==================================================================

    cl_shape = 0.d0
    !clone for each map switching from cosmomc order 
    ib = 0
    DO i=1,nmaps
       DO j =i,nmaps 
          cl_shape(2:lmax_xfaster,ib+1) = cl_tt(2:lmax_xfaster)
          IF(ncl_xfaster > 1) THEN
             cl_shape(2:lmax_xfaster,ib+2) = cl_ee(2:lmax_xfaster)
             cl_shape(2:lmax_xfaster,ib+3) = cl_bb(2:lmax_xfaster)
             cl_shape(2:lmax_xfaster,ib+4) = cl_te(2:lmax_xfaster)
             !clear out any noise in cross spectra
             cl_noise(2:lmax_xfaster,ib+4) = 0.d0
             !Set BB transfer == EE one
             transf(2:lmax_xfaster,ib+3) =  transf(2:lmax_xfaster,ib+2)
             IF(ncl_xfaster > 4) THEN
                !not computed by cosmomc
                cl_shape(2:lmax_xfaster,ib+5) = 0.d0
                cl_shape(2:lmax_xfaster,ib+6) = 0.d0
                cl_noise(2:lmax_xfaster,ib+5) = 0.d0
                cl_noise(2:lmax_xfaster,ib+6) = 0.d0
             ENDIF
          ENDIF
          ib = ib + ncl_xfaster
       ENDDO
    ENDDO

    CALL create_matrices_pol_multi

    !cl theory from shape passed from cosmomc
    cl_theory = 0.d0
    ix = 1
    iy = 0
    DO i = 1,nmaps
       DO j = i,nmaps
          DO l = 2, lmax_xfaster
             b1 = 0
             DO icl = 1,ncl_xfaster
                DO ib1 = 0, nbins_type(icl)-1
                   IF(icl==1) THEN
                      cl_theory(l,iy+1) = cl_theory(l,iy+1) + cbl(ib1,l,ix)
                   ELSEIF(icl==2) THEN
                      cl_theory(l,iy+2) = cl_theory(l,iy+2) + pcblee(ib1,l,ix) + &
                           mcblbb(ib1,l,ix)
                   ELSEIF(icl==3) THEN
                      cl_theory(l,iy+3) = cl_theory(l,iy+3) + pcblbb(ib1,l,ix) + &
                           mcblee(ib1,l,ix)
                   ELSEIF(icl==4) THEN
                      cl_theory(l,iy+4) = cl_theory(l,iy+4) + xcbl(ib1,l,ix)
                   ELSEIF(icl==5) THEN
                      cl_theory(l,iy+5) = cl_theory(l,iy+5) + xtbcbl(ib1,l,ix)
                   ELSEIF(icl==6) THEN
                      cl_theory(l,iy+6) = cl_theory(l,iy+6) + ebcbl(ib1,l,ix)
                   ENDIF
                   b1 = b1 + 1
                ENDDO
             ENDDO
          ENDDO
          ix = ix + 1
          iy = iy + ncl_xfaster
       ENDDO
    ENDDO

    trace = 0.d0
    DO l=low_l_cut,lmax_xfaster
       cl_mod(l,1:ncl_tot) = cl_theory(l,1:ncl_tot) + cl_noise(l,1:ncl_tot)

       cinv = 0.d0
       cmat = 0.d0
       a = 0.d0
       !make 3x3 block for each map
       ix = 1
       iy = 0
       DO i = 1,nmaps
          DO j = i,nmaps
             i_block = ndim*(i-1)
             j_block = ndim*(j-1)
             cmat(i_block + 1,j_block + 1) = cl_mod(l, iy+1)
             a(i_block + 1,j_block + 1) = cl_obs(l,iy+1)
             IF(ncl_xfaster>1) THEN
                cmat(i_block + 2,j_block + 2) = cl_mod(l, iy+2)
                cmat(i_block + 3,j_block + 3) = cl_mod(l, iy+3)
                cmat(i_block + 1,j_block + 2) = 0.d0
                a(i_block + 1,j_block + 2) = 0.d0
                !IF(l<te_lmax_xfaster) THEN
                   cmat(i_block + 1,j_block + 2) = cl_mod(l, iy+4)
                   a(i_block + 1,j_block + 2) = cl_obs(l,iy+4)
                !ENDIF
                cmat(i_block + 2,j_block + 1) = cmat(i_block + 1,j_block + 2)
                a(i_block + 2,j_block + 2) = cl_obs(l,iy+2)
                a(i_block + 3,j_block + 3) = cl_obs(l,iy+3)
                a(i_block + 2,j_block + 1) = a(i_block + 1,j_block + 2)
                IF(ncl_xfaster > 4) THEN
                   cmat(i_block + 1,j_block + 3) = cl_mod(l, iy+5)
                   cmat(i_block + 3,j_block + 1) = cmat(i_block + 1,j_block + 3)
                   cmat(i_block + 2,j_block + 3) = cl_mod(l, iy+6)
                   cmat(i_block + 3,j_block + 2) = cmat(i_block + 2,j_block + 3)
                   a(i_block + 1,j_block + 3) = cl_obs(l,iy+5) 
                   a(i_block + 3,j_block + 1) = a(i_block + 1,j_block + 3)
                   a(i_block + 2,j_block + 3) = cl_obs(l,iy+6) 
                   a(i_block + 3,j_block + 2) = a(i_block + 2,j_block + 3)
                ENDIF
             ENDIF
             ix = ix + 1
             iy = iy + ncl_xfaster
          ENDDO
       ENDDO

       !symmetrize
       DO i = 1,ndim*nmaps
          DO j = 1,i-1
             cmat(i,j) = cmat(j,i)
             a(i,j) = a(j,i)
          ENDDO
       ENDDO

       !a = cmat

       !call sPOTRF( 'L', ndim*nmaps, cmat, ndim*nmaps, INFO )
       !call SPOTRI( 'L', ndim*nmaps, cmat, ndim*nmaps, INFO )
       CALL SSYEV( 'V', 'L', ndim*nmaps, cmat, ndim*nmaps, W, WORK, 10*ndim*nmaps, INFO )

       !log and inverse of C^theory
       b = 0.d0
       d = 0.d0
       DO i=1,ndim*nmaps
          b(i,i) = 1.d0/w(i)
          d(i,i) = LOG(w(i))
       ENDDO
       cinv = MATMUL(MATMUL(cmat,b),TRANSPOSE(cmat))
       d =  MATMUL(MATMUL(cmat,d),TRANSPOSE(cmat))
       b = MATMUL(a,cinv)

       !log of C^obs
       CALL SSYEV( 'V', 'L', ndim*nmaps, a, ndim*nmaps, W, WORK, 10*ndim*nmaps, INFO )
       c = 0.d0
       DO i=1,ndim*nmaps
          c(i,i) = LOG(w(i))
       ENDDO
       c =  MATMUL(MATMUL(a,c),TRANSPOSE(a))

       trace_l = 0.d0
       ix = 1
       DO i = 1,nmaps
          i_block = ndim*(i-1)
          DO k = 1,ndim
             !calculate the likelihood for each l 3x3 matrix
             !likelihood is normalised to be 0.0 at max
             trace_l = trace_l + &
                  fsky_multi(i,k)*gdensity_multi(l,i) *(b(i_block+k,i_block+k) &
                  - c(i_block+k,i_block+k) &
                  + d(i_block+k,i_block+k)-1.d0) 
             ix = ix + 1
          ENDDO
       ENDDO
       trace  = trace + gdensity(l)*(2.d0*l+1.d0)*(trace_l)
    ENDDO

    like = trace

  CONTAINS

    SUBROUTINE create_matrices_pol_multi
      !==================================================================
      !
      ! create the masked projector for the shape spectrum
      !
      ! C_bl = \Sigma_{l'} K_ll' b(l') C_shape(l') X_b(l')
      !
      ! b(l) = transfer(l) * beam(l)^2 * pixel(l)^2
      !=====================================================================
      IMPLICIT NONE

      INTEGER :: ib, l0, l1, nl, lp, lmax_sum
      !==================================================================

      cbl = 0.d0
      IF(ncl_xfaster>1) THEN
         pcblee = 0.d0
         mcblee = 0.d0
         pcblbb = 0.d0
         mcblbb = 0.d0
         xcbl = 0.d0 
         IF(ncl_xfaster>4) THEN
            xtbcbl = 0.d0
            ebcbl = 0.d0
         ENDIF
      ENDIF
      ix = 1
      iy = 0
      DO i=1,nmaps
         DO j=i,nmaps
            b1 = 0
            DO icl = 1,ncl_xfaster
               DO ib = 0, nbins_type(icl)-1
                  l0 = ltab(0,b1)
                  l1 = ltab(1,b1)
                  !lmax_sum = MAX(l1+te_width,l1+(l1-l0+1)*5)
                  !DO lp = l0,l1
                  !max(l0-te_width,2),min(l1+te_width,lmax_xfaster) 
                  DO lp = 2,MIN(l1+(l1-l0+1)*3,lmax_xfaster)
                     IF(icl==1) THEN
                        cbl(ib,lp,ix) = SUM(kern(lp,l0:l1,ix)*transf(l0:l1,iy+1)* &
                             cl_shape(l0:l1,1))*kern_norm
                     ELSEIF(icl==2) THEN
                        pcblee(ib,lp,ix) = SUM(pkern(lp,l0:l1,ix)*transf(l0:l1,iy+2)* &
                             cl_shape(l0:l1,2))*kern_norm
                        mcblee(ib,lp,ix) = SUM(mkern(lp,l0:l1,ix)*transf(l0:l1,iy+3)* &
                             cl_shape(l0:l1,2))*kern_norm
                     ELSEIF(icl==3) THEN
                        pcblbb(ib,lp,ix) = SUM(pkern(lp,l0:l1,ix)*transf(l0:l1,iy+3)* &
                             cl_shape(l0:l1,3))*kern_norm
                        mcblbb(ib,lp,ix) = SUM(mkern(lp,l0:l1,ix)*transf(l0:l1,iy+2)* &
                             cl_shape(l0:l1,3))*kern_norm
                     ELSEIF(icl==4) THEN
                        xcbl(ib,lp,ix) = SUM(xkern(lp,l0:l1,ix)*transf(l0:l1,iy+4)* &
                             cl_shape(l0:l1,4))*kern_norm
                     ELSEIF(icl==5) THEN
                        xtbcbl(ib,lp,ix) = SUM(xkern(lp,l0:l1,ix)*transf(l0:l1,iy+5)* &
                             cl_shape(l0:l1,5))*kern_norm
                     ELSEIF(icl==6) THEN
                        ebcbl(ib,lp,ix) = SUM((pkern(lp,l0:l1,ix)-mkern(lp,l0:l1,ix))* &
                             transf(l0:l1,iy+6)* cl_shape(l0:l1,6))*kern_norm
                     ENDIF
                  ENDDO
                  b1 = b1 + 1
               ENDDO
            ENDDO
            ix = ix + 1
            iy = iy + ncl_xfaster
         ENDDO
      ENDDO
      RETURN
    END SUBROUTINE create_matrices_pol_multi

  END SUBROUTINE xfaster_lnL

END MODULE xfaster_likelihood_mod

