! -*- f90 -*-
! atmosphere-ocean carbon cycle box model
! mick follows, march 2015
! convert matlab to f90 - march/june 2016
! significant work by jonathan lauderale june 2016-oct 2019
! refresh, parallelization, and expansion by jonathan lauderdale 2020-
!
!         --------------------------------
!         |          ATMOSPHERE          |
!         |                              |
!         --------------------------------
!         | SUR- |                       |
!         | FACE ->   SURFACE  2         |
!         |  1   |                   |   |
!         |---^----------------------v---|
!         |   |                          | 
!         |          DEEP OCEAN  3       |
!         |                              | 
!         |                              | 
!         --------------------------------
!=======================================================================
       PROGRAM MICROCOSM_MODEL_ENSEMBLE
!=======================================================================

       USE MOD_PRECISION
       USE MOD_BOXES
       USE MOD_MODELMAIN_ENSEMBLE
       USE MOD_SUBROUTINES

       IMPLICIT NONE

       INTEGER :: outstepmax, id, id0
       
       REAL(KIND=wp) ::                                                &
            maxyears,                                                  &
            outputyears,                                               &
            m2deg,                                                     &
!            gaovla_opt,                                                &
            ligphi_in,                                                 &
            lt_lifein,                                                 &
            alpha_yr,                                                  &
            atpco2in,                                                  &
            psi_in,                                                    &
            dif_in,                                                    &
            lt_lifetime_factor

! Input arrays (nbox dimensionesix)
       REAL(KIND=wp), dimension (nbox) ::                              & 
            dx,                                                        &
            dy,                                                        &
            dz,                                                        &
            depth,                                                     &
            latitude,                                                  &
            thin,                                                      & 
            sain,                                                      &
            cain,                                                      & 
            alin,                                                      & 
            phin,                                                      & 
            niin,                                                      & 
            fein,                                                      & 
            ltin,                                                      & 
            fe_input,                                                  &
            dldz_in,                                                   &
            wind_in,                                                   &
            fopen_in,                                                  &
            eratio_in,                                                 &
            pbin,                                                      &
            ldocin,                                                    &
            rCLig_in,                                                  &
            pge_in

       REAL(KIND=wp), dimension (nbox, nbox) ::                        & 
            Kin,                                                       &
            Rin,                                                       &
            Pin
            
! Output arrays (nbox, by timestep dimensionesix)
       REAL(KIND=wp), dimension (:,:), allocatable ::                  &
            thout,                                                     &
            sout,                                                      &
            cout,                                                      &
            aout,                                                      &
            pout,                                                      &
            nout,                                                      &
            fout,                                                      &
            lout,                                                      &
            expout,                                                    &
            ocpco2out,                                                   &
            pbout,                                                     &
            ldocout                                                  
            
       REAL(KIND=wp), dimension (:), allocatable   ::                  &
            tout,                                                      &
            psout,                                                     &
            atpco2out

       INTEGER, dimension (:), allocatable   ::                        &
            nlout


       ! Define the containers for the parameters that are cycled through
       ! in the different ensemble runs
       ! Dependent on the type of ensemble one wants to run, this must
       ! be adjusted
       REAL(KIND=wp), dimension (:), allocatable ::                    &
            rFeC_pb_ens,                                               &
            mu0_ens,                                                   &
            m_l_ens,                                                   &
            m_q_ens,                                                   &
            kappa_ens,                                                 &
            kfe_p_ens,                                                 &
            kldoc_p_ens,                                               &
            pge_deep_ens,                                                   &
            phi_ens,                                                   &
            rCLig_ens,                                                 &
            lt_lifein_ens,                                             &
            ligphi_ens,                                                &
            lt_lifetime_factor_ens

       ! Define indices for the different ensemble parameters
       INTEGER ::                                                      &
            irFeC_pb,                                                  &
            imu0,                                                      &
            im_l,                                                      &
            im_q,                                                      &
            ikappa,                                                    &
            ikfe_p,                                                    &
            ikldoc_p,                                                  &
            ipge_deep,                                                 &
            iphi,                                                      &
            irCLig,                                                    &
            ilt_lifein,                                                &
            iligphi,                                                   &
            ilt_lifetime_factor

       ! csv file parameters
       CHARACTER(len=255) :: filename_csv
       CHARACTER(len=10), allocatable :: headers(:)
       INTEGER :: i

       ALLOCATE(headers(45))

              headers(1)  = 'id       '
              headers(2)  = 'dt(s)    '
              headers(3)  = 't(yr)    '
              headers(4)  = 'rFeC_pb  '
              headers(5)  = 'mu0      '
              headers(6)  = 'm_l      '
              headers(7)  = 'm_q      '
              headers(8)  = 'kappa    '
              headers(9)  = 'kfe_p    '
              headers(10) = 'kldoc_p  '
              headers(11) = 'pge_deep '
              headers(12) = 'phi      '
              headers(13) = 'rCLig    '
              headers(14) = 'lt_lifet '
              headers(15) = 'ligphi   '
              headers(16) = 'lt_deepf '
              headers(17) = 'PB(1)    '
              headers(18) = 'PB(2)    '
              headers(19) = 'PB(3)    '
              headers(20) = 'LDOC(1)  '
              headers(21) = 'LDOC(2)  '
              headers(22) = 'LDOC(3)  '
              headers(23) = 'Fe(1)    '
              headers(24) = 'Fe(2)    '
              headers(25) = 'Fe(3)    '
              headers(26) = 'PO4(1)   '
              headers(27) = 'PO4(2)   '
              headers(28) = 'PO4(3)   '
              headers(29) = 'NO3(1)   '
              headers(30) = 'NO3(2)   '
              headers(31) = 'NO3(3)   '
              headers(32) = 'Lig(1)   '
              headers(33) = 'Lig(2)   '
              headers(34) = 'Lig(3)   '
              headers(35) = 'DIC(1)   '
              headers(36) = 'DIC(2)   '
              headers(37) = 'DIC(3)   '
              headers(38) = 'ALK(1)   '
              headers(39) = 'ALK(2)   '
              headers(40) = 'ALK(3)   '
              headers(41) = 'OCPCO2(1)'
              headers(42) = 'OCPCO2(2)'
              headers(43) = 'OCPCO2(3)'
              headers(44) = 'ATPCO2   '
              headers(45) = 'Limit    '

!=======================================================================
       ! NEED TO CHANGE NAME OF OUTPUT FILE HERE
!=======================================================================
       ! Set the output filename
              filename_csv = 'Representative_member_Nr5_Ligandmorelabile.csv'

!=======================================================================
! Time parameters
       ! Input some initial parameters
       maxyears   = 1e4_wp
       outputyears= 1.0_wp
       outstepmax = int((maxyears/outputyears)+1)
!=======================================================================
       
! allocate memory
       allocate ( tout      (outstepmax) )
       allocate ( nlout     (outstepmax) )
       allocate ( psout     (outstepmax) )
       allocate ( atpco2out (outstepmax) )
       allocate ( thout     (outstepmax,nbox) )
       allocate ( sout      (outstepmax,nbox) )
       allocate ( cout      (outstepmax,nbox) )
       allocate ( aout      (outstepmax,nbox) )
       allocate ( pout      (outstepmax,nbox) )
       allocate ( nout      (outstepmax,nbox) )
       allocate ( fout      (outstepmax,nbox) )
       allocate ( lout      (outstepmax,nbox) )
       allocate ( expout    (outstepmax,nbox) )
       allocate ( ocpco2out   (outstepmax,nbox) )
       allocate ( pbout     (outstepmax,nbox) )
       allocate ( ldocout   (outstepmax,nbox) )
!===================================================================================
       ! NEED TO ALLOCATE AS MUCH MEMORY AS PARAMETER VALUES THAT ARE CYCLED THROUGH
       ! NEEDS TO SPECIFIED MANUALLY HERE
!===================================================================================
       allocate ( rFeC_pb_ens  (1:1) )
       allocate ( mu0_ens      (1:1) )
       allocate ( m_l_ens      (1:1) )
       allocate ( m_q_ens      (1:1) )
       allocate ( kappa_ens    (1:1) )
       allocate ( kfe_p_ens    (1:3) )
       allocate ( kldoc_p_ens  (1:3) )
       allocate ( pge_deep_ens (1:1) )
       allocate ( phi_ens      (1:1) )
       allocate ( rCLig_ens    (1:1) )
       allocate ( lt_lifein_ens(1:1) )
       allocate ( ligphi_ens   (1:1) )
       allocate ( lt_lifetime_factor_ens(1:1) )
!====================================================================================

! Initialize input arguements
       thin     =      0._wp
       sain     =     34._wp
       cain     =   2150._wp
       alin     =   2350._wp
!       phin     =      1._wp
       phin     =      2._wp
!       niin     =     16._wp
       niin     =     36._wp
       fein     =      0.0_wp
       ltin     =      0.0_wp
       atpco2in =    280._wp
       pbin     =      1._wp ! in cells µL-1
       ldocin   =      0.0_wp ! in mumol kg-1
       
! Overturning and mixing rates (m3/s)
!       psi_in = 20.e6_wp
       dif_in =  0.e6_wp
       
! Wind speed (m/s)for CO2 gas fluxes
       wind_in      =   0._wp
       
! Open surface fraction in contact with atmoshpere 
!  1 => fully open, <1 => flux impeded (e.g. by sea ice)
       fopen_in     =  0._wp
       
! Gamma over lambda for ligands "optimum" value (Lauderdale et al 2020)
!       gaovla_opt   = 0._wp !4398._wp
!       gaovla_opt   = 4398._wp

! Gamma ligand production rate (in phosphate, not carbon, units)
!       gamma_in     = 5.e-5_wp*106._wp
!       gamma_in     = 5.e-5_wp*106._wp

       ligphi_in = 5e-8_wp
       rCLig = 30.0_wp ! number of carbon atoms per ligand molecule
       
! Lambda ligand lifetime (s)
       lt_lifein    = 0.0_wp !(10._wp * 365._wp * 24._wp * 3600._wp) ! 1._wp/((5.e-5_wp*106._wp/106._wp)/4398._wp)
!       lt_lifein    = 1._wp/((gamma_in/106._wp)/gaovla_opt)
       
! Dust deposition in g Fe m-2 year-1
       fe_input     =  0._wp

! Biological production maximum rate (mol P/yr)
       alpha_yr      = 6.e-6_wp

! Deep ocean box lifetime modifier to capture the gradient due to
! photodegradation near the surface and slower loss in the deep
       dldz_in       =   0._wp

! Export ratio (export production / total production)
       eratio_in    = 0._wp

! Prokaryotic growth efficiency
       pge_in       = 0.15_wp

!=======================================================================
! Parameters to cycle through

! Prokaryotic parameters
! Prokaryotic biomass carbon to iron ratio

       ! rFeC_pb_ens(1) = 20.0_wp * 1.e-6_wp
       rFeC_pb_ens(1) = 40._wp * 1.e-6_wp
       ! rFeC_pb_ens(1) = 80._wp * 1.e-6_wp

! Prokaryotic maximum growth rate

       ! mu0_ens(1)  = 0.0_wp
       ! mu0_ens(1)  = 0.01_wp * 1.0_wp/sperd ! 
       ! mu0_ens(1)  = 0.03_wp  * 1.0_wp/sperd !
       ! mu0_ens(2)  = 0.1_wp  * 1.0_wp/sperd !
       ! mu0_ens(1)  = 0.1_wp  * 1.0_wp/sperd ! 
       ! mu0_ens(1)  = 1.2_wp  * 1.0_wp/sperd !
       ! mu0_ens(2)  = 1.5_wp  * 1.0_wp/sperd !
       ! mu0_ens(3)  = 2.0_wp  * 1.0_wp/sperd !
       mu0_ens(1)  = 6.94444444444444E-6_wp
       ! mu0_ens(2)  = mu0_ens(1) * 0.9_wp
       ! mu0_ens(3)  = mu0_ens(1) * 1.1_wp


       


! Prokaryotic linear mortality rate

       ! already tried m_l = 1 d-1 which is too high, always gives wrong numbers
       ! m_l_ens(1)  = 3.0e-3_wp * 1.0_wp/sperd ! in units of s-1
       ! m_l_ens(1)  = 3.0e-3_wp * 1.0_wp/sperd ! in units of s-1
       ! m_l_ens(2)  = 5.0e-3_wp * 1.0_wp/sperd ! in units of s-1
       ! m_l_ens(3)  = 1.0e-2_wp * 1.0_wp/sperd ! in units of s-1
       ! m_l_ens(1) = 5.0e-5_wp * 1.0_wp/sperd ! in units of s-1
       ! m_l_ens(2)  = 1.0e-4_wp * 1.0_wp/sperd ! in units of s-1
       ! m_l_ens(3)  = 3.0e-4_wp * 1.0_wp/sperd ! in units of s-1
       ! m_l_ens(4)  = 1.0e-3_wp * 1.0_wp/sperd ! in units of s-1
       ! m_l_ens(4)  = 1.0e-6_wp * 1.0_wp/sperd ! in units of s-1
       ! m_l_ens(3)  = 1.0e-7_wp * 1.0_wp/sperd ! in units of s-1
       ! m_l_ens(9) = 0                        ! in units of s-1
       m_l_ens(1) = 3.47222222222222E-9_wp
       ! m_l_ens(2) = m_l_ens(1) * 0.9_wp
       ! m_l_ens(3) = m_l_ens(1) * 1.1_wp





! Prokaryotic quadratic mortality rate
       m_q_ens(1) = 1.0e-20_wp ! treat pb in cells per m3, this is in units of m3 per cell per s
       ! m_q_ens(2) = m_q_ens(1) * 0.9_wp
       ! m_q_ens(3) = m_q_ens(1) * 1.1_wp
       ! m_q_ens(2) = 7.5e-19_wp ! treat pb in cells per m3, this is in units of m3 per cell per s
       ! m_q_ens(3) = 5.0e-19_wp ! treat pb in cells per m3, this is in units of m3 per cell per s
       ! m_q_ens(4) = 3.0e-19_wp ! treat pb in cells per m3, this is in units of m3 per cell per s
       ! m_q_ens(1)  = 1.0e-18_wp ! treat pb in cells per m3, this is in units of m3 per cell per s
       ! m_q_ens(2)  = 2.0e-18_wp ! treat pb in cells per m3, this is in units of m3 per cell per s
       ! m_q_ens(3)  = 5.0e-18_wp ! treat pb in cells per m3, this is in units of m3 per cell per s
       ! m_q_ens(4)  = 1.0e-17_wp ! treat pb in cells per m3, this is in units of m3 per cell per s
       ! m_q_ens(4)  = 1.0e-22_wp ! treat pb in cells per m3, this is in units of m3 per cell per s
       ! m_q_ens(3)  = 1.0e-23_wp ! treat pb in cells per m3, this is in units of m3 per cell per s


! fraction of dead prokaryotic biomass that released as LDOC

       kappa_ens(1) = 1.0_wp
       ! kappa_ens(2) = 0.0_wp
       


! Prokaryotic half saturation constant for iron

       ! kfe_p_ens(1) = 1.0e-10_wp *conv_molkg_molm3
       ! kfe_p_ens(2) = 3.0e-10_wp *conv_molkg_molm3
       kfe_p_ens(1) = 1.0e-9_wp  *conv_molkg_molm3
       ! kfe_p_ens(2) = kfe_p_ens(1) * 0.9_wp
       ! kfe_p_ens(3) = kfe_p_ens(1) * 1.1_wp
       ! kfe_p_ens(1) = 1.0e-10_wp  *conv_molkg_molm3

! Prokaryotic half saturation constant for LDOC

       ! kldoc_p_ens(1) = 3.0e-5_wp *conv_molkg_molm3
       ! kldoc_p_ens(2) = 1.0e-4_wp *conv_molkg_molm3
       ! kldoc_p_ens(1) = 5.0e-4_wp *conv_molkg_molm3
       kldoc_p_ens(1) = 1.0e-3_wp *conv_molkg_molm3
       ! kldoc_p_ens(2) = kldoc_p_ens(1) * 0.9_wp
       ! kldoc_p_ens(3) = kldoc_p_ens(1) * 1.1_wp
       ! kldoc_p_ens(2) = 1.0e-4_wp *conv_molkg_molm3
       ! kldoc_p_ens(4) = 1.0e-3_wp *conv_molkg_molm3
       ! kldoc_p_ens(5) = 1.0e-3_wp *conv_molkg_molm3
       ! kldoc_p_ens(6) = 3.0e-3_wp *conv_molkg_molm3
       ! kldoc_p_ens(7) = 1.0e-2_wp *conv_molkg_molm3

! Prokaryotic growth efficiency factor for multiplication
       pge_deep_ens(1) = 0.05_wp
       ! pge_deep_ens(2) = 0.2_wp
       ! pge_deep_ens(3) = 0.3_wp

! LDOC parameters
       ! phi_ens(1) = 1.4e-1_wp
       ! phi_ens(2) = 1.5e-1_wp
       ! phi_ens(3) = 1.6e-1_wp
       ! phi_ens(4) = 1.7e-1_wp
       ! phi_ens(4) = 1.75e-1_wp
       ! phi_ens(2) = 2.0e-1_wp
       ! phi_ens(1) = 0.0_wp
       ! phi_ens(1) = 1.0e-2_wp
       ! phi_ens(1) = 3.0e-2_wp
       ! phi_ens(2) = 1.0e-1_wp
       ! phi_ens(3) = 1.5e-1_wp
       ! phi_ens(1) = 1.0e-1_wp
       ! phi_ens(2) = 2.0e-1_wp
       phi_ens(1) = 5.0e-2_wp
       ! phi_ens(2) = phi_ens(1) * 0.9_wp
       ! phi_ens(3) = phi_ens(1) * 1.1_wp
       ! phi_ens(3) = 3.0e-1_wp
       ! phi_ens(4) = 4.0e-1_wp
       ! phi_ens(5) = 5.0e-1_wp
       ! phi_ens(2) = 1.0e-6_wp
       ! phi_ens(1) = 1.0e-5_wp
       ! phi_ens(2) = 1.0e-4_wp
       ! phi_ens(1) = 5.0e-3_wp
       ! phi_ens(2) = 2.5e-3_wp
       ! phi_ens(1) = 2.0e-2_wp
       

! Ligand parameters
       ! rCLig_ens(1) = 20.0_wp
       rCLig_ens(1) = 30.0_wp
       ! rCLig_ens(3) = 40.0_wp

! Ligand lifetimes for the surface ocean
       lt_lifein_ens(1) = 0.0_wp ! automatically sets the ad-hoc term to zero
       ! lt_lifein_ens(1) = 1._wp / 24._wp * 365._wp * 24._wp * 3600._wp ! half a month
       ! lt_lifein_ens(2) = 1._wp / 2._wp * 365._wp * 24._wp * 3600._wp ! half a year
       ! lt_lifein_ens(1) = 0.1_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(2) = 0.125892541179417_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(3) = 0.158489319246111_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(4) = 0.199526231496888_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(5) = 0.251188643150958_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(6) = 0.316227766016838_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(7) = 0.398107170553497_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(8) = 0.501187233627272_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(9) = 0.630957344480193_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(10) = 0.794328234724282_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(11) = 1.0_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(12) = 1.258925411794167_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(13) = 1.584893192461114_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(14) = 1.995262314968879_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(15) = 2.511886431509582_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(16) = 3.162277660168379_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(17) = 3.981071705534972_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(18) = 5.011872336272725_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(19) = 6.30957344480193_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(20) = 7.943282347242822_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(21) = 10.0_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(22) = 12.58925411794167_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(23) = 15.84893192461114_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(24) = 19.95262314968879_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(25) = 25.11886431509582_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(26) = 31.62277660168379_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(27) = 39.81071705534972_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(28) = 50.11872336272725_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(29) = 63.09573444801933_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(30) = 79.43282347242822_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(31) = 100.0_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(32) = 125.8925411794167_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(33) = 158.4893192461114_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(34) = 199.5262314968881_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(35) = 251.1886431509579_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(36) = 316.2277660168379_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(37) = 398.1071705534969_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(38) = 501.1872336272719_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(39) = 630.9573444801929_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(40) = 794.3282347242815_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(41) = 1000.0_wp * 365._wp * 24._wp * 3600._wp


! Ligand production rates as fraction of the primary production
       ! ligphi_ens(1) = 1.0e-5_wp
       ! ligphi_ens(1) = 7.5e-5_wp
       ligphi_ens(1) = 1.0e-3_wp
       ! ligphi_ens(2) = ligphi_ens(1) * 0.9_wp
       ! ligphi_ens(3) = ligphi_ens(1) * 1.1_wp
       ! ligphi_ens(1) = 5.0e-5_wp
       ! ligphi_ens(2) = 1.0e-4_wp
       ! ligphi_ens(3) = 1.5e-4_wp
       ! ligphi_ens(3) = 3.0e-4_wp
       ! ligphi_ens(4) = 1.0e-3_wp
       ! ligphi_ens(5) = 3.0e-3_wp
       ! ligphi_ens(5) = 5.0e-3_wp

       ! ligphi_ens(1) = 1.0e-1_wp
       ! ligphi_ens(2) = 7.9432823472e-2_wp
       ! ligphi_ens(3) = 6.3095734448e-2_wp
       ! ligphi_ens(4) = 5.0118723363e-2_wp
       ! ligphi_ens(5) = 3.9810717055e-2_wp
       ! ligphi_ens(6) = 3.1622776602e-2_wp
       ! ligphi_ens(7) = 2.5118864315e-2_wp
       ! ligphi_ens(8) = 1.995262315e-2_wp
       ! ligphi_ens(9) = 1.5848931925e-2_wp
       ! ligphi_ens(10) = 1.2589254118e-2_wp
       ! ligphi_ens(11) = 1.0e-2_wp
       ! ligphi_ens(12) = 7.9432823472e-3_wp
       ! ligphi_ens(13) = 6.3095734448e-3_wp
       ! ligphi_ens(14) = 5.0118723363e-3_wp
       ! ligphi_ens(15) = 3.9810717055e-3_wp
       ! ligphi_ens(16) = 3.1622776602e-3_wp
       ! ligphi_ens(17) = 2.5118864315e-3_wp
       ! ligphi_ens(18) = 1.995262315e-3_wp
       ! ligphi_ens(19) = 1.5848931925e-3_wp
       ! ligphi_ens(20) = 1.2589254118e-3_wp
       ! ligphi_ens(21) = 1.0e-3_wp
       ! ligphi_ens(22) = 7.9432823472e-4_wp
       ! ligphi_ens(23) = 6.3095734448e-4_wp
       ! ligphi_ens(24) = 5.0118723363e-4_wp
       ! ligphi_ens(25) = 3.9810717055e-4_wp
       ! ligphi_ens(26) = 3.1622776602e-4_wp
       ! ligphi_ens(27) = 2.5118864315e-4_wp
       ! ligphi_ens(28) = 1.995262315e-4_wp
       ! ligphi_ens(29) = 1.5848931925e-4_wp
       ! ligphi_ens(30) = 1.2589254118e-4_wp
       ! ligphi_ens(31) = 1.0e-4_wp
       ! ligphi_ens(32) = 7.9432823472e-5_wp
       ! ligphi_ens(33) = 6.3095734448e-5_wp
       ! ligphi_ens(34) = 5.0118723363e-5_wp
       ! ligphi_ens(35) = 3.9810717055e-5_wp
       ! ligphi_ens(36) = 3.1622776602e-5_wp
       ! ligphi_ens(37) = 2.5118864315e-5_wp
       ! ligphi_ens(38) = 1.995262315e-5_wp
       ! ligphi_ens(39) = 1.5848931925e-5_wp
       ! ligphi_ens(40) = 1.2589254118e-5_wp
       ! ligphi_ens(1) = 1.0e-5_wp
       ! ligphi_ens(2) = 7.9432823472e-6_wp
       ! ligphi_ens(3) = 6.3095734448e-6_wp
       ! ligphi_ens(4) = 5.0118723363e-6_wp
       ! ligphi_ens(5) = 3.9810717055e-6_wp
       ! ligphi_ens(6) = 3.1622776602e-6_wp
       ! ligphi_ens(7) = 2.5118864315e-6_wp
       ! ligphi_ens(8) = 1.995262315e-6_wp
       ! ligphi_ens(9) = 1.5848931925e-6_wp
       ! ligphi_ens(10) = 1.2589254118e-6_wp
       ! ligphi_ens(11) = 1.0e-6_wp
       ! ligphi_ens(12) = 7.9432823472e-7_wp
       ! ligphi_ens(13) = 6.3095734448e-7_wp
       ! ligphi_ens(14) = 5.0118723363e-7_wp
       ! ligphi_ens(15) = 3.9810717055e-7_wp
       ! ligphi_ens(16) = 3.1622776602e-7_wp
       ! ligphi_ens(17) = 2.5118864315e-7_wp
       ! ligphi_ens(18) = 1.995262315e-7_wp
       ! ligphi_ens(19) = 1.5848931925e-7_wp
       ! ligphi_ens(20) = 1.2589254118e-7_wp
       ! ligphi_ens(21) = 1.0e-7_wp



! Deep ocean box lifetime modifier, longer lifetime in the deep ocean
       lt_lifetime_factor_ens(1) = 1.0_wp
       ! lt_lifetime_factor_ens(1) = 1.0_wp
       ! lt_lifetime_factor_ens(2) = 5.0_wp
       ! lt_lifetime_factor_ens(1) = 10.0_wp
       ! lt_lifetime_factor_ens(2) = 20.0_wp
       ! lt_lifetime_factor_ens(3) = 30.0_wp
       ! lt_lifetime_factor_ens(4) = 50.0_wp
       ! lt_lifetime_factor_ens(5) = 75.0_wp
       ! lt_lifetime_factor_ens(1) = 100.0_wp
       ! lt_lifetime_factor_ens(7) = 150.0_wp
       ! lt_lifetime_factor_ens(8) = 200.0_wp
       ! lt_lifetime_factor_ens(9) = 300.0_wp
       ! lt_lifetime_factor_ens(10) = 500.0_wp

!===========================================================================
       ! SPECIFY FIRST FILE NUMBER FOR THE ENSEMBLE
       ! FILE NUMBER IS INCREMENTED BY 1 FOR EACH MEMBER
!===========================================================================
! File number identifier start for the ensemble
       id0            =     2040000

! Array inputs
#if defined(SIXBOX)
! For a 20SV AMOC, psi_in (i.e. Southern Ocean upwelling) needs to be 2x
       psi_in= 20.e6_wp

       dx    = [17.0e6_wp, 17.0e6_wp, 17.0e6_wp,                       &
                17.0e6_wp, 17.0e6_wp, 17.0e6_wp ]
       dy    = [ 4.0e6_wp,  4.0e6_wp,  2.0e6_wp,                       &
                 2.0e6_wp,  8.0e6_wp,  8.0e6_wp ]
       dz    = [ 100._wp, 3900._wp,  100._wp,                          & 
                3900._wp,  100._wp, 3900._wp ]
                
       depth = [       dz(1)/2._wp,                                    &
                 dz(1)+dz(2)/2._wp,                                    &
                       dz(3)/2._wp,                                    &
                 dz(3)+dz(4)/2._wp,                                    &
                       dz(5)/2._wp,                                    &
                 dz(5)+dz(6)/2._wp ]
                 
       m2deg    = 180._wp/(dy(1)+dy(3)+dy(5))  

       latitude = [((dy(1)/2._wp)+(dy(3)/2._wp)),                      &
                   ((dy(2)/2._wp)+(dy(4)/2._wp)),                      &
                   ((dy(3)/2._wp)              ),                      &
                   ((dy(4)/2._wp)              ),                      &
                   ((dy(5)/2._wp)+(dy(3)/2._wp)),                      &
                   ((dy(6)/2._wp)+(dy(4)/2._wp))                       &
                  ]
       latitude = -90._wp+(latitude*m2deg)

! define arrays (nbox*nbox long) of box connectivity for mixing and overturning (by rows)
! Box 1 mixes with box 2 and 3; 
! Box 2 mixes with box 1 and 4; 
! Box 3 mixes with box 1, 4 and 5; 
! Box 4 mixes with box 2, 3, and 6.
! Box 5 mixes with box 3 and 6.
! Box 4 mixes with box 4 and 5.
!                       Box1    Box2    Box3    Box4     Box5    Box6 
       Kin = RESHAPE([ 0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & ! Box1
                       1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, & ! Box2
                       1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, & ! Box3
                       0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, & ! Box4
                       0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, & ! Box5
                       0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp  & ! Box6
                     ], [ nbox, nbox ] )
       
! Box 1 is upstream of box 2; 
! Box 2 is upstream of box 1 and box 4; 
! Box 3 is upstream of box 1, 4, and 5; 
! Box 4 is upstream of box 2, 3, and 6;
! Box 5 is upstream of box 1 and box 6;
! Box 6 is upstream of box 4 and box 5.
!                       Box1     Box2     Box3     Box4     Box5     Box6  
       Pin = RESHAPE([ 0.00_wp, 1.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, & ! Box1
                       0.35_wp, 0.00_wp, 0.00_wp, 0.90_wp, 0.00_wp, 0.00_wp, & ! Box2
                       0.15_wp, 0.00_wp, 0.00_wp, 0.55_wp, 0.30_wp, 0.00_wp, & ! Box3
                       0.00_wp, 0.25_wp, 1.00_wp, 0.00_wp, 0.00_wp, 1.10_wp, & ! Box4
                       0.50_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.50_wp, & ! Box5
                       0.00_wp, 0.00_wp, 0.00_wp, 0.90_wp, 0.70_wp, 0.00_wp  & ! Box6
                     ], [ nbox, nbox ] )

! define array of remineralization coefficients (by rows)
! -1 indicates all of export is lost from cell, while 
! +1 indicates all of export is remineralized (gained) by cell
! Box 1 loses export from Box 1, which is completely remineralized in Box 2
! Box 3 loses export from Box 3, which is completely remineralized in Box 4 
! Box 5 loses export from Box 5, which is completely remineralized in Box 6 
!                       Box1    Box2    Box3    Box4     Box5    Box6 
       Rin = RESHAPE([-1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & ! Box1
                       0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & ! Box2
                       0.0_wp, 0.0_wp,-1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, & ! Box3
                       0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & ! Box4
                       0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,-1.0_wp, 1.0_wp, & ! Box5
                       0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp  & ! Box6
                     ], [ nbox, nbox ] )
                     
! Initial conditions
       thin(1:6)= [   20.0_wp,    4.0_wp,   -1.00_wp,   -1.00_wp,   20.0_wp,    4.0_wp ]
       sain(1:6)= [   35.5_wp,   35.5_wp,   34.75_wp,   34.75_wp,   35.0_wp,   35.0_wp ]
       cain(1:6)= [ 2100.0_wp, 2400.0_wp, 2100.00_wp, 2400.00_wp, 2100.0_wp, 2400.0_wp ]
       alin(1:6)= [ 2350.0_wp, 2400.0_wp, 2300.00_wp, 2400.00_wp, 2300.0_wp, 2400.0_wp ]
       phin(1:6)= [    0.0_wp,    2.5_wp,    2.50_wp,    2.50_wp,    0.0_wp,    2.5_wp ]
       niin(1:6)= [    0.0_wp,   36.0_wp,   30.00_wp,   36.00_wp,    0.0_wp,   36.0_wp ]

! Wind speed (m/s)for CO2 gas fluxes
       wind_in(1:6) = [ 5._wp, 0._wp, 10._wp, 0._wp, 5._wp, 0._wp ]
       
! Open surface fraction in contact with atmoshpere 
!  1 => fully open, <1 => flux impeded (e.g. by sea ice)
       fopen_in(1:6)= [ 1._wp, 0._wp, 1._wp, 0._wp, 1._wp, 0._wp ]

! Export ratio is smaller than 1 for the surface boxes
! Export ratio is 1 for the deep boxes
       eratio_in(1:6)= [ 0.1_wp, 0.5_wp, 0.1_wp, 0.5_wp, 0.1_wp, 0.5_wp ]
! ==============================================================================
! THESE SURFACE VALUES NEED TO BE ADJUSTED IF THE 6 BOX MODEL IS USED !!!
! ==============================================================================       
       
! Dust deposition in g Fe m-2 year-1
! Hydrothermal vent input of 1 Gmol/yr (Tagliabue et al., 2010)
! mol Fe/yr * g/mol * 1/area  == g Fe m-2 year-1....
!divide by 2.5e-3 because fe_sol is multiplied again within model.
       fe_input(1:6)= [ 1.5e-1_wp, (1.e9_wp*56._wp)/(2.5e-3_wp*dx(2)*dy(2)),   &
                        1.5e-3_wp, (1.e9_wp*56._wp)/(2.5e-3_wp*dx(4)*dy(4)),   &
                        7.5e-2_wp, (1.e9_wp*56._wp)/(2.5e-3_wp*dx(6)*dy(6)) ]

! Deep ocean box lifetime modifier to capture the gradient due to
! photodegradation near the surface and slower loss in the deep
! Modification: first order loss in the deep is set to 0
       dldz_in(1:6)  = [ 1._wp, 0._wp, 1._wp, 0._wp, 1._wp, 0._wp ]
       
#elif defined(FOURBOX)
! For a 20SV AMOC, psi_in (i.e. Southern Ocean upwelling) needs to be 2x
       psi_in= 40.e6_wp

       dx    = [17.0e6_wp, 17.0e6_wp, 17.0e6_wp, 17.0e6_wp]
       dy    = [ 1.0e6_wp,  3.0e6_wp, 12.0e6_wp, 16.0e6_wp]
       dz    = [50._wp,    50._wp,    50._wp,     5050._wp]
       depth = [       dz(1)/2._wp,                                    &
                       dz(2)/2._wp,                                    &
                       dz(3)/2._wp,                                    &
                 dz(1)+dz(4)/2._wp]
                 
       m2deg    = 180._wp/dy(4)  

       latitude = [(           +(dy(1)/2._wp)),                        &
                   (dy(1)      +(dy(2)/2._wp)),                        &
                   (dy(1)+dy(2)+(dy(3)/2._wp)),                        &
                   (           +(dy(4)/2._wp))                         &
                  ]
       latitude = -90._wp+(latitude*m2deg)

! define arrays (nbox*nbox long) of box connectivity for mixing and overturning (by rows)
! Box 1 mixes with box 2 and 4; 
! Box 2 mixes with box 1, 3 and 4; 
! Box 3 mixes with box 2 and 4; 
! Box 4 mixes with box 1, 2, and 3.
!                       Box1    Box2    Box3    Box4 
       Kin = RESHAPE([ 0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp,                 & ! Box1
                       1.0_wp, 0.0_wp, 1.0_wp, 1.0_wp,                 & ! Box2
                       0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp,                 & ! Box3
                       1.0_wp, 1.0_wp, 1.0_wp, 0.0_wp                  & ! Box4
                     ], [ nbox, nbox ] )
       
! Box 1 is upstream of box 4 (ie AABW downwelling); 
! Box 2 is upstream of box 1 and 3 (ie Antarctic divergence); 
! Box 3 is upstream of box 4 (ie NADW downwelling); 
! Box 4 is upstream of box 2 (ie SO upwelling).
!                       Box1    Box2    Box3    Box4 
       Pin = RESHAPE([ 0.0_wp, 0.0_wp, 0.0_wp, 0.5_wp,                 & ! Box1
                       0.5_wp, 0.0_wp, 0.5_wp, 0.0_wp,                 & ! Box2
                       0.0_wp, 0.0_wp, 0.0_wp, 0.5_wp,                 & ! Box3
                       0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp                  & ! Box4
                     ], [ nbox, nbox ] )

! define array of remineralization coefficients (by rows)
! -1 indicates all of export is lost from cell, while 
! +1 indicates all of export is remineralized (gained) by cell
! Box 1 loses export from Box 1, which is completely remineralized in Box 4 (Box 2 is adjacent)
! Box 2 loses export from Box 2, which is also completely remineralized in Box 4 (Box 1 is adjacent)
!                       Box1    Box2    Box3    Box4 
       Rin = RESHAPE([-1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp,                 & ! Box1
                       0.0_wp,-1.0_wp, 0.0_wp, 1.0_wp,                 & ! Box2
                       0.0_wp, 0.0_wp,-1.0_wp, 1.0_wp,                 & ! Box3
                       0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp                  & ! Box4
                     ], [ nbox, nbox ] )
                     
! Initial conditions
       thin(1:4)= [   -1._wp,    2._wp,   20._wp  ,    4._wp   ]
       sain(1:4)= [   35._wp,   34._wp,   35.50_wp,   34.75_wp ]
       cain(1:4)= [ 2100._wp, 2100._wp, 2100._wp  , 2400._wp   ]
       alin(1:4)= [ 2350._wp, 2300._wp, 2300._wp  , 2400._wp   ]
       phin(1:4)= [    2._wp,    2._wp,    0.0_wp ,    2.5_wp  ]
       niin(1:4)= [   32._wp,   32._wp,    0._wp  ,   36._wp   ]

! Wind speed (m/s)for CO2 gas fluxes
       wind_in(1:4) = [ 10._wp, 10._wp, 5._wp, 0._wp ]
       
! Open surface fraction in contact with atmoshpere 
!  1 => fully open, <1 => flux impeded (e.g. by sea ice)
       fopen_in(1:4)= [ 0.5_wp, 1._wp, 1._wp, 0._wp ]

! Export ratio is smaller than 1 for the surface boxes
! Export ratio is 1 for the deep boxes
       eratio_in(1:4)= [ 0.25_wp, 0.5_wp, 0.5_wp, 0.1_wp ]       
       
! Dust deposition in g Fe m-2 year-1
       fe_input(1:4)= [ 1.5e-3_wp, 1.5e-3_wp, 1.5e-1_wp,               &
! Hydrothermal vent input of 1 Gmol/yr (Tagliabue et al., 2010)
! mol Fe/yr * g/mol * 1/area  == g Fe m-2 year-1....
!divide by 2.5e-3 because fe_sol is multiplied again within model.
        (1.e9_wp*56._wp)/(2.5e-3_wp*dx(4)*dy(4)) ]

! Deep ocean box lifetime modifier to capture the gradient due to
! photodegradation near the surface and slower loss in the deep
! Modification: first order loss in the deep is set to 0
       dldz_in(1:4)  = [ 1._wp, 1._wp, 1._wp, 0._wp ]
#else

! Default to the three box model
       psi_in = 20.e6_wp
       
       dx    = [ 17.0e6_wp, 17.0e6_wp,   17.0e6_wp ]
       dy    = [  4.0e6_wp, 12.0e6_wp,   16.0e6_wp ]  
       dz    = [ 50.0_wp  , 50.0_wp  , 5050.0_wp   ]

       depth = [ dz(1)/2._wp, dz(2)/2._wp, dz(1)+dz(3)/2._wp]

       m2deg    = 180._wp/dy(3)  
       latitude = [(dy(1)       /2._wp) ,                              &
                   (dy(1)+(dy(2)/2._wp)),                              &
                   (dy(3)       /2._wp)                                &
                   ]
       latitude = -90._wp+(latitude*m2deg)

! Mixing mask array
       Kin = RESHAPE( [ 0._wp, 1._wp, 1._wp,                           &
                        1._wp, 0._wp, 1._wp,                           &
                        1._wp, 1._wp, 0._wp ],                         &
                      [ nbox, nbox ] )
       
! Overturning mask array
       Pin = RESHAPE([ 0._wp, 1._wp, 0._wp,                            &
                       0._wp, 0._wp, 1._wp,                            &
                       1._wp, 0._wp, 0._wp ],                          &
                        [ nbox, nbox ] )
       
! Remineralization coefficients 
! -1 indicates all of export is lost from cell, while 
! +1 indicates all of export is remineralized (gained) by cell
!  the first box (column one) loses export from box 1,
!  the second box (col two) loses export from box 2,
!  and the third box (col three) gains export from boxes 1 and 2 
       Rin = RESHAPE([ -1._wp, 0._wp, 1._wp,                           &
                        0._wp,-1._wp, 1._wp,                           &
                        0._wp, 0._wp, 0._wp ],                         &
                      [ nbox, nbox ] )

! Initial conditions
       thin(1:3)= [    2._wp,   20._wp  ,    4._wp   ]
       sain(1:3)= [   34._wp,   35.50_wp,   34.75_wp ]
       cain(1:3)= [ 2200._wp, 2100._wp  , 2400._wp   ]
       alin(1:3)= [ 2350._wp, 2350._wp  , 2400._wp   ]
       ! phin(1:3)= [    2._wp,    0._wp  ,    2.5_wp  ]
       ! niin(1:3)= [   25._wp,    0._wp  ,   35._wp   ]
       phin(1:3)= [    2.5_wp,    2.5_wp  ,    2.5_wp  ]
       niin(1:3)= [   35._wp,   35._wp  ,   35._wp   ]

! Initial concentrationesix in (u/n)mol/kg
! Here are some equilibrated values run for 100,000 yrs (round-off error notwithstanding)
! Make sure to compile without -DFIXATMPCO2
!       cain(1:3)= [ 2263.27105_wp, 2104.02729_wp, 2358.11830_wp ]
!       alin(1:3)= [ 2396.24755_wp, 2388.24068_wp, 2399.60156_wp ]
!       phin(1:3)= [ 1.85304_wp   , 0.31325_wp   , 2.49804_wp    ]
!       niin(1:3)= [ 24.68043_wp  , 0.04392_wp   , 35.00046_wp   ]
!       fein(1:3)= [ 0.01007_wp   , 0.37382_wp   , 0.55782_wp    ]
!       ltin(1:3)= [ 1.62217_wp   , 1.58451_wp   , 1.57992_wp    ]

! Wind speed (m/s)for CO2 gas fluxes
       wind_in(1:3) = [ 10._wp, 5._wp, 0._wp ]
       
! Open surface fraction in contact with atmoshpere 
!  1 => fully open, <1 => flux impeded (e.g. by sea ice)
       fopen_in(1:3)= [ 1._wp, 1._wp, 0._wp ]

! Export ratio is smaller than 1 for the surface boxes
! Export ratio is 1 for the deep boxes
       eratio_in(1:3)= [ 0.25_wp, 0.1_wp, 0.5_wp ]   


! Dust deposition in g Fe m-2 year-1
       fe_input(1:3)= [ 1.5e-3_wp, 1.5e-1_wp,                          &
! Hydrothermal vent input of 1 Gmol/yr (Tagliabue et al., 2010)
! mol Fe/yr * g/mol * 1/area  == g Fe m-2 year-1....
!divide by 2.5e-3 because fe_sol is multiplied again within model.
        (1.e9_wp*56._wp)/(2.5e-3_wp*dx(3)*dy(3)) ]
! Deep ocean box lifetime modifier to capture the gradient due to
! photodegradation near the surface and slower loss in the deep
! Modification: first order loss in the deep is set to 0
       dldz_in(1:3)  = [ 1._wp, 1._wp, 0.0_wp ]
#endif

!=======================================================================
! Initialize csv file
!=======================================================================
       CALL WRITE_CSV_HEADER(filename_csv, headers, SIZE(headers))

       DEALLOCATE(headers)
       
!=======================================================================
! Start loop over ensemble
!=======================================================================

id = id0

! Each ensemble member is defined by a unique combination of parameters
! Looping over all combinations of parameters that are specfied previously
DO irFeC_pb = 1, size(rFeC_pb_ens)
    rFeC_pb = rFeC_pb_ens(irFeC_pb)

    DO imu0 = 1, size(mu0_ens)
        mu0 = mu0_ens(imu0)

        DO im_l = 1, size(m_l_ens)
            m_l = m_l_ens(im_l)

            DO im_q = 1, size(m_q_ens)
                m_q = m_q_ens(im_q)

                DO ikappa = 1, size(kappa_ens)
                    kappa = kappa_ens(ikappa)

                    DO ikfe_p = 1, size(kfe_p_ens)
                        kfe_p = kfe_p_ens(ikfe_p)

                        DO ikldoc_p = 1, size(kldoc_p_ens)
                            kldoc_p = kldoc_p_ens(ikldoc_p)

                            DO ipge_deep = 1, size(pge_deep_ens)
                                pge_deep = pge_deep_ens(ipge_deep)

                                DO iphi = 1, size(phi_ens)
                                    phi = phi_ens(iphi)

                                    DO irCLig = 1, size(rCLig_ens)
                                        rCLig = rCLig_ens(irCLig)

                                        DO ilt_lifein = 1, size(lt_lifein_ens)
                                            lt_lifein = lt_lifein_ens(ilt_lifein)

                                          DO iligphi = 1, size(ligphi_ens)
                                            ligphi_in = ligphi_ens(iligphi)

                                            DO ilt_lifetime_factor = 1, size(lt_lifetime_factor_ens)
                                                lt_lifetime_factor = lt_lifetime_factor_ens(ilt_lifetime_factor)



                                                 WRITE(*,*) 'id = ', id

                                                        dldz_in(1:3)  = [ 1._wp, 1._wp, 1.0_wp / lt_lifetime_factor ]
                                                        !=================================================================================
                                                        ! Prokaryotic growth efficiency
                                                        pge_in(1:3)= [ 0.25_wp, 0.25_wp, 0.25_wp ] ! setting pge for the different boxes
                                                        !=================================================================================

                                                        !=================================================================================
                                                        ! part of the deep box ligands are refractory, assumed %
                                                        ! division by 2: 50% of the ligands are labile
                                                        ! division by 5: only 20% of the ligands are labile
                                                        ! division by 10: only 10% of the ligands are labile
                                                        ! ... and so on
                                                        rCLig_in(1:3) = [ rCLig, rCLig, rCLig / 5.0_wp ]
                                                        !=================================================================================

                                                        ! Initial conditions (random selection if desired)
                                                        ! call random_number(ltin) ! select random number between 0 and 1
                                                        ! ltin = ltin * 100.0_wp   ! scale to 0-100 nmol/kg
                                                        ! call random_number(fein) ! select random number between 0 and 1
                                                        ! fein = fein * 100.0_wp   ! scale to 0-100 nmol/kg

                                                 call model(                                                    &
                                                        id,                                                        &
                                                        filename_csv,                                              &
                                                        maxyears,                                                  &
                                                        outputyears,                                               &
                                                        outstepmax,                                                &
                                                        dx,                                                        &
                                                        dy,                                                        &
                                                        dz,                                                        &
                                                        depth,                                                     &
                                                        latitude,                                                  &
                                                        Kin,                                                       &
                                                        Rin,                                                       &
                                                        Pin,                                                       &
                                                        psi_in,                                                    &
                                                        dif_in,                                                    &    
                                                        alpha_yr,                                                  &
                                                        ligphi_in,                                                 &
                                                        lt_lifein,                                                 &
                                                        dldz_in,                                                   &
                                                        fe_input,                                                  &
                                                        wind_in,                                                   &
                                                        fopen_in,                                                  &
                                                        rCLig_in,                                                  &
                                                        thin,                                                      &
                                                        sain,                                                      &
                                                        cain,                                                      &
                                                        alin,                                                      &
                                                        phin,                                                      &
                                                        niin,                                                      &
                                                        fein,                                                      &
                                                        ltin,                                                      &
                                                        atpco2in,                                                  &
                                                        eratio_in,                                                 &
                                                        pge_in,                                                    &
                                                        pbin,                                                      &
                                                        ldocin,                                                    &
                                                        tout,                                                      &            
                                                        thout,                                                     &
                                                        sout,                                                      &
                                                        cout,                                                      &
                                                        aout,                                                      &
                                                        pout,                                                      &
                                                        nout,                                                      &
                                                        fout,                                                      &
                                                        lout,                                                      &
                                                        expout,                                                    &
                                                        nlout,                                                     &
                                                        psout,                                                     &
                                                        ocpco2out,                                                 &
                                                        atpco2out,                                                 &
                                                        pbout,                                                     &
                                                        ldocout                                                    &
                                                 )

                                                 id = id + 1

                                                 ! Closing all the loops that were opened
                                              END DO ! ilt_lifetime_factor
                                          END DO ! ligphi_ens
                                        END DO ! lt_lifein
                                    END DO ! rCLig_ens
                                END DO ! phi_ens
                            END DO ! pge_deep_ens
                        END DO ! kldoc_p_ens
                    END DO ! kfe_p_ens
                END DO ! ikappa
            END DO ! m_q_ens
        END DO ! m_l_ens
    END DO ! mu0_ens
END DO ! rFeC_pb_ens


!=======================================================================
       END PROGRAM MICROCOSM_MODEL_ENSEMBLE
!=======================================================================