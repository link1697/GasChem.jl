# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #MODULE: fullchem_mod.F90
#
# #DESCRIPTION: Contines arrays and routines for the GEOS_Chem "fullchem"
#  mechanism, which is implemented in KPP - generated Fortran code.
#\\
#\\
# #INTERFACE:
#
module fullchem_mod
#
# #USES:
#
  private
#
# #PUBLIC MEMBER FUNCTIONS:
#
  public  :: do_fullchem
  public  :: init_fullchem
  public  :: cleanup_fullchem
#
# #REVISION HISTORY:
#  14 Dec 2015 - M.S. Long - Initial version
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #PRIVATE TYPES:
#
  # Species ID flags (and logicals to denote if species are present)
#ifdef model_geos
#endif

  # Diagnostic flags
#ifdef model_geos
#endif

  # SAVEd scalars
  integer,  save        :: prevday   = - 1
  integer,  save        :: prevmonth = - 1

  # Arrays
  integer,  allocatable :: pl_kpp_id (:      )
  real(f4), allocatable :: jvcountday(:, :, :  )
  real(f4), allocatable :: jvcountmon(:, :, :  )
  real(f4), allocatable :: jvsumday  (:, :, :, :)
  real(f4), allocatable :: jvsummon  (:, :, :, :)

contains
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #ROUTINE: do_fullchem
#
# #DESCRIPTION: Driver subroutine for the KPP fullchem mechanism.
#\\
#\\
# #INTERFACE:
#
  subroutine do_fullchem( input_opt,  state_chm, state_diag, state_grid, state_met, rc                         )
#
# #USES:
#
#ifdef tomas
#ifdef bpch_diag
#endif
#endif
#
# #INPUT PARAMETERS:
#
    type(optinput), intent(in)    :: input_opt  # Input Options object
    type(grdstate), intent(in)    :: state_grid # Grid State object
#
# #INPUT / OUTPUT PARAMETERS:
#
    type(metstate), intent(inout) :: state_met  # Meteorology State object
    type(chmstate), intent(inout) :: state_chm  # Chemistry State object
    type(dgnstate), intent(inout) :: state_diag # Diagnostics State object
#
# #OUTPUT PARAMETERS:
#
    integer,        intent(out)   :: rc         # Success or failure
#
# #REVISION HISTORY:
#  14 Dec 2015 - M.S. Long - Initial version
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #LOCAL VARIABLES:
#
    # Scalars
    real(fp)               :: so4_frac,   t,         tin
    real(fp)               :: tout,       sr,        lwc

    # Strings
    character(len = 63)      :: origunit
    character(len = 255)     :: errmsg,   thisloc

    # SAVEd scalars
    logical,  save         :: firstchem = true
    integer,  save         :: ch4_year  = - 1

    # For

#ifdef model_classic
#ifndef no_omp
    integer, external      :: omp_get_thread_num
#endif
#endif

    # Arrays
    real(dp)               :: rcntrl (20)
    real(dp)               :: rstate (20)
    real(fp)               :: before(state_grid%nx, state_grid%ny, state_grid%nz, state_chm%nadvect       )

    # For tagged CO saving
    real(fp)               :: lch4, pco_tot, pco_ch4, pco_nmvoc

    # Objects
    type(species), pointer :: spcinfo

    # For testing purposes

    # OH reactivity and KPP reaction rate diagnostics
    real(fp)               :: ohreact
    real(dp)               :: vloc(nvar),     aout(nreact)
#ifdef model_geos
    real(f4)               :: noxtau, noxconc, nox_weight, nox_tau_weighted
    real(f4)               :: trop_nox_tau
    real(f4)               :: tropv_nox_tau(state_grid%nx, state_grid%ny)
    real(f4)               :: tropv_nox_mass(state_grid%nx, state_grid%ny)
    real(dp)               :: vdotout(nvar), localc(nspec)
#endif
#ifdef model_wrf
    real(dp)               :: localc(nspec)
#endif

    # Objects
    type(dgnmap), pointer :: mapdata = > null()

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # do_fullchem begins here#
    # NOTE: FlexChem timer is started in DO_CHEMISTRY (the calling routine)
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # Initialize
    rc        =  gc_success
    errmsg    =  ""
    thisloc   =  " - > at do_fullchem (in module geoscore / fullchem_mod.f90)"
    spcinfo   = > null()
    prtdebug  =  ( input_opt%lprt && input_opt%amiroot )
    day       =  get_day()    # Current day
    month     =  get_month()  # Current month
    year      =  get_year()   # Current year
    thread    =  1
    failed2x  = false

    # Set a switch that allows you to toggle off photolysis for testing
    # (default value : TRUE)
    do_photchem = true

    # Print information the first time that DO_FULLCHEM is called
    printfirsttimeinfo( input_opt, state_chm, firstchem, do_photchem )

    # Zero diagnostic archival arrays to make sure that we don"t have any
    # leftover values from the last timestep near the top of the chemgrid
    if (state_diag%archive_loss           ) state_diag%loss           = 0.0_f4 end
    if (state_diag%archive_prod           ) state_diag%prod           = 0.0_f4 end
    if (state_diag%archive_jval           ) state_diag%jval           = 0.0_f4 end
    if (state_diag%archive_jvalo3o1d      ) state_diag%jvalo3o1d      = 0.0_f4 end
    if (state_diag%archive_jvalo3o3p      ) state_diag%jvalo3o3p      = 0.0_f4 end
    if (state_diag%archive_jnoon          ) state_diag%jnoon          = 0.0_f4 end
    if (state_diag%archive_prodcofromch4  ) state_diag%prodcofromch4  = 0.0_f4 end
    if (state_diag%archive_prodcofromnmvoc) state_diag%prodcofromnmvoc = 0.0_f4 end
    if (state_diag%archive_ohreactivity   ) state_diag%ohreactivity   = 0.0_f4 end
    if (state_diag%archive_rxnrate        ) state_diag%rxnrate        = 0.0_f4 end
    if (state_diag%archive_kppdiags) 
       if (state_diag%archive_kppintcounts) state_diag%kppintcounts   = 0.0_f4 end
       if (state_diag%archive_kppjaccounts) state_diag%kppjaccounts   = 0.0_f4 end
       if (state_diag%archive_kpptotsteps ) state_diag%kpptotsteps    = 0.0_f4 end
       if (state_diag%archive_kppaccsteps ) state_diag%kppaccsteps    = 0.0_f4 end
       if (state_diag%archive_kpprejsteps ) state_diag%kpprejsteps    = 0.0_f4 end
       if (state_diag%archive_kppludecomps) state_diag%kppludecomps   = 0.0_f4 end
       if (state_diag%archive_kppsubsts   ) state_diag%kppsubsts      = 0.0_f4 end
       if (state_diag%archive_kppsmdecomps) state_diag%kppsmdecomps   = 0.0_f4 end
    end

    # Keep track of the boxes where it is local noon in the JNoonFrac
    # diagnostic. When time - averaged, this will be the fraction of time
    # that local noon occurred at a grid box. (bmy, 4 / 2 / 19)
    if ( state_diag%archive_jnoonfrac ) 
       where( state_met%islocalnoon )
          state_diag%jnoonfrac = 1.0_f4
       elsewhere
          state_diag%jnoonfrac = 0.0_f4
       end
    end

#if defined( model_geos )
    if ( state_diag%archive_noxtau     ) state_diag%noxtau(:, :, :) = 0.0_f4 end
    if ( state_diag%archive_tropnoxtau ) 
       state_diag%tropnoxtau(:, :) = 0.0_f4
       tropv_nox_mass(:, :) = 0.0_f4
       tropv_nox_tau(:, :)  = 0.0_f4
    end
#endif # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = end
    # Zero out certain species:
    # - isoprene oxidation counter species (dkh, bmy, 6 / 1 / 06)
    # - isoprene - NO3 oxidation counter species (hotp, 6 / 1 / 10)
    # - if SOA or SOA_SVPOA, aromatic oxidation counter species
    #      (dkh, 10 / 06 / 06)
    # - if SOA_SVPOA, LNRO2H and LNRO2N for NAP (hotp 6 / 25 / 09
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    for n = 1: state_chm%nspecies

       # Get info about this species from the species database
       spcinfo = > state_chm%spcdata(n)%info

       # isoprene oxidation counter species
       if ( trim( spcinfo%name ) = = "lisopoh" || trim( spcinfo%name ) = = "lisopno3" )  end
          state_chm%species(n)%conc(:, :, :) = 0.0_fp
       end

       # aromatic oxidation counter species
       if ( input_opt%lsoa || input_opt%lsvpoa ) 
          select case ( trim( spcinfo%name ) )
             case ( "lbro2h", "lbro2n", "ltro2h", "ltro2n", "lxro2h", "lxro2n", "lnro2h", "lnro2n" )
                state_chm%species(n)%conc(:, :, :) = 0.0_fp
          end
       end

       # Temporary fix for CO2
       # CO2 is a dead species and needs to be set to zero to
       # match the old SMVGEAR code (mps, 6 / 14 / 16)
       if ( trim( spcinfo%name ) = = "co2" )  end
          state_chm%species(n)%conc(:, :, :) = 0.0_fp
       end

       # Free pointer
       spcinfo = > null()

    end

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Convert species to [molec / cm3] (ewl, 8 / 16 / 16)
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    convert_spc_units( input_opt,            state_chm,   state_grid, state_met,           "molec / cm3", rc, origunit = origunit                               )
    if ( rc / = gc_success )  end
       errmsg = "unit conversion error#"
       gc_error( errmsg, rc, "fullchem_mod.f90")
       return
    end

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Call photolysis routine to compute J - Values
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    if ( input_opt%usetimers ) 
       timer_end  ( " = > flexchem",           rc )
       timer_start( " = > fast - jx photolysis", rc )
    end

    # Do Photolysis
    wavelength = 0
    fast_jx( wavelength, input_opt,  state_chm, state_diag, state_grid, state_met, rc                     )

    # Trap potential errors
    if ( rc / = gc_success )  end
       errmsg = "error encountered in "fast_jx"#"
       gc_error( errmsg, rc, thisloc )
       return
    end

    #### Debug
    if ( prtdebug ) 
       debug_msg( "### do_fullchem: after fast_jx" )
    end

    if ( input_opt%usetimers ) 
       timer_end  ( " = > fast - jx photolysis", rc )
       timer_start( " = > flexchem",           rc ) # end
    end

#if defined( model_geos ) || defined( model_wrf )
    # Init diagnostics
    if ( associated(state_diag%kpperror) ) 
       state_diag%kpperror(:, :, :) = 0.0
    end
#endif # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = end
    # Archive concentrations before chemistry
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    if ( state_diag%archive_concbeforechem ) 
       # Point to mapping obj specific to SpeciesConc diagnostic collection
       mapdata = > state_diag%map_concbeforechem

       #$OMP PARALLEL DO#$OMP DEFAULT( SHARED )#$OMP PRIVATE( N, S   )
       for s = 1: mapdata%nslots
          n = mapdata%slot2id(s)
          state_diag%concbeforechem(:, :, :, s) = state_chm%species(n)%conc(:, :, :)
       end
       #$OMP END PARALLEL DO

       # Free pointer
       mapdata = > null()
    end

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Set up integration convergence conditions and timesteps
    # (cf. M. J. Evans)
    #
    # NOTE: ATOL and RTOL are defined in gckpp_Global.F90 so they
    # are probably only used as INTENT(IN).  Therefore, it is
    # probably safe to define them here outside the OpenMP loop.
    # (bmy, 3 / 28 / 16)
    #
    # The ICNTRL vector specifies options for the solver.  We can
    # define ICNTRL outside of the "SOLVE CHEMISTRY" parallel loop
    # below because it is passed to KPP as INTENT(IN).
    # (bmy, 3 / 28 / 16)
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    #%%%%% TIMESTEPS %%%%%

    # mje Set up conditions for the integration
    # mje chemical timestep and convert it to seconds.
    dt        = get_ts_chem() # [s]
    t         = 0d0
    tin       = t
    tout      = t + dt

    #%%%%% CONVERGENCE CRITERIA %%%%%

    # Absolute tolerance
    atol      = 1e-2_dp

    # Relative tolerance
    rtol      = 0.5e-2_dp

    #%%%%% SOLVER OPTIONS %%%%%

    # Zero all slots of ICNTRL
    icntrl    = 0

    # 0 - non - autonomous, 1 - autonomous
    icntrl(1) = 1

    # 0 - vector tolerances, 1 - scalars
    icntrl(2) = 0

    # Select Integrator
    # ICNTRL(3) - > selection of a particular method.
    # For Rosenbrock, options are:
    # = 0 :  default method is Rodas3
    # = 1 :  method is  Ros2
    # = 2 :  method is  Ros3
    # = 3 :  method is  Ros4
    # = 4 :  method is  Rodas3
    # = 5:   method is  Rodas4
    icntrl(3) = 4

    # 0 - adjoint, 1 - no adjoint
    icntrl(7) = 1

    # Turn off calling Update_SUN, Update_RCONST, Update_PHOTO from within
    # the integrator.  Rate updates are done before calling KPP.
    icntrl(15) = - 1

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # %%%%% SOLVE CHEMISTRY -  - This is the main KPP solver loop %%%%%
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
100 format("no. of function calls:", i6, / , "no. of jacobian calls:", i6, / , "no. of steps:         ", i6, / , "no. of accepted steps:", i6, / , "no. of rejected steps ", i6, / , "       (except at very beginning)",          / , "no. of lu decompositions:             ", i6, / , "no. of forward / backward substitutions:", i6, / , "no. of singular matrix decompositions:", i6, / , / , "texit, the time corresponding to the      ",        / , "       computed y upon return:            ", f11.4, / , "hexit, last accepted step before exit:    ", f11.4, / , "hnew, last predicted step (not yet taken):", f11.4 )

    if ( input_opt%usetimers ) 
       timer_start( " - > flexchem loop", rc )
    end

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # MAIN LOOP: Compute reaction rates and chemical solver
    #
    # Variables not listed here are held THREADPRIVATE in gckpp_Global.F90
    # #$omp collapse(3) vectorizes the loop and #$OMP DYNAMIC(24) sends
    # 24 boxes at a time to each core...  when that core is finished,
    # it gets a nother chunk of 24 boxes.  This should lead to better
    # load balancing, and will spread the sunrise / sunset boxes across
    # more cores.
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #$OMP PARALLEL DO#$OMP DEFAULT( SHARED                                                   )#$OMP PRIVATE( I,        J,        L,       N                           )#$OMP PRIVATE( SO4_FRAC, IERR,     RCNTRL,  ISTATUS,   RSTATE           )#$OMP PRIVATE( SpcID,    KppID,    F,       P,         Vloc             )#$OMP PRIVATE( Aout,     Thread,   RC,      S,         LCH4             )#$OMP PRIVATE( OHreact,  PCO_TOT,  PCO_CH4, PCO_NMVOC, SR               )#$OMP PRIVATE( SIZE_RES, LWC                                            )#ifdef model_geos
    #$OMP PRIVATE( NOxTau,     NOxConc,  Vdotout, localC                    )#$OMP PRIVATE( NOx_weight, NOx_tau_weighted                             )#endif
#ifdef model_wrf
    #$OMP PRIVATE( localC                                                   )#endif
    #$OMP COLLAPSE( 3                                                       )#$OMP SCHEDULE( DYNAMIC, 24                                             )
    for l = 1: state_grid%nz
    for j = 1: state_grid%ny
    for i = 1: state_grid%nx

       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       # Initialize private loop variables for each (I, J, L)
       # Other private variables will be assigned in Set_Kpp_GridBox_Values
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       ierr      = 0                        # KPP success or failure flag
       istatus   = 0.0_dp                   # Rosenbrock output
       rcntrl    = 0.0_fp                   # Rosenbrock input
       rstate    = 0.0_dp                   # Rosenbrock output
       so4_frac  = 0.0_fp                   # Frac of SO4 avail for photolysis
       p         = 0                        # GEOS - Chem photolyis species ID
       lch4      = 0.0_fp                   # P / L diag: Methane loss rate
       pco_tot   = 0.0_fp                   # P / L diag: Total P(CO)
       pco_ch4   = 0.0_fp                   # P / L diag: P(CO) from CH4
       pco_nmvoc = 0.0_fp                   # P / L diag: P(CO) from NMVOC
       sr        = 0.0_fp                   # Enhancement to O2 catalysis rate
       lwc       = 0.0_fp                   # Liquid water content
       size_res  = false                  # Size resolved calculation?
       c         = 0.0_dp                   # KPP species conc"s
       rconst    = 0.0_dp                   # KPP rate constants
       photol    = 0.0_dp                   # Photolysis array for KPP
       k_cld     = 0.0_dp                   # Sulfur in - cloud rxn het rates
       k_mt      = 0.0_dp                   # Sulfur sea salt rxn het rates
       cfactor   = 1.0_dp                   # KPP conversion factor
#ifdef model_classic
#ifndef no_omp
       thread    = omp_get_thread_num() + 1 # OpenMP thread number
#endif
#endif
#if defined( model_geos ) || defined( model_wrf )
       localc    = 0.0_dp                   # Local backup array for C
#endif # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = end
       # Get photolysis rates (daytime only)
       #
       # NOTE: The ordering of the photolysis reactions here is
       # the order in the Fast - J definition file FJX_j2j.dat.
       # I"ve assumed that these are the same as in the text files
       # but this may have been changed.  This needs to be checked
       # through more thoroughly -  - M. Long (3 / 28 / 16)
       #
       # ALSO NOTE: We moved this section above the test to see if grid
       # box (I, J, L) is in the chemistry grid.  This will ensure that
       # J - value diagnostics are defined for all levels in the column.
       # This modification was validated by a geosfp_4x5_standard
       # difference test. (bmy, 1 / 18 / 18)
       #
       # Update SUNCOSmid threshold from 0 to cos(98 degrees) since
       # fast - jx allows for SZA down to 98 degrees. This is important in
       # the stratosphere - mesosphere where sunlight still illuminates at
       # high altitudes if the sun is below the horizon at the surface
       # (update submitted by E. Fleming (NASA), 10 / 11 / 2018)
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       if ( state_met%suncosmid(i, j) > - 0.1391731e + 0_fp ) 

          # Start timer
          if ( input_opt%usetimers ) 
             timer_start( timername = " - > photolysis rates", inloop    = true, threadnum = thread, rc        = rc                               )
          end

          # Get the fraction of H2SO4 that is available for photolysis
          # (this is only valid for UCX - enabled mechanisms)
          so4_frac = so4_photfrac( i, j, l )

          # Adjust certain photolysis rates:
          # (1) H2SO4 + hv - > SO2 + OH + OH   (UCX - based mechanisms)
          # (2) O3    + hv - > O2  + O         (UCX - based mechanisms)
          # (2) O3    + hv - > OH  + OH        (trop - only mechanisms)
          photrate_adj( input_opt, state_diag, state_met, i, j,         l,          so4_frac,  ierr         )

          # Loop over the FAST - JX photolysis species
          for n = 1: jvn_

             if ( do_photchem ) 

                # Copy photolysis rate from FAST_JX into KPP PHOTOL array
                photol(n) = zpj(l, n, i, j)

             end

             # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
             # HISTORY (aka netCDF diagnostics)
             #
             # Instantaneous photolysis rates [s - 1] (aka J - values)
             # and noontime photolysis rates [s - 1]
             #
             #    NOTE: Attach diagnostics here instead of in module
             #    fast_jx_mod.F90 so that we can get the adjusted photolysis
             #    rates (output from routne PHOTRATE_ADJ above).
             #
             # The mapping between the GEOS - Chem photolysis species and
             # the FAST - JX photolysis species is contained in the lookup
             # table in input file FJX_j2j.dat.

             # Some GEOS - Chem photolysis species may have multiple
             # branches for photolysis reactions.  These will be
             # represented by multiple entries in the FJX_j2j.dat
             # lookup table.
             #
             #    NOTE: For convenience, we have stored the GEOS - Chem
             #    photolysis species index (range: 1.0.State_Chm%nPhotol)
             #    for each of the FAST - JX photolysis species (range;
             #    1.0.JVN_) in the GC_PHOTO_ID array (located in module
             #    CMN_FJX_MOD.F90).
             # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

             # GC photolysis species index
             p = gc_photo_id(n)

             # If this FAST_JX photolysis species maps to a valid
             # GEOS - Chem photolysis species (for this simulation)...
             if ( p > 0 && p < = state_chm%nphotol )  end

                # Archive the instantaneous photolysis rate
                # (summing over all reaction branches)
                if ( state_diag%archive_jval ) 
                   s = state_diag%map_jval%id2slot(p)
                   if ( s > 0 ) 
                      state_diag%jval(i, j, l, s) = state_diag%jval(i, j, l, s) + photol(n)
                   end
                end

                # Archive the noontime photolysis rate
                # (summing over all reaction branches)
                if ( state_met%islocalnoon(i, j) ) 
                   if ( state_diag%archive_jnoon ) 
                      s = state_diag%map_jnoon%id2slot(p)
                      if ( s > 0 ) 
                         state_diag%jnoon(i, j, l, s) = state_diag%jnoon(i, j, l, s) + photol(n)
                      end
                   end
                end

             elseif ( p = = state_chm%nphotol + 1 )  end

                # J(O3_O1D).  This used to be stored as the nPhotol + 1st
                # diagnostic in Jval, but needed to be broken off
                # to facilitate cleaner diagnostic indexing (bmy, 6 / 3 / 20)
                if ( state_diag%archive_jvalo3o1d ) 
                   state_diag%jvalo3o1d(i, j, l) = state_diag%jvalo3o1d(i, j, l) + photol(n)
                end

             elseif ( p = = state_chm%nphotol + 2 )  end

                # J(O3_O3P).  This used to be stored as the nPhotol + 2nd
                # diagnostic in Jval, but needed to be broken off
                # to facilitate cleaner diagnostic indexing (bmy, 6 / 3 / 20)
                if ( state_diag%archive_jvalo3o3p ) 
                   state_diag%jvalo3o3p(i, j, l) = state_diag%jvalo3o3p(i, j, l) + photol(n)
                end

             end
          end

          # Stop timer
          if ( input_opt%usetimers ) 
             timer_end( timername = " - > photolysis rates", inloop    = true, threadnum = thread, rc        = rc                                 )
          end
       end

       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       # Test if we need to for the chemistry for box (I, J, L):
       # otherwise move onto the next box.
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

       # if we are not in the troposphere don"t do the chemistry#
       if ( ! state_met%inchemgrid(i, j, l) ) cycle

       # Skipping buffer zone (lzh, 08 / 10 / 2014)
       if ( state_grid%nestedgrid ) 
          if ( j < =                 state_grid%southbuffer ) cycle end
          if ( j >  state_grid%ny - state_grid%northbuffer ) cycle
          if ( i < =                 state_grid%eastbuffer  ) cycle end
          if ( i >  state_grid%nx - state_grid%westbuffer  ) cycle
       end

       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       # Initialize the KPP "C" vector of species concentrations [molec / cm3]
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       for n = 1: nspec
          spcid = state_chm%map_kppspc(n)
          c(n)  = 0.0_dp
          if ( spcid > 0 ) c(n) = state_chm%species(spcid)%conc(i, j, l) end
       end

       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       # Start KPP main timer
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       if ( input_opt%usetimers ) 
          timer_start( timername = " - > kpp", inloop    = true, threadnum = thread, rc        = rc                                  )
       end

       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       # CHEMISTRY MECHANISM INITIALIZATION (#1)
       #
       # Populate KPP global variables and arrays in gckpp_global.F90
       #
       # NOTE: This has to be done before Set_Sulfur_Chem_Rates, so that
       # the NUMDEN and SR_TEMP KPP variables will be populated first.
       # Otherwise this can lead to differences in output that are evident
       # when running with different numbers of OpenMP cores.
       # See https: / / github.com / geoschem / geos - chem / issues / 1157
       # -  - Bob Yantosca (08 Mar 2022)
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

       # Start timer
       if ( input_opt%usetimers ) 
          timer_start( timername = " - > init kpp", inloop    = true, threadnum = thread, rc        =  rc                                 )
       end

       # Copy values into the various KPP global variables
       set_kpp_gridbox_values( i          = i, j          = j, l          = l, input_opt  = input_opt, state_chm  = state_chm, state_grid = state_grid, state_met  = state_met, rc         = rc                         )

       # Stop timer
       if ( input_opt%usetimers ) 
          timer_end( timername =  " - > init kpp", inloop    = true, threadnum = thread, rc        =  rc                                   )
       end

       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       # CHEMISTRY MECHANISM INITIALIZATION (#2)
       #
       # Update reaction rates [1 / s] for sulfur chemistry in cloud and on
       # seasalt.  These will be passed to the KPP chemical solver.
       #
       # NOTE: This has to be done before fullchem_SetStateHet so that
       # State_Chm%HSO3_aq and State_Chm%SO3_aq will be populated first.
       # These are copied into State_Het%HSO3_aq and State_Het%SO3_aq.
       # See https: / / github.com / geoschem / geos - chem / issues / 1157
       # -  - Bob Yantosca (08 Mar 2022)
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

       # Start timer
       if ( input_opt%usetimers ) 
          timer_start( timername = "     rconst", inloop    = true, threadnum = thread, rc        = rc                                  )
       end

       # Compute sulfur chemistry reaction rates [1 / s]
       # If size_res = T, we"ll fullchem_HetDropChem below.
       set_sulfur_chem_rates( i          = i, j          = j, l          = l, input_opt  = input_opt, state_chm  = state_chm, state_diag = state_diag, state_grid = state_grid, state_met  = state_met, size_res   = size_res, rc         = rc                          )

       # Stop timer
       if ( input_opt%usetimers ) 
          timer_end( timername = "     rconst", inloop    = true, threadnum = thread, rc        = rc                                    )
       end

       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       # CHEMISTRY MECHANISM INITIALIZATION (#3)
       #
       # Populate the various fields of the State_Het object.
       #
       # NOTE: This has to be done after fullchem_SetStateHet so that
       # State_Chm%HSO3_aq and State_Chm%SO3_aq will be populated first.
       # These are copied into State_Het%HSO3_aq and State_Het%SO3_aq.
       # See https: / / github.com / geoschem / geos - chem / issues / 1157
       # -  - Bob Yantosca (08 Mar 2022)
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

       # Start timer
       if ( input_opt%usetimers ) 
          timer_start( timername = " - > init kpp", inloop    = true, threadnum = thread, rc        =  rc                                 )
       end

       # Populate fields of the State_Het object
       fullchem_setstatehet( i         = i, j         = j, l         = l, input_opt = input_opt, state_chm = state_chm, state_met = state_met, h         = state_het, rc        = rc                            )

       # Stop timer
       if ( input_opt%usetimers ) 
          timer_end( timername =  " - > init kpp", inloop    = true, threadnum = thread, rc        =  rc                                   )
       end

       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       # CHEMISTRY MECHANISM INITIALIZATION (#5)
       #
       # Call Het_Drop_Chem (formerly located in sulfate_mod.F90) to
       # estimate the in - cloud sulfate production rate in heterogeneous
       # cloud droplets based on the Yuen et al., 1996 parameterization.
       # Code by Becky Alexander (2011) with updates by Mike Long and Bob
       # Yantosca (2021).
       #
       # We will only Het_Drop_Chem if:
       # (1) It is at least 0.01% cloudy
       # (2) We are doing a size - resolved computation
       # (3) The grid box is over water
       # (4) The temperature is above - 5C
       # (5) Liquid water content is nonzero
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       if ( state_met%cldf(i, j, l) > 1.0e-4_fp ) 

          # Start timer
          if ( input_opt%usetimers ) 
             timer_start( timername = "     rconst", inloop    = true, threadnum = thread, rc        =  rc                              )
          end

          # Liquid water content (same formula from the old sulfate_mod.F90)
          lwc = ( state_met%ql(i, j, l) * state_met%airden(i, j, l) *   1.0e-3_fp           / state_met%cldf(i, j, l)               )


          # Eexecute fullchem_HetDropChem if criteria are satisfied
          # note: skip if lwc is very small, which will blow up equations#
          if ( ( size_res                                           ) && ( state_met%iswater(i, j)                             ) && ( temp                    > 268.15_fp                ) && ( lwc                     > 1.0e-20_fp               ) ) 

             fullchem_hetdropchem( i         = i, j         = j, l         = l, sr        = sr, input_opt = input_opt, state_met = state_met, state_chm = state_chm               )

             # Add result as an enhancement to O2 metal catalysis rate
             # as a 1st order reaction
             k_cld(3) = k_cld(3) + sr
          end

          # Stop timer
          if ( input_opt%usetimers ) 
             timer_end( timername = "     rconst", inloop    = true, threadnum = thread, rc        =  rc                                 )
          end
       end

       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       # Start KPP main timer and prepare arrays
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

       # Zero out dummy species index in KPP
       for f = 1: nfam
          kppid = pl_kpp_id(f)
          if ( kppid > 0 ) c(kppid) = 0.0_dp end
       end

       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       # Update reaction rates
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

       # Start timer
       if ( input_opt%usetimers ) 
          timer_start( timername = "     rconst", inloop    = true, threadnum = thread, rc        =  rc                                 )
       end

       # Update the array of rate constants
       update_rconst( )

       # Stop timer
       if ( input_opt%usetimers ) 
          timer_end( timername = "     rconst", inloop    = true, threadnum = thread, rc        =  rc                                   )
       end

       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       # HISTORY (aka netCDF diagnostics)
       #
       # Archive KPP reaction rates [s - 1]
       # See gckpp_Monitor.F90 for a list of chemical reactions
       #
       # NOTE: In KPP 2.5.0 + , VAR and FIX are now private to the integrator
       # and point to C.  Therefore, pass C(1:NVAR) instead of VAR and
       # C(NVAR + 1:NSPEC) instead of FIX to routine FUN.
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       if ( state_diag%archive_rxnrate ) 
          # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
          # Get equation rates (Aout)
          # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
          fun( v       = c(1:nvar), f       = c(nvar + 1:nspec), rct     = rconst, vdot    = vloc, aout    = aout                                          )

          for s = 1: state_diag%map_rxnrate%nslots
             n = state_diag%map_rxnrate%slot2id(s)
             state_diag%rxnrate(i, j, l, s) = aout(n)
          end
       end

       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       # Set options for the KPP Integrator (M. J. Evans)
       #
       # NOTE: Because RCNTRL(3) is set to an array value that
       # depends on (I, J, L), we must declare RCNTRL as PRIVATE
       # within the OpenMP parallel loop and define it just
       # before the to to Integrate. (bmy, 3 / 24 / 16)
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

       # Zero all slots of RCNTRL
       rcntrl    = 0.0_fp

       # Initialize Hstart (the starting value of the integration step
       # size with the value of Hnew (the last predicted but not yet 
       # taken timestep)  saved to the the restart file.
       rcntrl(3) = state_chm%kpphvalue(i, j, l)

       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       # Integrate the box forwards
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

       # Start timer
       if ( input_opt%usetimers ) 
          timer_start( timername = "     integrate 1", inloop    =  true, threadnum = thread, rc        = rc                                  )
       end

       # Call the KPP integrator
       integrate( tin,    tout,    icntrl, rcntrl, istatus, rstate, ierr                        )

       # Stop timer
       if ( input_opt%usetimers ) 
          timer_end( timername = "     integrate 1", inloop    = true, threadnum = thread, rc        = rc                                    )
       end

       # Print grid box indices to screen if integrate failed
       if ( ierr < 0 ) 
          #FIXME write(6, * ) "### integrate returned error at: ", i, j, l
       end

#if defined( model_geos ) || defined( model_wrf )
       # Print grid box indices to screen if integrate failed
       if ( ierr < 0 ) 
          #FIXME write(6, * ) "### integrate returned error at: ", i, j, l
          if ( associated(state_diag%kpperror) ) 
             state_diag%kpperror(i, j, l) = state_diag%kpperror(i, j, l) + 1.0
          end
       end
#endif # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = end
       # HISTORY: Archive KPP solver diagnostics
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       if ( state_diag%archive_kppdiags ) 

          # # of integrator calls
          if ( state_diag%archive_kppintcounts ) 
             state_diag%kppintcounts(i, j, l) = istatus(1)
          end

          # # of times Jacobian was constructed
          if ( state_diag%archive_kppjaccounts ) 
             state_diag%kppjaccounts(i, j, l) = istatus(2)
          end

          # # of internal timesteps
          if ( state_diag%archive_kpptotsteps ) 
             state_diag%kpptotsteps(i, j, l) = istatus(3)
          end

          # # of accepted internal timesteps
          if ( state_diag%archive_kpptotsteps ) 
             state_diag%kppaccsteps(i, j, l) = istatus(4)
          end

          # # of rejected internal timesteps
          if ( state_diag%archive_kpptotsteps ) 
             state_diag%kpprejsteps(i, j, l) = istatus(5)
          end

          # # of LU - decompositions
          if ( state_diag%archive_kppludecomps ) 
             state_diag%kppludecomps(i, j, l) = istatus(6)
          end

          # # of forward and backwards substitutions
          if ( state_diag%archive_kppsubsts ) 
             state_diag%kppsubsts(i, j, l) = istatus(7)
          end

          # # of singular - matrix decompositions
          if ( state_diag%archive_kppsmdecomps ) 
             state_diag%kppsmdecomps(i, j, l) = istatus(8)
          end
       end

       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       # Try another time if it failed
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       if ( ierr < 0 ) 

#if defined( model_geos ) || defined( model_wrf )
          # Save a copy of the C vector (GEOS and WRF only)
          localc    = c
#endif

          # Reset first time step and start concentrations
          # Retry the integration with non - optimized settings
          rcntrl(3) = 0.0_dp
          c         = 0.0_dp

          # Update rates again
          update_rconst( )

          # Start timer
          if ( input_opt%usetimers ) 
             timer_start( timername = "     integrate 2", inloop    =  true, threadnum = thread, rc        = rc                               )
          end

          # Integrate again
          integrate( tin,    tout,    icntrl, rcntrl, istatus, rstate, ierr                     )

          # Stop timer
          if ( input_opt%usetimers ) 
             timer_end( timername = "     integrate 2", inloop    =  true, threadnum = thread, rc        = rc                                 )
          end

          # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
          # HISTORY: Archive KPP solver diagnostics
          # This time, add to the existing value
          # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
          if ( state_diag%archive_kppdiags ) 

             # # of integrator calls
             if ( state_diag%archive_kppintcounts ) 
                state_diag%kppintcounts(i, j, l) = state_diag%kppintcounts(i, j, l) + istatus(1)
             end

             # # of times Jacobian was constructed
             if ( state_diag%archive_kppjaccounts ) 
                state_diag%kppjaccounts(i, j, l) = state_diag%kppjaccounts(i, j, l) + istatus(2)
             end

             # # of internal timesteps
             if ( state_diag%archive_kpptotsteps ) 
                state_diag%kpptotsteps(i, j, l) = state_diag%kpptotsteps(i, j, l) + istatus(3)
             end

             # # of accepted internal timesteps
             if ( state_diag%archive_kpptotsteps ) 
                state_diag%kppaccsteps(i, j, l) = state_diag%kppaccsteps(i, j, l) + istatus(4)
             end

             # # of rejected internal timesteps
             if ( state_diag%archive_kpptotsteps ) 
                state_diag%kpprejsteps(i, j, l) = state_diag%kpprejsteps(i, j, l) + istatus(5)
             end

             # # of LU - decompositions
             if ( state_diag%archive_kppludecomps ) 
                state_diag%kppludecomps(i, j, l) = state_diag%kppludecomps(i, j, l) + istatus(6)
             end

             # # of forward and backwards substitutions
             if ( state_diag%archive_kppsubsts ) 
                state_diag%kppsubsts(i, j, l) = state_diag%kppsubsts(i, j, l) + istatus(7)
             end

             # # of singular - matrix decompositions
#             IF ( State_Diag%Archive_KppSmDecomps ) THEN
#                State_Diag%KppSmDecomps(I, J, L) = #                State_Diag%KppSmDecomps(I, J, L) + ISTATUS(8)
#             ENDIF
          end

          # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
          # Exit upon the second failure
          # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
          if ( ierr < 0 ) 
             #FIXME write(6,     "(a   )" ) "## integrate failed twice ### "
             #FIXME write(errmsg, "(a, i3)" ) "integrator error code :", ierr
#if defined( model_geos ) || defined( model_wrf )
             if ( input_opt%kppstop ) 
                error_stop(errmsg, "integrate_kpp")
             # Revert to start values
             else
                c = localc
             end
             if ( associated(state_diag%kpperror) ) 
                state_diag%kpperror(i, j, l) = state_diag%kpperror(i, j, l) + 1.0
             end
#else
             # Set a flag to break out of loop gracefully
             # NOTE: You can set a GDB breakpoint here to examine the error
             #$OMP CRITICAL
             failed2x = true

             # Print concentrations at trouble box KPP error
             print * , repeat( "###", 79 )
             print * , "### kpp debug output#"
             print * , "### species concentrations at problem box ", i, j, l
             print * , repeat( "###", 79 )
             for n = 1: nspec
                print * , "### ", c(n), trim( adjustl( spc_names(n) ) )
             end

             # Print rate constants at trouble box KPP error
             print * , repeat( "###", 79 )
             print * , "### kpp debug output#"
             print * , "### reaction rates at problem box ", i, j, l
             print * , repeat( "###", 79 )
             for n = 1: nreact
                print * , rconst(n), trim( adjustl( eqn_names(n) ) )
             end
             #$OMP END CRITICAL
#endif
          end

       end


       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       # Continue upon successful return...
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

       # Revert Alkalinity (only when using sulfur chemistry in KPP)
       if ( ! state_chm%do_sulfatemod_seasalt ) 
          fullchem_convertequivtoalk()
       end

       # Save Hnew (the last predicted but not taken step) from the 3rd slot
       # of RSTATE into State_Chm so that it can be written to the restart
       # file.  For simulations that are broken into multiple stages,
       # Hstart will be initialized to the value of Hnew from the restart
       # file at startup (see above).
       state_chm%kpphvalue(i, j, l) = rstate(3)

       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       # Check we have no negative values and copy the concentrations
       # calculated from the C array back into State_Chm%Species%Conc
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

       # Loop over KPP species
       for n = 1: nspec

          # GEOS - Chem species ID
          spcid = state_chm%map_kppspc(n)

          # Skip if this is not a GEOS - Chem species
          if ( spcid < = 0 ) cycle end

          # Set negative concentrations to zero
          c(n) = max( c(n), 0.0_dp )

          # Copy concentrations back into State_Chm%Species
          state_chm%species(spcid)%conc(i, j, l) = real( c(n), kind = fp )

       end

       if ( input_opt%usetimers ) 

          # Stop main KPP timer
          timer_end( timername  = " - > kpp", inloop     = true, threadnum  = thread, rc         = rc                                   )

          # Start Prod / Loss timer
          timer_start( timername = " - > prod / loss diags", inloop    = true, threadnum = thread, rc        = rc                                  )
       end

#ifdef bpch_diag
#ifdef tomas
       #always calculate rate for TOMAS
       for f = 1: nfam

          # Determine dummy species index in KPP
          kppid =  pl_kpp_id(f)

          # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
          # FOR TOMAS MICROPHYSICS:
          #
          # Obtain P / L with a unit [kg S] for tracing
          # gas - phase sulfur species production (SO2, SO4, MSA)
          # (win, 8 / 4 / 09)
          # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 

          # Calculate H2SO4 production rate [kg s - 1] in each
          # time step (win, 8 / 4 / 09)
          if ( trim(fam_names(f)) = = "pso4" )  end
             # Hard - coded MW
             h2so4_rate(i, j, l) = c(kppid) / avo * 98.0e - 3_fp * state_met%airvol(i, j, l)    * 1.0e + 6_fp / dt

            if ( h2so4_rate(i, j, l) < 0.0d0) 
              #FIXME write( * , * ) "h2so4_rate negative in fullchem_mod.f90##", i, j, l, "was:", h2so4_rate(i, j, l), "  setting to 0.0d0"
              h2so4_rate(i, j, l) = 0.0d0
            end
          end
       end

#endif
#endif

#ifdef model_geos
       # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
       # Archive NOx lifetime [h]
       # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
       if ( state_diag%archive_noxtau || state_diag%archive_tropnoxtau ) 
          fun( c(1:nvar), c(nvar + 1:nspec), rconst, vloc,      aout = aout,       vdotout = vdotout )
          noxtau = vdotout(ind_no) + vdotout(ind_no2) + vdotout(ind_no3) + 2.0 * vdotout(ind_n2o5) + vdotout(ind_clno2) + vdotout(ind_hno2) + vdotout(ind_hno4)
          noxconc = c(ind_no) + c(ind_no2) + c(ind_no3) + 2.0 * c(ind_n2o5) + c(ind_clno2) + c(ind_hno2) + c(ind_hno4)
          # NOx chemical lifetime per grid cell
          if ( state_diag%archive_noxtau ) 
             noxtau = ( noxconc / ( - 1.0_f4 * noxtau) ) / 3600.0_f4
             if ( noxtau > 0.0_f4 ) 
                state_diag%noxtau(i, j, l) = min(1.0e10_f4, max(1.0e-10_f4, noxtau))
             else
                state_diag%noxtau(i, j, l) = max( - 1.0e10_f4, min( - 1.0e-10_f4, noxtau))
             end
          end
          # NOx chemical lifetime per trop. column
          if ( state_diag%archive_tropnoxtau ) 
             nox_weight = ( noxconc ) * state_met%airden(i, j, l) * state_met%delp_dry(i, j, l)
             nox_tau_weighted = ( noxconc / ( - 1.0_f4 * noxtau * 3600.0_f4 ) ) * nox_weight
             if ( abs(nox_tau_weighted) < 1.0e8 ) 
               nox_tau_weighted = ( nint(nox_tau_weighted) * 1.0e6 ) * 1.0e-6_f4
             else
                if ( nox_tau_weighted > 0.0 ) 
                   nox_tau_weighted = 1.0e8
                else
                   nox_tau_weighted = - 1.0e8
                end
             end
             if ( state_met%introposphere(i, j, l) ) 
               tropv_nox_mass(i, j) = tropv_nox_mass(i, j) + nox_weight
               tropv_nox_tau(i, j)  = tropv_nox_tau(i, j) + nox_tau_weighted
             end
          end
       end
#endif # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = end
       # HISTORY (aka netCDF diagnostics)
       #
       # Prod and loss of families or species [molec / cm3 / s]
       #
       # NOTE: KppId is the KPP ID # for each of the prod and loss
       # diagnostic species.  This is the value used to index the
       # KPP "C" array (in module gckpp_Global.F90).
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

       # Chemical loss of species or families [molec / cm3 / s]
       if ( state_diag%archive_loss ) 
          for s = 1: state_diag%map_loss%nslots
             kppid = state_diag%map_loss%slot2id(s)
             state_diag%loss(i, j, l, s) = c(kppid) / dt
          end
       end

       # Chemical production of species or families [molec / cm3 / s]
       if ( state_diag%archive_prod ) 
          for s = 1: state_diag%map_prod%nslots
             kppid = state_diag%map_prod%slot2id(s)
             state_diag%prod(i, j, l, s) = c(kppid) / dt
          end
       end

       # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
       # Archive prod / loss fields for the TagCO simulation [molec / cm3 / s]
       # (In practice, we only need to do this from benchmark simulations)
       # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
       if ( state_diag%archive_prodcofromch4 || state_diag%archive_prodcofromnmvoc ) 

          # Total production of CO
          pco_tot   = c(id_pco) / dt

          # Loss of CO from CH4
          lch4      = c(id_lch4) / dt

          # P(CO)_CH4 is LCH4. Cap so that it is never greater
          # than total P(CO) to prevent negative P(CO)_NMVOC.
          pco_ch4   = min( lch4, pco_tot )

          # P(CO) from NMVOC is the remaining P(CO)
          pco_nmvoc = pco_tot - pco_ch4

          # Archive P(CO) from CH4 for tagCO simulations
          if ( state_diag%archive_prodcofromch4 ) 
             state_diag%prodcofromch4(i, j, l) = pco_ch4
          end

          # Archive P(CO) from NMVOC for tagCO simulations
          if ( state_diag%archive_prodcofromnmvoc ) 
             state_diag%prodcofromnmvoc(i, j, l) = pco_nmvoc
          end

       end

       # Stop Prod / Loss timer
       if ( input_opt%usetimers ) 
          timer_end( timername =  " - > prod / loss diags", inloop    = true, threadnum = thread, rc        = rc                                    )
       end

       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       # HISTORY (aka netCDF diagnostics)
       #
       # Write out OH reactivity.  The OH reactivity is defined here as the
       # inverse of its life - time. In a crude ad - hoc approach, manually add
       # all OH reactants (ckeller, 9 / 20 / 2017)
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       if ( state_diag%archive_ohreactivity ) 

          # Start timer
          if ( input_opt%usetimers ) 
             timer_start( timername = " - > oh reactivity diag", inloop    = true, threadnum = thread, rc        = rc                               )
          end

          # Archive OH reactivity diagnostic
          get_ohreactivity ( c, rconst, ohreact )
          state_diag%ohreactivity(i, j, l) = ohreact

          # Stop timer
          if ( input_opt%usetimers ) 
             timer_end( timername = " - > oh reactivity diag", inloop    = true, threadnum = thread, rc        = rc                                 )
          end
       end
    end
    end
    end
    #$OMP END PARALLEL DO

    # Stop timer
    if ( input_opt%usetimers ) 
       timer_end( " - > flexchem loop", rc )
    end

    # Compute sum of in - loop timers
    if ( input_opt%usetimers ) 
       timer_sum_loop( " - > init kpp",            rc )
       timer_sum_loop( " - > het chem rates",      rc )
       timer_sum_loop( " - > photolysis rates",    rc )
       timer_sum_loop( " - > kpp",                 rc )
       timer_sum_loop( "     rconst",              rc )
       timer_sum_loop( "     integrate 1",         rc )
       timer_sum_loop( "     integrate 2",         rc )
       timer_sum_loop( " - > prod / loss diags",     rc )
       timer_sum_loop( " - > oh reactivity diag",  rc )
    end

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Return gracefully if integration failed 2x anywhere
    # (as we cannot break out of a parallel do loop#)
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    if ( failed2x ) 
       errmsg = "kpp failed to converge after 2 iterations#"
       gc_error( errmsg, rc, thisloc )
       return
    end

#if defined( model_geos )
    if ( state_diag%archive_tropnoxtau ) 
       where (tropv_nox_mass > 0.0_f4 )
          state_diag%tropnoxtau = tropv_nox_tau / tropv_nox_mass
       end
    end
#endif # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = end
    # Archive OH, HO2, O1D, O3P concentrations after solver
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    if ( do_diag_oh_ho2_o1d_o3p ) 

       # Save OH, HO2, O1D, O3P for the ND43 diagnostic
       # NOTE: These might not be needed for netCDF, as they will already
       # have been archived in State_Chm%Species%Conc output.
       diag_oh_ho2_o1d_o3p( input_opt,  state_chm, state_diag, state_grid, state_met, rc )

       # Trap potential errors
       if ( rc / = gc_success )  end
          errmsg = "error encountered in "diag_oh_ho2_o1d_o3p"#"
          gc_error( errmsg, rc, thisloc )
          return
       end

       if ( prtdebug ) 
          debug_msg( "### do_fullchem: after ohsave" )
       end
    end

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Archive quantities for computing OH metrics
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    diag_metrics( input_opt,  state_chm, state_diag, state_grid, state_met, rc                            )

    # Trap potential errors
    if ( rc / = gc_success )  end
       errmsg = "error encountered in "diag_mean_oh_and_ch4"
       gc_error( errmsg, rc, thisloc )
       return
    end

    if ( prtdebug ) 
       debug_msg( "### do_fullchem: after diag_metrics" )
    end

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Archive concentrations after chemistry
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    if ( state_diag%archive_concafterchem ) 
       # Point to mapping obj specific to SpeciesConc diagnostic collection
       mapdata = > state_diag%map_concafterchem

       #$OMP PARALLEL DO#$OMP DEFAULT( SHARED )#$OMP PRIVATE( N, S   )
       for s = 1: mapdata%nslots
          n = mapdata%slot2id(s)
          state_diag%concafterchem(:, :, :, s) = state_chm%species(n)%conc(:, :, :)
       end
       #$OMP END PARALLEL DO

       # Free pointer
       mapdata = > null()
    end

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Convert species back to original units (ewl, 8 / 16 / 16)
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    convert_spc_units( input_opt, state_chm,  state_grid, state_met, origunit,  rc )

    if ( rc / = gc_success )  end
       errmsg = "unit conversion error#"
       gc_error( errmsg, rc, "fullchem_mod.f90" )
       return
    end

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Apply high - altitude active nitrogen partitioning and H2SO4
    # photolysis approximations outside the chemistry grid
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    ucx_nox( input_opt, state_chm, state_grid, state_met )
    if ( prtdebug ) 
       debug_msg( "### chemdr: after ucx_nox" )
    end

    ucx_h2so4phot( input_opt, state_chm, state_grid, state_met )
    if ( prtdebug ) 
       debug_msg( "### chemdr: after ucx_h2so4phot" )
    end

    # Set State_Chm arrays for surface J - values used in HEMCO and
    # saved to restart file
    if ( rxn_o3_1 > = 0 )  end
       state_chm%joh(:, :) = zpj(1, rxn_o3_1, :, :)
    end
    if ( rxn_no2 > = 0 )  end
       state_chm%jno2(:, :) = zpj(1, rxn_no2, :, :)
    end

    # Set FIRSTCHEM = .FALSE. -  - we have gone thru one chem step
    firstchem = false

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: PrintFirstTimeInfo
#
# #DESCRIPTION: Prints information the first time DO_FULLCHEM is called
#\\
#\\
# #INTERFACE:
#
  subroutine printfirsttimeinfo( input_opt, state_chm, firstchem, do_photchem )
#
# #USES:
#
#
# #INPUT PARAMETERS:
#
    type(optinput), intent(in) :: input_opt     # Input Options object
    type(chmstate), intent(in) :: state_chm     # Chemistry State object
    logical,        intent(in) :: firstchem     # Is this the first call?
    logical,        intent(in) :: do_photchem   # Is photolysis turned on?
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC

    # Exit if we are not on the first chemistry timestep
    if ( ! firstchem ) return

    # Print on the root CPU only
    if ( input_opt%amiroot ) 

       # Write header
       #FIXME write( 6, "(a)" ) repeat( " = ", 79 )
       #FIXME write( 6, 100 )
 100   format("do_fullchem : interface between geos - chem and the solver")
       #FIXME write( 6, "(a)" )

       # Print where sulfur seasalt chemistry is done
       if ( state_chm%do_sulfatemod_seasalt ) 
          #FIXME write( 6, 110 )
 110      format( " * sulfur sea salt chemistry is computed in sulfate_mod" )
       else
          #FIXME write( 6, 120 )
 120      format( " * sulfur sea salt chemistry is computed in kpp" )
       end

       # Print where sulfur in - cloud chemistry is done
       if ( state_chm%do_sulfatemod_cld ) 
          #FIXME write( 6, 130 )
 130      format( " * sulfur in - cloud chemistry is computed in sulfate_mod" )
       else
          #FIXME write( 6, 140 )
 140      format( " * sulfur in - cloud chemistry is computed in kpp" )
       end

       # Print the status of photolysis: on or off
       if ( do_photchem ) 
          #FIXME write( 6, 150 )
 150      format(  " * photolysis is activated -  - rates computed by fast - jx" )
       else
          #FIXME write( 6, 160 )
 160      format(  " * photolysis has been deactivated for testing purposes"  )
       end

       # Write footer
       #FIXME write( 6, "(a)" ) repeat( " = ", 79 )
    end

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: set_sulfur_chem_rates
#
# #DESCRIPTION: Calls functions from the KPP rate - law library to compute
#  rates for sulfate chemistry reactions.  These are passed to the KPP
#  mechanism.
#\\
#\\
# #INTERFACE:
#
  subroutine set_sulfur_chem_rates( i,          j,          l, input_opt,  state_chm,  state_diag, state_grid, state_met,  size_res, rc                                      )
#
# #USES:
#
#
# #INPUT PARAMETERS:
#
    integer,        intent(in)    :: i, j, l      # X, Y, Z gridbox indices
    type(optinput), intent(in)    :: input_opt    # Input Options object
    type(grdstate), intent(in)    :: state_grid   # Grid State object
    type(metstate), intent(in)    :: state_met    # Meteorology State object
#
# #INPUT / OUTPUT PARAMETERS:
#
    type(chmstate), intent(inout) :: state_chm    # Chemistry State object
    type(dgnstate), intent(inout) :: state_diag   # Diagnostics State object
#
# #OUTPUT PARAMETERS:
#
    logical,        intent(out)   :: size_res     # Should we HetDropChem?
    integer,        intent(out)   :: rc           # success or failure#
#
# #REMARKS:
#  The routines below are based on or meant to replace reactions computed
#  outside of KPP within sulfate_mod.
#
#  Rates are set for the following:
#  1) Cloud sulfur chemistry
#  2) (If Dust acid uptake) Dust acid uptake reactions
# - MSL, 29 - Mar - 2021, 7 - May - 2021
#
#  Seasalt aerosol chemistry reaction rate - law functions are now contained
#  in module fullchem_RateLawFuncs.F90.
#
#  NOTE This routine defines the variables State_Chm%HSO3_aq and
#  State_Chm%SO3aq.  Therefore, we must Set_Sulfur_Chem_Rates after
#  Set_KPP_GridBox_Values, but before fullchem_SetStateHet.  Otherwise we
#  will not be able to copy State_Chm%HSO3_aq to State_Het%HSO3_aq and
#  State_Chm%SO3_aq to State_Het%SO3_aq properly.
# -  - Mike Long, Bob Yantosca (08 Mar 2022)
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC

    # Initialize
    rc       = gc_success
    size_res = false

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Do this when KPP is handling aqueous sulfur chemistry
    #
    # From Alexander et al., buffering capacity (or alkalinity) of
    # sea - salt aerosols is equal to 0.07 equivalents per kg dry sea salt
    # emitted Gurciullo et al., 1999.0 JGR 104(D17) 21, 719 - 21, 731.0
    # tdf; MSL
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    if ( ! state_chm%do_sulfatemod_seasalt ) 

       # Convert alkalinity from [molec / cm3] to equivalents
       fullchem_convertalktoequiv()

       # Compute reaction rates for aqueous sulfur chemistry
       # (i.e. S(IV) - >S(VI), HCl,  and HNO3)
       fullchem_sulfuraqchem( i          = i, j          = j, l          = l, input_opt  = input_opt, state_chm  = state_chm, state_grid = state_grid, state_met  = state_met, rc         = rc                          )
    end

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Do this when KPP is handling SO2 cloud chemistry ...
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    if ( ! state_chm%do_sulfatemod_cld ) 

       # Compute reaction rates [1 / s] for sulfate in cloud for KPP chem mech
       fullchem_sulfurcldchem(  i          = i, j          = j, l          = l, input_opt  = input_opt, state_chm  = state_chm, state_diag = state_diag, state_grid = state_grid, state_met  = state_met, size_res   = size_res, rc         = rc                        )

       # Update HSO3 - and SO3 -  - concentrations [molec / cm3]
       state_chm%fupdatehobr(i, j, l) = 1.0_fp
       state_chm%fupdatehocl(i, j, l) = 1.0_fp

    end

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: Set_Kpp_GridBox_Values
#
# #DESCRIPTION: Populates KPP variables in the gckpp_Global.F90 module
#  for a particular (I, J, L) grid box.
#\\
#\\
# #INTERFACE:
#
  subroutine set_kpp_gridbox_values( i,         j,         l, input_opt, state_chm, state_grid, state_met, rc                          )
#
# #USES:
#
#
# #INPUT PARAMETERS:
#
    integer,        intent(in)  :: i, j, l
    type(optinput), intent(in)  :: input_opt
    type(chmstate), intent(in)  :: state_chm
    type(grdstate), intent(in)  :: state_grid
    type(metstate), intent(in)  :: state_met
#
# #OUTPUT PARAMETERS:
#
    integer,        intent(out) :: rc
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #LOCAL VARIABLES:
#
    # Scalars
    real(f8)           :: consexp, vpresh2o

    # Strings
    character(len = 255) :: errmsg, thisloc

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # set_kpp_gridbox_values begins here#
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # Initialization
    rc      = gc_success
    na      = state_chm%naerotype
    errmsg  = ""
    thisloc = " - > at set_kpp_gridbox_values (in module geoscore / fullchem_mod.f90)"

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Populate global variables in gckpp_Global.F90
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # Solar quantities
    suncos          = state_met%suncosmid(i, j)

    # Pressure and density quantities
    numden          = state_met%airnumden(i, j, l)
    h2o             = state_met%avgw(i, j, l) * numden
    press           = get_pcenter( i, j, l )

    # Temperature quantities
    temp            = state_met%t(i, j, l)
    inv_temp        = 1.0_dp   / temp
    temp_over_k300  = temp     / 300.0_dp
    k300_over_temp  = 300.0_dp / temp
    sr_temp         = sqrt( temp )
    four_r_t        = 4.0_dp * con_r    * temp
    four_rgaslatm_t = 4.0_dp * rgaslatm * temp
    eight_rstarg_t  = 8.0_dp * rstarg   * temp

    # Relative humidity quantities
    consexp         = 17.2693882_dp * (temp - 273.16_dp) / (temp - 35.86_dp)
    vpresh2o        = consvap * exp( consexp ) / temp
    relhum          = ( h2o / vpresh2o ) * 100_dp

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: Diag_OH_HO2_O1D_O3P
#
# #DESCRIPTION: Archives the chemical production of OH, HO2, O1D, O3P.
#\\
#\\
# #INTERFACE:
#
  subroutine diag_oh_ho2_o1d_o3p( input_opt,  state_chm, state_diag, state_grid, state_met, rc )
#
# #USES:
#
#
# #INPUT PARAMETERS:
#
    type(optinput), intent(in)    :: input_opt   # Input Options object
    type(grdstate), intent(in)    :: state_grid  # Grid State object
#
# #INPUT / OUTPUT PARAMETERS:
#
    type(metstate), intent(inout) :: state_met   # Meteorology State object
    type(chmstate), intent(inout) :: state_chm   # Chemistry State object
    type(dgnstate), intent(inout) :: state_diag  # Diagnostics State object
#
# #OUTPUT PARAMETERS:
#
    integer,        intent(out)   :: rc          # Success or failure
#
# #REMARKS:
#  This routine replaces both OHSAVE and DIAGOH.  Those routines were needed
#  for SMVGEAR, when we had separate arrays for the non - advected species.
#  But now, all species are stored in State_Chm%Species%Conc, so the various
#  arrays (SAVEOH, SAVEHO2, etc.) are no longer necessary.  We can now just
#  just get values directly from State_Chm%Species%Conc.
#
#  Also note: for the netCDF diagnostics, we have removed multiplication by
#  LTOH etc arrays.  These are almost always set between 0 and 24.0
#
# #REVISION HISTORY:
#  06 Jan 2015 - R. Yantosca - Initial version
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #LOCAL VARIABLES:
#

    # Scalars

    # Strings
    character(len = 255) :: errmsg, thisloc

    # Pointers
    real(fp),      pointer  :: airnumden(:, :, :)
    type(spcconc), pointer  :: spc      (:    )

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # diag_oh_ho2_o1d_o3p begins here#
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # Assume success
    rc = gc_success

    # Point to the array of species concentrations
    airnumden = > state_met%airnumden
    spc       = > state_chm%species

    # Zero the netCDF diagnostic arrays (if activated) above the
    # tropopause or mesopause to avoid having leftover values
    # from previous timesteps
#ifdef model_geos
    if ( state_diag%archive_o3concafterchem  ) 
       state_diag%o3concafterchem  = 0.0_f4
    end
    if ( state_diag%archive_ro2concafterchem ) 
       state_diag%ro2concafterchem = 0.0_f4
    end
#endif
    if ( state_diag%archive_ohconcafterchem  ) 
       state_diag%ohconcafterchem  = 0.0_f4
    end
    if ( state_diag%archive_ho2concafterchem ) 
       state_diag%ho2concafterchem = 0.0_f4
    end
    if ( state_diag%archive_o1dconcafterchem ) 
       state_diag%o1dconcafterchem = 0.0_f4
    end
    if ( state_diag%archive_o3pconcafterchem ) 
       state_diag%o3pconcafterchem = 0.0_f4
    end

#$OMP PARALLEL DO#$OMP DEFAULT( SHARED )#$OMP PRIVATE( I, J, L )#$OMP SCHEDULE( DYNAMIC )
      for l = 1: state_grid%nz
      for j = 1: state_grid%ny
      for i = 1: state_grid%nx

         # Skip non - chemistry boxes
         if ( ! state_met%inchemgrid(i, j, l) ) 
            # Skip to next grid box
            cycle
         end

         # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
         # OH concentration [molec / cm3]
         # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
         if ( ok_oh ) 

            # HISTORY (aka netCDF diagnostics)
            if ( state_diag%archive_ohconcafterchem ) 
               state_diag%ohconcafterchem(i, j, l) = spc(id_oh)%conc(i, j, l)
            end
#ifdef model_geos
            if ( state_diag%archive_o3concafterchem ) 
               state_diag%o3concafterchem(i, j, l) = spc(id_o3)%conc(i, j, l)
            end
            if ( archive_ro2concafterchem ) 
               if ( id_a3o2 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_a3o2)%conc(i, j, l) end
               if ( id_ato2 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_ato2)%conc(i, j, l) end
               if ( id_b3o2 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_b3o2)%conc(i, j, l) end
               if ( id_bro2 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_bro2)%conc(i, j, l) end
               if ( id_eto2 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_eto2)%conc(i, j, l) end
               if ( id_ho2 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_ho2)%conc(i, j, l) end
               if ( id_limo2 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_limo2)%conc(i, j, l) end
               if ( id_mo2 > 0  )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_mo2)%conc(i, j, l) end
               if ( id_pio2 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_pio2)%conc(i, j, l) end
               if ( id_po2 > 0  )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_po2)%conc(i, j, l) end
               if ( id_prn1 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_prn1)%conc(i, j, l) end
               if ( id_r4n1 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_r4n1)%conc(i, j, l) end
               if ( id_r4o2 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_r4o2)%conc(i, j, l) end
               if ( id_tro2 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_tro2)%conc(i, j, l) end
               if ( id_xro2 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_xro2)%conc(i, j, l) end
               if ( id_ihoo1 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_ihoo1)%conc(i, j, l) end
               if ( id_ihoo4 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_ihoo4)%conc(i, j, l) end
               if ( id_ichoo > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_ichoo)%conc(i, j, l) end
               if ( id_ihpoo1 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_ihpoo1)%conc(i, j, l) end
               if ( id_ihpoo2 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_ihpoo2)%conc(i, j, l) end
               if ( id_ihpoo3 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_ihpoo3)%conc(i, j, l) end
               if ( id_iepoxaoo > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_iepoxaoo)%conc(i, j, l) end
               if ( id_iepoxboo > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_iepoxboo)%conc(i, j, l) end
               if ( id_c4hvp1 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_c4hvp1)%conc(i, j, l) end
               if ( id_c4hvp2 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_c4hvp2)%conc(i, j, l) end
               if ( id_hpald1oo > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_hpald1oo)%conc(i, j, l) end
               if ( id_hpald2oo > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_hpald2oo)%conc(i, j, l) end
               if ( id_isopnoo1 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_isopnoo1)%conc(i, j, l) end
               if ( id_isopnoo2 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_isopnoo2)%conc(i, j, l) end
               if ( id_ino2d > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_ino2d)%conc(i, j, l) end
               if ( id_ino2b > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_ino2b)%conc(i, j, l) end
               if ( id_idhnboo > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_idhnboo)%conc(i, j, l) end
               if ( id_idhndoo1 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_idhndoo1)%conc(i, j, l) end
               if ( id_idhndoo2 > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_idhndoo2)%conc(i, j, l) end
               if ( id_ihpnboo > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_ihpnboo)%conc(i, j, l) end
               if ( id_ihpndoo > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_ihpndoo)%conc(i, j, l) end
               if ( id_icnoo > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_icnoo)%conc(i, j, l) end
               if ( id_idnoo > 0 )state_diag%ro2concafterchem(i, j, l) = state_diag%ro2concafterchem(i, j, l) + spc(id_idnoo)%conc(i, j, l) end
            end
#endif

         end

         # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
         # HO2 concentration [v / v]
         # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
         if ( ok_ho2 ) 
            if ( state_diag%archive_ho2concafterchem ) 
               state_diag%ho2concafterchem(i, j, l) = ( spc(id_ho2)%conc(i, j, l) /   airnumden(i, j, l)      )
            end

         end

         # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
         # O1D concentration [molec / cm3]
         # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
         if ( ok_o1d ) 
            if ( state_diag%archive_o1dconcafterchem ) 
               state_diag%o1dconcafterchem(i, j, l) = spc(id_o1d)%conc(i, j, l)
            end
         end


         # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
         # O3P concentration [molec / cm3]
         # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
         if ( ok_o3p ) 
            if ( state_diag%archive_o3pconcafterchem ) 
               state_diag%o3pconcafterchem(i, j, l) = spc(id_o3p)%conc(i, j, l)
            end
         end

      end
      end
      end
#$OMP END PARALLEL DO

      # Free pointers
      airnumden = > null()
      spc       = > null()

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: Diag_Metrics
#
# #DESCRIPTION: Computes mass - weighted mean OH columns (full - atmosphere and
#  trop - only) that are needed to compute the overall mean OH concentration.
#  This is used as a metric as to how reactive, or "hot" the chemistry
#  mechanism is.
#\\
#\\
# #INTERFACE:
#
  subroutine diag_metrics( input_opt,  state_chm, state_diag, state_grid, state_met, rc                        )
#
# #USES:
#
#
# #INPUT PARAMETERS:
#
    type(optinput), intent(in)    :: input_opt    # Input Options object
    type(chmstate), intent(in)    :: state_chm    # Chemistry State object
    type(grdstate), intent(in)    :: state_grid   # Grid State object
    type(metstate), intent(in)    :: state_met    # Meteorology State object
#
# #INPUT / OUTPUT PARAMETERS:
#
    type(dgnstate), intent(inout) :: state_diag   # Diagnostics State object
#
# #OUTPUT PARAMETERS:
#
    integer,        intent(out)   :: rc           # Success or failure?
#
# #REMARKS:
#  References:
#  (1) Prather, M. and C. Spivakovsky, "Tropospheric OH and
#       the lifetimes of hydrochlorofluorocarbons", JGR,
#       Vol 95, No. D11, 18723 - 18729, 1990.0
#  (2) Lawrence, M.G, Joeckel, P, and von Kuhlmann, R., "What
#       does the global mean OH concentraton tell us?",
#       Atm. Chem. Phys, 1, 37 - 49, 2001.0
#  (3) WMO / UNEP Scientific Assessment of Ozone Depletion: 2010
#
# #REVISION HISTORY:
#  18 Aug 2020 - R. Yantosca - Initial version
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
    real(f8), parameter :: m3tocm3        = 1.0e + 6_f8
#
# #LOCAL VARIABLES:
#
    # SAVEd scalars
    logical,  save      :: first          = true
    integer,  save      :: id_oh          = - 1
    real(f8), save      :: mcm3tokgm3_oh  = - 1.0_f8

    # Scalars
    real(f8)            :: airmass_m,   airmass_kg,  airmassfull
    real(f8)            :: airmasstrop, ktrop,       lossohbych4
    real(f8)            :: lossohbymcf, ohconc_mcm3, ohmasswgt
    real(f8)            :: ohmassfull,  ohmasstrop,  volume

    # Strings
    character(len = 255)  :: errmsg,      thisloc

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # compute_mean_oh_and_ch4 begins here#
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # Initialize
    rc      = gc_success
    errmsg  = ""
    thisloc = " - > at compute_mean_oh (in module geoscore / diagnostics_mod.f90)"

    # Exit if we have not turned on the Metrics collection
    if ( ! state_diag%archive_metrics ) return

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # First - time setup
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    if ( first ) 

       # Get the species ID for OH
       id_oh  = ind_("oh")
       if ( id_oh < 0 ) 
          errmsg = "oh is not a defined species in this simulation###"
          gc_error( errmsg, rc, thisloc )
          return
       end

       # Convert [molec OH cm - 3] -  - > [kg OH m - 3]
       mcm3tokgm3_oh  = m3tocm3 * ( state_chm%spcdata(id_oh)%info%mw_g * 1.0e-3_f8 ) / avo


       # Reset first - time flag
       first  = false
    end

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Loop over surface boxes and compute mean OH in columns
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #$OMP PARALLEL DO#$OMP DEFAULT( SHARED                                                   )#$OMP PRIVATE( I,           J,           L,           airMass_kg        )#$OMP PRIVATE( airMass_m,   airMassFull, airMassTrop, Ktrop             )#$OMP PRIVATE( LossOHbyCH4, LossOHbyMCF, OHconc_mcm3, OHmassWgt         )#$OMP PRIVATE( OHmassFull,  OHmassTrop,  volume                         )#$OMP SCHEDULE( DYNAMIC, 4                                              )
    for j = 1: state_grid%ny
    for i = 1: state_grid%nx

       # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
       # Zero column - specific quantities
       # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
       airmass_kg  = 0.0_f8
       airmass_m   = 0.0_f8
       airmassfull = 0.0_f8
       airmasstrop = 0.0_f8
       ktrop       = 0.0_f8
       lossohbych4 = 0.0_f8
       lossohbymcf = 0.0_f8
       ohconc_mcm3 = 0.0_f8
       ohmasswgt   = 0.0_f8
       ohmassfull  = 0.0_f8
       ohmasstrop  = 0.0_f8
       volume      = 0.0_f8

       # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
       # Loop over the number of levels in the full - atmosphere column
       # Limit the computations to boxes in the chemistry grid
       # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
       for l = 1, state_met%chemgridlev(i:j)

          # Compute box volume [cm3], and air mass ([molec] and [kg])
          # Note: air mass in [molec] is also the atmospheric burden of
          # methyl chloroform (aka MCF, formula = CH3CCl3), since we assume
          # a unif orm mixing ratio ( = 1) of MCF in air. end
          volume      = state_met%airvol(i, j, l)    * m3tocm3
          airmass_m   = state_met%airnumden(i, j, l) * volume
          airmass_kg  = airmass_m / xnumolair

          # OH concentration [molec cm - 3]
          ohconc_mcm3 = state_chm%species(id_oh)%conc(i, j, l)

          # Airmass - weighted OH [kg air * (kg OH  m - 3)]
          ohmasswgt   = airmass_kg * ( ohconc_mcm3  * mcm3tokgm3_oh  )

          # Sum the air mass, mass - weighted CH4,
          # and mass - weighted OH in the full - atm column
          airmassfull = airmassfull + airmass_kg
          ohmassfull  = ohmassfull  + ohmasswgt

          # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
          # Only do the following for tropospheric boxes ...
          # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
          if ( state_met%introposphere(i, j, l) ) 

             # Sum the air mass, mass - weighted CH4,
             # and mass - weighted OH in the trop - only column
             airmasstrop = airmasstrop + airmass_kg
             ohmasstrop  = ohmasstrop  + ohmasswgt

             # Compute CH4 loss rate in troposphere
             # Ktrop (Arrhenius parameter) has units [cm3 / molec / s]
             # OHconc has units [molec / cm3]
             # AirMass has units [molec]
             # Resultant units of CH4 loss rate = [molec / s]
             ktrop = 2.45e - 12_f8 * exp( - 1775.0_f8 / state_met%t(i, j, l) )
             lossohbych4 = lossohbych4 + ( ktrop * ohconc_mcm3 * airmass_m )

             # Compute MCF loss rate in the troposphere
             # Ktrop (Arrhenius parameter) has units [cm3 / molec / s]
             # OHconc has units [molec / cm3]
             # AirMass has units [molec]
             # Resultant units of MCF loss rate = [molec / s]
             ktrop = 1.64e - 12_f8 * exp( - 1520.0_f8 / state_met%t(i, j, l) )
             lossohbymcf = lossohbymcf + ( ktrop * ohconc_mcm3 * airmass_m )

          end
       end

       # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
       # HISTORY (aka netCDF diagnostics)
       # Air mass [kg], full - atmosphere and trop - only column sums
       # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
       if ( state_diag%archive_airmasscolumnfull ) 
          state_diag%airmasscolumnfull(i, j) = airmassfull
       end

       if ( state_diag%archive_airmasscolumntrop ) 
          state_diag%airmasscolumntrop(i, j) = airmasstrop
       end

       # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
       # HISTORY (aka netCDF diagnostics)
       # Airmass - weighted mean OH [kg air * (kg OH m - 3)]
       # Full - atmosphere and trop - only column sums
       # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
       if ( state_diag%archive_ohwgtbyairmasscolumnfull ) 
          state_diag%ohwgtbyairmasscolumnfull(i, j) = ohmassfull
       end

       if ( state_diag%archive_ohwgtbyairmasscolumntrop ) 
          state_diag%ohwgtbyairmasscolumntrop(i, j) = ohmasstrop
       end

       # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
       # HISTORY (aka netCDF diagnostics)
       #
       # OH loss by CH4 + OH in troposphere [molec / s] and
       # OH loss by MCF + OH in troposphere [molec / s]
       # Full - atmosphere and trop - only column sums
       # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
       if ( state_diag%archive_lossohbych4columntrop ) 
          state_diag%lossohbych4columntrop(i, j) = lossohbych4
       end

       if ( state_diag%archive_lossohbymcfcolumntrop ) 
          state_diag%lossohbymcfcolumntrop(i, j) = lossohbymcf
       end

    end
    end
    #$OMP END PARALLEL DO

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: init_fullchem
#
# #DESCRIPTION: Subroutine Init\_FullChem is used to allocate arrays for the
#  KPP solver.
#\\
#\\
# #INTERFACE:
#
  subroutine init_fullchem( input_opt, state_chm, state_diag, rc )
#
# #USES:
#
#
# #INPUT PARAMETERS:
#
    type(optinput), intent(in)  :: input_opt   # Input Options object
    type(chmstate), intent(in)  :: state_chm   # Diagnostics State object
    type(dgnstate), intent(in)  :: state_diag  # Diagnostics State object
#
# #OUTPUT PARAMETERS:
#
    integer,        intent(out) :: rc          # Success or failure?
#
# #REVISION HISTORY:
#  14 Dec 2015 - M.S. Long - Initial version
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #LOCAL VARIABLES:
#
    # Scalars

    # Strings
    character(len = 255) :: errmsg,   thisloc

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # init_fullchem begins here#
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # Assume success
    rc       = gc_success

    # Do the following only if it is a full - chemistry simulation
    # NOTE: If future specialty simulations use the KPP solver,
    # modify the IF statement accordingly to allow initialization
    if ( ! input_opt%its_a_fullchem_sim ) return

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Initialize variables
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    errmsg   = ""
    thisloc  = " - > at init_fullchem (in module geoscore / fullchem_mod.f90)"
    prtdebug = ( input_opt%lprt && input_opt%amiroot )

    # Debug output
    if ( prtdebug ) 
       #FIXME write( 6, 100 )
100    format( " - init_fullchem: allocating arrays" )

       #FIXME write( 6 , "(a)" ) " kpp reaction reference "
       for n = 1: nreact
          #FIXME write( 6, "(i8, a3, a85)" ) n, " | ", eqn_names(n)
       end
    end

    # Initialize species flags
    id_ch4      = ind_( "ch4", "a"     ) # CH4 advected species
    id_ho2      = ind_( "ho2"          )
    id_nh3      = ind_( "nh3"          )
    id_o3p      = ind_( "o"            )
    id_o1d      = ind_( "o1d"          )
    id_oh       = ind_( "oh"           )
    id_salaal   = ind_( "salaal"       )
    id_salcal   = ind_( "salcal"       )
    id_so4      = ind_( "so4"          )
    id_salc     = ind_( "salc"         )

#ifdef model_geos
    # ckeller
    id_o3       = ind_( "o3"           )
    id_a3o2     = ind_( "a3o2"         )
    id_ato2     = ind_( "ato2"         )
    id_bro2     = ind_( "bro2"         )
    id_eto2     = ind_( "eto2"         )
    id_limo2    = ind_( "limo2"        )
    id_mo2      = ind_( "mo2"          )
    id_pio2     = ind_( "pio2"         )
    id_po2      = ind_( "po2"          )
    id_prn1     = ind_( "prn1"         )
    id_r4n1     = ind_( "r4n1"         )
    id_r4o2     = ind_( "r4o2"         )
    id_tro2     = ind_( "tro2"         )
    id_xro2     = ind_( "xro2"         )
    id_ihoo1    = ind_( "ihoo1"        )
    id_ihoo4    = ind_( "ihoo4"        )
    id_ihcoo    = ind_( "ihcoo"        )
    id_ihpoo1   = ind_( "ihpoo1"       )
    id_ihpoo2   = ind_( "ihpoo2"       )
    id_ihpoo3   = ind_( "ihpoo3"       )
    id_iepoxaoo = ind_( "iepoxaoo"     )
    id_iepoxboo = ind_( "iepoxboo"     )
    id_c4hvp1   = ind_( "c4hvp1"       )
    id_c4hvp2   = ind_( "c4hvp2"       )
    id_hpald1oo = ind_( "hpald1oo"     )
    id_hpald2oo = ind_( "hpald2oo"     )
    id_isopnoo1 = ind_( "isopnoo1"     )
    id_isopnoo2 = ind_( "isopnoo2"     )
    id_ino2b    = ind_( "ino2b"        )
    id_ino2d    = ind_( "ino2d"        )
    id_idhnboo  = ind_( "idhnboo"      )
    id_idhndoo1 = ind_( "idhndoo1"     )
    id_idhndoo2 = ind_( "idhndoo2"     )
    id_ihpnboo  = ind_( "ihpnboo"      )
    id_ihpndoo  = ind_( "ihpndoo"      )
    id_icnoo    = ind_( "icnoo"        )
    id_idnoo    = ind_( "idnoo"        )
#endif

    # Set flags to denote if each species is defined
    ok_ho2      = ( id_ho2 > 0         )
    ok_o1d      = ( id_o1d > 0         )
    ok_o3p      = ( id_o3p > 0         )
    ok_oh       = ( id_oh  > 0         )

    # Should we archive OH, HO2, O1D, O3P diagnostics?
    do_diag_oh_ho2_o1d_o3p = (#ifdef model_geos
                               state_diag%archive_o3concafterchem || state_diag%archive_ro2concafterchem || #endif
                               state_diag%archive_ohconcafterchem || state_diag%archive_ho2concafterchem || state_diag%archive_o1dconcafterchem || state_diag%archive_o3pconcafterchem          )

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Save physical parameters from the species database into KPP arrays
    # in gckpp_Global.F90.  These are for the hetchem routines.
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    for kppid = 1: state_chm%nkppspc + state_chm%nomitted
       n                  = state_chm%map_kppspc(kppid)
       if ( n > 0 ) 
          mw(kppid)       = state_chm%spcdata(n)%info%mw_g
          sr_mw(kppid)    = sqrt( mw(kppid ) )
          henry_k0(kppid) = state_chm%spcdata(n)%info%henry_k0
          henry_cr(kppid) = state_chm%spcdata(n)%info%henry_cr
       end
    end

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Allocate arrays
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # Initialize
    id_pco = - 1
    id_lch4 = - 1

    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    # Pre - store the KPP indices for each KPP prod / loss species or family
    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    if ( nfam > 0 ) 

       # Allocate mapping array for KPP Ids for ND65 bpch diagnostic
       allocate( pl_kpp_id( nfam ), stat = rc )
       gc_checkvar( "fullchem_mod.f90:pl_kpp_id", 0, rc )
       if ( rc / = gc_success ) return end

       # Loop over all KPP prod / loss species
       for n = 1: nfam

          # NOTE: KppId is the KPP ID # for each of the prod and loss
          # diagnostic species.  This is the value used to index the
          # KPP "VAR" array (in module gckpp_Global.F90).
          kppid = ind_( trim ( fam_names(n) ), "k" )

          # Find the KPP Id corresponding to PCO and LCH4
          # so that we can save output for tagged CO simulations
          if ( trim( fam_names(n) ) = = "pco"  ) id_pco  = kppid end
          if ( trim( fam_names(n) ) = = "lch4" ) id_lch4 = kppid end

          # Exit if an invalid ID is encountered
          if ( kppid < = 0 )  end
             errmsg = "invalid kpp id for prod / loss species: "            / / trim( fam_names(n) )
             gc_error( errmsg, rc, thisloc )
             return
          end

          # If the species ID is OK, save in ND65_Kpp_Id
          pl_kpp_id(n) = kppid
       end

    end

    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    # If we are archiving the P(CO) from CH4 and from NMVOC from a fullchem
    # simulation for the tagCO simulation, throw an error if we cannot find
    # the PCO or LCH4 prod / loss families in this KPP mechanism.
    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    if ( state_diag%archive_prodcofromch4 || state_diag%archive_prodcofromnmvoc ) 

       if ( id_pco < 1 ) 
          errmsg = "could not find pco in list of kpp families#  This   " / / "is needed to archive the prodcofromch4 diagnostic#"
          gc_error( errmsg, rc, thisloc )
          return
       end

       if ( id_lch4 < 1 ) 
          errmsg = "could not find lch4 in list of kpp families#  This "  / / "is needed to archive the prodcofromnmvoc diagnostic#"
          gc_error( errmsg, rc, thisloc )
          return
       end

    end

    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    # Initialize sulfate chemistry code (cf Mike Long)
    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    fullchem_initsulfurchem( rc )
    if ( rc / = gc_success )  end
       errmsg = "error encountered in "fullchem_initsulfurcldchem"#"
       gc_error( errmsg, rc, thisloc )
       return
    end

    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    # Initialize dust acid uptake code (Mike Long, Bob Yantosca)
    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    if ( input_opt%ldstup ) 
       aciduptake_initdustchem( rc )
       if ( rc / = gc_success )  end
          errmsg = "error encountered in "aciduptake_initdustchem"#"
          gc_error( errmsg, rc, thisloc )
          return
       end
    end

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: cleanup_fullchem
#
# #DESCRIPTION: Subroutine Cleanup\_FullChem deallocates module variables.
#\\
#\\
# #INTERFACE:
#
  subroutine cleanup_fullchem( rc )
#
# #USES:
#
#
# #OUTPUT PARAMETERS:
#
    integer, intent(out) :: rc          # Success or failure?
#
# #REVISION HISTORY:
#  24 Aug 2016 - M. Sulprizio - Initial version
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # cleanup_fullchem begins here#
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # Initialize
    rc = gc_success

    if ( allocated( pl_kpp_id ) ) 
       deallocate( pl_kpp_id, stat = rc  )
       gc_checkvar( "fullchem_mod.f90:pl_kpp_id", 2, rc )
       if ( rc / = gc_success ) return end
    end

    if ( allocated( jvcountday ) ) 
       deallocate( jvcountday, stat = rc  )
       gc_checkvar( "fullchem_mod.f90:jvcountday", 2, rc )
       if ( rc / = gc_success ) return end
    end

    if ( allocated( jvsumday ) ) 
       deallocate( jvsumday, stat = rc  )
       gc_checkvar( "fullchem_mod.f90:jvcountday", 2, rc )
       if ( rc / = gc_success ) return end
    end

    if ( allocated( jvcountmon ) ) 
       deallocate( jvcountmon, stat = rc  )
       gc_checkvar( "fullchem_mod.f90:jvcountmon", 2, rc )
       if ( rc / = gc_success ) return end
    end

    if ( allocated( jvsummon ) ) 
       deallocate( jvsummon, stat = rc  )
       gc_checkvar( "fullchem_mod.f90:jvcountmon", 2, rc )
       if ( rc / = gc_success ) return end
    end

  end
#EOC
end