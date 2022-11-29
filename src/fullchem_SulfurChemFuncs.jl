
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #MODULE: fullchem_SulfurChemFuncs
#
# #DESCRIPTION: FlexChem module for multiphase sulfate chemistry, via KPP.
#\\
#\\
# #INTERFACE:

# module fullchem_sulfurchemfuncs
#
# #USES:
#
#   private
#
# #PUBLIC MEMBER FUNCTIONS:
#
#   public :: fullchem_convertalktoequiv
#   public :: fullchem_convertequivtoalk
#   public :: fullchem_hetdropchem
#   public :: fullchem_initsulfurchem
#   public :: fullchem_sulfuraqchem
#   public :: fullchem_sulfurcldchem
#
# #PUBLIC TYPES:
#
  # Species ID flags
#
# #DEFINED_PARAMETERS
#
#   real(fp),  parameter   :: tcvv_s    = airmw / 32e + 0_fp # hard - coded MW
#   real(fp),  parameter   :: tcvv_n    = airmw / 14e + 0_fp # hard - coded MW
#   real(fp),  parameter   :: smallnum  = 1e-20_fp
#   real(fp),  parameter   :: cm3perm3  = 1.0e6_fp

tcvv_s    = airmw / 32e + 0 # hard - coded MW
tcvv_n    = airmw / 14e + 0 # hard - coded MW
smallnum  = 1e-20
cm3perm3  = 1.0e6

# contains
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: fullchem_ConvertAlkToEquiv
#
# #DESCRIPTION: Converts sea salt alkalinity to equivalents.  Abstracted
#  out from fullchem_mod.F90 to prevent compilation conflicts for other
#  KPP chemical mechanisms
#\\
#\\
# #INTERFACE:
#
function fullchem_convertalktoequiv()
#
# #USES:
#
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
    c(ind_salaal) = c(ind_salaal) * ( mw(ind_salaal) * 7.0e-5)
    c(ind_salcal) = c(ind_salcal) * ( mw(ind_salcal) * 7.0e-5)
end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: fullchem_ConvertEquivToAlk
#
# #DESCRIPTION: Converts sea salt alkalinity to equivalents.  Abstracted
#  out from fullchem_mod.F90 to prevent compilation conflicts for other
#  KPP chemical mechanisms
#\\
#\\
# #INTERFACE:
#
function fullchem_convertequivtoalk()
#
# #USES:
#
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
    c(ind_salaal) = c(ind_salaal) / ( mw(ind_salaal) * 7.0e-5_fp )
    c(ind_salcal) = c(ind_salcal) / ( mw(ind_salcal) * 7.0e-5_fp )
end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: fullchem_SulfurAqchem
#
# #DESCRIPTION: Main aqueous / aerosol chemistry driver routine.  Sets up the
#  vector of aqueous chemistry rates for the KPP chemistry solver.
#\\
#\\
# #INTERFACE:
#
function fullchem_sulfuraqchem( i,j,l, input_opt, state_chm,  state_grid, state_met, rc)
#
# #USES:
#
#
# #INPUT PARAMETERS:
#
#     integer,        intent(in)    :: i, j, l    # Lon, lat, level indices
#     type(metstate), intent(in)    :: state_met  # Meteorology State object
#     type(grdstate), intent(in)    :: state_grid  # Grid State object
#     type(optinput), intent(in)    :: input_opt  # Input Options object
#
# #INPUT / OUTPUT PARAMETERS:
#
#     type(chmstate), intent(in)    :: state_chm  # Chemistry State object
#
# OUTPUT PARAMETERS:
#
#     integer,        intent(out)   :: rc         # Success or failure
#
# #REMARKS:
#
#  # Reaction List (by K_MT() index)
#    1) SO2 + O3 + 2SALAAL -  - > SO4mm + O2 : From Sulfate_mod - 24 Mar 2021
#
# #REVISION HISTORY:
#  24 Mar 2021 - M. Long - Initial Version
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #LOCAL VARIABLES:

    # Scalars
#     real(fp)           :: k_ex

    # Strings
#     character(len = 255) :: errmsg
#     character(len = 255) :: thisloc

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # fullchem_sulfuraqchem begins here#
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # Initialize
    rc            = gc_success
    k_ex          = 0.0 
    k_mt          = 0.0 
    salaal_gt_0_1 = ( state_chm%species(id_salaal)%conc(i, j, l) > 0.1  )
    salcal_gt_0_1 = ( state_chm%species(id_salcal)%conc(i, j, l) > 0.1  )

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Reaction rates [1 / s] for fine sea salt alkalinity (aka SALAAL)
    #
    # K_MT(1) : SALAAL + SO2 + O3 = SO4 - SALAAL
    # K_MT(2) : SALAAL + HCl      = SALACL
    # K_MT(3) : SALAAL + HNO3     = NIT
    #
    # NOTE: SALAAL_gt_0_1 prevents div - by - zero errors
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    # SALAAL + SO2 + O3 = SO4 - SALAAL
    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    if ( salaal_gt_0_1 ) 

      # 1st order uptake
      k_ex = ars_l1k( area = state_chm%wetaeroarea(i, j, l, 11), radius = state_chm%aeroradi(i, j, l, 11), Γ  = 0.11 , srmw   = sr_mw(ind_so2)                              )

      # Assume SO2 is limiting, so recompute rxn rate accordingly
      k_mt1 = kiir1ltd( c(ind_so2), c(ind_o3), k_ex ) / state_chm%species(id_salaal)%conc(i, j, l)
    end

    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    # SALAAL + HCL = SALACL
    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    if ( salaal_gt_0_1 ) 

      # 1st order uptake
      k_ex = ars_l1k( area   = state_chm%wetaeroarea(i, j, l, 11), radius = state_chm%aeroradi(i, j, l, 11), Γ  = 0.07 , srmw   = sr_mw(ind_hcl)                              )

      # Assume HCl is limiting, so recompute reaction rate accordingly
      k_mt2 = kiir1ltd( c(ind_hcl), c(ind_salaal), k_ex )
    end

    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    # SALAAL + HNO3 = NIT
    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    if ( salaal_gt_0_1 ) 

      # 1st order uptake
      k_ex = ars_l1k( area   = state_chm%wetaeroarea(i, j, l, 11), radius = state_chm%aeroradi(i, j, l, 11), Γ  = 0.5 , srmw   = sr_mw(ind_hno3)                             )

      # Assume HNO3 is limiting, so recompute reaction rate accordingly
      k_mt3 = kiir1ltd( c(ind_hno3), c(ind_salaal), k_ex )
    end

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Reaction rates [1 / s] for coarse sea salt alkalinity (aka SALAAL)
    #
    # K_MT(4) : SALCAL + SO2 + O3 = SO4s - SALCAL
    # K_MT(5) : SALCAL + HCl      = SALCCL
    # K_MT(6) : SALCAL + HNO3     = NITs
    #
    # NOTE: SALCAL_gt_0_1 prevents div - by - zero errors
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    # SALCAL + SO2 + O3 = SO4s - SALCAL
    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    if ( salcal_gt_0_1 ) 

      # 1st order uptake
      k_ex = ars_l1k( area   = state_chm%wetaeroarea(i, j, l, 12), radius = state_chm%aeroradi(i, j, l, 12), Γ  = 0.11 , srmw   = sr_mw(ind_so2)                              )

      # Assume SO2 is limiting, so recompute rxn rate accordingly
      k_mt4 = kiir1ltd( c(ind_so2), c(ind_o3), k_ex ) / state_chm%species(id_salcal)%conc(i, j, l)
    end

    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    # SALCAL + HCl = SALCCL
    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    if ( salcal_gt_0_1 ) 

      # 1st order uptake
      k_ex = ars_l1k( area   = state_chm%wetaeroarea(i, j, l, 12), radius = state_chm%aeroradi(i, j, l, 12), Γ  = 0.07 , srmw   = sr_mw(ind_hcl)                              )

      # Assume HCl is limiting, so recompute rxn rate accordingly
      k_mt5 = kiir1ltd( c(ind_hcl), c(ind_salcal), k_ex )
    end

    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    # SALCAL + HNO3 = NITs
    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    if ( salcal_gt_0_1 ) 

      # 1st order uptake
      k_ex = ars_l1k( area   = state_chm%wetaeroarea(i, j, l, 12), radius = state_chm%aeroradi(i, j, l, 12), Γ  = 0.5 , srmw   = sr_mw(ind_hno3)                             )

      # Assume HNO3 is limiting, so recompute rxn rate accordingly
      k_mt(6) = kiir1ltd( c(ind_hno3), c(ind_salcal), k_ex )
    end
    return rc
end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: fullchem_SulfurCldChem
#
# #DESCRIPTION: Routine that compute reaction rates for sulfur chemistry
#  in cloud, so that these can be passed to the KPP chemical solver.
#\\
#\\
# #INTERFACE:
#
function fullchem_sulfurcldchem( i,         j,         l, input_opt,  state_chm, state_diag, state_grid, state_met, size_res, rc                                      )
#
# #USES:
#
#
# #INPUT PARAMETERS:
#
    type(optinput), intent(in)    :: input_opt   # Input Options object
    type(grdstate), intent(in)    :: state_grid  # Grid State object
    type(metstate), intent(in)    :: state_met   # Meteorology State object
    integer,        intent(in)    :: i, j, l
#
# #INPUT / OUTPUT PARAMETERS:
#
    type(chmstate), intent(inout) :: state_chm   # Chemistry State object
    type(dgnstate), intent(inout) :: state_diag  # Diagnostics State object
#
# #OUTPUT PARAMETERS:
#
    logical,        intent(out)   :: size_res    # Should we HetDropChem?
    integer,        intent(out)   :: rc          # Success or failure?
#
# #REVISION HISTORY:
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #LOCAL VARIABLES:
#
    character(len = 63)  :: origunit

    # Strings
    character(len = 255) :: errmsg, thisloc

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # fullchem_sulfurcldchem begins here#
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # Initialize
    rc       = gc_success
    size_res = false
    errmsg   = ""
    thisloc  = " - > at fullchem_sulfurcldchem (in kpp / fullchem / fullchem_sulfurchemfuncs.f90"

    # Should we print debug output?
    prtdebug             = ( input_opt%lprt && input_opt%amiroot )

    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    # SO2 chemistry
    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    set_so2( i          = i, j          = j, l          = l, input_opt  = input_opt, state_chm  = state_chm, state_diag = state_diag, state_grid = state_grid, state_met  = state_met, size_res   = size_res, rc         = rc                                           )

    # Trap potential errors
    if ( rc / = gc_success ) 
 end
       errmsg = "error encountered in "set_so2"#"
       gc_error( errmsg, rc, thisloc )
       return
    end

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: fullchem_HetDropChem
#
# #DESCRIPTION: Subroutine HET\_DROP\_CHEM estimates the in - cloud sulfate
#  production rate in heterogeneous cloud droplets based on the Yuen et al.,
#  1996 parameterization. (bec, 6 / 16 / 11)
#\\
#\\
# #INTERFACE:
#
function fullchem_hetdropchem( i,         j,         l,         sr, input_opt, state_met, state_chm          )
#
# #USES:
#
#
# #INPUT PARAMETERS:
#
    # integer,        intent(in)    :: i, j, l
    # type(optinput), intent(in)    :: input_opt   # Input Options object
    # type(metstate), intent(in)    :: state_met   # Meteorology State object
#
# #INPUT / OUTPUT PARAMETERS:
#
    # type(chmstate), intent(inout) :: state_chm   # Chemistry State object
#
# #OUTPUT PARAMETERS:
#
    # real(fp),       intent(out)   :: sr          # Sulfate production rate
#
# #REVISION HISTORY:
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
    # Dry sea - salt density [kg / m3]
    # real(dp), parameter :: ss_den        = 2200.0 

    # sigma of the size distribution for sea - salt (Jaegle et al., 2011)
    # real(fp), parameter :: sig_s         = 1.8e + 0 

    # geometric dry mean diameters [m] for computing lognormal size distribution
    # real(dp), parameter :: rg_s          = 0.4e-6  #(Jaegle et a., 2011)
    # real(dp), parameter :: rg_d2         = 1.5e-6  #(Ginoux et al., 2001)
    # real(dp), parameter :: rg_d3         = 2.5e-6 
    # real(dp), parameter :: rg_d4         = 4.0e - 6 

    # To prevent multiple divisions
    # real(dp), parameter   :: three_fourths = 3.0  / 4.0 
    # real(dp), parameter   :: nine_halves   = 9.0  / 2.0 

#
# #LOCAL VARIABLES:
#
    # real(dp)              :: α_nh3,  α_so2, α_h2o2
    # real(dp)              :: α_hno3, α_b,   α_cn
    # real(dp)              :: α_w,    α_so4, sum_gas
    # real(dp)              :: h,          ndss,      cn
    # real(dp)              :: w,          k,         arg
    # real(dp)              :: dtchem,     apv,       dsvi
    # real(dp)              :: b,          nh3,       so2
    # real(dp)              :: h2o2,       hno3,      so4
    # real(dp)              :: cnss,       mw_so4,    mw_salc
    # real(dp)              :: cvf,        r1,        r2
    # real(dp)              :: xx,         fc,        lst
    # real(dp)              :: xx1,        xx2,       xx3
    # real(dp)              :: xx4,        xx5,       gnh3

    # Pointers
    # real(fp), pointer     :: ad(:, :, :)
    # real(fp), pointer     :: airden(:, :, :)
    # real(fp), pointer     :: airvol(:, :, :)
    # real(fp), pointer     :: omega(:, :, :)
    # real(fp), pointer     :: u(:, :, :)
    # real(fp), pointer     :: v(:, :, :)

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # het_drop_chem begins here#
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # Initialize pointers
    ad     = > state_met%ad
    airden = > state_met%airden
    airvol = > state_met%airvol
    omega  = > state_met%omega
    u      = > state_met%u
    v      = > state_met%v

    # Zero output argument for safety"s sake
    sr     =  0.0 

    # Zero / initialize local variables for safety"s sake
    arg    =  0.0 
    b      =  0.0 
    cn     =  0.0 
    cvf    =  1.0e3_fp * airmw / ( airden(i, j, l) * avo ) # molec / cm3 - > v / v
    dsvi   =  0.0 
    dtchem =  get_ts_chem()                              # seconds
    gnh3   =  0.0 
    k      =  0.0 
    lst    =  0.0 
    r1     =  0.0 
    r2     =  0.0 
    w      =  0.0 
    xx     =  0.0 
    xx1    =  0.0 
    xx2    =  0.0 
    xx3    =  0.0 
    xx4    =  0.0 
    xx5    =  0.0 

    # FC is guaranteed to be > 1e-4, because HET_DROP_CHEM
    # is not called otherwise (bmy, 07 Oct 2021)
    fc     =  state_met%cldf(i, j, l)

    ## <<>> set the input units# EITHER CONVERT IN THE ROUTINE OR
    ## <<>> CONVERT BEFOREHAND. BUT EVERYTHING IS CURRENTLY mcl / cm3
    ## <<>> AND HET_DROP_CHEM EXPECTS V / V

    # XX * are calculated below to be consistent with
    # Sulfate_Mod(). Values are different when
    # computed with KPP - based variables. HET_DROP_CHEM()
    # could use some attention to make is consistent with
    # KPP.
    #
    # NOTE: Use function SafeExp, which will prevent the exponential from
    # blowing up.  Also if the entire expression will evaluate to zero
    #  skip the exponential, which is more computationally efficient.
    # -  - Bob Yantosca, 14 Oct 2021
    #

    # SO2 + H2O2
    r1  = c(ind_so2)  * cvf
    r2  = c(ind_h2o2) * cvf
    k   = k_cld(1)    / cvf / fc
    arg = ( r1 - r2 ) * ( k * dtchem )
    if ( issafeexp( arg ) && abs( arg ) > 0.0  ) 
       xx  = exp( arg )
       xx1 = ( r1 * r2 ) * ( xx - 1.0  ) / ( ( r1 * xx ) - r2 )
    else
       xx1 = whenexpcantbedone( r1, r2, k, dtchem )
    end

    # SO2 + O3
    r2  = c(ind_o3) * cvf
    k   = k_cld(2)  / cvf / fc
    arg = ( r1 - r2 ) * ( k * dtchem )
    if ( issafeexp( arg ) && abs( arg ) > 0.0  ) 
       xx = exp( arg )
       xx2 = ( r1 * r2 ) * ( xx - 1.0  ) / ( ( r1 * xx ) - r2 )
    else
       xx2 = whenexpcantbedone( r1, r2, k, dtchem )
    end

    # Metal catalyzed oxidation of SO2 pathway
    k   = - k_cld(3) / fc
    arg = k * dtchem
    xx3 = 0.0 
    if ( issafeexp( arg ) ) 
       xx  = exp( arg )
       xx3 = r1 * ( 1.0  - xx )
    end

    # HSO3 - + HOCl and SO3 -  - + HOCl
    r1  = c(ind_so2) * cvf * state_chm%hso3_aq(i, j, l)
    r2  = c(ind_hocl)  * cvf
    k   = hocluptkbyhso3m(state_het) / cvf
    arg = ( r1 - r2 ) * ( k * dtchem )
    if ( issafeexp( arg ) && abs( arg ) > 0.0  ) 
       xx  = exp( arg )
       xx4 = ( r1 * r2 )  * ( xx - 1.0  ) / ( ( r1 * xx ) - r2 )
    else
       xx4 = whenexpcantbedone( r1, r2, k, dtchem )
    end

    # SO3 -  - + HOCl (add to HSO3 - + HOCl rate)
    r1  = c(ind_so2) * cvf * state_chm%so3_aq(i, j, l)
    k   = hocluptkbyso3mm(state_het) / cvf
    arg = ( r1 - r2 ) * ( k * dtchem )
    if ( issafeexp( arg ) && abs( arg ) > 0.0  ) 
       xx  = exp( arg )
       xx4 = xx4 + ( ( r1 * r2 ) * ( xx - 1.0_fp ) / ( ( r1 * xx ) - r2 ) )
    else
       xx4 = xx4 + whenexpcantbedone( r1, r2, k, dtchem )
    end

    # HSO3 - + HOBr
    r1  = c(ind_so2) * cvf * state_chm%hso3_aq(i, j, l)
    r2  = c(ind_hobr)  * cvf
    k   = hobruptkbyhso3m(state_het) / cvf
    arg = ( r1 - r2 ) * ( k * dtchem )
    if ( issafeexp( arg ) && abs( arg ) > 0.0  ) 
       xx  = exp( arg )
       xx5 = ( r1 * r2 ) * ( xx - 1.0_fp ) / ( ( r1 * xx ) - r2 )
    else
       xx5 = whenexpcantbedone( r1, r2, k, dtchem )
    end

    # SO3 -  - + HOBr (add to HSO3 - + HOBr rate)
    r1  = c(ind_so2) * cvf * state_chm%so3_aq(i, j, l)
    k   = hobruptkbyso3mm(state_het) / cvf
    arg = ( r1 - r2 ) * ( k * dtchem )
    if ( issafeexp( arg ) && abs( arg ) > 0.0  ) 
       xx  = exp( arg )
       xx5 = xx5 + ( ( r1 * r2 ) * ( xx - 1.0  ) / ( ( r1 * xx ) - r2 ) )
    else
       xx5 = xx5 + whenexpcantbedone( r1, r2, k, dtchem )
    end

    # Sum of all rates
    lst = xx1 + xx2 + xx3 + xx4 + xx5

    #### Debug print
    #IF (I == 12 && J == 7 && L == 1) THEN
    #   #FIXME write( * , * ) "<<>> XX: ", XX1, XX2, XX3, XX4, XX5
    #ENDIF

    if ( lst > r1 ) 
       xx1 = ( r1 * xx1 ) / lst
       xx2 = ( r1 * xx2 ) / lst
       xx3 = ( r1 * xx3 ) / lst
       xx4 = ( r1 * xx4 ) / lst
       xx5 = ( r1 * xx5 ) / lst
       lst = xx1 + xx2 + xx3 + xx4 + xx5
     end

    # Convert gas phase concentrations from [v / v] to [pptv]
    nh3  = state_chm%species(id_nh3)%conc(i, j, l) * cvf * 1.0e + 12 
    so2  = max( c(ind_so2) * cvf - ( lst * fc ), 1.0e-20  ) * 1.0e + 12 
    h2o2 = c(ind_h2o2) * cvf * 1.0e12 
    hno3 = c(ind_hno3) * cvf * 1.0e12 

    # Set molecular weight local variables
    mw_so4  = state_chm%spcdata(id_so4)%info%mw_g
    mw_salc = state_chm%spcdata(id_salc)%info%mw_g

    # Convert sulfate aerosol concentrations from [v / v] to [ug / m3]
    so4 = ( c(ind_so4) * cvf * ad(i, j, l) * 1.0e + 9  ) / ( ( airmw / mw_so4 ) * airvol(i, j, l) )

    # Convert in cloud sulfate production rate from [v / v / timestep] to
    # [ug / m3 / timestep]
    b  = ( lst * ad(i, j, l) * 1.0e + 9  ) / ( ( airmw / mw_so4 ) * airvol(i, j, l) )

    # Convert coarse - mode aerosol concentrations from [v / v] to [# / cm3]
    # based on equation in Hofmann, Science, 1990.0
    # First convert from [v / v] to [kg / m3 air]
    cnss = state_chm%species(id_salc)%conc(i, j, l) * cvf * ad(i, j, l) / ( ( airmw / mw_salc ) * airvol(i, j, l) )

    # Now convert from [kg / m3 air] to [# / cm3 air]
    # Sea - salt
    arg  = nine_halves * ( log( sig_s ) )^2
    ndss = ( three_fourths * cnss                           ) / ( pi * ss_den * rg_s^3 * safeexp( arg, 0.0  ) ) * 1.0e - 6 

    # Total coarse mode number concentration [# / cm3]
    cn = ndss # sea - salt

    # Determine regression coefficients based on the local SO2 concentration
    if ( so2 < = 200.00_fp ) 
 end
       α_b    =  0.5318 
       α_nh3  = - 1.67e - 7 
       α_so2  =  2.59e - 6 
       α_h2o2 = - 1.77e - 7 
       α_hno3 = - 1.72e - 7 
       α_w    =  1.22e - 6 
       α_cn   =  4.58e - 6 
       α_so4  = - 1.00e - 5 
    elseif ( so2 > 200.00  && so2 < = 500.0  ) 
 end
       α_b    =  0.5591 
       α_nh3  =  3.62e - 6 
       α_so2  =  1.66e - 6 
       α_h2o2 =  1.06e - 7 
       α_hno3 = - 5.45e - 7 
       α_w    = - 5.79e - 7 
       α_cn   =  1.63e - 5 
       α_so4  = - 7.40e - 6 
    elseif ( so2 > 500.0  && so2 < 1000.0  ) 
       α_b    =  1.1547 
       α_nh3  = - 4.28e - 8 
       α_so2  = - 1.23e - 7 
       α_h2o2 = - 9.05e - 7 
       α_hno3 =  1.73e - 7 
       α_w    =  7.22e - 6 
       α_cn   =  2.44e - 5 
       α_so4  =  3.25e - 5 
    else                          # SO2 > 1000
       α_b    =  1.1795 
       α_nh3  =  2.57e - 7 
       α_so2  = - 5.54e - 7 
       α_h2o2 = - 1.08e - 6 
       α_hno3 =  1.95e - 6 
       α_w    =  6.14e - 6 
       α_cn   =  1.64e - 5 
       α_so4  =  2.48e - 6 
    end

    # Updraft velocity over the oceans [cm / s]
    # 500 cm / s is too high. Get W from the met field. (qjc, 04 / 10 / 16)
    #W = 500e + 0_fp
    w = - omega(i, j, l) / ( airden(i, j, l) * g0 ) * 100e + 0 

    # Compute H (integration time interval * air parcel velocity) [m]
    # DTCHEM is the chemistry timestep in seconds

    # Compute air parcel velocity [m / s]
    #APV = SQRT( (U(I, J, L) * U(I, J, L)) + (V(I, J, L) * V(I, J, L)) )
    #(qjc, 04 / 10 / 16)
    apv = sqrt( u(i, j, l)^2 + v(i, j, l)^2 ) + ( w^2 * 1.0e-4  )

    h   = dtchem * apv          #[m]

    sum_gas = ( α_nh3  * nh3  ) + ( α_so2  * so2  ) + ( α_h2o2 * h2o2 ) + ( α_hno3 * hno3 )

    dsvi = ( α_b * b ) + ( ( ( α_cn * cn) + ( α_w * w ) + ( α_so4 * so4 ) + sum_gas ) * h )

    # Only calculate SR when air parcel rises, in consistence with
    # Yuen et al. (1996) (qjc, 04 / 10 / 16)
    if ( w > 0.0  && c(ind_so2) > 0.0  ) 

       # additional sulfate production that can be
       # attributed to ozone [ug / m3 / timestep]
       # Don"t allow SR to be negative
       sr = max( ( dsvi - b ), 0.0  )

       # Skip further computation if SR = 0
 end
       if ( sr > 0.0  ) 

          # Convert SR from [ug / m3 / timestep] to [v / v / timestep]
          sr = sr * ( airmw / mw_so4 ) * 1.0e - 9  / airden(i, j, l)

          # Don"t produce more SO4 than SO2 available after AQCHEM_SO2
          # -  - SR is dSO4 / timestep (v / v) continue onvert
          #    to 1st order rate
          sr = min( sr, so2 / 1.0e12  ) / ( c(ind_so2) * cvf * dt )
       end
    end

    # Free pointers
    ad     = > null()
    airden = > null()
    airvol = > null()
    omega  = > null()
    u      = > null()
    v      = > null()

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: WhenExpCantBeDone
#
# #DESCRIPTION: Prevents floating point errors if exponential terms in routine
#  Het_Drop_Chem above can"t be done.  In the case of a negative XX, R should be
#  approximated as R1, instead of R2.  In other words,
#  R1 * R2 * ( XX - 1.0D0 ) / ( ( R1 * XX ) - R2 )
#  reaches different limits when XX reaches + Inf and - Inf.
#\\
#\\
# #INTERFACE:
#
  function whenexpcantbedone( r1, r2, k, dt ) result( r )
#
# #USES:
#
#
# #INPUT PARAMETERS:
#
    real(dp), intent(in) :: r1   # 1st term
    real(dp), intent(in) :: r2   # 2nd term
    real(dp), intent(in) :: k    # Rate [1 / s]
    real(dp), intent(in) :: dt   # timesetep [s]
#
# #RETURN VALUE:
#
    real(dp)             :: r    # new rate [1 / s]
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #LOCAL VARIABLES
#
    real(dp) :: diff

    dif f = r1 - r2
 end

    # R1 <  R2
    if ( diff < 0.0  ) 
       r = r1
       return
    end

    # R1 >  R2
    if ( diff > 0.0  ) 
       r = r2
       return
    end

    # R1 = = R2
    r = r1 - 1.0  / ( k * dt + ( 1.0  / r1 ) )

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: set_so2
#
# #DESCRIPTION: Subroutine SET\_SO2 is the SO2 chemistry subroutine.
#  (rjp, bmy, 11 / 26 / 02, 8 / 26 / 10) Adapted from CHEM_SO2() in SULFATE_MOD
#  (MSL - Spring 2021)
#\\
#\\
# #INTERFACE:
#
  subroutine set_so2( i,         j,          l,          input_opt, state_chm, state_diag, state_grid, state_met, size_res,  rc                                         )
#
# #USES:
#
#
# #INPUT PARAMETERS:
#
    integer,        intent(in)    :: i, j, l     # Grid box indices
    type(optinput), intent(in)    :: input_opt   # Input Options object
    type(grdstate), intent(in)    :: state_grid  # Grid State object
    type(metstate), intent(in)    :: state_met   # Meteorology State object
#
# #INPUT / OUTPUT PARAMETERS:
#
    type(chmstate), intent(inout) :: state_chm   # Chemistry State object
    type(dgnstate), intent(inout) :: state_diag  # Diagnostics State object
#
# #OUTPUT PARAMETERS:
#
    logical,        intent(out)   :: size_res    # Should we HetDropChem?
    integer,        intent(out)   :: rc          # Success or failure?
#
# #REMARKS:
#  Reaction List (by Rokjin Park)
#  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#  (1 ) SO2 production:
#       DMS + OH, DMS + NO3 (saved in CHEM_DMS)
#                                                                             .
#  (2 ) SO2 loss:
#       (a) SO2 + OH - > SO4
#       (b) SO2 - > drydep
#       (c) SO2 + H2O2 or O3 (aq) - > SO4
#                                                                             .
#  (3 ) SO2 = SO2_0 * exp( - bt) +  PSO2_DMS / bt * [1 - exp( - bt)]
#                                                                             .
#       where b is the sum of the reaction rate of SO2 + OH and the dry
#       deposition rate of SO2, PSO2_DMS is SO2 production from DMS in
#       MixingRatio / timestep.
#                                                                             .
#  If there is cloud in the gridbox (fraction = fc),  the aqueous
#  phase chemistry also takes place in cloud. The amount of SO2 oxidized
#  by H2O2 in cloud is limited by the available H2O2; the rest may be
#  oxidized due to additional chemistry, e.g, reaction with O3 or O2
#  (catalyzed by trace metal).
#
# #REVISION HISTORY:
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
    real(fp),  parameter  :: hplus_45  = 3.16227766016837953e - 5_fp  #pH = 4.5
    real(fp),  parameter  :: hplus_50  = 1.0e-5_fp  #pH = 5.0
    real(fp),  parameter  :: mindat    = 1.0e - 20_fp
#
# #LOCAL VARIABLES:
#
    # Scalars
    real(fp)              :: k0,     ki,      kk,     m,    l1
    real(fp)              :: l2,     l3,      ld,     f,    fc
    real(fp)              :: rk,     rkt,     dtchem, dt_t, tk
    real(fp)              :: f1,     rk1,     rk3,    so20, avo_over_lwc
    real(fp)              :: so2_cd, h2o20,   l2s,    l3s
    real(fp)              :: lwc,    kaqh2o2, kaqo3,  patm, rho, cnvfac
    real(fp)              :: alk,    alk1,    alk2,   so2_afterss
    real(fp)              :: alka,   alkc
    real(fp)              :: kt1,    kt2
    real(fp)              :: pso4e,  pso4f,   kt1n,    kt2n
    real(fp)              :: xx, kt1l, kt2l
    real(fp)              :: hplus,  so4nss, tnh3,   tno3,  gno3, anit
    real(fp)              :: lstot,  alkdst, alkss,  alkds, nh3, cl, tna
    real(fp)              :: sscvv,  aso4,   so2_sr, sr,    tanit
    real(fp)              :: tfa,  taa,   tdca    # (jmm, 12 / 03 / 2018)
    real(fp)              :: so2_gas,   ph2so4d_tot
    real(fp)              :: h2so4_cd,  h2so4_gas

    # (qjc, 04 / 10 / 16)
    real(fp)              :: l5, l5s, sro3, srhobr
    real(fp)              :: l5_1, l5s_1, l3_1, l3s_1, kaqo3_1
    real(fp)              :: hso3aq, so3aq
    real(fp)              :: so2_afterss0, rsiv, fupdatehobr_0
    real(fp)              :: hco3, hchobr, ko3, khobr, f_srhobr, hobr0
    real(fp)              :: tmp

    real(fp)              :: kaqo2, l4, l4s, mnii, feiii
    real(fp)              :: dust,  mn_ant,  mn_nat
    real(fp)              :: mn_tot, mn_d,    fe_d
    real(fp)              :: fe_ant, fe_nat,  fe_tot
    real(fp)              :: fe_d_ant, fe_d_nat

    real(fp)              :: l6, l6s, srhocl, l6_1, l6s_1      #XW
    real(fp)              :: fupdatehocl_0  #XW
    real(fp)              :: hchocl, khocl, f_srhocl, hocl0 #XW

    real(fp)              :: kaqhcho, kaqhms, kaqhms2, hmsc # JMM, MSL

    # Pointers
    type(spcconc), pointer :: spc(:)
    real(fp), pointer      :: ssalk(:)

    character(len = 255)     :: errmsg, thisloc

#ifdef luo_wetdep
    # For Luo et al wetdep scheme
#endif # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
 end
    # set_so2 begins here#
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    if ( id_h2o2 < 0 || id_so2 < 0  ) return

    # Initialize
    rc       = gc_success
    size_res = false
    errmsg   = ""
    thisloc  = " - > at set_so2 (in module kpp / fullchem / fullchem_sulfurchemfuncs.f90)"

    # Initialize variables
    spc                         = > state_chm%species
    ssalk                       = > state_chm%ssalk(i, j, l, :)
    state_chm%iscloud(i, j, l)    = 0.0_fp
    state_chm%phcloud(i, j, l)    = 0.0_fp
    state_chm%qlxphcloud(i, j, l) = 0.0_fp
    state_chm%hso3_aq(i, j, l)    = 1.0e-32_fp
    state_chm%so3_aq(i, j, l)     = 1.0e-32_fp
    dtchem                      = get_ts_chem()  # Timestep [s]
    is_fullchem                 = input_opt%its_a_fullchem_sim
    is_offline                  = ( ! is_fullchem )
    ld                          = 0.0_fp
    lstot                       = 0.0_fp
    rho                         = state_met%airden(i, j, l)
    cnvfac                      = 1.0e3_fp * airmw / ( rho * avo ) #mcl / cm3 - >v / v
    so20                        = spc(id_so2)%conc(i, j, l) * cnvfac
    so2_afterss                 = spc(id_so2)%conc(i, j, l) * cnvfac
    h2o20                       = spc(id_h2o2)%conc(i, j, l) * cnvfac
    kaqh2o2                     = 0.0_fp
    kaqo3                       = 0.0_fp
    kaqo3_1                     = 0.0_fp
    kaqo2                       = 0.0_fp
    k_cld                       = 0.0_fp
    hplus                       = 0.0_fp

    # Factor to convert AIRDEN from [kg air / m3] to [molec air / cm3]
    f                           = 1000.0e + 0_fp / airmw * avo * 1.0e - 6_fp

    # Meteorological data
    patm = state_met%pmid_dry( i, j, l ) / ( atm * 1.0e - 2_fp ) # Press, dry [atm]
    tk   = state_met%t(i, j, l)                                 # Temperature [K]
    fc   = state_met%cldf(i, j, l)                              # Cloud frac [1]

    # Get liquid water content [m3 H2O / m3 air] within cloud from met flds
    # Units: [kg H2O / kg air] * [kg air / m3 air] * [m3 H2O / 1e3 kg H2O]
#ifdef luo_wetdep
    # Luo et al wetdep scheme
    if ( is_qq3d ) 
       lwc = state_met%ql(i, j, l) * state_met%airden(i, j, l) * 1e-3_fp + max( 0.0_fp, state_chm%qq3d(i, j, l) * dtchem )
    else
       lwc = state_met%ql(i, j, l) * state_met%airden(i, j, l) * 1e-3_fp
    end
#else
    # Default scheme
    lwc = state_met%ql(i, j, l) * state_met%airden(i, j, l) * 1e-3_fp
#endif

    # QL can sometimes be negative, so force LWC to be positive
    lwc = max( 0.0_fp, lwc )

    # LWC is a grid - box averaged quantity. To improve the representation
    # of sulfate chemistry, we divide LWC by the cloud fraction and
    # compute sulfate chemistry based on the LWC within the cloud.  We
    # get the appropriate grid - box averaged mass of SO2 and sulfate by
    # multiplying these quantities by FC AFTER computing the aqueous
    # sulfur chemistry within the cloud. (lzh, jaf, bmy, 5 / 27 / 11)
    lwc = safediv( lwc, fc, 0.0_fp )


    # If (1) there is cloud, (2) there is SO2 present, (3) T > - 15 C, and
    # (4) liquid water content (LWC) is present (but not small enough to
    # make divisions blow up),  compute sulfate production in cloud.
    if ( ( fc          > 1.0e-4_fp  ) && ( so2_afterss > mindat     ) && #ifdef luo_wetdep
         ( tk          > 237.0_fp   ) && #else
         ( tk          > 258.0_fp   ) && #endif
         ( lwc         > 1.0e-20_fp ) ) 

       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       # NOTE...Sulfate production from aquatic reactions of SO2
       # with H2O2O3 is computed here and followings are
       # approximations or method used for analytical (integral)
       # solution of these computations.
       #
       # 1) with H2O2(aq)
       #      [HSO3 - ] + [H + ] + [H2O2(aq)] = > [SO4 = ]     (rxn)
       #      d[SO4 = ] / dt = k[H + ][HSO3 - ][H2O2(aq)] (M / s) (rate)
       #
       # we can re#FIXME write k[H + ][HSO3 - ] as K1 pSO2 hSO2,
       # where pSO2 is equilibrium vapor pressure of SO2(g)
       # in atm, and hSO2 is henry"s law constant for SO2
       #
       # Therefore, rate can be written as
       #
       #       k * K1 * pSO2 * hSO2 * pH2O2 * hH2O2,
       #
       # where pH2O2 is equilibrium vapor pressure of H2O2(g),
       # and hH2O2 is henry"s law constant for H2O2. Detailed
       # values are given in AQSET_SO2 routine.
       #
       # Let us define a fraction of gas phase of A species
       # in equilibrium with aqueous phase as
       #
       #        xA  = 1 / (1 + f),
       #
       # where  f   = hA * R * T * LWC,
       #        hA  = Henry"s constant,
       #        R   = gas constant,
       #        T   = temperature in kelvin,
       #        LWC = liquid water content [m3 / m3]
       #
       # As a result, the rate would become:
       #
       #    d[SO4 = ]
       # -  -  -  -  -  -  - = k K1 hSO2 hH2O2 xSO2 xH2O2 P P [SO2][H2O2]
       #      dt
       #      ^       ^                            ^   ^    ^
       #      |       |____________________________|   |    |
       #
       #   mole / l / s               mole / l / s            v / v  v / v
       #
       #
       # And we multiply rate by (LWC * R * T / P) in order to
       # convert unit from mole / l / s to v / v / s
       #
       # Finally we come to
       #
       #    d[SO4 = ]
       # -  -  -  -  -  -  - = KaqH2O2 [SO2][H2O2],
       #      dt
       #
       # where
       #
       #   KaqH2O2 = k K1 hSO2 hH2O2 xSO2 xH2O2 P LWC R T,
       #
       # this new rate corresponds to a typical second order
       # reaction of which analytical (integral) solution is
       #
       #   X  = A0 B0 ( exp[(A0 - B0) Ka t] - 1 )
       #      / ( A0 exp[(A0 - B0) Ka t] - B0 )
       #
       # inserting variables into solution  we get
       # [SO4 = ] =  [SO2][H2O2](exp[([SO2] - [H2O2]) KaqH2O2 t] - 1 )
       #        / ( [SO2] exp[([SO2] - [H2O2]) KaqH2O2 t] - [H2O2] )
       #
       # Note...Exactly same method can be applied to O3 reaction
       # in aqueous phase with different rate constants.
       # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

       # Get concentrations for cloud pH calculation (bec, 12 / 23 / 11)

       # <<>><<>><<>><<>>
       # HAVE TO DO SOME UNIT CONVERSION HERE. THIS ROUTINE IS CALLED
       # WITHIN FLEXCHEM_MOD WHERE SPECIES ARE IN MOLEC / CM3
       # <<>><<>><<>><<>>

       # Get sulfate concentration and convert from [v / v] to
       # [moles / liter]
       # Use a cloud scavenging ratio of 0.7

       so4nss = 1.0e + 3 * ( spc(id_so4)%conc(i, j, l) * 0.7e + 0_fp + spc(id_so4s)%conc(i, j, l) ) / ( lwc * avo ) # mcl / cm3 - > mol / L

       # Get HMS cloud concentration and convert from [v / v] to
       # [moles / liter] (jmm, 06 / 13 / 2018)
       # Use a cloud scavenging ratio of 0.7
       # assume nonvolatile like sulfate for realistic cloud pH
       hmsc = 0.0_fp
       if ( is_fullchem && id_hms > 0 ) 
          hmsc = 1.0e + 3 * spc(id_hms)%conc(i, j, l) * 0.7_fp / ( lwc * avo ) # mcl / cm3 - > mol / L
       end

       # Get total ammonia (NH3 + NH4 + ) concentration [v / v]
       # Use a cloud scavenging ratio of 0.7 for NH4 +
       tnh3 = ( ( spc(id_nh4)%conc(i, j, l) * 0.7e + 0_fp ) + spc(id_nh3)%conc(i, j, l) ) * cnvfac

       # Get total chloride (SALACL + HCL) concentration [v / v]
       # Use a cloud scavenging ratio of 0.7
       cl = ( spc(id_salacl)%conc(i, j, l) * 0.7e + 0_fp ) + spc(id_salccl)%conc(i, j, l)
       cl = ( cl + spc(id_hcl)%conc(i, j, l) ) * cnvfac

       # Get total formic acid concentration [v / v]
       # jmm (12 / 3 / 18)
       # no cloud scavenging because gases?
       tfa = spc(id_hcooh)%conc(i, j, l) * cnvfac

       # Get total acetic acid concentration [v / v]
       # jmm (12 / 3 / 18)
       # no cloud scavenging b / c gases?
       taa = spc(id_acta)%conc(i, j, l) * cnvfac

       # Get total sea salt NVC concentration expressed as NA + equivalents
       # and convert from [MND] to [moles / liter]
       # NVC is calculated to balance initial Cl - + alkalinity in
       # seas salt. Note that we should not consider SO4ss here.
       # Use a cloud scavenging ratio of 0.7 for fine aerosols
       tna      = 1.0e3_fp * ( spc(id_sala)%conc(i, j, l) * 0.7e + 0_fp + spc(id_salc)%conc(i, j, l) ) * ( 31.6e + 0_fp * 0.359e + 0_fp / 23.0e + 0_fp ) / ( lwc * avo ) # mcl / cm3 - > mol / L

       # Get total dust cation concentration [mol / L]
       # Use a cloud scavenging ratio of 1 for dust
       # to be consistent for how it was calculated for
       # metal catalyzed SO2 oxidation
       # Use asumption of dust being 3% soluble Ca2 + and
       # 0.6% soluble Mg2 + by mass (Fairlie et al., 2010)
       #
       # Dust treated at non - volatile cation and charge applied in
       # pH calculation
       #
       # Move dust calculation from SO2 Metal catalzyed oxidation
       # up here becasue needed for cloud pH
       # jmm (12 / 3 / 18)
       #
       # Get dust concentrations [MND - > ng / m3]

       dust = ( spc(id_dst1)%conc(i, j, l) * 0.7_fp + spc(id_dst2)%conc(i, j, l) + spc(id_dst3)%conc(i, j, l) + spc(id_dst4)%conc(i, j, l) ) * 1.0e + 15_fp * state_chm%spcdata(id_dst1)%info%mw_g / avo

       # Conversion from dust mass to Ca2 + and Mg2 + mol:
       #     0.071 * (1 / 40.08) + 0.011 * (1 / 24.31) = 2.22e - 3
       #     (Engelbrecht et al., 2016)
       #     1e-12_fp from m3 - >Lng - >g
       tdca     = dust * 2.22e - 15_fp / lwc

       # Get total nitrate (HNO3 + NIT) concentrations [v / v]
       # Use a cloud scavenging ratio of 0.7 for NIT
       tno3 = ( spc(id_hno3)%conc(i, j, l) + ( spc(id_nit)%conc(i, j, l)  * 0.7e + 0_fp ) + spc(id_nits)%conc(i, j, l) ) * cnvfac
       gno3 = spc(id_hno3)%conc(i, j, l) * cnvfac # For FaheyPandis decision algorithm

       # Calculate cloud pH
       get_hplus( so4nss, hmsc, tnh3, tno3,  so2_afterss,   cl, tna, tdca, tfa,    taa,  tk,   patm,  lwc, hplus_45, hplus  )

       # Store the cloud pH quantities
       state_chm%iscloud(i, j, l)    =  1.0_fp
       state_chm%phcloud(i, j, l)    = - 1.0_fp * log10(hplus)
       state_chm%qlxphcloud(i, j, l) = state_chm%phcloud(i, j, l) * state_met%ql(i, j, l)


       feiii = 0.0_fp
       mnii  = 0.0_fp
       if ( input_opt%lmetalcatso2 ) 

          # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
          # Metal catalyzed oxidation of SO2 pathway
          # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 

          # Get dust concentrations [v / v - > ng / m3]
#ifdef tomas
          # TOMAS uses its own dust tracers and does not
          # carry DST1 - 4.0  Set DUST to zero here. (mps, 2 / 2 / 18)
          dust = 0e + 0_fp
#endif
          # Calculate Fe and Mn natural [ng m - 3]
          # Assume that Fe is 3.5% of total dust mass based on
          # Taylor and McLennan [1985]
          fe_nat = dust * 35e - 3_fp
          # and Mn is 50 times less than Fe based on Desbouefs et al.[2005]
          mn_nat = fe_nat / 50e + 0_fp

          # Anthropogenic Fe concentrations [mcl / cm3 - > ng / m3]
          if ( id_pfe > 0 ) 
                fe_ant = spc(id_pfe)%conc(i, j, l) * cnvfac * 1.0e + 12_fp * state_met%ad(i, j, l) / ( airmw / state_chm%spcdata(id_pfe)%info%mw_g ) / state_met%airvol(i, j, l)
#             Fe_ant = Spc(id_pFe)%Conc(I, J, L) * 1.0e + 15_fp * #                  State_Chm%SpcData(id_DST1)%Info%MW_g / AVO
          else
             fe_ant = 0e + 0_fp
          end

          # Calculate Mn anthropogenic [ng m - 3]
          # assume anthropogenic Mn is 1 / 30 times anthropogenic Fe
          mn_ant = fe_ant / 10e + 0_fp

          # Calculate total Mn and Fe [ng m - 3]
          mn_tot = mn_ant + mn_nat
          fe_tot = fe_ant + fe_nat

          # Convert Mn and Fe [ng m - 3] to [mole l - 1]

          # Assume that 50% of Mn is dissolved [Spokes et al., 1994]
          # Hardcoded MW for Mn
          if ( lwc > 0e + 0_fp ) 
             # Units: ng / m3 * (g / ng) / (g / mol) / (m3 H2O / m3 air) * (m3 / L)
             mn_d = mn_tot * 1e-9_fp / 54.94e + 0_fp / lwc * 1e-3_fp
             mn_d = mn_d * 0.5e + 0_fp
          else
             mn_d = 0e + 0_fp
          end

          # Solubility of Fe is 10% for anthropogenic, and 1% for dust
          if ( lwc > 0e + 0_fp ) 
             fe_d_ant = fe_ant * 1e-9_fp / state_chm%spcdata(id_pfe)%info%mw_g / lwc * 1e-3_fp
             fe_d_nat = fe_nat * 1e-9_fp / state_chm%spcdata(id_pfe)%info%mw_g / lwc * 1e-3_fp
             fe_d     = fe_d_ant * 0.1e + 0_fp + fe_d_nat * 0.01e + 0_fp
          else
             fe_d     = 0e + 0_fp
          end

          # Impose a dependence of Fe speciation on sunlight
          if ( state_met%suncos(i, j) > 0e + 0_fp ) 
             # Assume 10% of dissolved Fe is in Fe(III)
             #oxidation state during the daytime
             feiii = fe_d * 0.1e + 0_fp
          else
             # Assume 90% of dissolved Fe is in Fe(III)
             # oxidation state during the nighttime
             feiii = fe_d * 0.9e + 0_fp
          end

          # Assume that dissolved Mn is in Mn(II) oxidation state all of
          # the time
          mnii = mn_d
       end

       # Compute aqueous rxn rates for SO2
       aqchem_so2( i       = i, j       = j, l       = l, lwc     = lwc, t       = tk, p       = patm, so2     = so2_afterss, h2o2    = spc(id_h2o2)%conc(i, j, l) * cnvfac, o3      = spc(id_o3)%conc(i, j, l)   * cnvfac, hcho    = spc(id_ch2o)%conc(i, j, l) * cnvfac, hplus   = hplus, mnii    = mnii, feiii   = feiii, kaqh2o2 = kaqh2o2, kaqo3   = kaqo3, kaqo3_1 = kaqo3_1, kaqo2   = kaqo2, hso3aq  = hso3aq, so3aq   = so3aq, kaqhcho = kaqhcho, kaqhms  = kaqhms, kaqhms2 = kaqhms2                                   )


       k_cld(1) = kaqh2o2 * fc * cnvfac   # v / v / s -  - > cm3 / mcl / s
       k_cld(2) = kaqo3   * fc * cnvfac   # v / v / s -  - > cm3 / mcl / s
       k_cld(3) = kaqo2   * fc           # 1 / s
       # vvvvvv Hold off using CloudHet2R until after initial S - chem benchmark
       # -  - MSL
       #K_CLD(1) = CloudHet2R( Spc(id_SO2)%Conc(I, J, L), #                       Spc(id_H2O2)%Conc(I, J, L), FC, KaqH2O2 * CNVFAC )
       #K_CLD(2) = CloudHet2R( Spc(id_SO2)%Conc(I, J, L), #                       Spc(id_O3)%Conc(I, J, L),   FC, KaqO3   * CNVFAC )
       #K_CLD(3) computed below

       # HMS reaction rates (skip if HMS isn"t defined)
       if ( is_fullchem && id_hms > 0 ) 
          k_cld(4) = kaqhcho * fc * cnvfac
          k_cld(5) = kaqhms  * fc
          k_cld(6) = kaqhms2 * fc
          # Leave comments here (bmy, 18 Jan 2022)
          #          CloudHet2R( Spc(id_HMS)%Conc(I, J, L), #                      Spc(id_CH2O)%Conc(I, J, L), FC, KaqHCHO * CNVFAC )
          #          cloudhet1r( fc, kaqhms ) # KaqHMS is pseudo - 1st order
          #          CloudHet2R( Spc(id_HMS)%Conc(I, J, L), #                      Spc(id_OH)%Conc(I, J, L), FC, ...)       ENDIF
       end

#ifdef tomas
       #%%%%%%%%%%%%%%%%% BUG FIX FOR TOMAS %%%%%%%%%%%%%%%%%%%%%%%
       # NOTE: TOMAS uses its own dust tracers and does not
       # carry ALKdst.  Set ALKdst to zero here. (bmy, 1 / 28 / 14)
       alkdst = 0e + 0_fp
#else

       # For other simulations, Sum up the contributions from
       # DST1 thru DST4 tracers into ALKdst. (bmy, 1 / 28 / 14)
       # mcl / cm3 - > ug / m3
       alkdst = ( spc(id_dst1)%conc(i, j, l) + spc(id_dst2)%conc(i, j, l) + spc(id_dst3)%conc(i, j, l) + spc(id_dst4)%conc(i, j, l) ) * cnvfac * 1.0e + 9_fp * state_met%ad(i, j, l) / ( airmw / state_chm%spcdata(id_dst1)%info%mw_g ) / state_met%airvol(i, j, l)
#endif

       # mcl / cm3 - > ug / m3
       alkss  = ( spc(id_sala)%conc(i, j, l) + spc(id_salc)%conc(i, j, l) ) * cnvfac * 1.0e + 9_fp * state_met%ad(i, j, l) / ( airmw / state_chm%spcdata(id_sala)%info%mw_g ) / state_met%airvol(i, j, l)

       alkds = alkdst + alkss

       # Get NH3 concentrations (v / v)
       nh3 = spc(id_nh3)%conc(i, j, l) * cnvfac

       # Initialize
       size_res = false

       # Fahey and Seinfeld decision algorithm
       # NOTE: This is ugly, needs refactoring.  For now, just added
       # whitespace to improve readability (bmy, 01 Oct 2021)
       if ( h2o20 > so2_afterss + 1e-9_fp ) 
          size_res = false

       elseif ( lwc < 0.1e-6_fp )  #10^ - 6 coversion from g / m3 -  - > m3 / m3
          size_res = true

       elseif ( gno3 > nh3 ) 

          if ( so2_afterss > = 5.0e - 9_fp && h2o20       > = so20              ) size_res = false
 end

          if ( lwc         > = 0.3e-6_fp && so2_afterss > = 3.0e - 9_fp && h2o20       > = so2_afterss       ) size_res = false
 end

          if ( alkds       > = 5.0e + 0_fp && lwc         > = 0.5e-6_fp && h2o20       > = so2_afterss       ) size_res = false
 end

          if ( lwc         > = 0.1e-6_fp && gno3        < = (nh3 + 2.0e - 9_fp)  ) size_res = false
 end

       elseif ( lwc > = 0.5e-6_fp ) 
 end

          if ( h2o20       > = ( 0.9_fp * so2_afterss )         ) size_res = false
 end

          if ( nh3         < = 1.0e - 9_fp && alkds       > = 5.0e + 0_fp && so2_afterss < = 10.0e - 9_fp         ) size_res = false
 end

       elseif ( lwc > = 0.3e-6_fp ) 
 end

          if ( nh3         > = (gno3 + 5.0e - 9_fp) && so2_afterss < = 10.0e - 9_fp         ) size_res = false
 end

          if ( gno3        < = 1.0e - 9_fp && nh3         > = (gno3 + 2.0e - 9_fp) ) size_res = false
 end

          if ( gno3        < = 7.0e - 9_fp && nh3         > = (gno3 + 3.0e - 9_fp) ) size_res = false
 end

          if ( alkds       > = 3.0e + 0_fp && nh3         < = 10e - 9_fp && so2_afterss < = 5e-9_fp           ) size_res = false
 end

          if ( alkds       > = 5.0e + 0_fp && nh3         < = 10.0e - 9_fp && so2_afterss < = 5.0e - 9_fp          ) size_res = false
 end

          if ( so2_afterss > = 1.5e-9_fp && h2o20       > = so2_afterss       ) size_res = false
 end

          if ( nh3         < = 12.0e - 9_fp && alkds       > = 10.0e + 0_fp         ) size_res = false
 end

          if ( nh3         < = 1.0e - 9_fp && alkds       > = 4.0e + 0_fp && so2_afterss < = 10.0e - 9_fp         ) size_res = false
 end

          if ( nh3         < = 5.0e - 9_fp && alkds       > = 6.0e + 0_fp && so2_afterss < = 10.0e - 9_fp         ) size_res = false
 end

          if ( nh3         < = 7.0e - 9_fp && alkds       > - 8.0e + 0_fp && so2_afterss < = 10.0e - 9_fp         ) size_res = false
 end

       elseif ( lwc > = 0.1e-6_fp ) 
 end

          if ( nh3         < = 1.0e - 9_fp && alkds       > = 5.0e + 0_fp          ) size_res = false
 end

          if ( nh3         < = 5.0e - 9_fp && alkds       > = 10.0e + 0_fp         ) size_res = false
 end

          if ( gno3        < = 1.0e - 9_fp && nh3         > = (gno3 + 2.0e - 9_fp) && so2_afterss < = 7.0e - 9_fp          ) size_res = false
 end

          if ( gno3        < = 1.0e - 9_fp && nh3         > = (gno3 + 2.0e - 9_fp) && alkds       > = 2.0e + 0_fp          ) size_res = false
 end

          if ( gno3        < = 3.0e - 9_fp && nh3         > = (gno3 + 4.0e - 9_fp) ) size_res = false
 end

          if ( gno3        < = 7.0e - 9_fp && nh3         > = (gno3 + 3.0e - 9_fp) && so2_afterss < = 5.0e - 9_fp          ) size_res = false
 end

          if ( gno3        < = 7.0e - 9_fp && nh3         > = (gno3 + 3.0e - 9_fp) && alkds       > = 4.0e + 0_fp && so2_afterss < = 9.0e - 9_fp          ) size_res = false
 end

          if ( alkds       > = 3.0e + 0_fp && nh3         < = 3.0e - 9_fp && so2_afterss < = 4.0e - 9_fp          ) size_res = false
 end

          if ( alkds       > = 5.0e + 0_fp && so2_afterss < = 5.0e - 9_fp && nh3         < = 7.0e - 9_fp          ) size_res = false
 end

          if ( nh3         > = (gno3 + 2.0e - 9_fp) && so2_afterss < = 5.0e - 9_fp          ) size_res = false
 end

          if ( nh3         > = (gno3 + 4.0e - 9_fp) && so2_afterss < = 10.0e - 9_fp         ) size_res = false
 end

          if ( alkds       > = 2.0e + 0_fp && nh3         < = 10.0e - 9_fp && h2o20       > = so2_afterss       ) size_res = false
 end

          if ( nh3         < = 1.0e - 9_fp && so2_afterss > = 3.0e - 9_fp && h2o20       > = so2_afterss       ) size_res = false
 end

       else

          size_res = true

       end


       state_chm%hso3_aq(i, j, l) = hso3aq
       state_chm%so3_aq(i, j, l)  = so3aq

    end
#>>    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#>>    # HISTORY (aka netCDF diagnostics)
#>>    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#>>
#>>    # P(SO4) from gas - phase oxidation [kg S / s]
#>>    IF ( State_Diag%Archive_ProdSO4fromGasPhase ) THEN
#>>       State_Diag%ProdSO4fromGasPhase(I, J, L) = #>>            ( L1  * State_Met%AD(I, J, L) / TCVV_S ) / DTCHEM
#>>    ENDIF
#>>
#>>    # P(SO4) from aqueous - phase oxidation with H2O2 [kg S / s]
#>>    IF ( State_Diag%Archive_ProdSO4fromH2O2inCloud ) THEN
#>>       State_Diag%ProdSO4fromH2O2inCloud(I, J, L) = #>>            ( L2S * State_Met%AD(I, J, L) / TCVV_S ) / DTCHEM
#>>    ENDIF
#>>
#>>    # P(SO4) from aqueous - phase oxidation with O3 [kg S / s]
#>>    IF ( State_Diag%Archive_ProdSO4fromO3InCloud ) THEN
#>>       State_Diag%ProdSO4fromO3InCloud(I, J, L) = #>>            ( ( L3S + SRo3 ) * State_Met%AD(I, J, L) / TCVV_S ) / DTCHEM
#>>    ENDIF
#>>
#>>    # P(SO4) from aqueous - phase oxidation with HOBr [kg S / s]
#>>    IF ( State_Diag%Archive_ProdSO4fromHOBrInCloud ) THEN
#>>       State_Diag%ProdSO4fromHOBrInCloud(I, J, L) = #>>            ( ( L5S + SRhobr ) * State_Met%AD(I, J, L) / TCVV_S ) / DTCHEM
#>>    ENDIF
#>>
#>>    # P(SO4) from aqueous - phase oxidation with O2 metal - catalyzed
#>>    # [kg S / s]
#>>    IF ( State_Diag%Archive_ProdSO4fromO2InCloudMetal ) THEN
#>>       State_Diag%ProdSO4fromO2InCloudMetal(I, J, L) = #>>            ( L4S * State_Met%AD(I, J, L) / TCVV_S ) / DTCHEM
#>>    ENDIF
#>>
#>>    # P(SO4) from O3 in sea salt aerosol [kg S / s]
#>>    IF ( State_Diag%Archive_ProdSO4fromO3inSeaSalt ) THEN
#>>       State_Diag%ProdSO4fromO3inSeaSalt(I, J, L) = #>>            ( ( PSO4E + PSO4F ) * State_Met%AD(I, J, L) / TCVV_S ) / DTCHEM
#>>    ENDIF
#>>
#>>    # P(SO4) by SRo3 [kg S / s]
#>>    IF ( State_Diag%Archive_ProdSO4fromSRO3 ) THEN
#>>       State_Diag%ProdSO4fromSRO3(I, J, L) = #>>            ( SRo3 * State_Met%AD(I, J, L) / TCVV_S ) / DTCHEM
#>>    ENDIF
#>>
#>>    # P(SO4) by SRhobr [kg S / s]
#>>    IF ( State_Diag%Archive_ProdSO4fromSRHOBr ) THEN
#>>       State_Diag%ProdSO4fromSRHOBr(I, J, L) = #>>            ( SRhobr * State_Met%AD(I, J, L) / TCVV_S ) / DTCHEM
#>>    ENDIF
#>>
#>>    # P(SO4) by o3s [kg S / s]
#>>    IF ( State_Diag%Archive_ProdSO4fromO3s ) THEN
#>>       State_Diag%ProdSO4fromO3s(I, J, L) = #>>            ( L3S_1 * State_Met%AD(I, J, L) / TCVV_S ) / DTCHEM
#>>    ENDIF
#>>

    # Free pointers
    spc     = > null()
    ssalk   = > null()

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: get_hplus
#
# #DESCRIPTION: Subroutine GET\_HPLUS computes H + concentrations in cloud
#  liquid water for pH dependent cloud chemistry. (bec, 4 / 11 / 11)
#\\
#\\
# #INTERFACE:
#
    subroutine get_hplus( so4nss, hmsc, tnh3, tno3, so2, cl, tna, tdca, tfa, taa,  t, pres, lwc,  ihplus, hplus )
#
# #USES:
#
#
# #INPUT PARAMETERS:
#
    real(fp),  intent(in)    :: so4nss # Total nss sulfate mixing ratio [M]
    real(fp),  intent(in)    :: hmsc # Total HMS mixing ratio [M]
    real(fp),  intent(in)    :: tno3   # Total nitrate (gas + particulate) mixing
                                       # ratio [v / v]
    real(fp),  intent(in)    :: tnh3   # NH3 mixing ratio [v / v]
    real(fp),  intent(in)    :: so2    # SO2 mixing ratio [v / v]
    real(fp),  intent(in)    :: cl     # Total chloride (gas + particulate) mixing
    real(fp),  intent(in)    :: tna    # Sodium (particulate) [v / v]
    real(fp),  intent(in)    :: tdca   # total ca2 + and mg2 + mixing ratio [m] # jmm 12 / 3 / 18
    real(fp),  intent(in)    :: taa    # acetic acid mixing ratio [v / v] # jmm 12 / 3 / 18
    real(fp),  intent(in)    :: tfa    # formic acid mixing ratio [v / v] # jmm 12 / 3 / 18
    real(fp),  intent(in)    :: t      # Temperature [K]
    real(fp),  intent(in)    :: pres   # Dry air partial ressure [atm]
    real(fp),  intent(in)    :: lwc    # Cloud liquid water content [m3 / m3]
    real(fp),  intent(in)    :: ihplus # Initial [H + ] [M]
#
# #OUTPUT PARAMETERS:
#
    real(fp),  intent(out)   :: hplus  # Calculated [H + ] [M]
# #REMARKS:
#  Calculation:
#  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#  Solve the following electroneutrality equation:
#  [h + ] = 2[so4 -  - ] + [cl - ] + [oh - ] + [hco3 - ] + 2[co3 -  - ] + [hso3 - ] + 2[so3 -  - ] + #
#          [NO3 - ] + [HCOO - ] + [CH3COO - ] - [Na] - 2[Ca] - [NH4]
#  Uses Newton"s method to solve the equation:
#     x_1 = x_0 - f(x_0) / f"(x_0)
#     iterate until converge
#
#  Let concentrations of [HCO3], [CO3], [HSO3], [SO3], [NO3] and [NH4] evolve
#  according to Henry"s law equilibrium.
#
#  To add new species:
# - Add species not affected by HPLUS to the "D" term
# - Add species that disassociate once using the kHNO3 and dkHNO3
#    functions
#      as a template
# - Add species that disassociate twice using the kSO21 and dkSO21
#    functions
#      as a template for the single charged ion and kSO22 and dkSO22
#      functions for
#      the double charged ion

#  Assume [S(VI)] = [SO4]nss (this applies for pH > 3)
#
# #REVISION HISTORY:
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
    # Water dissociation constants
    real(fp),  parameter   :: kw   = 1.0e-14_fp
    real(fp),  parameter   :: dhrkw = - 6710.0e + 0_fp
    real(fp),  parameter   :: minval = 0.01
#
# #LOCAL VARIABLES:
#
    real(fp)               :: d, kw_t, iph, newph, nhplus
    real(fp)               :: fhco3, fco3
    real(fp)               :: fhso3, fso3
    real(fp)               :: fhno3, fnh4, fhcl
    real(fp)               :: dhco3, dco3
    real(fp)               :: dhso3, dso3
    real(fp)               :: dhno3, dnh4, dhcl
    real(fp)               :: faa, ffa, daa, dfa
    real(fp)               :: f, df, nnhplus, fca, dca

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # get_hplus begins here#
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # Initial pH guess
    iph = - log10(ihplus)

    # Non - volatile aerosol concentration [M]
    # For now sulfate is the only non - volatile species
    d = (2.0e + 0_fp * so4nss) - tna - (2.0e + 0_fp * tdca) + (1.0e + 0_fp * hmsc)

    # Temperature dependent water equilibrium constant
    kw_t = kw * exp(dhrkw * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

    # Initialize
    newph   = 0.0
    count = 0

    do #FIXME while ( abs(iph - newph) > minval )

       count = count + 1

       if ( count == 1 ) 
          iph = iph
       else
          iph = newph
       end

       nhplus = 10.0e + 0_fp^( - iph)

       # Get f(x) terms
       fhco3  = kco21 ( pres, t, lwc, nhplus )

       fco3 = kco22 ( pres, t, lwc, nhplus )

       fhso3  = kso21 ( pres, t, lwc, nhplus, so2 )

       fso3 = kso22 ( pres, t, lwc, nhplus, so2 )

       fhno3 = khno3 ( pres, t, lwc, nhplus, tno3 )

       fnh4  = knh3  ( pres, t, lwc, nhplus, tnh3, kw_t )

       # include HCl in cloud pH calculations, xnw 10 / 17 / 17
       fhcl  = khcl  ( pres, t, lwc, nhplus, cl  )

       ffa   = kfa   ( pres, t, lwc, nhplus, tfa ) # jmm 12 / 3 / 18

       faa   = kaa   ( pres, t, lwc, nhplus, taa ) # jmm 12 / 3 / 18

       # Get f"(x) terms
       dhco3  = dkco21 ( pres, t, lwc, nhplus )

       dco3 = dkco22 ( pres, t, lwc, nhplus )

       dhso3  = dkso21 ( pres, t, lwc, nhplus, so2 )

       dso3 = dkso22 ( pres, t, lwc, nhplus, so2 )

       dhno3 = dkhno3 ( pres, t, lwc, nhplus, tno3 )

       dnh4  = dknh3  ( pres, t, lwc, nhplus, tnh3, kw_t )

       dhcl = dkhcl ( pres, t, lwc, nhplus, cl )

       dfa   = dkfa   ( pres, t, lwc, nhplus, tfa ) # jmm 12 / 3 / 18

       daa   = dkaa   ( pres, t, lwc, nhplus, taa ) # jmm 12 / 3 / 18
       # Calculate [Ca2 + ] in equilibrium with CaCO3(s)
       caco3_precip ( pres, t, nhplus, fca, dca )

       # if [Ca2 + ] in equilibrium with CacO3(s) is greater than total [Ca2 + ]
       #  all Ca is dissolved else [Ca2 + ] varies with [H + ]
       if ( fca >= tdca ) 
          # Non - volatile aerosol concentration [M]
          d = (2.0e + 0_fp * so4nss) - (tna + 2.0e + 0_fp * tdca)

          # Define f(x)
          f = d - nhplus + kw / nhplus + fhco3 + 2.0e + 0_fp * fco3 + fhso3 + 2.0e + 0_fp * fso3 + fhno3 - fnh4 + fhcl + ffa + faa

          # Define f"(x)
          df = - 1.0d0 - kw / nhplus / nhplus + dhco3 + 2.0e + 0_fp * dco3 + dhso3 + 2.0e + 0_fp * dso3 + dhno3 - dnh4 + dhcl + dfa + daa

       else
          # Non - volatile aerosol concentration [M]
          d = (2.0e + 0_fp * so4nss) - tna

          # Define f(x)
          f = d - nhplus + kw / nhplus + fhco3 + 2.0e + 0_fp * fco3 + fhso3 + 2.0e + 0_fp * fso3 + fhno3 - fnh4 + fhcl + ffa + faa - 2.0e + 0_fp * fca
          # Define f"(x)
          df = - 1.0d0 - kw / nhplus / nhplus + dhco3 + 2.0e + 0_fp * dco3 + dhso3 + 2.0e + 0_fp * dso3 + dhno3 - dnh4 + dhcl + dfa + daa - 2.0e + 0_fp * dca
       end

       # Apply Newton"s method
       nnhplus = nhplus - f / df

       # Set minimum [H + ] = 1.0d - 14 (pH = 14)
       nnhplus = max(nnhplus, 1.0e-14_fp)

       # Set maximum [H + ] = 1.0d - 1 (pH = 1)
       nnhplus = min(nnhplus, 1.0e-1_fp)

       # If solution does not converge after 50 iterations
       # average last 2 pH calculations
       if (count > 50) 
          newph = (( - log10(nnhplus)) + ( - log10(nhplus))) / 2.0e + 0_fp

          if (it_is_nan( newph )) 
             #FIXME write(6, * ) "newph = ", newph
             #FIXME write(6, * ) "nnhplus = ", nnhplus
             #FIXME write(6, * ) "nhplus = ", nhplus
             geos_chem_stop
          end

          exit
       else
          newph = - log10(nnhplus)

          if (it_is_nan( newph )) 
             #FIXME write(6, * ) "newph = ", newph
             #FIXME write(6, * ) "nnhplus = ", nnhplus
             geos_chem_stop
          end

       end

    end

    hplus = 10.0e + 0_fp^( - newph)

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: kCO21
#
# #DESCRIPTION: Function kCO21
#\\
#\\
# #INTERFACE:
#
  function kco21 ( p, t, lwc, hplus ) result ( kco2p )
#
# #INPUT PARAMETERS:
#
    real(fp),  intent(in) :: t, p, lwc, hplus
#
# #OUTPUT PARAMETERS:
#
    real(fp)              :: kco2p, kco2p2
#
# #REVISION HISTORY:
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
    # CO2 dissociation constants
    real(fp),  parameter  :: kc1 = 4.3e-7
    real(fp),  parameter  :: kc2 = 4.68e - 11
    real(fp),  parameter  :: dhrkc1 = - 1000.0
    real(fp),  parameter  :: dhrkc2 = - 1760.0
    real(fp),  parameter  :: hco2 = 3.4e-2
    real(fp),  parameter  :: dhco2 = 2.44e + 3_fp
    # CO2 concentration [v / v]
    real(fp),  parameter  :: co2 = 390.0e-6_fp
#
# #LOCAL VARIABLES:
#
    real(fp)              :: hco2_t, kc1_t, kc2_t
    real(fp)              :: hco2eff, xco2, pco2

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # kco21 begins here#
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    #CO2 dissolution constants
    hco2_t  = hco2 * exp(dhco2 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
    kc1_t   = kc1 * exp(dhrkc1 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
    kc2_t   = kc2 * exp(dhrkc2 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

    #CO2 dissolution
    hco2eff = hco2_t * (1.0e + 0_fp + (kc1_t / hplus) + ((kc1_t * kc2_t) / (hplus * hplus)))
    xco2    = 1.0e + 0_fp / ( 1.0e + 0_fp + ( hco2eff * 0.08205e + 0_fp * t * lwc ) )
    pco2    = co2 * p * xco2

    kco2p  = hco2_t / hplus * kc1_t * pco2

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: dkCO21
#
# #DESCRIPTION: Function dkCO21
#\\
#\\
# #INTERFACE:
#
      function dkco21 ( p, t, lwc, hplus ) result ( kco2p )
#
# #INPUT PARAMETERS:
#
      real(fp),  intent(in) :: t, p, lwc, hplus
#
# #OUTPUT PARAMETERS:
#
      real(fp)              :: kco2p, kco2p2
#
# #REVISION HISTORY:
#  25 Jan 2012 - M. Payer - Added ProTeX headers
#  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
#  pressure
#  22 Mar 2017 - M. Sulprizio - Dhco2 value is from Table 7.3 of Seinfeld
#  and
#                              Pandis (2006, pp 289) and should be
#                              positive for
#                              consistency with the way it is used here.
#                              Also,
#                              the value of R, in units of kcal mol - 1
#                              K - 1, is
#                              1.986x10^ - 3, not 0.04.0 (V. Shah)
#  15 Feb 2019 - J. Moch - updated function to make output
#  derivative of [HCO3 - ]
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
      # CO2 dissociation constants
      real(fp),  parameter  :: kc1 = 4.3e-7
      real(fp),  parameter  :: kc2 = 4.68e - 11
      real(fp),  parameter  :: dhrkc1 = - 1000.0
      real(fp),  parameter  :: dhrkc2 = - 1760.0
      real(fp),  parameter  :: hco2 = 3.4e-2
      real(fp),  parameter  :: dhco2 = 2.44e + 3_fp
      # CO2 concentration [v / v]
      real(fp),  parameter  :: co2 = 390.0e-6_fp
#
# #LOCAL VARIABLES:
#
      real(fp)              :: hco2_t, kc1_t, kc2_t

# #REMARKS:

      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      # dkco21 begins here#
      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

      #CO2 dissolution constants
      hco2_t = hco2 * exp(dhco2 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
      kc1_t = kc1 * exp(dhrkc1 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
      kc2_t = kc2 * exp(dhrkc2 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

      #CO2 dissolution

      kco2p  = kc1_t * hco2_t * co2 * p * ( kc1_t * kc2_t * hco2_t * 0.08205e + 0_fp * t * lwc - hco2_t * 0.08205e + 0_fp * t * lwc * hplus * hplus - hplus * hplus) / (kc1_t * kc2_t * hco2_t * 0.08205e + 0_fp * t * lwc + kc1_t * hco2_t * 0.08205e + 0_fp * t * lwc * hplus + hco2_t * 0.08205e + 0_fp * t * lwc * hplus * hplus + hplus * hplus)^2

      end
#EOC

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: kCO22
#
# #DESCRIPTION: Function kCO22
#\\
#\\
# #INTERFACE:
#
  function kco22 ( p, t, lwc, hplus ) result ( kco2p2 )
#
# #INPUT PARAMETERS:
#
    real(fp),  intent(in) :: t, p, lwc, hplus
#
# #OUTPUT PARAMETERS:
#
    real(fp)              :: kco2p, kco2p2
#
# #REVISION HISTORY:
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
    # CO2 dissociation constants
    real(fp),  parameter  :: kc1 = 4.3e-7
    real(fp),  parameter  :: kc2 = 4.68e - 11
    real(fp),  parameter  :: dhrkc1 = - 1000.0
    real(fp),  parameter  :: dhrkc2 = - 1760.0
    real(fp),  parameter  :: hco2 = 3.4e-2
    real(fp),  parameter  :: dhco2 = 2.44e + 3_fp
    # CO2 concentration [v / v]
    real(fp),  parameter  :: co2 = 390.0e-6_fp
#
# #LOCAL VARIABLES:
#
    real(fp)              :: hco2_t, kc1_t, kc2_t
    real(fp)              :: hco2eff, xco2, pco2

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # kco22 begins here#
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    #CO2 dissolution constants
    hco2_t  = hco2 * exp(dhco2 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
    kc1_t   = kc1 * exp(dhrkc1 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
    kc2_t   = kc2 * exp(dhrkc2 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

    #CO2 dissolution
    hco2eff = hco2_t * (1.0e + 0_fp + (kc1_t / hplus) + ((kc1_t * kc2_t) / (hplus * hplus)))
    xco2    = 1.0e + 0_fp / ( 1.0e + 0_fp  + ( hco2eff * 0.08205e + 0_fp * t * lwc ) )
    pco2    = co2 * p * xco2

    kco2p2 = kc1_t * kc2_t * hco2_t / hplus / hplus * pco2

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  GEOS - Chem Global Chemical Transport Model
#                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: dkCO22
#
# #DESCRIPTION: Function dkCO22
#\\
#\\
# #INTERFACE:
#
      function dkco22 ( p, t, lwc, hplus ) result ( kco2p2 )
#
# #INPUT PARAMETERS:
#
      real(fp),  intent(in) :: t, p, lwc, hplus
#
# #OUTPUT PARAMETERS:
#
      real(fp)              :: kco2p, kco2p2
#
# #REVISION HISTORY:
#  25 Jan 2012 - M. Payer - Added ProTeX headers
#  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
#  pressure
#  22 Mar 2017 - M. Sulprizio - Dhco2 value is from Table 7.3 of Seinfeld
#  and
#                              Pandis (2006, pp 289) and should be
#                              positive for
#                              consistency with the way it is used here.
#                              Also,
#                              the value of R, in units of kcal mol - 1
#                              K - 1, is
#                              1.986x10^ - 3, not 0.04.0 (V. Shah)
#  15 Feb 2019 - J. Moch - updated function to make output deriviate
#  of [CO3 -  - ]
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
      # CO2 dissociation constants
      real(fp),  parameter  :: kc1 = 4.3e-7
      real(fp),  parameter  :: kc2 = 4.68e - 11
      real(fp),  parameter  :: dhrkc1 = - 1000.0
      real(fp),  parameter  :: dhrkc2 = - 1760.0
      real(fp),  parameter  :: hco2 = 3.4e-2
      real(fp),  parameter  :: dhco2 = 2.44e + 3_fp
      # CO2 concentration [v / v]
      real(fp),  parameter  :: co2 = 390.0e-6_fp
#
# #LOCAL VARIABLES:
#
      real(fp)              :: hco2_t, kc1_t, kc2_t

      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      # dkco22 begins here#
      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

      #CO2 dissolution constants
      hco2_t = hco2 * exp(dhco2 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
      kc1_t = kc1 * exp(dhrkc1 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
      kc2_t = kc2 * exp(dhrkc2 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

      #CO2 dissolution

      kco2p2 = - 1.0e + 0_fp * kc1_t * kc2_t * hco2_t * co2 * p * ( kc1_t * hco2_t * 0.08205e + 0_fp * t * lwc + 2.0e + 0_fp * hco2_t * 0.08205e + 0_fp * t * lwc * hplus + 2.0e + 0_fp * hplus ) / ( kc1_t * kc2_t * hco2_t * 0.08205e + 0_fp * t * lwc + kc1_t * hco2_t * 0.08205e + 0_fp * t * lwc * hplus + hco2_t * 0.08205e + 0_fp * t * lwc * hplus * hplus + hplus * hplus )^2

      end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: kSO21
#
# #DESCRIPTION: Function kSO21
#\\
#\\
# #INTERFACE:
#
  function kso21 ( p, t, lwc, hplus, so2 ) result ( kso2p )
#
# #INPUT PARAMETERS:
#
    real(fp),  intent(in) :: t, p, lwc, hplus, so2
#
# #OUTPUT PARAMETERS:
#
    real(fp)              :: kso2p, kso2p2
#
# #REVISION HISTORY:
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
    # SO2 dissociation constants
    real(fp),  parameter  :: ks1 = 1.3e-2
    real(fp),  parameter  :: ks2 = 6.6e-8
    real(fp),  parameter  :: hso2 = 1.23
    real(fp),  parameter  :: dhso2 = 3.14e + 3_fp
    real(fp),  parameter  :: dhrkso21 = 1960.0
    real(fp),  parameter  :: dhrkso22 = 1500.0
#
# #LOCAL VARIABLES:
#
    real(fp)              :: hso2_t, ks1_t, ks2_t
    real(fp)              :: hso2eff, xso2, pso2

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # kso21 begins here#
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # SO2 dissolution constants
    hso2_t  = hso2 * exp(dhso2 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
    ks1_t   = ks1 * exp(dhrkso21 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
    ks2_t   = ks2 * exp(dhrkso22 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

    # SO2 dissolution
    hso2eff = hso2_t * (1.0e + 0_fp + (ks1_t / hplus) + ((ks1_t * ks2_t) / (hplus * hplus)))
    xso2    = 1.0e + 0_fp / ( 1.0e + 0_fp  + ( hso2eff * 0.08205e + 0_fp * t * lwc ) )
    pso2    = so2 * p * xso2

    kso2p   = hso2_t * ks1_t * pso2 / hplus

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  GEOS - Chem Global Chemical Transport Model
#                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: dkSO21
#
# #DESCRIPTION: Function dkSO21
#\\
#\\
# #INTERFACE:
#
      function dkso21 ( p, t, lwc, hplus, so2 ) result ( kso2p )
#
# #INPUT PARAMETERS:
#
      real(fp),  intent(in) :: t, p, lwc, hplus, so2
#
# #OUTPUT PARAMETERS:
#
      real(fp)              :: kso2p, kso2p2
#
# #REVISION HISTORY:
#  25 Jan 2012 - M. Payer - Added ProTeX headers
#  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
#  pressure
#  22 Mar 2017 - M. Sulprizio - Dhso2 value is from Table 7.3 of Seinfeld
#  and
#                              Pandis (2006, pp 289) and should be
#                              positive for
#                              consistency with the way it is used here.
#                              Also,
#                              the value of R, in units of kcal mol - 1
#                              K - 1, is
#                              1.986x10^ - 3, not 0.04.0 (V. Shah)
#  15 Feb 2019 - J. Moch - updated function to make output
#  derivative of [HSO3 - ]

#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
      # SO2 dissociation constants
      real(fp),  parameter  :: ks1 = 1.3e-2
      real(fp),  parameter  :: ks2 = 6.6e-8
      real(fp),  parameter  :: hso2 = 1.23
      real(fp),  parameter  :: dhso2 = 3.14e + 3_fp
      real(fp),  parameter  :: dhrkso21 = 1960.0
      real(fp),  parameter  :: dhrkso22 = 1500.0
#
# #LOCAL VARIABLES:
#
      real(fp)              :: hso2_t, ks1_t, ks2_t

      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      # dkso21 begins here#
      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



      # SO2 dissolution constants
      hso2_t = hso2 * exp(dhso2 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
      ks1_t = ks1 * exp(dhrkso21 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
      ks2_t = ks2 * exp(dhrkso22 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))


      kso2p  = ks1_t * hso2_t * so2 * p * ( ks1_t * ks2_t * hso2_t * 0.08205e + 0_fp * t * lwc - hso2_t * 0.08205e + 0_fp * t * lwc * hplus * hplus - hplus * hplus) / (ks1_t * ks2_t * hso2_t * 0.08205e + 0_fp * t * lwc + ks1_t * hso2_t * 0.08205e + 0_fp * t * lwc * hplus + hso2_t * 0.08205e + 0_fp * t * lwc * hplus * hplus + hplus * hplus)^2

      end
#EOC

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: kSO22
#
# #DESCRIPTION: Function kSO22
#\\
#\\
# #INTERFACE:
#
  function kso22 ( p, t, lwc, hplus, so2 ) result ( kso2p2 )
#
# #INPUT PARAMETERS:
#
    real(fp),  intent(in) :: t, p, lwc, hplus, so2
#
# #OUTPUT PARAMETERS:
#
    real(fp)              :: kso2p, kso2p2
#
# #REVISION HISTORY:
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
    # SO2 dissociation constants
    real(fp),  parameter  :: ks1 = 1.3e-2
    real(fp),  parameter  :: ks2 = 6.6e-8
    real(fp),  parameter  :: hso2 = 1.23
    real(fp),  parameter  :: dhso2 = 3.14e + 3_fp
    real(fp),  parameter  :: dhrkso21 = 1960.0
    real(fp),  parameter  :: dhrkso22 = 1500.0
#
# #LOCAL VARIABLES:
#
    real(fp)              :: hso2_t, ks1_t, ks2_t
    real(fp)              :: hso2eff, xso2, pso2

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # kso22 begins here#
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # SO2 dissolution constants
    hso2_t  = hso2 * exp(dhso2 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
    ks1_t   = ks1 * exp(dhrkso21 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
    ks2_t   = ks2 * exp(dhrkso22 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

    #SO2 dissolution
    hso2eff = hso2_t * (1.0e + 0_fp + (ks1_t / hplus) + ((ks1_t * ks2_t) / (hplus * hplus)))
    xso2    = 1.0e + 0_fp / ( 1.0e + 0_fp + ( hso2eff * 0.08205e + 0_fp * t * lwc ) )
    pso2    = so2 * p * xso2

    kso2p2 = ks1_t * ks2_t * hso2_t / hplus / hplus * pso2

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  GEOS - Chem Global Chemical Transport Model
#                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: dkSO22
#
# #DESCRIPTION: Function dkSO22
#\\
#\\
# #INTERFACE:
#
      function dkso22 ( p, t, lwc, hplus, so2 ) result ( kso2p2 )
#
# #INPUT PARAMETERS:
#
      real(fp),  intent(in) :: t, p, lwc, hplus, so2
#
# #OUTPUT PARAMETERS:
#
      real(fp)              :: kso2p, kso2p2
#
# #REVISION HISTORY:
#  25 Jan 2012 - M. Payer - Added ProTeX headers
#  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
#  pressure
#  22 Mar 2017 - M. Sulprizio - Dhso2 value is from Table 7.3 of Seinfeld
#  and
#                              Pandis (2006, pp 289) and should be
#                              positive for
#                              consistency with the way it is used here.
#                              Also,
#                              the value of R, in units of kcal mol - 1
#                              K - 1, is
#                              1.986x10^ - 3, not 0.04.0 (V. Shah)
#  15 Feb 2019 - J. Moch - updated function to make output
#  derivative [SO3 -  - ]

#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
      # SO2 dissociation constants
      real(fp),  parameter  :: ks1 = 1.3e-2
      real(fp),  parameter  :: ks2 = 6.6e-8
      real(fp),  parameter  :: hso2 = 1.23
      real(fp),  parameter  :: dhso2 = 3.14e + 3_fp
      real(fp),  parameter  :: dhrkso21 = 1960.0
      real(fp),  parameter  :: dhrkso22 = 1500.0
#
# #LOCAL VARIABLES:
#
      real(fp)              :: hso2_t, ks1_t, ks2_t

      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      # dkso22 begins here#
      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      # SO2 dissolution constants
      hso2_t = hso2 * exp(dhso2 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
      ks1_t  = ks1 * exp(dhrkso21 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
      ks2_t  = ks2 * exp(dhrkso22 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

      kso2p2 = - 1.0e + 0_fp * ks1_t * ks2_t * hso2_t * so2 * p * ( ks1_t * hso2_t * 0.08205e + 0_fp * t * lwc + 2.0e + 0_fp * hso2_t * 0.08205e + 0_fp * t * lwc * hplus + 2.0e + 0_fp * hplus ) / ( ks1_t * ks2_t * hso2_t * 0.08205e + 0_fp * t * lwc + ks1_t * hso2_t * 0.08205e + 0_fp * t * lwc * hplus + hso2_t * 0.08205e + 0_fp * t * lwc * hplus * hplus + hplus * hplus )^2

      end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: kHNO3
#
# #DESCRIPTION: Function kNO3
#\\
#\\
# #INTERFACE:
#
  function khno3 ( p, t, lwc, hplus, hno3 ) result ( khno3p )
#
# #INPUT PARAMETERS:
#
    real(fp),  intent(in) :: t, p, lwc, hplus, hno3
#
# #OUTPUT PARAMETERS:
#
    real(fp)              :: khno3p
#
# #REVISION HISTORY:
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
    # HNO3 dissociation constants
    real(fp),  parameter  :: kn1 = 15.4
    real(fp),  parameter  :: hhno3 = 2.1e5
    real(fp),  parameter  :: dhhno3 = 0.0
    real(fp),  parameter  :: dhrkn1 = 8700.0
#
# #LOCAL VARIABLES:
#
    real(fp)              :: hhno3_t, kn1_t
    real(fp)              :: hhno3eff, xhno3, phno3

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # khno3 begins here#
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # HNO3 dissolution constants
    hhno3_t  = hhno3 * exp(dhhno3 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
    kn1_t    = kn1 * exp(dhrkn1 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

    # HNO3 dissolution
    # The original Hhno3eff expression is valid for 298K (Seinfeld and Pandis
    # 2006, pp 299 - 301), and Kn1 has a strong temperature dependence. The
    # fix follows Eq. 7.59 of Seinfeld and Pandis (2006, pp 301).
    #Hhno3eff = 3.2e6 / HPLUS
    hhno3eff = hhno3_t * (1.0e + 0_fp + (kn1_t / hplus))
    xhno3    = 1.0e + 0_fp / ( 1.0e + 0_fp + ( hhno3eff * 0.08205e + 0_fp * t * lwc ) )
    phno3    = hno3 * p * xhno3

    khno3p = hhno3_t * kn1_t * phno3 / hplus

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  GEOS - Chem Global Chemical Transport Model
#                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: dkHNO3
#
# #DESCRIPTION: Function dkNO3
#\\
#\\
# #INTERFACE:
#
      function dkhno3 ( p, t, lwc, hplus, hno3 ) result ( khno3p )
#
# #INPUT PARAMETERS:
#
      real(fp),  intent(in) :: t, p, lwc, hplus, hno3
#
# #OUTPUT PARAMETERS:
#
      real(fp)              :: khno3p
#
# #REVISION HISTORY:
#  25 Jan 2012 - M. Payer - Added ProTeX headers
#  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
#  pressure
#  22 Mar 2017 - M. Sulprizio - Add fix for Hhno3eff from V. Shah
#  15 Feb 2019 - J. Moch - updated function to make output
#  derivative of [HNO3 - ]
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
      # HNO3 dissociation constants
      real(fp),  parameter  :: kn1 = 15.4
      real(fp),  parameter  :: hhno3 = 2.1e5
      real(fp),  parameter  :: dhhno3 = 0.0
      real(fp),  parameter  :: dhrkn1 = 8700.0
#
# #LOCAL VARIABLES:
#
      real(fp)              :: hhno3_t, kn1_t

      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      # dkhno3 begins here#
      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

      # HNO3 dissolution constants
      hhno3_t = hhno3 * exp(dhhno3 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
      kn1_t = kn1 * exp(dhrkn1 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

      # HNO3 dissolution
      # The original Hhno3eff expression is valid for 298K (Seinfeld and
      # Pandis
      # 2006, pp 299 - 301), and Kn1 has a strong temperature dependence.
      # The
      # fix follows Eq. 7.59 of Seinfeld and Pandis (2006, pp 301).

      khno3p = - 1.0e + 0_fp * kn1_t * hhno3_t * hno3 * p * ( 1.0e + 0_fp + hhno3_t * 0.08205e + 0_fp * t * lwc ) / ( kn1_t * hhno3_t * 0.08205e + 0_fp * t * lwc + hhno3_t * 0.08205e + 0_fp * t * lwc * hplus + hplus )^2

      end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: kHCl
#
# #DESCRIPTION: Function kHCl
#\\
#\\
# #INTERFACE:
#
  function khcl ( p, t, lwc, hplus, cl ) result ( khclp )
#
# #INPUT PARAMETERS:
#
    real(fp),  intent(in) :: t, p, lwc, hplus, cl
#
# #OUTPUT PARAMETERS:
#
    real(fp)              :: khclp
#
# #REVISION HISTORY:
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
    # HNO3 dissociation constants
    real(fp),  parameter  :: kcl = 1.74e + 6_fp
    real(fp),  parameter  :: hcl = 1.5e + 3_fp
    real(fp),  parameter  :: dhcl = 2.3e + 3_fp
    real(fp),  parameter  :: dhrkcl = 6900.0e + 0_fp
#
# #LOCAL VARIABLES:
#
    real(fp)              :: hcl_t, kcl_t
    real(fp)              :: hcleff, xcl, phcl

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # khcl begins here#
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # HCl dissolution constants
    hcl_t  = hcl * exp(dhcl * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
    kcl_t  = kcl * exp(dhrkcl * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

    #HCl dissolution
    hcleff = hcl_t * (1.0e + 0_fp + (kcl_t / hplus))
    xcl    = 1.0e + 0_fp / ( 1.0e + 0_fp + ( hcleff * 0.08205e + 0_fp * t * lwc ) )
    phcl   = cl * p * xcl

    khclp  = hcl_t * kcl_t * phcl / hplus

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  GEOS - Chem Global Chemical Transport Model
#                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: dkHCl
#
# #DESCRIPTION: Function dkHCl
#\\
#\\
# #INTERFACE:
#
      function dkhcl ( p, t, lwc, hplus, cl ) result ( khclp )
#
# #INPUT PARAMETERS:
#
      real(fp),  intent(in) :: t, p, lwc, hplus, cl
#
# #OUTPUT PARAMETERS:
#
      real(fp)              :: khclp
#
# #REVISION HISTORY:
#  03 Apr 2019 - X. Wang - Initial version
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
      # HCl dissociation constants
      real(fp),  parameter  :: kcl = 1.74e + 6_fp
      real(fp),  parameter  :: hcl = 1.5e + 3_fp
      real(fp),  parameter  :: dhcl = 2.3e + 3_fp
      real(fp),  parameter  :: dhrkcl = 6900.0e + 0_fp
#
# #LOCAL VARIABLES:
#
      real(fp)              :: hcl_t, kcl_t

      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      # dkhcl begins here#
      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

      # HCl dissolution constants
      hcl_t = hcl * exp(dhcl * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
      kcl_t = kcl * exp(dhrkcl * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

      # HCl dissolution
      # The fix follows Eq. 7.59 of Seinfeld and Pandis (2006, pp 301).

      khclp = - 1.0e + 0_fp * kcl_t * hcl_t * cl * p * ( 1.0e + 0_fp + hcl_t * 0.08205e + 0_fp * t * lwc ) / ( kcl_t * hcl_t * 0.08205e + 0_fp * t * lwc + hcl_t * 0.08205e + 0_fp * t * lwc * hplus + hplus )^2

      end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: kNH3
#
# #DESCRIPTION: Function kNH3
#\\
#\\
# #INTERFACE:
#
  function knh3 ( p, t, lwc, hplus, nh3, kw ) result ( knh3p )
#
# #INPUT PARAMETERS:
#
    real(fp),  intent(in) :: t, p, lwc, hplus, nh3, kw
#
# #OUTPUT PARAMETERS:
#
    real(fp)              :: knh3p
#
# #REVISION HISTORY:
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
    # NH3 dissociation contants
    real(fp),  parameter  :: ka1 = 1.7e-5
    real(fp),  parameter  :: hnh3 = 60.0
    real(fp),  parameter  :: dhnh3 = 4200e + 0_fp
    real(fp),  parameter  :: dhrka1 = - 450.0

    # Variables
    real(fp)              :: hnh3_t, ka1_t
    real(fp)              :: hnh3eff, xnh3, pnh3

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # knh3 begins here#
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    #NH3 dissolution constants
    hnh3_t  = hnh3 * exp(dhnh3 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
    ka1_t   = ka1 * exp(dhrka1 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

    #NH3 dissolution
    hnh3eff = hnh3_t * (1.0e + 0_fp + ((ka1_t * hplus) / kw))
    xnh3    = 1.0e + 0_fp / ( 1.0e + 0_fp + ( hnh3eff * 0.08205e + 0_fp * t * lwc ) )
    pnh3    = nh3 * p * xnh3

    knh3p   = hplus * hnh3_t * ka1_t * pnh3 / kw

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  GEOS - Chem Global Chemical Transport Model
#                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: dkNH3
#
# #DESCRIPTION: Function dkNH3
#\\
#\\
# #INTERFACE:
#
      function dknh3 ( p, t, lwc, hplus, nh3, kw ) result ( knh3p )
#

# #INPUT PARAMETERS:
#
      real(fp),  intent(in) :: t, p, lwc, hplus, nh3, kw
#
# #OUTPUT PARAMETERS:
#
      real(fp)              :: knh3p
#
# #REVISION HISTORY:
#  25 Jan 2012 - M. Payer - Added ProTeX headers
#  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
#  pressure
#  22 Mar 2017 - M. Sulprizio - Dhnh3 value is from Table 7.3 of Seinfeld
#  and
#                              Pandis (2006, pp 289) and should be
#                              positive for
#                              consistency with the way it is used here.
#                              Also,
#                              the value of R, in units of kcal mol - 1
#                              K - 1, is
#                              1.986x10^ - 3, not 0.04.0 (V. Shah)
#  15 Feb 2019 - J. Moch - updated function to make output
#  derivative of [NH4 + ]
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
      # NH3 dissociation contants
      real(fp),  parameter  :: ka1 = 1.7e-5
      real(fp),  parameter  :: hnh3 = 60.0
      real(fp),  parameter  :: dhnh3 = 4200e + 0_fp
      real(fp),  parameter  :: dhrka1 = - 450.0

      # Variables
      real(fp)              :: hnh3_t, ka1_t

      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      # dknh3 begins here#
      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

      #NH3 dissolution constants
      hnh3_t = hnh3 * exp(dhnh3 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
      ka1_t = ka1 * exp(dhrka1 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

      #NH3 dissolutionnyn

      knh3p = ka1_t * hnh3_t * nh3 * kw * p * ( 1.0e + 0_fp + hnh3_t * 0.08205e + 0_fp * t * lwc ) / ( hnh3_t * 0.08205e + 0_fp * t * lwc * ( kw + ka1_t * hplus ) + kw)^2

      end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  GEOS - Chem Global Chemical Transport Model
#                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: kFA
#
# #DESCRIPTION: Function kFA
#\\
#\\
# #INTERFACE:
#
      function kfa ( p, t, lwc, hplus, fa ) result ( kfap )
#

# #INPUT PARAMETERS:
#
      real(fp),  intent(in) :: t, p, lwc, hplus, fa
#
# #OUTPUT PARAMETERS:
#
      real(fp)              :: kfap
#
# #REVISION HISTORY:
#  25 Jan 2012 - M. Payer - Added ProTeX headers
#  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
#  pressure
#  17 Oct 2017 - M. Sulprizio - Dhck value is from Table 7.3 of Seinfeld
#  and
#                              Pandis (2006, pp 289) and should be
#                              positive for
#                              consistency with the way it is used here.
#                              Also,
#                              the value of R, in units of kcal mol - 1
#                              K - 1, is
#                              1.986x10^ - 3, not 0.04.0 (Qianjie Chen)
#  03 Dec 2018 - J. Moch - Modified for formic acid (HCOOH). Values
#                              taken from Sienfeld and Pandis. Made it
#                              to output is [FA]
#  01 May 2020 - V. Shah - Use correct equilibrium constants
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
      # HCOOH dissociation constants
      real(fp),  parameter  :: kformate = 1.8e-4_fp # equib const
      real(fp),  parameter  :: hfa = 8800e + 0_fp # henry const
      real(fp),  parameter  :: dhfa = 6100e + 0_fp # henry temp
      real(fp),  parameter  :: dhrkfa = 151.0e + 0_fp # equib temp
#
# #LOCAL VARIABLES:
#
      real(fp)              :: hfa_t, kfa_t
      real(fp)              :: hfaeff, xfa, pfa

      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      # kfa begins here#
      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

      # Formic acid dissolution constants
      hfa_t = hfa * exp(dhfa * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
      kfa_t = kformate * exp(dhrkfa * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

      #HCOOH  dissolution
      hfaeff = hfa_t * (1.0e + 0_fp + (kfa_t / hplus))
      xfa = 1.0e + 0_fp / ( 1.0e + 0_fp + ( hfaeff * 0.08205e + 0_fp * t * lwc ) )
      pfa = fa * p * xfa

      kfap = hfa_t * kfa_t * pfa / hplus

      end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  GEOS - Chem Global Chemical Transport Model
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: dkFA
#
# #DESCRIPTION: Function dkFA
#\\
#\\
# #INTERFACE:
#
      function dkfa ( p, t, lwc, hplus, fa ) result ( kfap )
#
# #INPUT PARAMETERS:
#
      real(fp),  intent(in) :: t, p, lwc, hplus, fa
#
# #OUTPUT PARAMETERS:
#
      real(fp)              :: kfap
#
# #REVISION HISTORY:
#  25 Jan 2012 - M. Payer - Added ProTeX headers
#  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
#  pressure
#  17 Oct 2017 - M. Sulprizio - Dhck value is from Table 7.3 of Seinfeld
#  and
#                              Pandis (2006, pp 289) and should be
#                              positive for
#                              consistency with the way it is used here.
#                              Also,
#                              the value of R, in units of kcal mol - 1
#                              K - 1, is
#                              1.986x10^ - 3, not 0.04.0 (Qianjie Chen)
#  03 Dec 2018 - J. Moch - Modified for formic acid (HCOOH). Values
#  taken from
#                              Sienfeld and Pandis. Made it to output is
#                              [FA]
#  01 May 2020 - V. Shah - Use correct equilibrium constants
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
      # HCOOH dissociation constants
      real(fp),  parameter  :: kformate = 1.8e-4_fp # equib const
      real(fp),  parameter  :: hfa = 8800e + 0_fp # henry const
      real(fp),  parameter  :: dhfa = 6100e + 0_fp # henry temp
      real(fp),  parameter  :: dhrkfa = 151.0e + 0_fp # equib temp
#
# #LOCAL VARIABLES:
#
      real(fp)              :: hfa_t, kfa_t

      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      # dkfa begins here#
      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

      # Formic acid dissolution constants
      hfa_t = hfa * exp(dhfa * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
      kfa_t = kformate * exp(dhrkfa * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

      #HCOOH  dissolution

      kfap = - 1.0e + 0_fp * kfa_t * hfa_t * fa * p * ( 1.0e + 0_fp + hfa_t * 0.08205e + 0_fp * t * lwc ) / ( kfa_t * hfa_t * 0.08205e + 0_fp * t * lwc + hfa_t * 0.08205e + 0_fp * t * lwc * hplus + hplus )^2

      end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  GEOS - Chem Global Chemical Transport Model
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: kAA
#
# #DESCRIPTION: Function kAA
#\\
#\\
# #INTERFACE:
#
      function kaa ( p, t, lwc, hplus, aa ) result ( kaap )
#
# #INPUT PARAMETERS:
#
      real(fp),  intent(in) :: t, p, lwc, hplus, aa
#
# #OUTPUT PARAMETERS:
#
      real(fp)              :: kaap
#
# #REVISION HISTORY:
#  25 Jan 2012 - M. Payer - Added ProTeX headers
#  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
#  pressure
#  17 Oct 2017 - M. Sulprizio - Dhck value is from Table 7.3 of Seinfeld
#  and
#                              Pandis (2006, pp 289) and should be
#                              positive for
#                              consistency with the way it is used here.
#                              Also,
#                              the value of R, in units of kcal mol - 1
#                              K - 1, is
#                              1.986x10^ - 3, not 0.04.0 (Qianjie Chen)
#  03 Dec 2018 - J. Moch - Modified for acetic acid (CH3COOH).
#  Values taken from
#                              Sienfeld and Pandis, value of [HCOOH]
#  01 May 2020 - V. Shah - Use correct equilibrium constants
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
      # CH3HCOOH dissociation constants
      real(fp),  parameter  :: kacetate = 1.75e - 5_fp
      real(fp),  parameter  :: haa = 4100e + 0_fp
      real(fp),  parameter  :: dhaa = 6200e + 0_fp
      real(fp),  parameter  :: dhrkaa = 50.0e + 0_fp
#
# #LOCAL VARIABLES:
#
      real(fp)              :: haa_t, kaa_t
      real(fp)              :: haaeff, xaa, paa
      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      # kaa begins here#
      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

      # Formic acid dissolution constants
      haa_t = haa * exp(dhaa * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
      kaa_t = kacetate * exp(dhrkaa * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

      #HCOOH  dissolution
      haaeff = haa_t * (1.0e + 0_fp + (kaa_t / hplus))
      xaa = 1.0e + 0_fp / ( 1.0e + 0_fp + ( haaeff * 0.08205e + 0_fp * t * lwc ) )
      paa = aa * p * xaa

      kaap = haa_t * kaa_t * paa / hplus

      end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  GEOS - Chem Global Chemical Transport Model
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: dkAA
#
# #DESCRIPTION: Function kdAA
#\\
#\\
# #INTERFACE:
#
      function dkaa ( p, t, lwc, hplus, aa ) result ( kaap )
#
# #INPUT PARAMETERS:
#
      real(fp),  intent(in) :: t, p, lwc, hplus, aa
#
# #OUTPUT PARAMETERS:
#
      real(fp)              :: kaap
#
# #REVISION HISTORY:
#  25 Jan 2012 - M. Payer - Added ProTeX headers
#  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
#  pressure
#  17 Oct 2017 - M. Sulprizio - Dhck value is from Table 7.3 of Seinfeld
#  and
#                              Pandis (2006, pp 289) and should be
#                              positive for
#                              consistency with the way it is used here.
#                              Also,
#                              the value of R, in units of kcal mol - 1
#                              K - 1, is
#                              1.986x10^ - 3, not 0.04.0 (Qianjie Chen)
#  03 Dec 2018 - J. Moch - Modified for acetic acid (CH3COOH).
#  Values taken from
#                              Sienfeld and Pandis. Output is
#                              derivative.
#  01 May 2020 - V. Shah - Use correct equilibrium constants
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
      # HCOOH dissociation constants
      real(fp),  parameter  :: kacetate = 1.75e - 5_fp
      real(fp),  parameter  :: haa = 4100e + 0_fp
      real(fp),  parameter  :: dhaa = 6200e + 0_fp
      real(fp),  parameter  :: dhrkaa = 50.0e + 0_fp
#
# #LOCAL VARIABLES:
#
      real(fp)              :: haa_t, kaa_t
      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      # kaa begins here#
      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

      # Formic acid dissolution constants
      haa_t = haa * exp(dhaa * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
      kaa_t = kacetate * exp(dhrkaa * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

      #HCOOH  dissolution
      kaap = - 1.0e + 0_fp * kaa_t * haa_t * aa * p * ( 1.0e + 0_fp + haa_t * 0.08205e + 0_fp * t * lwc ) / ( kaa_t * haa_t * 0.08205e + 0_fp * t * lwc + haa_t * 0.08205e + 0_fp * t * lwc * hplus + hplus )^2


      end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: CaCO3_PRECIP
#
# #DESCRIPTION: Subroutine CaCO3 to calculate [Ca + + ] in equilibrium with
# CaCO3(s) (dust particles) depending on [H + ]
#\\
#\\
# #INTERFACE:
#
      subroutine caco3_precip ( p,  t, hplus, fca, dca )
#
# #INPUT PARAMETERS:
#
      real(fp),        intent(in) :: t, p, hplus
#
# #OUTPUT PARAMETERS:
#
      real(fp),  intent(out):: fca, dca # [Ca2 + ] and d([Ca2 + ]) / d[H + ]
#
# #REVISION HISTORY:
#  25 Dec 2019 - V. Shah - Initial version
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
# #DEFINED PARAMETERS:
#
      real(fp),  parameter  :: kc1 = 4.3e-7_fp
      real(fp),  parameter  :: kc2 = 4.68e - 11_fp
      real(fp),  parameter  :: dhrkc1 = - 1000.0
      real(fp),  parameter  :: dhrkc2 = - 1760.0
      real(fp),  parameter  :: hco2 = 3.4e-2_fp
      real(fp),  parameter  :: dhco2 = 2.44e + 3_fp
      # CO2 concentration [v / v]
      real(fp),  parameter  :: co2 = 390.0e-6_fp
      real(fp),  parameter  :: ksp = 3.3e-9_fp
      real(fp),  parameter  :: dhrksp = - 1200e + 0_fp

# #LOCAL VARIABLES:
      real(fp)              :: hco2_t, kc1_t, kc2_t, ksp_t

# #REMARKS:

      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      # caco3_precip begins here#
      # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      #Temperature adjusted eq. constants
      #CO2 dissolution constants
      hco2_t = hco2 * exp(dhco2 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
      kc1_t = kc1 * exp(dhrkc1 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))
      kc2_t = kc2 * exp(dhrkc2 * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

      # CaCO3 eq constants
      ksp_t = ksp * exp(dhrksp * ((1.0e + 0_fp / t) - (1.0e + 0_fp / 298.0e + 0_fp)))

      #Ca concentrations [M]
      fca = ksp_t * hplus * hplus / (kc1_t * kc2_t * hco2_t * co2 * p)
      #derivative d[Ca2 + ] / dH +
      dca  = 2e + 0_fp * ksp_t * hplus / (kc1_t * kc2_t * hco2_t * co2 * p)

      end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: aqchem_so2
#
# #DESCRIPTION: Subroutine AQCHEM\_SO2 computes the reaction rates for aqueous
# SO2 chemistry. (rjp, bmy, 10 / 31 / 02, 12 / 12 / 02)
#\\
#\\
# #INTERFACE:
#
  subroutine aqchem_so2( i,      j,       l,      lwc,     t,      p, so2,    h2o2,    o3,     hcho,    hplus,  mnii, feiii,  kaqh2o2, kaqo3,  kaqo3_1, kaqo2,  hso3aq, so3aq,  kaqhcho, kaqhms, kaqhms2                   )
#
# #INPUT PARAMETERS:
#
    integer,  intent(in)  :: i, j, l # Coordinates, for diagnostic use -  - MSL

    real(fp), intent(in)  :: lwc     # Liq water content [m3 / m3] = 1.0E - 6 * L [g / m3]
    real(fp), intent(in)  :: t       # Temperature [K]
    real(fp), intent(in)  :: p       # Dry air partial pressure [atm]
    real(fp), intent(in)  :: so2     # SO2  mixing ratio [v / v]
    real(fp), intent(in)  :: h2o2    # H2O2 mixing ratio [v / v]
    real(fp), intent(in)  :: o3      # O3   mixing ratio [v / v]
    real(fp), intent(in)  :: hplus   # Concentration of H + ion (i.e. pH) [v / v]
    real(fp), intent(in)  :: mnii    # Concentration of MnII [mole / l]
    real(fp), intent(in)  :: feiii   # Concentration of FeIII [mole / l]
    real(fp), intent(in)  :: hcho    # HCHO   mixing ratio [v / v] (jmm, 06 / 13 / 18)

#
# #OUTPUT PARAMETERS:
#
    real(fp), intent(out) :: kaqh2o2 # Reaction rate for H2O2
    real(fp), intent(out) :: kaqo3   # Reaction rate for O3
    real(fp), intent(out) :: kaqo3_1 # only the SO3 -  - oxidation, (qjc, 04 / 10 / 16)
    real(fp), intent(out) :: kaqo2   # Reaction rate for O2 (metal cat)
    real(fp), intent(out) :: kaqhcho # Reaction rate for SO2 and HCHO (jmm, 06 / 13 / 18)
    real(fp), intent(out) :: kaqhms  # Reaction rate for HMS and OH - (jmm, 06 / 13 / 18)
    real(fp), intent(out) :: kaqhms2 # Reaction rate for HMS and OH(aq) (jmm, 06 / 28 / 18)
    real(fp), intent(out) :: hso3aq  # Cloud bisulfite [mol / l] (qjc, 06 / 10 / 16)
    real(fp), intent(out) :: so3aq   # Cloud sulfite   [mol / l] (qjc, 06 / 10 / 16)
#
# #REMARKS:
#  Chemical Reactions:
#  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#  (R1) HSO3 - + H2O2(aq) + H + = > SO4 -  - + 2H + + H2O [Jacob, 1986]
#                                                                             .
#      d[S(VI)] / dt = k[H + ][H2O2(aq)][HSO3 - ] / (1 + K[H + ])
#      [Seinfeld and Pandis, 1998, page 366]
#                                                                             .
#  (R2) SO2(aq) + O3(aq) = >
#       HSO3 - + O3(aq) = >
#       SO3 -  - + O3(aq) = >
#       [Jacob, 1986; Jacobson, 1999]
#                                                                             .
#       d[S(VI)] / dt = (k0[SO2(aq)] + k1[HSO3 - ] + K2[SO3 -  - ])[O3(aq)]
#       [Seinfeld and Pandis, 1998, page 363]
#                                                                             .
#  (R3) HSO3 - + HCHO(aq) = > HMS
#       SO3 -  - + HCHO(aq) = > HMS + OH - 
#       [Moch et al., 2018; Olson and Hoffman, 1986]
#                                                                             .
#       d[S(HMS)] / dt = (k1[HSO3 - ] + k2[SO3 -  - ])[HCHO(aq)]
#       [Seinfeld and Pandis, 2016, 309]
#
#  (R4) HMS + OH - = > HCHO(aq) + SO3 -  - 
#       [Moch et al., 2018; Deister et al., 1986]
#        (note treated as 1st order in contrast to other reactions here)
#
#  (R5) HMS + OH(aq) = (SO2, HO2, O2) = > HCHO + 2SO4 -  - + O2 + 3H + + 2H2O
#       [Jacob et al, 1986, Olson and Fessenden, 1992;
#        Seinfeld and Pandis, 2016, Table 7A.7]
#          Net reaction (R5):
#           HMS + OH(aq) = (O2) = > SO5 - + HCHO + H2O
#           HO2 < = > H + + O2 - 
#           SO5 - + O2 - = (H2O) = > HSO5 - + OH - + O2
#           SO2(aq) < = > HSO3 - + H +
#           H + + OH - < = > H2O
#           HSO5 - + HSO3 - = > 2SO4 -  - + 2H +
#
#  Reaction rates can be given as
#       Ra     = k [H2O2(ag)] [S(IV)]  [mole / liter * s]  OR
#       Krate  = Ra LWC R T / P        [1 / s]
#                                                                             .
#  Where:
#       LWC = Liquid water content(g / m3) * 10 - 6 [m3(water) / m3(gas)]
#       R   = 0.08205  (atm L / mol - K), Universal gas const.
#       T   = Temperature (K)
#       P   = Pressure (atm)
#                                                                             .
#  Procedure:
#  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#  (a ) Given [SO2] which is assumed to be total SO2 (gas + liquid) in
#        equilibrium between gas and liquid phase.
#                                                                             .
#  (b ) We can compute SO2(g) using Henry"s law
#          P(so2(g)) = Xg * [SO2]
#          Xg = 1 / (1 + Faq), Fraction of SO2 in gas
#       where:
#          Faq   = Kheff * R * T * LWC,
#          KHeff = Effective Henry"s constant
#                                                                             .
#  (c ) Then Calculate Aquous phase, S[IV] concentrations
#        S[IV] = Kheff * P(so2(g) in atm) [M]
#                                                                             .
#  (d ) The exact same procedure is applied to calculate H2O2(aq) and HCHO(aq)
#
# #REVISION HISTORY:
#  (1 ) Updated by Rokjin Park (rjp, bmy, 12 / 12 / 02)
#  See https: / / github.com / geoschem / geos - chem for complete history
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
    real(fp), parameter  :: r   = 0.08205e + 0_fp
    real(fp), parameter  :: doh = 1.0e-19_fp # [M cm^3 molec^ - 1]
#
# #LOCAL VARIABLES:
#
    real(fp)             :: kh2o2,   ra,     ks1,    ks2,    hcso2
    real(fp)             :: fhcso2,  xso2g,  siv,    hso3,   xso2aq
    real(fp)             :: xhso3,   xso3,   kh1,    hch2o2, fhch2o2
    real(fp)             :: xh2o2g,  h2o2aq, ko0,    ko1,    ko2
    real(fp)             :: hco3,    xo3g,   o3aq,   xhchog, hchcho
    real(fp)             :: fhchcho, khcho1, khcho2, khms,   kw1
    real(fp)             :: khc1,    khms2

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # aqchem_so2 begins here#
    #
    # Aqueous reaction rate
    # HSO3 - + H2O2 + H + = > SO4 -  - + 2H + + H2O [Jacob, 1986]
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # [Jacob, 1986]
    kh2o2 = 6.31e + 14_fp * exp( - 4.76e + 3_fp / t )

    ## [Jacobson, 1999]
    #KH2O2 = 7.45e + 0_fp7 * EXP( - 15.96e + 0_fp * ( (298.15 / T) - 1.0) ) / #        ( 1.0e + 0_fp + 13.0e + 0_fp * Hplus)

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Equilibrium reaction of SO2 - H2O
    #    SO2 + H2O = SO2(aq)        (s0)
    #    SO2(ag)   = HSO3 - + H +     (s1)
    #    HSO3 - = SO3 -  - + H +     (s2)
    #
    # Reaction constant for Aqueous chemistry -  - No big difference
    # between Jacob and Jacobson, choose one of them.
    #
    # Reaction rate dependent on Temperature is given
    #   H = A exp ( B (T. / T - 1) )
    #
    # For equilibrium reactions of SO2:
    #            As1      Bs1   As2      Bs2
    #  Seinfeld  1.30d - 2  7.02  6.60d - 8  3.76   [1998]
    #  Jacob     1.30d - 2  6.75  6.31d - 8  5.05   [1986]
    #  Jacobson  1.71d - 2  7.04  5.99d - 8  3.74   [1996]
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    ks1    = 1.30e - 2_fp * exp( 6.75e + 0_fp * ( 298.15e + 0_fp / t - 1.0e + 0_fp ) )
    ks2    = 6.31e - 8_fp * exp( 5.05e + 0_fp * ( 298.15e + 0_fp / t - 1.0e + 0_fp ) )

    # SIV Fraction
    xso2aq = 1.0e + 0_fp / (1.0e + 0_fp + ks1 / hplus + ks1 * ks2 / (hplus * hplus))
    xhso3  = 1.0e + 0_fp / (1.0e + 0_fp + hplus / ks1 + ks2 / hplus)
    xso3   = 1.0e + 0_fp / (1.0e + 0_fp + hplus / ks2 + hplus * hplus / (ks1 * ks2))

    # Henry"s constant [mol / l - atm] and Effective Henry"s constant for SO2
    hcso2  = 1.22e + 0_fp * exp( 10.55e + 0_fp * ( 298.15e + 0_fp / t - 1.0e + 0_fp) )
    fhcso2 = hcso2 * (1.0e + 0_fp + (ks1 / hplus) + (ks1 * ks2 / (hplus * hplus)))

    xso2g  = 1.0e + 0_fp / ( 1.0e + 0_fp + ( fhcso2 * r * t * lwc ) )
    siv    = fhcso2 * xso2g * so2 * p
    #HSO3   = Ks1 * HCSO2 * XSO2g * SO2 * P

    # Effective HSO3aq for HOBr + HSO3
    hso3aq = siv * xhso3           # unit: M (qjc, 06 / 10 / 16)

    # Effective SO3aq for HOBr + SO3
    so3aq  = siv * xso3            # unit: M (qjc, 06 / 10 / 16)

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # H2O2 equilibrium reaction
    # H2O2 + H2O = H2O2.H2O
    # H2O2.H2O   = HO2 - + H +   1)
    #
    # Reaction rate dependent on Temperature is given
    #   H = A exp ( B (T. / T - 1) )
    #
    # For equilibrium reactions of SO2
    #            Ah1       Bh1
    #  Jacob     1.58E - 12 - 12.49  [1986]
    #  Jacobson  2.20E - 12 - 12.52  [1996]
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    kh1 = 2.20e - 12_fp * exp( - 12.52e + 0_fp * ( 298.15e + 0_fp / t - 1.0e + 0_fp ) )

    # Henry"s constant [mol / l - atm] and Effective Henry"s constant for H2O2
    # [Seinfeld and Pandis, 1998]
    # HCH2O2  = 7.45D4 * EXP( 24.48e + 0_fp * ( 298.15e + 0_fp / T - 1.0e + 0_fp) )

    # [Jacobson, 1999]
    hch2o2  = 7.45e + 4_fp * exp( 22.21e + 0_fp * (298.15e + 0_fp / t - 1.0e + 0_fp) )
    fhch2o2 = hch2o2 * (1.0e + 0_fp + (kh1 / hplus))

    xh2o2g  = 1.0e + 0_fp / ( 1.0e + 0_fp + ( fhch2o2 * r * t * lwc ) )
    #H2O2aq  = FHCH2O2 * XH2O2g * H2O2 * P

    # Conversion rate from SO2 to SO4 via reaction with H2O2
    kaqh2o2  = kh2o2 * ks1 * fhch2o2 * hcso2 * xh2o2g * xso2g * p * lwc * r * t            # [v / v / s]

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #  Aqueous reactions of SO2 with O3
    #  SO2(aq) + O3 = >                       (0)
    #  HSO3 - + O3 = > SO4 -  - + H + + O2       (1)
    #  SO3 -  - + O3 = > SO4 -  - + O2            (2)
    #
    # NOTE
    # [Jacob, 1986]
    #    KO1  = 3.49E12 * EXP( - 4.83E3 / T )
    #    KO2  = 7.32E14 * EXP( - 4.03E3 / T )
    #
    # [Jacobson, 1999]
    #    KO0  = 2.40E + 4
    #    KO1  = 3.70E + 5 * EXP( - 18.56 * ((298.15 / T) - 1.0))
    #    KO2  = 1.50E + 9 * EXP( - 17.72 * ((298.15 / T) - 1.0))
    #
    # Rate constants from Jacobson is larger than those of Jacob
    # and results in faster conversion from S(IV) to S(VI)
    # We choose Jacob 1) 2) and Jacobson 0) here
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    ko0 = 2.40e + 4_fp
    ko1 = 3.49e + 12_fp * exp( - 4.83e + 3_fp / t )
    ko2 = 7.32e + 14_fp * exp( - 4.03e + 3_fp / t )

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # H2O2 equilibrium reaction
    # O3 + H2O = O3.H2O
    #  HCO3  = 1.13E - 2 * EXP( 8.51 * (298.15 / T - 1.0) ), SP
    #  HCO3  = 1.13E - 2 * EXP( 7.72 * (298.15 / T - 1.0) ), Jacobson
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # Calculate Henry"s Law constant for atmospheric temperature
    hco3  = 1.13e - 2_fp * exp( 8.51e + 0_fp * ( 298.15e + 0_fp / t - 1.0e + 0_fp ) )

    xo3g  = 1.0e + 0_fp / ( 1.0e + 0_fp + ( hco3 * r * t * lwc ) )
    #O3aq  = HCO3 * XO3g * O3 * P

    # Conversion rate from SO2 to SO4 via reaction with O3
    kaqo3 = (ko0 * xso2aq + ko1 * xhso3 + ko2 * xso3) * fhcso2 * xso2g * p * hco3 * xo3g * lwc * r * t   # [v / v / s]

    #(qjc, 04 / 10 / 16)
    kaqo3_1 = ko2 * xso3 * fhcso2 * xso2g * p * hco3 * xo3g * lwc * r * t   # [v / v / s]

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Metal (Fe, Mn) catalyzed O2 oxidation (bec, 7 / 12 / 04)
    # R = d[S(VI)] / dt = 750 * [Mn(II)] * [S(IV)] + 2600 * [Fe(III)] * [S(IV)] +
    #               1.0d10 * [Mn(II)] * [Fe(III)] * [S(IV)]
    # from Seinfeld and Pandis, 1998 pg. 371
    # S(IV) = HFCSO2 * XSO2 * P * [SO2]
    # R = KaqO2 * [SO2] (v / v / s)
    # KaqO2 = FHCSO2 * XSO2g * P *
    #        ((750 * [Mn(II)]) + (2600[Fe(III)]) + (1.0d10 * [Mn(II)] * [Fe(III)]))
    # in units of [M / s]
    # KaqO2 = FHCSO2 * XSO2g * P *
    #        ((750 * [Mn(II)]) + (2600[Fe(III)]) + (1.0d10 * [Mn(II)] * [Fe(III)])) *
    #        LWC * R * T / P
    # in units of [v / v / s]
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # Conversion rate from SO2 to SO4 via reaction with O2 (met cat)
    kaqo2 = fhcso2 * xso2g * ( (750e + 0_fp * mnii ) + ( 2600e + 0_fp * feiii ) + (1e + 10_fp * mnii * feiii ) ) * lwc * r * t   # [s - 1]

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #  Aqueous reactions of SO2 with HCHO
    #     HSO3 - + HCHO(aq) = > HMS + OH - (1)
    #     SO3 -  - + HCHO(aq) = > HMS                 (2)
    #
    #     NOTE:
    #     [Boyce and Hoffman, 1984]
    #        KHCHO1  = 7.9E2 * EXP( - 16.435 * ((298.15 / T) - 1.0))
    #        KHCHO2  = 2.5E7 * EXP( - 6.037 * ((298.15 / T) - 1.0))
    #
    #
    #  Aqueous reaction of HMS with OH - 
    #    HMS + OH - = > HCHO(aq) + SO3 -  - (3)
    #
    #     NOTE: unclear where B (E / R) value in Seinfeld and Pandis from,
    #     but close to Deister. Using Seinfeld and Pandis value for now
    #     [Deister et al., 1986]
    #        KHMS    = 3.6E3 * EXP( - 22.027 * ((298.15 / T) - 1.0))
    #     [Seinfeld and Pandis, 2016; Munger et al., 1986]
    #        KHMS    = 3.6E3 * EXP( - 15.09 * ((298.15 / T) - 1.0))
    #
    #
    #  Aqueous reaction of HMS with OH(aq)
    #    HMS + OH(aq) = (SO2, O2, HO2) = > 2SO4 -  - + HCHO + O2 + 3H + + 2H2O  (4)
    #
    #    NOTE: O2, SO2, and HO2 particpate in the stoichiometry but not kinetics.
    #          Assume steady state for sulfur radicals and the following reaction chain:
    #            HMS + OH(aq) = (O2) = > SO5 - + HCHO + H2O [Olsen and Fessenden, 1992]
    #            HO2 < = > H + + O2 - [Jacob, 1986]
    #            SO5 - + O2 - = (H2O) = > HSO5 - + OH - + O2   [Jacob, 1986]
    #            SO2(aq) < = > HSO3 - + H +
    #            H + + OH - < = > H2O
    #            HSO5 - + HSO3 - = > 2SO4 -  - + 2H +          [Jacob, 1986]
    #       Instead of assuming Henry"s law for OH, use the parameter from
    #       Jacob et al, 2005 that relates gas phase OH to aqueous phase OH
    #       accounting for the HO2(aq) / O2 - cylcing in cloud droplets:
    #        dOH = 1E - 19 [M cm^3 molec^ - 1]
    #     [Olson and Fessenden, 1992]
    #        KHMS2    = 6.2E8 * EXP( - 5.03 * ((298.15 / T) - 1.0))
    #
    #
    # (jmm, 06 / 28 / 18)
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    khcho1 = 7.9e + 2_fp * exp( - 16.44e + 0_fp# L / mol / s
         * ( 298.15e + 0_fp / t - 1.0e + 0_fp ) )
    khcho2 = 2.5e + 7_fp * exp( - 6.04e + 0_fp# L / mol / s
         * ( 298.15e + 0_fp / t - 1.0e + 0_fp ) )
    khms   = 3.6e + 3_fp * exp( - 15.09e + 0_fp# L / mol / s
         * ( 298.15e + 0_fp / t - 1.0e + 0_fp ) )
    khms2  = 2.65e + 8_fp * exp( - 5.03e + 0_fp# L / mol / s
         * ( 298.15e + 0_fp / t - 1.0e + 0_fp ) )

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # HCHO equilibrium reaction
    # HCHO(aq) + H2O   = HCH(OH)2
    #
    # Reaction rate dependent on Temperature is given
    #   H = A exp ( B (T. / T - 1) )
    #
    # For equilibrium reactions of HCHO
    #                             Ah1       Bh1
    #  Sienfeld and Pandis      2.53E3    13.48  [2016]
    #
    # (jmm, 06 / 15 / 18)
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    khc1 = 2.53e + 3_fp * exp( 13.48e + 0_fp * ( 298.15e + 0_fp / t - 1.0e + 0_fp ) )

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # H2O equilibrium reaction
    #  H2O   = H + + OH - 
    #
    # Reaction rate dependent on Temperature is given
    #   H = A exp ( B (T. / T - 1) )
    #
    # For equilibrium reactions of HCHO
    #                             Ah1       Bh1
    #  Sienfeld and Pandis       1E - 14 - 22.51  [2016]
    #
    # (jmm, 06 / 15 / 18)
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    kw1 = 1e-14_fp * exp( - 22.51e + 0_fp *  ( 298.15e + 0_fp / t - 1.0e + 0_fp ) )

    # Henry"s constant [mol / l - atm] and Effective Henry"s constant for HCHO
    # [Seinfeld and Pandis, 2016]
    # HCHCHO  = 2.5 * EXP( 21.6e + 0_fp * ( 298.15e + 0_fp / T - 1.0e + 0_fp) )
    # (jmm, - 6 / 15 / 18)
    hchcho  = 2.5e + 0_fp * exp( 21.6e + 0_fp *  (298.15e + 0_fp / t - 1.0e + 0_fp) )
    fhchcho = hchcho * (1.0e + 0_fp + khc1 )

    xhchog  = 1.0e + 0_fp / ( 1.0e + 0_fp + ( fhchcho * r * t * lwc ) )


    # Conversion rate from SO2 to HMS via reaction with HCHO
    # (jmm, 06 / 15 / 18)
    kaqhcho = (khcho1 * xhso3 + khcho2 * xso3) * fhcso2 * xso2g * p * hchcho * xhchog * lwc * r * t    # [v / v / s]

    # Conversion rate from HMS to SO2 via reaction with OH - 
    # (jmm, 06 / 15 / 18; MSL 1 / 18 / 22)
    kaqhms = khms * ( kw1 / hplus ) * 0.7e0_fp  # 70% scavenged by clouds; units [1 / s]

    # Conversion rate from HMS to SO42 - HCHO via reaction with OH(aq)
    # (jmm, 06 / 28 / 18; MSL, 01 / 14 / 22)
    kaqhms2 = khms2 * doh * 0.7e0_fp # 70% scavenged by clouds; units [cm3 / mcl / s]

  end
#EOC

  subroutine set_2r_cld( t, lwc, fc, hplus, cnvfac, p, so2, h2o2, kaqh2o2 )

    real(fp), parameter   :: r = 0.08205e + 0_fp
    real(fp) :: kh2o2, ks1, ks2, t, xso2aq, lwc, fc
    real(fp) :: xhso3, xso3, hcso2, hplus, fhcso2
    real(fp) :: xso2g, kh1, hch2o2, fhch2o2, xh2o2g
    real(fp) :: kaqh2o2
    real(fp) :: so2, h2o2, a, b, kab, cnvfac, p

    # [Jacob, 1986]
    kh2o2  = 6.31e + 14_fp * exp( - 4.76e + 3_fp / t )
    ks1    = 1.30e - 2_fp * exp( 6.75e + 0_fp * ( 298.15e + 0_fp / t - 1.0e + 0_fp ) )
    ks2    = 6.31e - 8_fp * exp( 5.05e + 0_fp * ( 298.15e + 0_fp / t - 1.0e + 0_fp ) )

    # SIV Fraction
    xso2aq = 1.0e + 0_fp / (1.0e + 0_fp + ks1 / hplus + ks1 * ks2 / (hplus * hplus))
    xhso3  = 1.0e + 0_fp / (1.0e + 0_fp + hplus / ks1 + ks2 / hplus)
    xso3   = 1.0e + 0_fp / (1.0e + 0_fp + hplus / ks2 + hplus * hplus / (ks1 * ks2))

    # Henry"s constant [mol / l - atm] and Effective Henry"s constant for SO2
    hcso2  = 1.22e + 0_fp * exp( 10.55e + 0_fp * ( 298.15e + 0_fp / t - 1.0e + 0_fp) )
    fhcso2 = hcso2 * (1.0e + 0_fp + (ks1 / hplus) + (ks1 * ks2 / (hplus * hplus)))

    xso2g  = 1.0e + 0_fp / ( 1.0e + 0_fp + ( fhcso2 * r * t * lwc ) )
    kh1 = 2.20e - 12_fp * exp( - 12.52e + 0_fp * ( 298.15e + 0_fp / t - 1.0e + 0_fp ) )

    # [Jacobson, 1999]
    hch2o2  = 7.45e + 4_fp * exp( 22.21e + 0_fp * (298.15e + 0_fp / t - 1.0e + 0_fp) )
    fhch2o2 = hch2o2 * (1.0e + 0_fp + (kh1 / hplus))

    xh2o2g  = 1.0e + 0_fp / ( 1.0e + 0_fp + ( fhch2o2 * r * t * lwc ) )

#    KaqH2O2  = kh2o2 * Ks1 * FHCH2O2 * HCSO2 * XH2O2g * XSO2ga = so2
    b = h2o2

    kab = kh2o2 * ks1 * fhch2o2 * hcso2 * xh2o2g * xso2g * p * lwc * r * t * cnvfac # cm2 / mcl / s

  end

  function cloudhet2r( a, b, fc, kab )  result( kx )
#
# #USES
#
#
# #INPUT PARAMETERS:
#
    real(fp), intent(in) :: fc         # Cloud Fraction [0 - 1]
    real(fp), intent(in) :: a, b       # Reactant Abundances
    real(fp), intent(in) :: kab        # Bimolecular Rate Constant
#
# #RETURN VALUE:
#
    real(fp)             :: kx # Grid - average loss frequency, cm3 / mcl / s
#
# #REVISION HISTORY:
#  23 Aug 2018 - C. D. Holmes - Initial version
#  15 May 2021 - M. Long - Revision for two reactants
#  See https: / / github.com / geoschem / geos - chem for complete history
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
    # Residence time of air in clouds, s
    real(fp), parameter :: tauc = 3600.0_fp
#
# #LOCAL VARIABLES:
#
    real(fp) :: kao, kbo, ff, denom, term1, term2
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#

    # Ratio of volume inside to outside cloud
    # FF has a range [0, + inf], so cap it at 1e30
    ff = safediv( fc, ( 1.0_fp - fc ), 1e30_fp )
    ff = min( ff, 1.0e30_fp )

    # Avoid div by zero for the TAUC / FF term
    term1 = 0.0_fp
    if ( ff > 0.0_fp ) term1 = tauc / ff
 end

    # Compute KAO and avoid div by zero
    #             term 1      term 2
    # KAO = 1 / ( (TAUC / FF) + ( 1 / (FC * KAB * B) ) ), units: 1 / s
    denom = fc * kab * b
    term2 = safediv( 1.0_fp, denom, 0.0_fp )
    denom = term1 + term2
    kao   = safediv( 1.0_fp, denom, 0.0_fp )

    # Compute KBO and avoid div by zero
    #             term 1      term 2
    # KBO = 1 / ( (TAUC / FF) + ( 1 / (FC * KAB * A) ) ), units: 1 / s
    denom = fc * kab * a
    term2 = safediv( 1.0_fp, denom, 0.0_fp )
    denom = term1 + term2
    kbo   = safediv( 1.0_fp, denom, 0.0_fp )

    if ( kao * a < = kbo * b ) 
 end
       #KX = KAO / B
       kx = safediv( kao, b, 0.0_fp )
    else
       #KX = KBO / A
       kx = safediv( kbo, a, 0.0_fp )
    end

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: CloudHet1R
#
# #DESCRIPTION: Function CloudHet calculates the loss frequency (1 / s) of gas
#  species due to heterogeneous chemistry on clouds in a partially cloudy grid
#  cell. The function uses the "entrainment limited uptake" equations of
#  Holmes et al. (2019). Both liquid and ice water clouds are treated.
#
#  For gasses that are that are consumed in multiple aqueous reactions with
#  different products, CloudHet can provide the loss frequency for each reaction
#  branch using the optional branch ratios (branchLiq, branchIce) as arguments.
#
#  Holmes, C.D., Bertram, T. H., Confer, K. L., Ronan, A. C., Wirks, C. K.,
#    Graham, K. A., Shah, V. (2019) The role of clouds in the tropospheric
#    NOx cycle: a new modeling approach for cloud chemistry and its global
#    implications, Geophys. Res. Lett. 46, 4980 - 4990,
#    https: / / doi.org / 10.1029 / 2019GL081990
#\\
#\\
# #INTERFACE:
#
  function cloudhet1r( fc, rate ) result( khet )
#
# #USES:
#
#
# #INPUT PARAMETERS:
#
    real(fp), intent(in) :: fc     # Cloud Fraction [0 - 1]
    real(fp), intent(in) :: rate   # 1st order reaction rate (1 / s)
#
# #RETURN VALUE:
#
    real(fp)             :: khet   # Grid - average loss frequency, 1 / s
#
# #REVISION HISTORY:
#  23 Aug 2018 - C. D. Holmes - Initial version
#  25 May 2021 - M. S. Long - Modified for 1st order aqueous reaction
#                             where diffusion is not limiting
#  See https: / / github.com / geoschem / geos - chem for complete history
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
#
# #DEFINED PARAMETERS:
#
    # Residence time of air in clouds, s
    real(fp), parameter :: tauc = 3600.0_fp
#
# #LOCAL VARIABLES:
#
    real(fp) :: ki, gam
    real(fp) :: kk, ff, xx
#
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#
    # If cloud fraction < 0.0001 (0.01%) or there is zero cloud surface area,
    #  return zero uptake
    if ( ( fc < 0.0001_fp ) || ( rate < = 0.0_fp ) ) 
 end
       khet = 0.0_fp
       return
    end

    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    # Loss frequency inside cloud
    #
    # Assume both water and ice phases are inside the same cloud, so mass
    # transport to both phases works in parallel (additive)
    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 

    # initialize loss, 1 / s
    ki  = rate # total loss rate of a gas in cloud

    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
    # Grid - average loss frequency
    #
    # EXACT expression for entrainment - limited uptake
    # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 

    # Ratio (in cloud) of heterogeneous loss to detrainment, s / s
    kk = ki * tauc

    # Ratio of volume inside to outside cloud
    # ff has a range [0, + inf], so cap it at 1e30
    ff = safediv( fc, (1e0_fp - fc), 1e30_fp )
    ff = min( ff, 1.0e30_fp )

    # Ratio of mass inside to outside cloud
    # xx has range [0, + inf], but ff is capped at 1e30, so shouldn"t overflow
    xx =     ( ff - kk - 1.0_fp       ) / 2.0_fp + sqrt( 1.0_fp    + ff * ff     + kk * kk                   + 2.0_fp * ff + 2.0_fp * kk - 2.0_fp * ff * kk ) / 2.0_fp

    # Overall heterogeneous loss rate, grid average, 1 / s
    # kHet = kI * xx / ( 1d0 + xx )
    #  Since the expression ( xx / (1 + xx) ) may behave badly when xx>>1,
    #  use the equivalent 1 / (1 + 1 / x) with an upper bound on 1 / x
    khet = ki / ( 1e0_fp + safediv( 1e0_fp, xx, 1e30_fp ) )

  end
#EOC
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: fullchem_InitSulfurChem
#
# #DESCRIPTION: Stores species indices in module variables for fast lookup.
#\\
#\\
# #INTERFACE:
#
  subroutine fullchem_initsulfurchem( rc )
#
# #USES:
#
#
# #INPUT / OUTPUT PARAMETERS:
#
    integer,        intent(out) :: rc          # Success or failure?
#
# #REVISION HISTORY:
#  02 Jun 2000 - R. Yantosca - Initial version
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

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # init_sulfate begins here#
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # Initialize
    rc       = gc_success
    errmsg   = ""
    thisloc  = " - > at fullchem_initsulfurchem (in kpp / fullchem / fullchem_sulfurchemfuncs.f90"

    # Define flags for species ID"s
    id_acta   = ind_( "acta"   )
    id_ch2o   = ind_( "ch2o"   )
    id_dms    = ind_( "dms"    )
    id_dst1   = ind_( "dst1"   )
    id_dst2   = ind_( "dst2"   )
    id_dst3   = ind_( "dst3"   )
    id_dst4   = ind_( "dst4"   )
    id_h2o2   = ind_( "h2o2"   )
    id_hcl    = ind_( "hcl"    )
    id_hcooh  = ind_( "hcooh"  )
    id_hms    = ind_( "hms"    )
    id_hno3   = ind_( "hno3"   )
    id_msa    = ind_( "msa"    )
    id_nh3    = ind_( "nh3"    )
    id_nh4    = ind_( "nh4"    )
    id_nit    = ind_( "nit"    )
    id_nits   = ind_( "nits"   )
    id_o3     = ind_( "o3"     )
    id_oh     = ind_( "oh"     )
    id_pfe    = ind_( "pfe"    )
    id_sala   = ind_( "sala"   )  # Sea salt aerosol     (fine mode  )
    id_salaal = ind_( "salaal" )  # Sea salt alkalinity  (fine mode  )
    id_salacl = ind_( "salacl" )  # Cl - on sea salt      (fine mode  )
    id_salc   = ind_( "salc"   )  # Sea salt aerosol     (coarse mode)
    id_salcal = ind_( "salcal" )  # Sea salt alkalinity  (coarse mode)
    id_salccl = ind_( "salccl" )  # Cl - on sea salt      (coarse mode)
    id_so2    = ind_( "so2"    )
    id_so4    = ind_( "so4"    )
    id_so4s   = ind_( "so4s"   )

  end
#EOC
end
