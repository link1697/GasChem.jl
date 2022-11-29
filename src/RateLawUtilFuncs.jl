temp = 288.15
#TEMP            = State_Met%T(I,J,L)
INV_TEMP        = 1.0    / temp
temp_over_K300  = temp     / 300.0 
k300_over_temp  = 300.0  / temp
SR_TEMP         = sqrt( temp )
NUMDEN = 1

# FOUR_R_T        = 4.0  * CON_R    * temp
# FOUR_RGASLATM_T = 4.0  * RGASLATM * temp
# EIGHT_RSTARG_T  = 8.0  * RSTARG   * temp

# #  Relative humidity quantities
# CONSEXP         = 17.2693882  * (temp - 273.16 ) / (temp - 35.86 )
# VPRESH2O        = CONSVAP * EXP( CONSEXP ) / temp
# RELHUM          = ( H2O / VPRESH2O ) * 100 

# define based on Set_Kpp_GridBox_Values
# https://github.com/geoschem/geos-chem/blob/0f85431c4beedb0af8325b54d789f910cd3d6357/GeosCore/fullchem_mod.F90

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#                  geos - chem global chemical transport model                  #
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOP
#
# #IROUTINE: rateLawUtilFuncs
#
# #DESCRIPTION: Provides common functions for computing reaction rates.
#\\
#\\
# #INTERFACE:
#
# module ratelawutilfuncs
#
# #USES:
#
#   public
#
# #PRIVATE MEMBER FUNCTIONS:
#
#   private :: coth
#
# #DEFINED PARAMETERS:
#
  # Minimum heterogeneous chemistry lifetime and reaction rate
#   real(dp), private, parameter :: het_min_lif e = 1.0e - 3 
#  end
#   real(dp), private, parameter :: het_min_rate = 1.0  / het_min_life
#EOP
# -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#BOC
# contains

  ##########################################################################
  ######                   ARRHENIUS FUNCTIONS                         #####
  ##########################################################################

function GCARR_ab( a0, b0 ) 
    # Arrhenius function, skipping computation of EXP( c0 / T ),
    # which evaluates to 1 when c0 = 0.0  This avoids excess CPU
    # cycles. (bmy, 12 / 18 / 20)
    #
   #  real(dp), intent(in) :: a0, b0
   #  real(dp)             :: k
    #
    k = a0 * k300_over_temp^b0
    return k
end

  function GCARR_ac( a0, c0 )  
    # Arrhenius function, skipping computation of ( 300 / T )^b0,
    # which evaluates to 1 when b0 = 0.0  This avoids excess CPU
    # cycles (bmy, 12 / 18 / 20)
    #
   #  real(dp), intent(in) :: a0, c0
   #  real(dp)             :: k
    #

    k = a0 * exp( c0 / temp )
    return k
  end

  function GCARR_abc( a0, b0, c0 )  
    # Arrhenius function, using all 3 terms.
    # Use this when a0, b0, c0 are all nonzero.
    #
    # real(dp), intent(in) :: a0, b0, c0
    # real(dp)             :: k
    #
    k = a0 * exp( c0 / temp ) * k300_over_temp^b0
    return k
  end

  ##########################################################################
  ######         COMMON FUNCTIONS FOR COMPUTING UPTAKE RATES           #####
  ##########################################################################

  function ars_l1k( area, radius, Γ, srmw )  
    #
    # Calculates the 1st - order loss rate of species on wet aerosol surface.
    #
    # real(dp), intent(in) :: area, radius, Γ, srmw
    # real(dp)             :: k,    dfkg
    #
    # If Γ or radius is very small, set rate to zero and return
    if ( Γ < 1.0e-30  || radius < 1.0e-30  ) 
       k = 0.0 
       return k
    end
    #
    # DFKG = Gas phase diffusion coeff [cm2 / s] (order of 0.1)
    dfkg = ( 9.45e + 17  / numden ) * sr_temp * sqrt( 3.472e - 2  + 1.0  / ( srmw * srmw ) )
    #
    # Compute ArsL1k according to the formula listed above
    k = area / ( (radius / dfkg) + 2.749064e - 4  * srmw / (Γ * sr_temp) )
    return k
  end

  function kiir1ltd( concgas, conceduct, kisource ) result( kii )
    #
    # Determine removal rates for both species in an uptake reaction.
    # - Assume that the 1st reactant (concGas) is limiting.
    # - Assume that the 2nd reactant (concEduct) is "abundant".
    # - Calculate the overall rate (kII) based only on the uptake
    #   rate of the first reactant.
    # NOTE: Rewritten for computational efficiency (bmy, 5 / 13 / 21)
    #
    # real(dp), intent(in) :: concgas, conceduct, kisource
    # real(dp)             :: kigas,   kieduct,   lifea,   lifeb, kii
    #
    kigas   = 0.0 
    kieduct = 0.0 
    kii     = 0.0 
    #
    # Prevent div by zero.  Now use 1.0 as the error trap for concEduct.
    # 100 and 1e-8 (the previous criteria) were too large and too small,
    # respectively.  See https: / / github.com / geoschem / geos - chem / issues / 1115.0
    # -  - Seb Eastham, Bob Yantosca (09 Feb 2022)
    if ( conceduct < 1.0                               ) return
    if ( ! is_safediv( concgas * kisource, conceduct ) ) return
    #
    # Compute rates
    kigas   = kisource
    kieduct = kigas * concgas / conceduct
    kii     = kigas           / conceduct
    #
    # Enforce a minimum lifetime?
    if ( kigas > 0.0  ) 
       #
       # Calculate lifetime of each reactant against removal
       lifea = safediv( 1.0 , kigas,   0.0  )
 end
       lifeb = safediv( 1.0 , kieduct, 0.0  )
 end
       #
       # Check if either lifetime is "too short"
       if ( ( lifea < lifeb ) && ( lifea < het_min_life ) ) 
          kii = safediv( het_min_rate, conceduct, 0.0  )
       elseif ( lifeb < het_min_life ) 
          kii = safediv( het_min_rate, concgas, 0.0  )
       end
    end
  end

#   function cloudhet( h, srmw, gamliq, gamice, brliq, brice ) result( khet )
#     #
#     # Function CloudHet calculates the loss frequency (1 / s) of gas species
#     # due to heterogeneous chemistry on clouds in a partially cloudy grid
#     # cell. The function uses the "entrainment limited uptake" equations of
#     # Holmes et al. (2019). Both liquid and ice water clouds are treated.
#     #
#     # For gasses that are that are consumed in multiple aqueous reactions
#     # with different products, CloudHet can provide the loss frequency for
#     # each reaction branch using the branch ratios (branchLiq, branchIce).
#     #
#     # Reference:
#     # Holmes, C.D., Bertram, T. H., Confer, K. L., Ronan, A. C., Wirks,
#     #   C. K., Graham, K. A., Shah, V. (2019) The role of clouds in the
#     #   tropospheric NOx cycle: a new modeling approach for cloud chemistry
#     #   and its global implications, Geophys. Res. Lett. 46, 4980 - 4990,
#     #   https: / / doi.org / 10.1029 / 2019GL081990
#     #
#     type(hetstate), intent(in) :: h              # Hetchem State object
#     real(dp),       intent(in) :: srmw           # SQRT( mol wt [g / mole] )
#     real(dp),       intent(in) :: gamliq         # Rxn prob, liquid [1]
#     real(dp),       intent(in) :: gamice         # Rxn prob, ice [1]
#     real(dp),       intent(in) :: brliq          # Frac of reactant consumed
#     real(dp),       intent(in) :: brice          #  in liqice branches [0 - 1]
#     real(dp)                   :: khet           # Grid - avg loss frequency [1 / s]
#     #
#     real(dp),       parameter  :: tauc = 3600.0 
#     real(dp)                   :: ki, gam, rd, area
#     real(dp)                   :: kk, ff, xx, branch, kib, ktmp

#     # If cloud fraction < 0.0001 (0.01%) or there is zero cloud surface
#     # area,  return zero uptake
#     if ( ( h%cldfr < 0.0001  ) || ( h%aliq + h%aice <= 0.0  ) ) 
#  end
#        khet = 0.0 
#        return
#     end

#     # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#     # Loss frequency inside cloud
#     #
#     # Assume both water and ice phases are inside the same cloud, so mass
#     # transport to both phases works in parallel (additive)
#     # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 

#     # initialize
#     ki   = 0.0 
#     kib  = 0.0 
#     ktmp = 0.0 
#     khet = 0.0 

#     # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#     # Liquid branch (skip if the liquid branching ratio is zero)
#     # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#     if ( brliq > 0.0  ) 

#        # Convert grid - average cloud condensate surface area density
#        # to in - cloud surface area density
#        area = safediv( h%aliq, h%cldfr, 0.0  )

#        # Skip if no area
#        if ( area > 0.0  ) 

#           # In - cloud loss frequency [1 / s]
#           ktmp = ars_l1k( area, h%rliq, gamliq, srmw )
#           ki   = ki + ktmp

#           # In - cloud loss frequency for liquid rxn branch [1 / s]
#           kib  = kib + ( ktmp * brliq )
#        end
#     end

#     # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#     # Ice branch (skip if the ice branching ratio is zero)
#     # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#     if ( brice > 0.0  ) 

#        # Convert grid - average cloud condensate surface area density
#        # to in - cloud surface area density
#        area = safediv( h%aice, h%cldfr, 0.0  )

#        # Skip if no area
#        if ( area > 0.0  ) 

#           # In - cloud loss frequency [1 / s]
#           ktmp = ars_l1k( area, h%rice, gamice, srmw )
#           ki   = ki + ktmp

#           # In - continue loud loss frequency for ice rxn branch [1 / s]
#           kib  = kib + ( ktmp * brice )
#        end
#     end

#     # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#     # Mean branch ratio for reaction of interest in cloud
#     # (averaged over ice and liquid)
#     #
#     # If the division can"t be done, set kHet = 0 and return
#     # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#     branch = safediv( kib, ki, 0.0  )
#     if ( ! branch > 0.0  ) 
#        khet = 0.0 
#        return
#     end

#     # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#     # Grid - average loss frequency
#     #
#     # EXACT expression for entrainment - limited uptake
#     # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 

#     # Ratio (in cloud) of heterogeneous loss to detrainment, s / s
#     kk = ki * tauc

#     # Ratio of volume inside to outside cloud
#     # ff has a range [0, + inf], so cap it at 1e30
#     ff = safediv( h%cldfr, h%clearfr, 1.0e + 30  )
#     ff = min( ff, 1.0e + 30  )

#     # Ratio of mass inside to outside cloud
#     # xx has range [0, + inf], but ff is capped at 1e30, so shouldn"t overflow.
#     xx =     ( ff - kk - 1.0        ) / 2.0  + sqrt( 1.0     + ff * ff     + kk * kk                   + 2.0  * ff + 2.0  * kk - 2.0  * ff * kk ) / 2.0 

#     # Do not let xx go negative, as this can cause numerical instability.
#     # See https: / / github.com / geoschem / geos - chem / issues / 1205
#     xx = max( xx, 0.0  )

#     # Overall heterogeneous loss rate, grid average, 1 / s
#     # kHet = kI * xx / ( 1d0 + xx )
#     #  Since the expression ( xx / (1 + xx) ) may behave badly when xx>>1,
#     #  use the equivalent 1 / (1 + 1 / x) with an upper bound on 1 / x
#     khet = ki / ( 1.0  + safediv( 1.0 , xx, 1.0e + 30  ) )

#     # Overall loss rate in a particular reaction branch, 1 / s
#     khet = khet * branch
#   end

# function cld_params( ad, cldf, frland, frocean, qi, ql, t, h )
#     #
#     # Returns ice and liquid cloud parameters (based on State_Met)
#     # for cloud particles.
#     #
#     # References:
#     #  Heymsfield, A. J., Winker, D., Avery, M., et al. (2014). Relationships
#     #   between ice water content and volume extinction coefficient from in
#     #   situ observations for temperatures from 0° to –86°C: implications
#     #   for spaceborne lidar retrievals. Journal of Applied Meteorology and
#     #   Climatology, 53(2), 479–505.0 https: / / doi.org / 10.1175 / JAMC - D - 13 - 087.1
#     #
#     #  Schmitt, C. G., Heymsfield, A. J. (2005). Total Surface Area Estimates
#     #   for Individual Ice Particles and Particle Populations. Journal of
#     #   Applied Meteorology, 44(4), 467–474.0 https: / / doi.org / 10.1175 / JAM2209.1
#     #
#    #  real(dp),       intent(in)    :: ad          # Air mass [kg]
#    #  real(dp),       intent(in)    :: cldf        # Cloud fraction [1]
#    #  real(dp),       intent(in)    :: frland      # Land fraction [1]
#    #  real(dp),       intent(in)    :: frocean     # Ocean fraction [1]
#    #  real(dp),       intent(in)    :: qi          # Ice mixing ratio [kg / kg]
#    #  real(dp),       intent(in)    :: ql          # Liquid mixing ratio [kg / kg]
#    #  real(dp),       intent(in)    :: t           # Temperature [K]
# #
# # #OUTPUT PARAMETERS:
# #
#    #  type(hetstate), intent(inout) :: h           # Hetchem State object
# #
# # #REMARKS:
# #EOP
# # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
# #BOC
# #
# # #DEFINED PARAMETERS:
# #
#     # Cloud droplet radius in continental warm clouds [cm]
#     real(dp), parameter :: cldr_cont = 6.0e-4 

#     # Cloud droplet radius in marine warm clouds [cm]
#     real(dp), parameter :: cldr_mari = 10.0e-4 

#     # Ice cloud droplet radius [cm]
#     real(dp), parameter :: cldr_ice  = 38.5e-4 

#     # Density of H2O liquid [kg / cm3]
#     real(dp), parameter :: dens_liq  = 0.001 

#     # Density of H2O ice [kg / cm3]
#     real(dp), parameter :: dens_ice  = 0.91e - 3 
# #
# # #LOCAL VARIABLES:
# #
#     real(dp) :: α, beta

#     # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#     # cld_params begins here#
#     # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#     # Exit if there is no cloud
#     if ( ( ql + qi  <= 0.0  ) || ( cldf  <= 0.0  ) ) 
#  end
#        h%rliq = cldr_cont
#        h%rice = cldr_ice
#        h%aliq = 0.0 
#        h%vliq = 0.0 
#        h%aice = 0.0 
#        h%vice = 0.0 
#        return
#     end

#     # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#     # In GC 12.0 and earlier, the liquid water volume was set to zero at
#     # temperatures colder than 258K and over land ice (Antarctica# Greenland). That was likely legacy code from GEOS - 4, which provided
#     # no information on cloud phase. As of GC 12.0, all met data sources
#     # provide explicit liquid and ice condensate amounts, so we use those
#     # as provided. (C.D. Holmes)
#     #
#     # Liquid water clouds
#     #
#     # Droplets are spheres, so
#     # Surface area = 3 * Volume / Radius
#     #
#     # Surface area density = Surface area / Grid volume
#     # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#     if ( frland > frocean ) 
#        h%rliq = cldr_cont      # Continental cloud droplet radius [cm]
#     else
#        h%rliq = cldr_mari      # Marine cloud droplet radius [cm]

#     end

#     # get the volume of cloud condensate [cm3(condensate) / cm3(air)]
#     # QL is [g / g]
#     h%vliq = ql * ad / dens_liq / h%vair
#     h%vice = qi * ad / dens_ice / h%vair
#     h%aliq = 3.0  * h%vliq / h%rliq

#     # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#     # Ice water clouds
#     #
#     # Surface area calculation requires information about ice crystal size
#     # and shape, which is a function of temperature. Use Heymsfield (2014)
#     # empirical relationships between temperature, effective radius,
#     # surface area and ice water content.
#     #
#     # Schmitt and Heymsfield (2005) found that ice surface area is about
#     # 9 times its cross - sectional area.
#     #
#     # For any shape,
#     #   Cross section area = pi * (Effective Radius)^2, so
#     #   Cross section area = 3 * Volume / ( 4 * Effective Radius ).
#     #
#     # Thus, for ice
#     #   Surface area = 9 * Cross section area
#     #                = 2.25 * 3 * Volume / Effective Radius
#     # (C.D. Holmes)
#     # -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 

#     # Heymsfield (2014) ice size parameters
#     if ( t < 202.0  )           # - 71 C
#        α = 83.3 
#        beta  = 0.0184 
#     elseif ( t < 217.0  )      # - 56 C
#        α = 9.1744e + 4 
#        beta  = 0.117 
#     else
#        α = 308.4 
#        beta  = 0.0152 
#     end

#     # Effective radius, cm
#     h%rice = 0.5  * α * exp( beta * ( t - 273.15  ) ) / 1e + 4 

#     # Ice surface area density, cm2 / cm3
#     h%aice = 3.0  * h%vice / h%rice * 2.25 

#   end

  ##########################################################################
  ######         COMMON FUNCTIONS FOR COMPUTING UPTAKE RATES           #####
  ##########################################################################

  function coth( x ) result( f_x )
    #
    # Hyperbolic cotangent = [1 + exp( - 2x)] / [1 - exp( - 2x)]
    #
    # real(dp), intent(in) :: x
    # real(dp)             :: y, f_x
    #
    y   = exp( - 2.0  * x )
    f_x = ( 1.0  + y ) / ( 1.0  - y )
  end

  function reactodiff_corr( radius, l ) result( corr )
    #
    # For x = radius / l, correction =  COTH( x ) - ( 1 / x )
    # Correction approaches 1 as x becomes large, corr(x>1000)~1
    # Correction approaches x / 3 as x goes towards 0
    #
    # real(dp), intent(in)  :: l, radius           # [cm] and [cm]
    # real(dp)              :: x, corr
    #
    x = radius / l
    if ( x > 1000.0  ) 
       corr = 1.0 
       return
    end
    if ( x < 0.1  ) 
       corr = x / 3.0 ;
       return
    end
    corr = coth(x) - ( 1.0  / x )
  end

  ##########################################################################
  ######   COMMON FUNCTIONS FOR ENFORCING SAFE NUMERICAL OPERATIONS    #####
  ##########################################################################

  function safediv( num, denom, alt ) result( quot )
    #
    # Performs "safe division", that is to prevent overflow, underlow,
    # NaN, or infinity errors.  An alternate value is returned if the
    # division cannot be performed.
    # real(dp), intent(in) :: num, denom, alt
    # real(dp)             :: ediff, quot
    #
    # Exponent difference (base 2)
    # For REAL * 8, max exponent = 1024 and min = - 1021
    ediff = exponent( num ) - exponent( denom )
#  end
    #
    if ( ediff > 1023 || denom == 0.0  ) 
 
      quot = alt
        
    elseif ( ediff < - 1020 ) 
       quot = 0.0 
    else
       quot = num / denom
    end
  end

  function is_safediv( num, denom )  
    #
    # Returns TRUE if a division can be performed safely.
    # real(dp), intent(in) :: num, denom
    # real(dp)             :: ediff
    #
    # Exponent difference (base 2)
    # For REAL * 8, max exponent = 1024 and min = - 1021
    safe  = true
    ediff = exponent( num ) - exponent( denom )
 
    #
    if ( ediff < - 1020 || ediff > 1023 || denom == 0.0  ) 
 
       safe = false
  
  end
  return safe

  function issafeexp( x ) result( safe )
    #
    # Returns TRUE if an exponential can be performed safely
    #
    # real(dp), intent(in) :: x
    #
    # Note EXP( 708 ) = 8.2e + 307 and EXP( - 708 ) = 3.3e-308, which are
    # very close to the maximum representable values at double precision.
    safe = ( abs( x ) < 709.0  )
  end

  function safeexp( x, alt ) result( y )
    #
    # Performs a "safe exponential", that is to prevent overflow, underflow,
    # underlow, NaN, or infinity errors when taking the value EXP( x ).  An
    # alternate value is returned if the exponential cannot be performed.
    #
    # real(dp), intent(in) :: x, alt
    # real(dp)             :: y
    #
    y = alt
    if ( abs( x ) < 709.0  ) y = exp( x )
 end
  end
#EOC
end