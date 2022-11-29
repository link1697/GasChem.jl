#USE gckpp_Parameters, ONLY: dp, NSPEC, NVAR, NFIX, NREACT


# Declaration of global variables

# ~~~ If you are using KPP within an OpenMP parallel environment,
# ~~~ then these variables must be declared THREADPRIVATE.  This means
# ~~~ that the compiler will make a private copy of these variables
# ~~~ (in stack memory) for each execution thread.  At the end of
# ~~~ the OpenMP parallel loop, these variables will be finalized,
# ~~~ and their memory deallocated.
# ~~~
# ~~~ NOTE: Because the OpenMP commands all begin with a comment
# ~~~ character, they will be ignored unless the code is compiled
# ~~~ with OpenMP parallelization turned on.
NSPEC = 1
NREACT = 1
TEMP = 1



# # C - Concentration of all species
#   REAL(kind=dp), TARGET :: C(NSPEC)
#   #$OMP THREADPRIVATE( C )
# # VAR - Concentrations of variable species (global)
#   REAL(kind=dp), POINTER :: VAR(:)
#   #$OMP THREADPRIVATE( VAR )
# # FIX - Concentrations of fixed species (global)
#   REAL(kind=dp), POINTER :: FIX(:)
#   #$OMP THREADPRIVATE( FIX )
# # RCONST - Rate constants (global)
#   REAL(kind=dp) :: RCONST(NREACT)
#   #$OMP THREADPRIVATE( RCONST )
# # TIME - Current integration time
#   REAL(kind=dp) :: TIME
#   #$OMP THREADPRIVATE( TIME )
# # SUN - Sunlight intensity between [0,1]
#   REAL(kind=dp) :: SUN
#   #$OMP THREADPRIVATE( SUN )
# # TEMP - Temperature
#   REAL(kind=dp) :: TEMP
#   #$OMP THREADPRIVATE( TEMP )

# # ~~~ If you are using KPP within an OpenMP parallel environment,
# # ~~~ these variables DO NOT need to be declared THREADPRIVATE.

# # TSTART - Integration start time
#   REAL(kind=dp) :: TSTART
# # TEND - Integration end time
#   REAL(kind=dp) :: TEND
# # DT - Integration step
#   REAL(kind=dp) :: DT
# # ATOL - Absolute tolerance
#   REAL(kind=dp) :: ATOL(NVAR)
# # RTOL - Relative tolerance
#   REAL(kind=dp) :: RTOL(NVAR)
# # STEPMIN - Lower bound for integration step
#   REAL(kind=dp) :: STEPMIN
# # STEPMAX - Upper bound for integration step
#   REAL(kind=dp) :: STEPMAX
# # CFACTOR - Conversion factor for concentration units
#   REAL(kind=dp) :: CFACTOR

# # INLINED global variable declarations

# # Inline common variables into gckpp_Global.F90
# #include "commonIncludeVars.H"

# # INLINED global variable declarations
 