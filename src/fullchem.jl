using Catalyst
using Unitful

export fullchem

# Add unit "ppb" to Unitful 
module MyUnits
using Unitful
@unit ppb "ppb" Number 1 / 1000000000 false
end
Unitful.register(MyUnits)


function fullchem()
    @variables A3O2(T)      = 0.22   #CH3CH2CH2OO; Primary RO2 from C3H8

    @variables ACET(T)      = 0.846   #CH3C(O)CH3; Acetone

    @variables ACTA(T)      = 0.138   #CH3C(O)OH; Acetic acid

    @variables AERI(T)      = 0.713   #I; Dissolved iodine

    @variables ALD2(T)      = 0.16   #CH3CHO; Acetaldehyde

    @variables ALK4(T)      = 0.602   #>= C4 alkanes

    @variables AONITA(T)    = 0.869   #Aerosol-phase organic nitrate from aromatic precursors

    @variables AROMRO2(T)   = 0.876   #generic peroxy radical from aromatic oxidation

    @variables AROMP4(T)    = 0.098   #Generic C4 product from aromatic oxidation

    @variables AROMP5(T)    = 0.107   #Generic C5 product from aromatic oxidation

    @variables ATO2(T)      = 0.768   #CH3C(O)CH2O2; RO2 from acetone

    @variables ATOOH(T)     = 0.994   #CH3C(O)CH2OOH; ATO2 peroxide

    @variables B3O2(T)      = 0.628   #CH3CH(OO)CH3; Secondary RO2 from C3H8

    @variables BALD(T)      = 0.987   #benzaldehyde and tolualdehyde

    @variables BENZ(T)      = 0.714   #C6H6; Benzene

    @variables BENZO(T)     = 0.889   #C6H5O radical

    @variables BENZO2(T)    = 0.879   #C6H5O2 radical

    @variables BENZP(T)     = 0.01   #hydroperoxide from BENZO2

    @variables Br(T)        = 0.396   #Br; Atomic bromine

    @variables Br2(T)       = 0.638   #Br2; Molecular bromine

    @variables BrCl(T)      = 0.368   #BrCl; Bromine chloride

    @variables BrNO2(T)     = 0.311   #BrNO2; Nitryl bromide

    @variables BrNO3(T)     = 0.684   #BrNO3; Bromine nitrate

    @variables BrO(T)       = 0.945   #BrO; Bromine monoxide

    @variables BRO2(T)      = 0.718   #C6H5O2 ; Peroxy radical from BENZ oxidation

    @variables BrSALA(T)    = 0.319   #Br; Fine sea salt bromine

    @variables BrSALC(T)    = 0.86   #Br; Course sea salt bromine

    @variables BZCO3(T)     = 0.579   #benzoylperoxy radical

    @variables BZCO3H(T)    = 0.829   #perbenzoic acid

    @variables BZPAN(T)     = 0.491   #peroxybenzoyl nitrate

    @variables C2H2(T)      = 0.51   #C2H2; Ethyne

    @variables C2H4(T)      = 0.714   #Ethylene

    @variables C2H6(T)      = 0.199   #C2H6; Ethane

    @variables C3H8(T)      = 0.654   #C3H8; Propane

    @variables C4HVP1(T)    = 0.481   #C4 hydroxy-vinyl-peroxy radicals from HPALDs

    @variables C4HVP2(T)    = 0.829   #C4 hydroxy-vinyl-peroxy radicals from HPALDs

    @variables CCl4(T)      = 0.55   #CCl4; Carbon tetrachloride

    @variables CFC11(T)     = 0.361   #CCl3F ; CFC-11, R-11, Freon 11

    @variables CFC12(T)     = 0.508   #CCl2F2; CFC-12, R-12, Freon 12

    @variables CFC113(T)    = 0.588   #C2Cl3F3; CFC-113, Freon 113

    @variables CFC114(T)    = 0.674   #C2Cl2F4; CFC-114, Freon 114

    @variables CFC115(T)    = 0.84   #C2ClF5; CFC-115, Freon 115

    @variables CH2Br2(T)    = 0.171   #CH3Br2; Dibromomethane

    @variables CH2Cl2(T)    = 0.414   #CH2Cl2; Dichloromethane

    @variables CH2I2(T)     = 0.32   #CH2I2; Diiodomethane

    @variables CH2IBr(T)    = 0.334   #CH2IBr; Bromoiodomethane

    @variables CH2ICl(T)    = 0.979   #CH2ICl; Chloroiodomethane

    @variables CH2O(T)      = 0.714   #CH2O; Formaldehyde

    @variables CH2OO(T)     = 0.516   #CH2OO; Criegee intermediate

    @variables CH3Br(T)     = 0.131   #CH3Br; Methyl bromide

    @variables CH3CCl3(T)   = 0.838   #CH3CCl3; Methyl chloroform

    @variables CH3CHOO(T)   = 0.93   #CH3CHOO; Criegee intermediate

    @variables CH3Cl(T)     = 0.937   #CH3Cl; Chloromethane

    @variables CH3I(T)      = 0.293   #CH3I; Methyl iodide

    @variables CH4(T)       = 0.218   #CH4; Methane

    @variables CHBr3(T)     = 0.271   #CHBr3; Tribromethane

    @variables CHCl3(T)     = 0.808   #CHCl3; Chloroform

    @variables Cl(T)        = 0.954   #Cl; Atomic chlorine

    @variables Cl2(T)       = 0.467   #Cl2; Molecular chlorine

    @variables Cl2O2(T)     = 0.73   #Cl2O2; Dichlorine dioxide

    @variables ClNO2(T)     = 0.323   #ClNO2; Nitryl chloride

    @variables ClNO3(T)     = 0.526   #ClONO2; Chlorine nitrate

    @variables ClO(T)       = 0.619   #ClO; Chlorine monoxide

    @variables ClOO(T)      = 0.323   #ClOO; Chlorine dioxide

    @variables CO(T)        = 0.224   #CO; Carbon monoxide

    @variables CO2(T)       = 0.17   #CO2; Carbon dioxide

    @variables CSL(T)       = 0.787   #cresols and xylols

    @variables DMS(T)       = 0.072   #(CH3)2S; Dimethylsulfide

    @variables EOH(T)       = 0.165   #C2H5OH; Ethanol

    @variables ETHLN(T)     = 0.793   #CHOCH2ONO2; Ethanal nitrate

    @variables ETHN(T)      = 0.446   #stable hydroxy-nitrooxy-ethane

    @variables ETHP(T)      = 0.66   #stable hydroxy-hydroperoxy-ethane

    @variables ETNO3(T)     = 0.086   #C2H5ONO2; Ethyl nitrate

    @variables ETO(T)       = 0.924   #hydroxy-alkoxy-ethane radical

    @variables ETOO(T)      = 0.185   #hydroxy-peroxy-ethane radical, formed from ethene + OH

    @variables ETO2(T)      = 0.568   #CH3CH2OO; Ethylperoxy radical

    @variables ETP(T)       = 0.893   #CH3CH2OOH; Ethylhydroperoxide

    @variables GLYC(T)      = 0.32   #HOCH2CHO; Glycoaldehyde

    @variables GLYX(T)      = 0.966   #CHOCHO; Glyoxal

    @variables H(T)         = 0.321   #H; Atomic hydrogen

    @variables H1211(T)     = 0.476   #CBrClF2; H-1211

    @variables H1301(T)     = 0.741   #CBrF3; H-1301

    @variables H2402(T)     = 0.789   #C2Br2F4; H-2402

    @variables H2O(T)       = 0.289   #H2O; Water vapor

    @variables H2O2(T)      = 0.308   #H2O2; Hydrogen peroxide

    @variables HAC(T)       = 0.243   #HOCH2C(O)CH3; Hydroxyacetone

    @variables HBr(T)       = 0.514   #HBr; Hypobromic acid

    @variables HC5A(T)      = 0.58   #C5H8O2; Isoprene-4,1-hydroxyaldehyde

    @variables HCFC123(T)   = 0.415   #C2HCl2F3; HCFC-123, R-123, Freon 123

    @variables HCFC141b(T)  = 0.702   #C(CH3)Cl2F; HCFC-141b, R-141b, Freon 141b

    @variables HCFC142b(T)  = 0.074   #C(CH3)ClF2; HCFC-142b, R-142b, Freon 142b

    @variables HCFC22(T)    = 0.831   #CHClF2 ; HCFC-22, R-22, Freon 22

    @variables HCl(T)       = 0.092   #HCl; Hydrochloric acid

    @variables HCOOH(T)     = 0.654   #HCOOH; Formic acid

    @variables HI(T)        = 0.232   #HI; Hydrogen iodide

    @variables HMHP(T)      = 0.203   #HOCH2OOH; Hydroxymethyl hydroperoxide

    @variables HMML(T)      = 0.346   #C4H6O3; Hydroxymethyl-methyl-a-lactone

    @variables HMS(T)       = 0.272   #HOCH2SO3-; hydroxymethanesulfonate

    @variables HNO2(T)      = 0.472   #HONO; Nitrous acid

    @variables HNO3(T)      = 0.121   #HNO3; Nitric acid

    @variables HNO4(T)      = 0.719   #HNO4; Pernitric acid

    @variables HO2(T)       = 0.281   #HO2; Hydroperoxyl radical

    @variables HOBr(T)      = 0.286   #HOBr; Hypobromous acid

    @variables HOCl(T)      = 0.3   #HOCl; Hypochlorous acid

    @variables HOI(T)       = 0.429   #HOI; Hypoiodous acid

    @variables HONIT(T)     = 0.589   #2nd gen monoterpene organic nitrate

    @variables HPALD1(T)    = 0.606   #O=CHC(CH3)=CHCH2OOH; d-4,1-C5-hydroperoxyaldehyde

    @variables HPALD1OO(T)  = 0.947   #peroxy radicals from HPALD1

    @variables HPALD2(T)    = 0.853   #HOOCH2C(CH3)=CHCH=O; d-1,4-C5-hydroperoxyaldehyde

    @variables HPALD2OO(T)  = 0.104   #peroxy radicals from HPALD2

    @variables HPALD3(T)    = 0.408   #O=CHC(CH3)OOHCH=CH2; b-2,1-C5-hydroperoxyaldehyde

    @variables HPALD4(T)    = 0.121   #CH2=C(CH3)CHOOHCH=O; b-3,4-C5-hydroperoxyaldehyde

    @variables HPETHNL(T)   = 0.918   #CHOCH2OOH; hydroperoxyethanal

    @variables I(T)         = 0.384   #I; Atmoic iodine

    @variables I2(T)        = 0.71   #I2; Molecular iodine

    @variables I2O2(T)      = 0.028   #I2O2; Diiodine dioxide

    @variables I2O3(T)      = 0.554   #I2O3; Diiodine trioxide

    @variables I2O4(T)      = 0.236   #I2O4; Diiodine tetraoxide

    @variables IBr(T)       = 0.311   #IBr; Iodine monobromide

    @variables ICHE(T)      = 0.905   #C5H8O3; Isoprene hydroxy-carbonyl-epoxides

    @variables ICHOO(T)     = 0.602   #peroxy radical from IEPOXD

    @variables ICl(T)       = 0.621   #ICl; Iodine monochloride

    @variables ICN(T)       = 0.984   #C5H7NO4; Lumped isoprene carbonyl nitrates

    @variables ICNOO(T)     = 0.859   #peroxy radicals from ICN

    @variables ICPDH(T)     = 0.068   #C5H10O5; Isoprene dihydroxy hydroperoxycarbonyl

    @variables IDC(T)       = 0.666   #C5H6O2; Lumped isoprene dicarbonyls

    @variables IDCHP(T)     = 0.855   #C5H8O5; Isoprene dicarbonyl hydroxy dihydroperoxide

    @variables IDHDP(T)     = 0.791   #C5H12O6; Isoprene dihydroxy dihydroperoxide

    @variables IDHNBOO(T)   = 0.551   #peroxy radicals from INPB

    @variables IDHNDOO1(T)  = 0.379   #peroxy radicals from INPD

    @variables IDHNDOO2(T)  = 0.534   #peroxy radicals from INPD

    @variables IDHPE(T)     = 0.505   #C5H10O5; Isoprene dihydroxy hydroperoxy epoxide

    @variables IDN(T)       = 0.93   #C5H8N2O6; Lumped isoprene dinitrates

    @variables IDNOO(T)     = 0.58   #peroxy radicals from IDN

    @variables IEPOXA(T)    = 0.758   #C5H10O3; trans-Beta isoprene epoxydiol

    @variables IEPOXAOO(T)  = 0.404   #peroxy radical from trans-Beta isoprene epoxydiol

    @variables IEPOXB(T)    = 0.165   #C5H10O3; cis-Beta isoprene epoxydiol

    @variables IEPOXBOO(T)  = 0.371   #peroxy radical from cis-Beta isoprene epoxydiol

    @variables IEPOXD(T)    = 0.507   #C5H10O3; Delta isoprene epoxydiol

    @variables IHN1(T)      = 0.595   #C5H9NO4; Isoprene-d-4-hydroxy-1-nitrate

    @variables IHN2(T)      = 0.64   #C5H9NO4; Isoprene-b-1-hydroxy-2-nitrate

    @variables IHN3(T)      = 0.691   #C5H9NO4; Isoprene-b-4-hydroxy-3-nitrate

    @variables IHN4(T)      = 0.111   #C5H9NO4; Isoprene-d-1-hydroxy-4-nitrate

    @variables IHOO1(T)     = 0.761   #peroxy radical from OH addition to isoprene at C1

    @variables IHOO4(T)     = 0.961   #peroxy radical from OH addition to isoprene at C4

    @variables IHPNBOO(T)   = 0.483   #peroxy radicals from INPB

    @variables IHPNDOO(T)   = 0.42   #peroxy radicals from INPD

    @variables IHPOO1(T)    = 0.96   #peroxy radical from ISOPOOH

    @variables IHPOO2(T)    = 0.776   #peroxy radical from ISOPOOH

    @variables IHPOO3(T)    = 0.62   #peroxy radical from ISOPOOH

    @variables INA(T)       = 0.693   #alkoxy radical from INO2D

    @variables INDIOL(T)    = 0.58   #Generic aerosol phase organonitrate hydrolysis product

    @variables INO(T)       = 0.326   #INO; Nitrosyl iodide

    @variables INO2B(T)     = 0.665   #beta-peroxy radicals from isoprene + NO3

    @variables INO2D(T)     = 0.938   #delta-peroxy radicals from isoprene + NO3

    @variables INPB(T)      = 0.206   #C5H9NO5; Lumped isoprene beta-hydroperoxy nitrates

    @variables INPD(T)      = 0.707   #C5H9NO5; Lumped isoprene delta-hydroperoxy nitrates

    @variables IO(T)        = 0.881   #IO; Iodine monoxide

    @variables IONITA(T)    = 0.587   #Aerosol-phase organic nitrate from isoprene precursors

    @variables IONO(T)      = 0.071   #IONO; Nitryl iodide

    @variables IONO2(T)     = 0.742   #IONO2; Iodine nitrate

    @variables IPRNO3(T)    = 0.757   #C3H8ONO2; Isopropyl nitrate

    @variables ISALA(T)     = 0.879   #I; Fine sea salt iodine

    @variables ISALC(T)     = 0.935   #I; Coarse sea salt iodine

    @variables ISOP(T)      = 0.187   #CH2=C(CH3)CH=CH2; Isoprene

    @variables ISOPNOO1(T)  = 0.054   #peroxy radicals from IHN2

    @variables ISOPNOO2(T)  = 0.717   #peroxy radicals from IHN3

    @variables ITCN(T)      = 0.628   #C5H9NO7; Lumped tetrafunctional isoprene carbonyl-nitrates

    @variables ITHN(T)      = 0.777   #C5H11NO7; Lumped tetrafunctional isoprene hydroxynitrates

    @variables KO2(T)       = 0.04   #RO2 from >3 ketones

    @variables LBRO2H(T)    = 0.188   #Dummy spc to track oxidation of BRO2 by HO2

    @variables LBRO2N(T)    = 0.665   #Dummy spc to track oxidation of BRO2 by NO

    @variables LIMO(T)      = 0.897   #C10H16; Limonene

    @variables LIMO2(T)     = 0.019   #RO2 from LIMO

    @variables LISOPOH(T)   = 0.739   #Dummy spc to track oxidation of ISOP by OH

    @variables LISOPNO3(T)  = 0.581   #Dummy spc to track oxidation of ISOP by NO3

    @variables LNRO2H(T)    = 0.514   #Dummy spc to track oxidation of NRO2 by HO2

    @variables LNRO2N(T)    = 0.998   #Dummy spc to track oxidation of NRO2 by NO

    @variables LTRO2H(T)    = 0.858   #Dummy spc to track oxidation of TRO2 by HO2

    @variables LTRO2N(T)    = 0.653   #Dummy spc to track oxidation of TRO2 by NO

    @variables LVOC(T)      = 0.876   #C5H14O5; Gas-phase low-volatility non-IEPOX product of ISOPOOH (RIP) oxidation

    @variables LVOCOA(T)    = 0.831   #C5H14O5; Aerosol-phase low-volatility non-IEPOX product of ISOPOOH (RIP) oxidation

    @variables LXRO2H(T)    = 0.82   #Dummy spc to track oxidation of XRO2 by HO2

    @variables LXRO2N(T)    = 0.266   #Dummy spc to track oxidation of XRO2 by NO

    @variables MACR(T)      = 0.935   #CH2=C(CH3)CHO; Methacrolein

    @variables MACR1OO(T)   = 0.054   #peroxyacyl radical from MACR + OH

    @variables MACR1OOH(T)  = 0.62   #CH2=C(CH3)C(O)OOH; Peracid from MACR

    @variables MACRNO2(T)   = 0.774   #Product of MCRHN + OH

    @variables MAP(T)       = 0.7   #CH3C(O)OOH; Peroxyacetic acid

    @variables MCO3(T)      = 0.326   #CH3C(O)OO; Peroxyacetyl radical

    @variables MCRDH(T)     = 0.69   #C4H8O3; Dihydroxy-MACR

    @variables MCRENOL(T)   = 0.555   #C4H6O2; Lumped enols from MVK/MACR

    @variables MCRHN(T)     = 0.809   #HOCH2C(ONO2)(CH3)CHO; Hydroxynitrate from MACR

    @variables MCRHNB(T)    = 0.705   #O2NOCH2C(OH)(CH3)CHO; Hydroxynitrate from MACR

    @variables MCRHP(T)     = 0.032   #HOCH2C(OOH)(CH3)CHO; Hydroxy-hydroperoxy-MACR

    @variables MCROHOO(T)   = 0.338   #peroxy radical from MACR + OH

    @variables MCT(T)       = 0.126   #methylcatechols

    @variables MEK(T)       = 0.238   #RC(O)R; Methyl ethyl ketone

    @variables MENO3(T)     = 0.669   #CH3ONO2; methyl nitrate

    @variables MGLY(T)      = 0.128   #CH3COCHO; Methylglyoxal

    @variables MO2(T)       = 0.195   #CH3O2; Methylperoxy radical

    @variables MOH(T)       = 0.061   #CH3OH; Methanol

    @variables MONITA(T)    = 0.265   #Aerosol-phase organic nitrate from monoterpene precursors

    @variables MONITS(T)    = 0.304   #Saturated 1st gen monoterpene organic nitrate

    @variables MONITU(T)    = 0.39   #Unsaturated 1st gen monoterpene organic nitrate

    @variables MP(T)        = 0.226   #CH3OOH; Methylhydroperoxide

    @variables MPAN(T)      = 0.137   #CH2=C(CH3)C(O)OONO2; Peroxymethacroyl nitrate (PMN)

    @variables MPN(T)       = 0.481   #CH3O2NO2; Methyl peroxy nitrate

    @variables MSA(T)       = 0.121   #CH4SO3; Methanesulfonic acid

    @variables MTPA(T)      = 0.282   #Lumped monoterpenes: a-pinene, b-pinene, sabinene, carene

    @variables MTPO(T)      = 0.051   #Other monoterpenes: Terpinene, terpinolene, myrcene, ocimene, other monoterpenes

    @variables MVK(T)       = 0.189   #CH2=CHC(=O)CH3; Methyl vinyl ketone

    @variables MVKDH(T)     = 0.367   #HOCH2CH2OHC(O)CH3; Dihydroxy-MVK

    @variables MVKHC(T)     = 0.866   #C4H6O3; MVK hydroxy-carbonyl

    @variables MVKHCB(T)    = 0.418   #C4H6O3; MVK hydroxy-carbonyl

    @variables MVKHP(T)     = 0.803   #C4H8O4; MVK hydroxy-hydroperoxide

    @variables MVKN(T)      = 0.57   #HOCH2CH(ONO2)C(=O)CH3; Hydroxynitrate from MVK

    @variables MVKOHOO(T)   = 0.446   #peroxy radical from MVK + OH

    @variables MVKPC(T)     = 0.685   #OCHCH(OOH)C(O)CH3; MVK hydroperoxy-carbonyl

    @variables N(T)         = 0.154   #N; Atomic nitrogen

    @variables N2O(T)       = 0.598   #N2O; Nitrous oxide

    @variables N2O5(T)      = 0.711   #N2O5; Dinitrogen pentoxide

    @variables NAP(T)       = 0.611   #C10H8; Naphthalene; IVOC surrogate

    @variables NIT(T)       = 0.296   #NIT; Fine mode inorganic nitrate

    @variables NITs(T)      = 0.364   #NITs; Coarse mode inorganic nitrate

    @variables NO(T)        = 0.687   #NO; Nitric oxide

    @variables NO2(T)       = 0.124   #NO2; Nitrogen dioxide

    @variables NO3(T)       = 0.205   #NO3; Nitrate radical

    @variables NPHEN(T)     = 0.589   #nitrophenols

    @variables NPRNO3(T)    = 0.503   #C3H8ONO2; n-propyl nitrate

    @variables NRO2(T)      = 0.198   #Peroxy radical from NAP oxidation

    @variables O(T)         = 0.52   #O(3P); Ground state atomic oxygen

    @variables O1D(T)       = 0.833   #O(1D); Excited atomic oxygen

    @variables O3(T)        = 0.312   #O3; Ozone

    @variables O3A(T)       = 0.473   #O3; Ozone in accum seasalt

    @variables O3C(T)       = 0.391   #O3; Ozone in coarse seasalt

    @variables OClO(T)      = 0.727   #OClO; Chlorine dioxide

    @variables OCS(T)       = 0.842   #COS; Carbonyl sulfide

    @variables OH(T)        = 0.439   #OH; Hydroxyl radical

    @variables OIO(T)       = 0.889   #OIO; Iodine dioxide

    @variables OLND(T)      = 0.266   #Monoterpene-derived NO3-alkene adduct

    @variables OLNN(T)      = 0.027   #Monoterpene-derived NO3 adduct

    @variables OTHRO2(T)    = 0.209   #Other C2 RO2 not from C2H6 oxidation

    @variables PAN(T)       = 0.704   #CH3C(O)OONO2; Peroxyacetylnitrate

    @variables PHEN(T)      = 0.494   #phenol

    @variables PIO2(T)      = 0.559   #RO2 from MTPA

    @variables PIP(T)       = 0.644   #Peroxides from MTPA

    @variables PO2(T)       = 0.728   #HOCH2CH(OO)CH3; RO2 from propene

    @variables PP(T)        = 0.697   #HOCH2CH(OOH)CH3; Peroxide from PO2

    @variables PPN(T)       = 0.443   #CH3CH2C(O)OONO2; Peroxypropionylnitrate

    @variables PRN1(T)      = 0.915   #O2NOCH2CH(OO)CH3; RO2 from propene + NO3

    @variables PROPNN(T)    = 0.421   #CH3C(=O)CH2ONO2; Propanone nitrate

    @variables PRPE(T)      = 0.219   #C3H6; >= C3 alkenes

    @variables PRPN(T)      = 0.276   #O2NOCH2CH(OOH)CH3; Peroxide from PRN1

    @variables PYAC(T)      = 0.54   #CH3COCOOH; Pyruvic acid

    @variables R4N1(T)      = 0.632   #RO2 from R4N2

    @variables R4N2(T)      = 0.89   #RO2NO; >= C4 alkylnitrates

    @variables R4O2(T)      = 0.02   #RO2 from ALK4

    @variables R4P(T)       = 0.118   #CH3CH2CH2CH2OOH; Peroxide from R4O2

    @variables RA3P(T)      = 0.822   #CH3CH2CH2OOH; Peroxide from A3O2

    @variables RB3P(T)      = 0.165   #CH3CH(OOH)CH3; Peroxide from B3O2

    @variables RCHO(T)      = 0.315   #CH3CH2CHO; >= C3 aldehydes

    @variables RCO3(T)      = 0.666   #CH3CH2C(O)OO; Peroxypropionyl radical

    @variables RIPA(T)      = 0.633   #HOCH2C(OOH)(CH3)CH=CH2; 1,2-ISOPOOH

    @variables RIPB(T)      = 0.152   #HOCH2C(OOH)(CH3)CH=CH2; 4,3-ISOPOOH

    @variables RIPC(T)      = 0.431   #C5H10O3; d(1,4)-ISOPOOH

    @variables RIPD(T)      = 0.823   #C5H10O3; d(4,1)-ISOPOOH

    @variables ROH(T)       = 0.283   #C3H7OH; > C2 alcohols

    @variables RP(T)        = 0.992   #CH3CH2C(O)OOH; Peroxide from RCO3

    @variables SALAAL(T)    = 0.964   #Accumulation mode seasalt aerosol alkalinity

    @variables SALCAL(T)    = 0.145   #Coarse mode seasalt aerosol alkalinity

    @variables SALACL(T)    = 0.676   #Cl; Fine chloride

    @variables SALCCL(T)    = 0.534   #Cl; Coarse chloride

    @variables SALASO2(T)   = 0.971   #SO2; Fine seasalt

    @variables SALCSO2(T)   = 0.027   #SO2; Coarse seasalt

    @variables SALASO3(T)   = 0.832   #SO3--; Fine seasalt

    @variables SALCSO3(T)   = 0.993   #SO3--; Coarse chloride

    @variables SO2(T)       = 0.764   #SO2; Sulfur dioxide

    @variables SO4(T)       = 0.098   #SO4; Sulfate

    @variables SO4s(T)      = 0.062   #SO4 on sea-salt; Sulfate

    @variables SOAGX(T)     = 0.649   #CHOCHO; Aerosol-phase glyoxal

    @variables SOAIE(T)     = 0.59   #C5H10O3; Aerosol-phase IEPOX

    @variables TOLU(T)      = 0.084   #C7H8; Toluene

    @variables TRO2(T)      = 0.817   #Peroxy radical from TOLU oxidation

    @variables XYLE(T)      = 0.436   #C8H10; Xylene

    @variables XRO2(T)      = 0.544   #Peroxy radical from XYLE oxidation

    @variables H2(T)        = 0.584   #H2; Molecular hydrogen

    @variables N2(T)        = 0.041   #N2; Molecular nitrogen

    @variables O2(T)        = 0.838   #O2; Molecular oxygen

    @variables RCOOH(T)     = 0.523   #C2H5C(O)OH; > C2 organic acids 

    rxs = [
        Reaction(GCARR_ac(3.00e-12, -1500.0e0), [O3, NO], [NO2, O2], [1, 1], [1, 1])                         GCARR_ac(3.00d-12, -1500.0d0);
        Reaction(GCARR_ac(1.70e-12, -940.0e0), [O3, OH], [HO2, O2], [1, 1], [1, 1])                         GCARR_ac(1.70d-12, -940.0d0);
        Reaction(GCARR_ac(1.00e-14, -490.0e0), [O3, HO2], [OH, O2, O2], [1, 1], [1, 1, 1])                    GCARR_ac(1.00d-14, -490.0d0);
        Reaction(GCARR_ac(1.20e-13, -2450.0e0), [O3, NO2], [O2, NO3], [1, 1], [1, 1])                        GCARR_ac(1.20d-13, -2450.0d0);
        Reaction(GCARR_ac(2.90e-16, -1000.0e0), [O3, MO2], [CH2O, HO2, O2], [1, 1], [1, 1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(1.80e-12, [OH, OH], [H2O, O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCJPLPR_aba(6.90e-31, 1.0e+00, 2.6e-11, 0.6e0), [#, M], [H2O2], [1, 1], [1])                        GCJPLPR_aba(6.90d-31, 1.0d+00, 2.6d-11, 0.6d0);
        Reaction(GCARR_ac(4.80e-11, 250.0e0), [OH, HO2], [H2O, O2], [1, 1], [1, 1])                        GCARR_ac(4.80d-11, 250.0d0);
        Reaction(1.80e-12, [OH, H2O2], [H2O, HO2], [1, 1], [1, 1])                      1.80d-12;
        Reaction(GCARR_ac(3.30e-12, 270.0e0), [HO2, NO], [OH, NO2], [1, 1], [1, 1])   #2013/02/12; JPL 10-6; BHH,JMAO,EAM
        Reaction(GC_HO2HO2_acac(3.00e-13, 460.0e0, 2.1e-33, 920.0e0), [HO2, HO2], [H2O2, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GC_OHCO_a(1.50e-13), [OH, CO], [HO2, CO2], [1, 1], [1, 1])   #2017/02/22; JPL 15-10; BHH,MJE
        Reaction(GCARR_ac(2.45e-12, -1775.0e0), [OH, CH4], [MO2, H2O], [1, 1], [1, 1])                       GCARR_ac(2.45d-12, -1775.0d0);
        Reaction(GC_RO2NO_B1_ac(2.80e-12, 300.0e0), [MO2, NO], [CH2O, HO2, NO2], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF
        Reaction(GC_RO2NO_A1_ac(2.80e-12, 300.0e0), [MO2, NO], [MENO3], [1, 1], [1])   #2019/05/10; Fisher2018; JAF
        Reaction(GCARR_abc(4.10e-13, 0.0e0, 750.0e0), [MO2, HO2], [MP, O2], [1, 1], [1, 1])                        GCARR_abc(4.10d-13, 0.0d0, 750.0d0);
        Reaction(GC_TBRANCH_1_acac(9.50e-14, 390.0e0, 2.62e1, -1130.0e0), [MO2, MO2], [MOH, CH2O, O2], [1, 1], [1, 1, 1])                GC_TBRANCH_1_acac(9.50d-14, 390.0d0, 2.62d1, -1130.0d0);
        Reaction(GC_TBRANCH_1_acac(9.50e-14, 390.0e0, 4.0e-2, 1130.0e0), [MO2, MO2], [.000CH2O, .000HO2], [1, 1], [2, 2])           GC_TBRANCH_1_acac(9.50d-14, 390.0d0, 4.0d-2, 1130.0d0);
        Reaction(1.60e-10 , [MO2, OH], [.13MOH, .87CH2O, .74HO2], [1, 1], [0, 0, 1])   #2021/09/22; Bates2021a; KHB,MSL
        Reaction(GCARR_ac(2.66e-12, 200.0e0), [MP, OH], [MO2, H2O], [1, 1], [1, 1])                        GCARR_ac(2.66d-12, 200.0d0);
        Reaction(GCARR_ac(1.14e-12, 200.0e0), [MP, OH], [CH2O, OH, H2O], [1, 1], [1, 1, 1])                  GCARR_ac(1.14d-12, 200.0d0);
        Reaction(GCARR_ac(2.66e-12, 200.0e0), [ATOOH, OH], [ATO2, H2O], [1, 1], [1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE
        Reaction(GCARR_ac(1.14e-12, 200.0e0), [ATOOH, OH], [MGLY, OH, H2O], [1, 1], [1, 1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE
        Reaction(GCARR_ac(5.50e-12, 125.0e0), [CH2O, OH], [CO, HO2, H2O], [1, 1], [1, 1, 1])                 GCARR_ac(5.50d-12, 125.0d0);
        Reaction(GC_OHHNO3_acacac(2.41e-14, 460.0e0, 2.69e-17, 2199.0e0, 6.51e-34, 1335.0e0), [HNO3, OH], [H2O, NO3], [1, 1], [1, 1])                      GC_OHHNO3_acacac(2.41d-14, 460.0d0, 2.69d-17, 2199.0d0, 6.51d-34, 1335.0d0);
        Reaction(GCARR_ac(1.80e-11, -390.0e0), [HNO2, OH], [H2O, NO2], [1, 1], [1, 1])                      GCARR_ac(1.80d-11, -390.0d0);
        Reaction(GCJPLPR_abcabc(9.05e-05, 3.4e0, -10900.0e0, 1.90e15, 0.3e0, -10900.0e0, 0.6e0), [#, M], [HO2, NO2], [1, 1], [1, 1])   #2017/02/22; JPL 15-10; BHH,MJE
        Reaction(GCARR_ac(1.30e-12, 380.0e0), [HNO4, OH], [H2O, NO2, O2], [1, 1], [1, 1, 1])                 GCARR_ac(1.30d-12, 380.0d0);
        Reaction(3.50e-12, [HO2, NO3], [OH, NO2, O2], [1, 1], [1, 1, 1])                  3.50d-12;
        Reaction(GCARR_ac(1.50e-11, 170.0e0), [NO, NO3], [.000NO2], [1, 1], [2])                        GCARR_ac(1.50d-11, 170.0d0);
        Reaction(2.20e-11, [OH, NO3], [HO2, NO2], [1, 1], [1, 1])                       2.20d-11;
        Reaction(GCJPLPR_abcabc(4.14e-04, 3.0e0, -10840.0e0, 2.76e14, -0.1e0, -10840.0e0, 0.6e0), [#, M], [NO2, NO3], [1, 1], [1, 1])   #2017/02/22; JPL 15-10; BHH,MJE
        Reaction(4.00e-13, [HCOOH, OH], [H2O, CO2, HO2], [1, 1], [1, 1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE
        Reaction(GCARR_ac(2.90e-12, -345.0e0), [MOH, OH], [HO2, CH2O], [1, 1], [1, 1])                      GCARR_ac(2.90d-12, -345.0d0);
        Reaction(GCARR_ac(4.50e-14, -1260.0e0), [NO2, NO3], [NO, NO2, O2], [1, 1], [1, 1, 1])                  GCARR_ac(4.50d-14, -1260.0d0);
        Reaction(5.80e-16, [NO3, CH2O], [HNO3, HO2, CO], [1, 1], [1, 1, 1])               5.80d-16;
        Reaction(GCARR_ac(1.40e-12, -1900.0e0), [ALD2, NO3], [HNO3, MCO3], [1, 1], [1, 1])                   GCARR_ac(1.40d-12, -1900.0d0);
        Reaction(GCJPLPR_abab(9.70e-29, 5.6e+00, 9.3e-12, 1.5e0, 0.6e0), [#, M], [PAN], [1, 1], [1])   #JPL Eval 17
        Reaction(GCJPLEQ_acabab(9.30e-29, 14000.0e0, 9.7e-29, 5.6e0, 9.3e-12, 1.5e0, 0.6e0), [PAN], [MCO3, NO2], [1], [1, 1])                           GCJPLEQ_acabab(9.30d-29, 14000.0d0, 9.7d-29, 5.6d0, 9.3d-12, 1.5d0, 0.6d0);
        Reaction(GCARR_ac(8.10e-12, 270.0e0), [MCO3, NO], [MO2, NO2, CO2], [1, 1], [1, 1, 1])                GCARR_ac(8.10d-12, 270.0d0);
        Reaction(GCARR_ac(7.66e-12, -1020.0e0), [C2H6, OH], [ETO2, H2O], [1, 1], [1, 1])   #2013/02/12; JPL 10-6; BHH,JMAO,EAM
        Reaction(GC_RO2NO_B2_aca(2.60e-12, 365.0e0, 2.0e0), [ETO2, NO], [ALD2, NO2, HO2], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF
        Reaction(GC_RO2NO_A2_aca(2.60e-12, 365.0e0, 2.0e0), [ETO2, NO], [ETNO3], [1, 1], [1])   #2019/05/10; Fisher2018; JAF
        Reaction(GCARR_ac(2.60e-12, 365.0e0), [OTHRO2, NO], [ALD2, NO2, HO2], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF
        Reaction(GC_TBRANCH_2_acabc(7.60e-12, -585.0e0, 5.87e0, 0.64e0, -816.0e0), [C3H8, OH], [B3O2], [1, 1], [1])                           GC_TBRANCH_2_acabc(7.60d-12, -585.0d0, 5.87d0, 0.64d0, -816.0d0);
        Reaction(GC_TBRANCH_2_acabc(7.60e-12, -585.0e0, 1.7e-1, -0.64e0, 816.0e0), [C3H8, OH], [A3O2], [1, 1], [1])                           GC_TBRANCH_2_acabc(7.60d-12, -585.0d0, 1.7d-1, -0.64d0, 816.0d0);
        Reaction(GC_RO2NO_B2_aca(2.90e-12, 350.0e0, 3.0e0), [A3O2, NO], [NO2, HO2, RCHO], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF
        Reaction(GC_RO2NO_A2_aca(2.90e-12, 350.0e0, 3.0e0), [A3O2, NO], [NPRNO3], [1, 1], [1])   #2019/05/10; Fisher2018; JAF
        Reaction(GCARR_ac(2.70e-12, 350.0e0), [PO2, NO], [NO2, HO2, CH2O, ALD2], [1, 1], [1, 1, 1, 1])         GCARR_ac(2.70d-12, 350.0d0);
        Reaction(GCARR_ac(9.10e-12, -405.0e0), [ALK4, OH], [R4O2], [1, 1], [1])                           GCARR_ac(9.10d-12, -405.0d0);
        Reaction(GC_RO2NO_A2_aca(2.70e-12, 350.0e0, 4.5e0), [R4O2, NO], [R4N2], [1, 1], [1])                           GC_RO2NO_A2_aca(2.70d-12, 350.0d0, 4.5d0);
        Reaction(GCARR_ac(2.80e-12, 300.0e0), [ATO2, NO], [NO2, CH2O, MCO3], [1, 1], [1, 1, 1])   #2017/07/27; Fix C creation; SAS,BHH,MJE
        Reaction(GC_RO2NO_B2_aca(2.70e-12, 360.0e0, 3.0e0), [B3O2, NO], [NO2, HO2, ACET], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF
        Reaction(GC_RO2NO_A2_aca(2.70e-12, 360.0e0, 3.0e0), [B3O2, NO], [IPRNO3], [1, 1], [1])   #2019/05/10; Fisher2018; JAF
        Reaction(GCARR_ac(2.70e-12, 350.0e0), [PRN1, NO], [.000NO2, CH2O, ALD2], [1, 1], [2, 1, 1])         GCARR_ac(2.70d-12, 350.0d0);
        Reaction(GCARR_ac(2.80e-12, -3280.0e0), [ALK4, NO3], [HNO3, R4O2], [1, 1], [1, 1])                   GCARR_ac(2.80d-12, -3280.0d0);
        Reaction(1.60e-12, [R4N2, OH], [R4N1, H2O], [1, 1], [1, 1])                     1.60d-12;
        Reaction(GCARR_ac(3.15e-14, 920.0e0), [ACTA, OH], [MO2, CO2, H2O], [1, 1], [1, 1, 1])   #2013/02/12; JPL 10-6; BHH,JMAO,EAM
        Reaction(GCARR_ac(6.00e-12, 410.0e0), [OH, RCHO], [RCO3, H2O], [1, 1], [1, 1])                     GCARR_ac(6.00d-12, 410.0d0);
        Reaction(GCJPLPR_abab(9.00e-28, 8.9e0, 7.7e-12, 0.2e0, 0.6e0), [#, M], [PPN], [1, 1], [1])   #JPL Eval 17
        Reaction(GCJPLEQ_acabab(9.00e-29, 14000.0e0, 9.00e-28, 8.9e0, 7.7e-12, 0.2e0, 0.6e0), [PPN], [RCO3, NO2], [1], [1, 1])                           GCJPLEQ_acabab(9.00d-29, 14000.0d0, 9.00d-28, 8.9d0, 7.7d-12, 0.2d0, 0.6d0);
        Reaction(6.50e-15, [RCHO, NO3], [HNO3, RCO3], [1, 1], [1, 1])                   6.50d-15;
        Reaction(1.33d-13 + 3.82d-11*exp(-2000.0e0/TEMP), [ACET, OH], [ATO2, H2O], [1, 1], [1, 1])   #JPL Eval 17, p1-62-D31; EVF
        Reaction(GCARR_ac(7.40e-13, 700.0e0), [R4O2, HO2], [R4P], [1, 1], [1])                           GCARR_ac(7.40d-13, 700.0d0);
        Reaction(GCARR_ac(7.40e-13, 700.0e0), [R4N1, HO2], [R4N2], [1, 1], [1])                          GCARR_ac(7.40d-13, 700.0d0);
        Reaction(GC_RO2HO2_aca(2.91e-13, 1300.0e0, 3.0e0), [B3O2, HO2], [RB3P], [1, 1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE
        Reaction(GC_RO2HO2_aca(2.91e-13, 1300.0e0, 3.0e0), [PRN1, HO2], [PRPN], [1, 1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE
        Reaction(GCARR_ac(1.30e-12, -25.0e0), [MEK, OH], [KO2, H2O], [1, 1], [1, 1])                       GCARR_ac(1.30d-12, -25.0d0);
        Reaction(8.00e-16, [MEK, NO3], [HNO3, KO2], [1, 1], [1, 1])                     8.00d-16;
        Reaction(3.35e-12, [EOH, OH], [HO2, ALD2], [1, 1], [1, 1])   #2013/02/12; JPL 10-6; BHH,JMAO,EAM
        Reaction(GCARR_ac(4.60e-12, 70.0e0), [ROH, OH], [HO2, RCHO], [1, 1], [1, 1])                      GCARR_ac(4.60d-12, 70.0d0);
        Reaction(4.10e-14, [ETO2, ETO2], [.000ALD2, .000HO2], [1, 1], [2, 2])         4.10d-14;
        Reaction(4.10e-14, [OTHRO2, OTHRO2], [.000ALD2, .000HO2], [1, 1], [2, 2])   #2019/05/10; Fisher2018; JAF
        Reaction(2.70e-14, [ETO2, ETO2], [EOH, ALD2], [1, 1], [1, 1])                   2.70d-14;
        Reaction(2.70e-14, [OTHRO2, OTHRO2], [EOH, ALD2], [1, 1], [1, 1])   #2019/05/10; Fisher2018; JAF
        Reaction(GCARR_ac(7.40e-13, 700.0e0), [HO2, ETO2], [ETP], [1, 1], [1])                           GCARR_ac(7.40d-13, 700.0d0);
        Reaction(GCARR_ac(7.40e-13, 700.0e0), [HO2, OTHRO2], [ETP], [1, 1], [1])   #2019/05/10; Fisher2018; JAF
        Reaction(GC_RO2HO2_aca(2.91e-13, 1300.0e0, 3.0e0), [A3O2, HO2], [RA3P], [1, 1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE
        Reaction(GC_RO2HO2_aca(2.91e-13, 1300.0e0, 3.0e0), [PO2, HO2], [PP], [1, 1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE
        Reaction(GCJPLPR_abab(4.60e-27, 4.0e0, 2.6e-11, 1.3e0, 0.5e0), [#, M], [PO2], [1, 1], [1])   #2017/02/22; JPL 15-10; BHH,MJE
        Reaction(GC_GLYCOH_B_a(8.00e-12), [GLYC, OH], [HCOOH, OH, CO], [1, 1], [1, 1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE
        Reaction(GCARR_ac(4.59e-13, -1156.0e0), [PRPE, NO3], [PRN1], [1, 1], [1])                          GCARR_ac(4.59d-13, -1156.0d0);
        Reaction(GCARR_ac(3.10e-12, 340.0e0), [GLYX, OH], [HO2, .000CO], [1, 1], [1, 2])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE
        Reaction(1.50e-11, [MGLY, OH], [MCO3, CO], [1, 1], [1, 1])                      1.50d-11;
        Reaction(GC_GLYXNO3_ac(1.40e-12, -1860.0e0), [GLYX, NO3], [HNO3, HO2, .000CO], [1, 1], [1, 1, 2])          GC_GLYXNO3_ac(1.40d-12, -1860.0d0);
        Reaction(GCARR_ac(3.36e-12, -1860.0e0), [MGLY, NO3], [HNO3, CO, MCO3], [1, 1], [1, 1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE
        Reaction(GC_HACOH_A_ac(2.15e-12, 305.0e0), [HAC, OH], [MGLY, HO2], [1, 1], [1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE
        Reaction(GCARR_ac(1.68e-12, 500.0e0), [MCO3, A3O2], [MO2, RCHO, HO2], [1, 1], [1, 1, 1])             GCARR_ac(1.68d-12, 500.0d0);
        Reaction(GCARR_ac(1.68e-12, 500.0e0), [MCO3, PO2], [MO2, ALD2, CH2O, HO2], [1, 1], [1, 1, 1, 1])       GCARR_ac(1.68d-12, 500.0d0);
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [MCO3, A3O2], [ACTA, RCHO], [1, 1], [1, 1])                  GCARR_ac(1.87d-13, 500.0d0);
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [MCO3, PO2], [ACTA, .350RCHO, .650HAC], [1, 1], [1, 0, 0])   GCARR_ac(1.87d-13, 500.0d0);
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [RCO3, MO2], [RCOOH, CH2O], [1, 1], [1, 1])                  GCARR_ac(1.87d-13, 500.0d0);
        Reaction(GCARR_ac(8.78e-12, 200.0e0), [R4P, OH], [.791OH, .209R4O2, .791RCHO], [1, 1], [0, 0, 0])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE
        Reaction(GCARR_ac(6.13e-13, 200.0e0), [RP, OH], [RCO3], [1, 1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE
        Reaction(GCARR_ac(8.78e-12, 200.0e0), [PP, OH], [.791OH, .209PO2, .791HAC], [1, 1], [0, 0, 0])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE
        Reaction(GCARR_ac(4.82e-11, -400.0e0), [LVOC, OH], [OH], [1, 1], [1])   #2017/06/14; Marais2016; EAM
        Reaction(GCARR_ac(6.13e-13, 200.0e0), [OH, MAP], [MCO3], [1, 1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE
        Reaction(1.40e-18, [C2H6, NO3], [ETO2, HNO3], [1, 1], [1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE
        Reaction(GCARR_ac(2.50e-12, 500.0e0), [MCO3, MCO3], [.000MO2], [1, 1], [2])                     GCARR_ac(2.50d-12, 500.0d0);
        Reaction(GCARR_ac(1.80e-12, 500.0e0), [MCO3, MO2], [CH2O, MO2, HO2], [1, 1], [1, 1, 1])              GCARR_ac(1.80d-12, 500.0d0);
        Reaction(GCARR_ac(2.00e-13, 500.0e0), [MCO3, MO2], [ACTA, CH2O], [1, 1], [1, 1])                   GCARR_ac(2.00d-13, 500.0d0);
        Reaction(GCARR_ac(1.68e-12, 500.0e0), [ATO2, MCO3], [MO2, MCO3, CH2O], [1, 1], [1, 1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE
        Reaction(GCARR_ac(1.68e-12, 500.0e0), [KO2, MCO3], [MO2, ALD2, MCO3], [1, 1], [1, 1, 1])             GCARR_ac(1.68d-12, 500.0d0);
        Reaction(GCARR_ac(1.68e-12, 500.0e0), [B3O2, MCO3], [MO2, HO2, ACET], [1, 1], [1, 1, 1])             GCARR_ac(1.68d-12, 500.0d0);
        Reaction(GCARR_ac(1.68e-12, 500.0e0), [PRN1, MCO3], [MO2, NO2, CH2O, ALD2], [1, 1], [1, 1, 1, 1])      GCARR_ac(1.68d-12, 500.0d0);
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [R4O2, MCO3], [MEK, ACTA], [1, 1], [1, 1])                   GCARR_ac(1.87d-13, 500.0d0);
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [ATO2, MCO3], [MGLY, ACTA], [1, 1], [1, 1])   #2017/07/27; Fix C creation; SAS,BHH,MJE
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [KO2, MCO3], [MEK, ACTA], [1, 1], [1, 1])                    GCARR_ac(1.87d-13, 500.0d0);
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [R4N1, MCO3], [RCHO, ACTA, NO2], [1, 1], [1, 1, 1])            GCARR_ac(1.87d-13, 500.0d0);
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [PRN1, MCO3], [RCHO, ACTA, NO2], [1, 1], [1, 1, 1])            GCARR_ac(1.87d-13, 500.0d0);
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [B3O2, MCO3], [ACET, ACTA], [1, 1], [1, 1])                  GCARR_ac(1.87d-13, 500.0d0);
        Reaction(GCARR_ac(1.68e-12, 500.0e0), [MCO3, ETO2], [MO2, ALD2, HO2], [1, 1], [1, 1, 1])             GCARR_ac(1.68d-12, 500.0d0);
        Reaction(GCARR_ac(1.68e-12, 500.0e0), [MCO3, OTHRO2], [MO2, ALD2, HO2], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [MCO3, ETO2], [ACTA, ALD2], [1, 1], [1, 1])                  GCARR_ac(1.87d-13, 500.0d0);
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [MCO3, OTHRO2], [ACTA, ALD2], [1, 1], [1, 1])   #2019/05/10; Fisher2018; JAF
        Reaction(GCARR_ac(8.50e-13, -2450.0e0), [NO3, NO3], [.000NO2, O2], [1, 1], [2, 1])                  GCARR_ac(8.50d-13, -2450.0d0);
        Reaction(GCJPLPR_abcabc(1.05e-02, 4.8e+00, -11234.0e0, 7.58e16, 2.1e0, -11234.0e0, 0.6e0), [#, M], [MO2, NO2], [1, 1], [1, 1])   #2012/02/12; Browne2011; ECB
        Reaction(GCARR_ac(1.20e-11, -280.0e0), [DMS, OH], [SO2, MO2, CH2O], [1, 1], [1, 1, 1])                GCARR_ac(1.20d-11, -280.0d0);
        Reaction(GC_DMSOH_acac(8.20e-39, 5376.0e0, 1.05e-5, 3644.0e0), [DMS, OH], [.750SO2, .250MSA, MO2], [1, 1], [0, 0, 1])       GC_DMSOH_acac(8.20d-39, 5376.0d0, 1.05d-5, 3644.0d0);
        Reaction(GCARR_ac(1.90e-13, 530.0e0), [DMS, NO3], [SO2, HNO3, MO2, CH2O], [1, 1], [1, 1, 1, 1])        GCARR_ac(1.90d-13, 530.0d0);
        Reaction(GCJPLPR_aba(3.30e-31, 4.3e+00, 1.6e-12, 0.6e0), [#, M], [SO4, HO2], [1, 1], [1, 1])                  GCJPLPR_aba(3.30d-31, 4.3d+00, 1.6d-12, 0.6d0);
        Reaction(GCARR_ac(1.60e-11, -780.0e0), [Br, O3], [BrO, O2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP
        Reaction(GCARR_ac(4.50e-12, 460.0e0), [BrO, HO2], [HOBr, O2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP
        Reaction(GCARR_ac(4.80e-12, -310.0e0), [Br, HO2], [HBr, O2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP
        Reaction(GCARR_ac(5.50e-12, 200.0e0), [HBr, OH], [Br, H2O], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP
        Reaction(GCARR_ac(2.40e-12,  40.0e0), [BrO, BrO], [.000Br, O2], [1, 1], [2, 1])   #2012/06/07; Parrella2012; JPP
        Reaction(GCARR_ac(2.80e-14, 860.0e0), [BrO, BrO], [Br2, O2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP
        Reaction(GCARR_ac(8.80e-12, 260.0e0), [BrO, NO], [Br, NO2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP
        Reaction(4.90e-11, [Br, BrNO3], [Br2, NO3], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP
        Reaction(GCARR_ac(2.10e-11, 240.0e0), [Br2, OH], [HOBr, Br], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP
        Reaction(GCARR_ac(1.20e-10, -430.0e0), [HOBr, O], [OH, BrO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(5.80e-12, -1500.0e0), [HBr, O], [OH, Br], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(1.70e-11, 250.0e0), [BrO, OH], [Br, HO2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP
        Reaction(1.60e-11, [Br, NO3], [BrO, NO2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP
        Reaction(GCARR_ac(1.70e-11, -800.0e0), [Br, CH2O], [HBr, HO2, CO], [1, 1], [1, 1, 1])   #2012/06/07; Parrella2012; JPP
        Reaction(GCARR_ac(1.80e-11, -460.0e0), [Br, ALD2], [HBr, MCO3], [1, 1], [1, 1])   #2017/07/27; Parrella2012,Fix C creation; SAS,BHH,MJE
        Reaction(GCARR_ac(1.66e-10, -7000.0e0), [Br, ACET], [HBr, ATO2], [1, 1], [1, 1])   #2017/07/27; Parrella2012,Fix C creation; SAS,BHH,MJE
        Reaction(GCARR_ac(2.36e-10, -6411.0e0), [Br, C2H6], [HBr, ETO2], [1, 1], [1, 1])   #2017/07/27; Parrella2012,Fix C creation; SAS,BHH,MJE
        Reaction(GCARR_ac(8.77e-11, -4330.0e0), [Br, C3H8], [HBr, A3O2], [1, 1], [1, 1])   #2017/07/27; Parrella2012,Fix C creation; SAS,BHH,MJE
        Reaction(GCARR_ac(9.00e-13, -360.0e0), [CHBr3, OH], [.000Br], [1, 1], [3])   #2017/02/22; JPL 15-10; BHH,MJE
        Reaction(GCARR_ac(2.00e-12, -840.0e0), [CH2Br2, OH], [.000Br], [1, 1], [2])   #2012/06/07; Parrella2012; JPP
        Reaction(GCARR_ac(1.42e-12, -1150.0e0), [CH3Br, OH], [Br, H2O, HO2], [1, 1], [1, 1, 1])   #2017/03/08; JPL 15-10; TS,BHH,MJE
        Reaction(GCARR_ac(1.63e-10, 60.0e0), [O1D, H2O], [.000OH], [1, 1], [2])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(2.15e-11, 110.0e0), [O1D, N2], [O, N2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(3.30e-11, 55.0e0), [O1D, O2], [O, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(1.20e-10, [O1D, H2], [H, OH], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(4.63e-11, 20.0e0), [O1D, N2O], [N2, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(7.25e-11, 20.0e0), [O1D, N2O], [.000NO], [1, 1], [2])   #2014/02/03; Eastham2014; SDE
        Reaction(1.31e-10, [O1D, CH4], [MO2, OH], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(0.09e-10, [O1D, CH4], [CH2O, H2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(0.35e-10, [O1D, CH4], [CH2O, H, HO2], [1, 1], [1, 1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(8.00e-12, -2060.0e0), [O, O3], [.000O2], [1, 1], [2])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(2.80e-12, -1800.0e0), [OH, H2], [H2O, H], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(1.80e-11, 180.0e0), [O, OH], [O2, H], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(3.00e-11, 200.0e0), [HO2, O], [OH, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(1.20e-10, [O1D, O3], [.000O2], [1, 1], [2])   #2014/02/03; Eastham2014; SDE
        Reaction(1.20e-10, [O1D, O3], [.000O, O2], [1, 1], [2, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(2.10e-11, -2200.0e0), [OCS, O], [CO, SO2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(1.10e-13, -1200.0e0), [OCS, OH], [CO2, SO2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(5.10e-12, 210.0e0), [NO2, O], [NO, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(1.00e-11, [NO3, O], [NO2, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(1.40e-12, -2000.0e0), [H2O2, O], [OH, HO2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(1.40e-10, -470.0e0), [H, O3], [OH, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(7.20e-11, [H, HO2], [.000OH], [1, 1], [2])   #2014/02/03; Eastham2014; SDE
        Reaction(1.60e-12, [H, HO2], [O, H2O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(6.90e-12, [H, HO2], [H2, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(1.50e-11, -3600.0e0), [N, O2], [NO, O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(2.10e-11, 100.0e0), [N, NO], [N2, O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(5.80e-12, 220.0e0), [N, NO2], [N2O, O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(1.90e-11, 230.0e0), [BrO, O], [Br, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(3.40e-11, -1600.0e0), [CH2O, O], [CO, HO2, OH], [1, 1], [1, 1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(1.80e-10, [O1D, CH3Br], [.440BrO, MO2, .560Br], [1, 1], [0, 1, 0])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(2.60e-12, -1100.0e0), [OH, Cl2], [HOCl, Cl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(1.80e-11, -600.0e0), [MO2, ClO], [ClOO, HO2, CH2O], [1, 1], [1, 1, 1])   #2017/03/20; JPL 15-10; TS,BHH,MJE
        Reaction(GCARR_ac(7.40e-12, 270.0e0), [OH, ClO], [HO2, Cl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(6.00e-13, 230.0e0), [OH, ClO], [HCl, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(1.40e-12, 600.0e0), [OH, OClO], [HOCl, O2], [1, 1], [1, 1])   #2017/02/22; JPL 15-10; BHH,MJE
        Reaction(GCARR_ac(6.00e-13, 670.0e0), [OH, Cl2O2], [HOCl, ClOO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(1.80e-12, -250.0e0), [OH, HCl], [H2O, Cl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(3.00e-12, -500.0e0), [OH, HOCl], [H2O, ClO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(2.40e-12, -1250.0e0), [OH, ClNO2], [HOCl, NO2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(1.20e-12, -330.0e0), [OH, ClNO3], [HOCl, NO3], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(1.96e-12, -1200.0e0), [OH, CH3Cl], [Cl, HO2, H2O], [1, 1], [1, 1, 1])   #2017/02/22; JPL 15-10; BHH,MJE
        Reaction(GCARR_ac(2.61e-12, -944.0e0), [OH, CH2Cl2], [.000Cl, HO2], [1, 1], [2, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(4.69e-12, -1134.0e0), [OH, CHCl3], [.000Cl, HO2], [1, 1], [3, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(1.64e-12, -1520.0e0), [OH, CH3CCl3], [.000Cl, H2O], [1, 1], [3, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(9.20e-13, -1560.0e0), [OH, HCFC22], [Cl, H2O], [1, 1], [1, 1])   #2017/02/22; JPL 15-10; BHH,MJE
        Reaction(GCARR_ac(1.25e-12, -1600.0e0), [OH, HCFC141b], [.000Cl, H2O], [1, 1], [2, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(1.30e-12, -1770.0e0), [OH, HCFC142b], [Cl, H2O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(7.40e-13, -900.0e0), [OH, HCFC123], [.000Cl, H2O], [1, 1], [2, 1])   #2017/02/22; JPL 15-10; BHH,MJE
        Reaction(GCARR_ac(7.10e-12, -1270.0e0), [CH4, Cl], [HCl, MO2], [1, 1], [1, 1])   #2017/03/08; JPL 15-10; TS,BHH,MJE
        Reaction(GCARR_ac(7.32e-11, -30.0e0), [CH2O, Cl], [CO, HCl, HO2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b; TS,JAS,SDE
        Reaction(GCARR_ac(2.30e-11, -200.0e0), [Cl, O3], [ClO, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(3.05e-11, -2270.0e0), [Cl, H2], [H, HCl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(1.10e-11, -980.0e0), [Cl, H2O2], [HO2, HCl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(1.40e-11, 270.0e0), [Cl, HO2], [O2, HCl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(3.60e-11, -375.0e0), [Cl, HO2], [OH, ClO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(2.80e-11, 85.0e0), [ClO, O], [Cl, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(2.60e-12, 290.0e0), [ClO, HO2], [O2, HOCl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(6.40e-12, 290.0e0), [ClO, NO], [Cl, NO2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(1.00e-12, -1590.0e0), [ClO, ClO], [Cl2, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(3.00e-11, -2450.0e0), [ClO, ClO], [Cl, ClOO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(3.50e-13, -1370.0e0), [ClO, ClO], [OClO, Cl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(2.30e-10, [ClOO, Cl], [Cl2, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(1.20e-11, [ClOO, Cl], [.000ClO], [1, 1], [2])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(9.50e-13, 550.0e0), [ClO, BrO], [Br, OClO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(2.30e-12, 260.0e0), [ClO, BrO], [Br, ClOO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(4.10e-13, 290.0e0), [ClO, BrO], [BrCl, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(3.60e-12, -840.0e0), [ClNO3, O], [ClO, NO3], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(6.50e-12, 135.0e0), [ClNO3, Cl], [Cl2, NO3], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(2.17e-11, -1130.0e0), [CH3Cl, Cl], [CO, .000HCl, HO2], [1, 1], [1, 2, 1])   #2014/02/03; Eastham2014; SDE
        Reaction(GCARR_ac(1.24e-12, -1070.0e0), [CH2Cl2, Cl], [CO, HCl, .000Cl, HO2], [1, 1], [1, 1, 2, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(3.77e-12, -1011.0e0), [CHCl3, Cl], [CO, HCl, .000Cl, HO2], [1, 1], [1, 1, 3, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(2.00e-13, [Cl, HCOOH], [HCl, CO2, H2O], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(1.60e-10, [Cl, MO2], [ClO, CH2O, HO2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(5.7e-11, [Cl, MP], [HCl, MO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(7.2e-11, -70.0e0), [Cl, C2H6], [HCl, ETO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(7.4e-11, [Cl, ETO2], [ClO, HO2, ALD2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(7.4e-11, [Cl, OTHRO2], [ClO, HO2, ALD2], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF
        Reaction(5.5e-11, [Cl, MOH], [HCl, CH2O, HO2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(9.6e-11, [Cl, EOH], [HCl, ALD2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(2.8e-14, [Cl, ACTA], [HCl, MO2, CO2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(6.54e-11, 60.0e0), [Cl, C3H8], [HCl, B3O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(8.12e-11, -90.0e0), [Cl, C3H8], [HCl, A3O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(7.70e-11, -1000.0e0), [Cl, ACET], [HCl, ATO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(7.60e-11, 500.0e0), [Cl, ISOP], [HCl, .5IHOO1, .5IHOO4], [1, 1], [1, 0, 0])   #2019/11/06; Sherwen2016b;KHB,TS,JAS,SDE
        Reaction(2.05e-10, [Cl, ALK4], [HCl, R4O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(3.60e-12, [Br, PRPE], [HBr, PO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(8.40e-11, -2620.0e0), [INO, INO], [I2, .000NO], [1, 1], [1, 2])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(2.90e-11, -2600.0e0), [IONO, IONO], [I2, .000NO2], [1, 1], [1, 2])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(1.50e-12, [I2, NO3], [I, IONO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(9.10e-11, -146.0e0), [IONO2, I], [I2, NO3], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(1.20e-11, [I, BrO], [IO, Br], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(3.00e-12, 510.0e0), [IO, BrO], [Br, I, O2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(1.20e-11, 510.0e0), [IO, BrO], [Br, OIO], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(1.50e-10, [OIO, OIO], [I2O4], [1, 1], [1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(1.10e-12, 542.0e0), [OIO, NO], [IO, NO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(5.10e-12, 280.0e0), [IO, ClO], [I, OClO], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(2.81e-12, 280.0e0), [IO, ClO], [I, Cl, O2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(1.02e-12, 280.0e0), [IO, ClO], [ICl, O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(2.30e-11, -870.0e0), [I, O3], [IO, O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2017;TS,JAS,SDE
        Reaction(GCARR_ac(1.50e-11, -1090.0e0), [I, HO2], [HI, O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(1.80e-10, [I2, OH], [HOI, I], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(3.00e-11, [HI, OH], [I, H2O], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(5.00e-12, [HOI, OH], [IO, H2O], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(1.30e-11, 570.0e0), [IO, HO2], [HOI, O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(9.10e-12, 240.0e0), [IO, NO], [I, NO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(6.00e-12, 500.0e0), [IO, IO], [I, OIO], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(GCARR_ac(2.90e-12, -1100.0e0), [CH3I, OH], [H2O, I, MO2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE
        Reaction(2.40e-12, [ETHLN, OH], [CH2O, CO2, NO2], [1, 1], [1, 1, 1])   #2017/06/15, Marais2016, EAM
        Reaction(6.70e-13, [PROPNN, OH], [NO2, MGLY], [1, 1], [1, 1])   #2017/07/14; MCMv3.3; KRT,JAF,CCM,EAM,KHB,RHS
        Reaction(1.20e-15, [CH2OO, CO], [CH2O], [1, 1], [1])   #2015/09/25; Millet2015; DBM,EAM
        Reaction(1.00e-14, [CH2OO, NO], [CH2O, NO2], [1, 1], [1, 1])   #2015/09/25; Millet2015; DBM,EAM
        Reaction(1.00e-15, [CH2OO, NO2], [CH2O, NO3], [1, 1], [1, 1])   #2015/09/25; Millet2015; DBM,EAM
        Reaction(1.40e-12, [CH2OO, O3], [CH2O], [1, 1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(3.70e-11, [CH2OO, SO2], [CH2O, SO4], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(1.20e-15, [CH3CHOO, CO], [ALD2], [1, 1], [1])   #2015/09/25; Millet2015; DBM,EAM
        Reaction(1.00e-14, [CH3CHOO, NO], [ALD2, NO2], [1, 1], [1, 1])   #2015/09/25; Millet2015; DBM,EAM
        Reaction(1.00e-15, [CH3CHOO, NO2], [ALD2, NO3], [1, 1], [1, 1])   #2015/09/25; Millet2015; DBM,EAM
        Reaction(7.00e-14, [CH3CHOO, SO2], [ALD2, SO4], [1, 1], [1, 1])   #2015/09/25; Millet2015; DBM,EAM
        Reaction(6.00e-18, [CH3CHOO, H2O], [ALD2, H2O2], [1, 1], [1, 1])   #2015/09/25; Millet2015; DBM,EAM
        Reaction(1.00e-17, [CH3CHOO, H2O], [ACTA], [1, 1], [1])   #2015/09/25; Millet2015; DBM,EAM
        Reaction(GCARR_ac(1.21e-11, 440.0e0), [MTPA, OH], [PIO2], [1, 1], [1])   #2017/03/23; IUPAC2010; EVF
        Reaction(GCARR_ac(1.21e-11, 440.0e0), [MTPO, OH], [PIO2], [1, 1], [1])   #2017/03/23; IUPAC2010; EVF
        Reaction(1.50e-11, [PIO2, HO2], [PIP], [1, 1], [1])   #2017/03/23; Roberts1992; EVF
        Reaction(1.20e-12, [PIO2, NO3], [HO2, NO2, RCHO, MEK], [1, 1], [1, 1, 1, 1])   #2017/03/23; Roberts1992; EVF
        Reaction(GCARR_ac(8.33e-13, 490.0e0), [MTPA, NO3], [.100OLNN, .900OLND], [1, 1], [0, 0])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS
        Reaction(GCARR_ac(8.33e-13, 490.0e0), [MTPO, NO3], [.100OLNN, .900OLND], [1, 1], [0, 0])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS
        Reaction(GCARR_ac(4.20e-11, 401.0e0), [LIMO, OH], [LIMO2], [1, 1], [1])   #2017/07/14; Gill2002; KRT,JAF,CCM,EAM,KHB,RHS
        Reaction(1.22e-11, [LIMO, NO3], [.500OLNN, .500OLND], [1, 1], [0, 0])   #2017/07/14; Fry2014,Atkinson2003; KRT,JAF,CCM,EAM,KHB,RHS
        Reaction(1.50e-11, [LIMO2, HO2], [PIP], [1, 1], [1])   #2017/07/14; Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS
        Reaction(4.00e-12, [OLNN, NO], [HO2, NO2, MONITS], [1, 1], [1, 1, 1])   #2017/07/14; Browne2014,Goliff2013; KRT,JAF,CCM,EAM,KHB,RHS
        Reaction(GCARR_ac(1.66e-13, 1300.0e0), [OLNN, HO2], [.700MONITS, .300MONITU], [1, 1], [0, 0])   #2017/07/14; Browne2014,Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS
        Reaction(GCARR_ac(1.66e-13, 1300.0e0), [OLND, HO2], [.700MONITS, .300MONITU], [1, 1], [0, 0])   #2017/07/14; Browne2014,Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS
        Reaction(4.80e-12, [MONITS, OH], [HONIT], [1, 1], [1])   #2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS
        Reaction(7.29e-11, [MONITU, OH], [HONIT], [1, 1], [1])   #2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS
        Reaction(1.67e-16, [MONITU, O3], [HONIT], [1, 1], [1])   #2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS
        Reaction(GCARR_ac(3.15e-13, -448.0e0), [MONITU, NO3], [HONIT], [1, 1], [1])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS
        Reaction(GCARR_ac(3.15e-13, -448.0e0), [MONITS, NO3], [HONIT], [1, 1], [1])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS
        Reaction(2.78e-04, [IONITA], [INDIOL, HNO3], [1], [1, 1])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS
        Reaction(2.78e-04, [MONITA], [INDIOL, HNO3], [1], [1, 1])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS
        Reaction(GC_OHHNO3_acacac(2.41e-14, 460.0e0, 2.69e-17, 2199.0e0, 6.51e-34, 1335.0e0), [HONIT, OH], [NO3, HAC], [1, 1], [1, 1])   #2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS
        Reaction(GCARR_ac(8.00e-13, -1000.0e0), [MENO3, OH], [CH2O, NO2], [1, 1], [1, 1])   #2019/05/16; JPL 15-10,Fisher2018; JAF
        Reaction(GCARR_ac(1.00e-12, -490.0e0), [ETNO3, OH], [ALD2, NO2], [1, 1], [1, 1])   #2019/05/16; JPL 15-10,Fisher2018; JAF
        Reaction(GCARR_ac(1.20e-12, -320.0e0), [IPRNO3, OH], [ACET, NO2], [1, 1], [1, 1])   #2019/05/16; JPL 15-10,Fisher2018; JAF
        Reaction(7.10e-13, [NPRNO3, OH], [RCHO, NO2], [1, 1], [1, 1])   #2019/05/16; JPL 15-10,Fisher2018; JAF
        Reaction(GC_ISO1(1.7e-11, 3.90e2, 9.33e-2, 5.05e15, -1.22e4, 1.79e14, -8.830e3), [ISOP, OH], [LISOPOH, IHOO1], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_ISO1(1.0e-11, 3.90e2, 2.26e-1, 2.22e9, -7.160e3, 1.75e14, -9.054e3), [ISOP, OH], [LISOPOH, IHOO4], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(ARRPLUS_abde(2.12e-13, -1300e0, -0.1644e0, 7.0485e-4), [IHOO1, HO2], [RIPC], [1, 1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(ARRPLUS_abde(2.12e-13, -1300e0, -0.2038e0, 9.0435e-4), [IHOO4, HO2], [RIPD], [1, 1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(ARRPLUS_abde(1.04e11, 9.746e3,  1.1644e0, -7.0485e-4), [IHOO1], [CH2O, OH, MVK], [1], [1, 1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(ARRPLUS_abde(1.88e11, 9.752e3, 1.2038e0, -9.0435e-4), [IHOO4], [MACR, OH, CH2O], [1], [1, 1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(ARRPLUS_ade(6.92e-14, 1.1644e0, -7.0485e-4), [IHOO1, IHOO1], [MVK, HO2, CH2O], [1, 1], [2, 2, 2])   #2019/11/06; Bates2019; KHB
        Reaction(ARRPLUS_ade(5.74e-12, 1.2038e0, -9.0435e-4), [IHOO4, IHOO4], [MACR, HO2, CH2O], [1, 1], [2, 2, 2])   #2019/11/06; Bates2019; KHB
        Reaction(ARRPLUS_ade(1.54e-12, 2.3682e0, -1.6092e-3), [IHOO1, IHOO4], [MACR, MVK, HO2, CH2O], [1, 1], [1, 1, 2, 2])   #2019/11/06; Bates2019; KHB
        Reaction(ARRPLUS_ade(2.0e-12, 1.1644e0, -7.0485e-4), [IHOO1, MO2], [MVK, HO2, CH2O], [1, 1], [1, 2, 2])   #2019/11/06; Bates2019; KHB
        Reaction(ARRPLUS_ade(2.0e-12, 1.2038e0, -9.0435e-4), [IHOO4, MO2], [MACR, HO2, CH2O], [1, 1], [1, 2, 2])   #2019/11/06; Bates2019; KHB
        Reaction(GC_NIT(2.7e-12, 3.50e2, 1.19e0,  6.0e0, 1.1644e0, 7.05e-4), [IHOO1, NO], [IHN2], [1, 1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_ALK(2.7e-12, 3.50e2, 1.19e0,  6.0e0, 1.1644e0, 7.05e-4), [IHOO1, NO], [NO2, MVK, HO2, CH2O], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_NIT(2.7e-12, 3.50e2, 1.421e0, 6.0e0, -0.1644e0, -7.05e-4), [IHOO1, NO], [IHN4], [1, 1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_NIT(2.7e-12, 3.50e2, 1.297e0, 6.0e0, 1.2038e0, 9.04e-4), [IHOO4, NO], [IHN3], [1, 1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_ALK(2.7e-12, 3.50e2, 1.297e0, 6.0e0, 1.2038e0, 9.04e-4), [IHOO4, NO], [NO2, MACR, HO2, CH2O], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_NIT(2.7e-12, 3.50e2, 1.421e0, 6.0e0, -0.2038e0, -9.04e-4), [IHOO4, NO], [IHN1], [1, 1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(3.00e-12, 650.0e0), [IDC, OH], [CO, HO2, MVKPC], [1, 1], [1, 1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(1.59e+13, -10000.0e0), [IHPOO1], [.176ICPDH, .824IDHPE, OH], [1], [0, 0, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_NIT(2.7e-12, 3.50e2, 2.1e0, 9.0e0, 1.0e0, 0.0e0), [IHPOO1, NO], [ITHN], [1, 1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(2.91e+13, -10000.0e0), [IHPOO2], [.548ICPDH, .452IDHPE, OH], [1], [0, 0, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_NIT(2.7e-12, 3.50e2, 2.315e0, 9.0e0, 1.0e0, 0.0e0), [IHPOO2, NO], [ITHN], [1, 1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(1.875e+13, -10000.0e0), [IHPOO3], [IDHPE], [1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_ALK(2.7e-12, 3.50e2, 3.079e0, 9.0e0, 1.0e0, 0.0e0), [IHPOO3, NO], [GLYC, HAC, NO2, OH], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_NIT(2.7e-12, 3.50e2, 3.079e0, 9.0e0, 1.0e0, 0.0e0), [IHPOO3, NO], [ITHN], [1, 1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(1.05e-11, -400.0e0), [IEPOXA, OH], [ICHE, HO2], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_EPO_a(5.82e-11, -4.00e2, 1.14e-20), [IEPOXA, OH], [.67IEPOXAOO, .33IEPOXBOO], [1, 1], [0, 0])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(8.25e-12, -400.0e0), [IEPOXB, OH], [ICHE, HO2], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_EPO_a(3.75e-11, -4.00e2, 8.91e-21), [IEPOXB, OH], [.81IEPOXAOO, .19IEPOXBOO], [1, 1], [0, 0])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(1.875e+13, -10000.0e0), [IEPOXAOO], [IDCHP, HO2], [1], [1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(1.0e+7, -5000.0e0), [IEPOXAOO], [OH, CO, MVKDH], [1], [1, 1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_NIT(2.7e-12, 3.50e2, 13.098e0, 8.0e0, 1.0e0, 0.0e0), [IEPOXAOO, NO], [ITCN], [1, 1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(1.875e+13, -10000.0e0), [IEPOXBOO], [IDCHP, HO2], [1], [1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(1.0e+7, -5000.0e0), [IEPOXBOO], [CO, OH, MCRDH], [1], [1, 1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_NIT(2.7e-12, 3.50e2, 16.463e0, 8.0e0, 1.0e0, 0.0e0), [IEPOXBOO, NO], [ITCN], [1, 1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_NIT(2.7e-12, 3.50e2, 13.098e0, 8.0e0, 1.0e0, 0.0e0), [ICHOO, NO], [ITCN], [1, 1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(1.875e+13, -10000.0e0), [ICHOO], [HO2, .000CO, HAC, OH], [1], [1, 2, 1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(2.70e-12, 350.0e0), [HPALD1OO, NO], [NO2, OH, CO2, MVK], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(2.38e-13, 1300.0e0), [HPALD1OO, HO2], [OH, OH, CO2, MVK], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(2.70e-12, 350.0e0), [HPALD2OO, NO], [NO2, OH, CO2, MACR], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(2.38e-13, 1300.0e0), [HPALD2OO, HO2], [OH, OH, CO2, MACR], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(7.14e-12, 390.0e0), [IHN2, OH], [ISOPNOO1], [1, 1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_EPO_a(6.30e-12, 390.0e0, 1.62e-19), [IHN2, OH], [.67IEPOXA, .33IEPOXB, NO2], [1, 1], [0, 0, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(1.02e-11, 390.0e0), [IHN3, OH], [ISOPNOO2], [1, 1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_EPO_a(1.05e-11, 390.0e0, 2.49e-19), [IHN3, OH], [.67IEPOXA, .33IEPOXB, NO2], [1, 1], [0, 0, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_EPO_a(1.55e-11, 390.0e0, 2.715e-19), [IHN1, OH], [IEPOXD, NO2], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(2.04e-11, 390.0e0), [IHN1, OH], [IDHNDOO1], [1, 1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_EPO_a(9.52e-12, 390.0e0, 2.715e-19), [IHN4, OH], [IEPOXD, NO2], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(2.95e-11, 390.0e0), [IHN4, OH], [IDHNDOO2], [1, 1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(1.875e+13, -10000.0e0), [ISOPNOO1], [ITCN, HO2], [1], [1, 1])   #2019/11/06; Bates2019; KHB
        Reaction(GC_NIT(2.7e-12, 350.0e0, 6.32e0, 11.0e0, 1.0e0, 0.0e0), [ISOPNOO1, NO], [IDN], [1, 1], [1])   #2019/11/06; Bates2019; KHB
        Reaction(GCARR_ac(1.875e+13, -10000.0e0), [ISOPNOO2], [ITCN, HO2], [1], [1, 1])   #2019/11/06; Bates2019; KHB


    
    
end
