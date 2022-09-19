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
        
function fullchem()

    @variables A3O2(T)      = 0.699   [unit = u"ppb"]   #CH3CH2CH2OO; Primary RO2 from C3H8
    @variables ACET(T)      = 0.272   [unit = u"ppb"]   #CH3C(O)CH3; Acetone
    @variables ACTA(T)      = 0.531   [unit = u"ppb"]   #CH3C(O)OH; Acetic acid
    @variables AERI(T)      = 0.374   [unit = u"ppb"]   #I; Dissolved iodine
    @variables ALD2(T)      = 0.402   [unit = u"ppb"]   #CH3CHO; Acetaldehyde
    @variables ALK4(T)      = 0.446   [unit = u"ppb"]   #>= C4 alkanes
    @variables AONITA(T)    = 0.369   [unit = u"ppb"]   #Aerosol-phase organic nitrate from aromatic precursors
    @variables AROMRO2(T)   = 0.146   [unit = u"ppb"]   #generic peroxy radical from aromatic oxidation
    @variables AROMP4(T)    = 0.846   [unit = u"ppb"]   #Generic C4 product from aromatic oxidation
    @variables AROMP5(T)    = 0.308   [unit = u"ppb"]   #Generic C5 product from aromatic oxidation
    @variables ATO2(T)      = 0.436   [unit = u"ppb"]   #CH3C(O)CH2O2; RO2 from acetone
    @variables ATOOH(T)     = 0.522   [unit = u"ppb"]   #CH3C(O)CH2OOH; ATO2 peroxide
    @variables B3O2(T)      = 0.994   [unit = u"ppb"]   #CH3CH(OO)CH3; Secondary RO2 from C3H8
    @variables BALD(T)      = 0.757   [unit = u"ppb"]   #benzaldehyde and tolualdehyde
    @variables BENZ(T)      = 0.955   [unit = u"ppb"]   #C6H6; Benzene
    @variables BENZO(T)     = 0.75   [unit = u"ppb"]   #C6H5O radical
    @variables BENZO2(T)    = 0.667   [unit = u"ppb"]   #C6H5O2 radical
    @variables BENZP(T)     = 0.134   [unit = u"ppb"]   #hydroperoxide from BENZO2
    @variables Br(T)        = 0.87   [unit = u"ppb"]   #Br; Atomic bromine
    @variables Br2(T)       = 0.143   [unit = u"ppb"]   #Br2; Molecular bromine
    @variables BrCl(T)      = 0.426   [unit = u"ppb"]   #BrCl; Bromine chloride
    @variables BrNO2(T)     = 0.153   [unit = u"ppb"]   #BrNO2; Nitryl bromide
    @variables BrNO3(T)     = 0.217   [unit = u"ppb"]   #BrNO3; Bromine nitrate
    @variables BrO(T)       = 0.65   [unit = u"ppb"]   #BrO; Bromine monoxide
    @variables BRO2(T)      = 0.739   [unit = u"ppb"]   #C6H5O2 ; Peroxy radical from BENZ oxidation
    @variables BrSALA(T)    = 0.58   [unit = u"ppb"]   #Br; Fine sea salt bromine
    @variables BrSALC(T)    = 0.716   [unit = u"ppb"]   #Br; Course sea salt bromine
    @variables BZCO3(T)     = 0.809   [unit = u"ppb"]   #benzoylperoxy radical
    @variables BZCO3H(T)    = 0.926   [unit = u"ppb"]   #perbenzoic acid
    @variables BZPAN(T)     = 0.6   [unit = u"ppb"]   #peroxybenzoyl nitrate
    @variables C2H2(T)      = 0.481   [unit = u"ppb"]   #C2H2; Ethyne
    @variables C2H4(T)      = 0.874   [unit = u"ppb"]   #Ethylene
    @variables C2H6(T)      = 0.919   [unit = u"ppb"]   #C2H6; Ethane
    @variables C3H8(T)      = 0.804   [unit = u"ppb"]   #C3H8; Propane
    @variables C4HVP1(T)    = 0.141   [unit = u"ppb"]   #C4 hydroxy-vinyl-peroxy radicals from HPALDs
    @variables C4HVP2(T)    = 0.65   [unit = u"ppb"]   #C4 hydroxy-vinyl-peroxy radicals from HPALDs
    @variables CCl4(T)      = 0.339   [unit = u"ppb"]   #CCl4; Carbon tetrachloride
    @variables CFC11(T)     = 0.752   [unit = u"ppb"]   #CCl3F ; CFC-11, R-11, Freon 11
    @variables CFC12(T)     = 0.421   [unit = u"ppb"]   #CCl2F2; CFC-12, R-12, Freon 12
    @variables CFC113(T)    = 0.575   [unit = u"ppb"]   #C2Cl3F3; CFC-113, Freon 113
    @variables CFC114(T)    = 0.4   [unit = u"ppb"]   #C2Cl2F4; CFC-114, Freon 114
    @variables CFC115(T)    = 0.743   [unit = u"ppb"]   #C2ClF5; CFC-115, Freon 115
    @variables CH2Br2(T)    = 0.976   [unit = u"ppb"]   #CH3Br2; Dibromomethane
    @variables CH2Cl2(T)    = 0.085   [unit = u"ppb"]   #CH2Cl2; Dichloromethane
    @variables CH2I2(T)     = 0.151   [unit = u"ppb"]   #CH2I2; Diiodomethane
    @variables CH2IBr(T)    = 0.412   [unit = u"ppb"]   #CH2IBr; Bromoiodomethane
    @variables CH2ICl(T)    = 0.912   [unit = u"ppb"]   #CH2ICl; Chloroiodomethane
    @variables CH2O(T)      = 0.416   [unit = u"ppb"]   #CH2O; Formaldehyde
    @variables CH2OO(T)     = 0.477   [unit = u"ppb"]   #CH2OO; Criegee intermediate
    @variables CH3Br(T)     = 0.356   [unit = u"ppb"]   #CH3Br; Methyl bromide
    @variables CH3CCl3(T)   = 0.028   [unit = u"ppb"]   #CH3CCl3; Methyl chloroform
    @variables CH3CHOO(T)   = 0.385   [unit = u"ppb"]   #CH3CHOO; Criegee intermediate
    @variables CH3Cl(T)     = 0.603   [unit = u"ppb"]   #CH3Cl; Chloromethane
    @variables CH3I(T)      = 0.449   [unit = u"ppb"]   #CH3I; Methyl iodide
    @variables CH4(T)       = 0.987   [unit = u"ppb"]   #CH4; Methane
    @variables CHBr3(T)     = 0.205   [unit = u"ppb"]   #CHBr3; Tribromethane
    @variables CHCl3(T)     = 0.365   [unit = u"ppb"]   #CHCl3; Chloroform
    @variables Cl(T)        = 0.838   [unit = u"ppb"]   #Cl; Atomic chlorine
    @variables Cl2(T)       = 0.639   [unit = u"ppb"]   #Cl2; Molecular chlorine
    @variables Cl2O2(T)     = 0.663   [unit = u"ppb"]   #Cl2O2; Dichlorine dioxide
    @variables ClNO2(T)     = 0.866   [unit = u"ppb"]   #ClNO2; Nitryl chloride
    @variables ClNO3(T)     = 0.598   [unit = u"ppb"]   #ClONO2; Chlorine nitrate
    @variables ClO(T)       = 0.475   [unit = u"ppb"]   #ClO; Chlorine monoxide
    @variables ClOO(T)      = 0.492   [unit = u"ppb"]   #ClOO; Chlorine dioxide
    @variables CO(T)        = 0.576   [unit = u"ppb"]   #CO; Carbon monoxide
    @variables CO2(T)       = 0.308   [unit = u"ppb"]   #CO2; Carbon dioxide
    @variables CSL(T)       = 0.341   [unit = u"ppb"]   #cresols and xylols
    @variables DMS(T)       = 0.465   [unit = u"ppb"]   #(CH3)2S; Dimethylsulfide
    @variables EOH(T)       = 0.215   [unit = u"ppb"]   #C2H5OH; Ethanol
    @variables ETHLN(T)     = 0.187   [unit = u"ppb"]   #CHOCH2ONO2; Ethanal nitrate
    @variables ETHN(T)      = 0.735   [unit = u"ppb"]   #stable hydroxy-nitrooxy-ethane
    @variables ETHP(T)      = 0.577   [unit = u"ppb"]   #stable hydroxy-hydroperoxy-ethane
    @variables ETNO3(T)     = 0.012   [unit = u"ppb"]   #C2H5ONO2; Ethyl nitrate
    @variables ETO(T)       = 0.246   [unit = u"ppb"]   #hydroxy-alkoxy-ethane radical
    @variables ETOO(T)      = 0.963   [unit = u"ppb"]   #hydroxy-peroxy-ethane radical, formed from ethene + OH
    @variables ETO2(T)      = 0.295   [unit = u"ppb"]   #CH3CH2OO; Ethylperoxy radical
    @variables ETP(T)       = 0.316   [unit = u"ppb"]   #CH3CH2OOH; Ethylhydroperoxide
    @variables GLYC(T)      = 0.773   [unit = u"ppb"]   #HOCH2CHO; Glycoaldehyde
    @variables GLYX(T)      = 0.148   [unit = u"ppb"]   #CHOCHO; Glyoxal
    @variables H(T)         = 0.109   [unit = u"ppb"]   #H; Atomic hydrogen
    @variables H1211(T)     = 0.703   [unit = u"ppb"]   #CBrClF2; H-1211
    @variables H1301(T)     = 0.277   [unit = u"ppb"]   #CBrF3; H-1301
    @variables H2402(T)     = 0.56   [unit = u"ppb"]   #C2Br2F4; H-2402
    @variables H2O(T)       = 0.633   [unit = u"ppb"]   #H2O; Water vapor
    @variables H2O2(T)      = 0.955   [unit = u"ppb"]   #H2O2; Hydrogen peroxide
    @variables HAC(T)       = 0.595   [unit = u"ppb"]   #HOCH2C(O)CH3; Hydroxyacetone
    @variables HBr(T)       = 0.19   [unit = u"ppb"]   #HBr; Hypobromic acid
    @variables HC5A(T)      = 0.723   [unit = u"ppb"]   #C5H8O2; Isoprene-4,1-hydroxyaldehyde
    @variables HCFC123(T)   = 0.797   [unit = u"ppb"]   #C2HCl2F3; HCFC-123, R-123, Freon 123
    @variables HCFC141b(T)  = 0.287   [unit = u"ppb"]   #C(CH3)Cl2F; HCFC-141b, R-141b, Freon 141b
    @variables HCFC142b(T)  = 0.847   [unit = u"ppb"]   #C(CH3)ClF2; HCFC-142b, R-142b, Freon 142b
    @variables HCFC22(T)    = 0.242   [unit = u"ppb"]   #CHClF2 ; HCFC-22, R-22, Freon 22
    @variables HCl(T)       = 0.445   [unit = u"ppb"]   #HCl; Hydrochloric acid
    @variables HCOOH(T)     = 0.826   [unit = u"ppb"]   #HCOOH; Formic acid
    @variables HI(T)        = 0.857   [unit = u"ppb"]   #HI; Hydrogen iodide
    @variables HMHP(T)      = 0.626   [unit = u"ppb"]   #HOCH2OOH; Hydroxymethyl hydroperoxide
    @variables HMML(T)      = 0.286   [unit = u"ppb"]   #C4H6O3; Hydroxymethyl-methyl-a-lactone
    @variables HMS(T)       = 0.51   [unit = u"ppb"]   #HOCH2SO3-; hydroxymethanesulfonate
    @variables HNO2(T)      = 0.743   [unit = u"ppb"]   #HONO; Nitrous acid
    @variables HNO3(T)      = 0.656   [unit = u"ppb"]   #HNO3; Nitric acid
    @variables HNO4(T)      = 0.691   [unit = u"ppb"]   #HNO4; Pernitric acid
    @variables HO2(T)       = 0.064   [unit = u"ppb"]   #HO2; Hydroperoxyl radical
    @variables HOBr(T)      = 0.034   [unit = u"ppb"]   #HOBr; Hypobromous acid
    @variables HOCl(T)      = 0.151   [unit = u"ppb"]   #HOCl; Hypochlorous acid
    @variables HOI(T)       = 0.358   [unit = u"ppb"]   #HOI; Hypoiodous acid
    @variables HONIT(T)     = 0.865   [unit = u"ppb"]   #2nd gen monoterpene organic nitrate
    @variables HPALD1(T)    = 0.604   [unit = u"ppb"]   #O=CHC(CH3)=CHCH2OOH; d-4,1-C5-hydroperoxyaldehyde
    @variables HPALD1OO(T)  = 0.486   [unit = u"ppb"]   #peroxy radicals from HPALD1
    @variables HPALD2(T)    = 0.928   [unit = u"ppb"]   #HOOCH2C(CH3)=CHCH=O; d-1,4-C5-hydroperoxyaldehyde
    @variables HPALD2OO(T)  = 0.317   [unit = u"ppb"]   #peroxy radicals from HPALD2
    @variables HPALD3(T)    = 0.887   [unit = u"ppb"]   #O=CHC(CH3)OOHCH=CH2; b-2,1-C5-hydroperoxyaldehyde
    @variables HPALD4(T)    = 0.639   [unit = u"ppb"]   #CH2=C(CH3)CHOOHCH=O; b-3,4-C5-hydroperoxyaldehyde
    @variables HPETHNL(T)   = 0.551   [unit = u"ppb"]   #CHOCH2OOH; hydroperoxyethanal
    @variables I(T)         = 0.416   [unit = u"ppb"]   #I; Atmoic iodine
    @variables I2(T)        = 0.169   [unit = u"ppb"]   #I2; Molecular iodine
    @variables I2O2(T)      = 0.604   [unit = u"ppb"]   #I2O2; Diiodine dioxide
    @variables I2O3(T)      = 0.976   [unit = u"ppb"]   #I2O3; Diiodine trioxide
    @variables I2O4(T)      = 0.512   [unit = u"ppb"]   #I2O4; Diiodine tetraoxide
    @variables IBr(T)       = 0.437   [unit = u"ppb"]   #IBr; Iodine monobromide
    @variables ICHE(T)      = 0.981   [unit = u"ppb"]   #C5H8O3; Isoprene hydroxy-carbonyl-epoxides
    @variables ICHOO(T)     = 0.26   [unit = u"ppb"]   #peroxy radical from IEPOXD
    @variables ICl(T)       = 0.479   [unit = u"ppb"]   #ICl; Iodine monochloride
    @variables ICN(T)       = 0.487   [unit = u"ppb"]   #C5H7NO4; Lumped isoprene carbonyl nitrates
    @variables ICNOO(T)     = 0.353   [unit = u"ppb"]   #peroxy radicals from ICN
    @variables ICPDH(T)     = 0.678   [unit = u"ppb"]   #C5H10O5; Isoprene dihydroxy hydroperoxycarbonyl
    @variables IDC(T)       = 0.99   [unit = u"ppb"]   #C5H6O2; Lumped isoprene dicarbonyls
    @variables IDCHP(T)     = 0.447   [unit = u"ppb"]   #C5H8O5; Isoprene dicarbonyl hydroxy dihydroperoxide
    @variables IDHDP(T)     = 0.985   [unit = u"ppb"]   #C5H12O6; Isoprene dihydroxy dihydroperoxide
    @variables IDHNBOO(T)   = 0.467   [unit = u"ppb"]   #peroxy radicals from INPB
    @variables IDHNDOO1(T)  = 0.678   [unit = u"ppb"]   #peroxy radicals from INPD
    @variables IDHNDOO2(T)  = 0.255   [unit = u"ppb"]   #peroxy radicals from INPD
    @variables IDHPE(T)     = 0.863   [unit = u"ppb"]   #C5H10O5; Isoprene dihydroxy hydroperoxy epoxide
    @variables IDN(T)       = 0.97   [unit = u"ppb"]   #C5H8N2O6; Lumped isoprene dinitrates
    @variables IDNOO(T)     = 0.601   [unit = u"ppb"]   #peroxy radicals from IDN
    @variables IEPOXA(T)    = 0.532   [unit = u"ppb"]   #C5H10O3; trans-Beta isoprene epoxydiol
    @variables IEPOXAOO(T)  = 0.935   [unit = u"ppb"]   #peroxy radical from trans-Beta isoprene epoxydiol
    @variables IEPOXB(T)    = 0.744   [unit = u"ppb"]   #C5H10O3; cis-Beta isoprene epoxydiol
    @variables IEPOXBOO(T)  = 0.885   [unit = u"ppb"]   #peroxy radical from cis-Beta isoprene epoxydiol
    @variables IEPOXD(T)    = 0.312   [unit = u"ppb"]   #C5H10O3; Delta isoprene epoxydiol
    @variables IHN1(T)      = 0.239   [unit = u"ppb"]   #C5H9NO4; Isoprene-d-4-hydroxy-1-nitrate
    @variables IHN2(T)      = 0.955   [unit = u"ppb"]   #C5H9NO4; Isoprene-b-1-hydroxy-2-nitrate
    @variables IHN3(T)      = 0.876   [unit = u"ppb"]   #C5H9NO4; Isoprene-b-4-hydroxy-3-nitrate
    @variables IHN4(T)      = 0.11   [unit = u"ppb"]   #C5H9NO4; Isoprene-d-1-hydroxy-4-nitrate
    @variables IHOO1(T)     = 0.009   [unit = u"ppb"]   #peroxy radical from OH addition to isoprene at C1
    @variables IHOO4(T)     = 0.419   [unit = u"ppb"]   #peroxy radical from OH addition to isoprene at C4
    @variables IHPNBOO(T)   = 0.299   [unit = u"ppb"]   #peroxy radicals from INPB
    @variables IHPNDOO(T)   = 0.835   [unit = u"ppb"]   #peroxy radicals from INPD
    @variables IHPOO1(T)    = 0.458   [unit = u"ppb"]   #peroxy radical from ISOPOOH
    @variables IHPOO2(T)    = 0.598   [unit = u"ppb"]   #peroxy radical from ISOPOOH
    @variables IHPOO3(T)    = 0.492   [unit = u"ppb"]   #peroxy radical from ISOPOOH
    @variables INA(T)       = 0.776   [unit = u"ppb"]   #alkoxy radical from INO2D
    @variables INDIOL(T)    = 0.073   [unit = u"ppb"]   #Generic aerosol phase organonitrate hydrolysis product
    @variables INO(T)       = 0.954   [unit = u"ppb"]   #INO; Nitrosyl iodide
    @variables INO2B(T)     = 0.904   [unit = u"ppb"]   #beta-peroxy radicals from isoprene + NO3
    @variables INO2D(T)     = 0.936   [unit = u"ppb"]   #delta-peroxy radicals from isoprene + NO3
    @variables INPB(T)      = 0.277   [unit = u"ppb"]   #C5H9NO5; Lumped isoprene beta-hydroperoxy nitrates
    @variables INPD(T)      = 0.929   [unit = u"ppb"]   #C5H9NO5; Lumped isoprene delta-hydroperoxy nitrates
    @variables IO(T)        = 0.185   [unit = u"ppb"]   #IO; Iodine monoxide
    @variables IONITA(T)    = 0.432   [unit = u"ppb"]   #Aerosol-phase organic nitrate from isoprene precursors
    @variables IONO(T)      = 0.548   [unit = u"ppb"]   #IONO; Nitryl iodide
    @variables IONO2(T)     = 0.475   [unit = u"ppb"]   #IONO2; Iodine nitrate
    @variables IPRNO3(T)    = 0.493   [unit = u"ppb"]   #C3H8ONO2; Isopropyl nitrate
    @variables ISALA(T)     = 0.048   [unit = u"ppb"]   #I; Fine sea salt iodine
    @variables ISALC(T)     = 0.511   [unit = u"ppb"]   #I; Coarse sea salt iodine
    @variables ISOP(T)      = 0.222   [unit = u"ppb"]   #CH2=C(CH3)CH=CH2; Isoprene
    @variables ISOPNOO1(T)  = 0.587   [unit = u"ppb"]   #peroxy radicals from IHN2
    @variables ISOPNOO2(T)  = 0.18   [unit = u"ppb"]   #peroxy radicals from IHN3
    @variables ITCN(T)      = 0.895   [unit = u"ppb"]   #C5H9NO7; Lumped tetrafunctional isoprene carbonyl-nitrates
    @variables ITHN(T)      = 0.652   [unit = u"ppb"]   #C5H11NO7; Lumped tetrafunctional isoprene hydroxynitrates
    @variables KO2(T)       = 0.237   [unit = u"ppb"]   #RO2 from >3 ketones
    @variables LBRO2H(T)    = 0.137   [unit = u"ppb"]   #Dummy spc to track oxidation of BRO2 by HO2
    @variables LBRO2N(T)    = 0.395   [unit = u"ppb"]   #Dummy spc to track oxidation of BRO2 by NO
    @variables LIMO(T)      = 0.506   [unit = u"ppb"]   #C10H16; Limonene
    @variables LIMO2(T)     = 0.531   [unit = u"ppb"]   #RO2 from LIMO
    @variables LISOPOH(T)   = 0.238   [unit = u"ppb"]   #Dummy spc to track oxidation of ISOP by OH
    @variables LISOPNO3(T)  = 0.488   [unit = u"ppb"]   #Dummy spc to track oxidation of ISOP by NO3
    @variables LNRO2H(T)    = 0.855   [unit = u"ppb"]   #Dummy spc to track oxidation of NRO2 by HO2
    @variables LNRO2N(T)    = 0.025   [unit = u"ppb"]   #Dummy spc to track oxidation of NRO2 by NO
    @variables LTRO2H(T)    = 0.142   [unit = u"ppb"]   #Dummy spc to track oxidation of TRO2 by HO2
    @variables LTRO2N(T)    = 0.398   [unit = u"ppb"]   #Dummy spc to track oxidation of TRO2 by NO
    @variables LVOC(T)      = 0.79   [unit = u"ppb"]   #C5H14O5; Gas-phase low-volatility non-IEPOX product of ISOPOOH (RIP) oxidation
    @variables LVOCOA(T)    = 0.751   [unit = u"ppb"]   #C5H14O5; Aerosol-phase low-volatility non-IEPOX product of ISOPOOH (RIP) oxidation
    @variables LXRO2H(T)    = 0.179   [unit = u"ppb"]   #Dummy spc to track oxidation of XRO2 by HO2
    @variables LXRO2N(T)    = 0.219   [unit = u"ppb"]   #Dummy spc to track oxidation of XRO2 by NO
    @variables MACR(T)      = 0.484   [unit = u"ppb"]   #CH2=C(CH3)CHO; Methacrolein
    @variables MACR1OO(T)   = 0.131   [unit = u"ppb"]   #peroxyacyl radical from MACR + OH
    @variables MACR1OOH(T)  = 0.498   [unit = u"ppb"]   #CH2=C(CH3)C(O)OOH; Peracid from MACR
    @variables MACRNO2(T)   = 0.96   [unit = u"ppb"]   #Product of MCRHN + OH
    @variables MAP(T)       = 0.134   [unit = u"ppb"]   #CH3C(O)OOH; Peroxyacetic acid
    @variables MCO3(T)      = 0.9   [unit = u"ppb"]   #CH3C(O)OO; Peroxyacetyl radical
    @variables MCRDH(T)     = 0.551   [unit = u"ppb"]   #C4H8O3; Dihydroxy-MACR
    @variables MCRENOL(T)   = 0.116   [unit = u"ppb"]   #C4H6O2; Lumped enols from MVK/MACR
    @variables MCRHN(T)     = 0.733   [unit = u"ppb"]   #HOCH2C(ONO2)(CH3)CHO; Hydroxynitrate from MACR
    @variables MCRHNB(T)    = 0.455   [unit = u"ppb"]   #O2NOCH2C(OH)(CH3)CHO; Hydroxynitrate from MACR
    @variables MCRHP(T)     = 0.373   [unit = u"ppb"]   #HOCH2C(OOH)(CH3)CHO; Hydroxy-hydroperoxy-MACR
    @variables MCROHOO(T)   = 0.132   [unit = u"ppb"]   #peroxy radical from MACR + OH
    @variables MCT(T)       = 0.809   [unit = u"ppb"]   #methylcatechols
    @variables MEK(T)       = 0.398   [unit = u"ppb"]   #RC(O)R; Methyl ethyl ketone
    @variables MENO3(T)     = 0.973   [unit = u"ppb"]   #CH3ONO2; methyl nitrate
    @variables MGLY(T)      = 0.757   [unit = u"ppb"]   #CH3COCHO; Methylglyoxal
    @variables MO2(T)       = 0.291   [unit = u"ppb"]   #CH3O2; Methylperoxy radical
    @variables MOH(T)       = 0.459   [unit = u"ppb"]   #CH3OH; Methanol
    @variables MONITA(T)    = 0.829   [unit = u"ppb"]   #Aerosol-phase organic nitrate from monoterpene precursors
    @variables MONITS(T)    = 0.264   [unit = u"ppb"]   #Saturated 1st gen monoterpene organic nitrate
    @variables MONITU(T)    = 0.432   [unit = u"ppb"]   #Unsaturated 1st gen monoterpene organic nitrate
    @variables MP(T)        = 0.142   [unit = u"ppb"]   #CH3OOH; Methylhydroperoxide
    @variables MPAN(T)      = 0.769   [unit = u"ppb"]   #CH2=C(CH3)C(O)OONO2; Peroxymethacroyl nitrate (PMN)
    @variables MPN(T)       = 0.497   [unit = u"ppb"]   #CH3O2NO2; Methyl peroxy nitrate
    @variables MSA(T)       = 0.795   [unit = u"ppb"]   #CH4SO3; Methanesulfonic acid
    @variables MTPA(T)      = 0.449   [unit = u"ppb"]   #Lumped monoterpenes: a-pinene, b-pinene, sabinene, carene
    @variables MTPO(T)      = 0.976   [unit = u"ppb"]   #Other monoterpenes: Terpinene, terpinolene, myrcene, ocimene, other monoterpenes
    @variables MVK(T)       = 0.433   [unit = u"ppb"]   #CH2=CHC(=O)CH3; Methyl vinyl ketone
    @variables MVKDH(T)     = 0.641   [unit = u"ppb"]   #HOCH2CH2OHC(O)CH3; Dihydroxy-MVK
    @variables MVKHC(T)     = 0.904   [unit = u"ppb"]   #C4H6O3; MVK hydroxy-carbonyl
    @variables MVKHCB(T)    = 0.875   [unit = u"ppb"]   #C4H6O3; MVK hydroxy-carbonyl
    @variables MVKHP(T)     = 0.234   [unit = u"ppb"]   #C4H8O4; MVK hydroxy-hydroperoxide
    @variables MVKN(T)      = 0.092   [unit = u"ppb"]   #HOCH2CH(ONO2)C(=O)CH3; Hydroxynitrate from MVK
    @variables MVKOHOO(T)   = 0.379   [unit = u"ppb"]   #peroxy radical from MVK + OH
    @variables MVKPC(T)     = 0.349   [unit = u"ppb"]   #OCHCH(OOH)C(O)CH3; MVK hydroperoxy-carbonyl
    @variables N(T)         = 0.934   [unit = u"ppb"]   #N; Atomic nitrogen
    @variables N2O(T)       = 0.673   [unit = u"ppb"]   #N2O; Nitrous oxide
    @variables N2O5(T)      = 0.556   [unit = u"ppb"]   #N2O5; Dinitrogen pentoxide
    @variables NAP(T)       = 0.934   [unit = u"ppb"]   #C10H8; Naphthalene; IVOC surrogate
    @variables NIT(T)       = 0.498   [unit = u"ppb"]   #NIT; Fine mode inorganic nitrate
    @variables NITs(T)      = 0.766   [unit = u"ppb"]   #NITs; Coarse mode inorganic nitrate
    @variables NO(T)        = 0.773   [unit = u"ppb"]   #NO; Nitric oxide
    @variables NO2(T)       = 0.377   [unit = u"ppb"]   #NO2; Nitrogen dioxide
    @variables NO3(T)       = 0.502   [unit = u"ppb"]   #NO3; Nitrate radical
    @variables NPHEN(T)     = 0.309   [unit = u"ppb"]   #nitrophenols
    @variables NPRNO3(T)    = 0.348   [unit = u"ppb"]   #C3H8ONO2; n-propyl nitrate
    @variables NRO2(T)      = 0.623   [unit = u"ppb"]   #Peroxy radical from NAP oxidation
    @variables O(T)         = 0.448   [unit = u"ppb"]   #O(3P); Ground state atomic oxygen
    @variables O1D(T)       = 0.212   [unit = u"ppb"]   #O(1D); Excited atomic oxygen
    @variables O3(T)        = 0.239   [unit = u"ppb"]   #O3; Ozone
    @variables O3A(T)       = 0.723   [unit = u"ppb"]   #O3; Ozone in accum seasalt
    @variables O3C(T)       = 0.086   [unit = u"ppb"]   #O3; Ozone in coarse seasalt
    @variables OClO(T)      = 0.616   [unit = u"ppb"]   #OClO; Chlorine dioxide
    @variables OCS(T)       = 0.189   [unit = u"ppb"]   #COS; Carbonyl sulfide
    @variables OH(T)        = 0.384   [unit = u"ppb"]   #OH; Hydroxyl radical
    @variables OIO(T)       = 0.732   [unit = u"ppb"]   #OIO; Iodine dioxide
    @variables OLND(T)      = 0.295   [unit = u"ppb"]   #Monoterpene-derived NO3-alkene adduct
    @variables OLNN(T)      = 0.843   [unit = u"ppb"]   #Monoterpene-derived NO3 adduct
    @variables OTHRO2(T)    = 0.43   [unit = u"ppb"]   #Other C2 RO2 not from C2H6 oxidation
    @variables PAN(T)       = 0.203   [unit = u"ppb"]   #CH3C(O)OONO2; Peroxyacetylnitrate
    @variables PHEN(T)      = 0.687   [unit = u"ppb"]   #phenol
    @variables PIO2(T)      = 0.185   [unit = u"ppb"]   #RO2 from MTPA
    @variables PIP(T)       = 0.144   [unit = u"ppb"]   #Peroxides from MTPA
    @variables PO2(T)       = 0.798   [unit = u"ppb"]   #HOCH2CH(OO)CH3; RO2 from propene
    @variables PP(T)        = 0.083   [unit = u"ppb"]   #HOCH2CH(OOH)CH3; Peroxide from PO2
    @variables PPN(T)       = 0.466   [unit = u"ppb"]   #CH3CH2C(O)OONO2; Peroxypropionylnitrate
    @variables PRN1(T)      = 0.038   [unit = u"ppb"]   #O2NOCH2CH(OO)CH3; RO2 from propene + NO3
    @variables PROPNN(T)    = 0.574   [unit = u"ppb"]   #CH3C(=O)CH2ONO2; Propanone nitrate
    @variables PRPE(T)      = 0.29   [unit = u"ppb"]   #C3H6; >= C3 alkenes
    @variables PRPN(T)      = 0.351   [unit = u"ppb"]   #O2NOCH2CH(OOH)CH3; Peroxide from PRN1
    @variables PYAC(T)      = 0.709   [unit = u"ppb"]   #CH3COCOOH; Pyruvic acid
    @variables R4N1(T)      = 0.564   [unit = u"ppb"]   #RO2 from R4N2
    @variables R4N2(T)      = 0.002   [unit = u"ppb"]   #RO2NO; >= C4 alkylnitrates
    @variables R4O2(T)      = 0.365   [unit = u"ppb"]   #RO2 from ALK4
    @variables R4P(T)       = 0.383   [unit = u"ppb"]   #CH3CH2CH2CH2OOH; Peroxide from R4O2
    @variables RA3P(T)      = 0.084   [unit = u"ppb"]   #CH3CH2CH2OOH; Peroxide from A3O2
    @variables RB3P(T)      = 0.392   [unit = u"ppb"]   #CH3CH(OOH)CH3; Peroxide from B3O2
    @variables RCHO(T)      = 0.93   [unit = u"ppb"]   #CH3CH2CHO; >= C3 aldehydes
    @variables RCO3(T)      = 0.144   [unit = u"ppb"]   #CH3CH2C(O)OO; Peroxypropionyl radical
    @variables RIPA(T)      = 0.786   [unit = u"ppb"]   #HOCH2C(OOH)(CH3)CH=CH2; 1,2-ISOPOOH
    @variables RIPB(T)      = 0.664   [unit = u"ppb"]   #HOCH2C(OOH)(CH3)CH=CH2; 4,3-ISOPOOH
    @variables RIPC(T)      = 0.69   [unit = u"ppb"]   #C5H10O3; d(1,4)-ISOPOOH
    @variables RIPD(T)      = 0.743   [unit = u"ppb"]   #C5H10O3; d(4,1)-ISOPOOH
    @variables ROH(T)       = 0.089   [unit = u"ppb"]   #C3H7OH; > C2 alcohols
    @variables RP(T)        = 0.787   [unit = u"ppb"]   #CH3CH2C(O)OOH; Peroxide from RCO3
    @variables SALAAL(T)    = 0.927   [unit = u"ppb"]   #Accumulation mode seasalt aerosol alkalinity
    @variables SALCAL(T)    = 0.472   [unit = u"ppb"]   #Coarse mode seasalt aerosol alkalinity
    @variables SALACL(T)    = 0.263   [unit = u"ppb"]   #Cl; Fine chloride
    @variables SALCCL(T)    = 0.188   [unit = u"ppb"]   #Cl; Coarse chloride
    @variables SALASO2(T)   = 0.469   [unit = u"ppb"]   #SO2; Fine seasalt
    @variables SALCSO2(T)   = 0.819   [unit = u"ppb"]   #SO2; Coarse seasalt
    @variables SALASO3(T)   = 0.816   [unit = u"ppb"]   #SO3--; Fine seasalt
    @variables SALCSO3(T)   = 0.014   [unit = u"ppb"]   #SO3--; Coarse chloride
    @variables SO2(T)       = 0.156   [unit = u"ppb"]   #SO2; Sulfur dioxide
    @variables SO4(T)       = 0.455   [unit = u"ppb"]   #SO4; Sulfate
    @variables SO4s(T)      = 0.804   [unit = u"ppb"]   #SO4 on sea-salt; Sulfate
    @variables SOAGX(T)     = 0.606   [unit = u"ppb"]   #CHOCHO; Aerosol-phase glyoxal
    @variables SOAIE(T)     = 0.605   [unit = u"ppb"]   #C5H10O3; Aerosol-phase IEPOX
    @variables TOLU(T)      = 0.931   [unit = u"ppb"]   #C7H8; Toluene
    @variables TRO2(T)      = 0.123   [unit = u"ppb"]   #Peroxy radical from TOLU oxidation
    @variables XYLE(T)      = 0.246   [unit = u"ppb"]   #C8H10; Xylene
    @variables XRO2(T)      = 0.524   [unit = u"ppb"]   #Peroxy radical from XYLE oxidation
    @variables H2(T)        = 0.613   [unit = u"ppb"]   #H2; Molecular hydrogen
    @variables N2(T)        = 0.012   [unit = u"ppb"]   #N2; Molecular nitrogen
    @variables O2(T)        = 0.937   [unit = u"ppb"]   #O2; Molecular oxygen
    @variables RCOOH(T)     = 0.558   [unit = u"ppb"]   #C2H5C(O)OH; > C2 organic acids 


    rxs = [

        Reaction(K_MT(1), [SO2, SALAAL, O3], [SO4 - SALAAL], [1, 1, 1], [1])  #  SO2  + SALAAL + O3  = SO4 - SALAAL         
        Reaction(K_MT(2), [HCl, SALAAL], [SALACL], [1, 1], [1])  #  HCl  + SALAAL       = SALACL               
        Reaction(K_MT(3), [HNO3, SALAAL], [NIT], [1, 1], [1])  #  HNO3 + SALAAL       = NIT                  
        Reaction(K_MT(4), [SO2, SALCAL, O3], [SO4s - SALCAL], [1, 1, 1], [1])  #  SO2  + SALCAL + O3  = SO4s - SALCAL        
        Reaction(K_MT(5), [HCl, SALCAL], [SALCCL], [1, 1], [1])  #  HCl  + SALCAL       = SALCCL               
        Reaction(K_MT(6), [HNO3, SALCAL], [NITs], [1, 1], [1])  #  HNO3 + SALCAL       = NITs                 
        Reaction(K_CLD(1), [SO2, H2O2], [SO4], [1, 1], [1])  #  SO2 + H2O2          = SO4                  
        Reaction(K_CLD(2), [SO2, O3], [SO4], [1, 1], [1])  #  SO2 + O3            = SO4                  
        Reaction(K_CLD(3), [#, O2], [SO4], [1, 1], [1])   #Mn & Fe catalysis + HET_DROP_CHEM()  #     #+O2           = SO4                  
        Reaction(K_CLD(4), [CH2O, SO2], [HMS], [1, 1], [1])   #Sep 2021; Moch2020; MSL  #  CH2O + SO2          = HMS                  
        Reaction(K_CLD(5), [HMS], [SO2, CH2O], [1], [1, 1])   #Sep 2021; Moch2020; MSL  #  HMS                 = SO2 + CH2O           
        Reaction(K_CLD(6), [HMS, OH], [SO4, CH2O - SO2], [1, 1], [2, 1])   #Sep 2021; Moch2020; MSL  #  HMS + OH            = 2SO4 + CH2O - SO2    


        Reaction(GCARR_ac(3.00e-12, -1500.0e0), [O3, NO], [NO2, O2], [1, 1], [1, 1])  #  O3 + NO = NO2 + O2   
        Reaction(GCARR_ac(1.70e-12, -940.0e0), [O3, OH], [HO2, O2], [1, 1], [1, 1])  #  O3 + OH = HO2 + O2   
        Reaction(GCARR_ac(1.00e-14, -490.0e0), [O3, HO2], [OH, O2, O2], [1, 1], [1, 1, 1])  #  O3 + HO2 = OH + O2 + O2   
        Reaction(GCARR_ac(1.20e-13, -2450.0e0), [O3, NO2], [O2, NO3], [1, 1], [1, 1])  #  O3 + NO2 = O2 + NO3   
        Reaction(GCARR_ac(2.90e-16, -1000.0e0), [O3, MO2], [CH2O, HO2, O2], [1, 1], [1, 1, 1])   #2014/02/03; Eastham2014; SDE  #  O3 + MO2 = CH2O + HO2 + O2   
        Reaction(1.80e-12, [OH, OH], [H2O, O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + OH = H2O + O   
        Reaction(GCJPLPR_aba(6.90e-31, 1.0e+00, 2.6e-11, 0.6e0), [#, M], [H2O2], [1, 1], [1])                        GCJPLPR_aba(6.90d-31, 1.0d+00, 2.6d-11, 0.6d0);  #     #+M = H2O2   
        Reaction(GCARR_ac(4.80e-11, 250.0e0), [OH, HO2], [H2O, O2], [1, 1], [1, 1])  #  OH + HO2 = H2O + O2   
        Reaction(1.80e-12, [OH, H2O2], [H2O, HO2], [1, 1], [1, 1])  #  OH + H2O2 = H2O + HO2   
        Reaction(GCARR_ac(3.30e-12, 270.0e0), [HO2, NO], [OH, NO2], [1, 1], [1, 1])   #2013/02/12; JPL 10-6; BHH,JMAO,EAM  #  HO2 + NO = OH + NO2   
        Reaction(GC_HO2HO2_acac(3.00e-13, 460.0e0, 2.1e-33, 920.0e0), [HO2, HO2], [H2O2, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  HO2 + HO2 = H2O2 + O2   
        Reaction(GC_OHCO_a(1.50e-13), [OH, CO], [HO2, CO2], [1, 1], [1, 1])   #2017/02/22; JPL 15-10; BHH,MJE  #  OH + CO = HO2 + CO2   
        Reaction(GCARR_ac(2.45e-12, -1775.0e0), [OH, CH4], [MO2, H2O], [1, 1], [1, 1])  #  OH + CH4 = MO2 + H2O   
        Reaction(GC_RO2NO_B1_ac(2.80e-12, 300.0e0), [MO2, NO], [CH2O, HO2, NO2], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF  #  MO2 + NO = CH2O + HO2 + NO2   
        Reaction(GC_RO2NO_A1_ac(2.80e-12, 300.0e0), [MO2, NO], [MENO3], [1, 1], [1])   #2019/05/10; Fisher2018; JAF  #  MO2 + NO = MENO3   
        Reaction(GCARR_abc(4.10e-13, 0.0e0, 750.0e0), [MO2, HO2], [MP, O2], [1, 1], [1, 1])  #  MO2 + HO2 = MP + O2   
        Reaction(GC_TBRANCH_1_acac(9.50e-14, 390.0e0, 2.62e1, -1130.0e0), [MO2, MO2], [MOH, CH2O, O2], [1, 1], [1, 1, 1])  #  MO2 + MO2 = MOH + CH2O + O2   
        Reaction(GC_TBRANCH_1_acac(9.50e-14, 390.0e0, 4.0e-2, 1130.0e0), [MO2, MO2], [CH2O, HO2], [1, 1], [2.000, 2.000])  #  MO2 + MO2 = 2.000CH2O + 2.000HO2   
        Reaction(1.60e-10 , [MO2, OH], [MOH, CH2O, HO2], [1, 1], [0.13, 0.87, 1.74])   #2021/09/22; Bates2021a; KHB,MSL  #  MO2 + OH = 0.13MOH + 0.87CH2O + 1.74HO2   
        Reaction(GCARR_ac(2.66e-12, 200.0e0), [MP, OH], [MO2, H2O], [1, 1], [1, 1])  #  MP + OH = MO2 + H2O   
        Reaction(GCARR_ac(1.14e-12, 200.0e0), [MP, OH], [CH2O, OH, H2O], [1, 1], [1, 1, 1])  #  MP + OH = CH2O + OH + H2O   
        Reaction(GCARR_ac(2.66e-12, 200.0e0), [ATOOH, OH], [ATO2, H2O], [1, 1], [1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  ATOOH + OH = ATO2 + H2O   
        Reaction(GCARR_ac(1.14e-12, 200.0e0), [ATOOH, OH], [MGLY, OH, H2O], [1, 1], [1, 1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  ATOOH + OH = MGLY + OH + H2O   
        Reaction(GCARR_ac(5.50e-12, 125.0e0), [CH2O, OH], [CO, HO2, H2O], [1, 1], [1, 1, 1])  #  CH2O + OH = CO + HO2 + H2O   
        Reaction(GC_OHHNO3_acacac(2.41e-14, 460.0e0, 2.69e-17, 2199.0e0, 6.51e-34, 1335.0e0), [HNO3, OH], [H2O, NO3], [1, 1], [1, 1])  #  HNO3 + OH = H2O + NO3   
        Reaction(GCARR_ac(1.80e-11, -390.0e0), [HNO2, OH], [H2O, NO2], [1, 1], [1, 1])  #  HNO2 + OH = H2O + NO2   
        Reaction(GCJPLPR_abcabc(9.05e-05, 3.4e0, -10900.0e0, 1.90e15, 0.3e0, -10900.0e0, 0.6e0), [#, M], [HO2, NO2], [1, 1], [1, 1])   #2017/02/22; JPL 15-10; BHH,MJE  #     #+M = HO2 + NO2   
        Reaction(GCARR_ac(1.30e-12, 380.0e0), [HNO4, OH], [H2O, NO2, O2], [1, 1], [1, 1, 1])  #  HNO4 + OH = H2O + NO2 + O2   
        Reaction(3.50e-12, [HO2, NO3], [OH, NO2, O2], [1, 1], [1, 1, 1])  #  HO2 + NO3 = OH + NO2 + O2   
        Reaction(GCARR_ac(1.50e-11, 170.0e0), [NO, NO3], [NO2], [1, 1], [2.000])  #  NO + NO3 = 2.000NO2   
        Reaction(2.20e-11, [OH, NO3], [HO2, NO2], [1, 1], [1, 1])  #  OH + NO3 = HO2 + NO2   
        Reaction(GCJPLPR_abcabc(4.14e-04, 3.0e0, -10840.0e0, 2.76e14, -0.1e0, -10840.0e0, 0.6e0), [#, M], [NO2, NO3], [1, 1], [1, 1])   #2017/02/22; JPL 15-10; BHH,MJE  #     #+M = NO2 + NO3   
        Reaction(4.00e-13, [HCOOH, OH], [H2O, CO2, HO2], [1, 1], [1, 1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  HCOOH + OH = H2O + CO2 + HO2   
        Reaction(GCARR_ac(2.90e-12, -345.0e0), [MOH, OH], [HO2, CH2O], [1, 1], [1, 1])  #  MOH + OH = HO2 + CH2O   
        Reaction(GCARR_ac(4.50e-14, -1260.0e0), [NO2, NO3], [NO, NO2, O2], [1, 1], [1, 1, 1])  #  NO2 + NO3 = NO + NO2 + O2   
        Reaction(5.80e-16, [NO3, CH2O], [HNO3, HO2, CO], [1, 1], [1, 1, 1])  #  NO3 + CH2O = HNO3 + HO2 + CO   
        Reaction(GCARR_ac(1.40e-12, -1900.0e0), [ALD2, NO3], [HNO3, MCO3], [1, 1], [1, 1])  #  ALD2 + NO3 = HNO3 + MCO3   
        Reaction(GCJPLPR_abab(9.70e-29, 5.6e+00, 9.3e-12, 1.5e0, 0.6e0), [#, M], [PAN], [1, 1], [1])   #JPL Eval 17  #     #+M = PAN   
        Reaction(GCJPLEQ_acabab(9.30e-29, 14000.0e0, 9.7e-29, 5.6e0, 9.3e-12, 1.5e0, 0.6e0), [PAN], [MCO3, NO2], [1], [1, 1])  #  PAN = MCO3 + NO2   
        Reaction(GCARR_ac(8.10e-12, 270.0e0), [MCO3, NO], [MO2, NO2, CO2], [1, 1], [1, 1, 1])  #  MCO3 + NO = MO2 + NO2 + CO2   
        Reaction(GCARR_ac(7.66e-12, -1020.0e0), [C2H6, OH], [ETO2, H2O], [1, 1], [1, 1])   #2013/02/12; JPL 10-6; BHH,JMAO,EAM  #  C2H6 + OH = ETO2 + H2O   
        Reaction(GC_RO2NO_B2_aca(2.60e-12, 365.0e0, 2.0e0), [ETO2, NO], [ALD2, NO2, HO2], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF  #  ETO2 + NO = ALD2 + NO2 + HO2   
        Reaction(GC_RO2NO_A2_aca(2.60e-12, 365.0e0, 2.0e0), [ETO2, NO], [ETNO3], [1, 1], [1])   #2019/05/10; Fisher2018; JAF  #  ETO2 + NO = ETNO3   
        Reaction(GCARR_ac(2.60e-12, 365.0e0), [OTHRO2, NO], [ALD2, NO2, HO2], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF  #  OTHRO2 + NO = ALD2 + NO2 + HO2   
        Reaction(GC_TBRANCH_2_acabc(7.60e-12, -585.0e0, 5.87e0, 0.64e0, -816.0e0), [C3H8, OH], [B3O2], [1, 1], [1])  #  C3H8 + OH = B3O2   
        Reaction(GC_TBRANCH_2_acabc(7.60e-12, -585.0e0, 1.7e-1, -0.64e0, 816.0e0), [C3H8, OH], [A3O2], [1, 1], [1])  #  C3H8 + OH = A3O2   
        Reaction(GC_RO2NO_B2_aca(2.90e-12, 350.0e0, 3.0e0), [A3O2, NO], [NO2, HO2, RCHO], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF  #  A3O2 + NO = NO2 + HO2 + RCHO   
        Reaction(GC_RO2NO_A2_aca(2.90e-12, 350.0e0, 3.0e0), [A3O2, NO], [NPRNO3], [1, 1], [1])   #2019/05/10; Fisher2018; JAF  #  A3O2 + NO = NPRNO3   
        Reaction(GCARR_ac(2.70e-12, 350.0e0), [PO2, NO], [NO2, HO2, CH2O, ALD2], [1, 1], [1, 1, 1, 1])  #  PO2 + NO = NO2 + HO2 + CH2O + ALD2   
        Reaction(GCARR_ac(9.10e-12, -405.0e0), [ALK4, OH], [R4O2], [1, 1], [1])  #  ALK4 + OH = R4O2   
        Reaction(GC_RO2NO_A2_aca(2.70e-12, 350.0e0, 4.5e0), [R4O2, NO], [R4N2], [1, 1], [1])  #  R4O2 + NO = R4N2   
        Reaction(GCARR_ac(2.80e-12, 300.0e0), [ATO2, NO], [NO2, CH2O, MCO3], [1, 1], [1, 1, 1])   #2017/07/27; Fix C creation; SAS,BHH,MJE  #  ATO2 + NO = NO2 + CH2O + MCO3   
        Reaction(GC_RO2NO_B2_aca(2.70e-12, 360.0e0, 3.0e0), [B3O2, NO], [NO2, HO2, ACET], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF  #  B3O2 + NO = NO2 + HO2 + ACET   
        Reaction(GC_RO2NO_A2_aca(2.70e-12, 360.0e0, 3.0e0), [B3O2, NO], [IPRNO3], [1, 1], [1])   #2019/05/10; Fisher2018; JAF  #  B3O2 + NO = IPRNO3   
        Reaction(GCARR_ac(2.70e-12, 350.0e0), [PRN1, NO], [NO2, CH2O, ALD2], [1, 1], [2.000, 1, 1])  #  PRN1 + NO = 2.000NO2 + CH2O + ALD2   
        Reaction(GCARR_ac(2.80e-12, -3280.0e0), [ALK4, NO3], [HNO3, R4O2], [1, 1], [1, 1])  #  ALK4 + NO3 = HNO3 + R4O2   
        Reaction(1.60e-12, [R4N2, OH], [R4N1, H2O], [1, 1], [1, 1])  #  R4N2 + OH = R4N1 + H2O   
        Reaction(GCARR_ac(3.15e-14, 920.0e0), [ACTA, OH], [MO2, CO2, H2O], [1, 1], [1, 1, 1])   #2013/02/12; JPL 10-6; BHH,JMAO,EAM  #  ACTA + OH = MO2 + CO2 + H2O   
        Reaction(GCARR_ac(6.00e-12, 410.0e0), [OH, RCHO], [RCO3, H2O], [1, 1], [1, 1])  #  OH + RCHO = RCO3 + H2O   
        Reaction(GCJPLPR_abab(9.00e-28, 8.9e0, 7.7e-12, 0.2e0, 0.6e0), [#, M], [PPN], [1, 1], [1])   #JPL Eval 17  #     #+M = PPN   
        Reaction(GCJPLEQ_acabab(9.00e-29, 14000.0e0, 9.00e-28, 8.9e0, 7.7e-12, 0.2e0, 0.6e0), [PPN], [RCO3, NO2], [1], [1, 1])  #  PPN = RCO3 + NO2   
        Reaction(6.50e-15, [RCHO, NO3], [HNO3, RCO3], [1, 1], [1, 1])  #  RCHO + NO3 = HNO3 + RCO3   
        Reaction(1.33d-13 + 3.82d-11*exp(-2000.0e0/TEMP), [ACET, OH], [ATO2, H2O], [1, 1], [1, 1])   #JPL Eval 17, p1-62-D31; EVF  #  ACET + OH = ATO2 + H2O   
        Reaction(GCARR_ac(7.40e-13, 700.0e0), [R4O2, HO2], [R4P], [1, 1], [1])  #  R4O2 + HO2 = R4P   
        Reaction(GCARR_ac(7.40e-13, 700.0e0), [R4N1, HO2], [R4N2], [1, 1], [1])  #  R4N1 + HO2 = R4N2   
        Reaction(GC_RO2HO2_aca(2.91e-13, 1300.0e0, 3.0e0), [B3O2, HO2], [RB3P], [1, 1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  B3O2 + HO2 = RB3P   
        Reaction(GC_RO2HO2_aca(2.91e-13, 1300.0e0, 3.0e0), [PRN1, HO2], [PRPN], [1, 1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  PRN1 + HO2 = PRPN   
        Reaction(GCARR_ac(1.30e-12, -25.0e0), [MEK, OH], [KO2, H2O], [1, 1], [1, 1])  #  MEK + OH = KO2 + H2O   
        Reaction(8.00e-16, [MEK, NO3], [HNO3, KO2], [1, 1], [1, 1])  #  MEK + NO3 = HNO3 + KO2   
        Reaction(3.35e-12, [EOH, OH], [HO2, ALD2], [1, 1], [1, 1])   #2013/02/12; JPL 10-6; BHH,JMAO,EAM  #  EOH + OH = HO2 + ALD2   
        Reaction(GCARR_ac(4.60e-12, 70.0e0), [ROH, OH], [HO2, RCHO], [1, 1], [1, 1])  #  ROH + OH = HO2 + RCHO   
        Reaction(4.10e-14, [ETO2, ETO2], [ALD2, HO2], [1, 1], [2.000, 2.000])  #  ETO2 + ETO2 = 2.000ALD2 + 2.000HO2   
        Reaction(4.10e-14, [OTHRO2, OTHRO2], [ALD2, HO2], [1, 1], [2.000, 2.000])   #2019/05/10; Fisher2018; JAF  #  OTHRO2 + OTHRO2 = 2.000ALD2 + 2.000HO2   
        Reaction(2.70e-14, [ETO2, ETO2], [EOH, ALD2], [1, 1], [1, 1])  #  ETO2 + ETO2 = EOH + ALD2   
        Reaction(2.70e-14, [OTHRO2, OTHRO2], [EOH, ALD2], [1, 1], [1, 1])   #2019/05/10; Fisher2018; JAF  #  OTHRO2 + OTHRO2 = EOH + ALD2   
        Reaction(GCARR_ac(7.40e-13, 700.0e0), [HO2, ETO2], [ETP], [1, 1], [1])  #  HO2 + ETO2 = ETP   
        Reaction(GCARR_ac(7.40e-13, 700.0e0), [HO2, OTHRO2], [ETP], [1, 1], [1])   #2019/05/10; Fisher2018; JAF  #  HO2 + OTHRO2 = ETP   
        Reaction(GC_RO2HO2_aca(2.91e-13, 1300.0e0, 3.0e0), [A3O2, HO2], [RA3P], [1, 1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  A3O2 + HO2 = RA3P   
        Reaction(GC_RO2HO2_aca(2.91e-13, 1300.0e0, 3.0e0), [PO2, HO2], [PP], [1, 1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  PO2 + HO2 = PP   
        Reaction(GCJPLPR_abab(4.60e-27, 4.0e0, 2.6e-11, 1.3e0, 0.5e0), [#, M], [PO2], [1, 1], [1])   #2017/02/22; JPL 15-10; BHH,MJE  #     #+M = PO2   
        Reaction(GC_GLYCOH_B_a(8.00e-12), [GLYC, OH], [HCOOH, OH, CO], [1, 1], [1, 1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  GLYC + OH = HCOOH + OH + CO   
        Reaction(GCARR_ac(4.59e-13, -1156.0e0), [PRPE, NO3], [PRN1], [1, 1], [1])  #  PRPE + NO3 = PRN1   
        Reaction(GCARR_ac(3.10e-12, 340.0e0), [GLYX, OH], [HO2, CO], [1, 1], [1, 2.000])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  GLYX + OH = HO2 + 2.000CO   
        Reaction(1.50e-11, [MGLY, OH], [MCO3, CO], [1, 1], [1, 1])  #  MGLY + OH = MCO3 + CO   
        Reaction(GC_GLYXNO3_ac(1.40e-12, -1860.0e0), [GLYX, NO3], [HNO3, HO2, CO], [1, 1], [1, 1, 2.000])  #  GLYX + NO3 = HNO3 + HO2 + 2.000CO   
        Reaction(GCARR_ac(3.36e-12, -1860.0e0), [MGLY, NO3], [HNO3, CO, MCO3], [1, 1], [1, 1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  MGLY + NO3 = HNO3 + CO + MCO3   
        Reaction(GC_HACOH_A_ac(2.15e-12, 305.0e0), [HAC, OH], [MGLY, HO2], [1, 1], [1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  HAC + OH = MGLY + HO2   
        Reaction(GCARR_ac(1.68e-12, 500.0e0), [MCO3, A3O2], [MO2, RCHO, HO2], [1, 1], [1, 1, 1])  #  MCO3 + A3O2 = MO2 + RCHO + HO2   
        Reaction(GCARR_ac(1.68e-12, 500.0e0), [MCO3, PO2], [MO2, ALD2, CH2O, HO2], [1, 1], [1, 1, 1, 1])  #  MCO3 + PO2 = MO2 + ALD2 + CH2O + HO2   
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [MCO3, A3O2], [ACTA, RCHO], [1, 1], [1, 1])  #  MCO3 + A3O2 = ACTA + RCHO   
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [MCO3, PO2], [ACTA, RCHO, HAC], [1, 1], [1, 0.350, 0.650])  #  MCO3 + PO2 = ACTA + 0.350RCHO + 0.650HAC   
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [RCO3, MO2], [RCOOH, CH2O], [1, 1], [1, 1])  #  RCO3 + MO2 = RCOOH + CH2O   
        Reaction(GCARR_ac(8.78e-12, 200.0e0), [R4P, OH], [OH, R4O2, RCHO], [1, 1], [0.791, 0.209, 0.791])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  R4P + OH = 0.791OH + 0.209R4O2 + 0.791RCHO   
        Reaction(GCARR_ac(6.13e-13, 200.0e0), [RP, OH], [RCO3], [1, 1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  RP + OH = RCO3   
        Reaction(GCARR_ac(8.78e-12, 200.0e0), [PP, OH], [OH, PO2, HAC], [1, 1], [0.791, 0.209, 0.791])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  PP + OH = 0.791OH + 0.209PO2 + 0.791HAC   
        Reaction(GCARR_ac(4.82e-11, -400.0e0), [LVOC, OH], [OH], [1, 1], [1])   #2017/06/14; Marais2016; EAM  #  LVOC + OH = OH   
        Reaction(GCARR_ac(6.13e-13, 200.0e0), [OH, MAP], [MCO3], [1, 1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  OH + MAP = MCO3   
        Reaction(1.40e-18, [C2H6, NO3], [ETO2, HNO3], [1, 1], [1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  C2H6 + NO3 = ETO2 + HNO3   
        Reaction(GCARR_ac(2.50e-12, 500.0e0), [MCO3, MCO3], [MO2], [1, 1], [2.000])  #  MCO3 + MCO3 = 2.000MO2   
        Reaction(GCARR_ac(1.80e-12, 500.0e0), [MCO3, MO2], [CH2O, MO2, HO2], [1, 1], [1, 1, 1])  #  MCO3 + MO2 = CH2O + MO2 + HO2   
        Reaction(GCARR_ac(2.00e-13, 500.0e0), [MCO3, MO2], [ACTA, CH2O], [1, 1], [1, 1])  #  MCO3 + MO2 = ACTA + CH2O   
        Reaction(GCARR_ac(1.68e-12, 500.0e0), [ATO2, MCO3], [MO2, MCO3, CH2O], [1, 1], [1, 1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  ATO2 + MCO3 = MO2 + MCO3 + CH2O   
        Reaction(GCARR_ac(1.68e-12, 500.0e0), [KO2, MCO3], [MO2, ALD2, MCO3], [1, 1], [1, 1, 1])  #  KO2 + MCO3 = MO2 + ALD2 + MCO3   
        Reaction(GCARR_ac(1.68e-12, 500.0e0), [B3O2, MCO3], [MO2, HO2, ACET], [1, 1], [1, 1, 1])  #  B3O2 + MCO3 = MO2 + HO2 + ACET   
        Reaction(GCARR_ac(1.68e-12, 500.0e0), [PRN1, MCO3], [MO2, NO2, CH2O, ALD2], [1, 1], [1, 1, 1, 1])  #  PRN1 + MCO3 = MO2 + NO2 + CH2O + ALD2   
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [R4O2, MCO3], [MEK, ACTA], [1, 1], [1, 1])  #  R4O2 + MCO3 = MEK + ACTA   
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [ATO2, MCO3], [MGLY, ACTA], [1, 1], [1, 1])   #2017/07/27; Fix C creation; SAS,BHH,MJE  #  ATO2 + MCO3 = MGLY + ACTA   
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [KO2, MCO3], [MEK, ACTA], [1, 1], [1, 1])  #  KO2 + MCO3 = MEK + ACTA   
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [R4N1, MCO3], [RCHO, ACTA, NO2], [1, 1], [1, 1, 1])  #  R4N1 + MCO3 = RCHO + ACTA + NO2   
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [PRN1, MCO3], [RCHO, ACTA, NO2], [1, 1], [1, 1, 1])  #  PRN1 + MCO3 = RCHO + ACTA + NO2   
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [B3O2, MCO3], [ACET, ACTA], [1, 1], [1, 1])  #  B3O2 + MCO3 = ACET + ACTA   
        Reaction(GCARR_ac(1.68e-12, 500.0e0), [MCO3, ETO2], [MO2, ALD2, HO2], [1, 1], [1, 1, 1])  #  MCO3 + ETO2 = MO2 + ALD2 + HO2   
        Reaction(GCARR_ac(1.68e-12, 500.0e0), [MCO3, OTHRO2], [MO2, ALD2, HO2], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF  #  MCO3 + OTHRO2 = MO2 + ALD2 + HO2   
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [MCO3, ETO2], [ACTA, ALD2], [1, 1], [1, 1])  #  MCO3 + ETO2 = ACTA + ALD2   
        Reaction(GCARR_ac(1.87e-13, 500.0e0), [MCO3, OTHRO2], [ACTA, ALD2], [1, 1], [1, 1])   #2019/05/10; Fisher2018; JAF  #  MCO3 + OTHRO2 = ACTA + ALD2   
        Reaction(GCARR_ac(8.50e-13, -2450.0e0), [NO3, NO3], [NO2, O2], [1, 1], [2.000, 1])  #  NO3 + NO3 = 2.000NO2 + O2   
        Reaction(GCJPLPR_abcabc(1.05e-02, 4.8e+00, -11234.0e0, 7.58e16, 2.1e0, -11234.0e0, 0.6e0), [#, M], [MO2, NO2], [1, 1], [1, 1])   #2012/02/12; Browne2011; ECB  #     #+M = MO2 + NO2   
        Reaction(GCARR_ac(1.20e-11, -280.0e0), [DMS, OH], [SO2, MO2, CH2O], [1, 1], [1, 1, 1])  #  DMS + OH = SO2 + MO2 + CH2O   
        Reaction(GC_DMSOH_acac(8.20e-39, 5376.0e0, 1.05e-5, 3644.0e0), [DMS, OH], [SO2, MSA, MO2], [1, 1], [0.750, 0.250, 1])  #  DMS + OH = 0.750SO2 + 0.250MSA + MO2   
        Reaction(GCARR_ac(1.90e-13, 530.0e0), [DMS, NO3], [SO2, HNO3, MO2, CH2O], [1, 1], [1, 1, 1, 1])  #  DMS + NO3 = SO2 + HNO3 + MO2 + CH2O   
        Reaction(GCJPLPR_aba(3.30e-31, 4.3e+00, 1.6e-12, 0.6e0), [#, M], [SO4, HO2], [1, 1], [1, 1])                  GCJPLPR_aba(3.30d-31, 4.3d+00, 1.6d-12, 0.6d0);  #     #+M = SO4 + HO2   
        Reaction(GCARR_ac(1.60e-11, -780.0e0), [Br, O3], [BrO, O2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  Br + O3 = BrO + O2   
        Reaction(GCARR_ac(4.50e-12, 460.0e0), [BrO, HO2], [HOBr, O2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  BrO + HO2 = HOBr + O2   
        Reaction(GCARR_ac(4.80e-12, -310.0e0), [Br, HO2], [HBr, O2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  Br + HO2 = HBr + O2   
        Reaction(GCARR_ac(5.50e-12, 200.0e0), [HBr, OH], [Br, H2O], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  HBr + OH = Br + H2O   
        Reaction(GCARR_ac(2.40e-12,  40.0e0), [BrO, BrO], [Br, O2], [1, 1], [2.000, 1])   #2012/06/07; Parrella2012; JPP  #  BrO + BrO = 2.000Br + O2   
        Reaction(GCARR_ac(2.80e-14, 860.0e0), [BrO, BrO], [Br2, O2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  BrO + BrO = Br2 + O2   
        Reaction(GCARR_ac(8.80e-12, 260.0e0), [BrO, NO], [Br, NO2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  BrO + NO = Br + NO2   
        Reaction(4.90e-11, [Br, BrNO3], [Br2, NO3], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  Br + BrNO3 = Br2 + NO3   
        Reaction(GCARR_ac(2.10e-11, 240.0e0), [Br2, OH], [HOBr, Br], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  Br2 + OH = HOBr + Br   
        Reaction(GCARR_ac(1.20e-10, -430.0e0), [HOBr, O], [OH, BrO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  HOBr + O = OH + BrO   
        Reaction(GCARR_ac(5.80e-12, -1500.0e0), [HBr, O], [OH, Br], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  HBr + O = OH + Br   
        Reaction(GCARR_ac(1.70e-11, 250.0e0), [BrO, OH], [Br, HO2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  BrO + OH = Br + HO2   
        Reaction(1.60e-11, [Br, NO3], [BrO, NO2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  Br + NO3 = BrO + NO2   
        Reaction(GCARR_ac(1.70e-11, -800.0e0), [Br, CH2O], [HBr, HO2, CO], [1, 1], [1, 1, 1])   #2012/06/07; Parrella2012; JPP  #  Br + CH2O = HBr + HO2 + CO   
        Reaction(GCARR_ac(1.80e-11, -460.0e0), [Br, ALD2], [HBr, MCO3], [1, 1], [1, 1])   #2017/07/27; Parrella2012,Fix C creation; SAS,BHH,MJE  #  Br + ALD2 = HBr + MCO3   
        Reaction(GCARR_ac(1.66e-10, -7000.0e0), [Br, ACET], [HBr, ATO2], [1, 1], [1, 1])   #2017/07/27; Parrella2012,Fix C creation; SAS,BHH,MJE  #  Br + ACET = HBr + ATO2   
        Reaction(GCARR_ac(2.36e-10, -6411.0e0), [Br, C2H6], [HBr, ETO2], [1, 1], [1, 1])   #2017/07/27; Parrella2012,Fix C creation; SAS,BHH,MJE  #  Br + C2H6 = HBr + ETO2   
        Reaction(GCARR_ac(8.77e-11, -4330.0e0), [Br, C3H8], [HBr, A3O2], [1, 1], [1, 1])   #2017/07/27; Parrella2012,Fix C creation; SAS,BHH,MJE  #  Br + C3H8 = HBr + A3O2   
        Reaction(GCARR_ac(9.00e-13, -360.0e0), [CHBr3, OH], [Br], [1, 1], [3.000])   #2017/02/22; JPL 15-10; BHH,MJE  #  CHBr3 + OH = 3.000Br   
        Reaction(GCARR_ac(2.00e-12, -840.0e0), [CH2Br2, OH], [Br], [1, 1], [2.000])   #2012/06/07; Parrella2012; JPP  #  CH2Br2 + OH = 2.000Br   
        Reaction(GCARR_ac(1.42e-12, -1150.0e0), [CH3Br, OH], [Br, H2O, HO2], [1, 1], [1, 1, 1])   #2017/03/08; JPL 15-10; TS,BHH,MJE  #  CH3Br + OH = Br + H2O + HO2   
        Reaction(GCARR_ac(1.63e-10, 60.0e0), [O1D, H2O], [OH], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  O1D + H2O = 2.000OH   
        Reaction(GCARR_ac(2.15e-11, 110.0e0), [O1D, N2], [O, N2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  O1D + N2 = O + N2   
        Reaction(GCARR_ac(3.30e-11, 55.0e0), [O1D, O2], [O, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  O1D + O2 = O + O2   
        Reaction(1.20e-10, [O1D, H2], [H, OH], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  O1D + H2 = H + OH   
        Reaction(GCARR_ac(4.63e-11, 20.0e0), [O1D, N2O], [N2, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  O1D + N2O = N2 + O2   
        Reaction(GCARR_ac(7.25e-11, 20.0e0), [O1D, N2O], [NO], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  O1D + N2O = 2.000NO   
        Reaction(1.31e-10, [O1D, CH4], [MO2, OH], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  O1D + CH4 = MO2 + OH   
        Reaction(0.09e-10, [O1D, CH4], [CH2O, H2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  O1D + CH4 = CH2O + H2   
        Reaction(0.35e-10, [O1D, CH4], [CH2O, H, HO2], [1, 1], [1, 1, 1])   #2014/02/03; Eastham2014; SDE  #  O1D + CH4 = CH2O + H + HO2   
        Reaction(GCARR_ac(8.00e-12, -2060.0e0), [O, O3], [O2], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  O + O3 = 2.000O2   
        Reaction(GCARR_ac(2.80e-12, -1800.0e0), [OH, H2], [H2O, H], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + H2 = H2O + H   
        Reaction(GCARR_ac(1.80e-11, 180.0e0), [O, OH], [O2, H], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  O + OH = O2 + H   
        Reaction(GCARR_ac(3.00e-11, 200.0e0), [HO2, O], [OH, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  HO2 + O = OH + O2   
        Reaction(1.20e-10, [O1D, O3], [O2], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  O1D + O3 = 2.000O2   
        Reaction(1.20e-10, [O1D, O3], [O, O2], [1, 1], [2.000, 1])   #2014/02/03; Eastham2014; SDE  #  O1D + O3 = 2.000O + O2   
        Reaction(GCARR_ac(2.10e-11, -2200.0e0), [OCS, O], [CO, SO2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OCS + O = CO + SO2   
        Reaction(GCARR_ac(1.10e-13, -1200.0e0), [OCS, OH], [CO2, SO2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OCS + OH = CO2 + SO2   
        Reaction(GCARR_ac(5.10e-12, 210.0e0), [NO2, O], [NO, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  NO2 + O = NO + O2   
        Reaction(1.00e-11, [NO3, O], [NO2, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  NO3 + O = NO2 + O2   
        Reaction(GCARR_ac(1.40e-12, -2000.0e0), [H2O2, O], [OH, HO2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  H2O2 + O = OH + HO2   
        Reaction(GCARR_ac(1.40e-10, -470.0e0), [H, O3], [OH, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  H + O3 = OH + O2   
        Reaction(7.20e-11, [H, HO2], [OH], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  H + HO2 = 2.000OH   
        Reaction(1.60e-12, [H, HO2], [O, H2O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  H + HO2 = O + H2O   
        Reaction(6.90e-12, [H, HO2], [H2, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  H + HO2 = H2 + O2   
        Reaction(GCARR_ac(1.50e-11, -3600.0e0), [N, O2], [NO, O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  N + O2 = NO + O   
        Reaction(GCARR_ac(2.10e-11, 100.0e0), [N, NO], [N2, O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  N + NO = N2 + O   
        Reaction(GCARR_ac(5.80e-12, 220.0e0), [N, NO2], [N2O, O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  N + NO2 = N2O + O   
        Reaction(GCARR_ac(1.90e-11, 230.0e0), [BrO, O], [Br, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  BrO + O = Br + O2   
        Reaction(GCARR_ac(3.40e-11, -1600.0e0), [CH2O, O], [CO, HO2, OH], [1, 1], [1, 1, 1])   #2014/02/03; Eastham2014; SDE  #  CH2O + O = CO + HO2 + OH   
        Reaction(1.80e-10, [O1D, CH3Br], [BrO, MO2, Br], [1, 1], [0.440, 1, 0.560])   #2014/02/03; Eastham2014; SDE  #  O1D + CH3Br = 0.440BrO + MO2 + 0.560Br   
        Reaction(GCARR_ac(2.60e-12, -1100.0e0), [OH, Cl2], [HOCl, Cl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + Cl2 = HOCl + Cl   
        Reaction(GCARR_ac(1.80e-11, -600.0e0), [MO2, ClO], [ClOO, HO2, CH2O], [1, 1], [1, 1, 1])   #2017/03/20; JPL 15-10; TS,BHH,MJE  #  MO2 + ClO = ClOO + HO2 + CH2O   
        Reaction(GCARR_ac(7.40e-12, 270.0e0), [OH, ClO], [HO2, Cl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + ClO = HO2 + Cl   
        Reaction(GCARR_ac(6.00e-13, 230.0e0), [OH, ClO], [HCl, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + ClO = HCl + O2   
        Reaction(GCARR_ac(1.40e-12, 600.0e0), [OH, OClO], [HOCl, O2], [1, 1], [1, 1])   #2017/02/22; JPL 15-10; BHH,MJE  #  OH + OClO = HOCl + O2   
        Reaction(GCARR_ac(6.00e-13, 670.0e0), [OH, Cl2O2], [HOCl, ClOO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + Cl2O2 = HOCl + ClOO   
        Reaction(GCARR_ac(1.80e-12, -250.0e0), [OH, HCl], [H2O, Cl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + HCl = H2O + Cl   
        Reaction(GCARR_ac(3.00e-12, -500.0e0), [OH, HOCl], [H2O, ClO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + HOCl = H2O + ClO   
        Reaction(GCARR_ac(2.40e-12, -1250.0e0), [OH, ClNO2], [HOCl, NO2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + ClNO2 = HOCl + NO2   
        Reaction(GCARR_ac(1.20e-12, -330.0e0), [OH, ClNO3], [HOCl, NO3], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + ClNO3 = HOCl + NO3   
        Reaction(GCARR_ac(1.96e-12, -1200.0e0), [OH, CH3Cl], [Cl, HO2, H2O], [1, 1], [1, 1, 1])   #2017/02/22; JPL 15-10; BHH,MJE  #  OH + CH3Cl = Cl + HO2 + H2O   
        Reaction(GCARR_ac(2.61e-12, -944.0e0), [OH, CH2Cl2], [Cl, HO2], [1, 1], [2.000, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  OH + CH2Cl2 = 2.000Cl + HO2   
        Reaction(GCARR_ac(4.69e-12, -1134.0e0), [OH, CHCl3], [Cl, HO2], [1, 1], [3.000, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  OH + CHCl3 = 3.000Cl + HO2   
        Reaction(GCARR_ac(1.64e-12, -1520.0e0), [OH, CH3CCl3], [Cl, H2O], [1, 1], [3.000, 1])   #2014/02/03; Eastham2014; SDE  #  OH + CH3CCl3 = 3.000Cl + H2O   
        Reaction(GCARR_ac(9.20e-13, -1560.0e0), [OH, HCFC22], [Cl, H2O], [1, 1], [1, 1])   #2017/02/22; JPL 15-10; BHH,MJE  #  OH + HCFC22 = Cl + H2O   
        Reaction(GCARR_ac(1.25e-12, -1600.0e0), [OH, HCFC141b], [Cl, H2O], [1, 1], [2.000, 1])   #2014/02/03; Eastham2014; SDE  #  OH + HCFC141b = 2.000Cl + H2O   
        Reaction(GCARR_ac(1.30e-12, -1770.0e0), [OH, HCFC142b], [Cl, H2O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + HCFC142b = Cl + H2O   
        Reaction(GCARR_ac(7.40e-13, -900.0e0), [OH, HCFC123], [Cl, H2O], [1, 1], [2.000, 1])   #2017/02/22; JPL 15-10; BHH,MJE  #  OH + HCFC123 = 2.000Cl + H2O   
        Reaction(GCARR_ac(7.10e-12, -1270.0e0), [CH4, Cl], [HCl, MO2], [1, 1], [1, 1])   #2017/03/08; JPL 15-10; TS,BHH,MJE  #  CH4 + Cl = HCl + MO2   
        Reaction(GCARR_ac(7.32e-11, -30.0e0), [CH2O, Cl], [CO, HCl, HO2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b; TS,JAS,SDE  #  CH2O + Cl = CO + HCl + HO2   
        Reaction(GCARR_ac(2.30e-11, -200.0e0), [Cl, O3], [ClO, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  Cl + O3 = ClO + O2   
        Reaction(GCARR_ac(3.05e-11, -2270.0e0), [Cl, H2], [H, HCl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  Cl + H2 = H + HCl   
        Reaction(GCARR_ac(1.10e-11, -980.0e0), [Cl, H2O2], [HO2, HCl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  Cl + H2O2 = HO2 + HCl   
        Reaction(GCARR_ac(1.40e-11, 270.0e0), [Cl, HO2], [O2, HCl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  Cl + HO2 = O2 + HCl   
        Reaction(GCARR_ac(3.60e-11, -375.0e0), [Cl, HO2], [OH, ClO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  Cl + HO2 = OH + ClO   
        Reaction(GCARR_ac(2.80e-11, 85.0e0), [ClO, O], [Cl, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + O = Cl + O2   
        Reaction(GCARR_ac(2.60e-12, 290.0e0), [ClO, HO2], [O2, HOCl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + HO2 = O2 + HOCl   
        Reaction(GCARR_ac(6.40e-12, 290.0e0), [ClO, NO], [Cl, NO2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + NO = Cl + NO2   
        Reaction(GCARR_ac(1.00e-12, -1590.0e0), [ClO, ClO], [Cl2, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + ClO = Cl2 + O2   
        Reaction(GCARR_ac(3.00e-11, -2450.0e0), [ClO, ClO], [Cl, ClOO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + ClO = Cl + ClOO   
        Reaction(GCARR_ac(3.50e-13, -1370.0e0), [ClO, ClO], [OClO, Cl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + ClO = OClO + Cl   
        Reaction(2.30e-10, [ClOO, Cl], [Cl2, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClOO + Cl = Cl2 + O2   
        Reaction(1.20e-11, [ClOO, Cl], [ClO], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  ClOO + Cl = 2.000ClO   
        Reaction(GCARR_ac(9.50e-13, 550.0e0), [ClO, BrO], [Br, OClO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + BrO = Br + OClO   
        Reaction(GCARR_ac(2.30e-12, 260.0e0), [ClO, BrO], [Br, ClOO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + BrO = Br + ClOO   
        Reaction(GCARR_ac(4.10e-13, 290.0e0), [ClO, BrO], [BrCl, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + BrO = BrCl + O2   
        Reaction(GCARR_ac(3.60e-12, -840.0e0), [ClNO3, O], [ClO, NO3], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClNO3 + O = ClO + NO3   
        Reaction(GCARR_ac(6.50e-12, 135.0e0), [ClNO3, Cl], [Cl2, NO3], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClNO3 + Cl = Cl2 + NO3   
        Reaction(GCARR_ac(2.17e-11, -1130.0e0), [CH3Cl, Cl], [CO, HCl, HO2], [1, 1], [1, 2.000, 1])   #2014/02/03; Eastham2014; SDE  #  CH3Cl + Cl = CO + 2.000HCl + HO2   
        Reaction(GCARR_ac(1.24e-12, -1070.0e0), [CH2Cl2, Cl], [CO, HCl, Cl, HO2], [1, 1], [1, 1, 2.000, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  CH2Cl2 + Cl = CO + HCl + 2.000Cl + HO2   
        Reaction(GCARR_ac(3.77e-12, -1011.0e0), [CHCl3, Cl], [CO, HCl, Cl, HO2], [1, 1], [1, 1, 3.000, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  CHCl3 + Cl = CO + HCl + 3.000Cl + HO2   
        Reaction(2.00e-13, [Cl, HCOOH], [HCl, CO2, H2O], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + HCOOH = HCl + CO2 + H2O   
        Reaction(1.60e-10, [Cl, MO2], [ClO, CH2O, HO2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + MO2 = ClO + CH2O + HO2   
        Reaction(5.7e-11, [Cl, MP], [HCl, MO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + MP = HCl + MO2   
        Reaction(GCARR_ac(7.2e-11, -70.0e0), [Cl, C2H6], [HCl, ETO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + C2H6 = HCl + ETO2   
        Reaction(7.4e-11, [Cl, ETO2], [ClO, HO2, ALD2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + ETO2 = ClO + HO2 + ALD2   
        Reaction(7.4e-11, [Cl, OTHRO2], [ClO, HO2, ALD2], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF  #  Cl + OTHRO2 = ClO + HO2 + ALD2   
        Reaction(5.5e-11, [Cl, MOH], [HCl, CH2O, HO2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + MOH = HCl + CH2O + HO2   
        Reaction(9.6e-11, [Cl, EOH], [HCl, ALD2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + EOH = HCl + ALD2   
        Reaction(2.8e-14, [Cl, ACTA], [HCl, MO2, CO2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + ACTA = HCl + MO2 + CO2   
        Reaction(GCARR_ac(6.54e-11, 60.0e0), [Cl, C3H8], [HCl, B3O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + C3H8 = HCl + B3O2   
        Reaction(GCARR_ac(8.12e-11, -90.0e0), [Cl, C3H8], [HCl, A3O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + C3H8 = HCl + A3O2   
        Reaction(GCARR_ac(7.70e-11, -1000.0e0), [Cl, ACET], [HCl, ATO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + ACET = HCl + ATO2   
        Reaction(GCARR_ac(7.60e-11, 500.0e0), [Cl, ISOP], [HCl, IHOO1, IHOO4], [1, 1], [1, 0.5, 0.5])   #2019/11/06; Sherwen2016b;KHB,TS,JAS,SDE  #  Cl + ISOP = HCl + 0.5IHOO1 + 0.5IHOO4   
        Reaction(2.05e-10, [Cl, ALK4], [HCl, R4O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + ALK4 = HCl + R4O2   
        Reaction(3.60e-12, [Br, PRPE], [HBr, PO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Br + PRPE = HBr + PO2   
        Reaction(GCARR_ac(8.40e-11, -2620.0e0), [INO, INO], [I2, NO], [1, 1], [1, 2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  INO + INO = I2 + 2.000NO   
        Reaction(GCARR_ac(2.90e-11, -2600.0e0), [IONO, IONO], [I2, NO2], [1, 1], [1, 2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IONO + IONO = I2 + 2.000NO2   
        Reaction(1.50e-12, [I2, NO3], [I, IONO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2 + NO3 = I + IONO2   
        Reaction(GCARR_ac(9.10e-11, -146.0e0), [IONO2, I], [I2, NO3], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IONO2 + I = I2 + NO3   
        Reaction(1.20e-11, [I, BrO], [IO, Br], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I + BrO = IO + Br   
        Reaction(GCARR_ac(3.00e-12, 510.0e0), [IO, BrO], [Br, I, O2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + BrO = Br + I + O2   
        Reaction(GCARR_ac(1.20e-11, 510.0e0), [IO, BrO], [Br, OIO], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + BrO = Br + OIO   
        Reaction(1.50e-10, [OIO, OIO], [I2O4], [1, 1], [1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  OIO + OIO = I2O4   
        Reaction(GCARR_ac(1.10e-12, 542.0e0), [OIO, NO], [IO, NO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  OIO + NO = IO + NO2   
        Reaction(GCARR_ac(5.10e-12, 280.0e0), [IO, ClO], [I, OClO], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + ClO = I + OClO   
        Reaction(GCARR_ac(2.81e-12, 280.0e0), [IO, ClO], [I, Cl, O2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + ClO = I + Cl + O2   
        Reaction(GCARR_ac(1.02e-12, 280.0e0), [IO, ClO], [ICl, O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + ClO = ICl + O2   
        Reaction(GCARR_ac(2.30e-11, -870.0e0), [I, O3], [IO, O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2017;TS,JAS,SDE  #  I + O3 = IO + O2   
        Reaction(GCARR_ac(1.50e-11, -1090.0e0), [I, HO2], [HI, O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I + HO2 = HI + O2   
        Reaction(1.80e-10, [I2, OH], [HOI, I], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2 + OH = HOI + I   
        Reaction(3.00e-11, [HI, OH], [I, H2O], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  HI + OH = I + H2O   
        Reaction(5.00e-12, [HOI, OH], [IO, H2O], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  HOI + OH = IO + H2O   
        Reaction(GCARR_ac(1.30e-11, 570.0e0), [IO, HO2], [HOI, O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + HO2 = HOI + O2   
        Reaction(GCARR_ac(9.10e-12, 240.0e0), [IO, NO], [I, NO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + NO = I + NO2   
        Reaction(GCARR_ac(6.00e-12, 500.0e0), [IO, IO], [I, OIO], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + IO = I + OIO   
        Reaction(GCARR_ac(2.90e-12, -1100.0e0), [CH3I, OH], [H2O, I, MO2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  CH3I + OH = H2O + I + MO2   
        Reaction(2.40e-12, [ETHLN, OH], [CH2O, CO2, NO2], [1, 1], [1, 1, 1])   #2017/06/15, Marais2016, EAM  #  ETHLN + OH = CH2O + CO2 + NO2   
        Reaction(6.70e-13, [PROPNN, OH], [NO2, MGLY], [1, 1], [1, 1])   #2017/07/14; MCMv3.3; KRT,JAF,CCM,EAM,KHB,RHS  #  PROPNN + OH = NO2 + MGLY   
        Reaction(1.20e-15, [CH2OO, CO], [CH2O], [1, 1], [1])   #2015/09/25; Millet2015; DBM,EAM  #  CH2OO + CO = CH2O   
        Reaction(1.00e-14, [CH2OO, NO], [CH2O, NO2], [1, 1], [1, 1])   #2015/09/25; Millet2015; DBM,EAM  #  CH2OO + NO = CH2O + NO2   
        Reaction(1.00e-15, [CH2OO, NO2], [CH2O, NO3], [1, 1], [1, 1])   #2015/09/25; Millet2015; DBM,EAM  #  CH2OO + NO2 = CH2O + NO3   
        Reaction(1.40e-12, [CH2OO, O3], [CH2O], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  CH2OO + O3 = CH2O   
        Reaction(3.70e-11, [CH2OO, SO2], [CH2O, SO4], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  CH2OO + SO2 = CH2O + SO4   
        Reaction(1.20e-15, [CH3CHOO, CO], [ALD2], [1, 1], [1])   #2015/09/25; Millet2015; DBM,EAM  #  CH3CHOO + CO = ALD2   
        Reaction(1.00e-14, [CH3CHOO, NO], [ALD2, NO2], [1, 1], [1, 1])   #2015/09/25; Millet2015; DBM,EAM  #  CH3CHOO + NO = ALD2 + NO2   
        Reaction(1.00e-15, [CH3CHOO, NO2], [ALD2, NO3], [1, 1], [1, 1])   #2015/09/25; Millet2015; DBM,EAM  #  CH3CHOO + NO2 = ALD2 + NO3   
        Reaction(7.00e-14, [CH3CHOO, SO2], [ALD2, SO4], [1, 1], [1, 1])   #2015/09/25; Millet2015; DBM,EAM  #  CH3CHOO + SO2 = ALD2 + SO4   
        Reaction(6.00e-18, [CH3CHOO, H2O], [ALD2, H2O2], [1, 1], [1, 1])   #2015/09/25; Millet2015; DBM,EAM  #  CH3CHOO + H2O = ALD2 + H2O2   
        Reaction(1.00e-17, [CH3CHOO, H2O], [ACTA], [1, 1], [1])   #2015/09/25; Millet2015; DBM,EAM  #  CH3CHOO + H2O = ACTA   
        Reaction(GCARR_ac(1.21e-11, 440.0e0), [MTPA, OH], [PIO2], [1, 1], [1])   #2017/03/23; IUPAC2010; EVF  #  MTPA + OH = PIO2   
        Reaction(GCARR_ac(1.21e-11, 440.0e0), [MTPO, OH], [PIO2], [1, 1], [1])   #2017/03/23; IUPAC2010; EVF  #  MTPO + OH = PIO2   
        Reaction(1.50e-11, [PIO2, HO2], [PIP], [1, 1], [1])   #2017/03/23; Roberts1992; EVF  #  PIO2 + HO2 = PIP   
        Reaction(1.20e-12, [PIO2, NO3], [HO2, NO2, RCHO, MEK], [1, 1], [1, 1, 1, 1])   #2017/03/23; Roberts1992; EVF  #  PIO2 + NO3 = HO2 + NO2 + RCHO + MEK   
        Reaction(GCARR_ac(8.33e-13, 490.0e0), [MTPA, NO3], [OLNN, OLND], [1, 1], [0.100, 0.900])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS  #  MTPA + NO3 = 0.100OLNN + 0.900OLND   
        Reaction(GCARR_ac(8.33e-13, 490.0e0), [MTPO, NO3], [OLNN, OLND], [1, 1], [0.100, 0.900])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS  #  MTPO + NO3 = 0.100OLNN + 0.900OLND   
        Reaction(GCARR_ac(4.20e-11, 401.0e0), [LIMO, OH], [LIMO2], [1, 1], [1])   #2017/07/14; Gill2002; KRT,JAF,CCM,EAM,KHB,RHS  #  LIMO + OH = LIMO2   
        Reaction(1.22e-11, [LIMO, NO3], [OLNN, OLND], [1, 1], [0.500, 0.500])   #2017/07/14; Fry2014,Atkinson2003; KRT,JAF,CCM,EAM,KHB,RHS  #  LIMO + NO3 = 0.500OLNN + 0.500OLND   
        Reaction(1.50e-11, [LIMO2, HO2], [PIP], [1, 1], [1])   #2017/07/14; Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS  #  LIMO2 + HO2 = PIP   
        Reaction(4.00e-12, [OLNN, NO], [HO2, NO2, MONITS], [1, 1], [1, 1, 1])   #2017/07/14; Browne2014,Goliff2013; KRT,JAF,CCM,EAM,KHB,RHS  #  OLNN + NO = HO2 + NO2 + MONITS   
        Reaction(GCARR_ac(1.66e-13, 1300.0e0), [OLNN, HO2], [MONITS, MONITU], [1, 1], [0.700, 0.300])   #2017/07/14; Browne2014,Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS  #  OLNN + HO2 = 0.700MONITS + 0.300MONITU   
        Reaction(GCARR_ac(1.66e-13, 1300.0e0), [OLND, HO2], [MONITS, MONITU], [1, 1], [0.700, 0.300])   #2017/07/14; Browne2014,Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS  #  OLND + HO2 = 0.700MONITS + 0.300MONITU   
        Reaction(4.80e-12, [MONITS, OH], [HONIT], [1, 1], [1])   #2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS  #  MONITS + OH = HONIT   
        Reaction(7.29e-11, [MONITU, OH], [HONIT], [1, 1], [1])   #2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS  #  MONITU + OH = HONIT   
        Reaction(1.67e-16, [MONITU, O3], [HONIT], [1, 1], [1])   #2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS  #  MONITU + O3 = HONIT   
        Reaction(GCARR_ac(3.15e-13, -448.0e0), [MONITU, NO3], [HONIT], [1, 1], [1])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS  #  MONITU + NO3 = HONIT   
        Reaction(GCARR_ac(3.15e-13, -448.0e0), [MONITS, NO3], [HONIT], [1, 1], [1])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS  #  MONITS + NO3 = HONIT   
        Reaction(2.78e-04, [IONITA], [INDIOL, HNO3], [1], [1, 1])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS  #  IONITA = INDIOL + HNO3   
        Reaction(2.78e-04, [MONITA], [INDIOL, HNO3], [1], [1, 1])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS  #  MONITA = INDIOL + HNO3   
        Reaction(GC_OHHNO3_acacac(2.41e-14, 460.0e0, 2.69e-17, 2199.0e0, 6.51e-34, 1335.0e0), [HONIT, OH], [NO3, HAC], [1, 1], [1, 1])   #2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS  #  HONIT + OH = NO3 + HAC   
        Reaction(GCARR_ac(8.00e-13, -1000.0e0), [MENO3, OH], [CH2O, NO2], [1, 1], [1, 1])   #2019/05/16; JPL 15-10,Fisher2018; JAF  #  MENO3 + OH = CH2O + NO2   
        Reaction(GCARR_ac(1.00e-12, -490.0e0), [ETNO3, OH], [ALD2, NO2], [1, 1], [1, 1])   #2019/05/16; JPL 15-10,Fisher2018; JAF  #  ETNO3 + OH = ALD2 + NO2   
        Reaction(GCARR_ac(1.20e-12, -320.0e0), [IPRNO3, OH], [ACET, NO2], [1, 1], [1, 1])   #2019/05/16; JPL 15-10,Fisher2018; JAF  #  IPRNO3 + OH = ACET + NO2   
        Reaction(7.10e-13, [NPRNO3, OH], [RCHO, NO2], [1, 1], [1, 1])   #2019/05/16; JPL 15-10,Fisher2018; JAF  #  NPRNO3 + OH = RCHO + NO2   
        Reaction(GC_ISO1(1.7e-11, 3.90e2, 9.33e-2, 5.05e15, -1.22e4, 1.79e14, -8.830e3), [ISOP, OH], [LISOPOH, IHOO1], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  ISOP + OH = LISOPOH + IHOO1   
        Reaction(GC_ISO1(1.0e-11, 3.90e2, 2.26e-1, 2.22e9, -7.160e3, 1.75e14, -9.054e3), [ISOP, OH], [LISOPOH, IHOO4], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  ISOP + OH = LISOPOH + IHOO4   
        Reaction(ARRPLUS_abde(2.12e-13, -1300e0, -0.1644e0, 7.0485e-4), [IHOO1, HO2], [RIPC], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHOO1 + HO2 = RIPC   
        Reaction(ARRPLUS_abde(2.12e-13, -1300e0, -0.2038e0, 9.0435e-4), [IHOO4, HO2], [RIPD], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHOO4 + HO2 = RIPD   
        Reaction(ARRPLUS_abde(1.04e11, 9.746e3,  1.1644e0, -7.0485e-4), [IHOO1], [CH2O, OH, MVK], [1], [1, 1, 1])   #2019/11/06; Bates2019; KHB  #  IHOO1 = CH2O + OH + MVK   
        Reaction(ARRPLUS_abde(1.88e11, 9.752e3, 1.2038e0, -9.0435e-4), [IHOO4], [MACR, OH, CH2O], [1], [1, 1, 1])   #2019/11/06; Bates2019; KHB  #  IHOO4 = MACR + OH + CH2O   
        Reaction(ARRPLUS_ade(6.92e-14, 1.1644e0, -7.0485e-4), [IHOO1, IHOO1], [MVK, HO2, CH2O], [1, 1], [2, 2, 2])   #2019/11/06; Bates2019; KHB  #  IHOO1 + IHOO1 = 2MVK + 2HO2 + 2CH2O   
        Reaction(ARRPLUS_ade(5.74e-12, 1.2038e0, -9.0435e-4), [IHOO4, IHOO4], [MACR, HO2, CH2O], [1, 1], [2, 2, 2])   #2019/11/06; Bates2019; KHB  #  IHOO4 + IHOO4 = 2MACR + 2HO2 + 2CH2O   
        Reaction(ARRPLUS_ade(1.54e-12, 2.3682e0, -1.6092e-3), [IHOO1, IHOO4], [MACR, MVK, HO2, CH2O], [1, 1], [1, 1, 2, 2])   #2019/11/06; Bates2019; KHB  #  IHOO1 + IHOO4 = MACR + MVK + 2HO2 + 2CH2O   
        Reaction(ARRPLUS_ade(2.0e-12, 1.1644e0, -7.0485e-4), [IHOO1, MO2], [MVK, HO2, CH2O], [1, 1], [1, 2, 2])   #2019/11/06; Bates2019; KHB  #  IHOO1 + MO2 = MVK + 2HO2 + 2CH2O   
        Reaction(ARRPLUS_ade(2.0e-12, 1.2038e0, -9.0435e-4), [IHOO4, MO2], [MACR, HO2, CH2O], [1, 1], [1, 2, 2])   #2019/11/06; Bates2019; KHB  #  IHOO4 + MO2 = MACR + 2HO2 + 2CH2O   
        Reaction(GC_NIT(2.7e-12, 3.50e2, 1.19e0,  6.0e0, 1.1644e0, 7.05e-4), [IHOO1, NO], [IHN2], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHOO1 + NO = IHN2   
        Reaction(GC_ALK(2.7e-12, 3.50e2, 1.19e0,  6.0e0, 1.1644e0, 7.05e-4), [IHOO1, NO], [NO2, MVK, HO2, CH2O], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  IHOO1 + NO = NO2 + MVK + HO2 + CH2O   
        Reaction(GC_NIT(2.7e-12, 3.50e2, 1.421e0, 6.0e0, -0.1644e0, -7.05e-4), [IHOO1, NO], [IHN4], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHOO1 + NO = IHN4   
        Reaction(GC_NIT(2.7e-12, 3.50e2, 1.297e0, 6.0e0, 1.2038e0, 9.04e-4), [IHOO4, NO], [IHN3], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHOO4 + NO = IHN3   
        Reaction(GC_ALK(2.7e-12, 3.50e2, 1.297e0, 6.0e0, 1.2038e0, 9.04e-4), [IHOO4, NO], [NO2, MACR, HO2, CH2O], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  IHOO4 + NO = NO2 + MACR + HO2 + CH2O   
        Reaction(GC_NIT(2.7e-12, 3.50e2, 1.421e0, 6.0e0, -0.2038e0, -9.04e-4), [IHOO4, NO], [IHN1], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHOO4 + NO = IHN1   
        Reaction(GCARR_ac(3.00e-12, 650.0e0), [IDC, OH], [CO, HO2, MVKPC], [1, 1], [1, 1, 1])   #2019/11/06; Bates2019; KHB  #  IDC + OH = CO + HO2 + MVKPC   
        Reaction(GCARR_ac(1.59e+13, -10000.0e0), [IHPOO1], [ICPDH, IDHPE, OH], [1], [0.176, 0.824, 1])   #2019/11/06; Bates2019; KHB  #  IHPOO1 = 0.176ICPDH + 0.824IDHPE + OH   
        Reaction(GC_NIT(2.7e-12, 3.50e2, 2.1e0, 9.0e0, 1.0e0, 0.0e0), [IHPOO1, NO], [ITHN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHPOO1 + NO = ITHN   
        Reaction(GCARR_ac(2.91e+13, -10000.0e0), [IHPOO2], [ICPDH, IDHPE, OH], [1], [0.548, 0.452, 1])   #2019/11/06; Bates2019; KHB  #  IHPOO2 = 0.548ICPDH + 0.452IDHPE + OH   
        Reaction(GC_NIT(2.7e-12, 3.50e2, 2.315e0, 9.0e0, 1.0e0, 0.0e0), [IHPOO2, NO], [ITHN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHPOO2 + NO = ITHN   
        Reaction(GCARR_ac(1.875e+13, -10000.0e0), [IHPOO3], [IDHPE], [1], [1])   #2019/11/06; Bates2019; KHB  #  IHPOO3 = IDHPE   
        Reaction(GC_ALK(2.7e-12, 3.50e2, 3.079e0, 9.0e0, 1.0e0, 0.0e0), [IHPOO3, NO], [GLYC, HAC, NO2, OH], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  IHPOO3 + NO = GLYC + HAC + NO2 + OH   
        Reaction(GC_NIT(2.7e-12, 3.50e2, 3.079e0, 9.0e0, 1.0e0, 0.0e0), [IHPOO3, NO], [ITHN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHPOO3 + NO = ITHN   
        Reaction(GCARR_ac(1.05e-11, -400.0e0), [IEPOXA, OH], [ICHE, HO2], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  IEPOXA + OH = ICHE + HO2   
        Reaction(GC_EPO_a(5.82e-11, -4.00e2, 1.14e-20), [IEPOXA, OH], [IEPOXAOO, IEPOXBOO], [1, 1], [0.67, 0.33])   #2019/11/06; Bates2019; KHB  #  IEPOXA + OH = 0.67IEPOXAOO + 0.33IEPOXBOO   
        Reaction(GCARR_ac(8.25e-12, -400.0e0), [IEPOXB, OH], [ICHE, HO2], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  IEPOXB + OH = ICHE + HO2   
        Reaction(GC_EPO_a(3.75e-11, -4.00e2, 8.91e-21), [IEPOXB, OH], [IEPOXAOO, IEPOXBOO], [1, 1], [0.81, 0.19])   #2019/11/06; Bates2019; KHB  #  IEPOXB + OH = 0.81IEPOXAOO + 0.19IEPOXBOO   
        Reaction(GCARR_ac(1.875e+13, -10000.0e0), [IEPOXAOO], [IDCHP, HO2], [1], [1, 1])   #2019/11/06; Bates2019; KHB  #  IEPOXAOO = IDCHP + HO2   
        Reaction(GCARR_ac(1.0e+7, -5000.0e0), [IEPOXAOO], [OH, CO, MVKDH], [1], [1, 1, 1])   #2019/11/06; Bates2019; KHB  #  IEPOXAOO = OH + CO + MVKDH   
        Reaction(GC_NIT(2.7e-12, 3.50e2, 13.098e0, 8.0e0, 1.0e0, 0.0e0), [IEPOXAOO, NO], [ITCN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IEPOXAOO + NO = ITCN   
        Reaction(GCARR_ac(1.875e+13, -10000.0e0), [IEPOXBOO], [IDCHP, HO2], [1], [1, 1])   #2019/11/06; Bates2019; KHB  #  IEPOXBOO = IDCHP + HO2   
        Reaction(GCARR_ac(1.0e+7, -5000.0e0), [IEPOXBOO], [CO, OH, MCRDH], [1], [1, 1, 1])   #2019/11/06; Bates2019; KHB  #  IEPOXBOO = CO + OH + MCRDH   
        Reaction(GC_NIT(2.7e-12, 3.50e2, 16.463e0, 8.0e0, 1.0e0, 0.0e0), [IEPOXBOO, NO], [ITCN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IEPOXBOO + NO = ITCN   
        Reaction(GC_NIT(2.7e-12, 3.50e2, 13.098e0, 8.0e0, 1.0e0, 0.0e0), [ICHOO, NO], [ITCN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  ICHOO + NO = ITCN   
        Reaction(GCARR_ac(1.875e+13, -10000.0e0), [ICHOO], [HO2, CO, HAC, OH], [1], [1, 2.000, 1, 1])   #2019/11/06; Bates2019; KHB  #  ICHOO = HO2 + 2.000CO + HAC + OH   
        Reaction(GCARR_ac(2.70e-12, 350.0e0), [HPALD1OO, NO], [NO2, OH, CO2, MVK], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  HPALD1OO + NO = NO2 + OH + CO2 + MVK   
        Reaction(GCARR_ac(2.38e-13, 1300.0e0), [HPALD1OO, HO2], [OH, OH, CO2, MVK], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  HPALD1OO + HO2 = OH + OH + CO2 + MVK   
        Reaction(GCARR_ac(2.70e-12, 350.0e0), [HPALD2OO, NO], [NO2, OH, CO2, MACR], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  HPALD2OO + NO = NO2 + OH + CO2 + MACR   
        Reaction(GCARR_ac(2.38e-13, 1300.0e0), [HPALD2OO, HO2], [OH, OH, CO2, MACR], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  HPALD2OO + HO2 = OH + OH + CO2 + MACR   
        Reaction(GCARR_ac(7.14e-12, 390.0e0), [IHN2, OH], [ISOPNOO1], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHN2 + OH = ISOPNOO1   
        Reaction(GC_EPO_a(6.30e-12, 390.0e0, 1.62e-19), [IHN2, OH], [IEPOXA, IEPOXB, NO2], [1, 1], [0.67, 0.33, 1])   #2019/11/06; Bates2019; KHB  #  IHN2 + OH = 0.67IEPOXA + 0.33IEPOXB + NO2   
        Reaction(GCARR_ac(1.02e-11, 390.0e0), [IHN3, OH], [ISOPNOO2], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHN3 + OH = ISOPNOO2   
        Reaction(GC_EPO_a(1.05e-11, 390.0e0, 2.49e-19), [IHN3, OH], [IEPOXA, IEPOXB, NO2], [1, 1], [0.67, 0.33, 1])   #2019/11/06; Bates2019; KHB  #  IHN3 + OH = 0.67IEPOXA + 0.33IEPOXB + NO2   
        Reaction(GC_EPO_a(1.55e-11, 390.0e0, 2.715e-19), [IHN1, OH], [IEPOXD, NO2], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  IHN1 + OH = IEPOXD + NO2   
        Reaction(GCARR_ac(2.04e-11, 390.0e0), [IHN1, OH], [IDHNDOO1], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHN1 + OH = IDHNDOO1   
        Reaction(GC_EPO_a(9.52e-12, 390.0e0, 2.715e-19), [IHN4, OH], [IEPOXD, NO2], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  IHN4 + OH = IEPOXD + NO2   
        Reaction(GCARR_ac(2.95e-11, 390.0e0), [IHN4, OH], [IDHNDOO2], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHN4 + OH = IDHNDOO2   
        Reaction(GCARR_ac(1.875e+13, -10000.0e0), [ISOPNOO1], [ITCN, HO2], [1], [1, 1])   #2019/11/06; Bates2019; KHB  #  ISOPNOO1 = ITCN + HO2   
        Reaction(GC_NIT(2.7e-12, 350.0e0, 6.32e0, 11.0e0, 1.0e0, 0.0e0), [ISOPNOO1, NO], [IDN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  ISOPNOO1 + NO = IDN   
        Reaction(GCARR_ac(1.875e+13, -10000.0e0), [ISOPNOO2], [ITCN, HO2], [1], [1, 1])   #2019/11/06; Bates2019; KHB  #  ISOPNOO2 = ITCN + HO2   
        Reaction(GC_ALK(2.7e-12, 350.0e0, 7.941e0, 11.0e0, 1.0e0, 0.0e0), [ISOPNOO2, NO], [MVKN, CH2O, HO2, NO2], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  ISOPNOO2 + NO = MVKN + CH2O + HO2 + NO2   
        Reaction(GC_NIT(2.7e-12, 350.0e0, 7.941e0, 11.0e0, 1.0e0, 0.0e0), [ISOPNOO2, NO], [IDN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  ISOPNOO2 + NO = IDN   
        Reaction(GCARR_ac(1.256e+13, -10000.0e0), [IDHNDOO1], [ITCN, HO2], [1], [1, 1])   #2019/11/06; Bates2019; KHB  #  IDHNDOO1 = ITCN + HO2   
        Reaction(GCARR_ac(5.092e+12, -10000.0e0), [IDHNDOO2], [ITCN, HO2], [1], [1, 1])   #2019/11/06; Bates2019; KHB  #  IDHNDOO2 = ITCN + HO2   
        Reaction(GC_NIT(2.7e-12, 350.0e0, 4.712e0, 11.0e0, 1.0e0, 0.0e0), [IDHNDOO1, NO], [IDN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IDHNDOO1 + NO = IDN   
        Reaction(GC_NIT(2.7e-12, 350.0e0, 2.258e0, 11.0e0, 1.0e0, 0.0e0), [IDHNDOO2, NO], [IDN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IDHNDOO2 + NO = IDN   


        Reaction(HO2uptk1stOrd( State_Het ), [HO2], [H2O], [1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  HO2 = H2O   
        Reaction(NO2uptk1stOrdAndCloud( State_Het ), [NO2], [HNO3, HNO2], [1], [0.500, 0.500])  #  NO2 = 0.500HNO3 + 0.500HNO2   
        Reaction(NO3uptk1stOrdAndCloud( State_Het ), [NO3], [HNO3], [1], [1])  #  NO3 = HNO3   
        Reaction(NO3hypsisClonSALA( State_Het ), [NO3], [NIT], [1], [1])   #2018/03/16; XW  #  NO3 = NIT   
        Reaction(NO3hypsisClonSALC( State_Het ), [NO3], [NITs], [1], [1])   #2018/03/16; XW  #  NO3 = NITs   
        Reaction(N2O5uptkByH2O( State_Het ), [N2O5, H2O], [HNO3], [1, 1], [2.000])  #  N2O5 + H2O = 2.000HNO3   
        Reaction(N2O5uptkByStratHCl( State_Het ), [N2O5, HCl], [ClNO2, HNO3], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  N2O5 + HCl = ClNO2 + HNO3   
        Reaction(N2O5uptkByCloud( State_Het ), [N2O5], [HNO3], [1], [2.000])   #2018/10/17; Cloud uptake, CDH  #  N2O5 = 2.000HNO3   
        Reaction(N2O5uptkBySALACl( State_Het ), [N2O5, SALACL], [ClNO2, HNO3], [1, 1], [1, 1])   #2018/01/19; Sherwen2017;TS,JAS,SDE,XW  #  N2O5 + SALACL = ClNO2 + HNO3   
        Reaction(N2O5uptkBySALCCl( State_Het ), [N2O5, SALCCL], [ClNO2, HNO3], [1, 1], [1, 1])   #2018/01/19; Sherwen2017;TS,JAS,SDE,XW  #  N2O5 + SALCCL = ClNO2 + HNO3   
        Reaction(OHuptkBySALACl( State_Het ), [OH, SALACL], [Cl2], [1, 1], [0.500])   #2018/03/12; XW  #  OH + SALACL = 0.500Cl2   
        Reaction(OHuptkBySALCCl( State_Het ), [OH, SALCCL], [Cl2], [1, 1], [0.500])   #2018/03/12; XW  #  OH + SALCCL = 0.500Cl2   
        Reaction(BrNO3uptkByH2O( State_Het ), [BrNO3, H2O], [HOBr, HNO3], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  BrNO3 + H2O = HOBr + HNO3   
        Reaction(BrNO3uptkByHCl( State_Het ), [BrNO3, HCl], [BrCl, HNO3], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  BrNO3 + HCl = BrCl + HNO3   
        Reaction(ClNO3uptkByH2O( State_Het ), [ClNO3, H2O], [HOCl, HNO3], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClNO3 + H2O = HOCl + HNO3   
        Reaction(ClNO3uptkByHCl( State_Het ), [ClNO3, HCl], [Cl2, HNO3], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClNO3 + HCl = Cl2 + HNO3   
        Reaction(ClNO3uptkByHBr( State_Het ), [ClNO3, HBr], [BrCl, HNO3], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClNO3 + HBr = BrCl + HNO3   
        Reaction(ClNO3uptkByBrSALA( State_Het ), [ClNO3, BrSALA], [BrCl, HNO3], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  ClNO3 + BrSALA = BrCl + HNO3   
        Reaction(ClNO3uptkByBrSALC( State_Het ), [ClNO3, BrSALC], [BrCl, HNO3], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  ClNO3 + BrSALC = BrCl + HNO3   
        Reaction(ClNO3uptkBySALACL( State_Het ), [ClNO3, SALACL], [Cl2, HNO3], [1, 1], [1, 1])   #2018/01/22; XW  #  ClNO3 + SALACL = Cl2 + HNO3   
        Reaction(ClNO3uptkBySALCCL( State_Het ), [ClNO3, SALCCL], [Cl2, HNO3], [1, 1], [1, 1])   #2018/01/22; XW  #  ClNO3 + SALCCL = Cl2 + HNO3   
        Reaction(ClNO2uptkBySALACL( State_Het ), [ClNO2, SALACL], [Cl2, HNO2], [1, 1], [1, 1])   #2018/01/22; XW  #  ClNO2 + SALACL = Cl2 + HNO2   
        Reaction(ClNO2uptkBySALCCL( State_Het ), [ClNO2, SALCCL], [Cl2, HNO2], [1, 1], [1, 1])   #2018/01/22; XW  #  ClNO2 + SALCCL = Cl2 + HNO2   
        Reaction(ClNO2uptkByHCl( State_Het ), [ClNO2, HCl], [Cl2, HNO2], [1, 1], [1, 1])   #2018/01/22; XW  #  ClNO2 + HCl = Cl2 + HNO2   
        Reaction(ClNO2uptkByBrSALA( State_Het ), [ClNO2, BrSALA], [BrCl, HNO2], [1, 1], [1, 1])   #2018/01/22; XW  #  ClNO2 + BrSALA = BrCl + HNO2   
        Reaction(ClNO2uptkByBrSALC( State_Het ), [ClNO2, BrSALC], [BrCl, HNO2], [1, 1], [1, 1])   #2018/01/22; XW  #  ClNO2 + BrSALC = BrCl + HNO2   
        Reaction(ClNO2uptkByHBr( State_Het ), [ClNO2, HBr], [BrCl, HNO2], [1, 1], [1, 1])   #2018/01/22; XW  #  ClNO2 + HBr = BrCl + HNO2   
        Reaction(HOClUptkByHCl( State_Het ), [HOCl, HCl], [Cl2, H2O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  HOCl + HCl = Cl2 + H2O   
        Reaction(HOClUptkByHBr( State_Het ), [HOCl, HBr], [BrCl, H2O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  HOCl + HBr = BrCl + H2O   
        Reaction(HOClUptkBySALACL( State_Het ), [HOCl, SALACL], [Cl2, H2O], [1, 1], [1, 1])   #2018/01/22; XW  #  HOCl + SALACL = Cl2 + H2O   
        Reaction(HOClUptkBySALCCL( State_Het ), [HOCl, SALCCL], [Cl2, H2O], [1, 1], [1, 1])   #2018/01/22; XW  #  HOCl + SALCCL = Cl2 + H2O   
        Reaction(HOClUptkByHSO3m( State_Het ) + HOClUptkBySO3mm( State_Het ), [HOCl, SO2], [SO4s, HCl], [1, 1], [1, 1])   #2018/11/08; XW; June 6, 2021, MSL  #  HOCl + SO2 = SO4s + HCl   
        Reaction(HOBrUptkByHBr( State_Het ), [HOBr, HBr], [Br2, H2O], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  HOBr + HBr = Br2 + H2O   
        Reaction(HOBrUptkByHCl( State_Het ), [HOBr, HCl], [BrCl, H2O], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  HOBr + HCl = BrCl + H2O   
        Reaction(HOBrUptkBySALACL( State_Het ), [HOBr, SALACL], [BrCl, H2O], [1, 1], [1, 1])   #2018/01/22; Sherwen2017;TS,JAS,SDE;XW  #  HOBr + SALACL = BrCl + H2O   
        Reaction(HOBrUptkBySALCCL( State_Het ), [HOBr, SALCCL], [BrCl, H2O], [1, 1], [1, 1])   #2018/01/22; Sherwen2017;TS,JAS,SDE,XW  #  HOBr + SALCCL = BrCl + H2O   
        Reaction(HOBrUptkByBrSALA( State_Het ), [HOBr, BrSALA], [Br2], [1, 1], [1])   #2017/09/22; Sherwen2017;TS,JAS,SDE  #  HOBr + BrSALA = Br2   
        Reaction(HOBrUptkByBrSALC( State_Het ), [HOBr, BrSALC], [Br2], [1, 1], [1])   #2017/09/22; Sherwen2017;TS,JAS,SDE  #  HOBr + BrSALC = Br2   
        Reaction(HOBrUptkByHSO3m( State_Het ) + HOBrUptkBySO3mm( State_Het ), [HOBr, SO2], [SO4s, HBr], [1, 1], [1, 1])   #2017/11/15; Chen2017; QJC; June 6, 2021, MSL  #  HOBr + SO2 = SO4s + HBr   
        Reaction(O3uptkByHBr( State_Het ), [O3, HBr], [HOBr], [1, 1], [1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  O3 + HBr = HOBr   
        Reaction(O3uptkByBrSALA( State_Het ), [O3, BrSALA], [HOBr], [1, 1], [1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  O3 + BrSALA = HOBr   
        Reaction(O3uptkByBrSALC( State_Het ), [O3, BrSALC], [HOBr], [1, 1], [1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  O3 + BrSALC = HOBr   
        Reaction(HBrUptkBySALA( State_Het ), [HBr], [BrSALA], [1], [1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  HBr = BrSALA   
        Reaction(HBrUptkBySALC( State_Het ), [HBr], [BrSALC], [1], [1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  HBr = BrSALC   
        Reaction(IuptkBySulf1stOrd( SR_MW(ine_HI), 0.10_ep, State_Het ), [HI], [AERI], [1], [1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  HI = AERI   
        Reaction(IuptkBySALA1stOrd( SR_MW(ine_HI), 0.10_ep, State_Het ), [HI], [ISALA], [1], [1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  HI = ISALA   
        Reaction(IuptkBySALC1stOrd( SR_MW(ine_HI), 0.10_ep, State_Het ), [HI], [ISALC], [1], [1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  HI = ISALC   
        Reaction(IuptkBySulf1stOrd( SR_MW(ine_I2O2), 0.02_ep, State_Het ), [I2O2], [AERI], [1], [2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2O2 = 2.000AERI   
        Reaction(IuptkBySALA1stOrd( SR_MW(ine_I2O2), 0.02_ep, State_Het ), [I2O2], [ISALA], [1], [2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2O2 = 2.000ISALA   
        Reaction(IuptkBySALC1stOrd( SR_MW(ine_I2O2), 0.02_ep, State_Het ), [I2O2], [ISALC], [1], [2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2O2 = 2.000ISALC   
        Reaction(IuptkBySulf1stOrd( SR_MW(ine_I2O3), 0.02_ep, State_Het ), [I2O3], [AERI], [1], [2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2O3 = 2.000AERI   
        Reaction(IuptkBySALA1stOrd( SR_MW(ine_I2O3), 0.02_ep, State_Het ), [I2O3], [ISALA], [1], [2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2O3 = 2.000ISALA   
        Reaction(IuptkBySALC1stOrd( SR_MW(ine_I2O3), 0.02_ep, State_Het ), [I2O3], [ISALC], [1], [2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2O3 = 2.000ISALC   
        Reaction(IuptkBySulf1stOrd( SR_MW(ine_I2O4), 0.02_ep, State_Het ), [I2O4], [AERI], [1], [2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2O4 = 2.000AERI   
        Reaction(IuptkBySALA1stOrd( SR_MW(ine_I2O4), 0.02_ep, State_Het ), [I2O4], [ISALA], [1], [2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2O4 = 2.000ISALA   
        Reaction(IuptkBySALC1stOrd( SR_MW(ine_I2O4), 0.02_ep, State_Het ), [I2O4], [ISALC], [1], [2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2O4 = 2.000ISALC   
        Reaction(IONO2uptkByH2O( State_Het ), [IONO2, H2O], [HOI, HNO3], [1, 1], [1, 1])   #2021/09/16 XW, TSherwen  #  IONO2 + H2O = HOI + HNO3   
        Reaction(IbrkdnByAcidBrSALA( SR_MW(ine_IONO), C(ine_IONO), 0.02_ep, State_Het ), [IONO, BrSALA], [IBr, HNO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2017;TS,JAS,SDE,XW  #  IONO + BrSALA = IBr + HNO2   
        Reaction(IbrkdnByAcidBrSALC( SR_MW(ine_IONO), C(ine_IONO), 0.02_ep, State_Het ), [IONO, BrSALC], [IBr, HNO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2017;TS,JAS,SDE,XW  #  IONO + BrSALC = IBr + HNO2   
        Reaction(IbrkdnByAcidSALACl( SR_MW(ine_IONO), C(ine_IONO), 0.02_ep, State_Het ), [IONO, SALACL], [ICl, HNO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2017;TS,JAS,SDE,XW  #  IONO + SALACL = ICl + HNO2   
        Reaction(IbrkdnByAcidSALCCl( SR_MW(ine_IONO), C(ine_IONO), 0.02_ep, State_Het ), [IONO, SALCCL], [ICl, HNO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2017;TS,JAS,SDE,XW  #  IONO + SALCCL = ICl + HNO2   
        Reaction(IbrkdnByAcidBrSALA( SR_MW(ine_IONO2), C(ine_IONO2), 0.01_ep, State_Het ), [IONO2, BrSALA], [IBr, HNO3], [1, 1], [1, 1])   #2017/09/22; Sherwen2017;TS,JAS,SDE,XW  #  IONO2 + BrSALA = IBr + HNO3   
        Reaction(IbrkdnByAcidBrSALC( SR_MW(ine_IONO2), C(ine_IONO2), 0.01_ep, State_Het ), [IONO2, BrSALC], [IBr, HNO3], [1, 1], [1, 1])   #2017/09/22; Sherwen2017;TS,JAS,SDE,XW  #  IONO2 + BrSALC = IBr + HNO3   
        Reaction(IbrkdnByAcidSALACl( SR_MW(ine_IONO2), C(ine_IONO2), 0.01_ep, State_Het ), [IONO2, SALACL], [ICl, HNO3], [1, 1], [1, 1])   #2017/09/22; Sherwen2017;TS,JAS,SDE,XW  #  IONO2 + SALACL = ICl + HNO3   
        Reaction(IbrkdnByAcidSALCCl( SR_MW(ine_IONO2), C(ine_IONO2), 0.01_ep, State_Het ), [IONO2, SALCCL], [ICl, HNO3], [1, 1], [1, 1])   #2017/09/22; Sherwen2017;TS,JAS,SDE,XW  #  IONO2 + SALCCL = ICl + HNO3   
        Reaction(IbrkdnByAcidBrSALA( SR_MW(ine_HOI), C(ine_HOI), 0.01_ep, State_Het ), [HOI, BrSALA], [IBr], [1, 1], [1])   #2017/09/22; Sherwen2017;TS,JAS,SDE,XW  #  HOI + BrSALA = IBr   
        Reaction(IbrkdnByAcidBrSALC( SR_MW(ine_HOI), C(ine_HOI), 0.01_ep, State_Het ), [HOI, BrSALC], [IBr], [1, 1], [1])   #2017/09/22; Sherwen2017;TS,JAS,SDE,XW  #  HOI + BrSALC = IBr   
        Reaction(IbrkdnByAcidSALACl( SR_MW(ine_HOI), C(ine_HOI), 0.01_ep, State_Het ), [HOI, SALACL], [ICl], [1, 1], [1])   #2017/09/22; Sherwen2017;TS,JAS,SDE,XW  #  HOI + SALACL = ICl   
        Reaction(IbrkdnByAcidSALCCl( SR_MW(ine_HOI), C(ine_HOI), 0.01_ep, State_Het ), [HOI, SALCCL], [ICl], [1, 1], [1])   #2017/09/22; Sherwen2017;TS,JAS,SDE,XW  #  HOI + SALCCL = ICl   
        Reaction(GLYXuptk1stOrd( SR_MW(ine_GLYX), State_Het), [GLYX], [SOAGX], [1], [1])   #2017/06/15; Marais2016, EAM  #  GLYX = SOAGX   
        Reaction(MGLYuptk1stOrd( SR_MW(ine_MGLY), State_Het), [MGLY], [SOAGX], [1], [1])   #2017/06/15; Marais2016, EAM  #  MGLY = SOAGX   
        Reaction(IEPOXuptk1stOrd( SR_MW(ine_IEPOXA), .FALSE., State_Het ), [IEPOXA], [SOAIE], [1], [1])   #2017/06/15; Marais2016, EAM  #  IEPOXA = SOAIE   
        Reaction(IEPOXuptk1stOrd( SR_MW(ine_IEPOXB), .FALSE., State_Het ), [IEPOXB], [SOAIE], [1], [1])   #2017/06/15; Marais2016, EAM  #  IEPOXB = SOAIE   
        Reaction(IEPOXuptk1stOrd( SR_MW(ine_IEPOXD), .FALSE., State_Het ), [IEPOXD], [SOAIE], [1], [1])   #2017/06/15; Marais2016, EAM  #  IEPOXD = SOAIE   
        Reaction(VOCuptk1stOrd( SR_MW(ine_LVOC), 1.0_ep, State_Het ), [LVOC], [LVOCOA], [1], [1])   #2017/06/15; Marais2016, EAM  #  LVOC = LVOCOA   
        Reaction(VOCuptk1stOrd( SR_MW(ine_MVKN), 5.0E-3_ep, State_Het ), [MVKN], [IONITA], [1], [1])   #2017/06/15; Marais2016, EAM  #  MVKN = IONITA   
        Reaction(VOCuptk1stOrd( SR_MW(ine_R4N2), 5.0E-3_ep, State_Het ), [R4N2], [IONITA], [1], [1])   #2017/06/15; Marais2016, EAM  #  R4N2 = IONITA   
        Reaction(VOCuptk1stOrd( SR_MW(ine_MONITS), 1.0E-2_ep, State_Het ), [MONITS], [MONITA], [1], [1])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS  #  MONITS = MONITA   
        Reaction(VOCuptk1stOrd( SR_MW(ine_MONITU), 1.0E-2_ep, State_Het ), [MONITU], [MONITA], [1], [1])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS  #  MONITU = MONITA   
        Reaction(VOCuptk1stOrd( SR_MW(ine_HONIT), 1.0E-2_ep, State_Het ), [HONIT], [MONITA], [1], [1])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS  #  HONIT = MONITA   
        Reaction(MGLYuptk1stOrd( SR_MW(ine_PYAC), State_Het ), [PYAC], [SOAGX], [1], [1])   #2019/11/06; Bates2019; KHB  #  PYAC = SOAGX   
        Reaction(IEPOXuptk1stOrd( SR_MW(ine_HMML), .TRUE., State_Het), [HMML], [SOAIE], [1], [1])   #2019/11/06; Bates2019; KHB  #  HMML = SOAIE   
        Reaction(VOCuptk1stOrd( SR_MW(ine_IHN1), 5.0E-3_ep, State_Het ), [IHN1], [IONITA], [1], [1])   #2019/11/06; Bates2019; KHB  #  IHN1 = IONITA   
        Reaction(VOCuptk1stOrd( SR_MW(ine_IHN2), 5.0E-2_ep, State_Het ), [IHN2], [IONITA], [1], [1])   #2019/11/06; Bates2019; KHB  #  IHN2 = IONITA   
        Reaction(VOCuptk1stOrd( SR_MW(ine_IHN3), 5.0E-3_ep, State_Het ), [IHN3], [IONITA], [1], [1])   #2019/11/06; Bates2019; KHB  #  IHN3 = IONITA   
        Reaction(VOCuptk1stOrd( SR_MW(ine_IHN4), 5.0E-3_ep, State_Het ), [IHN4], [IONITA], [1], [1])   #2019/11/06; Bates2019; KHB  #  IHN4 = IONITA   
        Reaction(IEPOXuptk1stOrd( SR_MW(ine_ICHE), .FALSE., State_Het ), [ICHE], [SOAIE], [1], [1])   #2019/11/06; Bates2019; KHB  #  ICHE = SOAIE   
        Reaction(VOCuptk1stOrd( SR_MW(ine_INPD), 5.0E-3_ep, State_Het ), [INPD], [IONITA], [1], [1])   #2019/11/06; Bates2019; KHB  #  INPD = IONITA   
        Reaction(VOCuptk1stOrd( SR_MW(ine_INPB), 5.0E-3_ep, State_Het ), [INPB], [IONITA], [1], [1])   #2019/11/06; Bates2019; KHB  #  INPB = IONITA   
        Reaction(VOCuptk1stOrd( SR_MW(ine_IDN), 5.0E-3_ep, State_Het ), [IDN], [IONITA], [1], [1])   #2019/11/06; Bates2019; KHB  #  IDN = IONITA   
        Reaction(VOCuptk1stOrd( SR_MW(ine_ITCN), 5.0E-3_ep, State_Het ), [ITCN], [IONITA], [1], [1])   #2019/11/06; Bates2019; KHB  #  ITCN = IONITA   
        Reaction(VOCuptk1stOrd( SR_MW(ine_ITHN), 5.0E-3_ep, State_Het ), [ITHN], [IONITA], [1], [1])   #2019/11/06; Bates2019; KHB  #  ITHN = IONITA   
        Reaction(VOCuptk1stOrd( SR_MW(ine_MCRHNB), 5.0E-3_ep, State_Het ), [MCRHNB], [IONITA], [1], [1])   #2019/11/06; Bates2019; KHB  #  MCRHNB = IONITA   
        Reaction(VOCuptk1stOrd( SR_MW(ine_MCRHN), 5.0E-3_ep, State_Het ), [MCRHN], [IONITA], [1], [1])   #2019/11/06; Bates2019; KHB  #  MCRHN = IONITA   
        Reaction(VOCuptk1stOrd( SR_MW(ine_NPHEN), 1.0E-2_ep, State_Het ), [NPHEN], [AONITA], [1], [1])  #  NPHEN = AONITA   


        Reaction(PHOTOL(2), [O3, hv], [O, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  O3 + hv = O + O2   
        Reaction(PHOTOL(3), [O3, hv], [O1D, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  O3 + hv = O1D + O2   
        Reaction(PHOTOL(1), [O2, hv], [O], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  O2 + hv = 2.000O   
        Reaction(PHOTOL(11), [NO2, hv], [NO, O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  NO2 + hv = NO + O   
        Reaction(PHOTOL(9), [H2O2, hv], [OH, OH], [1, 1], [1, 1])  #  H2O2 + hv = OH + OH   
        Reaction(PHOTOL(10), [MP, hv], [CH2O, HO2, OH], [1, 1], [1, 1, 1])  #  MP + hv = CH2O + HO2 + OH   
        Reaction(PHOTOL(7), [CH2O, hv], [HO2, H, CO], [1, 1], [1, 1, 1])   #2014/02/03; Eastham2014; SDE  #  CH2O + hv = HO2 + H + CO   
        Reaction(PHOTOL(8), [CH2O, hv], [H2, CO], [1, 1], [1, 1])  #  CH2O + hv = H2 + CO   
        Reaction(PHOTOL(16), [HNO3, hv], [OH, NO2], [1, 1], [1, 1])  #  HNO3 + hv = OH + NO2   
        Reaction(PHOTOL(15), [HNO2, hv], [OH, NO], [1, 1], [1, 1])  #  HNO2 + hv = OH + NO   
        Reaction(PHOTOL(17), [HNO4, hv], [OH, NO3], [1, 1], [1, 1])  #  HNO4 + hv = OH + NO3   
        Reaction(PHOTOL(18), [HNO4, hv], [HO2, NO2], [1, 1], [1, 1])  #  HNO4 + hv = HO2 + NO2   
        Reaction(PHOTOL(12), [NO3, hv], [NO2, O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  NO3 + hv = NO2 + O   
        Reaction(PHOTOL(13), [NO3, hv], [NO, O2], [1, 1], [1, 1])  #  NO3 + hv = NO + O2   
        Reaction(PHOTOL(14), [N2O5, hv], [NO3, NO2], [1, 1], [1, 1])  #  N2O5 + hv = NO3 + NO2   
        Reaction(PHOTOL(62), [ALD2, hv], [CH4, CO], [1, 1], [1, 1])  #  ALD2 + hv = CH4 + CO   
        Reaction(PHOTOL(76), [ACET, hv], [MCO3, MO2], [1, 1], [1, 1])  #  ACET + hv = MCO3 + MO2   
        Reaction(PHOTOL(77), [ACET, hv], [MO2, CO], [1, 1], [2.000, 1])  #  ACET + hv = 2.000MO2 + CO   
        Reaction(PHOTOL(72), [GLYX, hv], [HO2, CO], [1, 1], [2.000, 2.000])  #  GLYX + hv = 2.000HO2 + 2.000CO   
        Reaction(PHOTOL(73), [GLYX, hv], [H2, CO], [1, 1], [1, 2.000])  #  GLYX + hv = H2 + 2.000CO   
        Reaction(PHOTOL(74), [GLYX, hv], [CH2O, CO], [1, 1], [1, 1])  #  GLYX + hv = CH2O + CO   
        Reaction(PHOTOL(71), [MGLY, hv], [MCO3, CO, HO2], [1, 1], [1, 1, 1])  #  MGLY + hv = MCO3 + CO + HO2   
        Reaction(PHOTOL(63), [MVK, hv], [PRPE, CO], [1, 1], [1, 1])  #  MVK + hv = PRPE + CO   
        Reaction(PHOTOL(64), [MVK, hv], [MCO3, CH2O, CO, HO2], [1, 1], [1, 1, 1, 1])  #  MVK + hv = MCO3 + CH2O + CO + HO2   
        Reaction(PHOTOL(65), [MVK, hv], [MO2, RCO3], [1, 1], [1, 1])   #2014/05/23; Eastham2014; JMAO,SDE  #  MVK + hv = MO2 + RCO3   
        Reaction(PHOTOL(66), [MACR, hv], [CO, HO2, CH2O, MCO3], [1, 1], [1, 1, 1, 1])   #2014/05/23; Eastham2014; JMAO,SDE  #  MACR + hv = CO + HO2 + CH2O + MCO3   
        Reaction(PHOTOL(75), [HAC, hv], [MCO3, CH2O, HO2], [1, 1], [1, 1, 1])  #  HAC + hv = MCO3 + CH2O + HO2   
        Reaction(PHOTOL(79), [PRPN, hv], [OH, HO2, RCHO, NO2], [1, 1], [1, 1, 1, 1])  #  PRPN + hv = OH + HO2 + RCHO + NO2   
        Reaction(PHOTOL(80), [ETP, hv], [OH, HO2, ALD2], [1, 1], [1, 1, 1])  #  ETP + hv = OH + HO2 + ALD2   
        Reaction(PHOTOL(81), [RA3P, hv], [OH, HO2, RCHO], [1, 1], [1, 1, 1])  #  RA3P + hv = OH + HO2 + RCHO   
        Reaction(PHOTOL(82), [RB3P, hv], [OH, HO2, ACET], [1, 1], [1, 1, 1])  #  RB3P + hv = OH + HO2 + ACET   
        Reaction(PHOTOL(83), [R4P, hv], [OH, HO2, RCHO], [1, 1], [1, 1, 1])  #  R4P + hv = OH + HO2 + RCHO   
        Reaction(PHOTOL(84), [PP, hv], [OH, HO2, ALD2, CH2O], [1, 1], [1, 1, 1, 1])  #  PP + hv = OH + HO2 + ALD2 + CH2O   
        Reaction(PHOTOL(85), [RP, hv], [OH, HO2, ALD2], [1, 1], [1, 1, 1])  #  RP + hv = OH + HO2 + ALD2   
        Reaction(PHOTOL(99), [MAP, hv], [OH, MO2], [1, 1], [1, 1])  #  MAP + hv = OH + MO2   
        Reaction(PHOTOL(23), [Br2, hv], [Br], [1, 1], [2.000])   #2012/06/07; Parrella2012; JPP  #  Br2 + hv = 2.000Br   
        Reaction(PHOTOL(28), [BrO, hv], [Br, O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  BrO + hv = Br + O   
        Reaction(PHOTOL(32), [HOBr, hv], [Br, OH], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  HOBr + hv = Br + OH   
        Reaction(PHOTOL(29), [BrNO3, hv], [Br, NO3], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  BrNO3 + hv = Br + NO3   
        Reaction(PHOTOL(30), [BrNO3, hv], [BrO, NO2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  BrNO3 + hv = BrO + NO2   
        Reaction(PHOTOL(31), [BrNO2, hv], [Br, NO2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  BrNO2 + hv = Br + NO2   
        Reaction(PHOTOL(56), [CHBr3, hv], [Br], [1, 1], [3.000])   #2012/06/07; Parrella2012; JPP  #  CHBr3 + hv = 3.000Br   
        Reaction(PHOTOL(55), [CH2Br2, hv], [Br], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  CH2Br2 + hv = 2.000Br   
        Reaction(PHOTOL(50), [CH3Br, hv], [MO2, Br], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  CH3Br + hv = MO2 + Br   
        Reaction(PHOTOL(43), [CH3Cl, hv], [MO2, Cl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  CH3Cl + hv = MO2 + Cl   
        Reaction(PHOTOL(45), [CH2Cl2, hv], [Cl], [1, 1], [2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  CH2Cl2 + hv = 2.000Cl   
        Reaction(PHOTOL(33), [BrCl, hv], [Br, Cl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  BrCl + hv = Br + Cl   
        Reaction(PHOTOL(22), [Cl2, hv], [Cl], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  Cl2 + hv = 2.000Cl   
        Reaction(PHOTOL(27), [ClO, hv], [Cl, O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + hv = Cl + O   
        Reaction(PHOTOL(25), [OClO, hv], [ClO, O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OClO + hv = ClO + O   
        Reaction(PHOTOL(26), [Cl2O2, hv], [Cl, ClOO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  Cl2O2 + hv = Cl + ClOO   
        Reaction(PHOTOL(21), [ClNO2, hv], [Cl, NO2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClNO2 + hv = Cl + NO2   
        Reaction(PHOTOL(19), [ClNO3, hv], [Cl, NO3], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClNO3 + hv = Cl + NO3   
        Reaction(PHOTOL(20), [ClNO3, hv], [ClO, NO2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClNO3 + hv = ClO + NO2   
        Reaction(PHOTOL(24), [HOCl, hv], [Cl, OH], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  HOCl + hv = Cl + OH   
        Reaction(PHOTOL(44), [CH3CCl3, hv], [Cl], [1, 1], [3.000])   #2014/02/03; Eastham2014; SDE  #  CH3CCl3 + hv = 3.000Cl   
        Reaction(PHOTOL(42), [CCl4, hv], [Cl], [1, 1], [4.000])   #2014/02/03; Eastham2014; SDE  #  CCl4 + hv = 4.000Cl   
        Reaction(PHOTOL(37), [CFC11, hv], [Cl], [1, 1], [3.000])   #2014/02/03; Eastham2014; SDE  #  CFC11 + hv = 3.000Cl   
        Reaction(PHOTOL(38), [CFC12, hv], [Cl], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  CFC12 + hv = 2.000Cl   
        Reaction(PHOTOL(39), [CFC113, hv], [Cl], [1, 1], [3.000])   #2014/02/03; Eastham2014; SDE  #  CFC113 + hv = 3.000Cl   
        Reaction(PHOTOL(40), [CFC114, hv], [Cl], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  CFC114 + hv = 2.000Cl   
        Reaction(PHOTOL(41), [CFC115, hv], [Cl], [1, 1], [1])   #2014/02/03; Eastham2014; SDE  #  CFC115 + hv = Cl   
        Reaction(PHOTOL(47), [HCFC123, hv], [Cl], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  HCFC123 + hv = 2.000Cl   
        Reaction(PHOTOL(48), [HCFC141b, hv], [Cl], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  HCFC141b + hv = 2.000Cl   
        Reaction(PHOTOL(49), [HCFC142b, hv], [Cl], [1, 1], [1])   #2014/02/03; Eastham2014; SDE  #  HCFC142b + hv = Cl   
        Reaction(PHOTOL(46), [HCFC22, hv], [Cl], [1, 1], [1])   #2014/02/03; Eastham2014; SDE  #  HCFC22 + hv = Cl   
        Reaction(PHOTOL(53), [H1301, hv], [Br], [1, 1], [1])   #2014/02/03; Eastham2014; SDE  #  H1301 + hv = Br   
        Reaction(PHOTOL(51), [H1211, hv], [Cl, Br], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  H1211 + hv = Cl + Br   
        Reaction(PHOTOL(54), [H2402, hv], [Br], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  H2402 + hv = 2.000Br   
        Reaction(PHOTOL(101), [ClOO, hv], [Cl, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClOO + hv = Cl + O2   
        Reaction(PHOTOL(114), [I2, hv], [I], [1, 1], [2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2 + hv = 2.000I   
        Reaction(PHOTOL(115), [HOI, hv], [I, OH], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  HOI + hv = I + OH   
        Reaction(PHOTOL(116), [IO, hv], [I, O], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + hv = I + O   
        Reaction(PHOTOL(117), [OIO, hv], [I, O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  OIO + hv = I + O2   
        Reaction(PHOTOL(118), [INO, hv], [I, NO], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  INO + hv = I + NO   
        Reaction(PHOTOL(119), [IONO, hv], [I, NO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IONO + hv = I + NO2   
        Reaction(PHOTOL(120), [IONO2, hv], [I, NO3], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IONO2 + hv = I + NO3   
        Reaction(PHOTOL(121), [I2O2, hv], [I, OIO], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2O2 + hv = I + OIO   
        Reaction(PHOTOL(122), [CH3I, hv], [I], [1, 1], [1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  CH3I + hv = I   
        Reaction(PHOTOL(123), [CH2I2, hv], [I], [1, 1], [2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  CH2I2 + hv = 2.000I   
        Reaction(PHOTOL(124), [CH2ICl, hv], [I, Cl], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  CH2ICl + hv = I + Cl   
        Reaction(PHOTOL(125), [CH2IBr, hv], [I, Br], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  CH2IBr + hv = I + Br   
        Reaction(PHOTOL(126), [I2O4, hv], [OIO], [1, 1], [2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2O4 + hv = 2.000OIO   
        Reaction(PHOTOL(127), [I2O3, hv], [OIO, IO], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2O3 + hv = OIO + IO   
        Reaction(PHOTOL(128), [IBr, hv], [I, Br], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IBr + hv = I + Br   
        Reaction(PHOTOL(129), [ICl, hv], [I, Cl], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  ICl + hv = I + Cl   
        Reaction(PHOTOL(103), [MPN, hv], [CH2O, NO3, HO2], [1, 1], [1, 1, 1])   #2012/02/12; Browne2011; ECB  #  MPN + hv = CH2O + NO3 + HO2   
        Reaction(PHOTOL(104), [MPN, hv], [MO2, NO2], [1, 1], [1, 1])   #2012/02/12; Browne2011; ECB  #  MPN + hv = MO2 + NO2   
        Reaction(PHOTOL(97), [ATOOH, hv], [OH, CH2O, MCO3], [1, 1], [1, 1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  ATOOH + hv = OH + CH2O + MCO3   
        Reaction(PHOTOL(36), [N2O, hv], [N2, O1D], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  N2O + hv = N2 + O1D   
        Reaction(PHOTOL(34), [OCS, hv], [SO2, CO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OCS + hv = SO2 + CO   
        Reaction(PHOTOL(100), [SO4, hv], [SO2, OH], [1, 1], [1, 2.000])   #2014/02/03; Eastham2014; SDE  #  SO4 + hv = SO2 + 2.000OH   
        Reaction(PHOTOL(6), [NO, hv], [O, N], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  NO + hv = O + N   
        Reaction(PHOTOL(105), [PIP, hv], [RCHO, OH, HO2], [1, 1], [1, 1, 1])   #2017/03/23; Fischer2014; EVF  #  PIP + hv = RCHO + OH + HO2   
        Reaction(PHOTOL(107), [ETHLN, hv], [NO2, CH2O, CO, HO2], [1, 1], [1, 1, 1, 1])   #2017/06/15; Marais2016; EAM  #  ETHLN + hv = NO2 + CH2O + CO + HO2   
        Reaction(PHOTOL(111), [MONITS, hv], [MEK, NO2], [1, 1], [1, 1])   #2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS  #  MONITS + hv = MEK + NO2   
        Reaction(PHOTOL(112), [MONITU, hv], [RCHO, NO2], [1, 1], [1, 1])   #2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS  #  MONITU + hv = RCHO + NO2   
        Reaction(PHOTOL(113), [HONIT, hv], [HAC, NO2], [1, 1], [1, 1])   #2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS  #  HONIT + hv = HAC + NO2   
        Reaction(PHOTOL(130), [NITs, hv], [HNO2], [1, 1], [1])   #2018/07/19; Kasibhatla2018; PK, TMS  #  NITs + hv = HNO2   
        Reaction(PHOTOL(131), [NITs, hv], [NO2], [1, 1], [1])   #2018/07/19; Kasibhatla2018; PK, TMS  #  NITs + hv = NO2   
        Reaction(PHOTOL(132), [NIT, hv], [HNO2], [1, 1], [1])   #2018/07/19; Kasibhatla2018; PK, TMS  #  NIT + hv = HNO2   
        Reaction(PHOTOL(133), [NIT, hv], [NO2], [1, 1], [1])   #2018/07/19; Kasibhatla2018; PK, TMS  #  NIT + hv = NO2   
        Reaction(PHOTOL(134), [MENO3, hv], [NO2, HO2, CH2O], [1, 1], [1, 1, 1])   #2019/07/11; Fisher2018; JAF  #  MENO3 + hv = NO2 + HO2 + CH2O   
        Reaction(PHOTOL(135), [ETNO3, hv], [NO2, HO2, ALD2], [1, 1], [1, 1, 1])   #2019/07/11; Fisher2018; JAF  #  ETNO3 + hv = NO2 + HO2 + ALD2   
        Reaction(PHOTOL(136), [IPRNO3, hv], [NO2, HO2, ACET], [1, 1], [1, 1, 1])   #2019/07/11; Fisher2018; JAF  #  IPRNO3 + hv = NO2 + HO2 + ACET   
        Reaction(PHOTOL(137), [NPRNO3, hv], [NO2, HO2, RCHO], [1, 1], [1, 1, 1])   #2019/07/11; Fisher2018; JAF  #  NPRNO3 + hv = NO2 + HO2 + RCHO   
        Reaction(PHOTOL(86), [HMHP, hv], [OH, CH2O], [1, 1], [2, 1])   #2019/11/06; Bates2019; KHB  #  HMHP + hv = 2OH + CH2O   
        Reaction(PHOTOL(87), [HPETHNL, hv], [OH, CO, HO2, CH2O], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  HPETHNL + hv = OH + CO + HO2 + CH2O   
        Reaction(PHOTOL(88), [PYAC, hv], [MCO3, CO2, HO2], [1, 1], [1, 1, 1])   #2019/11/06; Bates2019; KHB  #  PYAC + hv = MCO3 + CO2 + HO2   
        Reaction(PHOTOL(89), [PROPNN, hv], [NO2, CH2O, MCO3], [1, 1], [1, 1, 1])   #2019/11/06; Bates2019; KHB  #  PROPNN + hv = NO2 + CH2O + MCO3   
        Reaction(PHOTOL(90), [MVKHC, hv], [CO, HO2, CH2O, MCO3], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  MVKHC + hv = CO + HO2 + CH2O + MCO3   
        Reaction(PHOTOL(109), [MCRHN, hv], [HAC, CO, HO2, NO2], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  MCRHN + hv = HAC + CO + HO2 + NO2   
        Reaction(PHOTOL(110), [MCRHNB, hv], [PROPNN, OH, CO, HO2], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  MCRHNB + hv = PROPNN + OH + CO + HO2   
        Reaction(PHOTOL(138), [RIPA, hv], [MVK, CH2O, HO2, OH], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  RIPA + hv = MVK + CH2O + HO2 + OH   
        Reaction(PHOTOL(139), [RIPB, hv], [MACR, CH2O, HO2, OH], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  RIPB + hv = MACR + CH2O + HO2 + OH   
        Reaction(PHOTOL(140), [RIPC, hv], [OH, HO2, HC5A], [1, 1], [1, 1, 1])   #2019/11/06; Bates2019; KHB  #  RIPC + hv = OH + HO2 + HC5A   
        Reaction(PHOTOL(141), [RIPD, hv], [OH, HO2, HC5A], [1, 1], [1, 1, 1])   #2019/11/06; Bates2019; KHB  #  RIPD + hv = OH + HO2 + HC5A   
        Reaction(PHOTOL(144), [HPALD3, hv], [CO, OH, HO2, MVK], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  HPALD3 + hv = CO + OH + HO2 + MVK   
        Reaction(PHOTOL(145), [HPALD4, hv], [CO, OH, HO2, MACR], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  HPALD4 + hv = CO + OH + HO2 + MACR   
        Reaction(PHOTOL(147), [IHN2, hv], [NO2, MVK, HO2, CH2O], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  IHN2 + hv = NO2 + MVK + HO2 + CH2O   
        Reaction(PHOTOL(148), [IHN3, hv], [NO2, MACR, HO2, CH2O], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  IHN3 + hv = NO2 + MACR + HO2 + CH2O   
        Reaction(PHOTOL(152), [INPD, hv], [NO2, IHOO1, IHOO4], [1, 1], [1, 0.841, 0.159])   #2019/11/06; Bates2019; KHB  #  INPD + hv = NO2 + 0.841IHOO1 + 0.159IHOO4   
        Reaction(PHOTOL(160), [ITCN, hv], [MGLY, OH, NO2, GLYC], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  ITCN + hv = MGLY + OH + NO2 + GLYC   
        Reaction(PHOTOL(162), [ETHP, hv], [ETO, OH], [1, 1], [1, 1])   #2021/09/22; Bates2021a; KHB,MSL  #  ETHP + hv = ETO + OH   
        Reaction(PHOTOL(163), [BALD, hv], [BENZO2, CO, HO2], [1, 1], [1, 1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  BALD + hv = BENZO2 + CO + HO2   
        Reaction(PHOTOL(164), [BZCO3H, hv], [BENZO2, OH, CO2], [1, 1], [1, 1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  BZCO3H + hv = BENZO2 + OH + CO2   
        Reaction(PHOTOL(165), [BENZP, hv], [BENZO], [1, 1], [1])   #2021/09/29; Bates2021b; KHB,MSL  #  BENZP + hv = BENZO   



    ]
end 
     
