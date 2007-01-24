"""
fields.py

Fields module corresponding to CPL6's mod.F90

This module is just a long list of the many constants defined for
the CCSM Coupler.
"""

# Component Names
atmname = 'atm'
ocnname = 'ocn'
icename = 'ice'
lndname = 'lnd'
rtmname = 'roff'
cplname = 'cpl'

component_names = [atmname, ocnname, icename, lndname, rtmname, cplname]

# The functions need oython analogues
#getField     #  returns string for nth aVect attribute
#getLongName  #  returns netCDF longname and unit strings
#getField     #  returns string for nth aVect attribute
#getLongName  #  returns netCDF longname and unit strings

# Define dictionaries for buffer descriptors

###
# integer buffer
###
ibuf = {}
ibuf["total"]    = 100  #  size of info-buffer
# These values are fortran array indices, they are
# all starting from 1.  However, this causes
# the values to be set incorrect
# ( i.e.: ibuf(1) in fortran is rcode,
# but rcode is ibuf[0] in python.
# After assigning these values I will
# subtract one from each.
#
# dict entry = fortran value - adjustment to python value
ibuf["rcode"]    =   1 - 1 #  error code
ibuf["cdate"]    =   2 - 1 #  current date: yymmdd
ibuf["sec"]      =   3 - 1 #  elapsed sec on date
ibuf["ncpl"]     =   4 - 1 #  cpl comm's per day
ibuf["nfields"]  =  10 - 1 
ibuf["gsize"]    =  11 - 1
ibuf["lsize"]    =  12 - 1
ibuf["gisize"]   =  13 - 1
ibuf["gjsize"]   =  14 - 1
ibuf["lisize"]   =  15 - 1
ibuf["ljsize"]   =  16 - 1
ibuf["stopeod"]  =  19 - 1
ibuf["stopnow"]  =  20 - 1
ibuf["resteod"]  =  21 - 1
ibuf["restnow"]  =  22 - 1
ibuf["histeod"]  =  23 - 1
ibuf["histnow"]  =  24 - 1
ibuf["histavg"]  =  25 - 1
ibuf["diageod"]  =  26 - 1
ibuf["diagnow"]  =  27 - 1
ibuf["infotim"]  =  28 - 1
ibuf["infobug"]  =  29 - 1
ibuf["precadj"]  =  31 - 1 #  precip adjustment factor (* 1.0e+6)
ibuf["ashift"]   =  32 - 1 #  albedo calculation time shift
ibuf["nbasins"]  =  33 - 1 #  number of active runoff basins
ibuf["xalbic"]   =  34 - 1 #  request extra albedo solar init msg
ibuf["inimask"]  =  36 - 1 #  flag cpl to send back domain spec for lnd
ibuf["dead"]     =  37 - 1 #  non-0 <=> dead model
ibuf["domain"]   =  40 - 1
ibuf["userest"]  =  41 - 1 #  non-0 <=> use restart data sent to cpl

###
# real buffer
###
rbuf = {}
rbuf["total"]    =  50  #  size of real info-buffer
# dict entry = fortran value - adjustment to python value
rbuf["spval"]    =   1 - 1 #  the special value
rbuf["eccen"]    =  10 - 1 #  Earth's eccentricity
rbuf["obliqr"]   =  11 - 1 #  Earth's Obliquity
rbuf["lambm0"]   =  12 - 1 #  longitude of perihelion at v-equinox
rbuf["mvelpp"]   =  13 - 1  #  Earth's Moving vernal equinox of orbit +pi

###
# grid fields
###
grid_total  = 7
grid_fields = ["lat","lon","area","aream","index","mask","pid"]

grid_indices = {}

i = 0
for field in grid_fields:
    grid_indices[field] = i
    i += 1
    
grid_lat    = 1     #  lat from component
grid_lon    = 2     #  lon from component
grid_area   = 3     #  area from component
grid_aream  = 4     #  area from mapping file
grid_index  = 5     #  global index
grid_mask   = 6     #  mask, 0 = inactive cell
grid_pid    = 7     #  proc id number

###
# a2c
###
a2c_total  = 19
a2c_states = ["Sa_z","Sa_u","Sa_v","Sa_tbot","Sa_ptem","Sa_shum",
              "Sa_dens","Sa_pbot","Sa_pslv"]
a2c_fluxes = ["Faxa_lwdn","Faxa_rainc","Faxa_rainl","Faxa_snowc",
              "Faxa_snowl","Faxa_swndr","Faxa_swvdr","Faxa_swndf",
              "Faxa_swvdf","Faxa_swnet"]
a2c_fields = a2c_states + a2c_fluxes

a2c_indices = {}
i = 0
for field in a2c_fields:
    a2c_indices[field] = i
    i += 1

a2c_z      =  1  #  bottom atm level height
a2c_u      =  2  #  bottom atm level zon wind
a2c_v      =  3  #  bottom atm level mer wind
a2c_tbot   =  4  #  bottom atm level temp
a2c_ptem   =  5  #  bottom atm level pot temp
a2c_shum   =  6  #  bottom atm level spec hum
a2c_dens   =  7  #  bottom atm level air den
a2c_pbot   =  8  #  bottom atm level pressure
a2c_pslv   =  9  #  sea level atm pressure
a2c_lwdn   = 10  #  downward lw heat flux
a2c_rainc  = 11  #  prec: liquid "convective"
a2c_rainl  = 12  #  prec: liquid "large scale"
a2c_snowc  = 13  #  prec: frozen "convective"
a2c_snowl  = 14  #  prec: frozen "large scale"
a2c_swndr  = 15  #  sw: nir direct  downward
a2c_swvdr  = 16  #  sw: vis direct  downward
a2c_swndf  = 17  #  sw: nir diffuse downward
a2c_swvdf  = 18  #  sw: vis diffuse downward
a2c_swnet  = 19  #  sw: net


###
# c2a
###
c2a_total  = 17
c2a_states = ["Sx_tref","Sx_qref","Sx_avsdr","Sx_anidr","Sx_avsdf","Sx_anidf",
              "Sx_t","So_t","Sx_snowh","Sx_ifrac","Sx_ofrac"]
c2a_fluxes = ["Faxx_taux","Faxx_tauy","Faxx_lat","Faxx_sen","Faxx_lwup","Faxx_evap"]
c2a_fields = c2a_states + c2a_fluxes

c2a_indices = {}
i = 0
for field in c2a_fields:
    c2a_indices[field] = i
    i += 1

c2a_tref  =  1  #  2m reference temperature
c2a_qref  =  2  #  2m reference specific humidity
c2a_avsdr =  3  #  albedo, visible, direct
c2a_anidr =  4  #  albedo, near-ir, direct
c2a_avsdf =  5  #  albedo, visible, diffuse
c2a_anidf =  6  #  albedo, near-ir, diffuse
c2a_t     =  7  #  surface temperature
c2a_sst   =  8  #  sea surface temperature
c2a_snowh =  9  #  surface snow depth
c2a_ifrac = 10  #  surface ice fraction
c2a_ofrac = 11  #  surface ocn fraction
c2a_taux  = 12  #  wind stress, zonal
c2a_tauy  = 13  #  wind stress, meridional
c2a_lat   = 14  #  latent          heat flux
c2a_sen   = 15  #  sensible        heat flux
c2a_lwup  = 16  #  upward longwave heat flux
c2a_evap  = 17  #  evaporation    water flux

###
# ice to coupler fields
###
i2c_total  = 22

i2c_states = ["Si_t","Si_tref","Si_qref","Si_ifrac","Si_avsdr",
              "Si_anidr","Si_avsdf","Si_anidf","index"]
i2c_fluxes = ["Faii_taux","Faii_tauy","Faii_lat","Faii_sen",
              "Faii_lwup","Faii_evap","Faii_swnet","Fioi_swpen",
              "Fioi_melth","Fioi_meltw","Fioi_salt","Fioi_taux",
              "Fioi_tauy"]
i2c_fields = i2c_states + i2c_fluxes

i2c_indices = {}
i = 0
for field in i2c_fields:
    i2c_indices[field] = i
    i += 1

i2c_t     =  1  #  temperature
i2c_tref  =  2  #  2m reference temperature
i2c_qref  =  3  #  2m reference specific humidity
i2c_ifrac =  4  #  fractional ice coverage
i2c_avsdr =  5  #  albedo: visible, direct
i2c_anidr =  6  #  albedo: near ir, direct
i2c_avsdf =  7  #  albedo: visible, diffuse
i2c_anidf =  8  #  albedo: near ir, diffuse
i2c_index =  9  #  global data compr index
i2c_taux  = 10  #  wind stress, zonal
i2c_tauy  = 11  #  wind stress, meridional
i2c_lat   = 12  #  latent          heat flux
i2c_sen   = 13  #  sensible        heat flux
i2c_lwup  = 14  #  upward longwave heat flux
i2c_evap  = 15  #  evaporation    water flux
i2c_swnet = 16  #  shortwave: net absorbed
i2c_swpen = 17  #  net SW penetrating ice
i2c_melth = 18  #  heat  flux from melting ice
i2c_meltw = 19  #  water flux from melting ice
i2c_salt  = 20  #  salt  flux from melting ice
i2c_otaux = 21  #  ice/ocn stress, zonal
i2c_otauy = 22  #  ice/ocn stress, meridional


###
# c2i
###
c2i_total  = 21
c2i_states = ["So_t","So_s","So_u","So_v","Sa_z","Sa_u","Sa_v",
              "Sa_ptem","Sa_tbot","Sa_shum","Sa_dens","Sa_dhdx","Sa_dhdy"]
c2i_fluxes = ["Fioo_q","Faxa_swndr","Faxa_swvdr","Faxa_swndf","Faxa_swvdf",
              "Faxa_lwdn","Faxc_rain","Faxc_snow"]
c2i_fields = c2i_states + c2i_fluxes

c2i_indices = {}
i = 0
for field in c2i_fields:
    c2i_indices[field] = i
    i += 1
    
c2i_ot    =  1  #  ocn temp
c2i_os    =  2  #  ocn salinity
c2i_ou    =  3  #  ocn u velocity
c2i_ov    =  4  #  ocn v velocity
c2i_dhdx  = 12  #  ocn surface slope, zonal
c2i_dhdy  = 13  #  ocn surface slope, merid
c2i_z     =  5  #  atm bottom layer height
c2i_u     =  6  #  atm u velocity
c2i_v     =  7  #  atm v velocity
c2i_ptem  =  8  #  atm potential temp
c2i_tbot  =  9  #  atm bottom temp
c2i_shum  = 10  #  atm specfic humidity
c2i_dens  = 11  #  atm air density
c2i_q     = 14  #  ocn freeze or melt heat
c2i_swndr = 15  #  atm sw near-ir, direct
c2i_swvdr = 16  #  atm sw visable, direct
c2i_swndf = 17  #  atm sw near-ir, diffuse
c2i_swvdf = 18  #  atm sw visable, diffuse
c2i_lwdn  = 19  #  long-wave down
c2i_rain  = 20  #  rain
c2i_snow  = 21  #  snow

###
# Land Fields
###

l2c_total  = 15
l2c_states = ["Sl_t","Sl_tref","Sl_qref","Sl_avsdr","Sl_anidr",
              "Sl_avsdf","Sl_anidf","Sl_snow"]
l2c_fluxes = ["Fall_taux","Fall_tauy","Fall_lat","Fall_sen","Fall_lwup",
              "Fall_evap","Fall_swnet"]
l2c_fields = l2c_states + l2c_fluxes

l2c_indices = {}
i = 0
for field in l2c_fields:
    l2c_indices[field] = i
    i+=1

l2c_t     =  1  #  temperature
l2c_tref  =  2  #  2m reference temperature
l2c_qref  =  3  #  2m reference specific humidity
l2c_avsdr =  4  #  albedo: direct , visible
l2c_anidr =  5  #  albedo: direct , near-ir
l2c_avsdf =  6  #  albedo: diffuse, visible
l2c_anidf =  7  #  albedo: diffuse, near-ir
l2c_snowh =  8  #  snow height
l2c_taux  =  9  #  wind stress, zonal
l2c_tauy  = 10  #  wind stress, meridional
l2c_lat   = 11  #  latent          heat flux
l2c_sen   = 12  #  sensible        heat flux
l2c_lwup  = 13  #  upward longwave heat flux
l2c_evap  = 14  #  evaporation    water flux
l2c_swnet = 15  #  2m reference temperature

c2l_total  = 18
c2l_states = ["Sa_z","Sa_u","Sa_v","Sa_tbot","Sa_ptem","Sa_shum","Sa_dens","Sa_pbot","Sa_pslv"]
c2l_fluxes = ["Faxa_lwdn","Faxa_rainc","Faxa_rainl","Faxa_snwoc","Faxa_snowl",
              "Faxa_swndr","Faxa_swvdr","Faxa_swndf","Faxa_swvdf"]
c2l_fields = c2l_states + c2l_fluxes

c2l_indices = {}
i = 0
for field in c2l_fields:
    c2l_indices[field] = i
    i+=1

# atm states
c2l_z     =  1  #  bottom atm level height
c2l_u     =  2  #  bottom atm level zon wind
c2l_v     =  3  #  bottom atm level mer wind
c2l_tbot = 4 # bottom atm level temp
c2l_ptem = 5 # bottom atm level pot temp
c2l_shum = 6 # bottom atm level spec hum
c2l_dens = 7 # bottom atm level air dens
c2l_pbot = 8 # bottom atm level pressure
c2l_pslv = 9 # sea level atm pressure
# Computed by atm
c2l_lwdn = 10 # downward longwave heat flux
c2l_rainc = 11 # precip: liquid, convective
c2l_rainl = 12 # precip: liquid, large-scale
c2l_snowc = 13 # precip: frozen, convective
c2l_snowl = 14 # precip: frozen, large-scale
c2l_swndr = 15 # shortwave: nir direct down
c2l_swvdr = 16 # shortwave: vis direct down
c2l_swndf = 17 # shortwave: nir diffuse down
c2l_swvdf = 18 # shortwave: vis diffuse down

# Special land grid initialization:
c2lg_total = 6
c2lg_fields = ["lon","lat","area","lfrac","maskl","maska"]
c2lg_indices = {}
i = 0
for field in c2lg_fields:
    c2lg_indices[field] = i
    i += 1

c2lg_alon = 1 # longitude
c2lg_alat = 2 # latitude
c2lg_aarea = 3 # cell area
c2lg_lfrac = 4 # lnd_fraction
c2lg_lmask = 5 # lnd mask
c2lg_amask = 6 # atm mask

###
# OCN Fields
###
o2c_total = 7
o2c_states = ["So_t","So_u","So_v","So_s","So_dhdx","So_dhdy"]
o2c_fluxes = ["Fioo_q"]
o2c_fields = o2c_states + o2c_fluxes

o2c_indices = {}
i = 0
for field in o2c_fields:
    o2c_indices[field] = i
    i += 1

# OCN States:
o2c_t = 1 # temperature
o2c_u = 2 # velocity, zonal
o2c_v = 3 # velocity, meridional
o2c_s = 4 # salinity
o2c_dhdx = 5 # surface slope, zonal
o2c_dhdy = 6 # surface slope, meridional
o2c_q = 7 # heat of fusion (q>0) melt pot

###
# c2o
###
c2o_total = 18
c2o_states = ["Si_ifrac","Sa_pslv","Faoc_duu10n"]
c2o_fluxes = ["Foxx_taux","Foxx_tauy","Foxx_swnet","Foxx_lat",
              "Foxx_sen","Foxx_lwup","Foxx_lwdn","Foxx_melth",
              "Foxx_salt","Foxx_prec","Foxx_snow","Foxx_rain",
              "Foxx_evap","Foxx_meltw","Forr_roff"]
c2o_fields = c2o_states + c2o_fluxes

c2o_indices = {}
i = 0
for field in c2o_fields:
    c2o_indices[field] = i
    i += 1

# OCN Model Input:
c2o_ifrac = 1 # state: ice fraction
c2o_press = 2 # state: sea level pressure
c2o_duu10 = 3 # state: 10m wind speed squared
c2o_taux = 4 # wind stress: zonal
c2o_tauy = 5 # wind stress: meridional
c2o_swnet = 6 # heat flux: shortwave net
c2o_lat = 7 # heat flux: latent
c2o_sen = 8 # heat flux: sensible
c2o_lwup = 9 # heat flux: long-wave up
c2o_lwdn = 10 # heat flux: long-wave down
c2o_melth = 11 # heat flux: melt
c2o_salt = 12 # salt flux
c2o_prec = 13 # water flux: rain + snow
c2o_snow = 14 # water flux: snow
c2o_rain = 15 # water flux: rain 
c2o_evap = 16 # water flux: evap
c2o_meltw = 17 # water flux: melt
c2o_roff = 18 # water flux: runoff

### 
# Run-off field
###

r2c_total = 1
r2c_states = []
r2c_fluxes = ["Forr_roff"]
r2c_fields = r2c_states + r2c_fluxes 
r2c_indices = {}
r2c_indices["Forr_roff"] = 1
r2c_runoff = 1

# Finished declarations, clean up variables we don't want
# visiable in the module
del i


