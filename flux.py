"""
flux.py

contains routines that do coupler's flux calculations
"""

# Standard Imports:
import sys,math

# Additional required modules:
import Numeric
import mpi

# CPL6 modules:
import cpl
import cpl.const
import cpl.control
import cpl.fields
import cpl.domain
import cpl.attributevector

import cpl.shr
import cpl.shr.cal
import cpl.shr.date
import cpl.shr.orb

# Definitions
FluxAlboFirstCall = True
FluxAlbiFirstCall = True
OceanAlbedoFirstCall = True
FluxEpbalFirstCall = True


AOFLUX_FIELDS = ["Faoc_sen",
                 "Faoc_lat",
                 "Faoc_lwup",
                 "Faoc_evap",
                 "Faoc_taux",
                 "Faoc_tauy",
                 "Faoc_tref",
                 "Faoc_qref",
                 "Faoc_duu10n",
                 "Faoc_swnet",
                 # Added these for fluxEpbal call
                 "aream",
                 "mask"]

# Functions:
def sign( x, y ):
    """
    The sign transfer function SIGN(X,Y)  takes the sign
    of the second argument and puts it on the first
    argument,
    RETURNS:
    ABS(X)
    if Y >= 0
    and -ABS(X)
    if Y < 0.
    """
    if (y >= 0):
        return abs(x)
    else:
        return (-1 * abs(x))

# subroutine flux_ao( bundle_ocn, bundle_atm, fabricate )
def fluxAtmOcn(av_ocn, dom_ocn, av_atm, dom_atm, fabricate=True):
    """
    av_flux = fluxAtmOcn(av_ocn, dom_ocn, av_atm, dom_atm[, fabricate=True])

    Input:
    av_ocn,dom_ocn -> ocn state fields on ocn domain
    av_atm,dom_atm -> atm state fields on ocn domain
    fabricate -> fabricate/clobber input data(Default: True)

    Output:
    av_flux -> flux fields on ocn grid
    
    Computes Atmosphere/Ocean Fluxes

    All data must be in the ocean domain
    (NOTE:  a domain includes a particular decomposition)
    """
    nloc, rcode, m, n = -1,0,0,0
    ## First call allocations and setup:
    if(nloc < 1):
        # set nloc
        nloc = av_ocn.lsize()
        # allocate input fields:
        uocn = Numeric.zeros((nloc),Numeric.Float64)
        vocn = Numeric.zeros((nloc),Numeric.Float64)
        tocn = Numeric.zeros((nloc),Numeric.Float64)
        rmask = Numeric.zeros((nloc),Numeric.Float64)
        mask = Numeric.zeros((nloc),Numeric.Int32)
        z = Numeric.zeros((nloc),Numeric.Float64)
        u = Numeric.zeros((nloc),Numeric.Float64)
        v = Numeric.zeros((nloc),Numeric.Float64)
        ptem = Numeric.zeros((nloc),Numeric.Float64)
        shum = Numeric.zeros((nloc),Numeric.Float64)
        dens = Numeric.zeros((nloc),Numeric.Float64)
        tbot = Numeric.zeros((nloc),Numeric.Float64)
        # allocate output fields:
        sen = Numeric.zeros((nloc),Numeric.Float64)
        lat = Numeric.zeros((nloc),Numeric.Float64)
        lwup = Numeric.zeros((nloc),Numeric.Float64)
        evap = Numeric.zeros((nloc),Numeric.Float64)
        taux = Numeric.zeros((nloc),Numeric.Float64)
        tauy = Numeric.zeros((nloc),Numeric.Float64)
        tref = Numeric.zeros((nloc),Numeric.Float64)
        qref = Numeric.zeros((nloc),Numeric.Float64)
        duu10n = Numeric.zeros((nloc),Numeric.Float64)
        # get mask, which is time-invariant:
        rmask,rmask_size = dom_ocn.lgrid.exportRAttr("mask")
        mask = rmask.astype(Numeric.Int32)

    
    if(fabricate):
        ## Fabricate "reasonable" data (using dead components):
        mask = Numeric.ones(nloc,'i')
        tocn = Numeric.ones(nloc,'d') * 290.0
        uocn = Numeric.zeros(nloc,'d') 
        vocn = Numeric.zeros(nloc,'d') 
        z = Numeric.ones(nloc,'d') * 55.0
        u = Numeric.ones(nloc,'d') * 0.0
        v = Numeric.ones(nloc,'d') * 2.0
        ptem = Numeric.ones(nloc,'d') * 301.0
        shum = Numeric.ones(nloc,'d') * 1.0e-2
        dens = Numeric.ones(nloc,'d') * 1.0
        tbot = Numeric.ones(nloc,'d') * 300.0
    else:
        # # Unpack Input State Fields:
        uocn,uocn_size = av_ocn.exportRAttr('So_u')
        vocn,vocn_size = av_ocn.exportRAttr('So_v')
        tocn,tocn_size = av_ocn.exportRAttr('So_t')
        z,z_size = av_atm.exportRAttr('Sa_z')
        u,u_size = av_atm.exportRAttr('Sa_u')
        v,v_size = av_atm.exportRAttr('Sa_v')
        ptem,ptem_size = av_atm.exportRAttr('Sa_ptem')
        shum,shum_size = av_atm.exportRAttr('Sa_shum')
        dens,dens_size = av_atm.exportRAttr('Sa_dens')
        tbot,tbot_size = av_atm.exportRAttr('Sa_tbot')
    ## Call the Physics Routine:
    fields = srfflx_ao( nloc, z, u, v, ptem,
               shum, dens, tbot, uocn, vocn,
               tocn, mask, sen, lat, lwup,
               evap, taux, tauy, tref, qref,
               duu10n )
    sen, lat, lwup = fields[0], fields[1], fields[2]
    evap, taux, tauy = fields[3], fields[4], fields[5]
    tref, qref, duu10n = fields[6], fields[7], fields[8]
    ## Pack Output Flux Fields:
    av_flux = cpl.attributevector.AttributeVector()
    av_flux.initialize([],AOFLUX_FIELDS,nloc)
    av_flux.importRAttr("Faoc_sen",sen)
    av_flux.importRAttr("Faoc_lat",lat)
    av_flux.importRAttr("Faoc_lwup",lwup)
    av_flux.importRAttr("Faoc_evap",evap)
    av_flux.importRAttr("Faoc_taux",taux)
    av_flux.importRAttr("Faoc_tauy",tauy)
    av_flux.importRAttr("Faoc_tref",tref)
    av_flux.importRAttr("Faoc_qref",qref)
    av_flux.importRAttr("Faoc_duu10n",duu10n)
    return av_flux

"""def srfflx_ao(imax, z, u, v, ptem,
               shum, dens, tbot, uocn, vocn,
               tocn, mask, sen, lat, lwup,
               evap, taux, tauy, tref, qref,
               duu10n):"""
def srfflx_ao( imax, zbot, ubot, vbot, thbot,
               qbot, rbot, tbot, us, vs,
               ts, mask, sen, lat, lwup,
               evap, taux, tauy, tref, qref,
               duu10n ):
    """
    Internal atm/ocn flux calculations.

    !------------------
    ! PURPOSE:
    !   computes atm/ocn surface fluxes
    !
    ! NOTES: 
    !   o all fluxes are positive downward
    !   o net heat flux = net sw + lw up + lw down + sen + lat
    !   o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
    !   o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
    ! 
    ! ASSUMPTIONS:
    !   o Neutral 10m drag coeff: cdn
    !   o Neutral 10m stanton number: ctn = .0327 sqrt(cdn), unstable
    !                                 ctn = .0180 sqrt(cdn), stable
    !   o Neutral 10m dalton number:  cen = .0346 sqrt(cdn)
    !   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
    !------------------
    """
    umin = 0.5    # Minimum wind speed         (m/s)
    zref = 10.0   # reference height             (m)
    ztref = 2.0   # reference height for air T   (m)

    # Local Functions:
    def qsat(Tk):
        """
        the saturation humidity of air (kg/m^3)
        """
        return (640380.0 / (math.exp(5107.4 / Tk)))

    def cdn(Umps):
        """
        neutral drag coeff at 10m

        NOTE:  There were no parenthesis in the fortran code.
        Does this do the same thing?  Should we make this
        less ambiguous by including parenthesis? 
        """
        return (0.0027 / Umps + 0.000143 + 0.0000764 * Umps )

    def psimhu(xd):
        return ( math.log((1.0+xd*(2.0+xd))*(1.0+xd*xd)/10) - 2.0*math.atan(xd) + 1.571 )

    def psixhu(xd):
        return ( 2.0 * math.log((1.0 + xd*xd)/2.0) ) 
    
    al2 = math.log( zref / ztref )

    for i in range(imax):
        if( mask[i] != 0 ):
            # --- Compute some needed quantities ---
            vmag = max(umin, math.sqrt( (ubot[i]-us[i])**2 + (vbot[i] - vs[i])**2 ))
            # thvbot: virtual temp (K)
            thvbot = thbot[i] * ( 1.0 + cpl.const.zvir * qbot[i])
            # sea surf hum (kg/kg)
            ssq = 0.98 * qsat( ts[i] ) / rbot[i]
            # pot temp diff (K)
            delt = thbot[i] - ts[i]
            # spec hum dif (kg/kg)
            delq = qbot[i] - ssq
            alz = math.log(zbot[i]/zref)
            cp = cpl.const.cpdair * (1.0 + cpl.const.cpvir * ssq )
            # ---
            # First estimate of Z/L and ustar, tstar, and qstar
            # ---
            # --- Neutral coefficients, z/L = 0.0 ---
            # sign?
            stable = 0.5 + sign(0.5, delt)
            rdn = math.sqrt( cdn(vmag) )
            rhn = (1.0-stable) * 0.0327 + stable * 0.018
            ren = 0.0346
            # --- ustar, tstar, qstar ----
            ustar = rdn * vmag
            tstar = rhn * delt
            qstar = ren * delq
            # --- compute stability & evalute all stability functions
            hol = cpl.const.karman * cpl.const.g * zbot[i] * \
                  (tstar / ( (thvbot + qstar)/(1.0/cpl.const.zvir+qbot[i])))
            hol = sign(min(abs(hol),10.0),hol)
            stable = 0.5 + sign(0.5, hol)
            xsq = max(math.sqrt(abs(1.0-16.0*hol)), 1.0)
            xqq = math.sqrt(xsq)
            psimh = -5.0*hol*stable + (1.0-stable)*psimhu(xqq)
            psixh = -5.0*hol*stable + (1.0-stable)*psixhu(xqq)
            # --- shift wind speed using old coefficient ---
            rd = rdn / (1.0 + rdn/cpl.const.karman*(alz-psimh) )
            u10n = vmag * rd / rdn
            # --- update transfer coefficients at 10m and neutral stability
            rdn = math.sqrt(cdn(u10n))
            ren = 0.0346
            rhn = (1.0-stable)*0.0327 + stable * 0.018
            # --- shift all coeffs to measurement height and stability
            # \/ - Unnecessary recomputing of rd? - \/
            rd = rdn / ( 1.0 + rdn / cpl.const.karman*(alz-psimh) )
            rh = rhn / (1.0 + rhn / cpl.const.karman * (alz-psixh) )
            re = ren / (1.0 + ren / cpl.const.karman * (alz-psixh) )
            # --- update ustar, tstar, qstar using updated, shifted coeffs
            ustar = rd * vmag
            tstar = rh * delt
            qstar = re * delq
            # ---
            # iterate to converge on Z/L, ustar, tstar, qstar
            # ---
            # --- compute stability & evaluate all stability functions
            hol = cpl.const.karman * cpl.const.g * zbot[i] * \
                  (tstar / ( (thvbot + qstar)/(1.0/cpl.const.zvir+qbot[i])))
            hol = sign(min(abs(hol),10.0),hol)
            stable = 0.5 + sign(0.5, hol)
            xsq = max(math.sqrt(abs(1.0-16.0*hol)), 1.0)
            xqq = math.sqrt(xsq)
            psimh = -5.0*hol*stable + (1.0-stable)*psimhu(xqq)
            psixh = -5.0*hol*stable + (1.0-stable)*psixhu(xqq)
            # --- shift wind speed using old coefficient ---
            rd = rdn / (1.0 + rdn/cpl.const.karman*(alz-psimh) )
            u10n = vmag * rd / rdn
            # --- update transfer coefficients at 10m and neutral stability
            rdn = math.sqrt(cdn(u10n))
            ren = 0.0346
            rhn = (1.0-stable)*0.0327 + stable * 0.018
            # --- shift all coeffs to measurement height and stability
            # \/ - Unnecessary recomputing of rd? - \/
            rd = rdn / ( 1.0 + rdn / cpl.const.karman*(alz-psimh) )
            rh = rhn / (1.0 + rhn / cpl.const.karman * (alz-psixh) )
            re = ren / (1.0 + ren / cpl.const.karman * (alz-psixh) )
            # --- update ustar, tstar, qstar using updated, shifted coeffs
            ustar = rd * vmag
            tstar = rh * delt
            qstar = re * delq

            # ---
            # Compute the fluxes
            # ---
            tau = rbot[i] * ustar * ustar
            # --- momentum flux
            taux[i] = tau * ( ubot[i] - us[i] ) / vmag
            tauy[i] = tau * ( vbot[i] - vs[i] ) / vmag
            # --- heat flux ---
            sen[i] =               cp * tau * tstar / ustar
            lat[i] = cpl.const.latvap * tau * qstar / ustar
            lwup[i] = -1 * cpl.const.stebol * ts[i] ** 4
            # --- water flux ---
            evap[i] = lat[i] / cpl.const.latvap
            # ---
            # compute diagnostics: 2m ref T & Q, 10m wind speed squared
            # ---
            hol = hol * ztref / zbot[i]
            zsq = max( 1.0, math.sqrt(abs(1.0-16.0*hol)) )
            zqq = math.sqrt(xsq)
            psix2 = -5.0*hol*stable + (1.0-stable) * psixhu(xqq)
            fac   = (rh/cpl.const.karman) * (alz + al2 - psixh + psix2 )
            tref[i] = thbot[i] - delt*fac
            tref[i] = tref[i] - 0.01 * ztref # pot temp to temp correction
            fac   = (re/cpl.const.karman) * (alz + al2 - psixh + psix2 )
            qref[i] = qbot[i] - delq * fac
            duu10n[i] = u10n*u10n # 10m wind speed squared
        else:
            # --- No valid data here -- out of domain ---
            sen[i] = cpl.const.spval
            lat[i] = cpl.const.spval
            lwup[i] = cpl.const.spval
            evap[i] = cpl.const.spval
            taux[i] = cpl.const.spval
            tauy[i] = cpl.const.spval
            tref[i] = cpl.const.spval
            qref[i] = cpl.const.spval
            duu10n[i] = cpl.const.spval
    return sen, lat, lwup, evap, taux, tauy, tref, qref, duu10n




def fluxSolar(av_atm, av_ocn, av_aoflux):
    """
    Computes Ocean Net Solar
    """

    swnet,swnet_size = av_aoflux.exportRAttr("Faoc_swnet")

    anidr,anidr_size = av_ocn.exportRAttr("So_anidr")
    avsdr,avsdr_size = av_ocn.exportRAttr("So_avsdr")
    anidf,anidf_size = av_ocn.exportRAttr("So_anidf")
    avsdf,avsdf_size = av_ocn.exportRAttr("So_avsdf")

    swndr,swndr_size = av_atm.exportRAttr("Faxa_swndr")
    swvdr,swvdr_size = av_atm.exportRAttr("Faxa_swvdr")
    swndf,swndf_size = av_atm.exportRAttr("Faxa_swndf")
    swvdf,swvdf_size = av_atm.exportRAttr("Faxa_swvdf")
    
    n_o = av_aoflux.lsize()
    for n in range(n_o):
        swnet[n] = (((1.0-anidr[n])*swndr[n]) + ((1.0-avsdr[n])*swvdr[n]) +
                    ((1.0-anidf[n])*swndf[n]) + ((1.0-avsdf[n])*swvdf[n]))
    return av_aoflux

def fluxEpbal(date,av_aoflux,av_i2c,av_prec,av_r2c,av_frac,
              dom_aoflux=None, dom_i2c=None, dom_prec=None,
              dom_r2c=None, dom_frac=None):
    """
    forces evap/precip/runoff balance
    
    ! !INPUT/OUTPUT PARAMETERS:
    
    type(shr_date)  ,intent(in   ) :: date       ! current date
    type(cpl_bundle),intent(in   ) :: bun_i2c    ! ice to cpl bundle: ice evap
    type(cpl_bundle),intent(in   ) :: bun_aoflux ! a/o flux bundle  : ocn evap
    type(cpl_bundle),intent(inout) :: bun_prec   ! a/x precip bundle: i+o prec
    type(cpl_bundle),intent(inout) :: bun_r2c    ! runoff bundle    : ocn roff
    type(cpl_bundle),intent(in   ) :: bun_frac   ! fractions on ocn domain
    
    !EOP
    
    !--- local ---
    integer(IN),parameter :: k_area = 1 ! index: area of ocn domain
    integer(IN),parameter :: k_prec = 2 ! index: water ~ precipitation
    integer(IN),parameter :: k_evap = 3 ! index: water ~ evaporation
    integer(IN),parameter :: k_roff = 4 ! index: water ~ runoff
    real(R8)              :: psum(4)    ! partial/local sum of area,prec,evap,roff
    real(R8)              :: gsum(4)    !   full/global sum of area,prec,evap,roff
    real(R8)              :: tprec      ! total precip
    real(R8)              :: tevap      ! total evap
    real(R8)              :: troff      ! total runoff
    real(R8)              :: tarea      ! total area
    real(R8)              :: da         ! area of one ocn grid cell = dth*dph
    real(R8)              :: dai        ! area of ocn grid covered by ice
    real(R8)              :: dao        ! area of ocn grid covered by atm
    real(R8)              :: factor     ! prec adjustment factor: evap/prec
    integer(IN),save      :: nloc       ! size of local data array
    integer(IN),save      :: k_xrain    ! i+o-rain aVect indicies
    integer(IN),save      :: k_xsnow    ! i+o-snow aVect indicies
    integer(IN),save      :: k_ievap    ! ice-evap aVect indicies
    integer(IN),save      :: k_oevap    ! ocn-evap aVect indicies
    integer(IN),save      :: k_oroff    ! ocn-roff aVect indicies
    integer(IN),save      :: k_ifrac    ! ice-frac aVect indicies
    integer(IN),save      :: k_afrac    ! atm-frac aVect indicies
    real(R8),allocatable,save ::  area(:)  ! cell area
    real(R8),allocatable,save :: imask(:)  ! domain mask
    real(R8),allocatable      ::  mask(:)  ! domain mask
    logical    ,save :: first_call =.true. ! flags 1st invocation of routine
    integer(IN)      :: n                  ! generic loop index
    integer(IN)      :: rcode              ! return code
    integer(IN)      :: year,month,day,sec ! date & time info
    """
    global FluxEpbalFirstCall
    
    if (cpl.control.fluxepbal == 'off'):
        return None
    
    
    print "apply adjustment factor to prec & runoff over ice & ocn"
    nloc = av_aoflux.size()
    xrain,xrain_size = av_prec.exportRAttr  ("Faxc_rain")
    xsnow,xsnow_size = av_prec.exportRAttr  ("Faxc_snow")
    ievap,ievap_size = av_i2c.exportRAttr   ("Faii_evap")
    oevap,oevap_size = av_aoflux.exportRAttr("Faoc_evap")
    oroff,oroff_size = av_r2c.exportRAttr("Forr_roff")
    ifrac,ifrac_size = av_frac.exportRAttr("ifrac")
    afrac,afrac_size = av_frac.exportRAttr("afrac")
    
    imask = Numeric.zeros((nloc),Numeric.Int32)
    area,area_size = av_aoflux.exportRAttr("aream")
    mask,mask_size = av_aoflux.exportRAttr("mask")
    
    # where (mask /= 0.0) imask = 1
    # Numeric.where( condition, x, y )
    # Where condition is true, we set the element to x, if false, y.
    """
    Shouldn't we avoid things like 'mask != 0.0' comparisons?  
    """
    #    imask[:] = 0
    imask = Numeric.where( (mask != 0.0), 1, 0)
        
    # --- compute local integrals ---
    psum = {} # local/partial sum
    #initialize psum:
    for i in ["k_area","k_prec","k_evap","k_roff"]:
        psum[i] = 0.0
    for n in xrange(nloc):
        if (imask[n] != 0):
            da  = area[n]
            dai = da * ifrac[n] 
            dao = da * afrac[n]
            psum["k_area"] += dai + dao
            psum["k_prec"] += dai * (xrain[n] + xsnow[n])
            psum["k_prec"] += dao * (xrain[n] + xsnow[n])
            psum["k_evap"] += dai * ievap[n]
            psum["k_evap"] += dao * oevap[n]
            psum["k_roff"] += da * oroff[n]
   
    # --- compute factor (on master process only) ---

    gsum = {} # global sum
    #call shr_mpi_sum(psum,gsum,cpl_comm_comp,subName//" psum")
    print "psum.keys():",psum.keys()
    for key in psum.keys():
        print "psum[%s]=%s"%(key,psum[key])
        gsum[key] = mpi.reduce(psum[key],1,mpi.MPI_DOUBLE,mpi.MPI_SUM,
                               0,cpl.comm.local_comm)
        print "gsum[%s]=%s"%(key,gsum[key])
    
    if (cpl.comm.component_pid == 0):
        tarea = gsum["k_area"]
        if( tarea != 0.0 ):
            tprec = gsum["k_prec"]/tarea
            tevap = gsum["k_evap"]/tarea
            troff = gsum["k_roff"]/tarea
        else:
            tprec = 0
            tevap = 0
            troff = 0
            print "WARNING: tarea = 0.0!  av_aoflux mask is suspect!"
        factor = cpl.control.fluxepfac
        if (cpl.control.fluxepbal == 'ocn'):
            # --- use factor supplied by ocn, NOTE: may not cause E+f(P+R)=0 
            if (factor <= 0.0):
                print 'WARNING: factor from ocn = ',factor
                print 'WARNING: resetting factor to 1.0'
                factor = 1.0
                   
        elif (cpl.control.fluxepbal == 'inst'):
            # --- compute factor st f(P+R)+E=0 at every timestep 
            if (FluxEpbalFirstCall):
                print "choosing factor st P'+E+R'=0"
                
            if ( (tprec+troff) > 0.0):
                factor = -tevap/(tprec+troff)
            else:
                factor=1.0
                print 'WARNING: avg  P,R,(P+R) = ',tprec,troff,tprec+troff
                print 'WARNING: setting factor = 1.0'
        else:
            # --- invalid option ---
            print 'ERROR: unknown epbal option: ',cpl.control.fluxepbal

    # --- document factor ---
    #call shr_date_getYMD(date,year,month,day,sec)
    """
    This function doesn't do anything with the year, month, day, or seconds...
    So i'm just not going to bother.
    """
    if ( factor < 0.75 or 1.25 < factor ):
        print '<E> <P> <R> =',tevap,tprec,troff,'  f =',factor
        print 'WARNING: erroneous adjustment factor?'
    elif ( day+sec == 1 ):
        print '<E> <P> <R> =',tevap,tprec,troff,'  f =',factor
    elif ( dbug > 1 and sec == 0 ):
        print '<E> <P> <R> =',tevap,tprec,troff,'  f =',factor
    else:
        pass
    
    # --- apply factor locally on all processes
    #call shr_mpi_bcast(factor,cpl_comm_comp,subName//" factor")
    factor = mpi.bcast(factor, 1, mpi.MPI_DOUBLE, 0, cpl.comm.local_comm)
    for n in range(nloc):
        if (imask[n] != 0):
            xrain[n] = xrain[n]*factor
            xsnow[n] = xsnow[n]*factor
            oroff[n] = oroff[n]*factor
            
    av_prec.importRAttr("Faxc_rain",xrain)
    av_prec.importRAttr("Faxc_snow",xsnow)
    #av_i2c.importRAttr("Faii_evap",ievap)
    #av_aoflux.importRAttr("Faoc_evap",oevap)
    av_r2c.importRAttr("Forr_roff",oroff)
    #av_frac.importRAttr("ifrac",ifrac)
    #av_frac.importRAttr("afrac",afrac)
    
    FluxEpbalFirstCall = False
    return av_prec, av_r2c

def fluxAlbo( date, av_ocn, dom_o ):
    """
    ROUTINE: flux_albo - ocean albedo calculation
    
    DESCRIPTION:
       if flux\_albav/=0 (ie. "on" or "true") \\
          Compute four effective daily avg surface albedos for all
          combinations of visible/near-infrared and direct/diffuse radiation
          without accounting for zenith angle (ie. a "daily average" albedo) \\
       else \\
          Compute four surface albedos for all combinations of visible/
          near-infrared and direct/diffuse radiation, accounting for
          instantaneous zenith angle
    
    REMARKS:
       o upon input, albedos are assumed to be a 60 degree reference albedo
       o albedos are computed by taking the 60 deg reference albedo
         and then adjusting this value based on zenith angle
       o Albedos are independent of spectral interval and other physical
         factors such as surface wind speed.
    
       For more details see Briegleb, Bruce P., 1992: Delta-Eddington
       Approximation for Solar Radiation in the NCAR Community Climate
       Model, Journal of Geophysical Research, Vol 97, D7, pp7603-7612.
    """
    albdif, albdir = 0.06, 0.07

    n_o = av_ocn.lsize()

    ### Allocate albedo, lat & lon arrays
    anidr = Numeric.zeros((n_o),Numeric.Float64)
    anidf = Numeric.zeros((n_o),Numeric.Float64)
    avsdr = Numeric.zeros((n_o),Numeric.Float64)
    avsdf = Numeric.zeros((n_o),Numeric.Float64)

    global FluxAlboFirstCall
    if( cpl.control.fluxalbav ):
        if( FluxAlboFirstCall ):
            print "ocn albedo is zenith angle independant"
        for n in xrange(n_o):
            anidr[n] = albdir
            avsdr[n] = albdir
            anidf[n] = albdif
            avsdf[n] = albdif
    else:
        if( FluxAlboFirstCall ):
            print "ocn albedo is zenith angle dependant"
        calday = date.getJulian( cpl.control.fluxashift )
        eccen = cpl.control.orbeccen
        mvelpp = cpl.control.orbmvelpp
        lambm0 = cpl.control.orblambm0
        obliqr = cpl.control.orbobliqr
        delta,eccf = cpl.shr.orb.decl( calday, eccen, mvelpp, lambm0, obliqr )

        lons, lons_size = dom_o.lgrid.exportRAttr("lon")
        lats, lats_size = dom_o.lgrid.exportRAttr("lat")

        for n in xrange(n_o):
            rlat = cpl.const.deg2rad * lats[n]
            rlon = cpl.const.deg2rad * lons[n]
            cosz = cpl.shr.orb.cosz( calday, rlat, rlon, delta )
            if (cosz > 0.0): # then sun hit
                anidr[n] = ( (0.026 / ((cosz**1.7)+0.065))+
                             ( (0.150 * (cosz - 0.10)) *
                               (cosz - 0.50) *
                               (cosz - 1.00) ) )
                avsdr[n] = anidr[n]
                anidf[n] = albdif
                avsdf[n] = albdif
            else: # dark side of the earth
                anidr[n] = 1.0
                avsdr[n] = 1.0
                anidf[n] = 1.0
                avsdf[n] = 1.0
        
    # Clean up and store results:
    print "av_ocn lsize:",av_ocn.lsize()
    print "av_ocn fields:",av_ocn.getFields()
    av_ocn.importRAttr("So_anidr",anidr)
    av_ocn.importRAttr("So_avsdr",avsdr)
    av_ocn.importRAttr("So_anidf",anidf)
    av_ocn.importRAttr("So_avsdf",avsdf)

    # Set bundle count to 1:
    # av_ocn.reset(1)
    FluxAlboFirstCall = False
    return av_ocn


def fluxAlbi( date, av_ice, dom_i ):
    """
    ROUTINE: flux_albi - ice albedo modification
    
    DESCRIPTION:
       if flux\_albav/=0 (ie. "on" or "true") \\
         Impose a zenith angle dependance on the ice model "reference albedo".
         Currently this only involves setting albedos to zero
         on the dark side of the earth. \\
       else \\
         do not alter ice albedos
    
    REMARKS:
       o upon input, albedos are assumed to be a 60 degree reference albedo
    
    REVISION HISTORY:
       199x-       - B. Kauffman -- original cpl5 version
       2002-Oct-26 - R. Jacob -- rewritten for cpl6
    """
    albdif, albdir = 0.06, 0.07

    n_i = av_ice.lsize()

    ### Allocate albedo, lat & lon arrays
    anidr = Numeric.zeros((n_i),Numeric.Float64)
    anidf = Numeric.zeros((n_i),Numeric.Float64)
    avsdr = Numeric.zeros((n_i),Numeric.Float64)
    avsdf = Numeric.zeros((n_i),Numeric.Float64)

    global FluxAlboFirstCall
    if( cpl.control.fluxalbav ):
        if( FluxAlboFirstCall ):
            print "ice albedo is zenith angle independant"
    else:
        if( FluxAlboFirstCall ):
            print "ice albedo is zenith angle dependant"
        calday = date.getJulian( cpl.control.fluxashift )
        eccen = cpl.control.orbeccen
        mvelpp = cpl.control.orbmvelpp
        lambm0 = cpl.control.orblambm0
        obliqr = cpl.control.orbobliqr
        delta,eccf = cpl.shr.orb.decl( calday, eccen, mvelpp, lambm0, obliqr )

        lons, lons_size = dom_i.lgrid.exportRAttr("lon")
        lats, lats_size = dom_i.lgrid.exportRAttr("lat")

        for n in xrange(n_i):
            rlat = cpl.const.deg2rad * lats[n]
            rlon = cpl.const.deg2rad * lons[n]
            cosz = cpl.shr.orb.cosz( calday, rlat, rlon, delta )
            if (cosz < 0.0): # dark side of the earth
                anidr[n] = 1.0
                avsdr[n] = 1.0
                anidf[n] = 1.0
                avsdf[n] = 1.0
        
    # Clean up and store results:
    av_ice.importRAttr("Si_anidr",anidr)
    av_ice.importRAttr("Si_avsdr",avsdr)
    av_ice.importRAttr("Si_anidf",anidf)
    av_ice.importRAttr("Si_avsdf",avsdf)

    FluxAlbiFirstCall = False
    return av_ice
    


