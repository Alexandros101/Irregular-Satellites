import numpy as np
import pandas as pd
import spiceypy as spice
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.constants import c
from astropy.time import Time
from astroquery.jplhorizons import Horizons
import rebound
import random

def Furnisher(k):
    '''
    This function is used to load all kernels needed in an operation.
    Comment out kernels not in use and add the ones in use.
    
    Arguments: NA
    Returns: NA
    
    '''
    spice.kclear()
    spice.furnsh('/Users/user/Downloads/naif0008.tls.txt')
    if k == '310+341+435':
            spice.furnsh('/Users/user/Downloads/jup310.bsp')
            spice.furnsh('/Users/user/Downloads/jup341.bsp')
            spice.furnsh('/Users/user/Downloads/de435.bsp')
    elif k == '310+341':
            spice.furnsh('/Users/user/Downloads/jup310.bsp')
            spice.furnsh('/Users/user/Downloads/jup341.bsp')
    elif k == '310+435':
            spice.furnsh('/Users/user/Downloads/de435.bsp')
            spice.furnsh('/Users/user/Downloads/jup310.bsp')
    elif k == '341+435':
            spice.furnsh('/Users/user/Downloads/jup341.bsp')
            spice.furnsh('/Users/user/Downloads/de435.bsp')
    elif k == '310':
            spice.furnsh('/Users/user/Downloads/jup310.bsp')
    elif k == '341':
            spice.furnsh('/Users/user/Downloads/jup341.bsp')
    elif k == '435':
            spice.furnsh('/Users/user/Downloads/de435.bsp')
    pass
    
def get_spice_function(name,cor,loc):
    """
    This wrapper function automates the creation of objects through the JPL Horizons database. 
    
    Arguments:
    
    name: str
    
    Stipulates the target object in Horizons. The major bodies in the Solar System have an id based on their position.
    Hence '5' refers to Jupiter and '3' to Earth. A single number designator refers to a barycenter and a designator
    such as '599' to the planetary center. For minor bodies in the Solar System, the id_type in the Horizons
    must be changed to "minorbody"
    
    cor: str
    
    Refers to the type of correction that the object has. Available arguments are 'NONE', 'LT','LT+S'
    
    loc: str
    
    Designates the location of observation. Names that start with "g@#" refer to barycenters where the number designates the 
    body that the observer is based at. Hence "g@0" refers to the Solar System barycenter. Also takes Earth location designators.
    Observatories are named after their code. Hence, Pan-Starrs observatory is referred as "f51"

    Returns:
    
    get_target_xyz function
    """    
    def get_target_xyz(t):
        """
        Returns the vectors of the Horizons body at a certain time t.
        
        Arguments:
        
        t: days
        
        Julian date of observation
        
        Returns:
    
        xyz: numpy array
        
        A position vector of the observed object
    
        uvw: numpy array
        
        An instantaneous velocity vector of the observed object
        
        radec: numpy array
        
        The right ascension and declination of the observed object
        """
        
        state,lighttime = spice.spkezr(name,t,'J2000',cor,loc)
        pos,lighttime = spice.spkpos(name,t,'J2000',cor,loc)
        range,ra,dec = spice.recrad(pos) 
        xyz = np.array([state[0],state[1],state[2]])/149597870.7#6.68459e-9
        uvw = np.array([state[3],state[4],state[5]])/149597870.7*24.*3600.#*6.68459e-9
        radec = np.array([ra,dec])
        return xyz,uvw,radec*180/np.pi
    return get_target_xyz
    
# From 432 (km^3/s^2)
GMEarth = 398600.435420
GMMoon = 27068703.151194

GMMercuryB = 22031.780000
GMVenusB = 324858.592000
GMEarthB = 403503.235502
GMMarsB = 42828.375214
GMJupiterB = 126712764.133446
GMSaturnB = 37940585.200000
GMUranusB = 5794556.465752
GMNeptuneB = 6836527.100580
GMSun = 132712440041.939400 # Just Sun
GMSunterr = GMSun + GMMercuryB + GMVenusB + GMEarthB + GMMarsB

# From 310

GMJupiter = 1.266865341960128E+08
GMIo = 5.959924010272514E+03
GMEuropa = 3.202739815114734E+03
GMGanymede = 9.887819980080976E+03
GMCallisto = 7.179304867611079E+03
GMAmalthea = 1.487604677404272E-01
GM10310 = 1.327132332639000E+11

GMSaturnB310 = 3.794058500000000E+07
GMUranusB310 = 5.794548600000000E+06
GMNeptuneB310 = 6.836527100580000E+06
GMJupiterB310 = 1.267127641334463E+08

MMercuryB = GMMercuryB/GMSun
MVenusB = GMVenusB/GMSun
MEarthB = GMEarthB/GMSun
MEarth = GMEarth/GMSun
MMoon = GMMoon/GMSun
MMarsB = GMMarsB/GMSun
MJupiter = (GMJupiterB-GMIo-GMEuropa-GMGanymede-GMCallisto)/GMSun
MSaturnB = GMSaturnB/GMSun
MUranusB = GMUranusB/GMSun
MNeptuneB = GMNeptuneB/GMSun
MIo = GMIo/GMSun
MEuropa = GMEuropa/GMSun
MGanymede = GMGanymede/GMSun
MCallisto = GMCallisto/GMSun
MAmalthea = GMAmalthea/GMSun

Furnisher("310+341+435")

def get_astroquery_function(name,cor,loc):
    """
    This wrapper function automates the creation of objects through the JPL Horizons database. 
    
    Arguments:
    
    name: str
    
    Stipulates the target object in Horizons. The major bodies in the Solar System have an id based on their position.
    Hence '5' refers to Jupiter and '3' to Earth. A single number designator refers to a barycenter and a designator
    such as '599' to the planetary center. For minor bodies in the Solar System, the id_type in the Horizons
    must be changed to "minorbody"
    
    cor: str
    
    Refers to the type of correction that the object has. Available arguments are "geometric","astrometric" and 
    "apparent"
    
    loc: str
    
    Designates the location of observation. Names that start with "g@#" refer to barycenters where the number designates the 
    body that the observer is based at. Hence "g@0" refers to the Solar System barycenter. Also takes Earth location designators.
    Observatories are named after their code. Hence, Pan-Starrs observatory is referred as "f51"

    Returns:
    
    get_target_xyz function
    """    
    def get_target_xyz(t):
        """
        Returns the vectors of the Horizons body at a certain time t.
        
        Arguments:
        
        t: days
        
        Julian date of observation
        
        Returns:
    
        xyz: numpy array
        
        A position vector of the observed object
    
        uvw: numpy array
        
        An instantaneous velocity vector of the observed object
        """
        
        obj = Horizons(id = name, location = loc, epochs = t, id_type = 'majorbody')
        obj1 = obj.vectors(aberrations = cor, refplane = 'earth')
        xyz = np.array([float(obj1['x']),float(obj1['y']),float(obj1['z'])])
        uvw = np.array([float(obj1['vx']),float(obj1['vy']),float(obj1['vz'])])
        obj2 = obj.ephemerides(refsystem = 'J2000', extra_precision = True)
        radec = np.array([float(obj2['RA']),float(obj2['DEC']),float(obj1['range'])])
        return xyz,uvw,radec
    return get_target_xyz 
    
def CordConv(xyz):
    '''
    This function takes in a position vector of a body relative to an observer and returns a radec.
    
    Arguments:
    
    xyz: numpy array
    
    Should be three values, the x,y,z of the position vector
    
    Returns:
    
    radec: numpy array
    
    Two values, the first for right ascension and the second for declination
    
    '''
    DEC = -(np.arccos(xyz[2]/np.linalg.norm(xyz))-np.pi/2)
    RA = (np.arctan2(xyz[1],xyz[0]))
    while (RA > 2*np.pi):
        RA -= 2*np.pi   
    while (RA < 0):
        RA += 2*np.pi   
    return np.array([RA*180/np.pi,DEC*180/np.pi])
    
def FinderSpicef51(name,t):
    '''
    
    Finds the radec of objects in the jovian system from the pansstars observatory. Used to test the correspondence
    of spice and Horizons.
    
    Arguments:
    
    name: str
    
    The name of the object in spice documentation. See Summary_Names.txt
    
    t: float
    
    A julian date
    
    '''
    thor = Time(t, format='jd', scale='utc')
    thor = thor.tt.value
    tstr = "jd " + str(t) + " utc"
    ts = spice.str2et(tstr)
    get_earth = get_astroquery_function('0','geometric','f51')
    get_jup = get_spice_function(name,"NONE","0")
    exyz = get_earth(thor)
    jxyz = get_jup(ts)
    fxyz = LT(ts, -exyz[0], get_jup)
    vectors = fxyz[0] + exyz[0]
    return CordConv(vectors)
    
def FinderHorizonsf51(name,t):
    '''
    
    Finds the radec of objects in the jovian system from the pansstars observatory from the Horizons data base.
    
    Arguments:
    
    name: str
    
    The name of the object in spice documentation. See Summary_Names.txt
    
    t: float
    
    A julian date
    
    '''
    thor = Time(t, format='jd', scale='utc')
    thor = thor.tt.value
    get_pans = get_astroquery_function(name,'astrometric','f51')
    vecs2 = get_pans(thor)
    return CordConv(vecs2[0])
    
def SimStart2(t,bsolar,msolar,iJovian,Ms,arr):
    '''
    Starts the simulation. Adds all particles
    
    Arguments:
    
    bsolar: list of strs
    
    The names of the bodies added. Typically: '1','2','3','4', 'Jupiter','Io','Ananke'. 'Jupiter' must be added
    last followed by the satellites in the system
    
    msolar: list of floats
    
    The masses of the bodies in bsolar. Must be at the order of bsolar
    
    iJovian: list, int
    
    The indexes (+1) that the satellites in the solar system have in bsolar
    
    '''
    
    sim = rebound.Simulation()
    k = np.sqrt(Ms*10**9)*(1.49597870700*10**11)**(-1.5) *86400.00000
    sim.G = k**2
    sim.add(m=1.)
    for i in range(0,len(bsolar)):
        get_planet = get_spice_function(bsolar[i],'NONE','SUN')
        xyz,uvw,radec = get_planet(t)
        xyz = xyz
        #uvw = (uvw*24*3600)
        sim.add(m=msolar[i],x=xyz[0],y=xyz[1],z=xyz[2],vx=uvw[0],vy=uvw[1],vz=uvw[2])
    ps = sim.particles
    len(ps)
    for i in range(0,len(arr)):
        sim.add(a=arr[i][0],e=arr[i][1],inc=arr[i][2],Omega=arr[i][3],omega=arr[i][4],f=arr[i][5],primary=ps[bsolar.index('Jupiter')+1])  
    sim.n_active = len(bsolar)+1
    sim.move_to_com()
    return sim
    
def Integrator(t,ten,bsolar,msolar,iJovian,Ms,arr):
    '''
    Integrates the system. If body1 or body2 refer to the Jupiter Barycenter, while there are bodies in the Jovian system
    teh function will compute the barycenter in the simulation
    
    Arguments:
    
    bsolar
    msolar
    iJovian
    
    body1 and body2: name in str
    
    The spice names of the bodies that we want to look into to check the error of our integrator
    
     '''
    sim = SimStart2(t,bsolar,msolar,iJovian,Ms,arr)
    year = 365.25 # days
    tmax = (ten - t)/86400 
    sim.integrate(tmax)
    return sim
    
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
df = pd.read_csv('wfilter.csv')

IDs = np.array(df['ID'])
times = np.array(df['JDUTC'])
smf_name = np.array(df['smf-name'])
smf_filepath = np.array(df['smf-filepath'])

def FinderReboundf51_2(t,years,arr,times,vecs,IDs):
    '''
    
    Finds the radec of objects in the jovian system from the pansstars observatory a number of years before
    the time of observation and integrates it to the time of observation using rebound.
    
    Arguments:
    
    name: str
    
    The name of the object in spice documentation. See Summary_Names.txt
    
    t: float
    
    A julian date
    
    '''
    thor = Time(t, format='jd', scale='utc')
    thor = thor.tt.value
    tstr = "jd " + str(t) + " utc"
    ts = spice.str2et(tstr)
    get_earth = get_astroquery_function('0','geometric','f51')
    sec_year = 86400*365.25
    sim = Integrator(ts - years*sec_year,ts,['1','2','3','4','6','7','8',
        'Jupiter','Io','Europa','Ganymede','Callisto'],[MMercuryB,MVenusB,
        MEarthB,MMarsB,MSaturnB,MUranusB,MNeptuneB,MJupiter,MIo,MEuropa,MGanymede,
        MCallisto],[9,10,11,12,13,14,15],GMSun,arr)
    ttemp = ts
    radec = []
    for j in range(0,len(times)):
        temp = "jd " + str(times[j]) + " utc"
        ttemp2 = spice.str2et(temp)
        k = sim.t + (ttemp2-ttemp)/86400
        ttemp = ttemp2
        sim.integrate(k)
        thor = Time(times[j], format='jd', scale='utc')
        thor = thor.tt.value
        ps = sim.particles
        dist = np.sqrt((ps[-1].x-ps[3].x)**2 + (ps[-1].y-ps[3].y)**2 + (ps[-1].z-ps[3].z)**2)
        ltime = float(((dist*u.AU).to(u.m)/c)/u.s)
        k = sim.t - ltime/(24*3600)
        sim.integrate(k)
        ps = sim.particles
        AnSun = np.array([float(ps[-1].x),float(ps[-1].y),float(ps[-1].z)])
        #F51Sun = get_earth(thor)
        #vecs.append(F51Sun[0])
        vectors = AnSun + vecs[j][0]
        radec.append(np.array([CordConv(vectors),times[j],IDs[j]]))
        k = sim.t + ltime/(24*3600)
        sim.integrate(k)
    return radec
    
xyzs = []
get_earth = get_astroquery_function('0','geometric','f51')
for j in range(0,len(times)):
    thor = Time(times[j], format='jd', scale='utc')
    thor = thor.tt.value
    F51Sun = get_earth(thor)
    xyzs.append(F51Sun)
    
temp = pd.read_csv('seedermain.csv', header = None)
temp = np.array(temp)
lst = list(temp)
k = []
for i in range(0, len(lst)):
    f = np.array([lst[i]])
    f = list(f)
    k.append(f)
    
radecmain = [] 
arrsmain = []
for i in range(0,10000):
    radecmain.append(FinderReboundf51_2((times[0]-1),1,k[i],times,xyzs,IDs))
    
def Mags(n):
    mags  = []
    for i in range(0,n):
        mags.append(random.uniform(18,28))   
    return mags
    
ID = []
SMFs = []
SMF_FLs = []
RAs = []
DECs = []
MAG = []
for i in range(0,743):
    for j in range(0,10000):
        ID.append(IDs[i])
        SMFs.append(smf_name[i])
        SMF_FLs.append(smf_filepath[i])
        RAs.append(radecmain[j][i][0][0])
        DECs.append(radecmain[j][i][0][1])
        MAG.append(k[j][-1][-1])
        
d = {'ID': ID, 'smf-name': SMFs, 'smf-filepath': SMF_FLs, 'RA': RAs, 'DEC': DECs, 'Mag': MAG}
df = pd.DataFrame(data=d)

df.to_csv('Satellites.csv', header=True, index=False)
