import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import math
import sys
import PyAstronomy
import math
#from __future__ import print_function,  division
from PyAstronomy import pyasl




# 2018 03 12
# 2018 04 11 -> corretta la definizione di SAA
# 2018 04 19 -> gestione SAA non presente

def atSidereal(mjd):


    a3 = -6.2e-6 / 3600.0;
    a2 = 0.093104 / 3600.0;
    a1 = 8640184.812866 / 3600.0 ;
    a0 = 24110.54841 / 3600.0;


    d = ((long) (mjd));
    x = (d - 51544.5) / 36525.;
    m = (mjd - d) * 24. * 60.;
    gsttod = (((a3 * x + a2) * x + a1) * x + a0)*15.0 + m * .25068447;
    while (gsttod >= 360.0): 
        gsttod = gsttod-360.0;
        #print(gsttod)
    
        #while (gsttod < 0.0 )   *gsttod += 360.0;
    #*gsttod = *gsttod * DEG2RAD;
    return gsttod;


def binary_search(a_list,  tm):
    """Performs iterative binary search to find the position of an integer in a given,  sorted,  list.
    a_list -- sorted list of integers
    item -- integer you are searching for the position of
    """

    first = 0
    last = len(a_list) - 1

    
    i = 0
    delta = 100000
    it = -1
    while first <= last:
        i = (first + last) / 2

        if (a_list[i] -tm) < delta:
            last = i - 1      
            if abs(a_list[i] -tm) < delta :
                delta = abs(a_list[i] -tm)
                it=i 
        elif (a_list[i] -tm) < delta:
            first = i + 1
            if abs(a_list[i] -tm) < delta :
                delta = abs(a_list[i] -tm)
                it=i 
    return it

##############INPUT####################################################################


if (len(sys.argv) < 4):
	print("python log_extractor.py <filename> <time> <deltaT>")
	print("Example: python log_extractor.py agql1802040218_1802040318.LOG.gz 444796813.8 2")	
	exit(-1)



filename = sys.argv[1]
timeS       = float(sys.argv[2])
deltaT       = float(sys.argv[3])


 
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################


hdulist = fits.open(filename,  ignore_missing_end=True)
tbdata = hdulist[1].data
cols = hdulist[1].columns
prihdr = hdulist[1].header



index = np.where((tbdata['time'] >= (timeS-deltaT)) & (tbdata['time'] <= (timeS+deltaT)))


timeOBT=np.array([tbdata['time']])
#SAA=np.array([tbdata['phase']])

if ((timeS < timeOBT[0][0]) | (timeS > timeOBT[0][-1])):
	print("Wrong LOG FILE")	
	print(timeS,timeOBT[0][0],timeOBT[0][-1])
	exit(-2)	


indSAA = np.where((tbdata['PHASE'] == 2) | (tbdata['PHASE'] == 1) )
print(len(indSAA))
if (len(indSAA)>1):
	saaStart=timeOBT[0][indSAA[0][0]]
	saaStop=timeOBT[0][indSAA[0][len(indSAA[0])-1]]


tbdata = tbdata[index]


timeOBT=np.array([tbdata['time']])


timeOBT=tbdata['time']



it =binary_search(timeOBT, timeS)


tbdata = tbdata[it]

position_x=tbdata['position_x']
position_y=tbdata['position_y']
position_z=tbdata['position_z']
mjd = 53005.000754444444444444 + (tbdata['time']/ 86400.0);


#print(mjd)

gsttod=atSidereal( mjd);

latitude = 90-math.degrees(math.acos(position_z/math.sqrt(position_x * position_x + position_y*position_y+ position_z*position_z)));				
longitude = math.degrees(math.atan2(position_y, position_x)) -gsttod;





if (longitude < 0):
    longitude = longitude+360;

last = hdulist[1].header['TFIELDS']-1
names=cols.names




print("TIME", round(tbdata['time'],6))


if(tbdata['mode'] == 0):
    print("MODE", "IDLE")    
if(tbdata['mode'] == 2):
    print("MODE", "OBSERVATION") 
    
    
    
if(tbdata['phase'] == 0):   
    print("PHASE", "NOMINAL")  
elif(tbdata['phase'] == 1):   
    print("PHASE", "SAA")       
elif(tbdata['phase'] == 2):   
    print("PHASE", "SAA")
elif(tbdata['phase'] == 3):   
    print("PHASE", "GROUND STATION CONTACT")       
elif(tbdata['phase'] == 4):   
    print("PHASE", "RECOVERY")       



if(tbdata['INSTR_STATUS'] == 15):
    print("AC", "ON")
    print("MCAL", "ON")
    print("SA", "ON")
    print("GRID", "ON")
if(tbdata['INSTR_STATUS'] == 11):
    print("AC", "ON")
    print("MCAL", "ON")
    print("SA", "OFF")
    print("GRID", "ON")
if(tbdata['INSTR_STATUS'] == 7):
    print("AC", "ON")
    print("MCAL", "OFF")
    print("SA", "OFF")
    print("GRID", "ON")  
if(tbdata['INSTR_STATUS'] == 3):
    print("AC", "ON")
    print("MCAL", "OFF")
    print("SA", "OFF")
    print("GRID", "ON") 
if(tbdata['INSTR_STATUS'] == 0):
    print("AC", "OFF")
    print("MCAL", "OFF")
    print("SA", "OFF")
    print("GRID", "OFF")    

	

if(len(indSAA)>1):
	if((tbdata['phase'] == 0) & (tbdata['time'] < saaStart)):
	    SAAdistance = saaStart-tbdata['time']
	    print("SAA_DISTANCE (PRE)", round(saaStart-tbdata['time']),"s")
	if((tbdata['phase'] == 0) & (tbdata['time'] > saaStart)):    
	    print(saaStop-tbdata['time'])
	    print("SAA_DISTANCE (POST)", round(saaStart-tbdata['time']),"s")
	if((tbdata['phase'] == 2)):     
	    print("SAA_DISTANCE (INTO)", round(0),"s")

else:
	 print("SAA_DISTANCE ", "NO SAA DATA","s") 

print("Q1", tbdata['Q1'])
print("Q2", tbdata['Q2'])
print("Q3", tbdata['Q3'])
print("Q4", tbdata['Q4'])
print("ATTITUDE_RA_X", round(tbdata['ATTITUDE_RA_X'],2),"degrees")
print("ATTITUDE_DEC_X",round( tbdata['ATTITUDE_DEC_X'],2),"degrees")
print("ATTITUDE_RA_Y", round(tbdata['ATTITUDE_RA_Y'],2),"degrees")
print("ATTITUDE_DEC_Y",round( tbdata['ATTITUDE_DEC_Y'],2),"degrees")
print("ATTITUDE_RA_Z", round(tbdata['ATTITUDE_RA_Z'],2),"degrees")
print("ATTITUDE_DEC_Z",round( tbdata['ATTITUDE_DEC_Z'],2),"degrees")

print("ATTITUDE_GLON_Y", round(tbdata['ATTITUDE_GLON_Y'],2),"degrees")
print("ATTITUDE_GLAT_Y", round(tbdata['ATTITUDE_GLAT_Y'],2),"degrees")

print("POSITION_X (ECI)", round(tbdata['POSITION_X'],2),"Km")
print("POSITION_Y (ECI)", round(tbdata['POSITION_Y'],2),"Km")
print("POSITION_Z (ECI)", round(tbdata['POSITION_Z'],2),"Km")


rarad=math.atan2(tbdata['POSITION_Y'], tbdata['POSITION_X'])

decrad=math.atan2(tbdata['POSITION_Z'], math.sqrt(tbdata['POSITION_X']*tbdata['POSITION_X'] + tbdata['POSITION_Y']*tbdata['POSITION_Y']))


radius=math.sqrt(tbdata['POSITION_X']*tbdata['POSITION_X']+tbdata['POSITION_Y']*tbdata['POSITION_Y']+tbdata['POSITION_Z']*tbdata['POSITION_Z'])



sat_ra =math.degrees(rarad)
if sat_ra < 0: sat_ra=sat_ra+360

print("SATELLITE_RA (ECI)", round(sat_ra,2),"degrees")
print("SATELLITE_DEC (ECI)", round(math.degrees(decrad),2),"degrees")

print("SATELLITE_RADIUS (Km)", round(radius,1),"Km")


print("EARTH_RA", round(tbdata['EARTH_RA'],2),"degrees")
print("EARTH_DEC", round(tbdata['EARTH_DEC'],2),"degrees")

#print("EARTH_THETA", tbdata['EARTH_THETA'])
#print("EARTH_PHI", tbdata['EARTH_PHI'])


print("EARTH_THETA", round(pyasl.getAngDist(tbdata['EARTH_RA'], tbdata['EARTH_DEC'], tbdata['ATTITUDE_RA_Y'], tbdata['ATTITUDE_DEC_Y']),2),"degrees")




print("RATEM_ST", tbdata['RATEM_ST'])
print("RATEM_ST_AC", tbdata['RATEM_ST_AC'])
print("RATEM_ST_MCAL_CT", tbdata['RATEM_ST_MCAL_CT'])
print("RATEM_ST_MCAL_AC_CT", tbdata['RATEM_ST_MCAL_AC_CT'])
print("RATEM_ST_MCAL", tbdata['RATEM_ST_MCAL'])
print("RATEM_ST_MCAL_AC", tbdata['RATEM_ST_MCAL_AC'])

print("LIVETIME", tbdata['LIVETIME'],"ms")


print("LATITUDE  (ECEF)",round(latitude,2),"degrees")
print("LONGITUDE (ECEF)",round(longitude,2),"degrees")


print("LOG_STATUS", 0)

hdulist.close()    










