# USAGE: python UPPER_LIMIT.py 472161854.643000 120.5 -40.2 570 3 1.5

import os, sys, math
import numpy as np
import time
from calendar import timegm
import datetime

ul_fluxes = [[0.0, 0.0, 1.4206e-07, 642.0, 1.104e-07, 72.0, 6.4838e-08, 10.0],
[105.0, 0.0, 1.4206e-07, 248.0, 1.104e-07, 23.0, 6.4838e-08, 5.0],
[105.0, 135.0, 1.4206e-07, 253.0, 1.104e-07, 29.0, 6.4838e-08, 4.0],
[105.0, 180.0, 1.4206e-07, 272.0, 1.104e-07, 29.0, 6.4838e-08, 4.0],
[105.0, 225.0, 1.4206e-07, 289.0, 1.104e-07, 38.0, 6.4838e-08, 6.0],
[105.0, 270.0, 1.4206e-07, 272.0, 1.104e-07, 25.0, 6.4838e-08, 7.0],
[105.0, 315.0, 1.4206e-07, 250.0, 1.104e-07, 25.0, 6.4838e-08, 4.0],
[105.0, 360.0, 1.4206e-07, 238.0, 1.104e-07, 26.0, 6.4838e-08, 7.0],
[105.0, 45.0, 1.4206e-07, 281.0, 1.104e-07, 32.0, 6.4838e-08, 8.0],
[105.0, 90.0, 1.4206e-07, 256.0, 1.104e-07, 31.0, 6.4838e-08, 10.0],
[110.0, 0.0, 1.4206e-07, 277.0, 1.104e-07, 26.0, 6.4838e-08, 5.0],
[120.0, 0.0, 1.4206e-07, 354.0, 1.104e-07, 47.0, 6.4838e-08, 4.0],
[120.0, 135.0, 1.4206e-07, 381.0, 1.104e-07, 46.0, 6.4838e-08, 12.0],
[120.0, 180.0, 1.4206e-07, 349.0, 1.104e-07, 47.0, 6.4838e-08, 12.0],
[120.0, 225.0, 1.4206e-07, 367.0, 1.104e-07, 35.0, 6.4838e-08, 6.0],
[120.0, 270.0, 1.4206e-07, 398.0, 1.104e-07, 32.0, 6.4838e-08, 4.0],
[120.0, 315.0, 1.4206e-07, 380.0, 1.104e-07, 38.0, 6.4838e-08, 3.0],
[120.0, 360.0, 1.4206e-07, 376.0, 1.104e-07, 33.0, 6.4838e-08, 4.0],
[120.0, 45.0, 1.4206e-07, 393.0, 1.104e-07, 34.0, 6.4838e-08, 7.0],
[120.0, 90.0, 1.4206e-07, 403.0, 1.104e-07, 36.0, 6.4838e-08, 8.0],
[130.0, 0.0, 1.4206e-07, 430.0, 1.104e-07, 55.0, 6.4838e-08, 9.0],
[135.0, 0.0, 1.4206e-07, 478.0, 1.104e-07, 64.0, 6.4838e-08, 13.0],
[135.0, 135.0, 1.4206e-07, 497.0, 1.104e-07, 49.0, 6.4838e-08, 13.0],
[135.0, 180.0, 1.4206e-07, 483.0, 1.104e-07, 54.0, 6.4838e-08, 7.0],
[135.0, 225.0, 1.4206e-07, 501.0, 1.104e-07, 49.0, 6.4838e-08, 9.0],
[135.0, 270.0, 1.4206e-07, 484.0, 1.104e-07, 49.0, 6.4838e-08, 9.0],
[135.0, 315.0, 1.4206e-07, 510.0, 1.104e-07, 51.0, 6.4838e-08, 14.0],
[135.0, 360.0, 1.4206e-07, 489.0, 1.104e-07, 45.0, 6.4838e-08, 7.0],
[135.0, 45.0, 1.4206e-07, 488.0, 1.104e-07, 48.0, 6.4838e-08, 11.0],
[135.0, 90.0, 1.4206e-07, 483.0, 1.104e-07, 48.0, 6.4838e-08, 11.0],
[150.0, 0.0, 1.4206e-07, 504.0, 1.104e-07, 54.0, 6.4838e-08, 9.0],
[150.0, 135.0, 1.4206e-07, 513.0, 1.104e-07, 54.0, 6.4838e-08, 10.0],
[150.0, 180.0, 1.4206e-07, 495.0, 1.104e-07, 44.0, 6.4838e-08, 10.0],
[150.0, 225.0, 1.4206e-07, 564.0, 1.104e-07, 66.0, 6.4838e-08, 8.0],
[150.0, 270.0, 1.4206e-07, 534.0, 1.104e-07, 66.0, 6.4838e-08, 8.0],
[150.0, 315.0, 1.4206e-07, 561.0, 1.104e-07, 51.0, 6.4838e-08, 8.0],
[150.0, 360.0, 1.4206e-07, 528.0, 1.104e-07, 48.0, 6.4838e-08, 9.0],
[150.0, 45.0, 1.4206e-07, 554.0, 1.104e-07, 48.0, 6.4838e-08, 19.0],
[150.0, 90.0, 1.4206e-07, 522.0, 1.104e-07, 58.0, 6.4838e-08, 18.0],
[15.0, 0.0, 1.4206e-07, 543.0, 1.104e-07, 69.0, 6.4838e-08, 17.0],
[15.0, 135.0, 1.4206e-07, 535.0, 1.104e-07, 69.0, 6.4838e-08, 15.0],
[15.0, 180.0, 1.4206e-07, 561.0, 1.104e-07, 61.0, 6.4838e-08, 15.0],
[15.0, 225.0, 1.4206e-07, 585.0, 1.104e-07, 71.0, 6.4838e-08, 6.0],
[15.0, 270.0, 1.4206e-07, 586.0, 1.104e-07, 70.0, 6.4838e-08, 9.0],
[15.0, 315.0, 1.4206e-07, 596.0, 1.104e-07, 60.0, 6.4838e-08, 10.0],
[15.0, 45.0, 1.4206e-07, 612.0, 1.104e-07, 60.0, 6.4838e-08, 13.0],
[15.0, 90.0, 1.4206e-07, 601.0, 1.104e-07, 59.0, 6.4838e-08, 12.0],
[165.0, 0.0, 1.4206e-07, 515.0, 1.104e-07, 45.0, 6.4838e-08, 6.0],
[165.0, 135.0, 1.4206e-07, 527.0, 1.104e-07, 44.0, 6.4838e-08, 5.0],
[165.0, 180.0, 1.4206e-07, 512.0, 1.104e-07, 37.0, 6.4838e-08, 8.0],
[165.0, 225.0, 1.4206e-07, 531.0, 1.104e-07, 37.0, 6.4838e-08, 6.0],
[165.0, 270.0, 1.4206e-07, 560.0, 1.104e-07, 42.0, 6.4838e-08, 5.0],
[165.0, 315.0, 1.4206e-07, 547.0, 1.104e-07, 50.0, 6.4838e-08, 6.0],
[165.0, 360.0, 1.4206e-07, 554.0, 1.104e-07, 51.0, 6.4838e-08, 6.0],
[165.0, 45.0, 1.4206e-07, 529.0, 1.104e-07, 43.0, 6.4838e-08, 5.0],
[165.0, 90.0, 1.4206e-07, 530.0, 1.104e-07, 47.0, 6.4838e-08, 6.0],
[180.0, 0.0, 1.4206e-07, 549.0, 1.104e-07, 50.0, 6.4838e-08, 8.0],
[30.0, 0.0, 1.4206e-07, 569.0, 1.104e-07, 62.0, 6.4838e-08, 14.0],
[30.0, 135.0, 1.4206e-07, 578.0, 1.104e-07, 65.0, 6.4838e-08, 12.0],
[30.0, 180.0, 1.4206e-07, 586.0, 1.104e-07, 63.0, 6.4838e-08, 13.0],
[30.0, 225.0, 1.4206e-07, 575.0, 1.104e-07, 52.0, 6.4838e-08, 7.0],
[30.0, 270.0, 1.4206e-07, 576.0, 1.104e-07, 61.0, 6.4838e-08, 15.0],
[30.0, 315.0, 1.4206e-07, 533.0, 1.104e-07, 61.0, 6.4838e-08, 20.0],
[30.0, 45.0, 1.4206e-07, 564.0, 1.104e-07, 61.0, 6.4838e-08, 18.0],
[30.0, 90.0, 1.4206e-07, 555.0, 1.104e-07, 54.0, 6.4838e-08, 17.0],
[45.0, 0.0, 1.4206e-07, 466.0, 1.104e-07, 41.0, 6.4838e-08, 14.0],
[45.0, 135.0, 1.4206e-07, 504.0, 1.104e-07, 67.0, 6.4838e-08, 19.0],
[45.0, 180.0, 1.4206e-07, 486.0, 1.104e-07, 50.0, 6.4838e-08, 8.0],
[45.0, 225.0, 1.4206e-07, 487.0, 1.104e-07, 51.0, 6.4838e-08, 7.0],
[45.0, 270.0, 1.4206e-07, 530.0, 1.104e-07, 51.0, 6.4838e-08, 7.0],
[45.0, 315.0, 1.4206e-07, 510.0, 1.104e-07, 55.0, 6.4838e-08, 12.0],
[45.0, 45.0, 1.4206e-07, 510.0, 1.104e-07, 55.0, 6.4838e-08, 16.0],
[45.0, 90.0, 1.4206e-07, 482.0, 1.104e-07, 57.0, 6.4838e-08, 15.0],
[60.0, 0.0, 1.4206e-07, 413.0, 1.104e-07, 62.0, 6.4838e-08, 11.0],
[60.0, 135.0, 1.4206e-07, 434.0, 1.104e-07, 75.0, 6.4838e-08, 12.0],
[60.0, 180.0, 1.4206e-07, 419.0, 1.104e-07, 66.0, 6.4838e-08, 12.0],
[60.0, 225.0, 1.4206e-07, 452.0, 1.104e-07, 49.0, 6.4838e-08, 4.0],
[60.0, 270.0, 1.4206e-07, 455.0, 1.104e-07, 52.0, 6.4838e-08, 9.0],
[60.0, 315.0, 1.4206e-07, 416.0, 1.104e-07, 52.0, 6.4838e-08, 8.0],
[60.0, 45.0, 1.4206e-07, 441.0, 1.104e-07, 54.0, 6.4838e-08, 8.0],
[60.0, 90.0, 1.4206e-07, 454.0, 1.104e-07, 55.0, 6.4838e-08, 8.0],
[75.0, 0.0, 1.4206e-07, 308.0, 1.104e-07, 32.0, 6.4838e-08, 5.0],
[75.0, 135.0, 1.4206e-07, 307.0, 1.104e-07, 37.0, 6.4838e-08, 15.0],
[75.0, 180.0, 1.4206e-07, 293.0, 1.104e-07, 36.0, 6.4838e-08, 15.0],
[75.0, 225.0, 1.4206e-07, 330.0, 1.104e-07, 30.0, 6.4838e-08, 6.0],
[75.0, 270.0, 1.4206e-07, 280.0, 1.104e-07, 32.0, 6.4838e-08, 7.0],
[75.0, 315.0, 1.4206e-07, 270.0, 1.104e-07, 34.0, 6.4838e-08, 5.0],
[75.0, 360.0, 1.4206e-07, 281.0, 1.104e-07, 46.0, 6.4838e-08, 11.0],
[75.0, 45.0, 1.4206e-07, 299.0, 1.104e-07, 46.0, 6.4838e-08, 7.0],
[75.0, 90.0, 1.4206e-07, 289.0, 1.104e-07, 40.0, 6.4838e-08, 9.0],
[90.0, 0.0, 1.4206e-07, 79.0, 1.104e-07, 10.0, 6.4838e-08, 2.0],
[90.0, 135.0, 1.4206e-07, 114.0, 1.104e-07, 12.0, 6.4838e-08, 5.0],
[90.0, 180.0, 1.4206e-07, 78.0, 1.104e-07, 9.0, 6.4838e-08, 6.0],
[90.0, 225.0, 1.4206e-07, 117.0, 1.104e-07, 11.0, 6.4838e-08, 3.0],
[90.0, 270.0, 1.4206e-07, 74.0, 1.104e-07, 11.0, 6.4838e-08, 6.0],
[90.0, 315.0, 1.4206e-07, 89.0, 1.104e-07, 15.0, 6.4838e-08, 3.0],
[90.0, 360.0, 1.4206e-07, 61.0, 1.104e-07, 12.0, 6.4838e-08, 1.0],
[90.0, 45.0, 1.4206e-07, 103.0, 1.104e-07, 20.0, 6.4838e-08, 4.0],
[90.0, 90.0, 1.4206e-07, 90.0, 1.104e-07, 12.0, 6.4838e-08, 4.0]]

TEMPO = str(sys.argv[1])
if int(len(TEMPO)) == 19:
        TTIME = timegm(time.strptime(TEMPO, '%Y-%m-%dT%H:%M:%S')) - 1072915200
else:
        TTIME = float(TEMPO)
LON_LB = float(sys.argv[2])
LAT_LB = float(sys.argv[3])
BKG = int(sys.argv[4])
N = int(sys.argv[5])
BETA = float(sys.argv[6])

if BETA != 1.0 and BETA != 1.5 and BETA != 2.0:
	print "\nOnly 1.0, 1.5, and 2.0 allowed!"
	exit(1)

THR = N*float(math.sqrt(BKG))

theta = float(os.popen("get_theta_phi.py %f %f %f" % (float(TTIME), float(LON_LB), float(LAT_LB))).read().split(' ')[0])
phi = float(os.popen("get_theta_phi.py %f %f %f" % (float(TTIME), float(LON_LB), float(LAT_LB))).read().split(' ')[1])

THETA = 0
PHI = 0

UL_ONE = []
UL_HALF = []
UL_TWO = []

if phi < 0:
	phi = phi + 360

for n in range(0,13):
	if theta >= (n*15)-7.5 and theta < (n*15)+7.5:
		THETA = n*15

for m in range(0,9):
	if phi >= (m*45)-22.5 and phi < (m*45)+22.5:
		PHI = m*45

if THETA == 0:
	THETA = 15
elif THETA == 110:
	THETA = 105
elif THETA == 130:
	THETA = 135
elif THETA == 180:
	THETA = 165
if PHI == 360:
	PHI = 0

for j in range(len(ul_fluxes)):

	if ul_fluxes[j][0] == THETA and ul_fluxes[j][1] == PHI:
		if BETA == 1.0:
			UL = [ float(LON_LB), float(LAT_LB), THETA, PHI, (ul_fluxes[j][2]/ul_fluxes[j][3])*THR]
		elif BETA == 1.5:
			UL = [float(LON_LB), float(LAT_LB), THETA, PHI, (ul_fluxes[j][4]/ul_fluxes[j][5])*THR]
		elif BETA == 2.0:
			UL = [ float(LON_LB), float(LAT_LB), THETA, PHI, (ul_fluxes[j][6]/ul_fluxes[j][7])*THR]

ul_fluxes.close()

print "\n\n%d sigma UL = %.2E erg cm^-2 at (l,b) = (%3.2f, %3.2f) = (theta,phi) = (%3.2f, %3.2f)" % (N, float(UL[4]), float(UL[0]), float(UL[1]), THETA, PHI)
print "for single power law model with photon index %2.1f\n\n" % BETA
