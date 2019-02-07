# NEW LCPLOT.py

# python NEW_LCPLOT.py 10553 168999780.444814 -30 30 0.032 -10 -5 0 168999779.000000 GRB090510 SOLO

import os, sys, math, os.path
from scipy.stats import poisson
import numpy as np
import matplotlib.pyplot as plt
import itertools
from itertools import groupby
import ROOT
from ROOT import TFile
from ROOT import gDirectory
import matplotlib.ticker as ticker
import time
from calendar import timegm
import datetime

PATH_DATA = str(os.environ['PATH_DATA'])
MCALPIPE = str(os.environ['MCALPIPE'])
AGILEPIPE = str(os.environ['AGILEPIPE'])
MCALSW = str(os.environ['MCALSW'])

CONT = int(sys.argv[1])
TEMPO = str(sys.argv[2])
if int(len(TEMPO)) == 19:
        TTIME = timegm(time.strptime(TEMPO, '%Y-%m-%dT%H:%M:%S')) - 1072915200
else:
        TTIME = float(TEMPO)
SX = float(sys.argv[3])
DX = float(sys.argv[4])
BIN = float(sys.argv[5])
BACK_START = float(sys.argv[6])
BACK_STOP = float(sys.argv[7])
SHIFT = int(sys.argv[8])
EXT_TEMPO = str(sys.argv[9])
if int(len(EXT_TEMPO)) == 19:
        EXT_TIME = timegm(time.strptime(EXT_TEMPO, '%Y-%m-%dT%H:%M:%S')) - 1072915200
else:
        EXT_TIME = float(EXT_TEMPO)
NAME = str(sys.argv[10])
FLAG = str(sys.argv[11])

OUTPATH_ROOT = "%s/%s" % (str(os.getcwd()), NAME)

os.system("mkdir -p %s" % OUTPATH_ROOT)

FITSPATH = "%s/DATA_2/COR" % PATH_DATA

os.chdir(OUTPATH_ROOT)

os.system("mkdir -p %s/RT" % OUTPATH_ROOT)

#####################
# CONVERT ROOT FILE #
#####################

print "\nFITS FILE IN:   %s/PKP%06.d_1_3908_000.lv1.cor.gz" % (FITSPATH, CONT)

if (os.path.isfile("%s/PKP%06.d_1_3908_000.lv1.cor.gz" % (FITSPATH, CONT)) == False) or (os.path.isfile("%s/PKP%06.d_1_3916_000.lv1.cor.gz" % (FITSPATH, CONT)) == False):
	print "\nNo file\n"

INPATH_3908_FITS = "%s/PKP%06.d_1_3908_000.lv1.cor.gz" % (FITSPATH, CONT)

INPATH_3916_FITS = "%s/PKP%06.d_1_3916_000.lv1.cor.gz" % (FITSPATH, CONT)

####################################################
# CONVERSION OF 3908/3916 FITS FILES TO ROOT FILES #
####################################################

if (os.path.isfile("%s/RT%06.d_3908.root" % (OUTPATH_ROOT, CONT)) == False):

	# creation of 3908 root file (no verbose)

	os.system("mcalanalyzer %s >& /dev/null" % INPATH_3908_FITS)

	# move the trigger list file into the OUTPUTH directory

	os.system("mv grboutput.txt %s/%06.d_triggers.txt" % (OUTPATH_ROOT, CONT))

	# creation of 3916 root file (already in the OUTPATH directory) (no verbose)

	os.system("%s/bin/convPKP3916toRoot %s %s/RT/RT%06.d_3916.root 1 >& /dev/null" % (MCALSW, INPATH_3916_FITS, OUTPATH_ROOT, CONT))

else:
	print "\nFILE GIA' PRESENTI"

print "\nROOT FILE IN:   %s/RT/%06.d/RT%06.d_3908.root" % (OUTPATH_ROOT, CONT, CONT)

##############################################################
# READ TRIGGER LIST FILE AND STORE T_START AND T_STOP VALUES #
##############################################################

TRG_T_START = []
TRG_T_STOP = []
TRG_PRE_BURST = []
TRG_POST_BURST = []
T_START = []
T_STOP = []

L_SUBMS = ''
L_1MS = ''
L_16MS = ''
L_64MS = ''
L_256MS = ''
L_1024MS = ''
L_8192MS = ''

ftrg = open("%s/%06.d_triggers.txt" % (OUTPATH_ROOT, CONT), "r")

M = 0

for line in ftrg:
	line = line.strip()
	columns = line.split()
	if int(columns[0]) == CONT:

		if float(columns[2])-5 <= TTIME <= float(columns[3])+10:
			M = int(columns[1])-1

		TRG_T_START.append(float(columns[2]))
		TRG_T_STOP.append(float(columns[3]))
		TRG_PRE_BURST.append(float(columns[2])+float(columns[7]))
		TRG_POST_BURST.append(float(columns[2])+float(columns[8]))
		if int(columns[12]) == 1:
			L_SUBMS = 'sub-ms'
		if int(columns[13]) == 1:
			L_1MS = '1ms'
		if int(columns[14]) == 1:
			L_16MS = '16ms'
		if int(columns[17]) == 2:
			L_64MS = '64ms'
		if int(columns[17]) == 3:
			L_256MS = '256ms'
		if int(columns[17]) == 4:
			L_1024MS = '1024ms'
		if int(columns[17]) == 5:
			L_8192MS = '8192ms'

print "\ndata are in trigger ", M+1, "between ", TRG_PRE_BURST[M], TRG_POST_BURST[M]

GRB_ROOTNAME  = ROOT.TFile.Open("%s/RT/RT%06.d_3908.root" % (OUTPATH_ROOT, CONT))
GRB_TREENAME  = GRB_ROOTNAME.Get("tdata")
GRB_CHAIN     = gDirectory.Get("tdata")
GRB_N_ENTRIES = GRB_CHAIN.GetEntriesFast()

EVTt = [] # <-- aggiunto per eliminare il ricircolo
EVT = []
EVT_RED = []

for atrg in GRB_ROOTNAME.tdata:
	EVTt.append([float(atrg.time), float(atrg.Etot), int(atrg.mult)])

GRB_ROOTNAME.Close()

# nel caso di un evento a cavallo tra due contatti (plotta due contatti insieme)

##
#GRB_ROOTNAME2  = ROOT.TFile.Open("%s/ARCHIVE/%06.d/RT%06.d_3908.root" % (OUTPATH_ROOT, CONT+2, CONT+2))
#GRB_TREENAME2  = GRB_ROOTNAME2.Get("tdata")
#GRB_CHAIN2     = gDirectory.Get("tdata")
#GRB_N_ENTRIES2 = GRB_CHAIN2.GetEntriesFast()
#for atrg in GRB_ROOTNAME2.tdata:
#	EVTt.append([float(atrg.time)-CORR, float(atrg.Etot), int(atrg.mult)])
#GRB_ROOTNAME2.Close()
##

EVT = dict((x[0], x) for x in EVTt).values() # <-- aggiunto per eliminare il ricircolo

EVT = sorted(EVT) # <-- aggiunto per eliminare il ricircolo

# merge events within 4 microseconds and introducing a last column to provide a flag for each trigger acquisition

q = 0
h = 0
trg = 0

while q < int(len(EVT)):
	if EVT[q][0]-EVT[q-1][0] <= 0.000004:
		EVT_RED.append([float(EVT[q-1][0]), float(EVT[q-1][1])+float(EVT[q][1]), int(EVT[q-1][2])+int(EVT[q][2])])
		q += 1
	else:
		if EVT[q][0] > EVT[q-1][0]:
			EVT_RED.append([float(EVT[q-1][0]), float(EVT[q-1][1]), int(EVT[q-1][2])])
	q += 1

# information on time

UTC_MCAL = str(time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(1072915200+TTIME)))
UTC_EXT = str(time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(1072915200+EXT_TIME)))

N_TRG = int(len(TRG_PRE_BURST))

if SX != 0.0 and DX != 0.0:
	TRG_PRE_BURST[M] = TTIME + SX
	TRG_POST_BURST[M] = TTIME + DX

N_BIN = int(((TRG_POST_BURST[M]+2.0-TRG_PRE_BURST[M])/BIN))

BS = np.arange((float)(TRG_PRE_BURST[M]-1.0)+(SHIFT*BIN/4.0), (float)(TRG_POST_BURST[M]+1.0)+(SHIFT*BIN/4.0)+2, BIN)
LIGHT = np.zeros((N_BIN, 2))
BACK = []
hist, bin_edges = np.histogram(EVT_RED, bins=BS)
lc = open("%s_lc.txt" % NAME, "w")
for c in range(N_BIN):
	LIGHT[c]=[BS[c+1], hist[c]]
	lc.write("\t%3.6f\t%9.6f\t%3.5f\n" % (BS[c+1]-TTIME, BS[c+1], hist[c]))
	if BACK_START != 0 and BACK_STOP != 0:
		if BACK_START <= BS[c+1]-TTIME <= BACK_STOP:
			BACK.append([BS[c+1], hist[c]])
lc.close()

if BACK_START != 0 and BACK_STOP != 0:
	BACK = np.array(BACK)
	BKG  = np.mean(BACK[ np.where( (BACK[:,1] != 0 ) ), 1 ])
else:
	BKG = 0.0

THR_15SIGMA = BKG + 15*math.sqrt(BKG)
THR_7SIGMA = BKG + 7*math.sqrt(BKG)
THR_6SIGMA = BKG + 6*math.sqrt(BKG)
THR_5SIGMA = BKG + 5*math.sqrt(BKG)
THR_4SIGMA = BKG + 4*math.sqrt(BKG)
THR_3SIGMA = BKG + 3*math.sqrt(BKG)
THR_2SIGMA = BKG + 2*math.sqrt(BKG)
THR_1SIGMA = BKG + 1*math.sqrt(BKG)

print "\nMCAL times: %s = %9.6f" % (UTC_MCAL, TTIME)
print "EXT times: %s = %9.6f" % (UTC_EXT, EXT_TIME)
print "\nBackground: %3.2f cts / %d ms  =  %3.2f cts/s\n" % (BKG, BIN*1000., float(BKG/BIN))

if FLAG == 'SOLO':

	LIGHTCURVE = np.zeros((N_BIN, 2))

	for c in range(0,N_BIN):
		LIGHTCURVE[c]=[BS[c+1], hist[c]]

	# plot GRB light curve in all energy bands and counts

	fig = plt.figure(figsize=(18,6))
	ax = fig.add_subplot(1,1,1)
	#ax.set_title('%s MCAL %s - IPN %s' % (NAME, UTC_MCAL, UTC_EXT), fontsize=22)
	ax.set_title('%s MCAL %s' % (NAME, UTC_MCAL), fontsize=22)
	ax.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,1], markersize=5, c='black', linewidth=1)
	ax.set_xlim(SX,DX)
	ax.set_ylim(0,max(hist)+(max(hist)/3.0)) # modifica dell'altezza massima +1/3
	plt.axvline(x=0, linewidth=7, color = 'magenta', alpha = 0.3)
	plt.axvline(x=-180, linewidth=4, linestyle='dotted', color = 'magenta', alpha = 0.3)
	plt.axvline(x=180, linewidth=4, linestyle='dotted', color = 'magenta', alpha = 0.3)
	plt.axvline(x=-120, linewidth=4, linestyle='dotted', color = 'magenta', alpha = 0.3)
	plt.axvline(x=120, linewidth=4, linestyle='dotted', color = 'magenta', alpha = 0.3)
	plt.axvline(x=-60, linewidth=4, linestyle='dotted', color = 'magenta', alpha = 0.3)
	plt.axvline(x=60, linewidth=4, linestyle='dotted', color = 'magenta', alpha = 0.3)
	if EXT_TIME != 0.0:
		plt.axvline(x=EXT_TIME-TTIME, linewidth=7, color = 'green', alpha = 0.3)
	if BKG != 0.0:
		plt.axhline(y=BKG, linewidth=2, color = 'red', alpha = 0.5, label='BKG = %3.2f cts/s' % float(BKG/BIN))
		plt.axhline(y=THR_5SIGMA, linewidth=2, color = 'green', linestyle = 'dashed', alpha = 0.4, label='$5\sigma$')
	#for J in range(0,len(TRG_T_START)):
		#plt.axvline(x=TRG_T_START[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		#plt.axvline(x=TRG_T_STOP[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		#plt.axvline(x=TRG_PRE_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
		#plt.axvline(x=TRG_POST_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
	ax.set_ylabel('MCAL High Res\n[counts / %.4f s]' % BIN, fontsize=20)
	plt.setp(ax.get_xticklabels(), fontsize=18)
	plt.setp(ax.get_yticklabels(), fontsize=18)
	#ax.set_xlabel('t - UT %s [s]' % UTC_MCAL, fontsize=20)
	ax.set_xlabel('t - %9.6f [s]' % TTIME, fontsize=20)
	plt.legend(loc='upper right', fontsize=20)
	plt.savefig('%s_%06.d_%9.6f.png' % (NAME, CONT, TTIME))
	plt.close(fig)

	os.system("rm -r %s/RT" % OUTPATH_ROOT)

if FLAG == 'MT':

        LIGHTCURVE = np.zeros((N_BIN, 4))

	EVT_RED_LE = [EVT_RED[x] for x in range(len(EVT_RED)) if EVT_RED[x][1] < 1.4]
	EVT_RED_HE = [EVT_RED[x] for x in range(len(EVT_RED)) if EVT_RED[x][1] >= 1.4]

	hist_LE, bin_edges = np.histogram(EVT_RED_LE, bins=BS)
	hist_HE, bin_edges = np.histogram(EVT_RED_HE, bins=BS)

	for c in range(0,N_BIN):
		LIGHTCURVE[c]=[BS[c+1], hist[c], hist_LE[c], hist_HE[c]]

        fig = plt.figure(figsize=(16,10))
	ax = fig.add_subplot(311)
	ax.set_title('Contact %06.d   Trigger %d\nBin %3.3f s -  T0 %f\n' % (CONT, M+1, BIN, TTIME), fontsize=16)
	ax.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,1], markersize=5, c='black', linewidth=1)
	ax.set_xlim(TRG_PRE_BURST[M]-TTIME-1.0,TRG_POST_BURST[M]-TTIME+1.0)
	ax.set_ylim(0,max(hist)+(max(hist)/3.0)) # modifica dell'altezza massima +1/3
	plt.axvline(x=0, linewidth=1, color = 'black', linestyle = 'dashed')
	ax.set_ylabel('[counts / %3.3f s]\n0.4 - 100 MeV\n' % BIN, fontsize=16)
	plt.setp(ax.get_xticklabels(), fontsize=16)
	plt.setp(ax.get_yticklabels(), fontsize=16)

	ax_LE = fig.add_subplot(312)
	ax_LE.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,2], markersize=5, c='black', linewidth=1)
	ax_LE.set_xlim(TRG_PRE_BURST[M]-TTIME-1.0,TRG_POST_BURST[M]-TTIME+1.0)
	plt.axvline(x=0, linewidth=1, color = 'black', linestyle = 'dashed')
	ax_LE.set_ylim(0,max(hist_LE)+(max(hist_LE)/3.0)) # modifica dell'altezza massima +1/3
	ax_LE.set_ylabel('[counts / %3.3f s]\n< 1.4 MeV\n' % BIN, fontsize=16)
	plt.setp(ax_LE.get_xticklabels(), fontsize=16)
	plt.setp(ax_LE.get_yticklabels(), fontsize=16)

	ax_HE = fig.add_subplot(313)
	ax_HE.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,3], markersize=5, c='black', linewidth=1)
	ax_HE.set_xlim(TRG_PRE_BURST[M]-TTIME-1.0,TRG_POST_BURST[M]-TTIME+1.0)
	plt.axvline(x=0, linewidth=1, color = 'black', linestyle = 'dashed')
	ax_HE.set_ylim(0,max(hist_HE)+(max(hist_HE)/3.0)) # modifica dell'altezza massima +1/3
	ax_HE.set_ylabel('[counts / %3.3f s]\n>= 1.4 MeV\n' % BIN, fontsize=16)
	plt.setp(ax_HE.get_xticklabels(), fontsize=16)
	plt.setp(ax_HE.get_yticklabels(), fontsize=16)

	ax_HE.set_xlabel('t - %s [s]' % UTC_MCAL, fontsize=16)
        plt.savefig('%s_%06.d_%9.6f.png' % (NAME, CONT, TTIME))
	plt.close(fig)

if FLAG == 'KW':

	# separate events with different energies and molteplicities PER KONUS-WIND

	EVT_RED_E1 = [EVT_RED[x] for x in range(len(EVT_RED)) if 0.02 < EVT_RED[x][1] <= 0.08]
	EVT_RED_E2 = [EVT_RED[x] for x in range(len(EVT_RED)) if 0.08 < EVT_RED[x][1] <= 0.31]
	EVT_RED_E3 = [EVT_RED[x] for x in range(len(EVT_RED)) if 0.31 < EVT_RED[x][1] <= 1.20]
	EVT_RED_E4 = [EVT_RED[x] for x in range(len(EVT_RED)) if 0.39 < EVT_RED[x][1] <= 1.14]
	EVT_RED_E5 = [EVT_RED[x] for x in range(len(EVT_RED)) if 1.14 < EVT_RED[x][1] <= 2.90]
	EVT_RED_E6 = [EVT_RED[x] for x in range(len(EVT_RED)) if 2.90 < EVT_RED[x][1] <= 7.20]
	EVT_RED_E7 = [EVT_RED[x] for x in range(len(EVT_RED)) if 7.20 < EVT_RED[x][1] <= 14.8]
	EVT_RED_E8 = [EVT_RED[x] for x in range(len(EVT_RED)) if EVT_RED[x][1] > 14.8]

	hist_E1, bin_edges = np.histogram(EVT_RED_E1, bins=BS)
	hist_E2, bin_edges = np.histogram(EVT_RED_E2, bins=BS)
	hist_E3, bin_edges = np.histogram(EVT_RED_E3, bins=BS)
	hist_E4, bin_edges = np.histogram(EVT_RED_E4, bins=BS)
	hist_E5, bin_edges = np.histogram(EVT_RED_E5, bins=BS)
	hist_E6, bin_edges = np.histogram(EVT_RED_E6, bins=BS)
	hist_E7, bin_edges = np.histogram(EVT_RED_E7, bins=BS)
	hist_E8, bin_edges = np.histogram(EVT_RED_E8, bins=BS)

	LIGHTCURVE = np.zeros((N_BIN, 10))

	for c in range(0,N_BIN):
		LIGHTCURVE[c]=[BS[c+1], hist[c], hist_E1[c], hist_E2[c], hist_E3[c], hist_E4[c], hist_E5[c], hist_E6[c], hist_E7[c], hist_E8[c]]

	# plot GRB light curve in all energy bands and counts

	fig = plt.figure(figsize=(12,18))

	ax = fig.add_subplot(9,1,1)
	ax.set_title('%s MCAL %s\n' % (NAME, UTC_MCAL), fontsize=22)
	ax.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,1], markersize=5, c='black', linewidth=1, label='0.3-100 MeV')
	ax.set_xlim(0,20)
	ax.set_ylabel('[counts / %3.3f s]\n' % BIN, fontsize=12)
	plt.legend(loc='upper right', fontsize=20)

	ax_E1 = fig.add_subplot(9,1,2)
	ax_E1.step(LIGHTCURVE[:,0]-TTIME-10000, LIGHTCURVE[:,2], markersize=5, c='blue', linewidth=1, label='20-80 keV')
	ax_E1.set_xlim(0,20)
	ax_E1.set_ylabel('[counts / %3.3f s]\n' % BIN, fontsize=12)
	plt.legend(loc='upper right', fontsize=20)

	ax_E2 = fig.add_subplot(9,1,3)
	ax_E2.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,3], markersize=5, c='blue', linewidth=1, label='80-310 keV')
	ax_E2.set_xlim(0,20)
	ax_E2.set_ylabel('[counts / %3.3f s]\n' % BIN, fontsize=12)
	plt.legend(loc='upper right', fontsize=20)

	ax_E3 = fig.add_subplot(9,1,4)
	ax_E3.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,4], markersize=5, c='blue', linewidth=1, label='310-1200 keV')
	ax_E3.set_xlim(0,20)
	ax_E3.set_ylabel('[counts / %3.3f s]\n' % BIN, fontsize=12)
	plt.legend(loc='upper right', fontsize=20)

	ax_E4 = fig.add_subplot(9,1,5)
	ax_E4.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,5], markersize=5, c='blue', linewidth=1, label='390-1140 keV')
	ax_E4.set_xlim(0,20)
	ax_E4.set_ylabel('[counts / %3.3f s]\n' % BIN, fontsize=12)
	plt.legend(loc='upper right', fontsize=20)

	ax_E5 = fig.add_subplot(9,1,6)
	ax_E5.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,6], markersize=5, c='blue', linewidth=1, label='1.0-2.9 MeV')
	ax_E5.set_xlim(0,20)
	ax_E5.set_ylabel('[counts / %3.3f s]\n' % BIN, fontsize=12)
	plt.legend(loc='upper right', fontsize=20)

	ax_E6 = fig.add_subplot(9,1,7)
	ax_E6.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,7], markersize=5, c='blue', linewidth=1, label='2.9-7.2 MeV')
	ax_E6.set_xlim(0,20)
	ax_E6.set_ylabel('[counts / %3.3f s]\n' % BIN, fontsize=12)
	plt.legend(loc='upper right', fontsize=20)

	ax_E7 = fig.add_subplot(9,1,8)
	ax_E7.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,8], markersize=5, c='blue', linewidth=1, label='7.2-14.8 MeV')
	ax_E7.set_xlim(0,20)
	ax_E7.set_ylabel('[counts / %3.3f s]\n' % BIN, fontsize=12)
	plt.legend(loc='upper right', fontsize=20)

	ax_E8 = fig.add_subplot(9,1,9)
	ax_E8.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,9], markersize=5, c='blue', linewidth=1, label='14.8-100 MeV')
	ax_E8.set_xlim(0,20)
	ax_E8.set_ylabel('[counts / %3.3f s]\n' % BIN, fontsize=12)
	plt.legend(loc='upper right', fontsize=20)

	plt.setp(ax_E8.get_xticklabels(), fontsize=16)
	ax_E8.set_xlabel('t - %9.6f [s]' % TTIME, fontsize=18)
	plt.savefig('%s_%06.d_%9.6f.png' % (NAME, CONT, TTIME))
	plt.close(fig)

if FLAG == 'AGILE':

	# separate events with different energies and molteplicities PER AGILE

	EVT_RED_E1 = [EVT_RED[x] for x in range(len(EVT_RED)) if EVT_RED[x][1] <= 0.7]
	EVT_RED_E2 = [EVT_RED[x] for x in range(len(EVT_RED)) if 0.7 < EVT_RED[x][1] <= 1.4]
	EVT_RED_E3 = [EVT_RED[x] for x in range(len(EVT_RED)) if 1.4 < EVT_RED[x][1] <= 2.8]
	EVT_RED_E4 = [EVT_RED[x] for x in range(len(EVT_RED)) if EVT_RED[x][1] > 2.8]
	EVT_RED_LE = [EVT_RED[x] for x in range(len(EVT_RED)) if EVT_RED[x][1] < 1.4]
	EVT_RED_HE = [EVT_RED[x] for x in range(len(EVT_RED)) if EVT_RED[x][1] >= 1.4]
	EVT_RED_RE = [EVT_RED[x] for x in range(len(EVT_RED)) if (EVT_RED[x][1] < 100.0 and EVT_RED[x][2] < 4)]
	EVT_RED_FE = [EVT_RED[x] for x in range(len(EVT_RED)) if (EVT_RED[x][1] >= 100.0 or EVT_RED[x][2] >= 4)]
	EVT_RED_MULT4 = [EVT_RED[x] for x in range(len(EVT_RED)) if EVT_RED[x][2] >= 4]
	EVT_RED_100MEV = [EVT_RED[x] for x in range(len(EVT_RED)) if EVT_RED[x][1] >= 100.0]

	hist_E1, bin_edges = np.histogram(EVT_RED_E1, bins=BS)
	hist_E2, bin_edges = np.histogram(EVT_RED_E2, bins=BS)
	hist_E3, bin_edges = np.histogram(EVT_RED_E3, bins=BS)
	hist_E4, bin_edges = np.histogram(EVT_RED_E4, bins=BS)
	hist_LE, bin_edges = np.histogram(EVT_RED_LE, bins=BS)
	hist_HE, bin_edges = np.histogram(EVT_RED_HE, bins=BS)
	hist_RE, bin_edges = np.histogram(EVT_RED_RE, bins=BS)
	hist_FE, bin_edges = np.histogram(EVT_RED_FE, bins=BS)
	hist_MULT4, bin_edges = np.histogram(EVT_RED_MULT4, bins=BS)
	hist_100MEV, bin_edges = np.histogram(EVT_RED_100MEV, bins=BS)

	LIGHTCURVE = np.zeros((N_BIN, 13))

	for c in range(0,N_BIN):
		LIGHTCURVE[c]=[BS[c+1], hist[c], hist_E1[c], hist_E2[c], hist_E3[c], hist_E4[c], hist_LE[c], hist_HE[c], hist_RE[c], hist_FE[c], float(int(hist_FE[c])/(int(hist[c]+1))), hist_MULT4[c], hist_100MEV[c]]

	# plot GRB light curve in all energy bands and counts

	fig = plt.figure(figsize=(12,18))
	ax = fig.add_subplot(12,1,1)
	ax.set_title('Contact %06.d   Trigger %d\nBin %3.3f s -  T0 %f\n' % (CONT, M+1, BIN, TTIME), fontsize=16)
	ax.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,1], markersize=5, c='black', linewidth=1)
	ax.set_xlim(TRG_PRE_BURST[M]-TTIME,TRG_POST_BURST[M]-TTIME)
	ax.set_ylim(0,max(hist)+(max(hist)/3.0)) # modifica dell'altezza massima +1/3
	plt.axvline(x=0, linewidth=5, color = 'magenta', alpha = 0.3)
	for J in range(0,len(TRG_T_START)):
		plt.axvline(x=TRG_T_START[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_T_STOP[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_PRE_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
		plt.axvline(x=TRG_POST_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
	ax.axhline(y=BKG, linewidth=1, c='red', alpha=0.5)
	ax.axvline(x=TTIME, linewidth=12, color = 'magenta', alpha = 0.2)
	plt.axhline(y=THR_6SIGMA, linewidth=2, color = 'black', linestyle = 'dashed', alpha = 0.4)
	plt.axhline(y=THR_5SIGMA, linewidth=2, color = 'red', linestyle = 'dashed', alpha = 0.4)
	plt.axhline(y=THR_4SIGMA, linewidth=2, color = 'orange', linestyle = 'dashed', alpha = 0.4)
	plt.axhline(y=THR_3SIGMA, linewidth=2, color = 'green', linestyle = 'dashed', alpha = 0.4)
	ax.set_ylabel('[counts / %3.3f s]\n0.4 - 100 MeV\n' % BIN, fontsize=12)

	ax_E1 = fig.add_subplot(12,1,2)
	ax_E1.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,2], markersize=5, c='blue', linewidth=1)
	ax_E1.set_xlim(TRG_PRE_BURST[M]-TTIME,TRG_POST_BURST[M]-TTIME)
	plt.axvline(x=0, linewidth=5, color = 'magenta', alpha = 0.3)
	for J in range(0,len(TRG_T_START)):
		plt.axvline(x=TRG_T_START[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_T_STOP[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_PRE_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
		plt.axvline(x=TRG_POST_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
	ax_E1.set_ylim(0,max(hist)+(max(hist)/3.0)) # modifica dell'altezza massima +1/3
	ax_E1.set_ylabel('[counts / %3.3f s]\n< 0.7 MeV\n' % BIN, fontsize=12)

	ax_E2 = fig.add_subplot(12,1,3)
	ax_E2.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,3], markersize=5, c='blue', linewidth=1)
	ax_E2.set_xlim(TRG_PRE_BURST[M]-TTIME,TRG_POST_BURST[M]-TTIME)
	plt.axvline(x=0, linewidth=5, color = 'magenta', alpha = 0.3)
	for J in range(0,len(TRG_T_START)):
		plt.axvline(x=TRG_T_START[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_T_STOP[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_PRE_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
		plt.axvline(x=TRG_POST_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
	ax_E2.set_ylim(0,max(hist)+(max(hist)/3.0)) # modifica dell'altezza massima +1/3
	ax_E2.set_ylabel('[counts / %3.3f s]\n0.7 - 1.4 MeV\n' % BIN, fontsize=12)

	ax_E3 = fig.add_subplot(12,1,4)
	ax_E3.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,4], markersize=5, c='blue', linewidth=1)
	ax_E3.set_xlim(TRG_PRE_BURST[M]-TTIME,TRG_POST_BURST[M]-TTIME)
	plt.axvline(x=0, linewidth=5, color = 'magenta', alpha = 0.3)
	for J in range(0,len(TRG_T_START)):
		plt.axvline(x=TRG_T_START[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_T_STOP[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_PRE_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
		plt.axvline(x=TRG_POST_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
	ax_E3.set_ylim(0,max(hist)+(max(hist)/3.0)) # modifica dell'altezza massima +1/3
	ax_E3.set_ylabel('[counts / %3.3f s]\n1.4 - 2.8\n' % BIN, fontsize=12)

	ax_E4 = fig.add_subplot(12,1,5)
	ax_E4.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,5], markersize=5, c='blue', linewidth=1)
	ax_E4.set_xlim(TRG_PRE_BURST[M]-TTIME,TRG_POST_BURST[M]-TTIME)
	plt.axvline(x=0, linewidth=5, color = 'magenta', alpha = 0.3)
	for J in range(0,len(TRG_T_START)):
		plt.axvline(x=TRG_T_START[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_T_STOP[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_PRE_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
		plt.axvline(x=TRG_POST_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
	ax_E4.set_ylim(0,max(hist)+(max(hist)/3.0)) # modifica dell'altezza massima +1/3
	ax_E4.set_ylabel('[counts / %3.3f s]\n> 2.8 MeV\n' % BIN, fontsize=12)

	ax_LE = fig.add_subplot(12,1,6)
	ax_LE.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,6], markersize=5, c='blue', linewidth=1)
	ax_LE.set_xlim(TRG_PRE_BURST[M]-TTIME,TRG_POST_BURST[M]-TTIME)
	plt.axvline(x=0, linewidth=5, color = 'magenta', alpha = 0.3)
	for J in range(0,len(TRG_T_START)):
		plt.axvline(x=TRG_T_START[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_T_STOP[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_PRE_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
		plt.axvline(x=TRG_POST_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
	ax_LE.set_ylim(0,max(hist)+(max(hist)/3.0)) # modifica dell'altezza massima +1/3
	ax_LE.set_ylabel('[counts / %3.3f s]\n< 1.4 MeV\n' % BIN, fontsize=12)

	ax_HE = fig.add_subplot(12,1,7)
	ax_HE.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,7], markersize=5, c='blue', linewidth=1)
	ax_HE.set_xlim(TRG_PRE_BURST[M]-TTIME,TRG_POST_BURST[M]-TTIME)
	plt.axvline(x=0, linewidth=5, color = 'magenta', alpha = 0.3)
	for J in range(0,len(TRG_T_START)):
		plt.axvline(x=TRG_T_START[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_T_STOP[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_PRE_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
		plt.axvline(x=TRG_POST_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
	ax_HE.set_ylim(0,max(hist)+(max(hist)/3.0)) # modifica dell'altezza massima +1/3
	ax_HE.set_ylabel('[counts / %3.3f s]\n>= 1.4 MeV\n' % BIN, fontsize=12)

	ax_MULT4 = fig.add_subplot(12,1,8)
	ax_MULT4.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,11], markersize=5, c='blue', linewidth=1)
	ax_MULT4.set_xlim(TRG_PRE_BURST[M]-TTIME,TRG_POST_BURST[M]-TTIME)
	plt.axvline(x=0, linewidth=5, color = 'magenta', alpha = 0.3)
	for J in range(0,len(TRG_T_START)):
		plt.axvline(x=TRG_T_START[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_T_STOP[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_PRE_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
		plt.axvline(x=TRG_POST_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
	ax_MULT4.set_ylim(0,max(hist)+(max(hist)/3.0)) # modifica dell'altezza massima +1/3
	ax_MULT4.set_ylabel('[counts / %3.3f s]\n mult >= 4\n' % BIN, fontsize=12)

	ax_100MEV = fig.add_subplot(12,1,9)
	ax_100MEV.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,12], markersize=5, c='blue', linewidth=1)
	ax_100MEV.set_xlim(TRG_PRE_BURST[M]-TTIME,TRG_POST_BURST[M]-TTIME)
	plt.axvline(x=0, linewidth=5, color = 'magenta', alpha = 0.3)
	for J in range(0,len(TRG_T_START)):
		plt.axvline(x=TRG_T_START[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_T_STOP[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_PRE_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
		plt.axvline(x=TRG_POST_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
	ax_100MEV.set_ylim(0,max(hist)+(max(hist)/3.0)) # modifica dell'altezza massima +1/3
	ax_100MEV.set_ylabel('[counts / %3.3f s]\n>= 100 MeV\n' % BIN, fontsize=12)

	ax_SE = fig.add_subplot(12,1,10)
	ax_SE.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,10], markersize=5, c='green', linewidth=1)
	ax_SE.set_xlim(TRG_PRE_BURST[M]-TTIME,TRG_POST_BURST[M]-TTIME)
	plt.axvline(x=0, linewidth=5, color = 'magenta', alpha = 0.3)
	for J in range(0,len(TRG_T_START)):
		plt.axvline(x=TRG_T_START[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_T_STOP[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_PRE_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
		plt.axvline(x=TRG_POST_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
	ax_SE.set_ylim(0,100) # modifica dell'altezza massima +1/3
	ax_SE.set_ylabel('[Scounts / %3.3f s]\n' % BIN, fontsize=12)

	ax_FE = fig.add_subplot(12,1,11)
	ax_FE.step(LIGHTCURVE[:,0]-TTIME, LIGHTCURVE[:,9], markersize=5, c='red', linewidth=1)
	ax_FE.set_xlim(TRG_PRE_BURST[M]-TTIME,TRG_POST_BURST[M]-TTIME)
	plt.axvline(x=0, linewidth=5, color = 'magenta', alpha = 0.3)
	for J in range(0,len(TRG_T_START)):
		plt.axvline(x=TRG_T_START[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_T_STOP[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_PRE_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
		plt.axvline(x=TRG_POST_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
	ax_FE.set_ylim(0,100) # modifica dell'altezza massima +1/3
	ax_FE.set_ylabel('[Scounts / %3.3f s]\n' % BIN, fontsize=12)

	ax_CTS = fig.add_subplot(12,1,12)
	ax_CTS.scatter(np.array(EVT_RED)[:,0]-TTIME, np.array(EVT_RED)[:,1], s=2, c='cyan')
	ax_CTS.set_xlim(TRG_PRE_BURST[M]-TTIME,TRG_POST_BURST[M]-TTIME)
	plt.axvline(x=0, linewidth=5, color = 'magenta', alpha = 0.3)
	for J in range(0,len(TRG_T_START)):
		plt.axvline(x=TRG_T_START[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_T_STOP[J]-TTIME, linewidth=2, linestyle='dotted', color = 'red', alpha = 0.5)
		plt.axvline(x=TRG_PRE_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
		plt.axvline(x=TRG_POST_BURST[J]-TTIME, linewidth=2, linestyle='dotted', color = 'blue', alpha = 0.5)
	ax_CTS.set_ylim(0,1)
	ax_CTS.set_ylabel('E [MeV]', fontsize=12)

	plt.setp(ax_CTS.get_xticklabels(), fontsize=16)
	ax_CTS.set_xlabel('t - %9.6f [s]' % TTIME, fontsize=18)
	plt.savefig('%s_%06.d_%9.6f.png' % (NAME, CONT, TTIME))
	plt.close(fig)

if FLAG != 'SOLO' and FLAG != 'AGILE' and FLAG != 'KW':
	print "\nNon capire cosa tu dire, badrone\n"

print "hai introdotto una flag =", FLAG





























sys.exit()









# reconstruct the geographic position of the satellite

t3916 = []
LON = []
LAT = []
LON3916 = 0.0
LAT3916 = 0.0

# read the 3916 file and select the nearest in time LON and LAT

ROOTNAME3916  = ROOT.TFile.Open("%s/%06.d/RT%06.d_3916.root" % (OUTPATH_ROOT, CONT, CONT))
TREENAME3916  = ROOTNAME3916.Get("tdata3916")
CHAIN3916     = gDirectory.Get("tdata3916")
N_ENTRIES3916 = CHAIN3916.GetEntriesFast()

for atrg in ROOTNAME3916.tdata3916:
	t3916.append(float(atrg.time))
	LON.append(float(atrg.lon))
	LAT.append(float(atrg.lat))
ROOTNAME3916.Close()



diff = []

for f in range(0, len(t3916)):

	diff.append(abs(t3916[f]-TTIME))

F = int(diff.index(min(diff)))

LON3916 = float(LON[F])
LAT3916 = float(LAT[F])

print LON3916, LAT3916, t3916[F], float(diff[F])

sys.exit()






















































#######
# END #
#######

