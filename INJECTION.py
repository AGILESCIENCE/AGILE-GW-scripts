# USAGE: python INJECTION.py INJ_TEST 123456789.123456 120.5 -40.5

import os, sys
import time
from calendar import timegm
import datetime

PATH_DATA = str(os.environ['PATH_DATA'])
MCALPIPE = str(os.environ['MCALPIPE'])
AGILEPIPE = str(os.environ['AGILEPIPE'])
MCALSW = str(os.environ['MCALSW'])

NAME = str(sys.argv[1])
TEMPO = str(sys.argv[2])
if int(len(TEMPO)) == 19:
        TTIME = timegm(time.strptime(TEMPO, '%Y-%m-%dT%H:%M:%S')) - 1072915200
else:
        TTIME = float(TEMPO)
LON = float(sys.argv[3])
LAT = float(sys.argv[4])

OUTPATH_ROOT = "%s/%s" % (str(os.getcwd()), NAME)

os.system("mkdir -p %s" % OUTPATH_ROOT)

print "OUTPATH_ROOT = ", OUTPATH_ROOT

os.chdir(OUTPATH_ROOT)

f_con = open("%s/%s_contour.con" % (OUTPATH_ROOT, NAME), "w")
f_con.write(	"%3.2f %3.2f\n" % (LON, LAT) +
		"%3.2f %3.2f\n" % (LON, LAT) +
		"%3.2f %3.2f\n" % (LON, LAT) +
		"%3.2f %3.2f\n" % (LON, LAT) +
		"%3.2f %3.2f\n" % (LON, LAT) +
		"%3.2f %3.2f\n" % (LON, LAT) +
		"%3.2f %3.2f\n" % (LON, LAT) +
		"%3.2f %3.2f\n" % (LON, LAT) +
		"%3.2f %3.2f\n" % (LON, LAT) +
		"%3.2f %3.2f\n" % (LON, LAT) +
		"%3.2f %3.2f\n" % (LON, LAT) +
		"%3.2f %3.2f\n" % (LON, LAT) +
		"%3.2f %3.2f\n" % (LON, LAT) +
		"%3.2f %3.2f\n" % (LON, LAT))
f_con.close()


f_rt = open("%s/%s_rt.ll" % (OUTPATH_ROOT, NAME), "w")
f_rt.write(	"source /opt/module/py_manual_27\n" +
		"source /opt/module/mcal_pipe_2.0\n" +
		"source activate py_mcal_27\n" +
		"source /opt/module/agile-B24-r5\n" +
		"source /opt/module/agile-mcal\n" +
		"source /opt/module/heasoft-6.25\n" +
		"xvfb-run -d -s \"-screen 0 2000x2000x24\" python %s/MCAL-ALERT.py %s/%s_run.xml\n" % (MCALPIPE, OUTPATH_ROOT, NAME) +
		"#rm -r %s/RESULTS" % OUTPATH_ROOT)
f_rt.close()


f_run = open("%s/%s_run.xml" % (OUTPATH_ROOT, NAME), "w")
f_run.write("<run id=\"000000\">\n" +
            "<parameter name=\"AlertInfo\" triggerid=\"0000000\" seqnum=\"1\" contname=\"%s_contour.con\" />\n" % NAME +
	    "<parameter name=\"TimeIntervals\" tmin=\"%9.6f\" tmax=\"%9.6f\" t0=\"%9.6f\" timeunit=\"s\" timesys=\"tt\"/>\n" % (float(TTIME)-100, float(TTIME)+100, float(TTIME)) +
            "<parameter name=\"Energy\" emin=\"0.0\" emax=\"10000.0\" energyBinID=\"\" />\n" +
	    "<parameter name=\"DirectoryList\" run=\"%s\" results=\"%s/RESULTS\" runprefix=\"000000\" fitspath=\"/data01/ASDC_PROC2/DATA_2/COR/\" />\n" % (OUTPATH_ROOT, OUTPATH_ROOT) +
            "<parameter name=\"DeleteRun\" value=\"0\" />\n" +
            "</run>")
f_run.close()


print "\n\nNOW LAUNCH THE FOLLOWING COMMAND:\nsource %s/%s_rt.ll\n\n" % (OUTPATH_ROOT, NAME)
