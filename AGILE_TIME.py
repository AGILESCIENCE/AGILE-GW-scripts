#############################################################################################################

# USAGE 1: UTC to On-Board Time ---> python AGILE_TIME.py UTC2OBT   2015-01-01 00:00:00
# USAGE 2: On-Board Time to UTC ---> python AGILE_TIME.py OBT2UTC   123456789
# USAGE 3: UTC to Epoch Time    ---> python AGILE_TIME.py UTC2EPOCH 2015-01-01 00:00:00
# USAGE 4: Epoch Time to UTC    ---> python AGILE_TIME.py EPOCH2UTC 123456789
# USAGE 5: UTC to Local Time    ---> python AGILE_TIME.py UTC2LT    2015-01-01 00:00:00

#############################################################################################################

import sys
import time
from calendar import timegm
import datetime

#################################
# USAGE 1: UTC to On-Board Time #
#################################

if str(sys.argv[1]) == 'UTC2OBT':

	# data e tempo in formato YYYY-MM-DD hh:mm:ss
	DATE = str(sys.argv[2])
	TIME = str(sys.argv[3])
	UTC = '%s %s' % (DATE, TIME)

	# epoch time di AGILE (detto anche OBT) (da 2004-01-01 00:00:00 = 1970-01-01 00:00:00 + 1072915200 s)
	OBT = timegm(time.strptime(UTC, '%Y-%m-%d %H:%M:%S')) - 1072915200
	print "OBT", OBT

#################################
# USAGE 2: On-Board Time to UTC #
#################################

elif str(sys.argv[1]) == 'OBT2UTC':

	# tempo in formato OBT time 123456789
	OBT = int(float(sys.argv[2])) + 1072915200

	# tempo in UTC
	UTC = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(OBT))
	print "UTC", UTC

##############################
# USAGE 3: UTC to Epoch Time #
##############################

elif str(sys.argv[1]) == 'UTC2EPOCH':

	# data e tempo in formato YYYY-MM-DD hh:mm:ss
	DATE = str(sys.argv[2])
	TIME = str(sys.argv[3])
	UTC = '%s %s' % (DATE, TIME)

	# epoch time (da 1970-01-01 00:00:00)
	EPOCH = timegm(time.strptime(UTC, '%Y-%m-%d %H:%M:%S'))
	print "EPOCH", EPOCH

##############################
# USAGE 4: Epoch Time to UTC #
##############################

elif str(sys.argv[1]) == 'EPOCH2UTC':

	# tempo in formato epoch time 123456789
	EPOCH = int(sys.argv[2])

	# tempo in UTC
	UTC = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(EPOCH))
	print "UTC", UTC

##############################
# USAGE 5: UTC to Local Time #
##############################

elif str(sys.argv[1]) == 'UTC2LT':

	# data e tempo in formato YYYY-MM-DD hh:mm:ss
	DATE = str(sys.argv[2])
	TIME = str(sys.argv[3])
	UTC = '%s %s' % (DATE, TIME)
	EPOCH = timegm(time.strptime(UTC, '%Y-%m-%d %H:%M:%S'))

	# local time dipendente dal fuso orario (attenzione: dipende dal fuso orario del pc)
	LT = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(EPOCH))
	print "LOCAL TIME", LT

#########
# ERROR #
#########

else:
	print "\nErrore: inserire \'UTC2OBT\' o \'OBT2UTC\' o \'UTC2LT\' \n"
