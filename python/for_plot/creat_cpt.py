import sys;
import math;
import string;

flag = int(sys.argv[3]);

C = [];
C.append([0,0,0])
C.append([155  ,0      ,0])
C.append([255     ,56     ,0])
C.append([255     ,255     ,0])
#C.append([255     ,255     ,255])
C.append([255     ,255     ,255])
C.append([52     ,255     ,19])
C.append([0      ,154     ,205])
C.append([0       ,0      ,225])
C.append([128       ,0       ,128])


bb = float(sys.argv[1]);
ee = float(sys.argv[2]);
dd = (ee-bb)/8.

k = 0;
for i in range (8):
	if (flag != -1):
		print "%4.3g %g %g %g %4.3g %g %g %g" % (bb+i*dd,C[i][0],C[i][1],C[i][2],    bb+(i+1)*dd,C[i+1][0],C[i+1][1],C[i+1][2]);
		k = i+1;
	else:
		ii = 8-i;
		k = ii-1;
		print "%4.3g %g %g %g %4.3g %g %g %g" % (bb+i*dd,C[ii][0],C[ii][1],C[ii][2],    bb+(i+1)*dd,C[ii-1][0],C[ii-1][1],C[ii-1][2]);

print "F %g %g %g" % (C[k][0],C[k][1],C[k][2]);
