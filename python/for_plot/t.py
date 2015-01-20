import sys
def cpt (m):
	y= m*m*m-8*m+6
	dy = 3*m*m-8
	v= m-y/dy
	if y/dy<0.01:
		print "haha"
	return y,dy,y/dy,v
m=float(sys.argv[1])

print cpt(m)

