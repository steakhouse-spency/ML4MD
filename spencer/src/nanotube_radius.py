from math import pi, sqrt, cos, sin, radians, acos, degrees

# global variables
# bl = 1.42
# a = bl * cos(radians(30)) * 2
# need = [105, 106, 115, 124, 144, 152]
# m_n = [(11,4), (13, 1), (13, 3), (9, 8), (14, 7), (15, 7)]

# function to compute nano tube diameter
def getDiam(n, m):
	bl = 1.42 # bond length for tube
	a = bl * cos(radians(30)) * 2
	return 2 * a * sqrt(m*m + m*n + n*n) / (2*pi)

def getC(n, m):
	return 246 * sqrt( ((n+m)**2) - n*m )

def getAngle(n, m, c):
	return degrees(acos((n+(m/2)) / c))


# def test(i):
# 	for n in range(i):
# 		for m in range(n, i):
# 			d = getDiam(n, m)
# 			if d >= 10 and d < 16 and int(d*10) in need:
# 				print("(%d, %d) = %f" % (m, n, d))


# print("a = %f" % (a))
# print("(m, n)")
# test(25)

c = getC(7, 14)
a = getAngle(7, 14, c)
print(a)