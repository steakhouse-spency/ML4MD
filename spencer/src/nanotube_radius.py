from math import pi, sqrt, cos, sin, radians

# global variables
bl = 1.42
a = bl * cos(radians(30)) * 2
need = [105, 106, 115, 124, 144, 152]
m_n = [(11,4), (13, 1), (13, 3), (9, 8), (14, 7), (15, 7)]

# get diameter of nano tube in angstrom
def getDiam(m, n):
	return 2 * a * sqrt(m*m + m*n + n*n) / (2*pi)

def test(i):
	for n in range(i):
		for m in range(n, i):
			d = getDiam(m, n)
			if d >= 10 and d < 16 and int(d*10) in need:
				print("(%d, %d) = %f" % (m, n, d))


print("a = %f" % (a))
print("(m, n)")
test(25)

