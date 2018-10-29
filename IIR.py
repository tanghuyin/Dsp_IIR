import math
import utils
from sympy import *
from sympy.plotting import plot
# Hz to rad/s
def f2Omega(f):
	Omega = f * math.pi * 2
	return Omega
def f2w(f, f_k):
	w = 2 * math.pi * f / f_k 
	return w

# Get g
def getG(As, Rp):
	tmp1 = pow(10, 0.1 * As) - 1
	tmp2 = pow(10, 0.1 * Rp) - 1
	g = pow(tmp1 / tmp2, 0.5)
	return g

# Get Lambda
def getLambda(f_st, f_p):
	Lambda = f2Omega(f_st) / f2Omega(f_p)
	return Lambda

# Get the order of a Butterworth filter
def getOrderOfButterworth(As, Rp, f_st, f_p):
	g = getG(As, Rp)
	Lambda = getLambda(f_st, f_p)
	N = math.ceil(math.log10(g) / math.log10(Lambda))
	N = int(N)
	return N

# Get the order of a ChebyshevI filter
def getOrderOfChebyshevI(As, Rp, f_st, f_p):
	g = getG(As, Rp)
	Lambda = getLambda(f_st, f_p)
	numerator = math.log10(g + pow(g*g - 1, 0.5))
	dominator = math.log10(Lambda + pow(Lambda*Lambda - 1, 0.5))
	N = math.ceil(numerator / dominator)
	#How to choose N is still a problem
	if (N - numerator/dominator) > 0.9:
		N -= 1
	N = int(N)
	return N

# Get epsilon of ChebyshevI filter
def getEpsilon(Rp):
	epsilon = pow((pow(10, 0.1*Rp) - 1), 0.5)
	return epsilon

# Get the 3dB cut off frequency in rad/s
def get3dBOmega(f_p, Rp, N):
	Omega_p = f2Omega(f_p)
	tmp1 = pow(10, 0.1 * Rp) - 1
	numerator = Omega_p
	dominator = pow(tmp1, 1.0 / (2 * N))
	Omega_c_3db = numerator / dominator
	return Omega_c_3db

'''
TO BE deleted
# Find the dominator of Normal ButterWorth from the table
def getNormalButterWorth(N):
	if N < 1 or N > 7:
		raise Exception('Order Error')
	else:
		return utils.Dominator_Butt[N -1]


# DeNormalzing the expression and get numerator and dominator
def deNormalButterWorth(Omega_c_3db, Dominator_Butt, N):
	numerator = pow(Omega_c_3db, N)
	dominator = []
	for i in range(N + 1):
		dominator.append(Dominator_Butt[N - i] * pow(Omega_c_3db, i))
	return numerator, dominator
'''
def getButterWorth(N, Omega_c_3dB):
	s, k = symbols('s, k')
	Hks = Omega_c_3dB**2 / (s**2 - 2*s*Omega_c_3dB*cos(pi/2 + pi*(2*k-1)/(2*N)) + Omega_c_3dB**2)
	if N % 2 == 0:
		Has = 1
	else:
		Has = Omega_c_3dB / (s + Omega_c_3dB)
		Has.n(5)
	for i in range(1, N / 2 +1):
		Has = Has * Hks.subs(k, i).evalf()

	return expand(Has)

def getChebyshevI(N, epsilon, Omega_p):
	s, k = symbols('s, k')
	a = sinh(asinh(1/epsilon)/N)
	b = cosh(acosh(1/epsilon)/N)
	deltak = -a * sin(pi*(2*k-1) / (2*N))
	omegak = b * cos(pi*(2*k-1) / (2*N))
	sk = deltak + I * omegak
	Hans = 1
	for i in range(1, N+1):
		Hans = Hans * (s - sk).subs(k, i)
	dominator = expand(Hans).evalf()
	numerator = 1 / (epsilon * 2**(N-1))
	Hans = (numerator / dominator)
	Has = Hans.subs(s, s / Omega_p)
	return Has

# # Find the dominator of Normal ButterWorth from the table
# def getNormalChebyshevI(N, dB):
# 	if N < 1 or N > 7:
# 		raise Exception("Order Error")
# 	elif dB != 0.5 and dB != 1 and dB != 2 and dB != 3:
# 		raise Exception("dB Error")

# 	else:
# 		if dB == 0.5:
# 			return utils.Dominator_Chebyshev_1_2dB[N - 1]
# 		if dB == 1:
# 			return utils.Dominator_Chebyshev_1dB[N - 1]
# 		if dB == 2:
# 			return utils.Dominator_Chebyshev_2dB[N - 1]
# 		if dB == 3:
# 			return utils.Dominator_Chebyshev_3dB[N - 1]	

# # Get numerator and dominator of Chebyshev
# def ChebyshevI(Rp, N, Dominator_Chebyshev):
# 	epsilon = getEpsilon(Rp)
# 	numerator = 1 / (epsilon * pow(2, N-1))
# 	dominator = Dominator_Chebyshev
# 	return numerator, dominator

# Analog to Analog/lp to bp
def A2ALp2Bp(f_p1, f_p2, f_st1, f_st2, As, Rp, isChebyOrButt):
	s = symbols('s')
	Omega_p1, Omega_p2, Omega_st1, Omega_st2 = f2Omega(f_p1), f2Omega(f_p2), f2Omega(f_st1), f2Omega(f_st2)
	Bp = Omega_p2 - Omega_p1
	Omega_p0 = pow(Omega_p1 * Omega_p2, 0.5)
	normal_Omega_st2 = (pow(Omega_st2, 2) - pow(Omega_p0, 2)) / (Omega_st2 * Bp)
	normal_Omega_st1 = (pow(Omega_st1, 2) - pow(Omega_p0, 2)) / (Omega_st1 * Bp)
	normal_Omega_st = min(abs(normal_Omega_st1), abs(normal_Omega_st2))
	if isChebyOrButt == 1:
		N = getOrderOfChebyshevI(As, Rp, normal_Omega_st, 1)
		epsilon = getEpsilon(Rp)
		Has = getChebyshevI(N, epsilon, 1)
		Hbp = Has.subs(s, (s**2 + Omega_p0**2) / (s * Bp))
		Hbp = simplify(Hbp).n(5)
		return Hbp
	else:
		N = getOrderOfButterworth(As, Rp, normal_Omega_st, 1)
		Omega_c_3dB = get3dBOmega(1 / (2 * math.pi), Rp, N)
		Has = getButterWorth(N, Omega_c_3dB)
		Hbp = Has.subs(s, (s**2 + Omega_p0**2) / (s * Bp))
		Hbp = simplify(Hbp).n(5)
		return Hbp

def A2ALp2Lp(f_p, f_st, As, Rp, isChebyOrButt):
	s = symbols('s')
	Omega_p, Omega_st = f2Omega(f_p), f2Omega(f_st)
	normal_Omega_st = Omega_st / Omega_p
	if isChebyOrButt == 1:
		N = getOrderOfChebyshevI(As, Rp, normal_Omega_st, 1)
		print N
		epsilon = getEpsilon(Rp)
		Has = getChebyshevI(N, epsilon, 1)
		Hlp = Has.subs(s, s / Omega_p)
		Hlp = simplify(Hlp).n(5)
		return Hlp
	else:
		N = getOrderOfButterworth(As, Rp, normal_Omega_st, 1)

		Omega_c_3dB = get3dBOmega(1 / (2 * math.pi), Rp, N)
		Has = getButterWorth(N, Omega_c_3dB)
		Hlp = Has.subs(s, s / Omega_p)
		Hlp = simplify(Hlp).n(5)
		return Hlp

def A2ALp2Hp(f_p, f_st, As, Rp, isChebyOrButt):
	s = symbols('s')
	Omega_p, Omega_st = f2Omega(f_p), f2Omega(f_st)
	normal_Omega_st =  -(Omega_p / Omega_st)
	if isChebyOrButt == 1:
		N = getOrderOfChebyshevI(As, Rp, normal_Omega_st, 1)
		epsilon = getEpsilon(Rp)
		Has = getChebyshevI(N, epsilon, 1)
		Hhp = Has.subs(s, Omega_p / s)
		Hhp = simplify(Hhp).n(5)
		return Hhp
	else:
		N = getOrderOfButterworth(As, Rp, abs(normal_Omega_st), 1)
		Omega_c_3dB = get3dBOmega(1 / (2 * math.pi), Rp, N)
		Has = getButterWorth(N, Omega_c_3dB)
		Hhp = Has.subs(s, Omega_p / s)
		Hhp = simplify(Hhp).n(5)
		return Hhp

# There is a BUG in this Function
def A2ALp2Bs(f_p1, f_p2, f_st1, f_st2, As, Rp, isChebyOrButt):
	s = symbols('s')
	nOmegast = symbols('nOmegast')
	Omega_p1, Omega_p2, Omega_st1, Omega_st2 = f2Omega(f_p1), f2Omega(f_p2), f2Omega(f_st1), f2Omega(f_st2)
	Bs = Omega_st2 - Omega_st1
	Omega_st0 = Omega_st1 * Omega_st2 # == Omega_st0**2
	normal_Omega_p2 = ((nOmegast * Bs * Omega_p2) / (Omega_st0 - Omega_p2**2))
	normal_Omega_p1 = ((nOmegast * Bs * Omega_p1) / (Omega_st0 - Omega_p1**2))
	a = normal_Omega_p1.subs(nOmegast, 1)
	b = normal_Omega_p2.subs(nOmegast, 1)
	normal_Omega_st = 1 / max(abs(a), abs(b))
	if isChebyOrButt == 1:
		N = getOrderOfChebyshevI(As, Rp, normal_Omega_st, 1)
		epsilon = getEpsilon(Rp)
		Has = getChebyshevI(N, epsilon, 1)
		Hbs = Has.subs(s, (s*Bs*normal_Omega_st) / (s**2 + Omega_st0))
		Hbs = simplify(Hbs).n(5)
		return Hbs
	else:
		N = getOrderOfButterworth(As, Rp, normal_Omega_st, 1)
		Omega_c_3dB = get3dBOmega(1 / (2 * math.pi), Rp, N)
		Has = getButterWorth(N, Omega_c_3dB)
		Hbs = Has.subs(s, (s*Bs*normal_Omega_st) / (s**2 + Omega_st0))
		Hbs = simplify(Hbs).n(5)
		return Hbs

def A2DLp2Bp(f_p1, f_p2, f_st1, f_st2, f_k, As, Rp, isChebyOrButt):
	w, z, s = symbols('w, z, s')
	w_p1 = f2Omega(f_p1) / f_k
	w_p2 = f2Omega(f_p2) / f_k
	w_st1 = f2Omega(f_st1) / f_k
	w_st2 = f2Omega(f_st2) / f_k
	w_0 = acos(cos((w_p2 + w_p1)/2) / cos((w_p2 - w_p1)/2))
	Omega = (cos(w_0) - cos(w)) / sin(w)
	Omega_p = Omega.subs(w, w_p2).evalf()
	Omega_st1 = Omega.subs(w, w_st1).evalf()
	Omega_st2 = Omega.subs(w, w_st2).evalf()
	Omega_st = min(abs(Omega_st1), abs(Omega_st2))
	if isChebyOrButt == 1:
		N = getOrderOfChebyshevI(As, Rp, Omega_st, Omega_p)
		epsilon = getEpsilon(Rp)
		Has = getChebyshevI(N, epsilon, Omega_p)
		Hbp = Hap.subs(s, (1 - 2 * (z**-1) * cos(w_0) + z**(-2)) / (1 - z**(-2)))
		Hbp = simplify(Hbp).n(5)
		return Hbp
	else:
		N = getOrderOfButterworth(As, Rp, Omega_st, Omega_p)
		Omega_c_3dB = get3dBOmega(Omega_p / (2 * math.pi), Rp, N)
		Has = getButterWorth(N, Omega_c_3dB)
		Hbp = Has.subs(s, (1 - 2 * (z**-1) * cos(w_0) + z**(-2)) / (1 - z**(-2)))
		# Hbp = expand(Hbp)
		Hbp = simplify(Hbp).n(5)
		amp = abs(Hbp.subs(z, exp(I*w)))
		# plot(20*log(amp), (w, 0, pi))
		return Hbp

def A2DLp2Bs(f_p1, f_p2, f_st1, f_st2, f_k, As, Rp, isChebyOrButt):
	w, z, s = symbols('w, z, s')
	w_p1 = f2Omega(f_p1) / f_k
	w_p2 = f2Omega(f_p2) / f_k
	w_st1 = f2Omega(f_st1) / f_k
	w_st2 = f2Omega(f_st2) / f_k
	w_0 = acos(cos((w_st2 + w_st1)/2) / cos((w_st2 - w_st1)/2))
	Omega = sin(w) / (cos(w) - cos(w_0)) 
	Omega_st = Omega.subs(w, w_st1).evalf()
	Omega_p1 = Omega.subs(w, w_p1).evalf()
	Omega_p2 = Omega.subs(w, w_p2).evalf()
	Omega_p = max(abs(Omega_p1), abs(Omega_p2))
	if isChebyOrButt == 1:
		N = getOrderOfChebyshevI(As, Rp, Omega_st, Omega_p)
		epsilon = getEpsilon(Rp)
		Has = getChebyshevI(N, epsilon, Omega_p)
		Hbs = Has.subs(s, ((1 - z**(-2)) / (1 - 2 * (z**-1) * cos(w_0) + z**(-2))))
		Hbs = simplify(Hbs).n(5)
		# amp = abs(Hbs.subs(z, exp(I*w)))
		# plot(20*log(amp), (w, 0, pi))
		return Hbs
	else:
		N = getOrderOfButterworth(As, Rp, Omega_st, Omega_p)
		Omega_c_3dB = get3dBOmega(Omega_p / (2 * math.pi), Rp, N)
		Has = getButterWorth(N, Omega_c_3dB)
		Hbs = Has.subs(s, ((1 - z**(-2)) / (1 - 2 * (z**-1) * cos(w_0) + z**(-2))))
		# Hbp = expand(Hbp)
		Hbs = simplify(Hbs).n(5)
		# amp = abs(Hbp.subs(z, exp(I*w)))
		# plot(20*log(amp), (w, 0, pi))
		return Hbs
		
def A2DLp2Hp(f_p, f_st, f_k, As, Rp, isChebyOrButt):
	w, z, s = symbols('w, z, s')
	w_p = f2Omega(f_p) / f_k
	w_st = f2Omega(f_st) / f_k
	Omega = cot(w / 2)
	Omega_p = Omega.subs(w, w_p).evalf()
	Omega_st = Omega.subs(w, w_st).evalf()
	if isChebyOrButt == 1:
		N = getOrderOfChebyshevI(As, Rp, Omega_st, Omega_p)
		epsilon = getEpsilon(Rp)
		Has = getChebyshevI(N, epsilon, Omega_p)
		Hhp = Has.subs(s, ((1 + z**(-1)) / (1 - z**(-1))))
		Hhp = simplify(Hhp).n(5)
		return Hhp
	else:
		N = getOrderOfButterworth(As, Rp, Omega_st, Omega_p)
		Omega_c_3dB = get3dBOmega(Omega_p / (2 * math.pi), Rp, N)
		Has = getButterWorth(N, Omega_c_3dB)
		Hhp = Has.subs(s, ((1 + z**(-1)) / (1 - z**(-1))))
		# Hbp = expand(Hbp)
		Hhp = simplify(Hhp).n(5)
		# amp = abs(Hbp.subs(z, exp(I*w)))
		# plot(20*log(amp), (w, 0, pi))
		return Hhp

def A2DLp2Lp(f_p, f_st, f_k, As, Rp, isChebyOrButt):
	w, z, s = symbols('w, z, s')
	w_p = f2Omega(f_p) / f_k
	w_st = f2Omega(f_st) / f_k
	Omega = tan(w / 2)
	Omega_p = Omega.subs(w, w_p).evalf()
	Omega_st = Omega.subs(w, w_st).evalf()
	if isChebyOrButt == 1:
		N = getOrderOfChebyshevI(As, Rp, Omega_st, Omega_p)
		epsilon = getEpsilon(Rp)
		Has = getChebyshevI(N, epsilon, Omega_p)
		Hlp = Has.subs(s, ((1 - z**(-1)) / (1 + z**(-1))))
		Hlp = simplify(Hlp).n(5)
		return Hlp
	else:
		N = getOrderOfButterworth(As, Rp, Omega_st, Omega_p)
		Omega_c_3dB = get3dBOmega(Omega_p / (2 * math.pi), Rp, N)
		Has = getButterWorth(N, Omega_c_3dB)
		Hlp = Has.subs(s, ((1 - z**(-1)) / (1 + z**(-1))))
		# Hbp = expand(Hbp)
		Hlp = simplify(Hlp).n(5)
		# amp = abs(Hlp.subs(z, exp(I*w)))
		# plot(20*log(amp), (w, 0, pi))
		return Hlp

def A2A2DLp2Lp(w_p, w_st, As, Rp, isChebyOrButt):
	z, s = symbols('z, s')
	Omega_p = tan(w_p / 2)
	Omega_st = tan(w_st / 2)
	f_p = Omega_p / (2 * math.pi)
	f_st = Omega_st / (2 * math.pi)
	Has = A2ALp2Lp(f_p, f_st, As, Rp, isChebyOrButt)
	Hlp = Has.subs(s, (1-z**(-1)) / (1+z**(-1)))
	Hlp = simplify(Hlp).n(5)
	return Hlp

def A2A2DLp2Hp(w_p, w_st, As, Rp, isChebyOrButt):
	z, s = symbols('z, s')
	Omega_p = tan(w_p / 2)
	Omega_st = tan(w_st / 2)
	f_p = Omega_p / (2 * math.pi)
	f_st = Omega_st / (2 * math.pi)
	Has = A2ALp2Hp(f_p, f_st, As, Rp, isChebyOrButt)
	Hhp = Has.subs(s, (1-z**(-1)) / (1+z**(-1)))
	Hhp = simplify(Hhp).n(5)
	return Hhp

def A2A2DLp2Bp(w_p1, w_p2, w_st1, w_st2, As, Rp, isChebyOrButt):
	z, s = symbols('z, s')
	Omega_p1 = tan(w_p1 / 2)
	Omega_st1 = tan(w_st1 / 2)
	Omega_p2 = tan(w_p2 / 2)
	Omega_st2 = tan(w_st2 / 2)
	f_p1 = Omega_p1 / (2 * math.pi)
	f_st1 = Omega_st1 / (2 * math.pi)
	f_p2 = Omega_p2 / (2 * math.pi)
	f_st2 = Omega_st2 / (2 * math.pi)
	Has = A2ALp2Bp(f_p1, f_p2, f_st1, f_st2, As, Rp, isChebyOrButt)
	Hbp = Has.subs(s, (1-z**(-1)) / (1+z**(-1)))
	Hbp = simplify(Hbp).n(5)
	return Hbp


def A2A2DLp2Bs(w_p1, w_p2, w_st1, w_st2, As, Rp, isChebyOrButt):
	z, s = symbols('z, s')
	Omega_p1 = tan(w_p1 / 2)
	Omega_st1 = tan(w_st1 / 2)
	Omega_p2 = tan(w_p2 / 2)
	Omega_st2 = tan(w_st2 / 2)
	f_p1 = Omega_p1 / (2 * math.pi)
	f_st1 = Omega_st1 / (2 * math.pi)
	f_p2 = Omega_p2 / (2 * math.pi)
	f_st2 = Omega_st2 / (2 * math.pi)
	Has = A2ALp2Bs(f_p1, f_p2, f_st1, f_st2, As, Rp, isChebyOrButt)
	Hbs = Has.subs(s, (1-z**(-1)) / (1+z**(-1)))
	Hbs = simplify(Hbs).n(5)
	return Hbs

def drawH(H):
	z, w = symbols('z, w')
	amp = abs(H.subs(z, exp(I*w)))
	plot(20*log(amp), (w, 0, pi))



