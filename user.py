# -*- coding: UTF-8 -*-
import IIR
from sympy import *
#isChebyOrButt = 0 巴特沃斯 / = 1 切比雪夫
#=====模拟低通->模拟低通=========
# H = IIR.A2ALp2Lp(f_p, f_st, As, Rp, isChebyOrButt)
#=====模拟低通->模拟高通=========
# H = IIR.A2ALp2Hp(f_p, f_st, As, Rp, isChebyOrButt)
#=====模拟低通->模拟带通=========
# H = IIR.A2ALp2Bp(f_p1, f_p2, f_st1, f_st2, As, Rp, isChebyOrButt)
#=====模拟低通->模拟带阻=========
# H = IIR.A2ALp2Bs(f_p1, f_p2, f_st1, f_st2, As, Rp, isChebyOrButt)
#===========双线性==============
#=====模拟低通->数字低通=========
# H = IIR.A2A2DLp2Lp(w_p, w_st, As, Rp, isChebyOrButt)
#=====模拟低通->数字高通=========
# H = IIR.A2A2DLp2Hp(w_p, w_st, As, Rp, isChebyOrButt)
#=====模拟低通->数字带通=========
# H = IIR.A2A2DLp2Bp(w_p1, w_p2, w_st1, w_st2, As, Rp, isChebyOrButt)
#=====模拟低通->数字带阻=========
# H = IIR.A2A2DLp2Bs(w_p1, w_p2, w_st1, w_st2, As, Rp, isChebyOrButt)
#===========直接型==============
#=====模拟低通->数字低通=========
# H = IIR.A2DLp2Lp(f_p, f_st, f_k, As, Rp, isChebyOrButt):
#=====模拟低通->数字高通=========
# H = IIR.A2DLp2Hp(f_p, f_st, f_k, As, Rp, isChebyOrButt)
#=====模拟低通->数字带通=========
# H = IIR.A2DLp2Bp(f_p1, f_p2, f_st1, f_st2, f_k, As, Rp, isChebyOrButt)
#=====模拟低通->数字带阻=========
# H = IIR.A2DLp2Bs(f_p1, f_p2, f_st1, f_st2, f_k, As, Rp, isChebyOrButt)
#============输出===============
# print pretty(H)
#============绘图===============
# print IIR.drawH(H)
#============转换===============
# IIR.f2Omega(f)
# IIR.f2w(f, f_k)
