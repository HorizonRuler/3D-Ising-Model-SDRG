'''
plot the scaling curves of the average magnetic moment
ln(m(2L)/m(L))/ln(L)=alpha if m(L)~L^(alpha)
graph alpha against 1/L where L is the one used to calc alpha
for the right st dev should get straight line
just do for all st dev
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

size8 = np.genfromtxt('first_order_SDRG_magnetization_moment_8.txt', delimiter=' : ')
size16 = np.genfromtxt('first_order_SDRG_magnetization_moment_16.txt', delimiter=' : ')
size32 = np.genfromtxt('first_order_SDRG_magnetization_moment_32.txt', delimiter=' : ')
alpha8 = np.log(size16[:,1] / size8[:,1]) / np.log(8)
alpha16 = np.log(size32[:,1] / size16[:,1]) / np.log(16)