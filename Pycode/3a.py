__author__ = 'VikiQiu'
import numpy as np
from read_data import getData2
from scipy.optimize import minimize

# get data
filename = '../data/S21_5.csv'
f, s21, Hf = getData2(filename)

# initial value
tn = 9.6e-9#1.201e-9
G0 = 1.8e6#8.486e5
N0 = 4.97e5#1.286e6
tp = 3.8e-12#2.884e-12
epsilon = 4.7e-8#3.888e-6
beta = 1e-5#2.68e-2
k = 1.5e-8#4.166e-8
P0 = 1.904172485e-3

q = 1.6e-19
ita = 3.645514e-01
I = 7.5e-3
Ith = I-P0/ita

#### functions ####
def Ns_pre(theta):
    tp, tn, G0, k, epsilon, beta, N0 = theta
    return (P0/(k*tp)+G0*N0*P0/(k+epsilon*P0))/(beta/tn+G0*P0/(k+epsilon*P0))

def Ss_pre(theta):
    tp, tn, G0, k, epsilon, beta, N0 = theta
    Ns = Ns_pre(theta)
    return (P0/q-Ns/tn)/(G0*(Ns-N0))

def y_pre(theta):
    tp, tn, G0, k, epsilon, beta, N0 = theta
    Ss = Ss_pre(theta)
    Ns = Ns_pre(theta)
    return 1/tp + 1/tn +G0*Ss/(1+epsilon*Ss) - G0*(Ns-N0)/(1+epsilon*Ss)

def z_pre(theta):
    tp, tn, G0, k, epsilon, beta, N0 = theta
    Ss = Ss_pre(theta)
    Ns = Ns_pre(theta)
    return (1/tp - G0*(Ns-N0)/(1+epsilon*Ss))*(1/tn+G0*Ss/(1+epsilon*Ss))+(beta/tn+G0*Ss/(1+epsilon*Ss))*G0*(Ns-N0)/(1+epsilon*Ss)

# stationary I
def I_pre(theta):
    tp, tn, G0, k, epsilon, beta, N0 = theta
    Ss = Ss_pre(theta)
    Ns = Ns_pre(theta)
    return q*(Ns/tn+G0*(Ns-N0)*P0/(k+epsilon*P0))/ita+Ith

# stationary P
def P_pre(theta):
    tp, tn, G0, k, epsilon, beta, N0 = theta
    Ns = Ns_pre(theta)
    return k*(P0/q-Ns/tn)/G0/(Ns-N0)

def Hf_loss(theta, f, Hf2):
    y = y_pre(theta)
    z = z_pre(theta)
    h = z**2/(4*np.pi**2*f**2*y**2+z**2+(4*np.pi**2*f**2)**2-8*np.pi**2*f**2*z)
    return (h-Hf2)**2

def loss_fn(theta, f, Hf2):
    I1 = I_pre(theta)
    P1 = P_pre(theta)
    return np.mean(Hf_loss(theta, f, Hf2)+(I1/I-1)**2)+(P1/P0-1)**2

#### optimization ####
print(loss_fn(np.array([tp, tn, G0, k, epsilon, beta, N0]), f, Hf))
xopt=minimize(loss_fn, np.array([tp, tn, G0, k, epsilon, beta, N0]), args=(f, Hf), tol=1e-5,
              options=dict(disp=True, maxiter = 1e20), method='Powell')
print(xopt)
print(xopt.success)
print(xopt.x)
