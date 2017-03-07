import numpy as np
import scipy.linalg as sp
import sys, ast
import xmltodict as x2d

def Aqp(omega_0, Ap):
    """Given the free particle Ap matrix and the frequency of the harmonic oscillator, computes the full drift matrix."""
    dAqp = np.zeros(np.asarray(Ap.shape) + 1)
    dAqp[0,1] = -np.power(omega_0,2)
    dAqp[1,0] = 1
    dAqp[1:,1:] = Ap.T
    return dAqp

def Dqp(omega_0, Dp):
    """Given the free particle Dp matrix and the frequency of the harmonic oscillator, computes the full D matrix."""
    dDqp = np.zeros(np.asarray(Dp.shape) + 1)
    dDqp[1:,1:] = Dp.T
    return dDqp

def Cqp(omega_0, Ap, Dp):
    """Given the free particle Ap and Dp matrices and the frequency of the harmonic oscillator, computes the full covariance matrix."""
    dAqp = Aqp(omega_0, Ap)
    dDqp = Dqp(omega_0, Dp)
    return sp.solve_continuous_are( -dAqp, np.zeros(dAqp.shape), dDqp, np.eye(dAqp.shape[-1]))

def Cvv(omega, omega_0, Ap, Dp, dw):
    """Given the Cp and Dp matrices for a harmonic oscillator of frequency omega_0, computes the value of the Fourier transform of the velocity velocity auto-correlation function."""
    try:
        #omega = np.maximum(omega, dw/100)
        #omega = np.maximum(omega_0, dw/100)
        dAqp = Aqp(omega_0, Ap)
        dAqp[1,1] = np.maximum(dAqp[1,1], dAqp[1,1] +2.0*dw)
        #dAqp[1,1] += 2*dw
        dDqp = Dqp(omega_0, Dp)
        dCqp = Cqp(omega_0, Ap, Dp)
        domega2 = np.eye(dAqp.shape[-1])*np.power(omega,2)
        dAqp2 = np.dot(dAqp,dAqp)
        dCvv = np.linalg.inv(dAqp2 + domega2)
        dCvv = np.dot(np.dot(dAqp, dCvv), dCqp)
        return dCvv[1,1]/(np.pi/2.0)
    except:
        return 0.0


def gleKernel(omega, Ap, Dp):
    """Given the Cp and Dp matrices for a harmonic oscillator of frequency omega_0, constructs the gle kernel for transformation of the velocity velocity autocorrelation function."""
    delta_omega = abs(omega[1]-omega[0])
    ngrid = len(omega) 
    dKer = np.zeros((ngrid,ngrid), float)
    for x in xrange(ngrid):
        for y in xrange(ngrid):
            dKer[x,y] = Cvv(omega[x], omega[y], Ap, Dp, delta_omega)
    return dKer*delta_omega

def ISRA(omega, ker, y, steps=500):
    """Given the thermostatted vvac spectrum and the range of frequencies, constructs the vibration density of states"""
    delta_omega = abs(omega[1]-omega[0])
    ngrid = len(omega)
    f = y
    cnvg= np.zeros((steps,2))
    for i in xrange(steps):
        h = np.dot(y,ker)
        D = np.dot(ker,ker.T)
        f = f*h/np.dot(f,D)
        cnvg[i] = np.asarray((np.linalg.norm((np.dot(f,ker.T) - y))**2, np.linalg.norm(np.gradient(np.gradient(f)))**2))
    return f, cnvg
