#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import gzip
import numpy as np
import struct
from scipy.optimize import curve_fit

def readgoftcube_txt(filename,compress=0):
    if (compress):
        f=gzip.open(filename,"r")
    else:
        f=open(filename,"r")
    alllines=f.readlines()
    chiantifile=alllines[3]
    abundfile=alllines[4]
    textlines=alllines[5:]
    data=[]
    for line in textlines:
        values=line.split()
        fval=[]
        for i in range(len(values)):
            fval.append(float(values[i]))
        data.append(fval)
    # now data contains
    # data[:][0] = x coordinate
    # data[:][1] = y coordinate
    # data[:][2] = wave coordinate
    # data[:][3] = intensity
    
    f.close()
    return np.array(data),chiantifile
    
def readgoftcube_dat(filename,compress=0):
    if (compress):
        f=gzip.open(filename,"rb")
    else:
        f=open(filename,"rb")
    
    int_in=f.read(4)
    dim=struct.unpack('i',int_in)[0]
    int_in=f.read(4)
    ng=struct.unpack('i',int_in)[0]
    int_in=f.read(4)
    nvars=struct.unpack('i',int_in)[0]
    int_in=f.read(8)
    chiantisize=struct.unpack('q',int_in)[0]
    chiantifile=f.read(chiantisize)
    int_in=f.read(8)
    abundsize=struct.unpack('q',int_in)[0]
    abundfile=f.read(abundsize)
    
    data=[]
    for i in range(ng):
        row=[]
        for j in range(dim+nvars):
            row.append(0.)
        data.append(row)
        
    for j in range(dim+nvars):
        for i in range(ng):
            float_in=f.read(4)
            data[i][j]=struct.unpack('f',float_in)[0]
    
    f.close()
    return np.array(data),chiantifile
    
def readgoftcubechianti(filename):
    parts=filename.split(".")
    compress=0
    if parts[-1]=="gz":
        compress=1
    
    if parts[-1-compress]=="txt":
        data,chiantifile=readgoftcube_txt(filename,compress=compress)
    
    if parts[-1-compress]=="dat":
        data,chiantifile=readgoftcube_dat(filename,compress=compress)
        
    # print len(data), len(data[:][0])
        
    return data,chiantifile

def readgoftcube(filename):
    data,chiantifile=readgoftcubechianti(filename)
    return data
    
def regulargoftcube(data):
    ng = len(data)
    nx = ng/np.shape(np.where(data[:,0]==data[0,0]))[1]
    ny = ng/np.shape(np.where(data[:,1]==data[0,1]))[1]
    if (ng/nx/ny == 1):
	    nl=1
	    dim=2
    else:
            nl = ng/np.shape(np.where(data[:,2]==data[0,2]))[1]
	    dim=3
    
    xvec=np.empty(nx)
    yvec=np.empty(ny)
    emiss=np.empty([nx,ny])
    lvec=1.

    if (dim == 3):
            lvec=np.empty(nl)
            emiss=np.empty([nx,ny,nl])

    for i in range(ng):
        j=i/nl
        k=i/nl/nx
        j-=k*nx
        l=i%nl
        xvec[j]=data[i,0]
        yvec[k]=data[i,1]
	if (dim == 3):
		lvec[l]=data[i,2]
		emiss[j,k,l]=data[i,dim]
	if (dim == 2):
		emiss[j,k]=data[i,dim]
    
    # there's a very weird ordering of the indices, it seems x needs to come last, in order to have an intuitive access with matplotlib
    emiss=np.swapaxes(emiss,0,dim-1)
    return emiss,xvec,yvec,lvec
    
def fit_gauss(x, a0, a1, a2):
#    print a0, a1, a2
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / 2)
    return y
    
def gaussfitgoftcube(emiss,lvec):
    shape=np.shape(emiss)
    nx=shape[2]
    ny=shape[1]
    nl=shape[0]
    peak=np.zeros([ny,nx])
    doppler=np.zeros([ny,nx])
    sigma=np.zeros([ny,nx])
    chisq=np.zeros([ny,nx])
    intens=np.zeros([ny,nx])
    
    l0=np.mean(lvec)
    errors=np.mean(emiss)
    for i in range(nx):
        for j in range(ny):
            spectralline=emiss[:,j,i]
            # localmem=np.max(spectralline)
            # p0=[localmem,l0,(np.amax(lvec)-np.amin(lvec))/4.]
            weights=spectralline/errors
            weights=spectralline/spectralline
            p=np.polyfit(lvec,np.log(spectralline),2,w=weights)
            sig=np.sqrt(-1./2/p[0])
            mu=-p[1]/2./p[0]
            a=np.exp(p[2]+mu**2/2/sig**2)
            p,pcov=curve_fit(fit_gauss,lvec,spectralline,p0=[a,mu,sig])
            a=p[0]
            mu=p[1]
            sig=p[2]
            peak[j,i]=a
            doppler[j,i]=(mu-l0)/l0*2.99e5 #in km/s
            sigma[j,i]=sig/l0*2.99e5 #in km/s
            intens[j,i]=a*sig*np.sqrt(2*np.pi)
            chisq[j,i]=np.sum((fit_gauss(lvec,a,mu,sig) - spectralline) ** 2)
        
    # p,pcov=curve_fit(fit_gauss,spectralline,lvec,p0=p0,sigma=np.full(nl,errors))
    # print a, mu, sigma
    # peak=p
    # chisq = sum((r / sigma) ** 2)
    
    return peak,doppler,sigma,chisq,intens
