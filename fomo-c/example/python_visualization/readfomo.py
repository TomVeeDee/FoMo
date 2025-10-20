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
    
def readgoftcube_dat(filename, compress=0):
    import gzip, struct, numpy as np

    if compress:
        f = gzip.open(filename, "rb")
    else:
        f = open(filename, "rb")

    # --- Try reading header until '#' (new format marker) ---
    header_bytes = b""
    while True:
        c = f.read(1)
        if not c:
            break
        header_bytes += c
        if c == b"#":
            break
        # Stop early if header gets too long
        if len(header_bytes) > 60:
            break

    header = header_bytes.decode(errors="ignore")

    # ------------------------------------------------------------------
    # Case 1: New format (FoMo-v...#)
    # ------------------------------------------------------------------
    if header.startswith("FoMo-") and header.strip().endswith("#"):
        dim = struct.unpack("i", f.read(4))[0]
        ng = struct.unpack("i", f.read(4))[0]
        nvars = struct.unpack("i", f.read(4))[0]

        # read units
        units = []
        for _ in range(dim + nvars):
            size_buf = f.read(8)
            if len(size_buf) < 8:
                raise EOFError("Unexpected EOF while reading unit size.")
            size = struct.unpack("q", size_buf)[0]
            unitstring = f.read(size) if size > 0 else b""
            units.append(unitstring.decode(errors="ignore"))

        chiantisize = struct.unpack("q", f.read(8))[0]
        chiantifile = f.read(chiantisize)
        abundsize = struct.unpack("q", f.read(8))[0]
        abundfile = f.read(abundsize)

    # ------------------------------------------------------------------
    # Case 2: Old format (no FoMo header)
    # ------------------------------------------------------------------
    else:
        # rewind to start
        f.seek(0)
        dim = struct.unpack("i", f.read(4))[0]
        ng = struct.unpack("i", f.read(4))[0]
        nvars = struct.unpack("i", f.read(4))[0]
        chiantisize = struct.unpack("q", f.read(8))[0]
        chiantifile = f.read(chiantisize)
        abundsize = struct.unpack("q", f.read(8))[0]
        abundfile = f.read(abundsize)
        units = []

    # ------------------------------------------------------------------
    # Common data read section
    # ------------------------------------------------------------------
    ncols = dim + nvars
    data = np.zeros((ng, ncols), dtype=np.float32)

    for j in range(ncols):
        buf = f.read(ng * 4)
        if len(buf) < ng * 4:
            raise EOFError("Unexpected EOF in data section.")
        data[:, j] = np.frombuffer(buf, dtype=np.float32, count=ng)

    f.close()
    return data, units, chiantifile

def readgoftcube(filename):
    parts = filename.split(".")
    compress = 0
    if parts[-1] == "gz":
        compress = 1

    if parts[-1 - compress] == "txt":
        data, chiantifile = readgoftcube_txt(filename, compress=compress)
    elif parts[-1 - compress] == "dat":
        data, _, chiantifile = readgoftcube_dat(filename, compress=compress)
    else:
        raise ValueError("Unknown file format: " + filename)

    return data
    
def regulargoftcube(data):
    print('data',data)
    ng = len(data)
    nx = round(ng/np.shape(np.where(data[:,0]==data[0,0]))[1])
    ny = round(ng/np.shape(np.where(data[:,1]==data[0,1]))[1])
    if (ng/nx/ny == 1):
	    nl=1
	    dim=2
    else:
            nl = round(ng/np.shape(np.where(data[:,2]==data[0,2]))[1])
            dim=3
#    print('dim',dim)
    xvec=np.empty(nx)
    yvec=np.empty(ny)
    emiss=np.empty([nx,ny])
    lvec=1.

    if (dim == 3):
            lvec=np.empty(nl)
            emiss=np.empty([nx,ny,nl])

    for i in range(ng):
        j=i/nl
        #begin modification by DS
        #changing values of j
        k=i/nl/nx
        k=int(k)
        j-=k*nx
        l=i%nl
        j=int(j)
        xvec[j]=data[i,0]
        yvec[k]=data[i,1]
        
        if (dim == 3):
                lvec[l]=data[i,2]
                emiss[j,k,l]=data[i,dim]
        if (dim == 2):
                emiss[j,k]=data[i,dim]
#                print('j',j,'k',k,'emis',emiss[j,k],'xvec',xvec[j],'yvec',yvec[k])
    
    # there's a very weird ordering of the indices, it seems x needs to come last, in order to have an intuitive access with matplotlib
    emiss=emiss[:,::-1]#DS change the yaxis 
    emiss=np.swapaxes(emiss,0,dim-1)
 #   print('emiss2',emiss)
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
            sigma[j,i]=sig
            intens[j,i]=a*sig*np.sqrt(2*np.pi)
            chisq[j,i]=np.sum((fit_gauss(lvec,a,mu,sig) - spectralline) ** 2)
        
    # p,pcov=curve_fit(fit_gauss,spectralline,lvec,p0=p0,sigma=np.full(nl,errors))
    # print a, mu, sigma
    # peak=p
    # chisq = sum((r / sigma) ** 2)
    
    return peak,doppler,sigma,chisq,intens
