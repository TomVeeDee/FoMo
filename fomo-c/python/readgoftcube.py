#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import gzip

def readgoftcube_txt(filename,compress=0):
    if (compress):
        f=gzip.open(filename,"r")
    else:
        f=open(filename,"r")
    textlines=f.readlines()[5:]
    data=[]
    for line in textlines:
        data.append(line.split())
    # now data contains
    # data[:][0] = x coordinate
    # data[:][1] = y coordinate
    # data[:][2] = wave coordinate
    # data[:][3] = intensity
    print len(data), len(data[:][0])
    
    return data
    
def readgoftcube(filename):
    parts=filename.split(".")
    compress=0
    if parts[-1]=="gz":
        compress=1
    
    if parts[-1-compress]=="txt":
        data=readgoftcube_txt(filename,compress=compress)
    
    if parts[-1-compress]=="dat":
        print "not implemented yet"
        # data=readgoftcube_dat(filename,compress=compress)
        
    return data