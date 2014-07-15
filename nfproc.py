import os
import sys
import fileinput
import re
import random
import math
from operator import itemgetter, attrgetter
import subprocess
from optparse import OptionParser
import copy

import time
import argparse
from dateutil import parser as dparser
import calendar

import pylab as P
import numpy as np
from scipy.stats import binom
from scipy.stats import nbinom
from scipy.stats import norm
from scipy.stats import poisson
from scipy.stats import chisquare

# Multiple plots
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt

#
# Data processing here
#
def graph(plt, x='Sample', y='Port number', loc=1):
    if loc!=-1:
        plt.legend(loc=loc)
    plt.xlabel(x)
    plt.ylabel(y) #,rotation='horizontal')
    plt.grid(True)
    plt.show()
    plt.close()
    
def nfloat(x):
    if x=='nan': return 'nan'
    return float(x)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='NAT netflow data processor.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o','--output',    help='Output file name from finder', required=False, default='graph.txt')
    parser.add_argument('-t','--space',     help='Time in ms to wait between packet send', required=False, default=10, type=int)
    parser.add_argument('-l','--lmbd_start',help='On which lambda to start', required=False, default=-1, type=float)
    parser.add_argument('-s','--strategy',  help='Strategy to use (poisson, i2j, fibo, their, ij, binom, simple)', required=False, default='poisson')
    parser.add_argument('-r','--rounds',    help='Simulation rounds', required=False, type=int, default=1000)
    parser.add_argument('-e','--errors',    help='Maximum steps by algorithm', required=False, type=int, default=1000)
    parser.add_argument('-d','--dot',       help='Graphviz dot illustration', required=False, type=int, default=0)
    parser.add_argument('-a','--ascii',     help='Ascii illustration', required=False, type=int, default=0)
    parser.add_argument('-n','--nfdump',    help='NFdump file', required=False, default=None)
    parser.add_argument('-m','--nfdump_sorted',help='NFdump sorted file', required=False, default=None)
    parser.add_argument('-f','--filter',    help='NFdump filter', required=False, default=None)
    parser.add_argument('-g','--hostnet',   help='NFdump host address', required=False, default="147.250.")
    parser.add_argument('--lmbd',           help='Default Poisson lambda for simulations', required=False, type=float, default=0.1)
    parser.add_argument('--pval',           help='Graph pvalue', required=False, default=False, action='store_true')
    parser.add_argument('file', action="store", nargs='+')
    args = parser.parse_args()
    
    keys  = []
    succ  = []
    mean  = []
    styles = ['--bx', '-.g2', ':.r', '--|k', ':m+', '--1c']
    
    #for i, fname in enumerate(args.file):
    fname=args.file[0]
    fh    = open(fname)
    dat   = fh.readlines()
        
    n    = 0
    last = -1
    
    kk = [] # keys
    ex = [] # E[X]
    vx = [] # V[X]
    ss = [] # ssum
    
    pk = [[] for i in range(0,6)] # p-value, key
    pv = [[] for i in range(0,6)] # p-value
    ck = [[] for i in range(0,6)] # chi-square, key
    cv = [[] for i in range(0,6)] # chi-square
    
    for d in dat:
        d = str(d).strip()
        if d.startswith('#') or d.startswith('New'): continue
        arr = [nfloat(x) for x in filter(None, d.split('|'))]
        if len(arr)==0: continue
        
        # Process the file, E[X], V[X], SUM
        if arr[0] > n: n=int(arr[0])
        if last != arr[0]:
            kk.append(arr[0])
            ex.append(arr[2])
            vx.append(arr[3])
            ss.append(arr[4])
            
        idx = int(arr[1])
        
        # Distribution
        if arr[5] != 'nan': 
            pk[idx].append(arr[0])
            ck[idx].append(arr[0])
            pv[idx].append(arr[5])
            cv[idx].append(arr[6])
        
        # last
        last = arr[0]
        
    # Process output to nicely looking graph
    x = np.array(range(0, n))
    
    # hypothesis tests results
    hypo_0 = len(filter(lambda x: x >= 0.05, pv[0]))
    hypo_4 = len(filter(lambda x: x >= 0.05, pv[4]))
    
    print "Statistical data"
    print "Mean EX %03.4f; Median EX %03.4f; V[Mean] %03.4f; Mean VX %03.4f; Median VX %03.4f;" % (np.mean(ex), np.median(ex), np.var(ex), np.mean(vx), np.median(vx))
    
    
    print "Hypothesis testing result"
    print "Poisson: %01.5f; median p-value: %01.8f; median chi-square: %01.8f" % (hypo_0/float(len(pk[0])), np.median(pv[0]), np.median(cv[0]))
    print "NBinom:  %01.5f; median p-value: %01.8f; median chi-square: %01.8f" % (hypo_4/float(len(pk[4])), np.median(pv[4]), np.median(cv[4]))
    
    print "%01.3f & %01.3f & %03.4f & %03.4f & %03.4f" % ((1 - hypo_0/float(len(pk[0]))) * 100, (1 - hypo_4/float(len(pk[4]))) * 100, np.mean(ex), np.var(ex), np.mean(vx))
    
    # e,x
    ex_np = np.array(ex)
    vx_np = np.array(vx)
    exvx  = np.abs(vx_np / ex_np)
    plt.plot(x, ex_np, 'b+', label="E[X]")
    #plt.plot(x, vx_np, 'r3', label="V[X]")
    graph(plt)
    
    #
    # Histogram for E[X]
    #
    P.grid(True)
    P.Figure()
    P.hist(ex_np, 20, normed=0, histtype='bar')
    P.legend()
    graph(plt, x='$E[X]$', y='Count', loc=-1)
    
    #
    # Histogram for dispersion
    #
     
    #P.grid(True)
    #P.Figure()
    #P.hist(exvx, 20, normed=1, histtype='bar')
    #P.legend()
    #graph(plt, x='$E[X] / V[X]$', y='Count', loc=-1)
    
    
    #
    # P-value histogram
    #
    
    # p-value with critical region
    pk_p = np.array(pk[0]) # poisson, key
    pk_n = np.array(pk[4]) # nbin, key
    pv_p = np.array(pv[0]) # poisson, value
    pv_n = np.array(pv[4]) # nbin, value
    
    # critical region
    plt.axvspan(0.0, 0.05, facecolor='r', alpha=0.6)
    
    # Stacked histogram for poisson and negative binomial 
    P.grid(True)
    P.Figure()
    P.hist([pv_p, pv_n], 20, normed=0, histtype='bar', color=['green', 'blue'], label=['Po', 'NB'])
    P.legend()
    graph(plt, x='p-value', y='Count', loc=-1)
    
    #
    # P-value scatter plot
    #
    plt.plot(pk_p, pv_p, 'g+', label="Po")
    plt.plot(pk_n, pv_n, 'b3', label="NB")
    plt.axhspan(0.0, 0.05, facecolor='r', alpha=0.5) # p-value reqion
    graph(plt, y='p-value', loc=-1)
    
          