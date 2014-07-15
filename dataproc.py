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

from scipy.stats import binom
from scipy.stats import nbinom
from scipy.stats import norm
from scipy.stats import poisson
from scipy.stats import chisquare
import numpy as np
import matplotlib.pyplot as plt

#
# Data processing here
#


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='NAT data processor.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    parser.add_argument('--mean',           help='Graph main', required=False, default=False, action='store_true')
    parser.add_argument('file', action="store", nargs='+')
    args = parser.parse_args()
    
    keys  = []
    succ  = []
    mean  = []
    styles = ['--bx', '-.g2', ':.r', '--|k', ':m+', '--1c']
    
    for i, fname in enumerate(args.file):
        fh    = open(fname)
        dat   = fh.readlines()
        
        k, s, m = [], [], []
        for d in dat:
            d = str(d).strip()
            if d.startswith('#') or d.startswith('New'): continue
            
            arr = [float(x) for x in filter(None, d.split('|'))]
            if len(arr)==0: continue
            
            k.append(arr[0]) 
            s.append(arr[2])
            m.append(arr[3])
        keys.append(k)
        succ.append(s)
        mean.append(m)
        
        x = np.array(k)
        y = np.array(m if args.mean else s)
        
        tt = plt.plot(x, y, styles[i], label=chr(ord("A")+i))
    #plt.plot(xp1, pxp1, '--')
    #plt.plot(xp2, pxp2, 'g-')
    #plt.plot(xp3, pxp3, 'k-.')
    
    if args.mean: plt.legend(loc=1)
    else:         plt.legend(loc=3)
    
    plt.xlim(-0.01, max(max(keys)) * 1.1)
    if args.mean: pass #plt.ylim(0.0,max(y)*1.1)
    else:         plt.ylim(0.0,1.1)
    
    
    plt.xlabel('$\lambda$')
    plt.ylabel('Mean step success' if args.mean else 'success rate [%]') #,rotation='horizontal')
    plt.grid(True)
    plt.show()
        
        
        
        
        
        