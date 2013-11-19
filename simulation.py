#
# NAT protocol simulation
#
#
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
from scipy.stats import poisson
import numpy as np
import matplotlib.pyplot as plt
import time
import argparse

# Fibonacchi list generator
def fibGenerator():
    a, b = 0, 1
    yield 0
    while True:
        a, b = b, a + b
        yield a
            
def f7(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]

class Strategy:
    '''
    Base abstract class for NAT strategy
    '''
    def init(self, params=None):
        raise Exception("Not implemented yet...")
    def reset(self, nats=[], sim=None, params=None):
        self.init(None)
    def next(self, party, step):
        raise Exception("Not implemented yet...")   # return (srcPort, dstPort)
    def silent(self, time1, time2, lmbd):
        '''
        Tells to strategy how long did silent period take just before start and guessed lambda.
        '''
        pass

class Nat:
    '''
    Base abstract class for NAT allocation
    '''
    
    # timeout of a created connection in milli seconds
    timeout = 3*60*1000
    
    # port pool available for new allocations
    pool = None
    poolLen = 0
    
    def init(self, params=None):
        raise Exception("Not implemented yet...") 
    def reset(self):
        raise Exception("Not implemented yet...")
    def alloc(self, srcIP, srcPort, dstIP, dstPort, timeNow):
        raise Exception("Not implemented yet...") 
    def occupy(self, num, timeNow):
        raise Exception("Not implemented yet...")
    def freePorts(self):
        raise Exception("Not implemented yet...")

class Quartet:
    '''
    SrcIP, srcPort, DstIP, DstPort
    '''
    srcIP=0
    srcPort=0
    dstIP=0
    dstPort=0
    
    def __init__(self, srcIP, srcPort, dstIP, dstPort):
        self.srcIP = srcIP
        self.srcPort = srcPort
        self.dstIP = dstIP
        self.dstPort = dstPort
    def __cmp__(self, other):
        if self.srcIP == other.srcIP          \
            and self.srcPort == other.srcPort \
            and self.dstIP == other.dstIP     \
            and self.dstPort == other.dstPort: return 0
        else: return 1
    def __eq__(self, other):
        return self.__cmp__(other) == 0
    def __ne__(self, other):
        return self.__cmp__(other) != 0
    def __str__(self):
        return "%s:%05d --> %s:%05d" % (self.srcIP, self.srcPort, self.dstIP, self.dstPort)
    def __hash__(self):
        prime=31
        result=1
        result = prime * result + self.srcIP
        result = prime * result + self.srcPort
        result = prime * result + self.dstIP
        result = prime * result + self.dstPort
        return result
    
class SymmetricNat(Nat):
    '''
    Base class for symmetric NAT. 
    '''
    # allocation table for NAT; key = quartet; value = external port
    allocations = None 
    # port -> (quartet, expire time). Quartet may be null
    allocatedPorts = None
    
    # port is the key
    # port -> quartet, expire
    # quartet -> port
    
    def init(self, params=None):
        self.allocations = {}
        self.pool = range(0, 65536)
        self.poolLen = len(self.pool)
        self.allocatedPorts = {}
        
    def reset(self):
        self.allocations = {}
        self.allocatedPorts = {}
    
    def nextPort(self):
        '''
        Uses port pool array and pointer to last allocated port to obtain next in the sequence.
        In case of random allocation randomly generates index to a pool and returns a port on the index.
        '''
        raise Exception("Not implemented yet... This class is abstract, you have to override this method in subclass")
    
    def nextFreePort(self, timeNow):
        '''
        Returns next free port in the sequence, takes existing associations into account and their expiration
        '''
        tries=0                                                # pool exhaustion check
        port=-1
        while tries <= self.poolLen:
            port   = self.nextPort()                           # linear port allocation rule here
            tries += 1                                         # check pool exhaustion
            if port in self.allocatedPorts:
                # next port is already allocated, what about timeout?
                tup = self.allocatedPorts[port]
                # check expiration first
                if (tup[1] + self.timeout) < timeNow:
                    if (tup[0] != None):
                            del self.allocations[tup[0]]       # expired -> delete from allocation table
                    del self.allocatedPorts[port]              # delete from allocation set
                    break                                      # slot is free now, can allocate
                else: continue                       # slot is in use, continue with search
            else: break                              # slot is free, assign
        # check if pool is exhausted - all ports are allocated currently
        if tries >= self.poolLen or port==-1:
            raise Exception("Port pool exhausted")
        # return resulting port, should not be -1
        return port
    
    def alloc(self, srcIP, srcPort, dstIP, dstPort, timeNow):
        '''
        Basic allocation method for new connection
        '''
        q = Quartet(srcIP, srcPort, dstIP, dstPort)
        
        # Check for existing allocation for a given quartet
        if q in self.allocations:
            port = self.allocations.get(q)
            # check expiration time, if a record is too old, it has to be removed from alloc table.
            tup = self.allocatedPorts.get(port)
            if (tup[1] + self.timeout) < timeNow:
                del self.allocatedPorts[port]     # delete from allocation set
                del self.allocations[q]           # expired -> delete from allocation table
            else:
                self.allocatedPorts[port] = (q, timeNow)    # update last query access
                return port                                 # external port returned
        
        # if here -> not in allocation list, create a new allocation
        port=self.nextFreePort(timeNow)
        # create a new allocation
        self.allocatedPorts[port] = (q, timeNow)
        self.allocations[q] = port
        
        return port
    
    def occupy(self, num, timeNow):
        '''
        Simulates another connections created randomly
        '''
        for i in range(0, num):
            port = self.nextFreePort(timeNow)
            self.allocatedPorts[port] = (None, timeNow)
        return 1
    
    def freePorts(self):
        return (self.poolLen - len(self.allocatedPorts))
    
    def trulyFreePorts(self, timeNow):
        cp = copy.deepcopy(self.allocatedPorts)
        for port in cp:
            # next port is already allocated, what about timeout?
            tup = self.allocatedPorts[port]
            # check expiration first
            if (tup[1] + self.timeout) < timeNow:
                if (tup[0] != None):
                        del self.allocations[tup[0]]       # expired -> delete from allocation table
                del self.allocatedPorts[port]              # delete from allocation set
        return self.freePorts()

class SymmetricRandomNat(SymmetricNat):
    '''
    Symmetric NAT with random allocation function 
    '''
    def nextPort(self):
        '''
        Randomly generates index to a pool and returns a port on the index.
        '''
        return self.pool[random.randint(0, self.poolLen-1)]
    
class SymmetricIncrementalNat(SymmetricNat):
    # index of last allocated port. Index to pool[]
    lastPort = 0
    
    def reset(self):
        self.allocatedPorts={}
        self.allocations={}
        self.lastPort = 0
        
    def nextPort(self):
        '''
        Uses port pool array and pointer to last allocated port to obtain next in the sequence.
        ''' 
        self.lastPort = (self.lastPort + 1) % self.poolLen # linear port allocation rule here
        return self.pool[self.lastPort]                  # just a shortcut

class TheirStragegy(Strategy):
    delta = [1000,1000]
    def init(self, params=None):
        pass    
    def next(self, party, step):
        if party==0: return (step,self.delta[0])
        if party==1: return (step,self.delta[1])
    
class I2JStragegy(Strategy):
    startPos=[0,0]
    def init(self, params=None):
        pass
    def silent(self,  time1, time2, lmbd):
        # Sets different starting point for algorithm than zero. Takes silent period
        # duration into account together with predicted workload to start algorithm
        # on a right place. It there are too many errors in the prediction (new connections)
        # small-step-big-step can have problems with cathing them up.
        #
        # Use expected value instead of a random sample as a starting point. E(X) = lmbd, X ~ Po(lmbd)
        # should be the central.
        #self.startPos=[int(lmbd * time1), int(lmbd * time2)]
        #self.startPos=[NatSimulation.poisson(lmbd, time1), NatSimulation.poisson(lmbd, time2)]
        return self.startPos
        
    def next(self, party, step):
        if party==0: return (0, int(self.startPos[0]+step))
        if party==1: return (0, int(self.startPos[1]+2*step))
        
class FiboStrategy(Strategy):
    fibn = []
    b    = []
    startPos=[0,0]
    def init(self, params=None):
        self.fibn = []
        fib = fibGenerator()
        for n in range(22):
            self.fibn.append(next(fib))        
        for i in range(1, len(self.fibn)-1):
            for j in range(0, self.fibn[i-1]):
                #int(NatSimulation.poisson(0.1, 10 * (1+self.fibn[i+1] + j)    ))
                self.b.append(self.fibn[i+1] + j)
        #sys.exit(1)
    
    def silent(self,  time1, time2, lmbd):
        #self.startPos=[int(lmbd * time1), int(lmbd * time2)]
        #self.startPos=[NatSimulation.poisson(lmbd, time1), NatSimulation.poisson(lmbd, time2)]
        return self.startPos
    def next(self, party, step):
        if party==0: return (0, self.startPos[party] +   self.b[step])
        if party==1: return (0, self.startPos[party] + 2*self.b[step])
        return

class PoissonStrategy(Strategy):
    startPos=[0,0]
    nats = None
    sim  = None
    lmbd = 0.1
    dupl = False
    coef = 1.4773
    
    b = [[],[]]
    def init(self, params=None): 
        #self.reset()
        self.gen()
        #print self.b[0], "len=", len(self.b[0])
    
    def reset(self, nats=[], sim=None, params=[]):
        self.sim = sim
        
        if len(nats)==2: self.nats = nats
        if self.sim!=None: self.lmbd = sim.lmbd
        
        self.gen()
        pass
    
    def genPart(self, party):
        # lambda on both sides
        lmbd = self.sim.lmbd if self.sim!=None else self.lmbd
        
        # port scan interval from simulation
        t = self.sim.portScanInterval if self.sim != None else 10
        
        b    = []             # local array
        seen = set()          # duplicity check
        seen_add = seen.add
        
        for step in range(0, 3001):     
            #self.b[0].append(int(NatSimulation.poisson(self.lmbd, 10 * (1+step*1.875)    )))
            #self.b[1].append(int(NatSimulation.poisson(self.lmbd, 10 * (1+step*1.875)    )))
            
            #self.b[0].append(int(NatSimulation.poisson(self.lmbd, 10 * (1+step*1.778)    )))
            #self.b[1].append(int(NatSimulation.poisson(self.lmbd, 10 * (1+step*1.778)    )))
            
            # dupl = False, lmbd=0.01, t=10
            #x = int(NatSimulation.poisson(lmbd, t * (1+step*4.5)    ))
            
            # dupl = False, lmbd=0.01, t=10
            #x = int(NatSimulation.poisson(self.lmbd, 10 * (1+step*5.5)    ))
            
            # dupl=True, lmbd=0.1, t=10
            #x = int(NatSimulation.poisson(self.lmbd, 10 * (1+step*1.778)    ))
            
            # dupl=False, lmbd=0.1, t=10
            #x = int(NatSimulation.poisson(self.lmbd, t * (1+step*1.485)    ))
            
            x  = round(step * 2.03) 
            #x = int(NatSimulation.poisson(lmbd, t * (1+step*self.coef)    ))
            if self.dupl or (x not in seen and not seen_add(x)): 
                b.append(x)
        return b
                
    def gen(self):
        self.b = [self.genPart(0), self.genPart(1)]
            
        #self.b[0] = list(set(self.b[0]))
        #self.b[1] = list(set(self.b[1]))
         
        #self.b[0] = f7(self.b[0])
        #self.b[1] = f7(self.b[1])
        
    def silent(self,  time1, time2, lmbd):
        self.lmbd = lmbd
        self.startPos=[int(lmbd * time1), int(lmbd * time2)] # expected value
        #self.startPos=[NatSimulation.poisson(lmbd, time1), NatSimulation.poisson(lmbd, time2)]        
        return self.startPos
    
    def next(self, party, step):
        #return (0, int(self.startPos[party] + NatSimulation.poisson(self.lmbd, 10 * (1+step*1.77)    )))
        return (0, int(self.startPos[party] + self.b[party][min(step, len(self.b[party])-1)]))
        
        #self.startPos[party] += 1+NatSimulation.poisson(self.lmbd, 10)#*(1+step*0.77))
        #return (0, int(self.startPos[party]))

class NatSimulation:
    
    # Lambda for Poisson process generator. Time unit = 1 ms
    # Intuitively: represents rate of creation of new events in given period.
    #  Average number of arrivals per unit time is lambda.        [analysis&synthesis]
    #  The expected length of interarrival intervals is 1/lambda. [analysis&synthesis]
    #  Interarrival intervals are independent and distributed exponentially with parameter lambda. [analysis&synthesis]
    #
    # @see http://www.columbia.edu/~ks20/4703-Sigman/4703-07-Notes-PP-NSPP.pdf
    # @see http://filebox.vt.edu/users/pasupath/papers/poisson_streams.pdf
    # @see http://www.math.wsu.edu/faculty/genz/416/lect/l05-45.pdf
    # @see http://preshing.com/20111007/how-to-generate-random-timings-for-a-poisson-process/
    lmbd = 0.1
    
    # Number of miliseconds for silent period to take [ms].
    # Based on basic ping / round trip time it takes to communicate 
    # IP with another peer 
    silentPeriodBase = 0#1000
    
    # Lambda for Pois(lmbd) for silent period variability.
    # Silent period time = silentPeriodBase + Pois(lmbd) [ms]
    silentPeriodlmbd = 0#100
    
    # Number of rounds for simulation
    simulationRounds = 1000
    
    # how many rounds has fast simulation? Fast is extended to deep simulation if it has
    # more than 50% probability of success
    simulationRoundsFast = 100  
    
    # Number of errors that are handled by algorithm 
    errors = 1000
    
    # number of milliseconds between consecutive port scans
    portScanInterval = 10
    
    # number of connections to establish
    numCon = 1
    
    # very compact output
    compact=True
    
    @staticmethod
    def poisson(lmbd, t):
        '''
        Uses Numpy package to take sample from poisson distribution
        '''
        return int(np.random.poisson(lmbd*t, 1)[0])
    
    @staticmethod
    def uniform(lmbd, t):
        return random.random() * lmbd * t
    
    @staticmethod
    def poissonSample(lmbd, t):
        '''
        Generates number of events in Poisson process in time [0, t]
        source: http://www.math.wsu.edu/faculty/genz/416/lect/l05-45.pdf
        '''
        u = random.random()
        N = 0
        p = math.exp(-lmbd * t)
        F = p
        while u > F:
            N = N+1
            p = lmbd*t*p/N 
            F = F + p
        return N
    
    def poissonCDF(self, lmbd, x):
        '''
        Returns P(X <= x), X~Poisson(lmbd)
        
        P(X <= x) = e^{-lambda} * \sum_{i=0}^{k}{ \frac{lambda^i}{i!} }
        '''
        F = 1
        res = 0
        for i in range(0, x+1):
            res = res + (math.pow(lmbd, i) / F)
            F = F * (i+1)
        return (math.exp(-lmbd) * res) 
    
    def getNumOfNewConnections(self, tim):
        '''
        Simple wrapper for poission. Returns number of new connections
        created. It is assumed they are distributed according to Poisson distribution.
        '''
        return int(np.random.poisson(self.lmbd*tim, 1)[0])
    
    def poissonSimulate(self, T):
        '''
        Simulates Poisson process with arrival times
        source: http://www.columbia.edu/~ks20/4703-Sigman/4703-07-Notes-PP-NSPP.pdf
        '''
        t = 0.0
        N = 0
        while t <= T:
            # U ~ U(0,1), uniform distribution
            U = random.random()
            
            # next time of the event, exponential distribution
            t = t + (-(1/self.lmbd) * math.log(U))
            if (t > T): return N
        
            # increment the event counter
            N = N + 1
            print "New event will occur: " + str(t) + ("; now events: %02d" % N) 
        
        return N
  
    def coefFinder(self, natA, natB, strategy, baseStep=0.10, start=0.1, maxc=15.0, epsdiff=0.25, maxp_acc=0.1, depth=0):
        '''
        Finds coefficient that maximizes probability of establishing connection for Poisson strategy
        '''
        
        probs = {}
        curc = start
        
        maxp   = 0.0
        maxcur = curc
        resm   = None
        
        minp   = 1.0
        mincur = curc
        
        bpoint = curc
        
        # Scanning with small step
        # TODO: if scanning step is too big and it is not possible to find the peak
        # in probability function step has to be smaller and interval re-scanned to
        # find desired maximum. 
        while curc < maxc:
            strategy.coef = curc    # sets current coefficient to strategy
            
            sys.stdout.write(" curc=%03.4f; " % curc)
            res = self.simulation(natA, natB, strategy)
            probs[curc] = res[0]
            
            if res[0] > maxp:
                maxp = res[0]
                maxcur = curc
                resm   = res
                
            if res[0] < minp:
                minp = res[0]
                mincur = curc
            
            # check if we are performing worse than before
            if maxp >= maxp_acc and maxp > res[0] and abs(maxp-res[0]) >= epsdiff:
                print "; We are getting worse...; "
                bpoint = curc
                break
            # too good solution
            if maxp >= 0.99: break
            
            curc += baseStep    # increment current coefficient in strategy to the next round        
        print probs
        
        # Coefficient binary finding, maximum should be somewhere in the middle
        #print self.coefFinderInterval(natA, natB, strategy, bpoint-3*baseStep, bpoint, 0)
        
        # recursive call on this function - try finer step
        if maxp > 0.99 or depth >= 3:
            print "Ending recursion; max=%03.4f bp=%03.4f" % (maxp, maxcur)
            return (maxp, curc, resm[2] if resm != None else '0') 
        else:
            return self.coefFinder(natA, natB, strategy, baseStep/10.0, maxcur-2*baseStep, maxcur+2*baseStep, epsdiff/1.0, maxp, depth+1)
        
        pass
  
    def coefFinderInterval(self, natA, natB, strategy, cl, cr, step=0, prec=0.001):
        '''
        Recursive binary search for finding coefficient that maximizes probability of successful connection establishment
        ''' 
        eps = 100 # search epsilon, accuracy, 100ms
        t   = 0
        
        stepNull = step==0
        cc       = 0
        while cl < cr and (cr - cl) > 0.0001:
            if stepNull: step=(cr-cl) / 20.0
            cc = cl + (cr-cl) / 2.0

            # mid-1
            strategy.coef = cc-step
            sys.stdout.write("c-1 [%02.03f, %02.03f] curc=%03.4f; " % (cl, cr, strategy.coef))
            res_cm = self.simulation(natA, natB, strategy)
            
            # mid
            strategy.coef = cc
            sys.stdout.write("c   [%02.03f, %02.03f] curc=%03.4f; " % (cl, cr, strategy.coef))
            res_c = self.simulation(natA, natB, strategy)
            
            #mid+1
            strategy.coef = cc+step
            sys.stdout.write("c+1 [%02.03f, %02.03f] curc=%03.4f; " % (cl, cr, strategy.coef))
            res_cp = self.simulation(natA, natB, strategy)
            
            print ""
            sys.stdout.flush()
            
            # decision which path to take
            if res_cm[0] >= 0.99:
                return cc-step
            elif res_cp[0] >= 0.99:
                return cc+step
            elif res_c[0] >= 0.99:
                return cc
            elif res_cm[0] <= res_c[0] and res_c[0] >= res_cp[0]: # middle is peak; could return already but we might obtain better peak by "zooming"
                print "shrinking interval"
                cl = cl + (cr-cl) / 4.0
                cr = cr - (cr-cl) / 4.0
            elif res_cm[0] >= res_c[0] and res_c[0] <= res_cp[0]: # middle is low-peak;
                if res_cm[0] >= res_cp[0]:   # if left side is bigger
                    print "lowpeak, going left..."
                    cr = cc-step
                else:                       # right side is bigger
                    print "lowpeak, going right..." 
                    cl = cc+step
            elif res_cm[0] >= res_c[0]:
                print "going left..."
                cr = cc-step
            elif res_cp[0] >= res_c[0]:
                print "going right..."
                cl = cc+step
            else:
                print "dafuq?"
                return cc   
        return cc
    
    def simulation(self, natA, natB, strategy):
        '''
        Simple simulation of NAT traversal algorithm.
        '''
        
        nats = [natA, natB]
        successCnt = 0.0
        stopOnFirstMatch = self.simulationRounds != 1
        getTime = lambda: int(round(time.time() * 1000))
        simStart = getTime()
        
        successAcc = [0,0]              # accumulator for steps needed to connect if successfully
        realRounds = self.simulationRounds
        for sn in range(0, self.simulationRounds):
            # reset NATs
            nats[0].reset()
            nats[1].reset()
            strategy.reset(nats, self)
            
            # generate silent period time
            curSilentA = self.silentPeriodBase + self.poisson(self.silentPeriodlmbd, 1)
            curSilentB = self.silentPeriodBase + self.poisson(self.silentPeriodlmbd, 1)
            
            if not self.compact:
                print "\n##%03d. Current silent period time: [%03.3f, %03.3f]" % (sn, curSilentA, curSilentB) 
            
            # generate new TCP connections for silent period on both sides, same lambda.
            kA = self.poisson(self.lmbd, curSilentA)
            kB = self.poisson(self.lmbd, curSilentB)
            
            # reflect errors to NAT allocation
            nats[0].occupy(kA, 0)
            nats[1].occupy(kB, 0)
            
            # assume we are always starting from port 0 on both sides
            #i = kA + 2*kB
            #j = kA +   kB
            #print "Offsets that would match without any further errors i=%02d; j=%02d" % (i, j)
            
            # now simulate the protocol, phase with port scanning
            #targetA = 2 * self.numCon * self.errors
            #targetB =     self.numCon * self.errors
            
            # set silent period duration to the strategy
            sData = []
            sData = strategy.silent(curSilentA, curSilentB, self.lmbd)
            
            mapA  = [{}, {}]                # mapping of the current port to index
            scanA = [set([]), set([])]      # list of a tuple (assigned port, destination port)
            portsA = [set([]), set([])]     # set of an allocated ports
            totalLagA = [0, 0]              # total number of errors during protocol
            foundSomething = False
            for i in range(0, self.errors):
                
                # A scan
                #dstA  = b[i]#1*i #- stageChange*(stageNumA)/10.0# destination of scan o the other side
                for party in [0,1]:
                    # Obtain next tuple (source port, destination port) from strategy
                    nextA = strategy.next(party, i)
                    dstA  = nextA[1]
                    # Obtain external NAT port by querying NAT for allocation a new connection
                    curA  = nats[party].alloc(party, nextA[0], party ^ 0x1, dstA, i*self.portScanInterval)
                    
                    # Waiting between consecutive scans, compute number of new connections by 
                    # using Poisson process. Now generating new allocations to the new round/step of the protocol.
                    curLag = self.getNumOfNewConnections(self.portScanInterval)
                    totalLagA[party] += curLag
                    
                    # Reflect allocations meanwhile to the NAT
                    nats[party].occupy(curLag, i*self.portScanInterval)
                    
                    # Add protocol to the maps.
                    toAdd  = (curA, dstA) if party==0 else (dstA, curA)     # swap pair for other party - in order to find set intersection
                    scanA[party].add(toAdd)
                    portsA[party].add(curA)
                    mapA[party][curA] = i
                    #print "A scan: %d [%03d] --> [%03d] lag=%02d i=%03d toAdd=%s" % (party, curA, dstA, curLag, i, str(toAdd))
                    
                    if stopOnFirstMatch and toAdd in scanA[party ^ 0x1]: 
                        foundSomething = True
                if foundSomething: break
            
            if not self.compact:
                print "New connections meanwhile silent period [%02d, %02d]; " % (kA, kB), "start=", sData, "; totalLags [%02d %02d]" % (totalLagA[0], totalLagA[1])
            
            # OK is there any intersection in both sets?
            res = list(scanA[0].intersection(scanA[1]))
            # sort by minimum element in tuple
            res.sort(key=lambda tup: min(tup[0], tup[1]))
            
            # Generate DOT graph
            if self.simulationRounds==1:
                self.generateDot(portsA[0], portsA[1], scanA[0], scanA[1], res)
            
            if (len(res) == 0): 
                if not self.compact:
                    print "Algorithm failed, no intersecting points"
                else:
                    sys.stdout.write('.')
                    sys.stdout.flush()
                    
            # Stop early if poor performance
            if False and sn == self.simulationRoundsFast and 2.0*successCnt < self.simulationRoundsFast:
                sys.stdout.write('Z')
                sys.stdout.flush()
                realRounds = sn+1
                break
            
            # fail -> nothing to do now
            if (len(res) == 0): 
                continue
            
            if not self.compact: 
                print "RES: ", res, "i=%02d" % mapA[0][res[0][0]], "; j=%02d" % mapA[1][res[0][1]]
            else:
                sys.stdout.write('X')
                sys.stdout.flush()
            
            successCnt += 1.0
            successAcc[0] += mapA[0][res[0][0]]
            successAcc[1] += mapA[1][res[0][1]]
        
        simEnd = getTime()
        simTotal = simEnd - simStart    
            
        # Report results after simulation is done
        print "\nSuccess count: %02.3f ; cnt=%03d; lmbd=%01.3f; scanInterval=%04d ms; base sleep=%04d; average steps: %04.3f %04.3f; time elapsed=%04.3f s" % \
            (successCnt / realRounds    if realRounds > 0 else 0, 
             successCnt, 
             self.lmbd, 
             self.portScanInterval, 
             self.silentPeriodBase,
             successAcc[0] / successCnt if successCnt > 0 else 0,
             successAcc[1] / successCnt if successCnt > 0 else 0,
             simTotal/1000.0)
        
        return (successCnt / realRounds    if realRounds > 0 else 0, 
                successCnt, 
                successAcc[0] / successCnt if successCnt > 0 else 0,
                successAcc[1] / successCnt if successCnt > 0 else 0,)
    
    def generateDot(self, portsA, portsB, scanA, scanB, res):
        '''
        Generate DOT image for protocol run.
        '''
        
        dot = "digraph finite_state_machine {\n"
        maxport = max(max(portsA), max(portsB))
        for p in range(0, maxport):
            ps = str(p)
            inA = len([i for i in res if i[0] == p]) > 0
            inB = len([i for i in res if i[1] == p]) > 0
            
            line = "node [shape=circle, fixedsize=true, width=1, height=1, style=filled, colorscheme=orrd9, fillcolor=\"%s\" pos=\"%f,%f!\" label=\"%s\"] P%s;\n"
            desc = "node [shape=plaintext, width=2, pos=\"%f,%f!\" label=\"%s\"] DSC%s;\n"
            
            # A
            whoping  = [i for i in scanB if i[0] == p]
            
            color = 7
            if p in portsA:
                color = 3 if len(whoping) > 0 else 1
            if inA:
                color = "#00ff005f"
            dot = dot + line % (color, 10, 1.5*p, p, "A"+ps)
            
            # A desc
            if p in portsA:
                bcounter = [i for i in scanA if i[0] == p][0][1]
                acounter = [i for i in scanB if i[1] == bcounter] if bcounter in portsB else []
                acounter = acounter[0][0] if len(acounter)>0 else -1
                dot = dot + desc % (8, 1.5*p, "%03d -> %03d -> %03d\\n%s" % (p, bcounter, acounter, str(whoping)), "A"+ps)
            else:
                dot = dot + desc % (8, 1.5*p, "%s" % (str(whoping)), "A"+ps)
                
            
            # B
            whoping  = [i for i in scanA if i[1] == p]
            
            color = 7
            if p in portsB:
                color = 3 if len(whoping) > 0 else 1
            if inB:
                color = "#00ff005f"
            dot = dot + line % (color, 130, 1.5*p, p, "B"+ps)
            
            # B desc
            if p in portsB:
                acounter = [i for i in scanB if i[1] == p]
                acounter = acounter[0][0] if (len(acounter)>0) else -1
                bcounter = [i for i in scanA if i[0] == acounter][0][1] if acounter in portsA else -1
                dot = dot + desc % (132, 1.5*p, "%03d -> %03d -> %03d\\n%s" % (p, acounter, bcounter, str(whoping)), "B"+ps)
            else:
                dot = dot + desc % (132, 1.5*p, "%s" % (str(whoping)), "B"+ps)
        dot = dot + "\n\n"
        
        # add connections representing scan
        extraArrow =  "[penwidth=\"3.0\", arrowsize=\"2.5\"]"
        for tup in scanA:
            dot = dot + "PA%d -> PB%d %s\n" % (tup[0], tup[1], extraArrow if tup in res else "")
        for tup in scanB:
            dot = dot + "PB%d -> PA%d %s\n" % (tup[1], tup[0], extraArrow if tup in res else "")
        
        # generate graphviz image only for 1 round - illustrative run only
        dot = dot + "fontsize=32;}"
        f = open('dotfile.dot', 'w')
        f.write(dot)
        f.close()
        
        # generate SVG file
        print "GraphViz output: ", subprocess.Popen('neato -Tsvg < dotfile.dot > dotfile.svg', shell=True).communicate()[0]
    
    def poolExhaustionNat(self, natA, timeout):
        return self.poolExhaustion(timeout, natA.poolLen, self.lmbd)
    
    def poolExhaustion(self, timeout, poolsize, lmbd):
        '''
        Computes how long does it take to exhaust port pool given the new connection creation rate
        
        Related:
            Simulates Poisson process with arrival times
            source: http://www.columbia.edu/~ks20/4703-Sigman/4703-07-Notes-PP-NSPP.pdf
        '''
        t = 0.0
        N = 0
        while N <= poolsize:
            # U ~ U(0,1), uniform distribution
            U = random.random()
            
            # next time of the event, exponential distribution
            t = t + (-(1/lmbd) * math.log(U))
            if (N > poolsize): return t
        
            # increment the event counter
            N = N + 1
            #print "New event will occur: " + str(t) + ("; now events: %02d" % N)
        print "Port pool will be exhausted in %05.3f ms = %05.3f s = %05.3f min = %05.3f h" % (t, t/1000.0, t/1000.0/60, t/1000.0/60/60)
        print "P(X > portPoolSize) = %02.3f where X~Poisson(timeout * lamda)" % (1.0-poisson.cdf(poolsize, lmbd * timeout))
        return t 
    
    def poolExhaustionEx(self, natA, timeout):
        '''
        Computes a simulation of port pool exhaustion of a NAT taking into consideration timeout.
        It is the same thing like poolExhaustion in this way: if result from poolExhaustion is below timeout of the NAT, 
        it will be exhausted with given setting. Otherwise some ports will timeout and NAT will never exhaust its pool size.
        
        Corresponds to the sample of P(X >= poolSize), X~Poisson(timeout * lambda). X is number of new connections in 
        the given timeout interval. This gives us probability that NAT will be exceeded.
        '''
        t = 0.0
        N = 0
        try:
            while True:
                # U ~ U(0,1), uniform distribution
                U = random.random()
                
                # next time of the event, exponential distribution
                nextEvt = (-(1/self.lmbd) * math.log(U)) 
                t = t + nextEvt 
                
                # add new port allocation by that time
                natA.occupy(1, t)
            
                # increment the event counter
                N = N + 1
                
                if (N % 10000) == 0:
                    freePorts = natA.trulyFreePorts(t)
                    freeRatio = freePorts / float(natA.poolLen)
                    print "New event will occur: " + str(t) + ("; now events: %02d; evt=%05.3f; freePorts=%d, %02.2f %%" % (N, nextEvt, freePorts, freeRatio))
            
        except Exception, e:
            print "Port pool will be exhausted in %05.3f ms = %05.3f s = %05.3f min = %05.3f h" % (t, t/1000.0, t/1000.0/60, t/1000.0/60/60)
            print "Exception: ", e
            pass
        
        return 0
    
    def getLambdaExhaustionCDF(self, natA, prob):
        '''
        Gets lambda such that:
        
        P(X > poolsize) >= prob, X ~ Poisson(lambda * timeout)
        '''
        # at first we have to find proper interval where to find by binary search
        timeout  = natA.timeout
        poolsize = natA.poolLen
        
        lmbd = 1.0
        lmbdL=-1
        lmbdR=-1
        
        while (lmbdL==-1) or (lmbdR==-1):
            probc = 1.0 - poisson.cdf(poolsize, lmbd * timeout)
            print "current lambda: %02.3f ; prob=%02.3f" % (lmbd, probc)
            
            # left side of the interval. If fits, set and continue to find right side, otherwise keep left 
            if lmbdL==-1:
                if probc <= prob: 
                    lmbdL = lmbd
                    lmbd  = lmbd * 2.0
                    continue
                else: 
                    lmbd = lmbd/2.0
                    continue
            
            # right side of the interval, if here, we have left side and finding the right side
            if probc > prob: 
                lmbdR = lmbd
                break
            else:
                lmbd = lmbd * 2.0
        
        print "Interval found: [%02.03f, %02.03f]" % (lmbdL, lmbdR)
        return self.getLambdaExhaustionCDFinterval(timeout, poolsize, prob, lmbdL, lmbdR)
        
    def getLambdaExhaustionCDFinterval(self, timeout, poolsize, prob, l, r):
        eps = 0.0001
        while l < r:
            c = (l+r)/2
            probc = 1.0 - poisson.cdf(poolsize, c * timeout)
            
            print "\nNew iteration [%02.03f, %02.03f]; c=%02.3f probc=%02.3f vs. prob=%02.3f" % (l, r, c, probc, prob)
            
            if probc >= (prob-eps) and probc <= (prob+eps): break 
            if probc < prob: l = c
            if probc > prob: r = c
        return l       
    
    def getLambdaExhaustion(self, natA):
        '''
        Get a lambda that will cause exhaustion for a given NAT
        '''
        # at first we have to find proper interval where to find by binary search
        timeout  = natA.timeout
        poolsize = natA.poolLen
        
        lmbd = 1.0
        lmbdL=-1
        lmbdR=-1
        
        while (lmbdL==-1) or (lmbdR==-1):
            print "current lambda: %02.3f" % lmbd
            t = self.poolExhaustion(timeout, poolsize, lmbd)
            
            # left side of the interval. If fits, set and continue to find right side, otherwise keep left 
            if lmbdL==-1:
                if t > timeout: 
                    lmbdL = lmbd
                    lmbd  = lmbd * 2.0
                    continue
                else: 
                    lmbd = lmbd/2.0
                    continue
            
            # right side of the interval, if here, we have left side and finding the right side
            if t < timeout:
                lmbdR = lmbd
                break
            else:
                lmbd = lmbd * 2.0
        
        print "Interval found: [%02.03f, %02.03f]" % (lmbdL, lmbdR)
        return self.getLambdaExhaustionInterval(timeout, poolsize, lmbdL, lmbdR)
        
    def getLambdaExhaustionInterval(self, timeout, poolsize, lmbdL, lmbdR):
        '''
        Recursive binary search for finding lambda that will exhaust port pool size
        ''' 
        eps = 100 # search epsilon, accuracy, 100ms
        t   = 0
        while lmbdL < lmbdR and (lmbdR - lmbdL) > 0.00001:
            
            print "\nNew iteration [%02.03f, %02.03f]" % (lmbdL, lmbdR)
            cLmbd = (lmbdL+lmbdR) / 2.0
            t = self.poolExhaustion(timeout, poolsize, cLmbd)
            
            if t >= (timeout-eps) and t <= (timeout+eps): break 
            if t < timeout: lmbdR = cLmbd
            if t > timeout: lmbdL = cLmbd
        return lmbdL
    
    def portDistributionFunction(self, lmbd, t, pstep):
        '''
        Measures distribution function of the ports on NAT with Poisson process.
        This matters since the whole nature of NAT is incremental. Number of 
        new connections in different time windows should hold distribution also
        considered together, Po(lmbd*(t1+t2)), but considering port numbers it 
        makes a difference, also taking my port allocation into account. 
        
        For instance port 6 can be reached by 2,2,2 or 3,3. 
        '''
        ports = 1000
        portDistrib = []
        for i in range(0, ports): portDistrib.append(0)
        for sn in range(0, 10000): # self.simulationRounds):
            
            ports = []
            step    = 0
            curPort = 0
            while curPort < ports:
                if(curPort!=0 and pstep!=-1 and step==pstep):
                    portDistrib[curPort] += 1        # port is used
                
                curPort += self.poisson(lmbd, t) # add new connections by Poisson process
                curPort += 1                     # add my port, I made it by a new connection
                step    += 1
                
                ports.append(curPort)
        pass
    
        print "Here comes a histogram..."
        #hist = np.histogram(portHist, bins=range(0, ports))
        
        pos = np.arange(ports)
        width = 1.0     # gives histogram aspect to the bar diagram
        
        ax = plt.axes()
        ax.set_xticks(pos + (width / 2))
        ax.set_xticklabels(range(0, ports))
        
        plt.bar(pos, portDistrib, width, color='r')
        plt.show()
    
        
    
# main executable code    



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='NAT simulator.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o','--output',help='Output file name from finder', required=False, default='graph.txt')
    parser.add_argument('-t','--space',help='Time in ms to wait between packet send', required=False, default=10, type=int)
    parser.add_argument('-l','--lmbd_start',help='On which lambda to start', required=False, default=-1, type=float)
    parser.add_argument('-s','--strategy',help='Strategy to use (poisson, i2j, fibo, their)', required=False, default='poisson')
    args = parser.parse_args()
    
    ns = NatSimulation()
    
    # create a symmetric nat both for Alice and Bob
    natA = SymmetricIncrementalNat()
    natB = SymmetricIncrementalNat()
    
    natA.init(None)
    natB.init(None)
    
 #   print ns.poolExhaustionNat(natA, 3*60*1000)
#    
#    print "Lambda that will exhaust given NAT: "
#    lmbd = ns.getLambdaExhaustion(natA)
#    print "\nLambda that will exhaust given NAT: ", lmbd
#    
#    print ns.getLambdaExhaustionCDF(natA, 0.999)
#    sys.exit()
    
    if args.strategy == 'i2j':
        print "I2J Strategy: "
        strategy = I2JStragegy()
    elif args.strategy == 'fibo':
        print "Fibonacci strategy"
        strategy = FiboStrategy()
    elif args.strategy == 'their':
        print "Their strategy"
        strategy = TheirStragegy()
    elif args.strategy == 'poisson':
        print "Poisson strategy"    
        strategy  = PoissonStrategy()
    strategy.init(None)
    
    # Their
    #strategy.delta = [200, 200]

    ns.portDistributionFunction(0.2, 10, 60)
    sys.exit(3)
    
    #strategy.dupl = True
    #strategy.coef = 1.8
    ns.simulation(natA, natB, strategy)
    sys.exit(3)
    
    # generating graph for moving lambda
    f = open(args.output, 'a+')
    ns.portScanInterval = args.space
    
    f.write("New start at %s; scanInterval=%d; strategy=%s\n" % (time.time(), ns.portScanInterval, args.strategy))
    print "Scanning port interval: %d" % ns.portScanInterval
    
    lmbdArr = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, \
               0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, \
               0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, \
               0.2, 0.21, 0.22, 0.23, 0.24, 0.25]
    for clmb in lmbdArr:
        print "# Current lambda: %03.4f" % clmb
        if args.lmbd_start!=-1 and clmb < args.lmbd_start: continue
        
        ns.lmbd = clmb
        try:
            if args.strategy == 'poisson':
                res = ns.coefFinder(natA, natB, strategy, 0.10, 0.1)
                f.write("%03.4f|%03.4f|%03.4f|%03.4f\n" % (ns.lmbd, ns.portScanInterval, res[0], res[1])) # python will convert \n to os.linesep
            else:
                res = ns.simulation(natA, natB, strategy)
                f.write("%03.4f|%03.4f|%03.4f|%03.4f\n" % (ns.lmbd, ns.portScanInterval, res[0], res[2])) # python will convert \n to os.linesep
            f.flush()
        except Exception, e:
            print "Exception!", e
    f.close()
    #ns.simulateThem()
    

        