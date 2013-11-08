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

# Fibonacchi list generator
def fibGenerator():
    a, b = 0, 1
    yield 0
    while True:
        a, b = b, a + b
        yield a
            
class Strategy:
    '''
    Base abstract class for NAT strategy
    '''
    def init(self, params=None):
        raise Exception("Not implemented yet...")
    def reset(self, params=None):
        self.init(None)
    def next(self, party, step):
        raise Exception("Not implemented yet...")   # return (srcPort, dstPort)

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
    delta = []
    def init(self, params=None):
        pass
    def next(self, party, step):
        if party==0: return (step,self.delta[0])
        if party==1: return (step,self.delta[1])
    
class I2JStragegy(Strategy):
    def init(self, params=None):
        pass
    def next(self, party, step):
        if party==0: return (0,step)
        if party==1: return (0,2*step)
        
class FiboStrategy(Strategy):
    fibn = []
    b    = []
    def init(self, params=None):
        self.fibn = []
        fib = fibGenerator()
        for n in range(22):
            self.fibn.append(next(fib))
        print self.fibn
        
        for i in range(1, len(self.fibn)-1):
            for j in range(0, self.fibn[i-1]):
                self.b.append(self.fibn[i+1] + j)
        print self.b
    def reset(self):
        pass
    def next(self, party, step):
        return (0, self.b[step])

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
    lmbd = 0.367
    
    # Number of miliseconds for silent period to take [ms].
    # Based on basic ping / round trip time it takes to communicate 
    # IP with another peer 
    silentPeriodBase = 500
    
    # Lambda for Pois(lmbd) for silent period variability.
    # Silent period time = silentPeriodBase + Pois(lmbd) [ms]
    silentPeriodlmbd = 100
    
    # Number of rounds for simulation
    simulationRounds = 1
    
    # Number of errors that are handled by algorithm 
    errors = 1000
    
    # number of milliseconds between consecutive port scans
    portScanInterval = 10
    
    # number of connections to establish
    numCon = 1
    
    def poisson(self, lmbd, t):
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
        return self.poisson(self.lmbd, tim)
    
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
    
    def simulation(self, natA, natB, strategy):
        '''
        Simple simulation of NAT traversal algorithm.
        '''
        
        nats = [natA, natB]
        successCnt = 0.0
        stopOnFirstMatch = self.simulationRounds != 1
        for sn in range(0, self.simulationRounds):
            # reset NATs
            nats[0].reset()
            nats[1].reset()
            strategy.reset()
            
            # generate silent period time
            curSilentA = self.silentPeriodBase + self.poisson(self.silentPeriodlmbd, 1)
            curSilentB = self.silentPeriodBase + self.poisson(self.silentPeriodlmbd, 1)
            print "\n##%03d. Current silent period time: [%03.3f, %03.3f]" % (sn, curSilentA, curSilentB) 
            
            # generate new TCP connections for silent period on both sides, same lambda.
            kA = self.poisson(self.lmbd, curSilentA)
            kB = self.poisson(self.lmbd, curSilentB)
            print "New connections meanwhile silent period [%02d, %02d]" % (kA, kB)
            
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
            
            print "totalLags [%02d %02d]" % (totalLagA[0], totalLagA[1])
            
            # OK is there any intersection in both sets?
            res = list(scanA[0].intersection(scanA[1]))
            # sort by minimum element in tuple
            res.sort(key=lambda tup: min(tup[0], tup[1]))
            
            # Generate DOT graph
            if self.simulationRounds==1:
                self.generateDot(portsA[0], portsA[1], scanA[0], scanA[1], res)
            
            if (len(res) == 0): 
                print "Algorithm failed, no intersecting points"
                continue
            
            print "RES: ", res, "i=%02d" % mapA[0][res[0][0]], "; j=%02d" % mapA[1][res[0][1]]
            successCnt += 1.0
            
            
        print "Success count: %02.3f ; cnt=%d" % (successCnt / self.simulationRounds, successCnt)
        
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
        while True: #lmbdL < lmbdR and (lmbdR - lmbdL) > 0.0001:
            
            print "\nNew iteration [%02.03f, %02.03f]" % (lmbdL, lmbdR)
            cLmbd = (lmbdL+lmbdR) / 2.0
            t = self.poolExhaustion(timeout, poolsize, cLmbd)
            
            if t >= (timeout-eps) and t <= (timeout+eps): break 
            if t < timeout: lmbdR = cLmbd
            if t > timeout: lmbdL = cLmbd
        return lmbdL
    
# main executable code    
if __name__ == "__main__":    
    ns = NatSimulation()
    
    # create a symmetric nat both for Alice and Bob
    natA = SymmetricIncrementalNat()
    natB = SymmetricIncrementalNat()
    
    
    natA.init(None)
    natB.init(None)
    
    print ns.poolExhaustionNat(natA, 3*60*1000)
    
    print "Lambda that will exhaust given NAT: "
    lmbd = ns.getLambdaExhaustion(natA)
    print "\nLambda that will exhaust given NAT: ", lmbd
    
    print ns.getLambdaExhaustionCDF(natA, 0.999)
    sys.exit()
    
    print "I2J Strategy: "
    strategy = FiboStrategy()
    #strategy = I2JStragegy()
    #strategy = TheirStragegy()
    strategy.init(None)
    
    # Their
    #strategy.delta = [100, 100]
    
    ns.simulation(natA, natB, strategy)
    #ns.simulateThem()
    

        