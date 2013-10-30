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
from portNums import lastPort

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
    def next(self, party, step):
        raise Exception("Not implemented yet...")   # return (srcPort, dstPort)

class Nat:
    '''
    Base abstract class for NAT allocation
    '''
    def init(self, params=None):
        raise Exception("Not implemented yet...")    
    def get(self, srcIP, srcPort, dstIP, dstPort):
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
    
class SymmetricIncrementalNat(Nat):
    # timeout of a created connection in seconds
    timeout = 3*60
    # allocation table for NAT; key = quartet; value = (external port, time added)
    allocations = None
    # port pool available for new allocations
    pool = None
    poolLen = 0
    # reversed mapping from allocations (duplicit structure, could be infered from allocations)
    allocatedPorts = None
    # index of last allocated port. Index to pool[]
    lastPort = 0
    
    def init(self, params=None):
        self.allocations = {}
        self.pool = range(0, 65536)
        self.poolLen = len(self.pool)
        self.allocatedPorts = {}
        
    def alloc(self, srcIP, srcPort, dstIP, dstPort, timeNow):
        q = Quartet(srcIP, srcPort, dstIP, dstPort)
        
        # Check for existing allocation for a given quartet
        if q in self.allocations:
            tup = self.allocations.get(q)
            # check expiration time, if a record is too old, it has to be removed from alloc table.
            if (tup[1] + self.timeout) < timeNow:
                del self.allocatedPorts[tup[0]]  # delete from allocation set
                del self.allocations[q]             # expired -> delete from allocation table
            else:
                self.allocations[q][1] = timeNow    # update last query access
                return tup[0]                       # external port returned
        
        # if here -> not in allocation list, create a new allocation
        tries=0
        while tries <= self.poolLen:
            lastPort = (lastPort + 1) % self.poolLen # linear port allocation rule here
            tries    += 1                            # check pool exhaustion
            if self.pool[lastPort] in self.allocatedPorts:
                # next port is already allocated, find quartet for it
                q2  = self.allocatedPorts[self.pool[lastPort]]
                tup = self.allocations.get(q2)
                # check expiration first
                if (tup[1] + self.timeout) < timeNow:
                    del self.allocatedPorts[tup[0]]  # delete from allocation set
                    del self.allocations[q2]          # expired -> delete from allocation table
                    break                            # slot is free now, can allocate
                else: continue                       # slot is in use, continue with search
            else: break                              # slot is free, assign
        # check if pool is exhausted - all ports are allocated currently
        if tries == self.poolLen:
            raise Exception("Port pool exhausted")
        # create a new allocation
        port = self.pool[lastPort]
        self.allocations[q] = (port, timeNow)
        self.allocatedPorts[port] = q
        
        return port

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
    lmbd = 0.05
    
    # Number of miliseconds for silent period to take [ms].
    # Based on basic ping / round trip time it takes to communicate 
    # IP with another peer 
    silentPeriodBase = 500
    
    # Lambda for Pois(lmbd) for silent period variability.
    # Silent period time = silentPeriodBase + Pois(lmbd) [ms]
    silentPeriodlmbd = 100
    
    # Number of rounds for simulation
    simulationRounds = 1000
    
    # Number of errors that are handled by algorithm 
    errors = 3000
    
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
        
    def simulateThem(self):
        '''
        Simulates delta guessing algorithm
        '''
        successCnt = 0.0
        for sn in range(0, self.simulationRounds):
            # generate silent period time
            curSilentA = self.silentPeriodBase + self.poisson(self.silentPeriodlmbd, 1)
            curSilentB = self.silentPeriodBase + self.poisson(self.silentPeriodlmbd, 1)
            print "\n##%03d. Current silent period time: [%03.3f, %03.3f]" % (sn, curSilentA, curSilentB) 
            
            # generate new TCP connections for silent period on both sides, same lambda.
            kA = self.poisson(self.lmbd, curSilentA)
            kB = self.poisson(self.lmbd, curSilentB)
            print "New connections meanwhile silent period [%02d, %02d]" % (kA, kB)
            
            # set delta ports
            deltaA = 30
            deltaB = 30
            
            # now simulate the protocol, phase with port scanning
            targetA = 30
            targetB = 30
            
            curA = kA
            curB = kB
            
            mapA  = {} # mapping of current port to i
            mapB  = {} # mapping of current port to j
            scanA = set([])
            scanB = set([])
            portsA = set([])
            portsB = set([])
            totalLagA, totalLagB = 0, 0
            
            dot = "digraph finite_state_machine {\n"
            for i in range(0, max(targetA, targetB)):
                # A scan
                dstA  = deltaA
                
                # waiting between consecutive scans
                curLag = self.poisson(self.lmbd, self.portScanInterval)
                scanA.add((curA, dstA))
                portsA.add(curA)
                mapA[curA] = i
                #print "A scan: [%03d] --> [%03d] lag=%02d i=%03d" % (curA, dstA, curLag, i)
                
                curA += curLag + 1
                totalLagA += curLag
                
                # B scan
                dstB  = deltaB
                
                # waiting between consecutive scans
                curLag = self.poisson(self.lmbd, self.portScanInterval)
                toAdd  = (dstB, curB)
                scanB.add(toAdd) # reversed order here is OK, want to find intersection of lists
                portsB.add(curB)
                mapB[curB] = i
                #print "B scan: [%03d] --> [%03d] lag=%02d i=%03d" % (curB, dstB, curLag, i)
                                
                curB += curLag + 1
                totalLagB += curLag
            
            print "totalLags [%02d %02d]" % (totalLagA, totalLagB)
            
            # OK is there any intersection in both sets?
            res = list(scanA.intersection(scanB))
            # sort by minimum element in tuple
            res.sort(key=lambda tup: min(tup[0], tup[1]))
            # Generate DOT graph
            if self.simulationRounds==1:
                self.generateDot(portsA, portsB, scanA, scanB, res)
            
            if (len(res) == 0): 
                print "Algorithm failed, no intersecting points"
                continue
            
            print "RES: ", res, "i=%02d" % mapA[res[0][0]], "; j=%02d" % mapB[res[0][1]]
            successCnt += 1.0
            
            
        print "Success count: %02.3f ; cnt=%d" % (successCnt / self.simulationRounds, successCnt)
        
    def simulation(self):
        '''
        Simple simulation of NAT traversal algorithm.
        '''
        
        # fib
        fibn = []
        fib = self.fibGenerator()
        for n in range(22):
            fibn.append(next(fib))
            
        print fibn
        b = []
        for i in range(1, len(fibn)-1):
            for j in range(0, fibn[i-1]):
                b.append(fibn[i+1] + j)
                
        print b
        #return 1
        
        successCnt = 0.0
        stopOnFirstMatch = self.simulationRounds != 1
        for sn in range(0, self.simulationRounds):
            # generate silent period time
            curSilentA = self.silentPeriodBase + self.poisson(self.silentPeriodlmbd, 1)
            curSilentB = self.silentPeriodBase + self.poisson(self.silentPeriodlmbd, 1)
            print "\n##%03d. Current silent period time: [%03.3f, %03.3f]" % (sn, curSilentA, curSilentB) 
            
            # generate new TCP connections for silent period on both sides, same lambda.
            kA = self.poisson(self.lmbd, curSilentA)
            kB = self.poisson(self.lmbd, curSilentB)
            print "New connections meanwhile silent period [%02d, %02d]" % (kA, kB)
            
            # assume we are always starting from port 0 on both sides
            i = kA + 2*kB
            j = kA +   kB
            print "Offsets that would match without any further errors i=%02d; j=%02d" % (i, j)
            
            # now simulate the protocol, phase with port scanning
            targetA = 2 * self.numCon * self.errors
            targetB =     self.numCon * self.errors
            
            curA = kA
            curB = kB
            
            mapA  = {} # mapping of current port to i
            mapB  = {} # mapping of current port to j
            scanA = set([])
            scanB = set([])
            portsA = set([])
            portsB = set([])
            totalLagA, totalLagB = 0, 0
            
            stageNumA = 0
            stageNumB = 0
            stageChange = 200 #self.errors / 10
            
            for i in range(0, max(targetA, targetB)):
                # A scan
                dstA  = b[i]#1*i #- stageChange*(stageNumA)/10.0# destination of scan o the other side
                #dstA = i
                
                # waiting between consecutive scans
                curLag = self.getNumOfNewConnections(self.portScanInterval)
                toAdd  = (curA, dstA)
                scanA.add(toAdd)
                portsA.add(curA)
                mapA[curA] = i
                #print "A scan: [%03d] --> [%03d] lag=%02d i=%03d" % (curA, dstA, curLag, i)
                
                if stopOnFirstMatch and toAdd in scanB:    
                    break
                
                curA += curLag + 1
                totalLagA += curLag
                
                # B scan
                if i <= targetB:
                    #dstB  = 2*(i/self.numCon) + i%self.numCon + 20*(stageNumB)
                    #dstB  = (2+stageNumB)*i
                    dstB  = b[i] #2*i
                    #dstB  = 2*i
                    
                    # waiting between consecutive scans
                    curLag = self.getNumOfNewConnections(self.portScanInterval)
                    toAdd  = (dstB, curB)
                    scanB.add(toAdd) # reversed order here is OK, want to find intersection of lists
                    portsB.add(curB)
                    mapB[curB] = i
                    #print "B scan: [%03d] --> [%03d] lag=%02d i=%03d" % (curB, dstB, curLag, i)
                    
                    if stopOnFirstMatch and toAdd in scanA:
                        break
                    
                    curB += curLag + 1
                    totalLagB += curLag
                
                #if (i/2.0) >= stageChange * (stageNumA+1):
                #    stageNumA+=1
                #    print "StageA++  i=%03d" % i
                if i       >= stageChange * (stageNumB+1):
                    stageNumB+=1
                    print "StageB++  i=%03d" % i
            
            print "totalLags [%02d %02d]" % (totalLagA, totalLagB)
            
            # OK is there any intersection in both sets?
            res = list(scanA.intersection(scanB))
            # sort by minimum element in tuple
            res.sort(key=lambda tup: min(tup[0], tup[1]))
            
            # Generate DOT graph
            if self.simulationRounds==1:
                self.generateDot(portsA, portsB, scanA, scanB, res)
            
            if (len(res) == 0): 
                print "Algorithm failed, no intersecting points"
                continue
            
            print "RES: ", res, "i=%02d" % mapA[res[0][0]], "; j=%02d" % mapB[res[0][1]]
            successCnt += 1.0
            
            
        print "Success count: %02.3f ; cnt=%d" % (successCnt / self.simulationRounds, successCnt)
        
    def generateDot(self, portsA, portsB, scanA, scanB, res):
        '''
        Generate DOT image
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
                acounter = [i for i in scanB if i[1] == bcounter][0][0] if bcounter in portsB else -1
                dot = dot + desc % (8, 1.5*p, "%03d -> %03d -> %03d\\n%s" % (p, bcounter, acounter, str(whoping)), "A"+ps)
            else:
                dot = dot + desc % (8, 1.5*p, "%s" % (str(whoping)), "A"+ps)
                
            
            # B
            whoping  = [i for i in scanB if i[0] == p]
            
            color = 7
            if p in portsB:
                color = 3 if len(whoping) > 0 else 1
            if inB:
                color = "#00ff005f"
            dot = dot + line % (color, 130, 1.5*p, p, "B"+ps)
            
            # B desc
            if p in portsB:
                acounter = [i for i in scanB if i[1] == p][0][0]
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
    
# main executable code    
if __name__ == "__main__":
    ns = NatSimulation()
    

    
    ns.simulation()
    #ns.simulateThem()
    

        