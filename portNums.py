import os, sys, fileinput, re

if sys.argv[1] == None or os.path.exists(sys.argv[1])==False:
    print "Usage: %s file_to_analyze" % sys.argv[0]
    sys.exit()


class portRec:
    srcIP=0
    srcPort=0
    dstPort=0
    timeRecv=0
    order=0
    
    timeSent=0
    intPort=0
    
    def __init__(self):
        pass
    
    def __str__(self):
        return "[int]:%05d -> %s:%05d -> [srv]:%05d  order:%05d  sent: %s  recv: %s" % (int(self.intPort), self.srcIP, int(self.srcPort), int(self.dstPort), int(self.order), self.timeSent, self.timeRecv)
    
    def csv(self):
        return "%05d;%s;%05d;%05d;%05d;%s;%s" % (int(self.intPort), self.srcIP, int(self.srcPort), int(self.dstPort), int(self.order), self.timeSent, self.timeRecv)
    
# compiled regular expressions
flineRe = re.compile(r"^[\d]+:[\d]+:[\d]+\.[\d]+ IP \(tos 0x0, ttl")  # regex for first line record
addrsRe = re.compile(r"^([\d]+\.[\d]+\.[\d]+\.[\d]+)\.([\d]+) > ([\d]+\.[\d]+\.[\d]+\.[\d]+)\.([\d]+)")  # regex for addresses
recRe = re.compile(r"\|\|t=([\d]+);s=([\d]+);d=([\d]+)\|\|")     # regex for record

def parseLastRecord(lastRec, order):
    toRet = portRec()
    lstIdx = len(lastRec)-1
    
    toRet.order = order
    line1split = lastRec[0].strip().split(" ", 2)
    toRet.timeRecv = line1split[0]
    
    line2strip = lastRec[1].strip()
    m = addrsRe.match(line2strip)
    if(m==None):
        print "Warning, unrecognized record: ", lastRec
        return None
    toRet.srcIP = m.group(1)
    toRet.srcPort = int(m.group(2))
    toRet.dstPort = int(m.group(4))
    
    m = recRe.search(lastRec[lstIdx])
    if(m==None):
        #print "Record not detected in the last line: ", lastRec[lstIdx]
        return toRet
    toRet.timeSent = int(m.group(1))
    toRet.intPort = int(m.group(2))
    
    return toRet
    

# filename to analyze    
fname=sys.argv[1]

# DB
portdb = []

portCounts = []
for i in range(0,65536): portCounts.append(0)

dstPorts = set()

# stores last parsed record
curOrder = 0
lastRecord=[]
for line in fileinput.input(fname):
    if flineRe.match(line):             # new line - finish last record processing, start new one
        #print "Last record: ", lastRecord
        
        # PROCESS LAST RECORD HERE
        if len(lastRecord)>0:
            objRec = parseLastRecord(lastRecord, curOrder)
            curOrder+=1
        
            #print objRec
            portdb.append(objRec)
            
            portCounts[objRec.srcPort]+=1
            dstPorts.add(objRec.dstPort)
            
            if objRec.srcPort==10000 or objRec.srcPort==10001:
                print "TargetPort: %0d5, data: %s" % (objRec.dstPort, str(objRec))  
                
        
        # new last record
        lastRecord = [line]
    else:
        lastRecord.append(line)

print "DONE reading data"
print "Distinct destination ports: ", dstPorts

firstPort=-1
lastPort=-1
for i in range(0,65536): 
    if portCounts[i] > 0 and firstPort==-1:
        firstPort=i
    if portCounts[i] > 0:
        lastPort=i

print "Port min:max interval: [%d, %d]" % (firstPort, lastPort)

# look for gaps in interval
for i in range(firstPort,lastPort+1): 
    if portCounts[i]==0:
        print "PortGap: %05d" % i
        
#for i in portdb:
    #print i.csv()
 

