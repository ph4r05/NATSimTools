
import os, sys, socket, threading, SocketServer, time, traceback, copy
from threading import Thread, Lock
from datetime import datetime

def utc():
    return int(datetime.utcnow().strftime("%s"))

# transaction record
class txobj:
    """Transaction object in transaction database"""
    txname=''          # Identifier of transaction
    participants=0     # number of participants in transaction. 0=empty transaction; 1=half/open; 2=established
    startTime=0        # utc when participants -> 1
    fullTime=0         # utc when participants -> 2
    plist=[]           # participant list. Example: [['127.0.0.1', 88, id, type], ['192.168.1.1', 45, id, type]]
    pSockets=[]        # [(socket, client_addr)]
    params=[]
    
    def __init__(self):
        self.plist=[]
        self.pSockets=[]
        self.params=[]

# not used anymore, just for demonstration purposes
class ThreadingUDPServer(SocketServer.ThreadingMixIn, SocketServer.UDPServer): 
    """
    This class works similar to the TCP handler class, except that
    self.request consists of a pair of data and client socket, and since
    there is no connection the client address must be given explicitly
    when sending data back via sendto().
    """
    def handle(self):
        data = self.request[0].strip()
        socket = self.request[1]
        print "{} wrote:".format(self.client_address[0])
        print data
        socket.sendto(data.upper(), self.client_address)

# just simple request handler for threaded UDP server
class ThreadedUDPRequestHandler(SocketServer.BaseRequestHandler):
    def handle(self):
        #data = self.request.recv(1024)          # TCP variant - receive 1024 bytes from byte stream
        data = self.request[0].strip()           # Receive UDP message, stripped from white space on both ends
        socket = self.request[1]                 # Socket to send message back
        cur_thread = threading.current_thread()
        response = "%s: %s Addr: %s Time: %s" % (cur_thread.name, data, str(self.client_address), utc())
        print response
        print data
        
        # do some serious stuff
        self.server.txman.txdbLock.acquire()      # acquire mutex for transactions
        try:
            action, tail = data.split("|", 1)
            print "Action: [", action, "]; tail: [", tail
            
            payload = tail.split("|")
            if action=='txbegin' and len(payload)>=3:
                txid=payload[0].strip()
                myid=payload[1].strip()
                mytype=payload[2].strip()
                tx=None
                
                print "Transaction start; txid=%s; pId=%s" % (txid, myid)
                print self.server.txman.txdb
                if (txid in self.server.txman.txdb)==False:
                    tx = txobj()
                    tx.txname = txid
                    tx.startTime = utc()
                    tx.participants=1
                    tx.pSockets=[(socket, self.client_address)]
                    tx.plist=[(self.client_address[0], self.client_address[1], myid, mytype)]
                    tx.params=[]
                else:
                    tx = self.server.txman.txdb[txid]
                    notInTransaction=True
                    for i,p in enumerate(copy.deepcopy(tx.plist)):  # check if is already present in transaction
                        if myid==p[2]:  # if present, just update records for him 
                            print "Already in transaction: [%s]" % myid, "list: ", tx.plist
                            notInTransaction=False
                            tx.pSockets[i]= (socket, self.client_address)
                            tx.plist[i] = (self.client_address[0], self.client_address[1], myid, mytype)
                            break
                    if notInTransaction:    
                        tx.participants = 2
                        tx.fullTime = utc()
                        tx.pSockets.append((socket, self.client_address))
                        tx.plist.append((self.client_address[0], self.client_address[1], myid, mytype))
                
                
                if len(payload)>3:    # analyze some parameters for transaction (delay for example)
                    for param in payload[3:]:
                        print "Analyzing parameter: ", param
                        try:
                            pname,pval = param.strip().split("=")
                            tx.params.append((pname,pval))
                        except Exception,e:
                            print "Parameter error [%s]" % param, e
                
                self.server.txman.txdb[txid] = tx
                if tx.participants>=2: # if transaction is saturated, do some stuff
                    peers = ["%s;%s;%s;%s" % (p[0],p[1],p[2],p[3]) for p in tx.plist]
                    params = ["%s=%s" % (p[0], str(p[1])) for p in tx.params]
                    
                    
                    fullDelay=0
                    for pname,pval in tx.params:
                        if "fullDelay" == pname:
                            fullDelay=int(pval)
                    
                    if fullDelay!=0:
                        print "Some full delay here, going to sleep %s ms" % fullDelay
                        time.sleep(fullDelay/1000.0)
                    
                    txDescription="txstarted|"+txid+"|"+str(tx.startTime)+"|"+str(utc())+"|PEERS|"+("|".join(peers))+"|PARAMETERS|"+("|".join(params))
                    print "Going to broadcast transaction to peers [[ %s ]]" % txDescription
                    for sock,addr in tx.pSockets:
                        sock.sendto(txDescription, addr)
                        print "Sending something to: ", addr
                    print "Transaction broadcasted"
                    
                    self.server.txman.txdb[txid] = None
                    del self.server.txman.txdb[txid]
                    del tx
                    print "Transaction deleted"
                else:
                    print "Waiting to saturate transaction [%s]" % txid                
            elif action=='txsetup' and len(payload)>=2:
                print "Transaction setup"
                txid=payload[0].strip()
                for params in payload[1:]:
                    print "Param: %s" % params
                socket.sendto(data.upper(), self.client_address)
            else:
                print "Unknown command"
                socket.sendto(data.upper(), self.client_address)
        except Exception,e:
            print "Some exception happened:", e
            traceback.print_exc()
        finally:
            self.server.txman.txdbLock.release()
        #self.request.sendall(response)        # TCP variant

# transaction manager shared among threaded servers
class TXManager():
    txdb = {}           # Transaction database
    txdbLock = None     # mutex for transaction database
    
    def __init__(self):
        self.txdbLock = Lock()
    
    def cleanOldTxs(self):
        """Cleans old transactions"""
        self.txdbLock.acquire()
        try:
            curTime = utc()
            for i,p in enumerate(copy.deepcopy(self.txdb)):
                tx = self.txdb[p]
                if (tx.participants==1 and (curTime-tx.startTime)>60*10) \
                or (tx.participants==2 and (curTime-tx.fullTime)>60*10):
                    print "Removing old transaction: %d" % i
                    del self.txdb[p]
        except Exception,e:
            print "Some exception during cleaning old transactions: ", e
            traceback.print_exc()
        finally:
            self.txdbLock.release()
    pass

# threaded UDP server with central transaction database
class ThreadedUDPServer(SocketServer.ThreadingMixIn, SocketServer.UDPServer):
    txman = None
    def __init__(self, server_address, RequestHandlerClass, bind_and_activate=True):
        SocketServer.UDPServer.__init__(self, server_address, RequestHandlerClass, bind_and_activate=bind_and_activate)
    def printMe(self):
        print "Hello from server!"

# simple client routine to connect to server and communicate something
def client(ip, port, message):
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM) #socket.SOCK_STREAM)
    sock.connect((ip, port))
    try:
        # As you can see, there is no connect() call; UDP has no connections.
        # Instead, data is directly sent to the recipient via sendto().
        sock.sendto(message + "\n", (HOST, PORT))
        received = sock.recv(1024)
        #sock.sendall(message)
        #response = sock.recv(1024)
        print "Received: {}".format(received)
    finally:
        sock.close()

if __name__ == "__main__":
    # Port 0 means to select an arbitrary unused port
    HOST, PORT, numServers = "0.0.0.0", 9999, 100
    
    # initialize common transaction manager for servers
    txman = TXManager()
    servers = []
    for i in range(0, 100):
        server = ThreadedUDPServer((HOST, PORT+i), ThreadedUDPRequestHandler)
        server.txman = txman
        ip, port = server.server_address
        # Start a thread with the server -- that thread will then start one
        # more thread for each request
        server_thread = threading.Thread(target=server.serve_forever)
        # Exit the server thread when the main thread terminates
        server_thread.daemon = True
        server_thread.start()
        print "Server loop running in thread:", server_thread.name, "IP:port %s:%s" % (ip, port)
    #client(ip, port, "Hello World 1")

    #server.shutdown()
    #server.serve_forever()
    i=0;
    try:
        while True:
            time.sleep(1)
            # clean old transactions
            if (i % 5) == 0:
                txman.cleanOldTxs()
            i = (i + 1) % 65535
    except Exception,e:
        print "Exception during cleaning old transactions; e: ", e
    print "Finishing process..."
    
    #server = SocketServer.UDPServer((HOST, PORT), ThreadingUDPServer)
    #server.serve_forever()
