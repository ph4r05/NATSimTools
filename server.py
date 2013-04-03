import os, sys, socket, threading, SocketServer, time, traceback, copy
from threading import Thread, Lock
from datetime import datetime

sDataArr=[]
def utc():
    return int(datetime.utcnow().strftime("%s"))

class txobj:
    """Transaction object in transaction database"""
    txname=''          # Identifier of transaction
    participants=0     # number of participants in transaction. 0=empty transaction; 1=half/open; 2=established
    startTime=0        # utc when participants -> 1
    fullTime=0         # utc when participants -> 2
    plist=[]           # participant list. Example: [['127.0.0.1', 88, type], ['192.168.1.1', 45, type]]
    p1IP="0.0.0.0"     # participant #1 IP address
    p1Port=0           # participant #1 port number
    p2IP="0.0.0.0"     # participant #2 IP address

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
        response = "{}: {} Addr: {} Time: {}".format(cur_thread.name, data, str(self.client_address), utc())
        print response
        print data
        
        # do some serious stuff
        self.server.txdbLock.acquire()      # acquire mutex for transactions
        try:
            action, tail = data.split("|", 1)
            print "Action: [", action, "]; tail: [", tail
            
            payload = tail.split("|")
            if action=='txbegin' and len(payload)>=3:
                txid=payload[0].strip()
                myid=payload[1].strip()
                mytype=payload[2].strip()
                print "Transaction start; txid=%s; pId=%s" % (txid, myid)
                if (txid in self.server.txdb)==False:
                    tx = txobj()
                    tx.startTime = utc()
                    tx.participants=1
                    tx.plist.append((self.client_address[0], self.client_address[1],))
                    
                    
                else:
                    print "Existing transaction"
                
                
            elif action=='txsetup':
                print "Transaction setup"
                
                
            else:
                print "Unknown command"
            socket.sendto(data.upper(), self.client_address)
        except Exception,e:
            print "Some exception happened:", e
            traceback.print_exc()
        finally:
            self.server.txdbLock.release()
                            
        #self.request.sendall(response)        # TCP variant

# threaded UDP server with central database
class ThreadedUDPServer(SocketServer.ThreadingMixIn, SocketServer.UDPServer):
    txdb = {}           # Transaction database
    txdbLock = None     # mutex for transaction database
    def __init__(self, server_address, RequestHandlerClass, bind_and_activate=True):
        SocketServer.UDPServer.__init__(self, server_address, RequestHandlerClass, bind_and_activate=bind_and_activate)
        self.txdbLock = Lock()
        
    def printMe(self):
        print "Hello from server!"
    
    def cleanOldTxs(self):
        """Cleans old transactions"""
        self.txdbLock.acquire()
        try:
            curTime = utc()
            for i,p in enumerate(copy.deepcopy(self.txdb)):
                if (p.participants==1 and (curTime-p.startTime)>60*5) \
                or (p.participants==2 and (curTime-p.fullTime)>60*5):
                    print "Removing old transaction: %d" % i
                    del self.txdb[i]
        except Exception,e:
            print "Some exception during cleaning old transactions: ", e
        finally:
            self.txdbLock.release()
    pass

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
    HOST, PORT = "localhost", 9999

    server = ThreadedUDPServer((HOST, PORT), ThreadedUDPRequestHandler)
    ip, port = server.server_address

    # Start a thread with the server -- that thread will then start one
    # more thread for each request
    server_thread = threading.Thread(target=server.serve_forever)
    # Exit the server thread when the main thread terminates
    server_thread.daemon = True
    server_thread.start()
    print "Server loop running in thread:", server_thread.name

    #client(ip, port, "Hello World 1")
    #client(ip, port, "Hello World 2")
    #client(ip, port, "Hello World 3")

    #server.shutdown()
    #server.serve_forever()
    i=0;
    try:
        while True:
            time.sleep(1)
            # clean old transactions
            if (i % 5) == 0:
                server.cleanOldTxs()
            i = (i + 1) % 65535
    except Exception,e:
        pass
    print "Finishing process..."
    
    #server = SocketServer.UDPServer((HOST, PORT), ThreadingUDPServer)
    #server.serve_forever()
