import time

class Stopwatch:
    def __init__(self, name="unnamed"):
        self.name = name
        self.starttime=0
        self.stoptime=0
        self.isrunning = False
        self.isreset = True
        
    def start(self):
        if( not self.isrunning ):
            if( self.isreset ):
                self.starttime=time.time()
                self.stopttime=self.starttime
                self.isrunning = True
                self.isreset = False
            else:
                self.isrunning = True
        
    def __str__(self):
        if( self.isrunning ):
            currenttime = time.time()
            return "elapsed time on timer %s: %s"%( self.name, currenttime-self.starttime )
        else:
            return "elapsed time on timer %s: %s"%( self.name, self.stoptime-self.starttime)

    def __repr__(self):
        return self.__str__()

    def stop(self):
        if( self.isrunning ):
            self.stoptime=time.time()
            self.isrunning = False
            
    def getElapsedTime(self):
        if( self.isrunning ):
            currenttime = time.time()
            return (currenttime - self.starttime)
        else:
            return (self.stoptime-self.starttime)

    def isRunning(self):
        return self.isrunning

    def reset( self ):
        self.starttime = time.time()
        self.stoptime = time.time()
        self.isreset = True

class Timer:
    def __init__( self, names=[] ):
        self.timers = {}
        self.names = names
        for name in self.names:
            self.timers[name] = Stopwatch(name)

    def start(self, *names):
        for name in names:
            self.timers[name].start()

    def stop(self, *names):
        for name in names:
            self.timers[name].stop()

    def reset(self, *names):
        for name in names:
            self.timers[name].reset()

    def __str__(self):
        s = ""
        s += "Data from timer for:\n"
        s += "*"+str(self.timers.keys())+"\n"
        for name in self.names:
            s += "\t*%s"%(self.timers[name]) + "\n"
        return s
    
    def __repr__(self):
        return self.__str__()

