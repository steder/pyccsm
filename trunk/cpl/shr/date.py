import cal

secsPerDay = 86400.0

class Date:
    def __init__(self, initdate=None):
        self.y = 0
        self.m = 0
        self.d = 0
        self.s = 0
        self.cdate = 0
        self.eday = 0
        self.stepInDay = 0
        self.stepsPerDay = 1
        return

    def initYMD(self, y,m,d, ns):
        """init date given Year, month, day, and steps per day"""
        if( cal.validYMD( y,m,d ) ):
            self.y = y
            self.m = m
            self.d = d
            self.s = 0
            self.stepInDay = 0
            self.stepsPerDay = ns
            self.cdate = cal.ymd2date( y,m,d )
            self.eday = cal.ymd2eday( y,m,d )
        else:
            raise ValueError,"%s:invalid y,m,d = %s,%s,%s"%("shr.date.initYMD",y,m,d)
        return self

    def initEday(self,edays, ns):
        """init date given elapsed days"""
        y,m,d = cal.eday2ymd( edays )
        if( cal.validYMD( y,m,d )):
            self.initYMD(y,m,d,ns)
        else:
            raise ValueError,"%s:invalid eday = %s"%("shr.date.initEday",edays)
        return self

    def initCDate(self, cdate, ns):
        """init date given coded date(yymmdd)"""
        if( cal.validDate(cdate)):
            y,m,d = cal.date2ymd(cdate)
            self.initYMD(y,m,d,ns)
        else:
            raise ValueError,"%s:invalid date = %s"%("shr.date.initCDate",cdate)
        return self

    def adv1step(self):
        """
        advance the date one time step
        """
        self.stepInDay += 1
        if( self.stepInDay < self.stepsPerDay ):
            self.s = int( (secsPerDay * self.stepInDay)/self.stepsPerDay )
        else:
            self.advNextDay()
        return

    def advNextDay(self):
        """ advance date to next day, 0 seconds """
        self.eday += 1
        self.stepInDay = 0
        self.s = 0
        self.y,self.m,self.d = cal.eday2ymd( self.eday )
        self.cdate = cal.eday2date( self.eday )
        return


    def getymd( self ):
        """
        returns integers yy,mm,dd,ss
        """
        return self.y, self.m, self.d, self.s

    def getCDate(self):
        """
        returns integers cdate, seconds
        """
        return self.cdate,self.s

    def getEday(self):
        """
        returns integers eday,seconds
        """
        return self.eday,self.s

    def getStepsPerDay(self):
        return self.stepsPerDay

    def getStepInDay(self):
        return self.stepInDay

    def getJulian(self, shift=None):
        """
        return julian day (as a real number)
        """
        totsec = self.s
        if( shift ):
            totsec += shift

        julian = cal.elapsDaysStrtMonth(self.y, self.m)
        julian += self.d
        julian += totsec/secsPerDay
        return julian

    def __equal__(self, source):
        return

    def __cmp__(self, source):
        return
