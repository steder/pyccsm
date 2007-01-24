"""
MODULE: shr_cal_mod -- calendar module, relates elapsed days to calendar dat\e.

DESCRIPTION:
   These calendar routines do conversions between...
   \begin{itemize}
   \item the integer number of elapsed days
   \item the integers year, month, day (three inter-related integers)
   \item the integer coded calendar date (yyyymmdd)
   \end{itemize}
   Possible uses include: a calling routine can increment the elapsed days
   integer and use this module to determine what the corresponding calendar
   date is;  this module can be used to determine how many days apart two
   arbitrary calendar dates are.

REMARKS:
Following are some internal assumptions.  These assumptions are somewhat
arbitrary -- they were chosen because they result in the simplest code give\n
the requirements of this module.  These assumptions can be relaxed as
necessary:
o the valid range of years is [0,9999]
o elapsed days = 0 <=> January 1st, year 0000
o all years have 365 days (no leap years)
  This module is hard-coded to implement a 365-day calendar, ie. there
  are no leap years.  This module can be modified to implement a calendar
  with leap years if this becomes desireable.  This would make the internal
  logic of this module more complex, but would not require any change to th\e
  module API or the calling code because the module API hides these details
  from all external routines.
"""

# Elapsed *D*ays on the *S*tart of the *M*onth
dsm = [0,31,59,90,120,151,181,212,243,273,304,334]
# Days Per Month
dpm = [31,28,31,30,31,30,31,31,30,31,30,31]

def eday2date( eday ):
    """
    Converts elapsed days to coded-date
    """
    eday = int(eday)
    year = eday / 365
    day = eday % 365
    for k in xrange(12):
        if day >= dsm[k]:
            month = k
    day = day - dsm[month] + 1

    return ( (year*10000) + (month*100) + day )

def eday2ymd( eday ):
    """
    Converts elapsed days to year/month/day 

    Assumptions:
      This calendar has a year zero (but no day or month zero)
    """
    eday = int(eday)
    year = eday/365
    day = eday % 365
    for k in xrange(12):
        if( day >= dsm[k] ):
            month = k
    day = day - dsm[k] + 1
    return year,month,day
    
def date2ymd( cdate ):
    """
    Converts coded-date to year/month/day
    """
    cdate = int(cdate)
    if (not validDate( cdate )):
        raise ValueError,'%s:Invalid date: %s'%("shr.cal.date2ymd",cdate)
    year =  int( cdate / 10000 )
    month = int( (cdate%10000)/100 )
    day = cdate % 100
    return year,month,day

def date2eday( cdate ):
    """
    Converts coded-date to elapsed days
    """
    cdate = int(cdate)
    if( not validDate( cdate )):
        raise ValueError,'%s:Invalid date: %s'%("shr.cal.date2eday",cdate)
    y,m,d = date2ymd(cdate)
    return (y * 365) + dsm[m] + (day - 1)

def ymd2date( year, month, day ):
    """
    converts year, month, day to coded-date
    """
    if( not validYMD( year, month, day ) ):
        raise ValueError,'%s:Invalid date: %s/%s/%s'%("shr.cal.ymd2date",year,month,day)
    date = (year * 10000) + (month * 100) + day
    return date

def ymd2eday( year, month, day ):
    """
    Converts year, month, day to elapsed days
    """
    if( not validYMD( year, month, day ) ):
        raise ValueError,'%s:Invalid date: %s/%s/%s'%("shr.cal.ymd2eday",year,month,day)
    eday = (year * 365) + dsm[month] + (day - 1)
    return eday

def validDate( date ):
    """
    Determines if the given coded-date is a valid date
    """
    y =  int( date / 10000 )
    m = int( (date%10000)/100 )
    d = date % 100
    return validYMD(y,m,d)
    
def validYMD( y, m, d ):
    """
    Determines if the given year, month, and day indicates a valid date
    """
    valid = True
    if not ( 0 <= y <= 9999 ):
        valid = False
    if not ( 1 <= m <= 12 ):
        valid = False
    if not ( 1 <= d <= dpm[m] ):
        valid = False
    return valid
    
def numDaysInMonth( year, month ):
    """
    Return the number of days in a month.

    ! Appears to ignore leap years !
    """
    return dpm[month]

def elapsDaysStrtMonth( year, month ):
    """
    return the number of elapsed days at the start of a month

    !Again, appears to ignore leap years!
    """
    return dsm[month]
