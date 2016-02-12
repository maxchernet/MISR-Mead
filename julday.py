#Julian day conversions

import matplotlib.dates as dates
import datetime
import numpy as np

#Convert from julian day to a date as string
def jul2date(jul_date):
    date1 = dates.num2date(dates.julian2num(jul_date))
    return '%d.%d.%d'%(date1.year, date1.month, date1.day)

#Convert from julian day to a date as date
def jul2date2(jul_date):
    date1 = dates.num2date(dates.julian2num(jul_date))
    return date1

#Convert day of year to date if year is known
def doy2date(year, doy):
    date_0 = datetime
    jul_0 = dates.num2julian( dates.date2num(date_0.date(year,1,1)) )
    date = dates.num2date(dates.julian2num(jul_0+doy-1)) #jul2date(jul_0+doy-1)
    return date

#Convert date to julian day
def date2jul(date):
        return dates.num2julian( dates.date2num(date) )

# Convert date to day of year
def jul2doy(jul_date):
        y = jul2date2(jul_date).year
        d = datetime
        jul_date_0 = date2jul(d.date(y, 1, 1))
        #print 'jul2doy: ', jul_date - jul_date_0 + 1
        return int(np.round(jul_date - jul_date_0)) + 1

