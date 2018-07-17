# -*- coding: utf-8 -*-

from astropy.coordinates import SkyCoord, EarthLocation, get_sun, get_moon, \
AltAz
from astroplan import FixedTarget, Observer
from astropy.table import QTable
from astropy.time import Time
from plots import visibility_plot
import astropy.units as u
import numpy as np

Time.FORMATS.setdefault('mjd')

def visibility_table(observer, target, obstime, dt=1*u.hour, dt_number=24):
    """Given a target and a date for observation returns....
    Inputs:
        target:
        date:
        observer (optional):
        dt (optional):
        dt_number (optional:
    Outpus
        table
    """
    utcoffset = observer.timezone.utcoffset(obstime.datetime).seconds*u.s
    delta_midnight = np.linspace(0, dt_number, dt_number + 1)*u.hour
    times = obstime  + delta_midnight
    frame = AltAz(obstime=times, location=observer.location)
    moon_altazs = get_moon(times).transform_to(frame)
    sun_altazs = get_sun(times).transform_to(frame)
    target_altazs = target.coord.transform_to(frame)
    night = observer.is_night(times)
    up = observer.target_is_up(times, target)
    table = QTable()
    table['date'] = times
    # table['night'] = night
    table['target'] = target_altazs
    table['visible'] = up & night
    table['airmass'] = target_altazs.secz
    table['moon_separation'] = np.sqrt((target_altazs.alt - moon_altazs.alt)**2
                                       + (target_altazs.az - moon_altazs.az)**2)
    table['moon_iluminationn'] = observer.moon_illumination(times)
    import pdb; pdb.set_trace()
    return(table)

def visibility(observer, coord, obstime, **kargs):
    """Function that calculates the visibility of an object in
    """
    # import pdb; pdb.set_trace()
    try:
        coord = SkyCoord(coord, obstime=obstime, **kargs)
    except:
        raise ValueError("Invalid coordinates")
    target = FixedTarget(coord)
    # import pdb; pdb.set_trace()
    table = visibility_table(observer, target, obstime)
    if sum(table['visible'] == True) == 0:
        print(f'Target {coord} is not visible in {date}')
        visible = False
    else:
        cond = table['airmass'] > 0
        min_airmass = min(table[cond]['airmass'])
        min_airmass_table = table['airmass'] == min_airmass
        min_airmass_time = table[min_airmass_table]['date']
        print(table[table['visible'] == True] )
        print(f'\n min air mass at {min_airmass_time}')
        visible = True
    return visible

def get_coordinates():
    """
    """
    coord = input('\nInput target coordinates (ra dec): ')
    coord_units = {'1':(u.hourangle, u.deg), '2':(u.deg, u.deg),
                   '3':(u.rad, u.rad)}
    coord_unit = 0
    while coord_unit not in  coord_units.keys():
        coord_unit = input("\nSelect coordinates units (1, 2 or 3): \n\n" +
                       ' 1 - (hour,deg) \n' + ' 2 - (deg,deg) \n' +
                       ' 3 - (rad,rad) \n\n')
    return coord, coord_units[coord_unit]

def get_date():
    f"""
    Function that ask the user for the dates of the observations nihts.
    If more than one observation nigth dates must be separated by ','.
    Date format must be {settings.time_format}. This can be modified in the
    settings.py file.
    """
    dates = input("\nInput observation date. If more thatn one use ',' between dates. ")
    dates = dates.split(',')
    times = []
    for date in dates:
        try:
            time = Time(float(date), format=settings.time_format)
        except:
            raise ValueError(f'All dates must be in {settings.time_format}' )
    return times

def main():
    #define observer located at ckoirama observatory
    location = EarthLocation.from_geodetic(-69.930544*u.deg,-24.089385*u.deg,
                                           980*u.m)
    ckoirama = Observer(location=location, name="Ckoirama",
                        timezone="Chile/Continental")
    #ask user for target coordinates
    coord, coord_unit = get_coordinates()
    #ask user for observation date
    time = get_date()
    visibility(ckoirama, coord, time, unit=coord_unit)

if __name__ == '__main__':
    main()
