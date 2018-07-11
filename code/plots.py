# -*- coding: utf-8 -*-
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.visualization import astropy_mpl_style
from astropy.coordinates import get_sun, get_moon, AltAz
plt.style.use(astropy_mpl_style)

def visibility_plot(observer, time, target):
    # import pdb; pdb.set_trace()
    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    utcoffset = observer.timezone.utcoffset(time.datetime).seconds*u.s
    utcoffset = utcoffset.to(u.hour)
    times = time - utcoffset + delta_midnight
    frame = AltAz(obstime=times, location=observer.location)
    sun_altazs = get_sun(times).transform_to(frame)
    moon_altazs = get_moon(times).transform_to(frame)
    # import pdb; pdb.set_trace()
    target_altazs = target.coord.transform_to(frame)

    plt.plot(delta_midnight, sun_altazs.alt, color='r', label='Sun')
    plt.plot(delta_midnight, moon_altazs.alt, color=[0.75]*3, ls='--', label='Moon')
    plt.scatter(delta_midnight, target_altazs.alt,
                c=target_altazs.az, label='target', lw=0, s=8,
                cmap='viridis')
    plt.fill_between(delta_midnight.to('hr').value, 0, 90,
                     sun_altazs.alt < -0*u.deg, color='0.5', zorder=0)
    plt.fill_between(delta_midnight.to('hr').value, 0, 90,
                     sun_altazs.alt < -18*u.deg, color='k', zorder=0)
    plt.colorbar().set_label('Azimuth [deg]')
    plt.legend(loc='upper left')
    plt.xlim(-12, 12)
    plt.xticks(np.arange(13)*2 -12)
    plt.ylim(0, 90)
    plt.xlabel('Hours from EDT Midnight')
    plt.ylabel('Altitude [deg]')
    plt.show()
