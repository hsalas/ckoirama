# -*- coding: utf-8 -*-
import settings
import numpy as np
import astropy.units as u

from astropy.time import Time
from astroplan import Observer
from astropy.table import Table
from visibility import visibility, get_targets
from astropy.coordinates import EarthLocation


def main():
    #define observer located at ckoirama observatory
    location = EarthLocation.from_geodetic(-69.930544*u.deg,-24.089385*u.deg,
    980*u.m)
    ckoirama = Observer(location=location, name="Ckoirama",
                        timezone="Chile/Continental")
    time = Time('2018-05-06 12:00:00')


if __name__ == '__main__':
    main()
