# Response functions

These files describe the response functions of a variety of filters.
In each file, the first column is the wavelength (in nm) in vacuum (\*1),
and the second column is the response function.
The response functions include miscellaneous effects
such as those of the telescope, the CCD and the atmosphere.

(\*1): Originally, wavelengths in all files in this directory were in air.
The wavelengths have been converted to those in vacuum by solving
`wavelength_in_air = vacuum_to_air(wavelength_in_vacuum)` with an iterative method.
The function `vacuum_to_air(wavelength_in_vacuum)` is the conversion
given by Peck & Reeder (1972):
```
def vacuum_to_air(wavelength_in_vacuum):
    ss = (1e3 / wavelength_in_vacuum)**2
    n = 1.0000806051 + 0.02480990 / (132.274 - ss) + 0.000174557 / (39.32957 - ss)
    return wavelength_in_vacuum / n
```
This vacuum-to-air conversion is used by
US National Institute of Standards and Technology
according to Morton (2003) (Bibcode: 2003ApJS..149..205M).

## HSC

Computed by Masayuki Tanaka.

The response functions include
the primary mirror,
prime focus optics,
filter (averaged over the surface area),
CCD response,
and atmosphere (airmass=1.2 PVM=1.5mm).

## SDSS

Computed by Masayuki Tanaka.

The response functions include
the telescope (2 aluminum surfaces),
camera,
filter (average of the 6 columns),
CCD response,
and atmosphere (airmass=1.3).

## PS1

Quoted from Table 3 in J. L. Tonry et al 2012 ApJ 750 99.

## Gaia

"Gaia DR2 revised passbands and zeropoints"
https://www.cosmos.esa.int/web/gaia/iow_20180316
