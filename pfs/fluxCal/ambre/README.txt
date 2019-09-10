Stellar spectra intended for use as in-field flux standards
===========================================================

Provided by Masayuki Tanaka, who writes:

> I am not sure how many spectra you need, but I put Teff=6500, 7000, 7500K
> spectra with Z=0 and -1.  File names are in this format:
>
> > Spectrum files are named as follows:
> > p2500:g+3.5:m0.0:t01:z+0.00:a+0.00.AMBRE
> > where
> > 'p' (or 's') means that radiative transfer has been solved in plan-parallel (or spherical) geometry
> > '2500' : Teff in K
> > 'g+3.5': log of surface gravity (g in units of cm/s2)
> > 'm0.0' : mass of the corresponding MARCS model atmosphere (in Msun)
> > 't0.1' : microturbulence velocity adopted in the MARCS model atmosphere and the spectral synthesis (in km/s)
> > 'z+0.00' : mean stellar metallicity (equivalent to [Metal/H] or [Fe/H])
> > 'a+0.00' : [alpha/Fe] chemical index
>
> In the fits files,
>
> > 1- wavelength (in \AA)
> > 2- relative flux (normalized to the local continuum)
> > 3- absolute flux (in erg/cm^2/s/A)
>
> Use the columns 1 and 3, and scale the fluxes to the input magnitudes
> to the simulator.  The spectra are of high resolution (R~15000), and
> so you want to smooth them to the PFS resolution.
