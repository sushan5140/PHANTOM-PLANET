README

Each Keck/HIRES spectrum consists of three individual files. Each file is 
prefixed with an letter indicating the HIRES chip. Each chip covers the 
following range of wavelengths:

'i' (6540-8000 Ang)
'r' (4975-6420 Ang)
'b' (3640-4800 Ang)
The CPS radial velocity analysis uses only the central ('r') spectrum, where
the iodine absorption lines exisit (roughly 4000-6000 Angstroms.)

Within each .fits file, there is a single extension. This is the value of dn
per pixel. In order to obtain the flux or photons per pixel, multiply each 
pixel by the gain, which is equal to 2.0 for all three HIRES CCDs.

The blaze function has NOT been removed from each spectral order.


The wavelength solution and explanation, accurate to ~1 pixel can be found here:
https://www.astro.caltech.edu/~howard/cps/hires_cps.html

The CKS spectra are defined with respect to the rest frame of HIRES. 
For precision spectroscopy, it is more convenient to register (shift) spectra 
onto the rest frame of the star, i.e. H-alpha at 6563 angstroms. The open source
code Empirical-SpecMatch, developed by Yee, Petigura, and von Braun (2017) 
contains a module to facilitate the registration of HIRES spectra. For example,
shifting registering KOI-7 (systemic RV = −60.8 km/s).

The spectra contained here have not been shifted to the stars' rest frames.

Associated data calculated with these spectra includes precise relative
radial velocity measurements of the target stars.

Documentation about calculating radial velocities can be found here:
https://ui.adsabs.harvard.edu/abs/1996PASP..108..500B/abstract

And information on the precise radial velocity observing setup can be found here:
https://ui.adsabs.harvard.edu/abs/2010ApJ...721.1467H/abstract

TABLE DATA:
Columns and definitions:
index               Incremental
star                Starname
ra                  Right_Ascension(Hours)
dec                 Declinatin(Degrees)
observation_id      Filename
bjd                 Barycentric Julian Date(days)
mnvel               Precise Radial Velocity (m/s)
errvel              Precise Radial Velocity Error (m/s)
svalue              Activity metric: Mt. Wilson Style S-value
svalue_err          S-value error(absolute)
