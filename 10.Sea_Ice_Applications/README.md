# ICESAT-2HackWeek-sea-ice-tutorials

Authors: Ellen Buckley and Alek Petty 

Date: June 2020

ICESat-2 hackweek repository for the sea ice tutorials. Primarily hosts the Jupyter notebooks which contain the individual tutorials. 


![icesat2_seaice](./Images/icesat2_seaice.png?raw=true "ICESat-2 profiling the sea ice surface, taken from the ICESat-2 website (Satellite image courtesy of Orbital
Earth image illustrating AMSR-E sea ice courtesy of the NASA Scientific Visualization Studio)")



## Setup

This GitHub repository primarily hosts the Jupyter Notebooks needed for the hackweek tutorials. The notebooks should work without any extra steps if you're working in the ICESat-2 Pangeo environment that has been created for this hackweek. Just clone this repository into your Pangeo user account (git clone https://github.com/ICESAT-2HackWeek/sea-ice-tutorials when logged in). 

The example data files are being stored in the shared folder on Pangeo '/home/jovyan/tutorial-data/sea-ice/'. 


## Notebooks

0. DataAccess.ipynb
* run this notebook to download the ICESat-2 granules required for the tutorials

1. ATL03.ipynb
* General understanding of the data included in a typical ATL03 file.
* Reading in, plotting and basic analysis of ATL03 data.
* Understanding along-track distance and time variables.
* Understanding of differences between weak and strong beams

2. ATL07.ipynb
* General understanding of the data included in a typical ATL07 file.
* Reading in, plotting and basic analysis of ATL07 data.
* How variables change with different surface types.
* How clouds affect surface returns

3. ATL10.ipynb
* General understanding of the difference between ATL07 and ATL10
* Reading in, plotting and basic analysis of ATL10 data.
* How you too can calculate sea ice thickness from freeboard!

4. GriddingDemo.ipynb
* demo binning of ICESat-2 along-track data to the 25 km NSIDC grid.


## Background

NASA's ICESat-2 launched succesfully in September 2018 from Vadenberg Airforce Base, California. ICESat-2 carries onboard a single instrument – the Advanced Topographic Laser Altimeter System (ATLAS). ATLAS measures the travel times of laser pulses to derive the distance between the spacecraft and Earth’s surface and thus an estiamte of the Earth's surface elevation with respect to a reference surface. The orbit pattern results in a higher density of shots within the polar regions. One of ICESat-2's primary mission objectives is to: *Estimate sea-ice thickness to examine ice/ocean/atmosphere exchanges of energy, mass and moisture.*

![icesat2_profiling](./Images/icesat2_profiling.png?raw=true "ICESat-2 profiling the sea ice surface, figure taken from the ATL07/10 ATBD document")

ICESat-2 employs a photon-counting system to obtain high measurement sensitivity with lower resource (power) demands on the satellite platform compared to the analog waveform approach of the ICESat laser altimeter GLAS. The ATLAS instrument transmits laser pulses at 532 nm with a high pulse repetition frequency of 10 kHz. The ICESat-2 nominal orbit altitude of ~500 km results in laser footprints of ~17 m on the ground, which are separated by only ~0.7 m along-track (resulting in substantial overlap between the shots). This relatively small footprint size combined with the high pulse repetition frequency and  precise time of flight estimates were partly chosen to enable the measurements of sea surface height within leads to carry a precision of 3 cm or less. 

The laser is split into 6 beams (three pairs of strong and weak beams) which provide individual profiles of elevation. The multiple beams address the need for unambiguous separation of ice sheet slope from height changes. For sea ice, this provides multiple profiles of sea ice and sea surface heights, increasing overall profiling coverage and enabling assessments of beam reliability. 

The beam configuration and their separation are shown above: the beams within each pair have different transmit energies (‘weak’ and‘strong’, with an energy ratio between them of approximately 1:4) and are separated by 90 m in the across-track direction. The beam pairs are separated by ~3.3 km in the across-track direction, and the strong and weak beams are separated by ~2.5 km in the along-track direction. The observatory orientation is an important consideration, as this changes the labelling of the beams - i.e. 'gt1r' refers to the left-side strong beam when in the forward direction (beam 5) but the left-side weak beam when in the backward direction (beam 2). 

The ICESat-2 products of most interest to the sea ice community are:

* ATL03: Along-track photon heights (1. ATL03.ipynb tutorial) 
* ATL07: Along-track segment surace heights (2. ATL07.ipynb tutorial) 
* ATL09: Cloud products (no direct tutorial provided, used mainly in ATL07 production for cloud filtering)
* ATL10: Along-track segment (and 10 km swath) freeboards (3. ATL10.ipynb tutorial) 
* ATL20: Gridded monthly sea ice freeboard (expected release summer 2020).

We provide in the notebooks a brief summary of these data products, but encourage the user to read the ATBD or references provided above and at the start of the Jupyter Notebooks for the more complete (and probably accurate) descriptions.

The ICESat-2 data products are provided in the Hierarchical Data Format – version 5 (HDF-5) format and have recently been made publicly available through the National Snow and Ice Data Center (NSIDC - https://nsidc.org/data/icesat-2). See the hdf5 tutorial (https://github.com/ICESAT-2HackWeek/intro-hdf5) for more information on this data format.


## References

Kwok, R., Markus, T., Kurtz, N. T., Petty, A. A., Neumann, T. A., Farrell, S. L., et al. (2019). Surface height and sea ice freeboard of the Arctic Ocean from ICESat-2: Characteristics and early results. Journal of Geophysical Research: Oceans, 124, doi: 10.1029/2019JC015486.

Markus, T., Neumann, T., Martino, A., Abdalati, W., Brunt, K., Csatho, B., et al. (2017). The Ice, Cloud and land Elevation Satellite-2 (ICESat-2): Science requirements, concept, and implementation. Remote Sensing of the Environment, 190, 260-273, doi: 10.1016/j.rse.2016.12.029.

Neumann, T A., Martino, A. J., Markus, T., Bae, S., Bock, M. R., Brenner, A. C., et al. (2019). The Ice, Cloud, and Land Elevation Satellite – 2 mission: A global geolocated photon product derived from the Advanced Topographic Laser Altimeter System. Remote Sensing of Environment, 233, 111325, doi: 10.1016/j.rse.2019.111325.

Petty, A. A., N. T. Kurtz, R. Kwok, T. Markus, T. A. Neumann (2020), Winter Arctic sea ice thickness from ICESat‐2 freeboards, Journal of Geophysical Research: Oceans, 125, e2019JC015764. doi: 10.1029/2019JC015764.