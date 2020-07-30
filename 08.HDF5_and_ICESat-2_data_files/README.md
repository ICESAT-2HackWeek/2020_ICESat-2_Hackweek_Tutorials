# Intro to HDF5 and ICESat-2 Data Files

Instructor: Fernando Paolo (paolofer@jpl.nasa.gov)  

## Tutorials:  

Part 1: Introduction to HDF5 data model  
Part 2: Reduction of ICESat-2 data files  

## Goals:  

- Familiarize with HDF5 data model  
- Familiarize with HDF5 basic tools  
- Have a quick overview of cloud options  
- Download IS2 files according region of interest  
- Extract variables of interest and filter in time and space  
- Prepare data for large-scale processing  
- Learn simple/generic data parallelization strategy  
- Inspect data files from the command line  
- Plot millions of points efficiently

## Schedule

Monday June 15: Hack-day #1 (10:30 - 11:30 AM)

## Setup

1. Login to JupyterHub: https://icesat-2.hackweek.io/hub/login
2. Open a Terminal: Find Terminal icon at the bottom (or click the `+` sign on the top left)
3. Clone GitHub repo (type on the terminal): `git clone https://github.com/ICESAT-2HackWeek/intro-hdf5.git`

## Questions

Post your question in the [#questions](https://icesat2hackweek.slack.com/archives/C014V14KA3G) channel on Slack

## Credit

The algorithms used in this tutorial are downscaled versions of an open-source Python altimetry package, currently being developed at JPL/Caltech: [captoolkit](https://github.com/fspaolo/captoolkit). This package provides a range of algorithms for common tasks in altimetry data processing. We also use the convenient data download interface from [icepyx](https://github.com/icesat2py/icepyx). 

## License

The content of this project itself is licensed under the [Creative Commons Attribution 3.0 Unported license](https://creativecommons.org/licenses/by/3.0/), and the underlying source code used to format and display that content is licensed under the [MIT licenseo](https://github.com/github/choosealicense.com/blob/gh-pages/LICENSE.md).
