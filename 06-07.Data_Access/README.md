# Access and Customize ICESat-2 Data Tutorials

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/ICESAT-2HackWeek/data-access/blob/master/LICENSE)

## Presenters and Authors

Jessica Scheick ([@JessicaS11](https://github.com/JessicaS11)), University of Maine, lead `icepyx` developer

Amy Steiker ([@asteiker](https://github.com/asteiker)), National Snow and Ice Data Center Distributed Active Archive Center (NSIDC DAAC)


## Summary

These tutorials will walk you though how to access ICESat-2 data, with an emphasis on programmatic access from the NASA National Snow and Ice Data Center Distributed Active Archive Center (NSIDC DAAC) using the [`icepyx`](https://github.com/icesat2py/icepyx) Python library. Access parameters include spatial and temporal filters, as well as customization services including subsetting and reformatting.

This tutorial is hosted on the [Github ICESat-2 HackWeek organization](https://github.com/ICESAT-2HackWeek/data-access).


## Key Learning Objectives

1. Understand ICESat-2 data access methods, using OpenAltimetry for data discovery and encouraging programmatic access to NSIDC resources.
2. Become familiar with the `icepyx` library, which provides a wrapper to NSIDC's Application Programming Interface (API).
3. Create an icesat2data object (the primary class of `icepyx`) and use it to search, order (with and without subsetting), and download ICESat-2 data.


## Schedule (all times in PDT)

*Date:* Friday 11 June

- 12:10-12:30 [Amy Steiker]

[ICESat-2 Data Access Resource Intro](https://github.com/ICESAT-2HackWeek/data-access/blob/master/notebooks/ICESat-2_NSIDC_DAAC_Data_Resource_Intro.ipynb)

*Date:* Monday 14 June

- 9:00-9:10 [Jessica Scheick]: Intro to icepyx
- 9:10-9:25 [Amy Steiker]: Search Parameters
- 9:25-9:40 [Jessica Scheick]: Query and login
- 9:40-9:50 [Amy Steiker]: Data ordering
- 9:50-9:55 [Jessica Scheick]: Data download
- 9:55-10:00 [Jessica Scheick]: Wrap up

[ICESat-2 Data Programmatic Data Access](https://github.com/ICESAT-2HackWeek/data-access/blob/master/notebooks/ICESat-2_NSIDC-DAAC_DataAccess.ipynb)

[Extra tutorial notebook on subsetting](https://github.com/ICESAT-2HackWeek/data-access/blob/master/notebooks/ICESat-2_NSIDC_DataAccess2_Subsetting.ipynb)


## Setup
1. Login to JupyterHub:
    https://icesat-2.hackweek.io
2. Open a Terminal:
    Find Terminal icon at the bottom (or click the + sign on the top left)
3. Clone GitHub repo (type on the terminal):
    `git clone https://github.com/ICESAT-2HackWeek/data-access.git`

## Questions

Post your question in the [#questions](https://icesat2hackweek.slack.com/archives/C014V14KA3G) channel on Slack.

## Credit

This collection of tutorials is conceptually based on the [2019 ICESat-2 Hackweek Data Access tutorial](https://github.com/ICESAT-2HackWeek/data-access/tree/v2019) by [@asteiker](https://github.com/asteiker) and [@wallinb](https://github.com/wallinb), which were modified into [example notebooks](https://github.com/icesat2py/icepyx/tree/master/doc/examples) to showcase the usage of the `icepyx` library. The materials for 2020 were co-developed from those examples by the presenters.

## License

All code and content in this repository is free software: you can redistribute it and/or modify it under the terms of the BSD License. A copy of this license is provided in [LICENSE](LICENSE).
