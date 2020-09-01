[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3966463.svg)](https://doi.org/10.5281/zenodo.3966463)

# 2020 ICESat-2 Hackweek Tutorials
Combined repository for the final tutorial material presented at the 2020 ICESat-2 Cryosphere-themed Hackweek presented virtually by the University of Washington from 8-18 June 2020.

## Background
The [ICESat-2 Cryospheric Science Hackweek](https://icesat-2hackweek.github.io/learning-resources/) was the first virtual Hackweek held by the University of Washington. While originally planned as a five-day, in-person workshop, the event was shifted to a fully virtual/remote setting in light of stay-at-home orders and travel restrictions in place to curb the spread of COVID-19. To accomodate multiple time zones and limit the daily duration of online tutorial sessions, the event was spread out over the course of ten days. The first week had three half-days of interactive tutorials/lectures. The second week had four days that included some interactive tutorials/lectures and scheduled times where instructors were available to help participants with a facilitated exploration of datasets and hands-on software development.

Participants learned about the ICESat-2 satellite, sensors, and datasets as well as technologies and tools for accessing and processing ICESat-2 data with a focus on the cryosphere. They also self-organized into teams for hacking projects involving ICESat-2 data.

These tutorials were largely developed by volunteer instructors. Each tutorial was prepared and distributed via a topical repository under the [ICESat-2 Hackweek Github organization](https://github.com/ICESAT-2HackWeek). Participants followed the live or recorded video of each tutorial, and had the option to clone the repository and interactively run the examples on the Pangeo JupyterHub environment created explicitly for the event. 

This `2020_ICESat-2_Hackweek_Tutorials` repository contains the final collection of Jupyter notebooks and slides presented at the virtual event. It centralizes the final content from the individual tutorial repositories and provides a tagged "release" of the material presented during the hackweek with a DOI for distribution to the larger community. Most notebooks were rendered to include all output (including embedded plots). These notebooks can be identified by the "\_rendered" at the end of the filename. Some of these tutorials may continue to evolve within their respective repositories (links are below). 

## Running these tutorials
During the Hackweek participants worked in a [Pangeo](https://pangeo.io/) environment specifically created for the event and hosted on AWS us-west-2. The JupyterHub included a shared volume with datasets used during some tutorials that are too large to be included within this repository. However, these data are publicly available from [NSIDC](https://nsidc.org/data/icesat-2) and can be easily obtained using the [icepyx library](https://icepyx.readthedocs.io/en/latest/). Where possible, code to download the needed data using `icepyx` has been included within each tutorial notebook.

**Please note that the tutorials presented here used version 0.2.0 of `icepyx`. These tutorials are set up to use that version in the provided Binder link, but they will require modification for more recent versions (>= v0.3.0) of `icepyx`. Up-to-date data access tutorials/examples are available from the `icepyx`  [repository](https://github.com/icesat2py/icepyx) and [associated documentation](https://icepyx.readthedocs.io/en/latest/getting_started/example_link.html).**

### Re-create the ICEsat-2 Hackweek JupyterLab environment with Binder
Clicking this button will launch a [binder](https://mybinder.org/) replica of the [JupyterLab computing environment](https://github.com/ICESAT-2HackWeek/jupyterhub-2020) described above. With the exception of those tutorials denoted with an asterisk(\*), this will allow you to run the tutorials presented during the Hackweek. Be aware the session is ephemeral. **Your home directory will not persist, so use this binder only for running tutorials or other short-lived demos!**

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ICESAT-2HackWeek/2020_ICESat-2_Hackweek_Tutorials/binder?urlpath=git-pull%3Frepo%3Dhttps%253A%252F%252Fgithub.com%252FICESAT-2HackWeek%252F2020_ICESat-2_Hackweek_Tutorials%26urlpath%3Dlab%252Ftree%252F2020_ICESat-2_Hackweek_Tutorials%252F%26branch%3Dbinder)

## Tutorials
\* Tutorial filenames within this repository denoted below with an asterisk(\*) cannot be run in Binder due to large dataset requirements.

### 01. Introductory Session (slides)
*Anthony Arendt and Charley Haley*
* [Slides](https://github.com/ICESAT-2HackWeek/2020_ICESat-2_Hackweek_Tutorials/blob/master/01.IntroductorySession_slides.pdf)

### 02. ICESat-2 Mission: Satellite, Sensor, and Data
*Axel Schweiger*
* Welcome by the NASA Cryosphere Program Manager (*Thorsten Markus*)
* intro to the ICESat-2 mission (*Tom Neumann*)
    * [Slides](https://github.com/ICESAT-2HackWeek/intro_ICESat2/blob/master/HackWeekIntroNeumann2020.pptx)
* ICESat-2 data products (*Ben Smith*)
    * [Slides](https://github.com/ICESAT-2HackWeek/intro_ICESat2/blob/master/ICESat-2_data_products_Hackweek2020.pptx)

### 03. Git and GitHub
*Fernando Perez*
* [intro-git repository](https://github.com/ICESAT-2HackWeek/intro-git)
* [Slides](https://docs.google.com/presentation/d/1pOWte7V5UbnVBvRktvLbLTRluDwrGbXtIdAZhzAd1AE/edit?usp=sharing)
* [Video](https://youtu.be/O2lLC_s_d20)

### 04. Jupyter and iPython
*Fernando Perez*
* [intro-jupyter repository](https://github.com/ICESAT-2HackWeek/intro-jupyter)
* [Slides](https://docs.google.com/presentation/d/1TfY7rnCuGQDrlvsf2-P9lNADT2vwiJsdb7vmgZ3SDmA/edit?usp=sharing)
* [Video](https://youtu.be/Jft9-RnmH1Y)

### 05. Geospatial Analysis with Python
*David Shean*
* [geospatial-analysis repository](https://github.com/ICESAT-2HackWeek/geospatial-analysis)
* [Video](https://youtu.be/46vxJYqUMsM)

### 06. Introduction to ICESat-2 Sea Ice and Land Ice Products and Data Access
* Sea ice products: overview of products, algorithms, and parameters for sea ice investigations (*Alek Petty*)
    * [sea ice lesson slides](https://drive.google.com/file/d/1e3VFvBRBHcY5_gjEyWVjA-l7tL2K4HfQ/view?usp=sharing)
* Land ice products: overview of products, algorithms, and parameters for land ice investigations (*Ben Smith*)
* ICEsat-2 data access: basic data explore and visualization in OpenAltimetry (*Jessica Scheick and Amy Steiker*)
    * [ICESat-2 data access repository](https://github.com/ICESAT-2HackWeek/data-access)
* [Video](https://youtu.be/UH9t3Fn7lN0)

### 07. Programmatic ICESat-2 data access
*Jessica Scheick and Amy Steiker*
* [ICESat-2 data access repository](https://github.com/ICESAT-2HackWeek/data-access)
* [Video](https://youtu.be/RJNSDVnfGvU)

### 08. Introduction to HDF5 and ICESat-2 data files
*Fernando Paolo*
* [Introduction to HDF5 and ICESat-2 data files repository](https://github.com/ICESAT-2HackWeek/intro-hdf5)
* [Video](https://youtu.be/-sEDiMEne3k)

### 09. Land ice applications\*
*Ben Smith*
* [Land ice applications repository](https://github.com/ICESAT-2HackWeek/Land_Ice_Applications)
* [Video](https://youtu.be/qtkVd2xc-U8)

### 10. Sea ice applications\*
*Ellen Buckley*
* [Sea ice applications repository](https://github.com/ICESAT-2HackWeek/sea-ice-tutorials)
* [Video](https://youtu.be/zl-qifVcPV4)

 ### 11. Science data generation\*
 *Johan Nilsson*
 * [Science data generation repository](https://github.com/ICESAT-2HackWeek/ScienceDataGeneration)
 * [Video](https://www.youtube.com/watch?v=PZlGcI3SMts)
 
### 12. Machine learning
*Yara Mohajerani and Shane Grigsby*
* [Machine learning repository](https://github.com/ICESAT-2HackWeek/Machine-Learning)
* [Video](https://www.youtube.com/watch?v=GQ2RSOSOmdU)

## Citation
This content is original material prepared and/or modified for the Hackweek by a dedicated team of volunteer instructors. We released these materials with a digital object identifier (DOI) to provide an easy way for both contributors and users to cite this material, as it is not necessarily appropriate for a peer-reviewed journal article publication. If you find these tutorials useful and/or adapt some of the underlying source code for your research (whether or not you attended the Hackweek), we request that you *Star* the repository (clicking the button in upper right corner) and cite as:

Anthony Arendt, Jessica Scheick, David Shean, Ellen Buckley, Shane Grigsby, Charley Haley, Lindsey Heagy, Yara Mohajerani, Tom Neumann, Johan Nilsson, Thorsten Markus, Fernando Paolo, Fernando Perez, Alek Petty, Axel Schweiger, Ben Smith, Amy Steiker, Sebastian Alvis, Scott Henderson, Nick Holschuh, Zheng Liu, Tyler Sutterley. (2020). ICESAT-2HackWeek/2020_ICESat-2_Hackweek_Tutorials (Version 1.0). Zenodo. http://doi.org/10.5281/zenodo.3966463.

With such a large number of contributors at varying career stages (with many building publication records) collaboratively working on a diverse array of tasks for this event, author order is somewhat subjective. The author order listed above was determined by:
1. Involvement in planning the Hackweek and in preparing this repository for release
2. Alphabetical order by last name for tutorial leads and presenters
3. Alphabetical order by last name for other planning team members who did not present

Special thanks to Yu-Chan Chao (UW APL) and Jane Koh (UW eScience) for all of their technical expertise and help in planning and executing the event!

Please click on the Zenodo badge [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3966463.svg)](https://doi.org/10.5281/zenodo.3966463) for the latest citation information and export options.

## License
The content of this project is licensed under the [Creative Commons Attribution 3.0 Unported license](https://creativecommons.org/licenses/by/3.0/), and the underlying source code used to format and display that content is licensed under the [MIT license](LICENSE.md).
