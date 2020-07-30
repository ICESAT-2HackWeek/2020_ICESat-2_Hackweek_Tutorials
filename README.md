[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3966463.svg)](https://doi.org/10.5281/zenodo.3966463)

# 2020 ICESat-2 Hackweek Tutorials
Combined repository for the final tutorial material presented at the 2020 ICESat-2 Cryosphere-themed Hackweek presented virtually by the University of Washington from 8-18 June 2020.


## Background
The [ICESat-2 Cryospheric Science Hackweek]() was the first virtual Hackweek held by the University of Washington. Originally meant to be a five-day, in-person workshop, organizers quickly regrouped to make the event virtual in light of stay-at-home orders and travel restrictions in place to curb the spread of COVID-19. To accomodate multiple time zones and limit the daily duration of online tutorial sessions, the event was spread out over the course of ten days. The first week had three half-days of interactive tutorials/lectures. The second week had four days that included some interactive tutorials/lectures and scheduled times where instructors were available to help participants with a facilitated exploration of datasets and hands-on software development.

Participants learned about the ICESat-2 satellite, sensors, and datasets as well as technologies and tools for accessing and processing ICESat-2 data with a focus on the cryosphere. They also participated in hacking projects.

These tutorials were largely developed by volunteer instructors. Each tutorial was prepared and distributed via a topical repository under the [ICESat-2 Hackweek Github organization](https://github.com/ICESAT-2HackWeek). Participants were welcomed to either follow along on the video presentation or clone the repository and run it locally on the provided Pangeo JupyterHub environment created explicitly for the event. This 2020_ICESat-2_Hackweek_Tutorials repository contains the final collection of tutorials presented at the virtual event. It centralizes the final content from the individual tutorial repositories and provides a tagged "release" of the material presented during the hackweek with a DOI for distribution to the larger community. Some of these tutorials are may continue to evolve within their respective repositories (links are below).


## Re-create the icesat2 hackweek JupyterLab environment with Pangeo Binder

## Tutorials
### 01. Introductory Session (slides)
*Anthony Arendt and Charley Haley*

* [Slides](https://docs.google.com/presentation/d/1kNc6u4mz9qt5TI-DCosSL6jZ5M7q_k3godlBrOE891c/edit?usp=sharing)

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

### 04. Jupyter and iPython
*Fernando Perez*

* [intro-jupyter repository](https://github.com/ICESAT-2HackWeek/intro-jupyter)
* [Slides](https://docs.google.com/presentation/d/1TfY7rnCuGQDrlvsf2-P9lNADT2vwiJsdb7vmgZ3SDmA/edit?usp=sharing)


### 05. overview of python / numpy / pandas / matplotlib / geospatial data processing
*David Shean*

* [geospatial-analysis repository](https://github.com/ICESAT-2HackWeek/geospatial-analysis)

### 06. Introduction to ICESat-2 Sea Ice and Land Ice Products and Data Access

* Sea ice products: overview of products, algorithms, and parameters for sea ice investigations (*Alek Petty*)

    * [sea ice lesson slides](https://drive.google.com/file/d/1e3VFvBRBHcY5_gjEyWVjA-l7tL2K4HfQ/view?usp=sharing)

* Land ice products: overview of products, algorithms, and parameters for land ice investigations (*Ben Smith*)

* ICEsat-2 data access: basic data explore and visualization in OpenAltimetry (*Jessica Scheick and Amy Steiker*)

    * [ICESat-2 data access repository](https://github.com/ICESAT-2HackWeek/data-access)

### 07. Programmatic ICESat-2 data access
*Jessica Scheick and Amy Steiker*

* [ICESat-2 data access repository](https://github.com/ICESAT-2HackWeek/data-access)
 
### 08. Introduction to HDF5 and ICESat-2 data files
*Fernando Paolo*

* [Introduction to HDF5 and ICESat-2 data files repository](https://github.com/ICESAT-2HackWeek/intro-hdf5)

### 09. Land ice applications
*Ben Smith*

* [Land ice applications repository](https://github.com/ICESAT-2HackWeek/Land_Ice_Applications)

### 10. Sea ice applications
*Ellen Buckley*
 
* [Sea ice applications repository](https://github.com/ICESAT-2HackWeek/sea-ice-tutorials)
 
 ### 11. Science data generation
 *Johan Nilsson*
 
* [Science data generation repository](https://github.com/ICESAT-2HackWeek/ScienceDataGeneration)

### 12. Machine learning
*Yara Mohajerani and Shane Grigsby*
 
* [Machine learning repository](https://github.com/ICESAT-2HackWeek/Machine-Learning)


## Citation and License
Tutorial content is made up of original material modified specifically for the Hackweek by a dedicated team of volunteer instructors. We release these materials with a digital object identifier (DOI) to provide an easy way for both contributors and users to cite their work, given it is not necessarily appropriate for a peer-reviewed journal article publication. If you find these tutorials useful or adapt some of the underlying source code for your research (whether or not you attended the Hackweek), we request that you cite it as:

Anthony Arendt, Jessica Scheick, David Shean, Ellen Buckley, Shane Grigsby, Charley Haley, Lindsey Heagy, Yara Mohajerani, Tom Neumann, Johan Nilsson, Thorsten Markus, Fernando Paolo, Fernando Perez, Alek Petty, Axel Schweiger, Ben Smith, Amy Steiker, Sebastian Alvis, Scott Henderson, Nick Holschuh, Zheng Liu, Tyler Sutterley. (2020). ICESAT-2HackWeek/2020_ICESat-2_Hackweek_Tutorials (Version 1.0). Zenodo. http://doi.org/10.5281/zenodo.3966463.

With such a large number of contributors completing a diverse array of tasks and at varying career stages (and therefore with a varying level of need for publications), author order was challenging to establish. Here, author order was determined by:
1. involvement in planning the Hackweek and in preparing this summary repository for release on Zenodo
2. alphabetical order by last name for tutorial leads and presenters
3. alphabetical order by last name for other planning team members who did not present

Special thanks to Yu-Chan and Jane for all of their technical expertise and help in planning and executing the event!

Please click on the Zenodo badge [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3966463.svg)](https://doi.org/10.5281/zenodo.3966463) for the latest citation information and export options.


The content of this project is licensed under the [Creative Commons Attribution 3.0 Unported license](https://creativecommons.org/licenses/by/3.0/), and the underlying source code used to format and display that content is licensed under the [MIT license](LICENSE.md).
