# SigMT: An open-source python package for magnetotelluric data processing

**NOTE:** A major update to version 2, which includes a graphical user interface (GUI), will be released in a few days (in September 2024). Please stay tuned for updates to this repository.
##
SigMT is a python package designed for the processing of the raw magnetotelluric (MT) data to obtain the MT impedance and tipper estimates. It works in an automated way, so that manual time series inspection and editing are not required. Mahalanobis based data selection tool is implemented in the package to avoid the manual editing of time series. The final impedance estimation is done using the robust estimation method. Different data selection tools such as coherency threshold, polarization direction are included in this package.

## How to cite
Ajithabh, K.S., Patro, P.K., 2023. SigMT: An open-source Python package for magnetotelluric data processing. Computers & Geosciences, 171, 105270. https://doi.org/10.1016/j.cageo.2022.105270

## How to use test data

Test data is given in `data` folder in the repository. Give path to `data` folder in the `main_script.py` file for the testing purpose. 
For example give, `project_path = 'C:/Users/Username/Documents/GitHub/SigMT/data/'`

Complete site data can be downloaded using the link https://devapps.ngri.res.in/patro/TEST.zip 

## Downloads
* Download sample MT data for testing: https://devapps.ngri.res.in/patro/TEST.zip
* Download user manual: https://github.com/ajithabhks/SigMT/blob/main/docs/User_Manual_for_SigMT.pdf

## Contact details
* K. S. Ajithabh

  `Email:` ajithabhks@gmail.com

* Prasanta K. Patro

  `Email:` patrobpk@ngri.res.in


