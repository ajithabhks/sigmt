# SigMT: An open-source python package for magnetotelluric data processing

SigMT is a python package designed for the processing of the raw magnetotelluric (MT) data to obtain the MT impedance and tipper estimates. It works in an automated way, so that manual time series inspection and editing are not required. Mahalanobis based data selection tool is implemented in the package to avoid the manual editing of time series. The final impedance estimation is done using the robust estimation method. Different data selection tools such as coherency threshold, polarization direction are included in this package.

## How to install
Please note that SigMT currently supports only Metronix data format (.ats).

Open anaconda prompt and type:

```
pip install sigmt
```

After installation, type:

```
sigmt
```

## How to cite
If you use SigMT for publication, please cite the following paper:

Ajithabh, K.S., Patro, P.K., 2023. SigMT: An open-source Python package for magnetotelluric data processing. Computers & Geosciences, 171, 105270. https://doi.org/10.1016/j.cageo.2022.105270

## Discussions
If you have any questions, feedback, or suggestions regarding SigMT, please feel free to post them in the [discussions section](https://github.com/ajithabhks/sigmt/discussions).

## Downloads
* Read user guide: https://github.com/ajithabhks/sigmt/blob/main/docs/user_guide_sigmt.ipynb
* Metronix test data link: https://cloud.geo-metronix.de/s/GcigJA3Zp8zTAif

  Download Northern_Mining.zip, use Sarıçam site in the ts folder.

## Contact details
* K. S. Ajithabh

  `Email:` ajithabhks@gmail.com


