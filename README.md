sacramentoR
-----------

`sacramentoR` builds out the SAC-SMA distributed model within the R environment.
In addition, `sacramentoR` contains functions and tools to translate
[prmsR](https://github.com/dongkwanp/prmsR) hydrologic
models to SAC-SMA models as easily and as effectively as possible. This
library is intended to be used by [University of Massachusetts Amherst
Hydrosystems Research Group](http://blogs.umass.edu/hydrosystems/). The
function of this library is to:

1.  Build a distributed Sacramento Soil Moisture Accounting Model
2.  Translate PRMS Models to SAC-SMA

Certain parts of the code is inspired from [Dr. Sungwook Wi](https://scholar.google.com/citations?user=G6H6U8AAAAAJ&hl=en)'s Sacramento model in MATLAB.

Questions? Comments? Concerns? Submit an issue!

Dependencies
------------

`sacramentoR` requires the following libraries to run.

    install.packages(c('tidyverse', 'xts', 'zoo'))

In addition, the following packages are recommended to use the
translation functions from `prmsR`:

    install.packages(c('sp', 'rgdal'))

Installation Instructions
-------------------------

`sacramentoR` is in constant development. Therefore, using the most
up-to-date version of the library is recommended. This library can be
accessed via HTTPS or SSH. Both options are outlined below.

### Access via HTTPS

Within R, execute the following lines in order to install `sacramentoR`.

    install.packages(c('devtools', 'git2r'))
    devtools::install_git('https://github.com/dongkwanp/sacramentoR.git')

### Access via SSH

    install.packages(c('git2r', 'devtools'))
    devtools::install_git("git@github.com:dongkwanp/sacramentoR.git")

How to Use
----------

Call the library using the following command.

    library('sacramentoR')

Available Models List
---------------------

`sacramentoR` is still very early in development. However, the following
modules are available for use:

**Daylight Models**

1.  N method
2.  CBM Model

**Potential Evapotranspiration Models**

1.  Hamon Method
2.  Hargreaves Method (Planned)

**Snow Models**

1.  SNOW17
2.  Isnobal (Planned)

**Hydrology Models**

1.  Sacramento Soil Moisture Accounting Model
2.  Conversion from PRMS (via
    [prmsR](https://github.com/dongkwanp/prmsR))

**Routing Models**

1.  Lohmann Routing Model

**Utility Functions**

1.  prmsR to sacramentoR translation function
2.  A\_v calculation function
3.  Area Ratio calculation function
4.  Leap Year calculation function

Troubleshooting
---------------

`sacramentoR` has been tested and is in regular use on the Linux
environment. However, the code should be platform neutral and work with
Windows or MacOS systems as well. If there are any issues or problems,
feel free to open an issue in our [issue
tracker](https://github.com/dongkwanp/sacramentoR/issues).
