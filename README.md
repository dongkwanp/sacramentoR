sacramentoR
-----------

`sacramentoR` builds out the SAC-SMA-DS model within the R environment.
In addition, `sacramentoR` contains functions and tools to translate
[prmsR](https://projects.cloudwaterlab.com/UMass/prmsR) hydrologic
models to SAC-SMA models as easily and as effectively as possible. This
library is intended to be used by [University of Massachusetts Amherst
Hydrosystems Research Group](http://blogs.umass.edu/hydrosystems/). The
function of this library is to:

1.  Build a distributed Sacramento Soil Moisture Accounting Model
2.  Translate PRMS Models to SAC-SMA

Questions? Comments? Concerns? Contact [Don](mailto:donpark@umass.edu).

Dependencies
------------

`sacramentoR` does not require any libraries to run. However, to use the
translation functions from `prmsR`, the following library is required:

    install.packages(c('sp', 'rgdal'))

Installation Instructions
-------------------------

`sacramentoR` is in constant development. Therefore, using the most
up-to-date version of the library is recommended. This library can be
accessed via HTTPS or SSH. Both options are outlined below.

### Access via HTTPS

Within R, execute the following lines in order to install `sacramentoR`.

    install.packages(c('devtools', 'git2r'))
    devtools::install_git('https://projects.cloudwaterlab.com/UMass/sacramentoR.git')

### Access via SSH

    install.packages(c('git2r', 'devtools'))
    devtools::install_git("git@projects.cloudwaterlab.com:UMass/sacramentoR.git")

How to Use
----------

Call the library using the following command.

    library('sacramentoR')

Troubleshooting
---------------

`sacramentoR` has been tested and is in regular use on the Linux
environment. However, the code should be platform neutral and work with
Windows or MacOS systems as well. If there are any issues or problems,
feel free to open an issue in our [issue
tracker](https://projects.cloudwaterlab.com/UMass/sacramentoR/issues).
