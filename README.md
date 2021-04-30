# mRNAmarkup - a workflow for annotating transcript sequences


## Installation

mRNAmarkup is implemented as a bash script which should work on any Linux or
UNIX system.  Please see [INSTALL](./INSTALL) for requirements of other software
and specific set-up instructions.


## Quick Start

If you want to avoid any installation hussles, the following will do nicely,
using our mRNAmarkup [Singularity](https://www.sylabs.io/docs/) container that
encapsulates all scripts, programs, and system packages.


```
git clone https://github.com/BrendelGroup/mRNAmarkup
cd mRNAmarkup/
singularity pull http://BrendelGroup.org/SingularityHub/mRNAmarkup.sif
./xsetup
cd db
singularity exec -e -B ${PWD}/.. ../mRNAmarkup.sif bash 0README
cd ..
cd data
singularity exec -e -B ${PWD}/.. ../mRNAmarkup.sif ./xtest
xdiff
```


## Reference

manuscript to be submitted


## Contact

Please direct all comments and suggestions to
[Volker Brendel](<mailto:vbrendel@indiana.edu>)
at [Indiana University](http://brendelgroup.org/).
