## Test environments
* local OS X install, R-release
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Ubuntu Linux 16.04 LTS, R-release, GCC
* Fedora Linux, R-devel, clang, gfortran
* Debian Linux, R-devel, GCC ASAN/UBSAN

## R CMD check results
There were no ERRORs or WARNINGs. 
Notes in each platform:

1. Platform:   Windows Server 2008 R2 SP1, R-devel, 32/64 bit
 *Note 1:
    checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Rui Costa <rui.costa@ebi.ac.uk>'
    New submission
 *Comment to note 1:
    Seems to be an unavoidable comment for first submissions.

> checking examples ...
  ** running examples for arch 'i386' ... OK
  ** running examples for arch 'x64' ... NOTE
  Examples with CPU (user + system) or elapsed time > 5s
                   user system elapsed
  probtrans_ebsurv  5.2      0     5.2

0 errors ✓ | 0 warnings ✓ | 2 notes x


## Downstream dependencies
There are currently no downstream dependencies for this package.
