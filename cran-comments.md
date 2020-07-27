## Test environments
* Local macOS Catalina install, R-release
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Ubuntu Linux 16.04 LTS, R-release, GCC
* Fedora Linux, R-devel, clang, gfortran
* Debian Linux, R-devel, GCC ASAN/UBSAN

## R CMD check results
There were no ERRORs or WARNINGs for any platform tested.

Notes in each platform:

1. Local macOS Catalina install, R-release: 0 notes.

2. Windows Server 2008 R2 SP1, R-devel, 32/64 bit: 2 notes:   

  * Note 1:  
   
      >checking CRAN incoming feasibility ... NOTE  
      >Maintainer: 'Rui Costa <rui.costa@ebi.ac.uk>'  
      >New submission
    
      Comment to note 1:  
      This seems to be an unavoidable comment for first submissions.  

  * Note 2: 

      > checking examples ...     
      >  ** running examples for arch 'i386' ... OK  
      >  ** running examples for arch 'x64' ... NOTE  
      >  Examples with CPU (user + system) or elapsed time > 5s  
      >                   user system elapsed  
      >  probtrans_ebsurv  5.2      0     5.2  

      Comment to note 2:  
      Negligible overrun of the 5 sec threshold?

3. Ubuntu Linux 16.04 LTS, R-release, GCC: 1 note:   

  * Note:  
   
      >checking CRAN incoming feasibility ... NOTE  
      >Maintainer: 'Rui Costa <rui.costa@ebi.ac.uk>'  
      >New submission
    
      Comment to note:  
      This seems to be an unavoidable comment for first submissions. 
      
4. Fedora Linux, R-devel, clang, gfortran: 2 notes

  * Note 1:  
   
      >checking CRAN incoming feasibility ... NOTE  
      >Maintainer: 'Rui Costa <rui.costa@ebi.ac.uk>'  
      >New submission
    
      Comment to note 1:  
      This seems to be an unavoidable comment for first submissions. 

  * Note 2: 

      > checking examples ...     
      >  Examples with CPU (user + system) or elapsed time > 5s  
      >                   user system elapsed  
      >  probtrans_ebsurv  8.666  0.009   8.676  

      Comment to note 2:  
      Exceeding the 5 sec threshold by 3.6 secs on a single example - can it be considered a benign flaw?
      
5. Debian Linux, R-devel, GCC ASAN/UBSAN: 0 notes.


## Downstream dependencies
There are currently no downstream dependencies for this package.
