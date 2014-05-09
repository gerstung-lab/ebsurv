all: roxygen build

roxygen: CoxHD
	R -e 'library(roxygen2); roxygenize("CoxHD", unlink.target=TRUE)'
 
check: CoxHD
	R CMD check CoxHD

build: CoxHD
	R CMD build --no-vignettes CoxHD

install: CoxHD
	R CMD install CoxHD