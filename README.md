# CellDifferentiationNetworks
This repository includes the R code that supports the results of the article "*Scalable inference of cell differentiation networks in gene therapy clonal tracking studies of haematopoiesis*", https://doi.org/10.1093/bioinformatics/btad605 .

The R code can be found in the *source* folder of the *master* branch. To reproduce the results reported in the paper it is sufficient to run all the *.sh* files of the *master* branch.The syntax of the sh file names is explained as follows.

	▪	Karen-C1_x-y-z-n —> runs the comparative study 1 with measurement noise parameters equal to x, sampling frequency equal to y, false-negative rate equal to z, and n distinct clones
	▪	Karen-C2_x-y-z-n —> runs the comparative study 2 with measurement noise parameters equal to x, sampling frequency equal to y, false-negative rate equal to z, and n distinct clones
	▪	Karen-sim-MS-XY —> runs the simulation study with true model X and candidate model Y
	▪	Karen-scaling.sh —> runs the scalability study
	▪	Karen-resources-x_y_z_w —> runs a simulation study with sampling frequency equal to x,  false-negative rate equal to y, z distinct clones, using w distinct processors
	▪	Karen-MS-Y.sh —> runs Karen to fit model Y on clonal tracking data from the genotoxicity study.
	▪	Karen-RM-Y.sh —> runs Karen to fit model Y on clonal tracking data from the rhesus macaque study.
	▪	Karen-CT-X-Y.sh —> runs Karen to fit model Y on clonal tracking data from the gene therapy clinical trial X.

The R package used for the analyses is publicy available at https://cran.r-project.org/package=Karen.
