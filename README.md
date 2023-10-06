# CellDifferentiationNetworks
This repository includes the R code that supports the findings of the article "*Scalable inference of cell differentiation networks in gene therapy clonal tracking studies of haematopoiesis*", https://doi.org/10.1093/bioinformatics/btad605 .

The R package used for the analyses is publicy available at https://cran.r-project.org/package=Karen.
This can be installed by the usual R command

```R
install.packages("Karen")
```

The package includes the datasets being analysed in the paper.
Furthermore, the vignette file https://cran.r-project.org/web/packages/Karen/vignettes/Karen.pdf includes executable examples that allow the user to apply our method to their own clonal tracking datasets.

To reproduce the figures and findings of our paper, it is sufficient to navigate to the main folder containing all the *.sh* files, and run all the *.sh* files by using the bash command

```bash
sbatch ./FILE_NAME.sh
```

The *.sh* files can be found in the main folder of the *master* branch. The syntax of the *.sh* file names is explained as follows.

	▪	Karen-C1_x-y-z-n —> runs the comparative study 1 with measurement noise parameters equal to x, sampling frequency equal to y, false-negative rate equal to z, and n distinct clones
	▪	Karen-C2_x-y-z-n —> runs the comparative study 2 with measurement noise parameters equal to x, sampling frequency equal to y, false-negative rate equal to z, and n distinct clones
	▪	Karen-sim-MS-XY —> runs the simulation study with true model X and candidate model Y
	▪	Karen-scaling.sh —> runs the scalability study
	▪	Karen-resources-x_y_z_w —> runs a simulation study with sampling frequency equal to x,  false-negative rate equal to y, z distinct clones, using w distinct processors
	▪	Karen-MS-Y.sh —> runs Karen to fit model Y on clonal tracking data from the genotoxicity study.
	▪	Karen-RM-Y.sh —> runs Karen to fit model Y on clonal tracking data from the rhesus macaque study.
	▪	Karen-CT-X-Y.sh —> runs Karen to fit model Y on clonal tracking data from the gene therapy clinical trial X.

Each *.sh* file runs a specific *.R* file, whose code can be found in the *source* folder of the *master* branch.
