SHELL=/bin/bash
current_date := $(shell date +'%Y-%m-%d_%H-%M')
SINGULARITY=singularity exec /beegfs/scratch/bruening_scratch/pklemm/singularity/singularity-images/mytidyverse-4.3.1-2.simg

.PHONY: plots
plots:
	-mkdir -p docs/plots
	$(SINGULARITY) Rscript -e 'rmarkdown::render("analysis/pnoc/plots.Rmd", params = list("dev" = FALSE))'
	cp analysis/pnoc/plots.html docs/analysis/${current_date}_custom.html
