#!/bin/sh
#$ -S /bin/sh
#$ -cwd

R CMD BATCH ${arg1}
