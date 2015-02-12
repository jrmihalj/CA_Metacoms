#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#PBS -q short

R CMD BATCH ${arg1}
