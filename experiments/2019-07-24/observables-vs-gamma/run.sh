#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --account=def-arutenbe
#SBATCH --output=slurmoutput/run_%A_%a.out


module restore standard_modules

python observables-vs-gamma.py $1 $2