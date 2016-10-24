#!/bin/bash
# Lines starting with #PBS are treated by bash as comments, but interpreted by qsub
# as arguments.  For more details about usage of these arguments see "man qsub"

#
# Name the job.
#SBATCH -J smEtaFreeVp

#
# Set a walltime for the job. The time format is HH:MM:SS

# Run for 30 seconds:
#SBATCH --time=5:00:00

#
# Select one node
#

#SBATCH -N 1
# Select one task per node (similar to one processor per node)
#SBATCH --ntasks-per-node 12

# Set output file name with job number
#SBATCH -o inv_smEtaFreeVp_run.0-%j.out
# Use the janus-debug(1hr maximum)  or use 'himem'
#SBATCH --qos=himem


# Execute the program.
#csh do_run.AZcm.angc.20disc.XYdip.AZmat4.5perc.LT10.csh  pointA.txt 100  AD_test.largerEtaRange
#csh do_run.AZcm.angc.20disc.XYdip.AZmat4.5perc.LT10.smEta.csh pointA.txt 100  AD_test.largerEtaRange
csh do_run.AZcm.angc.20disc.XYdip.AZmat4.5perc.LT10.smEta.csh pointA.txt 100  AD_testFreeVpAni.largerEtaRange

