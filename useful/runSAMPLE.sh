#!/bin/bash
# Lines starting with #PBS are treated by bash as comments, but interpreted by qsub
# as arguments.  For more details about usage of these arguments see "man qsub"

#
# Name the job.
#SBATCH -J touch

#
# Set a walltime for the job. The time format is HH:MM:SS

# Run for 30 seconds:
#SBATCH --time=12:00:00

#
# Select one node
#

#SBATCH -N 1
# Select one task per node (similar to one processor per node)
#SBATCH --ntasks-per-node 12

# Set output file name with job number
#SBATCH -o testjob-%j.out
# Use the janus-debug(1hr maximum)  or use 'himem'
#SBATCH --qos=janus


# Execute the program.
lfs find /lustre/janus_scratch/jixi7887/ -type f -ctime +60 -print0 | xargs -0 touch
