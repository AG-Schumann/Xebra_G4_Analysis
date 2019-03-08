#!/usr/bin/python

import sys
import os
import glob
import string
import time
from subprocess import *

#  Write your Username  #
username = "ybiondi"
queue="15:00:00"
#  Check input command line arguments  #

arguments_list = sys.argv[1:]

if len(arguments_list)==2:
	print 'Starting the processing of the MC events for DARWIN'
else:
	sys.exit('The number of command-line arguments does not match the necessary number to execute this program. Exiting program')

#  Assign homedir and binary file from command-line arguments  #

homedir =  arguments_list[0]
print 'This is your home folder for Simulations:', homedir 
binary = arguments_list[1]
print 'This is the current binary that you are using to process events', binary

#  Show the list of Isotopes to be analyzed  #

with open('isotopes_file') as f:
	isotopes_list = f.read().splitlines()
if len(isotopes_list)<1:
	print 'You need to set the isotopes/folders to be process!'
	sys.exit()
else:
	print 'The following isotopes/folders will be processed', isotopes_list

#  Setting directories and creating folders for log files ##
logdir=   homedir + "/" +  "logs" 

#  Nodes management in parallel processing  #
maxnodes = 100 # multiple isotopes
joboffset = 0 #To submit more series of jobs (number of nodes already submitted)
counter = 0 #when it reaches maxnodes no more jobs are submitted

#  Process files for different isotopes, allocate nodes  #

for isotope in isotopes_list:

	working_dir =  homedir + "/" + isotope + "/"
	file_list = os.listdir(working_dir)
	counter = 0
	proc_dir = working_dir + "proc"

	if not os.path.exists(proc_dir):
		os.makedirs(proc_dir)

	for filename in file_list:

		dataset = filename.strip('.root')
		job_prename= "DARWIN" + "_" + isotope + "_"+ dataset
		logfile_prename = logdir + "/" + dataset

		if counter + joboffset < maxnodes :
			proc1 = Popen(["squeue"],stdout=PIPE)
			proc2 = Popen(["grep",username], stdin=proc1.stdout, stdout=PIPE)
			proc3 = Popen(["wc","-l"], stdin=proc2.stdout, stdout=PIPE)
			allocated = int(proc3.communicate()[0])
			proc3.stdout.close()
			proc2.stdout.close()
			proc1.stdout.close()
			
			if (maxnodes > allocated):

				logfile = logfile_prename + ".log"
				jobname = job_prename + ".job"

				os.system('sbatch --job-name=%s --time=%s --export=binary_file=\'%s\',work_dir=%s,inputfile=%s ./proc_many_files.slurm' %(jobname,queue,binary,working_dir,dataset))

			print "\nsubmitting.."	
			print "jobname: ", jobname
			print "queue: ", queue
			print "binary: ", binary
			print "dataset: ", dataset
			print "working folder: ", working_dir
			print "Log file name: ", logfile + "\n"

			print "Submitted job: " + str(joboffset+counter+1)
			counter += 1
			time.sleep(2)

          
		else:
			print "\n\nAllocated nodes: " + str(allocated)
			print "You can't submit! Sleeping 2.5 mins"
			time.sleep(150)
