#!/usr/bin/python
import random
from datetime import datetime
import os
import sys
import argparse
import time


parser = argparse.ArgumentParser("Specify the number of jobs, and the name of working folder.")
parser.add_argument("jobs_number")
parser.add_argument("folder_name")
parser.add_argument("-sc", type=bool, default=False,
        help="If the code need to be self-consistent, the argument should be set as -sc=True.")
args = parser.parse_args()

jobs_number = args.jobs_number
folder_name = args.folder_name
selfConsistent = args.sc

selfConsistent = True

##### Modify parameters here  ###############
# Cluster="Rutgers"
# Cluster="PBS"
Cluster = "local"
# Cluster = "condor"
############################################

Number = int(jobs_number)
print "Creating {0} jobs ...".format(Number)
PIDList = range(Number)


def CreateFolder(path):
    if not os.path.exists(path):
        os.system("mkdir "+path)


rootdir = os.getcwd()
execute = "feyncalc.exe"
random.seed(datetime.now())

homedir = os.path.join(rootdir, folder_name)
CreateFolder(homedir)

os.system("cp {0} {1}".format(execute, homedir))
os.system("cp {0} {1}".format("parameter", homedir))
# os.system("cp {0} {1}".format("green.data", homedir))
# os.system("cp {0} {1}".format("dispersion.data", homedir))

if Cluster != "Rutgers":
    outfilepath = os.path.join(homedir, "outfile")
    CreateFolder(outfilepath)
    jobfilepath = os.path.join(homedir, "jobfile")
    CreateFolder(jobfilepath)
else:
    infilepath = homedir

for pid in PIDList:
    time.sleep(0.3)
    seed = random.randint(0, 2**31-1)
    # print pid, seed
    outfile = os.path.join(outfilepath, "_out{0}".format(pid))  # output files
    jobfile = os.path.join(jobfilepath, "_job{0}.sh".format(pid))  # job files

    if Cluster == "local":
        os.chdir(homedir)
        os.system("./{0} {1} {2} > {3} &".format(execute, pid, seed, outfile))
        os.chdir("..")

    elif Cluster == "condor":
        with open(jobfile, "w") as fjob:
            fjob.write("executable = {0}\n".format(execute))
            fjob.write("arguments = {0} {1}\n".format(pid, seed))
            fjob.write("output ={0}\n".format(outfile))
            fjob.write("initialdir ={0}\n".format(homedir))
            fjob.write("queue")

        os.chdir(homedir)
        os.system("condor_submit {0}".format(jobfile))
        os.system("rm {0}".format(jobfile))
        os.chdir("..")
    elif Cluster == "PBS":
        with open(jobfile, "w") as fjob:
            fjob.write("#!/bin/sh\n"+"#PBS -N " + folder_name+"_"+jobfile+"\n")
            fjob.write("#PBS -o "+homedir+"/Output\n")
            fjob.write("#PBS -e "+homedir+"/Error\n")
            fjob.write("#PBS -l walltime=2000:00:00\n")
            fjob.write("echo $PBS_JOBID >>"+homedir+"/id_job.log\n")
            fjob.write("cd "+homedir+"\n")
            fjob.write(
                "./{0} {1} {2} > {3}".format(execute, pid, seed, outfile))

        os.chdir(homedir)
        os.system("qsub {0}".format(jobfile))
        os.system("rm {0}".format(jobfile))
        os.chdir("..")
    else:
        print("{0} means no submission.".format(Cluster))

print("\nJobs has submitted.")

if selfConsistent:
    os.chdir(rootdir)
    os.system("cp  -r  selfconsistent/  {0}".format(homedir))
    # os.system("./merge_sigma.py  " + folder_name + " > merge_sigma.log &")

sys.exit(0)
