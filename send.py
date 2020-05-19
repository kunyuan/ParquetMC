#!/usr/bin/python
import random
from datetime import datetime
import os
import sys

##### Modify parameters here  ###############
# Cluster="Rutgers"
# Cluster="PBS"
Cluster = "local"
# Cluster = "condor"
############################################

assert len(sys.argv) == 2, "Number of jobs is needed as a parameter!"
Number = int(sys.argv[1])
print "Creating {0} jobs ...".format(Number)
PIDList = range(Number)


def CreateFolder(path):
    if not os.path.exists(path):
        os.system("mkdir "+path)


rootdir = os.getcwd()
execute = "feyncalc.exe"
random.seed(datetime.now())


homedir = os.path.join(rootdir, "Data")
CreateFolder(homedir)

os.system("cp {0} {1}".format(execute, homedir))
os.system("cp {0} {1}".format("parameter", homedir))

if Cluster != "Rutgers":
    outfilepath = os.path.join(homedir, "outfile")
    CreateFolder(outfilepath)
    jobfilepath = os.path.join(homedir, "jobfile")
    CreateFolder(jobfilepath)
else:
    infilepath = homedir

for pid in PIDList:
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
            fjob.write("#!/bin/sh\n"+"#PBS -N "+jobfile+"\n")
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
sys.exit(0)
