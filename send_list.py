#!/usr/bin/python

import os
import sys
import argparse


parser = argparse.ArgumentParser("Specify the number of jobs")
parser.add_argument("jobs_number")
parser.add_argument("folder_prefix")
args = parser.parse_args()

jobs_number = args.jobs_number
folder_prefix = args.folder_prefix


if "sigma" in folder_prefix:
    fprefix = "sigma_"
elif "polar" in folder_prefix:
    fprefix = "polar_"
# elif "A" in folder_prefix:
#     fprefix = "A_"
# elif "F" in folder_prefix:
#     fprefix = "F_"
else:
    fprefix = folder_prefix


def FileName(parameter_line):
    pl = parameter_line.strip().split(",")
    pl = [p.strip() for p in pl]
    name = os.path.join("data", fprefix + "B{0}_R{1}_M{2}_L{3}_O{4}".format(pl[1],pl[2],pl[3],pl[4],pl[0]) )
    return name



with open("parameter_list", "r") as f:
    lines = f.readlines()
        
num = 0
for i, l in enumerate(lines):
    if l.strip() == "":
        num = i - num
        break
    elif l.strip()[0] == "#":
        num += 1
parameterList = lines[1:num+1]
# num represents how many parameters.


for parameter in parameterList:
    foldername = FileName(parameter)
    print("--------- send_list.py ------------\nNew Parameter: " + foldername)
    with open("parameter", "w") as f1:
        content = lines[0] + parameter + "".join(lines[num+1:])
        f1.write(content)
    
    os.system("python send.py " + jobs_number + " " + foldername )

print("All parameters have been send.")
sys.exit(0)
