#this is just scripting
#This script find the expected limits, then executes another script that makes the toys, does the limits

import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d","--datacard",type=str,help="data card to input")
parser.add_argument("-i","--info",type=str,help="extra info to append to name")
parser.add_argument("-t","--ntoys",type=int,help="number of toys to generate")
args = parser.parse_args()

print "------Running Asymptotic Limits to Get the r values------"

dcard = args.datacard
info  = args.info
ntoys = args.ntoys

#Generate deliniations
#This just takes my naming convention to get the signal point out
#Run2_161718_mumu_Zp2000ND800NS200_Zpt100_Hpt300_MET75.txt --> example
splitnz = dcard.split("mumu_")
protoname = splitnz[-1]
name = protoname.split("_Zpt")[0]
print "The name signifier for the output is ",name+"_"+info


#Run the limits
print "Running the limit, might print a warning"
asympout = subprocess.check_output(["combine","-M","AsymptoticLimits",dcard,"--noFitAsimov","--run","blind","-n",name])

#Gather the limits for injection
limout = asympout.split(" -- AsymptoticLimits ( CLs ) --")[-1]
lines = limout.split('\n')[1:6]#just the lines printed with the limit
limits = {}

for line in lines:
    key = line.split(":")[0]
    limit = float(line.split(" r < ")[-1])
    limits[key] = limit#Expected quantile is the key, limit stored as float

print(limits)

limits["bkgonly"] = 0.0

rvals = ["Expected 50.0%","Expected 84.0%","Expected 16.0%","bkgonly"]

#Calls Another script for doing the limits
for r in rvals:
    subprocess.call(["python","generateSignalInjectionFiles.py","-d",dcard,"-i",info,"-exps",str(limits[r]),"-t",str(ntoys)])


