#This is the script that takes a datacard and an injected signal strenght
#this makes the toys, and does the Fit Diagnostics

import glob
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d","--datacard",type=str,help="data card to input")
parser.add_argument("-i","--info",type=str,help="extra info to append to name")
parser.add_argument("-exps","--expectedsignal",type=float,help="expected signal to inject")
parser.add_argument("-t","--ntoys",type=int,help="number of toys to generate")
args = parser.parse_args()

print "------Doing signal injection test------"

#Gather command line options
dcard = args.datacard
info  = args.info
exps  = args.expectedsignal
ntoys = args.ntoys

#Generate deliniations
#This just takes my naming convention to get the signal point out
#Run2_161718_mumu_Zp2000ND800NS200_Zpt100_Hpt300_MET75.txt --> example
splitnz = dcard.split("mumu_")
protoname = splitnz[-1]
name = protoname.split("_Zpt")[0]
print "The name signifier for the output is ",name+"_expectedsignal"+str(exps)+"_"+info
nametocombine = name+"_expectedsignal"+str(exps)+"_"+info

#Generate the toys
print "   ---Generating Injected Toys---   "
print "   injecting signal r = ",exps
print "   generating toys: ",ntoys

#The combine incantation
subprocess.call("combine -M GenerateOnly --saveToys -d "+dcard+" -n "+nametocombine+" --toysFrequentist --bypassFrequentist -t "+str(ntoys)+" --expectSignal "+str(exps),shell=True)


print "     ---Doing Fit Diagnositcs---   "
    
#This is the way we have been doing it with mixed understanding of how the --toyFrequentists works. I think this takes the toys as the prefit model, because --bypassFrequentist only acts when the toys are generated
#subprocess.call("combine -M FitDiagnostics -d "+dcard+" -n "+nametocombine+" -t "+str(ntoys)+" --expectSignal "+str(exps)+" --rMin -100.0 --rMax 100.0 --toysFile=higgsCombine"+nametocombine+".GenerateOnly.mH120.123456.root --toysFrequentist --bypassFrequentist",shell=True)

#test without bypass frequentist, and toys frequentist turned off.
#toys frequentists takes the data_obs as the input model, giving toys directly takes the toys
subprocess.call("combine -M FitDiagnostics -d "+dcard+" -n "+nametocombine+" -t "+str(ntoys)+" --expectSignal "+str(exps)+" --rMin -1000.0 --rMax 1000.0 --toysFile=higgsCombine"+nametocombine+".GenerateOnly.mH120.123456.root",shell=True)

#Different ranges
#subprocess.call("combine -M FitDiagnostics -d "+dcard+" -n "+nametocombine+" -t "+str(ntoys)+" --expectSignal "+str(exps)+" --rMin -50.0 --rMax 50.0 --toysFile=higgsCombine"+nametocombine+".GenerateOnly.mH120.123456.root",shell=True)
#subprocess.call("combine -M FitDiagnostics -d "+dcard+" -n "+nametocombine+" --toysFrequentist --bypassFrequentist -t "+str(ntoys)+" --expectSignal "+str(exps)+" --toysFile higgsCombine"+nametocombine+".GenerateOnly.mH120.123456.root --rMin -100.0 --rMax 50.0",shell=True)
