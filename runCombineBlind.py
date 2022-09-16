import glob
import subprocess

datacards = glob.glob("Run2_161718_ZllHbbMET_datacard*")

for card in datacards:
    #splitnz = card.split("_")
    #name = splitnz[6]
    splitnz = card.split("mumu_")
    protoname = splitnz[-1]
    name = protoname.split(".root")[0]
    print "About to run combine for: ",name
    subprocess.call(["combine","-M","AsymptoticLimits",card,"--noFitAsimov","--run","blind","-n",name,"--rMin","-100.0","--rMax","100.0"])
