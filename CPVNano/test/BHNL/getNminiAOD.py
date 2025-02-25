import os
import sys 

filename = '/BToHNLEMuX_HNLToEMuPi_SoftQCD_b_mHNL4p5_ctau0p1mm_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-Custom_RDStar_BParking_102X_upgrade2018_realistic_v15-v2/MINIAODSIM'

os.system('rm tmp.txt')

command = 'dasgoclient --query="file dataset={ds} | grep file.nevents" > tmp.txt'.format(ds=filename)
os.system(command)

f = open('tmp.txt')
lines = f.readlines()
n_events = 0
for line in lines:
  n_events = n_events + float(line)

print 'filename: {}'.format(filename)
print 'n_events: {}'.format(n_events)

  


