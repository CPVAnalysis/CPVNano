import os
import sys
import ROOT



def analyse():

  inputfile = 'bparknano.root'
  f = ROOT.TFile.Open(inputfile)
  tree = f.Get('Events')

  n_events = tree.GetEntries()
  print 'number of events ({}): {}'.format(inputfile, n_events)

  for entry in tree:
    count_phi = 0
    phi_idx = []
    for icand in range(0, entry.nGenPart):
      if abs(entry.GenPart_pdgId[icand]) == 333:
        mother_idx = entry.GenPart_genPartIdxMother[icand]
        mother_pdgid = entry.GenPart_pdgId[mother_idx]
        if abs(mother_pdgid) == 531:
          count_phi += 1
          phi_idx.append(icand)

    if count_phi != 2:
      print 'WARNING - did not find exactly two phi daughters'
      continue

    print 'number of phis: {}'.format(count_phi) 

    phi1_idx = phi_idx[0]
    phi2_idx = phi_idx[1]

    count_phi1_daughters = 0
    count_phi2_daughters = 0

    for icand in range(0, entry.nGenPart):
      if abs(entry.GenPart_pdgId[icand]) == 321 and entry.GenPart_genPartIdxMother[icand] == phi1_idx:
        count_phi1_daughters += 1

      if abs(entry.GenPart_pdgId[icand]) == 321 and entry.GenPart_genPartIdxMother[icand] == phi2_idx:
        count_phi2_daughters += 1

    print 'number of daughters of phi1: {}'.format(count_phi1_daughters) 
    print 'number of daughters of phi2: {}'.format(count_phi2_daughters) 


if __name__ == '__main__':
  ROOT.gROOT.SetBatch(True)
  analyse()



