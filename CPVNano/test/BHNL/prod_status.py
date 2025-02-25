import os
import sys
import glob
from os import path
import ROOT

user = 'anlyon'
version_label = 'V42_Bc'
tree_name = 'signal_tree'
#tag = '06Feb23'
#tag = '06Feb23_partial_v3'
#tag = '06Feb23_partial'
#tag = '06Feb23_norm'
tag = '06Feb23_15Jun23'

def getMass(point):
  mass_str = point[point.rfind('mass')+4:point.find('ctau')-1]
  mass = mass_str.replace('p', '.')
  return float(mass)

indirectory = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}'.format(user, version_label)

# get all the subdirectories (signal points)
pointdirs = [f for f in glob.glob('{}/*'.format(indirectory))]
for pointdir in pointdirs:
  print pointdir

# sort them by mass
pointdirs.sort(key=getMass)

#for pointdir in pointdirs:
#  point = pointdir[pointdir.rfind('mass'):]

count = 0
percentage_tot = 0
for pointdir in pointdirs:
  count = count + 1
  point = pointdir[pointdir.rfind('mass'):]
  #if 'mass1p95_ctau10p0' not in point: continue
  flat_filename = pointdir + '/nanoFiles/merged/flat_bparknano_{}.root'.format(tag)

  # check the file
  if not path.exists(flat_filename):
    continue
  rootfile = ROOT.TNetXNGFile.Open(flat_filename, 'r')
  if not rootfile: continue
  if not rootfile.GetListOfKeys().Contains('signal_tree'): continue

  flat_file = ROOT.TFile.Open(flat_filename, 'READ')
  tree = flat_file.Get(tree_name)

  hist_name = 'hist_{}'.format(point)
  hist = ROOT.TH1D(hist_name, hist_name, 100, 0, 10)
  tree.Project(hist_name, 'hnl_mass', 'ismatched==1')
  n_events = int(hist.Integral())

  if n_events < 3000: colour = 'red'
  if n_events >= 3000 and n_events < 8500: colour = 'blue'
  if n_events >= 8501: colour = 'darkgreen'

  percentage_perpoint = round(n_events / 8500. * 100, 1)
  if percentage_perpoint > 100: percentage_perpoint = 100
  percentage_tot += percentage_perpoint

  #print '{} & {} & {}\% \\\ '.format(point.replace('_', '\_'), '\color{' + colour + '}{' + str(n_events) + '}', percentage_perpoint)
  print '{} & {} \\\ '.format(point.replace('_', '\_'), '\color{' + colour + '}{' + str(n_events) + '}')

print 'total percentage: {}\%'.format(round(percentage_tot / float(count), 1))





