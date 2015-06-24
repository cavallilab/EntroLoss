# entroloss.py - Quantification of entropy-loss in replica-averaged modelling.
# Copyright (C) 2015 SIMON OLSSON and ANDREA CAVALLI
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#
try: 
  from numpy import *
except:
  print "numpy not found, please install and try again"
  quit()

from sys import argv

def sum_probs(f_samples, N, n_samples=5000):
 """
 computes sum probabilities
 """
 num_f_samples, dim_f_samples = shape(f_samples)
 samples = []
 for s in range(n_samples):
   indeces = random.randint(0,num_f_samples-1, N)
   samples.append(f_samples[indeces].sum(axis=0))
 return array(samples) 

def hist_prod(hista, histb):
  """
  returns product of input histograms
  """
  if hista[1] is not histb[1]:
    print "Histograms have different number of bins"
    quit()
  else:
    return (hista[0]*histb[0],hista[1])

def hist_norm(hist):
  """
  returns normalized histogram
  """
  h0s = float(sum(hist[0]))
  if h0s==0:
    return (ones(hist[0].shape)/float(sum(ones(hist[0].shape))*abs(hist[1][1]-hist[1][0])),hist[1])
  else:
    return (hist[0]/float(sum(hist[0])*abs(hist[1][1]-hist[1][0])),hist[1])

def hist_norm_const(hist):
  """
  returns normalization constant of input histogram
  """
  return float(sum(hist[0])*abs(hist[1][1]-hist[1][0]))

def hist_mean(hist):
  return average([(hist[1][i+1]+hist[1][i])/2. for i in range(len(hist[1])-1)], weights=hist[0])

def estimate_missing_bins(hist, nsamples, eval_at_edge=False):
  """ 
  Infers missing bins using a first order polynomial fit to the log of the observed bins
  hist: input sparse histogram    
  """
  nzero_counts = where(hist[0]>0.)
  x_vals = array([(hist[1][i+1]+hist[1][i])/2. for i in range(len(hist[1])-1)])[nzero_counts]
  linfit = polyfit(x_vals, log(hist[0][nzero_counts]/nsamples), 1)
  if eval_at_edge:
    #evaluates extrapolation at lower bin edge
    return (exp(polyval(linfit,[hist[1][i] for i in range(len(hist[1])-1)])),hist[1])
  else:
    #evaluates extrapolation at bin-center
    return (exp(polyval(linfit,[(hist[1][i+1]+hist[1][i])/2. for i in range(len(hist[1])-1)])),hist[1])

def estimate_missing_bins_fancy(samples, hist, eval_at_edge=False):
  """ 
  Infers missing bins using a gaussian fit to original samples 
  hist: input sparse histogram    
  """
  fancy_hist = histogram(samples, density=True, bins=128)
  nzero_counts = where(fancy_hist[0]>0.)
  x_vals = array([(fancy_hist[1][i+1]+fancy_hist[1][i])/2. for i in range(len(fancy_hist[1])-1)])[nzero_counts]
  linfit = polyfit(x_vals, log(fancy_hist[0][nzero_counts]), 3)
  h2=(exp(polyval(linfit,[(fancy_hist[1][i+1]+fancy_hist[1][i])/2. for i in range(len(fancy_hist[1])-1)])),fancy_hist[1])
  h2c = hist_norm_const(h2)
  if eval_at_edge:
    #evaluates extrapolation at lower bin edge
    return (exp(polyval(linfit,[hist[1][i] for i in range(len(hist[1])-1)])),hist[1])
  else:
    #evaluates extrapolation at bin-center
    return (exp(polyval(linfit,[(hist[1][i+1]+hist[1][i])/2. for i in range(len(hist[1])-1)]))/h2c,hist[1])


def hist_kl(hista,histb):
  """
  evaluates KL-divergence
  """
  if hista[1] is not histb[1]:
    print "Histograms have different number of bins"
    quit()
  else:
    return sum((hista[0]+1e-6)*(log(hista[0]+1e-6)-log(histb[0]+1e-6)))
 

if __name__=="__main__":
  hist_bins = 24
  random.seed(303808)
  if len(argv) == 1:
    print "input please"
    quit()
  

  distances = loadtxt(argv[1]).T
  noes = loadtxt(argv[2]).tolist()
  Ns = map(int, argv[3:])

  num_dists,num_noes = shape(distances)
  distances = distances 
  madist = distances.max(axis=0)
  minimum_dists = distances.min(axis=0)
  midist = array([mid for mid in minimum_dists])
  dbins = {}
  dhist = {}
  
  for j,noe in enumerate(noes):
    dbins[j] = linspace(midist[j],madist[j], hist_bins)
    dhist[j] = hist_norm(histogram(distances[:,j], range=(midist[j], madist[j]), bins=dbins[j]))
 
    #print "histmean for distances", hist_mean(dhist[j]), " and samplemean:", distances[:,j].mean()
    if isnan(hist_mean(dhist[j])):
      print "NaN histogram mean of index",j
      print "exiting..."
      quit()
  
  for N in Ns: 
    sprobs = sum_probs(distances, N, 64*num_dists)
    kls=[]
    for j,noe in enumerate(noes):
        if noe<minimum_dists[j]:
          continue
        # normalization factor
        ssamples = (N+1)*noe-sprobs[:,j]
        shists = histogram(ssamples,range=(midist[j], madist[j]), bins=dbins[j])
        if(mean(ssamples)<midist[j] or mean(ssamples)>madist[j]):
          print "#trying to extrapolate bins, not using this datapoint to evaluate entropy loss"
          shists = estimate_missing_bins_fancy(ssamples, shists)
          lpf ="#"
        else:
          shists = hist_norm(shists)
          lpf =""
        nprodhist=hist_norm(hist_prod(shists, dhist[j]))
        if lpf=="":
          kls.append(hist_kl(dhist[j], nprodhist))
    print N+1, mean(kls)
