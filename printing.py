# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 14:46:01 2017

@author: Aryn Harmon
This is a module which contains printFigures(); this function plots a comparison of the FDM and the ROM (or any two models). It should plot side by side, titles, axis titles, and in the correct typeface. It defaults to outputting a single svg file, but 'png' can be specified instead. 
"""
import matplotlib.pyplot as plt
from matplotlib import rc

def printFigures( positionVector,tstep,model1,model1Title,model2,model2Title,nLines,filetype='svg' ):
  rc('font',**{'family':'serif','serif':['CMU Serif']})
  rc('text', usetex=True)
  plt.rcParams['font.size'] = 24
  plt.figure(figsize=(15, 7))
  plt.subplot( 1,2,1 )
  plt.title( model1Title )
  plt.xlabel( 'Position' )
  plt.ylabel( 'Amplitude' )
  plt.plot( positionVector, model1[:,0] )
  for line in range(1,nLines+1):
    plt.plot( positionVector, model1[:,line*int(tstep/nLines)-1] )
  plt.subplot( 1,2,2 )
  plt.title( model2Title )
  plt.xlabel( 'Position' )
  plt.ylabel( 'Amplitude' )
  plt.plot( positionVector, model2[:,0] )
  for line in range(1,nLines+1):
    plt.plot( positionVector, model2[:,line*int(tstep/nLines)-1] )
  plt.tight_layout()
  if filetype == 'svg':
    plt.savefig('fig.svg')
  if filetype == 'png':
    plt.savefig('fig.png',dpi=300)
    plt.show()