'''
Script for the computation of the expected D-D* correlation function
'''

import os
import sys
import argparse
import pandas as pd
import numpy as np
import yaml
from ROOT import TFile, TGaxis, TCanvas, TF1, TGraph, TH1F, TSpline3, TLegend, TLatex, TLine # pylint: disable=import-error,no-name-in-module
from ROOT import kAzure, kGreen, kOrange, kRed, kBlack, kGray, kFullCircle, kOpenCircle, kFullDiamond, gStyle, kRainBow, gROOT # pylint: disable=import-error,no-name-in-module


def WeightedAverage(graph, weights):
    '''
    Compute the weighted average of a graph with an histogram (TH1)
    '''
    smeared = 0
    counts = weights.Integral(1, weights.GetNbinsX())
    for iBin in range(weights.GetNbinsX()):
        freq = weights.GetBinContent(iBin + 1) / counts
        y = graph.Eval(weights.GetBinCenter(iBin+1))
        smeared += freq * y
    return smeared


def SmearGraph(graph, matrix):
    '''
    Smear a graph with a smearing matrix that has:
     x axis: true variable
     y axis: reconstructed variable.
    '''
    gSmeared = TGraph(1)

    iPoint = 0
    for iBin in range(matrix.GetNbinsX()):
        hProj = matrix.ProjectionY(f'hProj_{iBin+1}', iBin+1, iBin+1)
        counts = hProj.Integral(1, hProj.GetNbinsX())

        if counts < 1:
            continue

        x = matrix.GetXaxis().GetBinCenter(iBin+1)
        ySmear = WeightedAverage(graph, hProj)
        gSmeared.SetPoint(iPoint, x, ySmear)
        iPoint += 1

    return gSmeared


gStyle.SetPadBottomMargin(0.12)
gStyle.SetPadTopMargin(0.05)
gStyle.SetPadLeftMargin(0.12)
gStyle.SetPadRightMargin(0.035)
gStyle.SetTitleOffset(1.2, 'yz')
gStyle.SetTitleSize(0.045, 'xy')
gStyle.SetLabelSize(0.045, 'xy')
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPalette(kRainBow)
TGaxis.SetMaxDigits(3)

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('cfgFileName', metavar='text', default='config.yml')
parser.add_argument('--yMax', type=float, default=8.)
parser.add_argument('--xMax', type=float, default=0.5)
parser.add_argument('--source', type=int, default=1)
parser.add_argument('-b', action='store_true', default=False)
args = parser.parse_args()

gROOT.SetBatch(args.b)

with open(args.cfgFileName, 'r') as ymlcfgFileName:
    cfg = yaml.load(ymlcfgFileName, yaml.FullLoader)

Dspecie = cfg['species']['D']
Dstarspecie = cfg['species']['Dstar']

if Dspecie == 411:
    Dtitle = 'D^{#plus}'
    brD = 0.0938
elif Dspecie == -411:
    Dtitle = 'D^{#minus}'
    brD = 0.0938
elif Dspecie == 421:
    Dtitle = 'D^{0}'
    brD = 0.03951
elif Dspecie == -421:
    Dtitle = '#bar{D}^{0}'
    brD = 0.03951
else:
    print(f'ERROR: D specie {Dspecie} not supported! Exit')
    sys.exit()

if Dstarspecie == 413:
    Dstartitle = 'D*^{#plus}'
    brDstar = 0.677 * 0.03951
elif Dstarspecie == -413:
    Dstartitle = 'D*^{#minus}'
    brDstar = 0.677 * 0.03951
elif Dstarspecie == 423:
    Dstartitle = 'D*^{0}'
    brDstar = 0.353 * 0.03951
elif Dstarspecie == -423:
    Dstartitle = '#bar{D}*^{0}'
    brDstar = 0.353 * 0.03951
else:
    print(f'ERROR: D* specie {Dstarspecie} not supported! Exit')
    sys.exit()


inFiles = {}
hResoSE = {}
for magField in ['5kG', '10kG', '20kG']:
    inFiles[magField] = TFile.Open(cfg['smearfiles'][magField])

    hResoSE[magField] = inFiles[magField].Get(f'hResoSE_{Dstarspecie}_{Dspecie}')
    hResoSE[magField].SetDirectory(0)
    hResoSE[magField].RebinX(cfg['predictions']['smear']['rebin'])
    hResoSE[magField].RebinY(cfg['predictions']['smear']['rebin'])
    inFiles[magField].Close()

    cSmearingMatrix = TCanvas('cSmearingMatrix', '', 600, 600)
    hResoSE[magField].Draw('colz')
    cSmearingMatrix.SaveAs(magField + cfg['output']['file'])

predictions = {'nosmear': pd.read_csv(cfg['predictions'][f'{args.source}fm'], names=['kstar', 'cf', 'idk'], sep=' ')}

predColors = {'nosmear': kAzure+4,
              '5kG': kGreen+2,
              '10kG': kOrange+7,
              '20kG': kRed+1}

gPred = {}
gPred['nosmear'] = TGraph(1)
for iP, (kStar, cf) in enumerate(
    zip(predictions['nosmear']['kstar'].to_numpy(), predictions['nosmear']['cf'].to_numpy())):
    gPred['nosmear'].SetPoint(iP, kStar/1000, cf)
gPred['nosmear'].SetLineColor(predColors['nosmear'])
gPred['nosmear'].SetLineWidth(2)

gPred['5kG'] = SmearGraph(gPred['nosmear'], hResoSE['5kG'])
gPred['10kG'] = SmearGraph(gPred['nosmear'], hResoSE['10kG'])
gPred['20kG'] = SmearGraph(gPred['nosmear'], hResoSE['20kG'])

gPred['5kG'].SetLineColor(predColors['5kG'])
gPred['10kG'].SetLineColor(predColors['10kG'])
gPred['20kG'].SetLineColor(predColors['20kG'])

gPred['5kG'].SetMarkerColor(predColors['5kG'])
gPred['10kG'].SetMarkerColor(predColors['10kG'])
gPred['20kG'].SetMarkerColor(predColors['20kG'])

gPred['5kG'].SetLineWidth(2)
gPred['10kG'].SetLineWidth(2)
gPred['20kG'].SetLineWidth(2)

cPred = TCanvas('cPred', '', 600, 600)
cPred.DrawFrame(0, 0, 0.3, 10, ';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
for key, pred in gPred.items():
    pred.Draw('lc same')
leg = TLegend(0.3, 0.55, 0.9, 0.9, 'Expected CF')
leg.AddEntry(gPred['nosmear'], 'No smearing')
leg.AddEntry(gPred['20kG'], '#it{p}_{T} resolution for B = 2 T')
leg.AddEntry(gPred['10kG'], '#it{p}_{T} resolution for B = 1 T')
leg.AddEntry(gPred['5kG'], '#it{p}_{T} resolution for B = 0.5 T')
leg.Draw('same')
cPred.SaveAs('cSmearedPred.pdf')

for magField in cfg['smearfiles']:
    print(magField)
    oFile = TFile(f'predictions/Tcc/{args.source}fm_smeared_{magField}.root', 'recreate')
    gPred[magField].SetTitle(';#it{k}* (GeV/#it{c});#it{C}(#it{k}*)')
    gPred[magField].Write('gCorrelationFunction')
    oFile.Close()

input('Press enter to exit')
