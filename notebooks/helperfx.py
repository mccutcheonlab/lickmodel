# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 08:36:38 2019

@author: James Rig
"""
import numpy as np
import xlrd
import csv


def medfilereader(filename, varsToExtract = 'all',
                  sessionToExtract = 1,
                  verbose = False,
                  remove_var_header = False):
    if varsToExtract == 'all':
        numVarsToExtract = np.arange(0,26)
    else:
        numVarsToExtract = [ord(x)-97 for x in varsToExtract]
    
    f = open(filename, 'r')
    f.seek(0)
    filerows = f.readlines()[8:]
    datarows = [isnumeric(x) for x in filerows]
    matches = [i for i,x in enumerate(datarows) if x == 0.3]
    if sessionToExtract > len(matches):
        print('Session ' + str(sessionToExtract) + ' does not exist.')
    if verbose == True:
        print('There are ' + str(len(matches)) + ' sessions in ' + filename)
        print('Analyzing session ' + str(sessionToExtract))
    
    varstart = matches[sessionToExtract - 1]
    medvars = [[] for n in range(26)]
    
    k = int(varstart + 27)
    for i in range(26):
        medvarsN = int(datarows[varstart + i + 1])
        
        medvars[i] = datarows[k:k + int(medvarsN)]
        k = k + medvarsN
        
    if remove_var_header == True:
        varsToReturn = [medvars[i][1:] for i in numVarsToExtract]
    else:
        varsToReturn = [medvars[i] for i in numVarsToExtract]

    if np.shape(varsToReturn)[0] == 1:
        varsToReturn = varsToReturn[0]
    return varsToReturn

def metafilemaker(xlfile, metafilename, sheetname='metafile', fileformat='csv'):
    with xlrd.open_workbook(xlfile) as wb:
        sh = wb.sheet_by_name(sheetname)  # or wb.sheet_by_name('name_of_the_sheet_here')
        
        if fileformat == 'csv':
            with open(metafilename+'.csv', 'w', newline="") as f:
                c = csv.writer(f)
                for r in range(sh.nrows):
                    c.writerow(sh.row_values(r))
        if fileformat == 'txt':
            with open(metafilename+'.txt', 'w', newline="") as f:
                c = csv.writer(f, delimiter="\t")
                for r in range(sh.nrows):
                    c.writerow(sh.row_values(r))
    
def metafilereader(filename):
    
    f = open(filename, 'r')
    f.seek(0)
    header = f.readlines()[0]
    f.seek(0)
    filerows = f.readlines()[1:]
    
    tablerows = []
    
    for i in filerows:
        tablerows.append(i.split('\t'))
        
    header = header.split('\t')
    # need to find a way to strip end of line \n from last column - work-around is to add extra dummy column at end of metafile
    return tablerows, header

def isnumeric(s):
    try:
        x = float(s)
        return x
    except ValueError:
        return float('nan')

def extractlicks(self, substance):
    licks = medfilereader(self.medfile,
                              varsToExtract = sub2var(self, substance),
                                                remove_var_header = True)
    lickData = lickCalc(licks, burstThreshold=0.5, binsize=120)        

    return lickData