#!/usr/bin/python
""" 
Script for submitting jobs
Stephen Inglis
"""

import random
import numpy as np
import glob
import os
import subprocess as sp

def countFolders(path):
    folders = glob.glob(path + 'r*')
    if folders == []:
        return 0
    fmax = 0
    for F in folders:
        t = int(F.split('r')[1])
        fmax = max(t,fmax)
    return fmax+1

def replace(d):
    text = ''
    for i in d:
        text += " -e's/%s/%s/'" % (i,d[i])
    return text

def main():
    # Dictionary of values for modifying the input template
    input_dict = { '__Lx':16,
    '__Ly':16,
    '__betaLow':1.0,
    '__betaHigh':1.0,
    '__JLow':1.0,
    '__JHigh':1.0,
    '__Eq':10000,
    '__MCS':10000,
    '__bins':100,
    '__seed':12345
    }

    # Folder structure will be of the form
    # L%03d/R%.5f/r%04d/bins.txt
    # L = size, B = beta, r = random realization number
    folder_fmt = 'L%03d/B%.6f/r%04d'
    oldRep = 0
    numRep = 1

    fout = open('runme3','w')
    # Allows for changing from normal and modified simulation
    normal = True

    Tc = 1.0/np.log(1+np.sqrt(2))
    Bc = 1.0/Tc
    Lset = [8*i for i in range(1,2)]
    Bset = [1.0/(1.5)**i for i in range(1,101)]
    Bset += [0.02*Bc*i for i in range(1,101)]
    Bset = set(Bset)
    Bset = list(Bset)
    Bset.sort()
    repset = range(oldRep,numRep+oldRep)
    job_min = 0
    job_max = 1
    c_job = 0
    for iL in Lset:
        input_dict['__Lx'] = iL
        input_dict['__Ly'] = iL
        for iB in Bset:
            input_dict['__betaLow'] = iB
            if normal:
                input_dict['__betaHigh'] = iB
            else:
                input_dict['__betaHigh'] = 2*iB
            for ir in repset:
                input_dict['__seed'] = random.randint(1,999999)
                dir_name = folder_fmt % (iL,iB,ir)
                runtime = 0.
                runtime = int(runtime) + 1
                if (c_job >= job_min) and (c_job < job_max):
                    try:
                        os.makedirs(dir_name)
                    except OSError:
                        pass
                    sp.call('sed' + replace(input_dict) + ' <../../../param_template.txt >param.txt',shell=True,cwd=dir_name)
                    sp.call('cp ../../../gradient .',shell=True,cwd=dir_name)
                    sp.call('sqsub --mpp 800M -r %dh -e run.err -o run.out ./gradient' % runtime,shell=True,cwd=dir_name)
                c_job += 1

if __name__ == '__main__':
    main()
