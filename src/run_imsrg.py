#!/usr/local/bin/python3
import os,sys
import subprocess
HOME=os.path.expanduser('../testwork/')

noshell = [] # no-valence-shell
pshell = ['0p1','0p3'] # p-shell
sdshell = ['0d3','0d5','1s1'] # sd-shell
pfshell = ['1p1','1p3','0f7','0f5'] # pf-shell
sdgshell = ['0g9','0g7','1d5','1d3','2s1'] # sdg-shell

exe = './imsrg++'
args = {}
hwlist = [10,15,20,25,30,35,40,45,50]
elist = [2,4,6,8]

vplist = noshell
vnlist = vplist

def SetDefaultArgs():
    hw = 20
    args['hw'] = hw
    args['BetaCM'] = 0.0
    args['A'] = 16
    args['reference'] = 'O16'
    args['custom_valence_space'] = ''
    args['smax'] = 500
    args['emax'] = 3
    args['e3max'] = 0
    args['method'] = 'magnus'
    args['omega_norm_max'] = 0.25
    args['file3e1max'] = 12
    args['file3e2max'] = 28
    args['file3e3max'] = 12
    args['Operators'] = 'Rp2,Rm2'
    args['scratch'] = 'SCRATCH'
    args['use_brueckner_bch'] = False
    args['denominator_delta_orbit'] = 'all'
    args['outputdir'] = './output'
def SetHamilFile(hw, emax2=0, e2max2=0, \
        emax3=0, e2max3=0, e3max3=0):
    args['hw'] = hw
    #args['2bme'] = HOME + 'TwBME-HO_NN-only_' + \
    #        'N3LO_EM500_srg2.00_hw' + \
    #        str(hw) + '_emax' + str(emax2) + \
    #        '_e2max' + str(e2max2) + '.kshell.snt'
    args['2bme'] = HOME + "tbme_em500n3lo_srg2.0hw"+str(hw)+"emax"+str(emax2)+".snt"
    args['2bme'] = "/Users/sotauu/Desktop/local_chiEFTint.jl/tbme_em500n3lo_srg2.0hw20emax12.snt"
    args['2bme'] = "/Users/sotauu/Desktop/local_chiEFTint.jl/tbme_em500n3lo2n3n_srg2.0hw20emax8.snt"
    
    #args['3bme'] = HOME + 'ThBME-srg2.00_cD-0.20cE0.098_lam400_e3max' + \
    #        str(e3max3) + '_hw'+ str(hw) + \
    #        '_NNN-full.txt'
    args['3bme'] = 'none'
    args['file3e1max'] = emax3
    args['file3e2max'] = e2max3
    args['file3e3max'] = e3max3
def SetModelSpace(emax=0, e3max=0):
    args['emax'] = emax
    args['e3max'] = e3max
def SetModelSpaceHF(emax=0, e2max=0, e3max=0):
    args['emax_hf'] = emax
    args['e2max_hf'] = e2max
    args['e3max'] = e3max
def SetModelSpaceIMSRG(emax=0, e2max=0):
    args['emax_imsrg'] = emax
    args['e2max_imsrg'] = e2max

def SetArg(key, value):
    if(key == 'other'):
        arg = args[key]
        arg += ' ' + value
        args[key] = arg
        return
    args[key] = value
def GetCommand():
    if(not os.path.isdir(args['outputdir'])):
        subprocess.call('mkdir -p ' + args['outputdir'], shell=True)
    if(not os.path.isdir(args['scratch'])):
        subprocess.call('mkdir -p ' + args['scratch'], shell=True)
    cmd = exe + ' '
    for key in args.keys():
        if(key == 'other'): cmd += args[key] + ' '
        else:
            cmd += key + '=' + str(args[key]) + ' '
    return cmd

valence=''
for orb in vplist:
    porb = 'p' + orb + ','
    valence += porb
for orb in vnlist:
    norb = 'n' + orb + ','
    valence += norb

SetDefaultArgs()

elist = [8];hwlist = [20]
Nucl = 'He4';A = 4;Core = 'He4'
#Nucl = 'O16';A = 16;Core = 'O16'
Nucl = 'Ca40';A = 40;Core = 'Ca40'

SetArg('fmt2','tokyo')
SetArg('fmt3','tokyo')
#SetArg('fmt2','myg-ascii')
#SetArg('fmt3','myg-ascii')
SetArg('outputdir', Nucl)
SetArg('basis','HF')
#SetArg('basis','HO') # what is the difference?
SetArg('valence_file_format','tokyo')
SetArg('custom_valence_space', Core+','+valence)
SetArg('A', A)
SetArg('reference', Nucl)
SetArg('Operators', '')
SetArg('dsmax', 0.5)

for emax in elist:
    for hw in hwlist:
        SetArg('denominator_delta', hw)

        SetModelSpace(emax)
        SetHamilFile(hw, emax2=emax, e2max2=2*emax)

        cmd = GetCommand()
        subprocess.call(cmd, shell=True)

