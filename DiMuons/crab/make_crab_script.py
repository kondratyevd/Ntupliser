#! /usr/bin/env python

import sys, os
#adding the working dir path
sys.path.append(os.getcwd())

from python.Samples import *

import time  ## For timestamp on crab jobs
import os  ## For executable permissions on scripts
import subprocess 


## Read the code tag for the production.
def get_prod_version():
  try: ## try with an annotated tag
    print('Getting annotated tag...')
    return subprocess.check_output(["git","describe"]).strip()
  except: 
    print('No annotated tag. Trying with lightweighted tag.')
    pass 
  try: ## otherwise try with a lightweighted tag
    print('Getting lightweighted tag...')
    return subprocess.check_output(["git","describe","--tag"]).strip()
  except:  
    print('No lightweighted tag. Using latest committ.')
    pass  
  ## otherwise pick the latest commit
  return subprocess.check_output(["git","rev-parse","--short"]).strip()

prod_version = get_prod_version()
print("Production using code version {0} starting" .format(prod_version))

samps = []

## Get the samples you want to make a crab config file for 
test_run = False
test_str = '_prod2017_{0}'.format(prod_version)
 
# samps.extend(SingleMu)
# samps.extend(Signal)
# samps.extend(Background)
samps.extend(DataAndMC)

# test_run = True
# test_str = '_v1'
# samps.append(SingleMu_2016C)
# samps.append(H2Mu_gg)
# samps.append(ZJets_MG_HT_2500_inf)

# samps.append(SingleMu_2017B)

#samps.append(Zd2Mu_125)
#samps.append(Zd2Mu_150)

crab_prod_dir = 'crab_%s-%s'%(time.strftime('%Y_%m_%d_%H_%M'),prod_version)
crab_analyzers_dir = crab_prod_dir+'/analyzers'
crab_configs_dir = crab_prod_dir+'/configs'
print('crab production directory = ' + crab_prod_dir)
os.mkdir(crab_prod_dir)
os.mkdir(crab_analyzers_dir)
os.mkdir(crab_configs_dir)

for samp in samps:
    print '\nCreating analyzer and crab config for %s' % samp.name
    ## Open a file for writing
    out_file = open('%s/%s.py' % (crab_analyzers_dir,samp.name), 'w')
    
    ## Open the analyzer template to change
    in_file = open('crab/templates/EDAnalyzer.py', 'r')
    
    ## Read in the template and replace the parameters to
    ## make a python file which we will submit using crab
    for line in in_file:
        line = line.replace('samp.isData', str(samp.isData))
        line = line.replace('samp.name', "'%s'" % samp.name)
        line = line.replace('samp.DAS', "'%s'" % samp.DAS)
        line = line.replace('samp.GT', "'%s'" % samp.GT)
        line = line.replace('samp.files', str(samp.files))
        line = line.replace('samp.JSON', "'%s'" % samp.JSON)
        line = line.replace('samp.inputDBS','%s' %samp.inputDBS)

        out_file.write(line)
    
    ## Close the file
    print '  * Wrote %s' % out_file.name
    out_file.close()
    in_file.close()
    
    ## Open the crab submission template to change
    in_file = open('crab/templates/crab_config.py', 'r')
    out_file = open('%s/%s.py' %(crab_configs_dir,samp.name), 'w')
    
    # Read in the template and replace the parameters to make a
    # crab submission file that uses the above CMSSW analyzer
    for line in in_file:
        if 'requestName' in line:
            line = line.replace("= 'STR'", "= '%s_%s%s'" % (samp.name, time.strftime('%Y_%m_%d_%H_%M'), test_str.replace(".","p")) ) 

        if 'psetName' in line: 
            line = line.replace("= 'STR'", "= '%s/%s.py'" % (crab_analyzers_dir, samp.name) )

        if 'inputDataset' in line:
            line = line.replace("= 'STR'", "= '%s'" % samp.DAS)

        if samp.isData and 'lumiMask' in line: 
            line = line.replace('# config.Data.lumiMask', 'config.Data.lumiMask')
            line = line.replace("= 'STR'", "= '%s'" % samp.JSON)

        if 'splitting' in line:
            if samp.isData:
                line = line.replace("= 'STR'", "= 'LumiBased'")
            else:
                line = line.replace("= 'STR'", "= 'FileBased'")

        if 'unitsPerJob' in line:
            if test_run:
                line = line.replace('= NUM', '= 1')
            elif samp.isData:
                line = line.replace('= NUM', '= 250')  ## 100
            # elif samp.name == 'ZJets_MG' or ('ZJets_MG' in samp.name and '_B' in samp.name) or samp.name == 'ZZ_4l_AMC':
            #     line = line.replace('= NUM', '= 3')  ## 10-file jobs fail with too much RAM
            else:
                line = line.replace('= NUM', '= 5')  ## 5

        if 'inputDBS' in line:
            line = line.replace("= 'DBS'", "= '%s'"  % samp.inputDBS)

        if 'outputDatasetTag' in line:
            line = line.replace("= 'STR'", "= '%s'" % samp.name)

        out_file.write(line)
    
    print '  * Wrote %s' % out_file.name
    out_file.close()
    in_file.close()


print '\nCreating submit_all.sh and check_all.sh scripts'

out_file = open('%s/submit_all.sh' % crab_prod_dir, 'w')
out_file.write('#!/bin/bash\n')
out_file.write('\n')
#out_file.write('source /cvmfs/cms.cern.ch/crab3/crab.csh\n')
out_file.write('voms-proxy-init --voms cms --valid 168:00\n')
out_file.write('\n')
for samp in samps:
    out_file.write('crab submit -c %s/%s.py\n' % (crab_configs_dir, samp.name))
out_file.close()
os.chmod('%s/submit_all.sh' % crab_prod_dir, 0o744)

out_file = open('%s/check_all.sh' % crab_prod_dir, 'w')
out_file.write('#!/bin/bash\n')
out_file.write('\n')
# out_file.write('source /cvmfs/cms.cern.ch/crab3/crab.csh\n')
# out_file.write('voms-proxy-init --voms cms --valid 168:00\n')
out_file.write('\n')
for samp in samps:
    out_file.write('crab status -d %s/logs/crab_%s*\n' % (crab_prod_dir,samp.name))
out_file.close()
os.chmod('%s/check_all.sh' % crab_prod_dir, 0o744)
