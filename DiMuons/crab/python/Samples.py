
class sample:
    def __init__(self, name='', DAS='', inputDBS='global', nEvt=0, files=[], GT='94X_mc2017_realistic_v12', 
                 JEC='Fall17_17Nov2017_V4_MC', runs=[], JSON=[], isData=False):
        self.name   = name   ## User-assigned dataset name
        self.DAS    = DAS    ## DAS directory
        self.inputDBS = inputDBS  # to be used in crab in case of private production. config.Data.inputDBS = 'global' or 'phys03'
        self.nEvt   = nEvt   ## Number of events in dataset
        self.files  = files  ## Local files for testing
        self.GT     = GT     ## Global tag
        self.JEC    = JEC    ## Jet energy corrections global tag
        self.runs   = runs   ## Run range
        self.JSON   = JSON   ## JSON file
        self.isData = isData ## Is data

# =======================================================================================================
# ------------------------------- DATA ------------------------------------------------------------------
# =======================================================================================================

# The JSON file details the valid lumi sections
## JSON files: https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/
JSON_2017 = ['data/JSON/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt']  ## 41.37 /fb

PB_wdir = '/afs/cern.ch/work/b/bortigno/hmm_940p1_ufl/src/Ntupliser/DiMuons/'

#####################################
###  Global tag and dataset info  ###
#####################################

# 2017 Data
# More info at https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2017Analysis
# Suggested GT (see frontier conditions twiki for more info)
# 94X_dataRun2_ReReco_EOY17_v2

#SingleMuon dataset list
#/SingleMuon/Run2017A-PromptReco-v2/MINIAOD
#/SingleMuon/Run2017A-PromptReco-v3/MINIAOD
#/SingleMuon/Run2017B-04Jul2017-v2/MINIAOD
#/SingleMuon/Run2017B-06Jul2017-v2/MINIAOD
#/SingleMuon/Run2017B-12Sep2017-v1/MINIAOD
#/SingleMuon/Run2017B-17Nov2017-v1/MINIAOD
#/SingleMuon/Run2017B-22Jun2017-v1/MINIAOD
#/SingleMuon/Run2017B-23Jun2017-v1/MINIAOD
#/SingleMuon/Run2017B-PromptReco-v1/MINIAOD
#/SingleMuon/Run2017B-PromptReco-v2/MINIAOD
#/SingleMuon/Run2017C-12Sep2017-v1/MINIAOD
#/SingleMuon/Run2017C-17Nov2017-v1/MINIAOD
#/SingleMuon/Run2017C-PromptReco-v1/MINIAOD
#/SingleMuon/Run2017C-PromptReco-v2/MINIAOD
#/SingleMuon/Run2017C-PromptReco-v3/MINIAOD
#/SingleMuon/Run2017D-17Nov2017-v1/MINIAOD
#/SingleMuon/Run2017D-PromptReco-v1/MINIAOD
#/SingleMuon/Run2017E-17Nov2017-v1/MINIAOD
#/SingleMuon/Run2017E-PromptReco-v1/MINIAOD
#/SingleMuon/Run2017F-17Nov2017-v1/MINIAOD
#/SingleMuon/Run2017F-PromptReco-v1/MINIAOD
#/SingleMuon/Run2017G-PromptReco-v1/MINIAOD
#/SingleMuon/Run2017H-PromptReco-v1/MINIAOD


#ReReco MINIAOD 17Nov2017

SingleMu_2017B = sample ( name   = 'SingleMu_2017B',
                          DAS    = '/SingleMuon/Run2017B-17Nov2017-v1/MINIAOD',
                          GT     = '94X_dataRun2_ReReco_EOY17_v2',
#                          JEC    = '',
                          JSON   = JSON_2017[0],
                          isData = True )

SingleMu_2017C = sample ( name   = 'SingleMu_2017C',
                          DAS    = '/SingleMuon/Run2017C-17Nov2017-v1/MINIAOD',
                          GT     = '94X_dataRun2_ReReco_EOY17_v2',
#                          JEC    = '',
                          JSON   = JSON_2017[0],
                          isData = True )


SingleMu_2017D = sample ( name   = 'SingleMu_2017D',
                          DAS    = '/SingleMuon/Run2017D-17Nov2017-v1/MINIAOD',
                          GT     = '94X_dataRun2_ReReco_EOY17_v2',
#                          JEC    = '',
                          JSON   = JSON_2017[0],
                          isData = True )


SingleMu_2017E = sample ( name   = 'SingleMu_2017E',
                          DAS    = '/SingleMuon/Run2017E-17Nov2017-v1/MINIAOD',
                          GT     = '94X_dataRun2_ReReco_EOY17_v2',
#                          JEC    = '',
                          JSON   = JSON_2017[0],
                          isData = True )


SingleMu_2017F = sample ( name   = 'SingleMu_2017F',
                          DAS    = '/SingleMuon/Run2017F-17Nov2017-v1/MINIAOD',
                          GT     = '94X_dataRun2_ReReco_EOY17_v2',
#                          JEC    = '',
                          JSON   = JSON_2017[0],
                          isData = True )


## Up-to-date GT info: https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2017Analysis

## Jet energy correction info
## https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
## https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections

SingleMu = []  ## All SingleMuon datasets

SingleMu.append(SingleMu_2017B)
SingleMu.append(SingleMu_2017C)
SingleMu.append(SingleMu_2017D)
SingleMu.append(SingleMu_2017E)
SingleMu.append(SingleMu_2017F)

# =======================================================================================================
# ------------------------------- SIGNAL ----------------------------------------------------------------
# =======================================================================================================

## DAS location: dataset=/*HToMuMu_M*_13TeV*/*Fall17*/MINIAODSIM
DAS_era_sig = 'RunIIFall17DRPremix-94X_mc2017_realistic*/MINIAODSIM'

#For the moment there are no signal samples ready. Using the A. Marini private prod
#/GluGlu_HToMuMu_M125_13TeV_amcatnloFXFX_pythia8/amarini-Fall17_94X-MINIAODSIM-65c6b29ab922da94b788da955c09b417/USER
#/VBFH_HToMuMu_M125_13TeV_amcatnloFXFX_pythia8/amarini-Fall17_94X-MINIAODSIM-65c6b29ab922da94b788da955c09b417/USER
#/ttH_HToMuMu_M125_13TeV_powheg_pythia8/amarini-Fall17_94X-MINIAODSIM-65c6b29ab922da94b788da955c09b417/USER
#/WPlusH_HToMuMu_M125_13TeV_powheg_pythia8/amarini-Fall17_94X-MINIAODSIM-65c6b29ab922da94b788da955c09b417/USER
#/WMinusH_HToMuMu_M125_13TeV_powheg_pythia8/amarini-Fall17_94X-MINIAODSIM-65c6b29ab922da94b788da955c09b417/USER
#/ZH_HToMuMu_M125_13TeV_powheg_pythia8/amarini-Fall17_94X-MINIAODSIM-65c6b29ab922da94b788da955c09b417/USER


## Gluon-gluon fusion
H2Mu_gg = sample( name  = 'H2Mu_gg',
                  DAS   = '/GluGlu_HToMuMu_M125_13TeV_amcatnloFXFX_pythia8/amarini-Fall17_94X-MINIAODSIM-65c6b29ab922da94b788da955c09b417/USER',
                  inputDBS = 'phys03',
                  nEvt  = -1
                )

#H2Mu_gg_120 = sample( name  = 'H2Mu_gg_120',
#                      DAS   = '/GluGlu_HToMuMu_M120_13TeV_powheg_pythia8/'+DAS_era_sig,
#                      nEvt  = 250000, ## 250 k
#                      files = [ AWB_dir+H2Mu_gg_dir+'C0801715-85C0-E611-97A8-001E67396A18.root' ] )
#
#H2Mu_gg_130 = sample( name  = 'H2Mu_gg_130',
#                      DAS   = '/GluGlu_HToMuMu_M130_13TeV_powheg_pythia8/'+DAS_era_sig,
#                      nEvt  = 250000, ## 250 k
#                      files = [ AWB_dir+H2Mu_gg_dir+'C0801715-85C0-E611-97A8-001E67396A18.root' ] )

## Vector boson fusion
H2Mu_VBF = sample( name = 'H2Mu_VBF',
                   DAS  = '/VBFH_HToMuMu_M125_13TeV_amcatnloFXFX_pythia8/amarini-Fall17_94X-MINIAODSIM-65c6b29ab922da94b788da955c09b417/USER',
                   inputDBS = 'phys03',
                   nEvt = -1 )

#H2Mu_VBF_120 = sample( name = 'H2Mu_VBF_120',
#                       DAS  = '/VBF_HToMuMu_M120_13TeV_powheg_pythia8/'+DAS_era_sig,
#                       nEvt = 249200 ) ## 250 k
#
#H2Mu_VBF_130 = sample( name = 'H2Mu_VBF_130',
#                       DAS  = '/VBF_HToMuMu_M130_13TeV_powheg_pythia8/'+DAS_era_sig,
#                       nEvt = 249200 ) ## 250 k

## WH (+)
H2Mu_WH_pos = sample( name = 'H2Mu_WH_pos',
                      DAS  = '/WPlusH_HToMuMu_M125_13TeV_powheg_pythia8/amarini-Fall17_94X-MINIAODSIM-65c6b29ab922da94b788da955c09b417/USER',
                      inputDBS = 'phys03',
                      nEvt = -1 )

#H2Mu_WH_pos_120 = sample( name = 'H2Mu_WH_pos_120',
#                          DAS  = '/WPlusH_HToMuMu_M120_13TeV_powheg_pythia8/'+DAS_era_sig,
#                          nEvt = 124547 ) ## 125 k
#
#H2Mu_WH_pos_130 = sample( name = 'H2Mu_WH_pos_130',
#                          DAS  = '/WPlusH_HToMuMu_M130_13TeV_powheg_pythia8/'+DAS_era_sig,
#                          nEvt = 124547 ) ## 125 k

## WH (-)
H2Mu_WH_neg = sample( name = 'H2Mu_WH_neg',
                      DAS  = '/WMinusH_HToMuMu_M125_13TeV_powheg_pythia8/amarini-Fall17_94X-MINIAODSIM-65c6b29ab922da94b788da955c09b417/USER',
                      inputDBS = 'phys03',
                      nEvt = -1 )

#H2Mu_WH_neg_120 = sample( name = 'H2Mu_WH_neg_120',
#                          DAS  = '/WMinusH_HToMuMu_M120_13TeV_powheg_pythia8/'+DAS_era_sig,
#                          nEvt = 125000 ) ## 125 k
#
#H2Mu_WH_neg_130 = sample( name = 'H2Mu_WH_neg_130',
#                          DAS  = '/WMinusH_HToMuMu_M130_13TeV_powheg_pythia8/'+DAS_era_sig,
#                          nEvt = 125000 ) ## 125 k

## ZH
H2Mu_ZH = sample( name = 'H2Mu_ZH',
                  DAS  = '/ZH_HToMuMu_M125_13TeV_powheg_pythia8/amarini-Fall17_94X-MINIAODSIM-65c6b29ab922da94b788da955c09b417/USER',
                  inputDBS = 'phys03',
                  nEvt = 249748 ) ## 250 k

#H2Mu_ZH_120 = sample( name = 'H2Mu_ZH_120',
#                  DAS  = '/ZH_HToMuMu_M120_13TeV_powheg_pythia8/'+DAS_era_sig,
#                  nEvt = 249748 ) ## 250 k
#
#H2Mu_ZH_130 = sample( name = 'H2Mu_ZH_130',
#                      DAS  = '/ZH_HToMuMu_M130_13TeV_powheg_pythia8/'+DAS_era_sig,
#                      nEvt = 249748 ) ## 250 k


## ttH
H2Mu_ttH = sample( name = 'H2Mu_ttH',
                   DAS  = '/ttH_HToMuMu_M125_13TeV_powheg_pythia8/amarini-Fall17_94X-MINIAODSIM-65c6b29ab922da94b788da955c09b417/USER',
                   inputDBS = 'phys03',
                   nEvt = -1 )

###############
# Dark Photon #
###############

#ZdToMuMu_M20_eps0p02_eta2p6 = sample( name  = 'Zd2Mu_20',
#                                      DAS   = '/DarkPhoton/avartak-ZdToMuMu-M20-eps0p02-etal2p6_SCOUT-d069282eca76b6fde70ab5192475fd6a/USER',
#                                      nEvt  = 250000,
#                                      files = [ PB_wdir + 'Zd2Mu_M20_test.root',
#                                      '/store/user/avartak/DarkPhoton/ZdToMuMu-M20-eps0p02-etal2p6_SCOUT/170901_214525/0000/genscout_1.root' ])
#
#Zd150 = sample( name = 'Zd150',
#                      DAS  = '/DarkPhoton/avartak-ZdToMuMu-M150-eps0p02_MINIAOD-230b6435bde4f6030b269a9cb8e2b63c/USER',
#                      nEvt = 249748 ) ## 250 k




Signal = []  ## All H2Mu signal samples
Signal.append(H2Mu_gg)
Signal.append(H2Mu_VBF)
Signal.append(H2Mu_WH_pos)
Signal.append(H2Mu_WH_neg)
Signal.append(H2Mu_ZH)
# Signal.append(H2Mu_ttH)

#Signal.append(H2Mu_gg_120)
#Signal.append(H2Mu_VBF_120)
#Signal.append(H2Mu_WH_pos_120)
#Signal.append(H2Mu_WH_neg_120)
#Signal.append(H2Mu_ZH_120)
#
#Signal.append(H2Mu_gg_130)
#Signal.append(H2Mu_VBF_130)
#Signal.append(H2Mu_WH_pos_130)
#Signal.append(H2Mu_WH_neg_130)
#Signal.append(H2Mu_ZH_130)



# =======================================================================================================
# ------------------------------- BACKGROUND ------------------------------------------------------------
# =======================================================================================================

Background = []  ## All H2Mu background samples

###################
###  Drell-Yan  ###
###################

# Monte Carlo Fall17 status of production
# https://cms-pdmv.cern.ch/pmp/historical?r=RunIIFall17DRPremix,RunIIFall17DRStdmix
# Recommended global tags
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
# for 944 and higher: Analysis level GTs with updates of Jet probability calibration and Jet energy corrections on the top of 94X_mc2017_realistic_v10 GT
# 94X_mc2017_realistic_v12
# for 940 and higher: final GT containing the latest MC conditions for 2017
# 94X_mc2017_realistic_v10

# DYJetsToLL Fall 17 list MINIAODSIM - pb 21.02.2017
# DAS query dataset=/DYJetsToLL*/RunIIFall17*/MINIAODSIM
#/DYJetsToLL_M-4to50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
#/DYJetsToLL_M-4to50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
#/DYJetsToLL_M-4to50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10_ext1-v1/MINIAODSIM
#/DYJetsToLL_M-4to50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2/MINIAODSIM
#/DYJetsToLL_M-4to50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10_ext1-v1/MINIAODSIM
#/DYJetsToLL_M-4to50_HT-600toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
#/DYJetsToLL_M-4to50_HT-70to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
#/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
#/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10_ext1-v1/MINIAODSIM
#/DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
#/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
#/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
#/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10_ext1-v1/MINIAODSIM
#/DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
#/DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
#/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
#  Creation time: 2017-12-05 21:12:33, Dataset size: 1.1TB, Number of blocks: 71, Number of events: 26923935, Number of files: 372, Physics group: NoGroup, Status: VALID, Type: mc
#/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10_ext1-v1/MINIAODSIM
#  Creation time: 2018-01-04 06:26:07, Dataset size: 7.4TB, Number of blocks: 203, Number of events: 185998625, Number of files: 2752, Physics group: NoGroup, Status: VALID, Type: mc
#/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-RECOSIMstep_94X_mc2017_realistic_v10-v1/MINIAODSIM
#/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-RECOSIMstep_94X_mc2017_realistic_v10_ext1-v1/MINIAODSIM
#/DYJetsToLL_M-50_Zpt-150toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-RECOSIMstep_94X_mc2017_realistic_v10-v1/MINIAODSIM


# DYJets inclusive Fall17
ZJets_AMC = sample( name = 'ZJets_AMC',
                    DAS  = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10_ext1-v1/MINIAODSIM',
                    nEvt = -1
                  )

# DYJets HT bins Fall17
#
#/DYJetsToLL_M-4to50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10_ext1-v1/ AODSIM;
#SUS-RunIIFall17DRPremix-00008
#94%
#/DYJetsToLL_M-4to50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10_ext1-v1/ AODSIM;
#SUS-RunIIFall17DRPremix-00009
#97%
#/DYJetsToLL_M-4to50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10_ext1-v1/ AODSIM;
#SUS-RunIIFall17DRPremix-00010
#80%
#/DYJetsToLL_M-4to50_HT-600toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10_ext1-v1/ AODSIM;
#SUS-RunIIFall17DRPremix-00011
#100%
#/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10_ext1-v1/ AODSIM;
#SUS-RunIIFall17DRPremix-00013
#100%
#/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10_ext1-v1/ AODSIM;
#SUS-RunIIFall17DRPremix-00014
#83%
#/DYJetsToLL_M-4to50_HT-70to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10_ext1-v1/ AODSIM;



## Mass bins
#/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10-v2/AODSIM; # NO MINIAOD for the moment - pb 21.02.2018
#/ZToMuMu_NNPDF31_13TeV-powheg_M_120_200/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2/MINIAODSIM
#/ZToMuMu_NNPDF31_13TeV-powheg_M_1400_2300/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
#/ZToMuMu_NNPDF31_13TeV-powheg_M_200_400/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
#/ZToMuMu_NNPDF31_13TeV-powheg_M_3500_4500/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
#/ZToMuMu_NNPDF31_13TeV-powheg_M_400_800/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2/MINIAODSIM
#/ZToMuMu_NNPDF31_13TeV-powheg_M_4500_6000/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
#/ZToMuMu_NNPDF31_13TeV-powheg_M_6000_Inf/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
#/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM


#Background.append(ZJets_MG)
#Background.append(ZJets_MG_HER)
Background.append(ZJets_AMC)
#Background.append(ZJets_hiM)
#Background.append(ZJets_MG_HT_70_100)
#Background.append(ZJets_MG_HT_100_200_A)
#Background.append(ZJets_MG_HT_100_200_B)
#Background.append(ZJets_MG_HT_200_400_A)
#Background.append(ZJets_MG_HT_200_400_B)
#Background.append(ZJets_MG_HT_400_600_A)
#Background.append(ZJets_MG_HT_400_600_B)
#Background.append(ZJets_MG_HT_600_800)
#Background.append(ZJets_MG_HT_800_1200)
#Background.append(ZJets_MG_HT_1200_2500)
#Background.append(ZJets_MG_HT_2500_inf)

####################
###  Single top  ###
####################

#Background.append(tW_pos_1)
#Background.append(tW_pos_2)
#Background.append(tW_neg_1)
#Background.append(tW_neg_2)

###############
###  ttbar  ###
###############

# TTbar inclusive
# Dataset: /TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
# Creation time: 2018-01-08 16:44:57, Dataset size: 8.3TB, Number of blocks: 203, Number of events: 153596015, Number of files: 2849, Physics group: NoGroup, Status: VALID, Type: mc

tt = sample( name = 'tt',
             DAS  = '/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM',
             nEvt = -1
           )

Background.append(tt)
#Background.append(tt_ll_MG_1)
#Background.append(tt_ll_MG_2)
#Background.append(tt_ll_AMC)

#################
###  Diboson  ###
#################

# Background.append(WW)
# Background.append(WW_HER)
# Background.append(WW_up)
# Background.append(WW_down)
# Background.append(WZ_2l)
# Background.append(WZ_3l_AMC)
# Background.append(WZ_3l_POW)
# Background.append(ZZ_2l_2v)
# Background.append(ZZ_2l_2q)
# Background.append(ZZ_4l_AMC)
# Background.append(ZZ_4l_POW)

#################
###  Triboson  ##
#################

#Background.append(WWW)
#Background.append(WWZ)
#Background.append(WZZ)
#Background.append(ZZZ)

#####################
###  Single top+X  ##
#####################

# Background.append(tZq)
# Background.append(tZq_HER)
# Background.append(tZW)

################
###  ttbar+X  ##
################

# Background.append(ttW_1)
# Background.append(ttW_2)
# Background.append(ttZ)
# Background.append(ttH)


DataAndMC = []
DataAndMC.extend(SingleMu)
DataAndMC.extend(Signal)
DataAndMC.extend(Background)

MC = []
MC.extend(Signal)
MC.extend(Background)
