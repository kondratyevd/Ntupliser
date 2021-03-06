
#ifndef MU_PAIR_HELPER
#define MU_PAIR_HELPER

#include "Ntupliser/DiMuons/interface/CommonIncludes.h"
#include "Ntupliser/DiMuons/interface/MuPairInfo.h"
#include "Ntupliser/DiMuons/interface/MuonInfo.h"

typedef struct {
  TLorentzVector nom;
  TLorentzVector up1;
  TLorentzVector down1;
  TLorentzVector up2;
  TLorentzVector down2;
} pairVecSys;

typedef struct {
  TLorentzVector nom;
  TLorentzVector up;
  TLorentzVector down;
} muVecSys;

void FillMuPairInfos( MuPairInfos& _pairInfos, const MuonInfos _muonInfos, pat::MuonCollection muonsSelected);

bool pair_is_OS( std::pair< bool, std::pair<int, int> > i,
		 std::pair< bool, std::pair<int, int> > j );

void FillMuPairMasses( muVecSys& mu1_vec, muVecSys& mu2_vec, pairVecSys& pair_vec, 
		       double& massErr, const double MASS_MUON,
		       const MuonInfo _mu1, const MuonInfo _mu2,
		       const double _mu1_pt, const double _mu2_pt,
		       const double _mu1_ptErr, const double _mu2_ptErr );

double cosThetaCS(TLorentzVector *v1,TLorentzVector *v2, double sqrtS);
double phiCS(TLorentzVector *v1,TLorentzVector *v2, double sqrtS);
double phiCS_Zprime(double px_mum, double py_mum, double px_mup, 
		    double py_mup, double pt_dil, double eta_dil,
		    double phi_dil, double mass_dil);

#endif  // #ifndef MU_PAIR_HELPER
