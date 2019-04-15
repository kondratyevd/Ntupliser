#include "Ntupliser/DiMuons/interface/MuPairHelper.h"
#include "Ntupliser/DiMuons/interface/CompositeCandMassResolution.h"

void FillMuPairInfos( MuPairInfos& _pairInfos, const MuonInfos _muonInfos, pat::MuonCollection muonsSelected) {

  _pairInfos.clear();
  if (_muonInfos.size() < 2)
    return;

  double const MASS_MUON  = 0.105658367; // GeV/c^2
  double const PI         = 3.14159265359;

  // 4-vectors: nominal, mu1_up, mu1_down, mu2_up, mu2_down
  muVecSys mu1_vec;
  muVecSys mu2_vec;
  pairVecSys pair_vec;
  double massErr;

  std::vector< std::pair< bool, std::pair<int, int> > > isOS;
  isOS.clear();

  // Sort pairs by OS/SS, then highest mu1 pT, then highest mu2 pT
  // Muons come sorted by pT, so only need to stable sort by OS/SS 
  for (int i = 0; i < int(_muonInfos.size()); i++) {
    for (int j = i+1; j < int(_muonInfos.size()); j++) {
      bool osPair = (_muonInfos.at(i).charge + _muonInfos.at(j).charge == 0);
      isOS.push_back(std::make_pair(osPair, std::make_pair(i, j)));
    }
  }

  int nPairs = isOS.size();
  std::stable_sort( isOS.begin(), isOS.end(), pair_is_OS );
  
  for (int i = 0; i < int(isOS.size()); i++) {
    
    MuPairInfo _pairInfo;
    
    int iMu1 = isOS.at(i).second.first;
    int iMu2 = isOS.at(i).second.second;
    
    _pairInfo.iMu1 = iMu1;
    _pairInfo.iMu2 = iMu2;
    _pairInfo.charge = _muonInfos.at(iMu1).charge + _muonInfos.at(iMu2).charge; 

    FillMuPairMasses( mu1_vec, mu2_vec, pair_vec, massErr, MASS_MUON, 
		    _muonInfos.at(iMu1), _muonInfos.at(iMu2),
		    _muonInfos.at(iMu1).pt, _muonInfos.at(iMu2).pt,
		    _muonInfos.at(iMu1).ptErr, _muonInfos.at(iMu2).ptErr );

    _pairInfo.mass    = pair_vec.nom.M();
    _pairInfo.massErr = massErr;
    _pairInfo.pt      = pair_vec.nom.Pt();
    _pairInfo.eta     = pair_vec.nom.PseudoRapidity();
    _pairInfo.rapid   = pair_vec.nom.Rapidity();
    _pairInfo.phi     = pair_vec.nom.Phi();

    double _dR   = mu1_vec.nom.DeltaR(mu2_vec.nom);
    double _dEta = mu1_vec.nom.PseudoRapidity() - mu2_vec.nom.PseudoRapidity();
    double _dPhi = mu1_vec.nom.DeltaPhi(mu2_vec.nom);

    double _dThetaStarEta = acos( tanh(_dEta/2) );
    double _dPhiStar      = tan( (PI - fabs(_dPhi)) / 2) * sin(_dThetaStarEta);

    _pairInfo.dR   = _dR;
    _pairInfo.dEta = _dEta;
    _pairInfo.dPhi = _dPhi;
    _pairInfo.dPhiStar = _dPhiStar;
    _pairInfo.dThetaStarEta = _dThetaStarEta;

    if ( _muonInfos.at(iMu1).pt_PF > 0 && _muonInfos.at(iMu2).pt_PF > 0 ) {
      FillMuPairMasses( mu1_vec, mu2_vec, pair_vec, massErr, MASS_MUON, 
		      _muonInfos.at(iMu1), _muonInfos.at(iMu2),
		      _muonInfos.at(iMu1).pt_PF, _muonInfos.at(iMu2).pt_PF,
		      _muonInfos.at(iMu1).ptErr, _muonInfos.at(iMu2).ptErr );
      
      _pairInfo.mass_PF    = pair_vec.nom.M();
      _pairInfo.massErr_PF = massErr;
      _pairInfo.pt_PF      = pair_vec.nom.Pt();
    }

    if ( _muonInfos.at(iMu1).pt_trk > 0 && _muonInfos.at(iMu2).pt_trk > 0 ) {
      FillMuPairMasses( mu1_vec, mu2_vec, pair_vec, massErr, MASS_MUON, 
		      _muonInfos.at(iMu1), _muonInfos.at(iMu2),
		      _muonInfos.at(iMu1).pt_trk, _muonInfos.at(iMu2).pt_trk,
		      _muonInfos.at(iMu1).ptErr_trk, _muonInfos.at(iMu2).ptErr_trk );
      
      _pairInfo.mass_trk    = pair_vec.nom.M();
      _pairInfo.massErr_trk = massErr;
      _pairInfo.pt_trk      = pair_vec.nom.Pt();
    }
    
    if ( _muonInfos.at(iMu1).pt_KaMu > 0 && _muonInfos.at(iMu2).pt_KaMu > 0 ) {
      FillMuPairMasses( mu1_vec, mu2_vec, pair_vec, massErr, MASS_MUON, 
		      _muonInfos.at(iMu1), _muonInfos.at(iMu2),
		      _muonInfos.at(iMu1).pt_KaMu, _muonInfos.at(iMu2).pt_KaMu,
		      _muonInfos.at(iMu1).ptErr_KaMu, _muonInfos.at(iMu2).ptErr_KaMu );
      
      _pairInfo.mass_KaMu    = pair_vec.nom.M();
      _pairInfo.massErr_KaMu = massErr;
      _pairInfo.pt_KaMu      = pair_vec.nom.Pt();
    }
    
    if ( _muonInfos.at(iMu1).pt_KaMu_clos_up > 0 && _muonInfos.at(iMu2).pt_KaMu_clos_up > 0 ) {
      FillMuPairMasses( mu1_vec, mu2_vec, pair_vec, massErr, MASS_MUON, 
		      _muonInfos.at(iMu1), _muonInfos.at(iMu2),
		      _muonInfos.at(iMu1).pt_KaMu_clos_up, _muonInfos.at(iMu2).pt_KaMu_clos_up,
		      _muonInfos.at(iMu1).ptErr_KaMu, _muonInfos.at(iMu2).ptErr_KaMu );
      
      _pairInfo.mass_KaMu_clos_up = pair_vec.nom.M();
      _pairInfo.pt_KaMu_clos_up   = pair_vec.nom.Pt();
    }

    if ( _muonInfos.at(iMu1).pt_KaMu_clos_down > 0 && _muonInfos.at(iMu2).pt_KaMu_clos_down > 0 ) {
      FillMuPairMasses( mu1_vec, mu2_vec, pair_vec, massErr, MASS_MUON, 
		      _muonInfos.at(iMu1), _muonInfos.at(iMu2),
		      _muonInfos.at(iMu1).pt_KaMu_clos_down, _muonInfos.at(iMu2).pt_KaMu_clos_down,
		      _muonInfos.at(iMu1).ptErr_KaMu, _muonInfos.at(iMu2).ptErr_KaMu );
      
      _pairInfo.mass_KaMu_clos_down = pair_vec.nom.M();
      _pairInfo.pt_KaMu_clos_down   = pair_vec.nom.Pt();
    }

    if ( _muonInfos.at(iMu1).pt_KaMu_sys_up > 0 && _muonInfos.at(iMu2).pt_KaMu_sys_up > 0 ) {
      FillMuPairMasses( mu1_vec, mu2_vec, pair_vec, massErr, MASS_MUON, 
		      _muonInfos.at(iMu1), _muonInfos.at(iMu2),
		      _muonInfos.at(iMu1).pt_KaMu_sys_up, _muonInfos.at(iMu2).pt_KaMu_sys_up,
		      _muonInfos.at(iMu1).ptErr_KaMu, _muonInfos.at(iMu2).ptErr_KaMu );
      
      _pairInfo.mass_KaMu_sys_up = pair_vec.nom.M();
      _pairInfo.pt_KaMu_sys_up   = pair_vec.nom.Pt();
    }

    if ( _muonInfos.at(iMu1).pt_KaMu_sys_down > 0 && _muonInfos.at(iMu2).pt_KaMu_sys_down > 0 ) {
      FillMuPairMasses( mu1_vec, mu2_vec, pair_vec, massErr, MASS_MUON, 
		      _muonInfos.at(iMu1), _muonInfos.at(iMu2),
		      _muonInfos.at(iMu1).pt_KaMu_sys_down, _muonInfos.at(iMu2).pt_KaMu_sys_down,
		      _muonInfos.at(iMu1).ptErr_KaMu, _muonInfos.at(iMu2).ptErr_KaMu );
      
      _pairInfo.mass_KaMu_sys_down = pair_vec.nom.M();
      _pairInfo.pt_KaMu_sys_down   = pair_vec.nom.Pt();
    }


    // Kinemtic Fit
    if ( _muonInfos.at(iMu1).pt_kinfit > 0 && _muonInfos.at(iMu2).pt_kinfit > 0 ) {
      FillMuPairMasses( mu1_vec, mu2_vec, pair_vec, massErr, MASS_MUON, 
		      _muonInfos.at(iMu1), _muonInfos.at(iMu2),
		      _muonInfos.at(iMu1).pt_kinfit, _muonInfos.at(iMu2).pt_kinfit,
		      _muonInfos.at(iMu1).ptErr_kinfit, _muonInfos.at(iMu2).ptErr_kinfit );
      
      _pairInfo.mass_kinfit    = pair_vec.nom.M();
      _pairInfo.massErr_kinfit = massErr;
      _pairInfo.pt_kinfit      = pair_vec.nom.Pt();
    }

    // TODO: Implement systematics for kinfit
 
    // Rochester
    if ( _muonInfos.at(iMu1).pt_Roch > 0 && _muonInfos.at(iMu2).pt_Roch > 0 ) {
      FillMuPairMasses( mu1_vec, mu2_vec, pair_vec, massErr, MASS_MUON, 
		      _muonInfos.at(iMu1), _muonInfos.at(iMu2),
		      _muonInfos.at(iMu1).pt_Roch, _muonInfos.at(iMu2).pt_Roch,
		      _muonInfos.at(iMu1).ptErr_Roch, _muonInfos.at(iMu2).ptErr_Roch );
      
      _pairInfo.mass_Roch    = pair_vec.nom.M();
      _pairInfo.massErr_Roch = massErr;
      _pairInfo.pt_Roch      = pair_vec.nom.Pt();
    }
    
    if ( _muonInfos.at(iMu1).pt_Roch_sys_up > 0 && _muonInfos.at(iMu2).pt_Roch_sys_up > 0 ) {
      FillMuPairMasses( mu1_vec, mu2_vec, pair_vec, massErr, MASS_MUON, 
		      _muonInfos.at(iMu1), _muonInfos.at(iMu2),
		      _muonInfos.at(iMu1).pt_Roch_sys_up, _muonInfos.at(iMu2).pt_Roch_sys_up,
		      _muonInfos.at(iMu1).ptErr_Roch, _muonInfos.at(iMu2).ptErr_Roch );
      
      _pairInfo.mass_Roch_sys_up = pair_vec.nom.M();
      _pairInfo.pt_Roch_sys_up   = pair_vec.nom.Pt();
    }
    
    if ( _muonInfos.at(iMu1).pt_Roch_sys_down > 0 && _muonInfos.at(iMu2).pt_Roch_sys_down > 0 ) {
      FillMuPairMasses( mu1_vec, mu2_vec, pair_vec, massErr, MASS_MUON, 
		      _muonInfos.at(iMu1), _muonInfos.at(iMu2),
		      _muonInfos.at(iMu1).pt_Roch_sys_down, _muonInfos.at(iMu2).pt_Roch_sys_down,
		      _muonInfos.at(iMu1).ptErr_Roch, _muonInfos.at(iMu2).ptErr_Roch );
      
      _pairInfo.mass_Roch_sys_down = pair_vec.nom.M();
      _pairInfo.pt_Roch_sys_down   = pair_vec.nom.Pt();
    }

    // Add event-by-event mass resolution
    pat::Muon mu1 = muonsSelected.at(iMu1);
    pat::Muon mu2 = muonsSelected.at(iMu2);
    pat::CompositeCandidate mumu;
    mumu.addDaughter(mu1, "mu1");
    mumu.addDaughter(mu2, "mu2");
    AddFourMomenta addP4;
    addP4.set(mumu);
    CompositeCandMassResolution *res = new CompositeCandMassResolution();
    double mass_res = res->getMassResolution(mumu);
    //std::cout << "mu1 eta = " << mu1.eta() << "mu2 eta = "<< mu2.eta() << " resolution = "<<mass_res<<std::endl;
    _pairInfo.mass_res = mass_res;


    // Add angular variables in Collins-Soper frame: implementation from ChargedHiggs
    _pairInfo.cosThetaCS = cosThetaCS(&mu1_vec.nom, &mu2_vec.nom, 13.);
    //_pairInfo.phiCS = phiCS(&mu1_vec.nom, &mu2_vec.nom, 13.);
    
    //implementation of phiCS from Zprime analysis
    double phicszprime;
    if (_muonInfos.at(iMu1).charge < 0 ){
      phicszprime = phiCS_Zprime(mu1_vec.nom.Px(), mu1_vec.nom.Py(), mu2_vec.nom.Px(), mu2_vec.nom.Py(), _pairInfo.pt_Roch, _pairInfo.eta,  _pairInfo.phi, _pairInfo.mass_Roch);

    } else {
      phicszprime    = phiCS_Zprime(mu2_vec.nom.Px(), mu2_vec.nom.Py(), mu1_vec.nom.Px(), mu1_vec.nom.Py(), _pairInfo.pt_Roch, _pairInfo.eta,  _pairInfo.phi, _pairInfo.mass_Roch);
    }
    _pairInfo.phiCS = phicszprime;
    //    std::cout << "mu1_pt = " << mu1_vec.nom.Pt()<<std::endl;
    //std::cout << "cosThetaCS = "<<_pairInfo.cosThetaCS<<", phiCS = "<<_pairInfo.phiCS<< ", phiCS as in Zprime = "<< phicszprime <<std::endl;

    _pairInfos.push_back( _pairInfo );
  } // End loop: for (int i = 0; i < isOS.size(); i++)
  
  if ( int(_pairInfos.size()) != nPairs )
    std::cout << "Bizzare error: muon _pairInfos.size() = " << _pairInfos.size()
	      << ", nPairs = " << nPairs << std::endl;
  
} // End void FillMuPairInfos( MuPairInfos& _pairInfos, const MuonInfos _muonInfos )

bool pair_is_OS( std::pair< bool, std::pair<int, int> > i, 
		 std::pair< bool, std::pair<int, int> > j) {
  return (i.first || !j.first);
}

void FillMuPairMasses( muVecSys& mu1_vec, muVecSys& mu2_vec, pairVecSys& pair_vec, 
		     double& massErr, const double MASS_MUON,
		     const MuonInfo _mu1, const MuonInfo _mu2, 
		     const double _mu1_pt, const double _mu2_pt,
		     const double _mu1_ptErr, const double _mu2_ptErr ) {

  mu1_vec.nom.SetPtEtaPhiM(_mu1_pt, _mu1.eta, _mu1.phi, MASS_MUON);
  mu2_vec.nom.SetPtEtaPhiM(_mu2_pt, _mu2.eta, _mu2.phi, MASS_MUON);
  
  mu1_vec.up.SetPtEtaPhiM(_mu1_pt + _mu1_ptErr, _mu1.eta, _mu1.phi, MASS_MUON);
  mu2_vec.up.SetPtEtaPhiM(_mu2_pt + _mu2_ptErr, _mu2.eta, _mu2.phi, MASS_MUON);
  
  mu1_vec.down.SetPtEtaPhiM(_mu1_pt - _mu1_ptErr, _mu1.eta, _mu1.phi, MASS_MUON);
  mu2_vec.down.SetPtEtaPhiM(_mu2_pt - _mu2_ptErr, _mu2.eta, _mu2.phi, MASS_MUON);
  
  pair_vec.nom   = mu1_vec.nom  + mu2_vec.nom;
  pair_vec.up1   = mu1_vec.up   + mu2_vec.nom;
  pair_vec.down1 = mu1_vec.down + mu2_vec.nom;
  pair_vec.up2   = mu1_vec.nom  + mu2_vec.up;
  pair_vec.down2 = mu1_vec.nom  + mu2_vec.down;
  
  massErr = sqrt( pow(pair_vec.nom.M() - pair_vec.up1.M(),   2) +
		  pow(pair_vec.nom.M() - pair_vec.down1.M(), 2) + 
		  pow(pair_vec.nom.M() - pair_vec.up2.M(),   2) + 
		  pow(pair_vec.nom.M() - pair_vec.down2.M(), 2) ) / 2.;
  
} // End void FillMuPairMasses()

double cosThetaCS(TLorentzVector *v1,TLorentzVector *v2, double sqrtS) {
  double pMass = 0.938272081; // mass from PDG
  TLorentzVector b1,b2,v;
  b1.SetPx(0); b1.SetPy(0);
  b2.SetPx(0); b2.SetPy(0);
  double beamE = 500.*sqrtS; // 1/2 sqrtS in GeV
  b1.SetPz( std::hypot(beamE,pMass)); b1.SetE(beamE);
  b2.SetPz(-std::hypot(beamE,pMass)); b2.SetE(beamE);

  v=*v1+*v2;
  TVector3 boostToVFrame = -v.BoostVector();

  // Boost to higgs frame
  TLorentzVector refV_v1 = *v1; refV_v1.Boost(boostToVFrame);
  TLorentzVector refV_b1 = b1; refV_b1.Boost(boostToVFrame);
  TLorentzVector refV_b2 = b2; refV_b2.Boost(boostToVFrame);

  // Getting beam 3-vector from 4-vectors
  TVector3 refV_vb1_direction = refV_b1.Vect().Unit();
  TVector3 refV_vb2_direction = refV_b2.Vect().Unit();

  // Definition of zz directions
  TVector3 direction_cs = (refV_vb1_direction - refV_vb2_direction).Unit(); // CS direction

  return TMath::Cos(direction_cs.Angle(refV_v1.Vect()));
}


double phiCS(TLorentzVector *v1,TLorentzVector *v2, double sqrtS) {
  double pMass = 0.938272081; // mass from PDG
  TLorentzVector b1,b2,v;
  b1.SetPx(0); b1.SetPy(0);
  b2.SetPx(0); b2.SetPy(0);
  double beamE = 500.*sqrtS; // 1/2 sqrtS in GeV
  b1.SetPz( std::hypot(beamE,pMass)); b1.SetE(beamE);
  b2.SetPz(-std::hypot(beamE,pMass)); b2.SetE(beamE);

  v=*v1+*v2;
  TVector3 boostToVFrame = -v.BoostVector();

  // Boost to higgs frame
  TLorentzVector refV_v1 = *v1; refV_v1.Boost(boostToVFrame);
  TLorentzVector refV_b1 = b1; refV_b1.Boost(boostToVFrame);
  TLorentzVector refV_b2 = b2; refV_b2.Boost(boostToVFrame);

  // Getting beam 3-vector from 4-vectors
  TVector3 refV_vb1_direction = refV_b1.Vect().Unit();
  TVector3 refV_vb2_direction = refV_b2.Vect().Unit();

  // Definition of zz directions
  TVector3 direction_cs = (refV_vb1_direction - refV_vb2_direction).Unit(); // CS direction
  //return (xAxis,yAxis,CSAxis)
  auto yAxis=(refV_b1.Vect().Unit()-refV_b2.Vect().Unit()).Unit();
  auto xAxis=yAxis.Cross(direction_cs).Unit();
  double phi=std::atan2(refV_v1.Vect()*yAxis , refV_v1.Vect() *xAxis) ;
  if(phi<0) return phi + 2*TMath::Pi();
  return phi;
}

double phiCS_Zprime(double px_mum, double py_mum, double px_mup, 
		    double py_mup, double pt_dil, double eta_dil,
		    double phi_dil, double mass_dil) {
  // Analytical calculation of collins-soper phi using 4-vector of dilepton 
  // and the px,py components of mu+, and mu- in the lab frame (can also use 
  // dilepton CM frame).

  TLorentzVector v_dil;
  TVector3 v3_beam, v3_R_T;

  // Create 4-vector out of Z' compononets.
  v_dil.SetPtEtaPhiM(pt_dil, eta_dil, phi_dil, mass_dil);

  // Approximate the longitudinal momentum of the quark to be in the same 
  // direction as that of the dilepton (here the beam is defined as the 
  // direction of the quark).
  double beamEnergy = 6500.;
  if (v_dil.Pz() > 0.) v3_beam.SetXYZ(0., 0., beamEnergy);
  else  v3_beam.SetXYZ(0., 0., -beamEnergy);

  // Make a transverse unit vector in the direction of beam x dilepton
  v3_R_T = (v3_beam.Cross(v_dil.Vect())).Unit();

  // Store transverse components of vectors in appropriate containers.
  TVector2 v2_delta_T(px_mum-px_mup, py_mum-py_mup);
  TVector2 v2_Q_T(v_dil.X(), v_dil.Y());
  v2_Q_T = v2_Q_T.Unit();
  TVector2 v2_R_T(v3_R_T.X(), v3_R_T.Y());

  double Q_term = (sqrt((mass_dil*mass_dil + (pt_dil*pt_dil))))/mass_dil;
  double delta_R_term = v2_delta_T*v2_R_T;
  double delta_Q_term = v2_delta_T*v2_Q_T;
  double phi_cs = atan2(Q_term*delta_R_term, delta_Q_term);
  if (phi_cs < 0.) phi_cs += 2*TMath::Pi();

  return phi_cs;
}
