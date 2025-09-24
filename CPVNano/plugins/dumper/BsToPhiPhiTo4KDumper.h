//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar  5 09:51:52 2025 by ROOT version 6.12/07
// from TTree BsToPhiPhiTo4KDumper/BsToPhiPhiTo4KDumper
// found on file: bparknano.root
//////////////////////////////////////////////////////////

#ifndef BsToPhiPhiTo4KDumper_h
#define BsToPhiPhiTo4KDumper_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class BsToPhiPhiTo4KDumper : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<UInt_t> run = {fReader, "run"};
   TTreeReaderValue<UInt_t> luminosityBlock = {fReader, "luminosityBlock"};
   TTreeReaderValue<ULong64_t> event = {fReader, "event"};
   TTreeReaderValue<UInt_t> nBsToPhiPhiTo4K = {fReader, "nBsToPhiPhiTo4K"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_beta = {fReader, "BsToPhiPhiTo4K_Bs_beta"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_gamma = {fReader, "BsToPhiPhiTo4K_Bs_gamma"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_charge = {fReader, "BsToPhiPhiTo4K_Bs_charge"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_cos2D = {fReader, "BsToPhiPhiTo4K_Bs_cos2D"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_cxx = {fReader, "BsToPhiPhiTo4K_Bs_cxx"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_cyx = {fReader, "BsToPhiPhiTo4K_Bs_cyx"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_cyy = {fReader, "BsToPhiPhiTo4K_Bs_cyy"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_czx = {fReader, "BsToPhiPhiTo4K_Bs_czx"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_czy = {fReader, "BsToPhiPhiTo4K_Bs_czy"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_czz = {fReader, "BsToPhiPhiTo4K_Bs_czz"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_eta = {fReader, "BsToPhiPhiTo4K_Bs_eta"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxy = {fReader, "BsToPhiPhiTo4K_Bs_lxy"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxyErr = {fReader, "BsToPhiPhiTo4K_Bs_lxyErr"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxy_sig = {fReader, "BsToPhiPhiTo4K_Bs_lxy_sig"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxy_corr = {fReader, "BsToPhiPhiTo4K_Bs_lxy_corr"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxyErr_corr = {fReader, "BsToPhiPhiTo4K_Bs_lxyErr_corr"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxy_sig_corr = {fReader, "BsToPhiPhiTo4K_Bs_lxy_sig_corr"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxyz = {fReader, "BsToPhiPhiTo4K_Bs_lxyz"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxyzErr = {fReader, "BsToPhiPhiTo4K_Bs_lxyzErr"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxyz_corr = {fReader, "BsToPhiPhiTo4K_Bs_lxyz_corr"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxyzErr_corr = {fReader, "BsToPhiPhiTo4K_Bs_lxyzErr_corr"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lx = {fReader, "BsToPhiPhiTo4K_Bs_lx"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_ly = {fReader, "BsToPhiPhiTo4K_Bs_ly"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lx_posbsz = {fReader, "BsToPhiPhiTo4K_Bs_lx_posbsz"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_ly_posbsz = {fReader, "BsToPhiPhiTo4K_Bs_ly_posbsz"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lx_posbspv = {fReader, "BsToPhiPhiTo4K_Bs_lx_posbspv"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_ly_posbspv = {fReader, "BsToPhiPhiTo4K_Bs_ly_posbspv"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_ct_2D_cm = {fReader, "BsToPhiPhiTo4K_Bs_ct_2D_cm"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_ct_2D_cm_posbsz = {fReader, "BsToPhiPhiTo4K_Bs_ct_2D_cm_posbsz"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_ct_2D_cm_posbspv = {fReader, "BsToPhiPhiTo4K_Bs_ct_2D_cm_posbspv"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_ct_3D_cm = {fReader, "BsToPhiPhiTo4K_Bs_ct_3D_cm"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_mass = {fReader, "BsToPhiPhiTo4K_Bs_mass"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_massErr = {fReader, "BsToPhiPhiTo4K_Bs_massErr"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_mass_corr = {fReader, "BsToPhiPhiTo4K_Bs_mass_corr"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_phi = {fReader, "BsToPhiPhiTo4K_Bs_phi"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_pt = {fReader, "BsToPhiPhiTo4K_Bs_pt"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_ptErr = {fReader, "BsToPhiPhiTo4K_Bs_ptErr"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_px = {fReader, "BsToPhiPhiTo4K_Bs_px"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_py = {fReader, "BsToPhiPhiTo4K_Bs_py"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_sv_chi2 = {fReader, "BsToPhiPhiTo4K_Bs_sv_chi2"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_sv_ndof = {fReader, "BsToPhiPhiTo4K_Bs_sv_ndof"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_sv_prob = {fReader, "BsToPhiPhiTo4K_Bs_sv_prob"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_vx = {fReader, "BsToPhiPhiTo4K_Bs_vx"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_vy = {fReader, "BsToPhiPhiTo4K_Bs_vy"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_vz = {fReader, "BsToPhiPhiTo4K_Bs_vz"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_beamspot_x = {fReader, "BsToPhiPhiTo4K_beamspot_x"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_beamspot_y = {fReader, "BsToPhiPhiTo4K_beamspot_y"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_beamspot_z = {fReader, "BsToPhiPhiTo4K_beamspot_z"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_cos_theta_star_phi1 = {fReader, "BsToPhiPhiTo4K_cos_theta_star_phi1"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_cos_theta_star_phi2 = {fReader, "BsToPhiPhiTo4K_cos_theta_star_phi2"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_deltaR_k1k2 = {fReader, "BsToPhiPhiTo4K_deltaR_k1k2"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_deltaR_k1k3 = {fReader, "BsToPhiPhiTo4K_deltaR_k1k3"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_deltaR_k1k4 = {fReader, "BsToPhiPhiTo4K_deltaR_k1k4"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_deltaR_k2k3 = {fReader, "BsToPhiPhiTo4K_deltaR_k2k3"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_deltaR_k2k4 = {fReader, "BsToPhiPhiTo4K_deltaR_k2k4"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_deltaR_k3k4 = {fReader, "BsToPhiPhiTo4K_deltaR_k3k4"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_deltaR_max = {fReader, "BsToPhiPhiTo4K_deltaR_max"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_deltaR_min = {fReader, "BsToPhiPhiTo4K_deltaR_min"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_deltaR_phi1phi2 = {fReader, "BsToPhiPhiTo4K_deltaR_phi1phi2"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1_eta = {fReader, "BsToPhiPhiTo4K_k1_eta"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1_mass = {fReader, "BsToPhiPhiTo4K_k1_mass"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1_phi = {fReader, "BsToPhiPhiTo4K_k1_phi"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1_pt = {fReader, "BsToPhiPhiTo4K_k1_pt"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1k3_mass = {fReader, "BsToPhiPhiTo4K_k1k3_mass"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1k3_pt = {fReader, "BsToPhiPhiTo4K_k1k3_pt"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1k4_mass = {fReader, "BsToPhiPhiTo4K_k1k4_mass"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1k4_pt = {fReader, "BsToPhiPhiTo4K_k1k4_pt"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2_eta = {fReader, "BsToPhiPhiTo4K_k2_eta"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2_mass = {fReader, "BsToPhiPhiTo4K_k2_mass"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2_phi = {fReader, "BsToPhiPhiTo4K_k2_phi"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2_pt = {fReader, "BsToPhiPhiTo4K_k2_pt"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2k3_mass = {fReader, "BsToPhiPhiTo4K_k2k3_mass"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2k3_pt = {fReader, "BsToPhiPhiTo4K_k2k3_pt"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2k4_mass = {fReader, "BsToPhiPhiTo4K_k2k4_mass"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2k4_pt = {fReader, "BsToPhiPhiTo4K_k2k4_pt"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k3_eta = {fReader, "BsToPhiPhiTo4K_k3_eta"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k3_mass = {fReader, "BsToPhiPhiTo4K_k3_mass"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k3_phi = {fReader, "BsToPhiPhiTo4K_k3_phi"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k3_pt = {fReader, "BsToPhiPhiTo4K_k3_pt"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k4_eta = {fReader, "BsToPhiPhiTo4K_k4_eta"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k4_mass = {fReader, "BsToPhiPhiTo4K_k4_mass"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k4_phi = {fReader, "BsToPhiPhiTo4K_k4_phi"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k4_pt = {fReader, "BsToPhiPhiTo4K_k4_pt"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_phi1_eta = {fReader, "BsToPhiPhiTo4K_phi1_eta"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_phi1_mass = {fReader, "BsToPhiPhiTo4K_phi1_mass"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_phi1_phi = {fReader, "BsToPhiPhiTo4K_phi1_phi"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_phi1_pt = {fReader, "BsToPhiPhiTo4K_phi1_pt"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_phi2_eta = {fReader, "BsToPhiPhiTo4K_phi2_eta"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_phi2_mass = {fReader, "BsToPhiPhiTo4K_phi2_mass"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_phi2_phi = {fReader, "BsToPhiPhiTo4K_phi2_phi"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_phi2_pt = {fReader, "BsToPhiPhiTo4K_phi2_pt"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_isMatched = {fReader, "BsToPhiPhiTo4K_isMatched"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_k1_idx = {fReader, "BsToPhiPhiTo4K_k1_idx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_k2_idx = {fReader, "BsToPhiPhiTo4K_k2_idx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_k3_idx = {fReader, "BsToPhiPhiTo4K_k3_idx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_k4_idx = {fReader, "BsToPhiPhiTo4K_k4_idx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_phi1_idx = {fReader, "BsToPhiPhiTo4K_phi1_idx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_phi2_idx = {fReader, "BsToPhiPhiTo4K_phi2_idx"};
   TTreeReaderValue<UInt_t> nPhiToKK = {fReader, "nPhiToKK"};
   TTreeReaderArray<Float_t> PhiToKK_cos_theta_star_k1 = {fReader, "PhiToKK_cos_theta_star_k1"};
   TTreeReaderArray<Float_t> PhiToKK_cos_theta_star_k2 = {fReader, "PhiToKK_cos_theta_star_k2"};
   TTreeReaderArray<Float_t> PhiToKK_deltaR_postfit = {fReader, "PhiToKK_deltaR_postfit"};
   TTreeReaderArray<Float_t> PhiToKK_deltaR_prefit = {fReader, "PhiToKK_deltaR_prefit"};
   TTreeReaderArray<Float_t> PhiToKK_energy_diff_phi_daughters_cm = {fReader, "PhiToKK_energy_diff_phi_daughters_cm"};
   TTreeReaderArray<Float_t> PhiToKK_energy_diff_prefitphi_daughters_lab = {fReader, "PhiToKK_energy_diff_prefitphi_daughters_lab"};
   TTreeReaderArray<Float_t> PhiToKK_p_daughters_cm = {fReader, "PhiToKK_p_daughters_cm"};
   TTreeReaderArray<Float_t> PhiToKK_phi_charge = {fReader, "PhiToKK_phi_charge"};
   TTreeReaderArray<Float_t> PhiToKK_phi_cos2D = {fReader, "PhiToKK_phi_cos2D"};
   TTreeReaderArray<Float_t> PhiToKK_phi_cxx = {fReader, "PhiToKK_phi_cxx"};
   TTreeReaderArray<Float_t> PhiToKK_phi_cyx = {fReader, "PhiToKK_phi_cyx"};
   TTreeReaderArray<Float_t> PhiToKK_phi_cyy = {fReader, "PhiToKK_phi_cyy"};
   TTreeReaderArray<Float_t> PhiToKK_phi_czx = {fReader, "PhiToKK_phi_czx"};
   TTreeReaderArray<Float_t> PhiToKK_phi_czy = {fReader, "PhiToKK_phi_czy"};
   TTreeReaderArray<Float_t> PhiToKK_phi_czz = {fReader, "PhiToKK_phi_czz"};
   TTreeReaderArray<Float_t> PhiToKK_eta = {fReader, "PhiToKK_eta"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_DCASig = {fReader, "PhiToKK_phi_k1_DCASig"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_charge = {fReader, "PhiToKK_phi_k1_charge"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_chi2 = {fReader, "PhiToKK_phi_k1_chi2"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_covLamLam = {fReader, "PhiToKK_phi_k1_covLamLam"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_covLamPhi = {fReader, "PhiToKK_phi_k1_covLamPhi"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_covPhiPhi = {fReader, "PhiToKK_phi_k1_covPhiPhi"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_covQopLam = {fReader, "PhiToKK_phi_k1_covQopLam"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_covQopPhi = {fReader, "PhiToKK_phi_k1_covQopPhi"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_covQopQop = {fReader, "PhiToKK_phi_k1_covQopQop"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_dxy = {fReader, "PhiToKK_phi_k1_dxy"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_dxyS = {fReader, "PhiToKK_phi_k1_dxyS"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_dz = {fReader, "PhiToKK_phi_k1_dz"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_dzS = {fReader, "PhiToKK_phi_k1_dzS"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_eta = {fReader, "PhiToKK_phi_k1_eta"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_mass = {fReader, "PhiToKK_phi_k1_mass"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_normalisedChi2 = {fReader, "PhiToKK_phi_k1_normalisedChi2"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_pdgid = {fReader, "PhiToKK_phi_k1_pdgid"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_phi = {fReader, "PhiToKK_phi_k1_phi"};
   TTreeReaderArray<Float_t> PhiToKK_k1_pt = {fReader, "PhiToKK_k1_pt"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_ptErr = {fReader, "PhiToKK_phi_k1_ptErr"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_validFraction = {fReader, "PhiToKK_phi_k1_validFraction"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_vx = {fReader, "PhiToKK_phi_k1_vx"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_vy = {fReader, "PhiToKK_phi_k1_vy"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_vz = {fReader, "PhiToKK_phi_k1_vz"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_DCASig = {fReader, "PhiToKK_phi_k2_DCASig"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_charge = {fReader, "PhiToKK_phi_k2_charge"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_chi2 = {fReader, "PhiToKK_phi_k2_chi2"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_covLamLam = {fReader, "PhiToKK_phi_k2_covLamLam"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_covLamPhi = {fReader, "PhiToKK_phi_k2_covLamPhi"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_covPhiPhi = {fReader, "PhiToKK_phi_k2_covPhiPhi"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_covQopLam = {fReader, "PhiToKK_phi_k2_covQopLam"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_covQopPhi = {fReader, "PhiToKK_phi_k2_covQopPhi"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_covQopQop = {fReader, "PhiToKK_phi_k2_covQopQop"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_dxy = {fReader, "PhiToKK_phi_k2_dxy"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_dxyS = {fReader, "PhiToKK_phi_k2_dxyS"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_dz = {fReader, "PhiToKK_phi_k2_dz"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_dzS = {fReader, "PhiToKK_phi_k2_dzS"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_eta = {fReader, "PhiToKK_phi_k2_eta"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_mass = {fReader, "PhiToKK_phi_k2_mass"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_normalisedChi2 = {fReader, "PhiToKK_phi_k2_normalisedChi2"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_pdgid = {fReader, "PhiToKK_phi_k2_pdgid"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_phi = {fReader, "PhiToKK_phi_k2_phi"};
   TTreeReaderArray<Float_t> PhiToKK_k2_pt = {fReader, "PhiToKK_k2_pt"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_ptErr = {fReader, "PhiToKK_phi_k2_ptErr"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_validFraction = {fReader, "PhiToKK_phi_k2_validFraction"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_vx = {fReader, "PhiToKK_phi_k2_vx"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_vy = {fReader, "PhiToKK_phi_k2_vy"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_vz = {fReader, "PhiToKK_phi_k2_vz"};
   TTreeReaderArray<Float_t> PhiToKK_phi_lxy = {fReader, "PhiToKK_phi_lxy"};
   TTreeReaderArray<Float_t> PhiToKK_phi_lxy_sig = {fReader, "PhiToKK_phi_lxy_sig"};
   TTreeReaderArray<Float_t> PhiToKK_mass = {fReader, "PhiToKK_mass"};
   TTreeReaderArray<Float_t> PhiToKK_phi_masserr = {fReader, "PhiToKK_phi_masserr"};
   TTreeReaderArray<Float_t> PhiToKK_phi = {fReader, "PhiToKK_phi"};
   TTreeReaderArray<Float_t> PhiToKK_pt = {fReader, "PhiToKK_pt"};
   TTreeReaderArray<Float_t> PhiToKK_sv_chi2 = {fReader, "PhiToKK_sv_chi2"};
   TTreeReaderArray<Float_t> PhiToKK_sv_ndof = {fReader, "PhiToKK_sv_ndof"};
   TTreeReaderArray<Float_t> PhiToKK_sv_prob = {fReader, "PhiToKK_sv_prob"};
   TTreeReaderArray<Float_t> PhiToKK_phi_vx = {fReader, "PhiToKK_phi_vx"};
   TTreeReaderArray<Float_t> PhiToKK_phi_vy = {fReader, "PhiToKK_phi_vy"};
   TTreeReaderArray<Float_t> PhiToKK_phi_vz = {fReader, "PhiToKK_phi_vz"};
   TTreeReaderArray<Float_t> PhiToKK_px_diff_phi_daughters_lab = {fReader, "PhiToKK_px_diff_phi_daughters_lab"};
   TTreeReaderArray<Float_t> PhiToKK_px_diff_prefitphi_daughters_lab = {fReader, "PhiToKK_px_diff_prefitphi_daughters_lab"};
   TTreeReaderArray<Float_t> PhiToKK_py_diff_phi_daughters_lab = {fReader, "PhiToKK_py_diff_phi_daughters_lab"};
   TTreeReaderArray<Float_t> PhiToKK_py_diff_prefitphi_daughters_lab = {fReader, "PhiToKK_py_diff_prefitphi_daughters_lab"};
   TTreeReaderArray<Float_t> PhiToKK_pz_diff_phi_daughters_lab = {fReader, "PhiToKK_pz_diff_phi_daughters_lab"};
   TTreeReaderArray<Float_t> PhiToKK_pz_diff_prefitphi_daughters_lab = {fReader, "PhiToKK_pz_diff_prefitphi_daughters_lab"};
   TTreeReaderArray<Int_t> PhiToKK_isMatched = {fReader, "PhiToKK_isMatched"};
   TTreeReaderArray<Int_t> PhiToKK_k1_idx = {fReader, "PhiToKK_k1_idx"};
   TTreeReaderArray<Int_t> PhiToKK_k2_idx = {fReader, "PhiToKK_k2_idx"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k1_highPurityFlag = {fReader, "PhiToKK_phi_k1_highPurityFlag"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k1_islost = {fReader, "PhiToKK_phi_k1_islost"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k1_ispacked = {fReader, "PhiToKK_phi_k1_ispacked"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k1_ndof = {fReader, "PhiToKK_phi_k1_ndof"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k1_numberOfLostHits = {fReader, "PhiToKK_phi_k1_numberOfLostHits"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k1_numberOfPixelHits = {fReader, "PhiToKK_phi_k1_numberOfPixelHits"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k1_numberOfPixelLayers = {fReader, "PhiToKK_phi_k1_numberOfPixelLayers"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k1_numberOfTrackerLayers = {fReader, "PhiToKK_phi_k1_numberOfTrackerLayers"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k1_numberOfValidHits = {fReader, "PhiToKK_phi_k1_numberOfValidHits"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k1_numberOfValidPixelHits = {fReader, "PhiToKK_phi_k1_numberOfValidPixelHits"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k1_qualityIndex = {fReader, "PhiToKK_phi_k1_qualityIndex"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k2_highPurityFlag = {fReader, "PhiToKK_phi_k2_highPurityFlag"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k2_islost = {fReader, "PhiToKK_phi_k2_islost"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k2_ispacked = {fReader, "PhiToKK_phi_k2_ispacked"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k2_ndof = {fReader, "PhiToKK_phi_k2_ndof"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k2_numberOfLostHits = {fReader, "PhiToKK_phi_k2_numberOfLostHits"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k2_numberOfPixelHits = {fReader, "PhiToKK_phi_k2_numberOfPixelHits"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k2_numberOfPixelLayers = {fReader, "PhiToKK_phi_k2_numberOfPixelLayers"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k2_numberOfTrackerLayers = {fReader, "PhiToKK_phi_k2_numberOfTrackerLayers"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k2_numberOfValidHits = {fReader, "PhiToKK_phi_k2_numberOfValidHits"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k2_numberOfValidPixelHits = {fReader, "PhiToKK_phi_k2_numberOfValidPixelHits"};
   TTreeReaderArray<Int_t> PhiToKK_phi_k2_qualityIndex = {fReader, "PhiToKK_phi_k2_qualityIndex"};
   TTreeReaderValue<UInt_t> nMuon = {fReader, "nMuon"};
   TTreeReaderArray<Float_t> Muon_caloCompatibility = {fReader, "Muon_caloCompatibility"};
   TTreeReaderArray<Float_t> Muon_dsaToSlimmedMatching_deltaPtRel = {fReader, "Muon_dsaToSlimmedMatching_deltaPtRel"};
   TTreeReaderArray<Float_t> Muon_dsaToSlimmedMatching_deltaR = {fReader, "Muon_dsaToSlimmedMatching_deltaR"};
   TTreeReaderArray<Float_t> Muon_dsaToSlimmedMatching_deltadxyRel = {fReader, "Muon_dsaToSlimmedMatching_deltadxyRel"};
   TTreeReaderArray<Float_t> Muon_dsaToSlimmedMatching_deltadzRel = {fReader, "Muon_dsaToSlimmedMatching_deltadzRel"};
   TTreeReaderArray<Float_t> Muon_dxy = {fReader, "Muon_dxy"};
   TTreeReaderArray<Float_t> Muon_dxyS = {fReader, "Muon_dxyS"};
   TTreeReaderArray<Float_t> Muon_dxyS_BS = {fReader, "Muon_dxyS_BS"};
   TTreeReaderArray<Float_t> Muon_dxyS_BS_alaRdst = {fReader, "Muon_dxyS_BS_alaRdst"};
   TTreeReaderArray<Float_t> Muon_dxy_BS = {fReader, "Muon_dxy_BS"};
   TTreeReaderArray<Float_t> Muon_dxy_BS_alaRdst = {fReader, "Muon_dxy_BS_alaRdst"};
   TTreeReaderArray<Float_t> Muon_dz = {fReader, "Muon_dz"};
   TTreeReaderArray<Float_t> Muon_dzS = {fReader, "Muon_dzS"};
   TTreeReaderArray<Float_t> Muon_dz_alaRdst = {fReader, "Muon_dz_alaRdst"};
   TTreeReaderArray<Float_t> Muon_eta = {fReader, "Muon_eta"};
   TTreeReaderArray<Float_t> Muon_globalNormalisedChi2 = {fReader, "Muon_globalNormalisedChi2"};
   TTreeReaderArray<Float_t> Muon_kinkFinderChi2 = {fReader, "Muon_kinkFinderChi2"};
   TTreeReaderArray<Float_t> Muon_localPositionChi2 = {fReader, "Muon_localPositionChi2"};
   TTreeReaderArray<Float_t> Muon_mass = {fReader, "Muon_mass"};
   TTreeReaderArray<Float_t> Muon_matched_dpt = {fReader, "Muon_matched_dpt"};
   TTreeReaderArray<Float_t> Muon_matched_dr = {fReader, "Muon_matched_dr"};
   TTreeReaderArray<Float_t> Muon_pfiso03Rel_all = {fReader, "Muon_pfiso03Rel_all"};
   TTreeReaderArray<Float_t> Muon_pfiso03_all = {fReader, "Muon_pfiso03_all"};
   TTreeReaderArray<Float_t> Muon_pfiso04Rel_all = {fReader, "Muon_pfiso04Rel_all"};
   TTreeReaderArray<Float_t> Muon_pfiso04_all = {fReader, "Muon_pfiso04_all"};
   TTreeReaderArray<Float_t> Muon_phi = {fReader, "Muon_phi"};
   TTreeReaderArray<Float_t> Muon_pt = {fReader, "Muon_pt"};
   TTreeReaderArray<Float_t> Muon_segmentCompatibility = {fReader, "Muon_segmentCompatibility"};
   TTreeReaderArray<Float_t> Muon_validHitFraction = {fReader, "Muon_validHitFraction"};
   TTreeReaderArray<Float_t> Muon_vx = {fReader, "Muon_vx"};
   TTreeReaderArray<Float_t> Muon_vy = {fReader, "Muon_vy"};
   TTreeReaderArray<Float_t> Muon_vz = {fReader, "Muon_vz"};
   TTreeReaderArray<Int_t> Muon_charge = {fReader, "Muon_charge"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu10p5_IP3p5 = {fReader, "Muon_fired_HLT_Mu10p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu12_IP6 = {fReader, "Muon_fired_HLT_Mu12_IP6"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu7_IP4 = {fReader, "Muon_fired_HLT_Mu7_IP4"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_IP3 = {fReader, "Muon_fired_HLT_Mu8_IP3"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_IP5 = {fReader, "Muon_fired_HLT_Mu8_IP5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_IP6 = {fReader, "Muon_fired_HLT_Mu8_IP6"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8p5_IP3p5 = {fReader, "Muon_fired_HLT_Mu8p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu9_IP4 = {fReader, "Muon_fired_HLT_Mu9_IP4"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu9_IP5 = {fReader, "Muon_fired_HLT_Mu9_IP5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu9_IP6 = {fReader, "Muon_fired_HLT_Mu9_IP6"};
   TTreeReaderArray<Int_t> Muon_inTimeMuon = {fReader, "Muon_inTimeMuon"};
   TTreeReaderArray<Int_t> Muon_indexMatchedSlimmedMuon = {fReader, "Muon_indexMatchedSlimmedMuon"};
   TTreeReaderArray<Int_t> Muon_isDSAMuon = {fReader, "Muon_isDSAMuon"};
   TTreeReaderArray<Int_t> Muon_isGlobalMuon = {fReader, "Muon_isGlobalMuon"};
   TTreeReaderArray<Int_t> Muon_isGlobalNotTrackerMuon = {fReader, "Muon_isGlobalNotTrackerMuon"};
   TTreeReaderArray<Int_t> Muon_isGlobalOrTrackerMuon = {fReader, "Muon_isGlobalOrTrackerMuon"};
   TTreeReaderArray<Int_t> Muon_isMatchedToSlimmedMuon = {fReader, "Muon_isMatchedToSlimmedMuon"};
   TTreeReaderArray<Int_t> Muon_isPF = {fReader, "Muon_isPF"};
   TTreeReaderArray<Int_t> Muon_isSlimmedMuon = {fReader, "Muon_isSlimmedMuon"};
   TTreeReaderArray<Int_t> Muon_isTrackerMuon = {fReader, "Muon_isTrackerMuon"};
   TTreeReaderArray<Int_t> Muon_isTrackerNotGlobalMuon = {fReader, "Muon_isTrackerNotGlobalMuon"};
   TTreeReaderArray<Int_t> Muon_isTriggering = {fReader, "Muon_isTriggering"};
   TTreeReaderArray<Int_t> Muon_isTriggeringBPark = {fReader, "Muon_isTriggeringBPark"};
   TTreeReaderArray<Int_t> Muon_looseId = {fReader, "Muon_looseId"};
   TTreeReaderArray<Int_t> Muon_mediumId = {fReader, "Muon_mediumId"};
   TTreeReaderArray<Int_t> Muon_numberOfPixelLayers = {fReader, "Muon_numberOfPixelLayers"};
   TTreeReaderArray<Int_t> Muon_numberOfStations = {fReader, "Muon_numberOfStations"};
   TTreeReaderArray<Int_t> Muon_numberOfTrackerLayers = {fReader, "Muon_numberOfTrackerLayers"};
   TTreeReaderArray<Int_t> Muon_numberOfValidMuonHits = {fReader, "Muon_numberOfValidMuonHits"};
   TTreeReaderArray<Int_t> Muon_numberOfValidPixelHits = {fReader, "Muon_numberOfValidPixelHits"};
   TTreeReaderArray<Int_t> Muon_passDSAMuonID = {fReader, "Muon_passDSAMuonID"};
   TTreeReaderArray<Int_t> Muon_pdgId = {fReader, "Muon_pdgId"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu10p5_IP3p5 = {fReader, "Muon_prescale_HLT_Mu10p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu12_IP6 = {fReader, "Muon_prescale_HLT_Mu12_IP6"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu7_IP4 = {fReader, "Muon_prescale_HLT_Mu7_IP4"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_IP3 = {fReader, "Muon_prescale_HLT_Mu8_IP3"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_IP5 = {fReader, "Muon_prescale_HLT_Mu8_IP5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_IP6 = {fReader, "Muon_prescale_HLT_Mu8_IP6"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8p5_IP3p5 = {fReader, "Muon_prescale_HLT_Mu8p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu9_IP4 = {fReader, "Muon_prescale_HLT_Mu9_IP4"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu9_IP5 = {fReader, "Muon_prescale_HLT_Mu9_IP5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu9_IP6 = {fReader, "Muon_prescale_HLT_Mu9_IP6"};
   TTreeReaderArray<Int_t> Muon_softId = {fReader, "Muon_softId"};
   TTreeReaderArray<Int_t> Muon_tightId = {fReader, "Muon_tightId"};
   TTreeReaderArray<Int_t> Muon_trackerHighPurityFlag = {fReader, "Muon_trackerHighPurityFlag"};
   TTreeReaderArray<Int_t> Muon_triggerIdLoose = {fReader, "Muon_triggerIdLoose"};
   TTreeReaderArray<UChar_t> Muon_pfIsoId = {fReader, "Muon_pfIsoId"};
   TTreeReaderArray<UChar_t> Muon_tkIsoId = {fReader, "Muon_tkIsoId"};
   //TTreeReaderValue<UInt_t> nTriggerMuon = {fReader, "nTriggerMuon"};
   //TTreeReaderArray<Float_t> TriggerMuon_eta = {fReader, "TriggerMuon_eta"};
   //TTreeReaderArray<Float_t> TriggerMuon_mass = {fReader, "TriggerMuon_mass"};
   //TTreeReaderArray<Float_t> TriggerMuon_phi = {fReader, "TriggerMuon_phi"};
   //TTreeReaderArray<Float_t> TriggerMuon_pt = {fReader, "TriggerMuon_pt"};
   //TTreeReaderArray<Float_t> TriggerMuon_vx = {fReader, "TriggerMuon_vx"};
   //TTreeReaderArray<Float_t> TriggerMuon_vy = {fReader, "TriggerMuon_vy"};
   //TTreeReaderArray<Float_t> TriggerMuon_vz = {fReader, "TriggerMuon_vz"};
   //TTreeReaderArray<Int_t> TriggerMuon_charge = {fReader, "TriggerMuon_charge"};
   //TTreeReaderArray<Int_t> TriggerMuon_pdgId = {fReader, "TriggerMuon_pdgId"};
   TTreeReaderValue<Float_t> fixedGridRhoFastjetAll = {fReader, "fixedGridRhoFastjetAll"};
   TTreeReaderValue<Float_t> fixedGridRhoFastjetCentral = {fReader, "fixedGridRhoFastjetCentral"};
   TTreeReaderValue<Float_t> fixedGridRhoFastjetCentralCalo = {fReader, "fixedGridRhoFastjetCentralCalo"};
   TTreeReaderValue<Float_t> fixedGridRhoFastjetCentralChargedPileUp = {fReader, "fixedGridRhoFastjetCentralChargedPileUp"};
   TTreeReaderValue<Float_t> fixedGridRhoFastjetCentralNeutral = {fReader, "fixedGridRhoFastjetCentralNeutral"};
   TTreeReaderValue<UChar_t> HLT_Mu7_IP4 = {fReader, "HLT_Mu7_IP4"};
   TTreeReaderValue<UChar_t> HLT_Mu8_IP6 = {fReader, "HLT_Mu8_IP6"};
   TTreeReaderValue<UChar_t> HLT_Mu8_IP5 = {fReader, "HLT_Mu8_IP5"};
   TTreeReaderValue<UChar_t> HLT_Mu8_IP3 = {fReader, "HLT_Mu8_IP3"};
   TTreeReaderValue<UChar_t> HLT_Mu8p5_IP3p5 = {fReader, "HLT_Mu8p5_IP3p5"};
   TTreeReaderValue<UChar_t> HLT_Mu9_IP6 = {fReader, "HLT_Mu9_IP6"};
   TTreeReaderValue<UChar_t> HLT_Mu9_IP5 = {fReader, "HLT_Mu9_IP5"};
   TTreeReaderValue<UChar_t> HLT_Mu9_IP4 = {fReader, "HLT_Mu9_IP4"};
   TTreeReaderValue<UChar_t> HLT_Mu10p5_IP3p5 = {fReader, "HLT_Mu10p5_IP3p5"};
   TTreeReaderValue<UChar_t> HLT_Mu12_IP6 = {fReader, "HLT_Mu12_IP6"};
   TTreeReaderValue<UChar_t> L1_SingleMu7er1p5 = {fReader, "L1_SingleMu7er1p5"};
   TTreeReaderValue<UChar_t> L1_SingleMu8er1p5 = {fReader, "L1_SingleMu8er1p5"};
   TTreeReaderValue<UChar_t> L1_SingleMu9er1p5 = {fReader, "L1_SingleMu9er1p5"};
   TTreeReaderValue<UChar_t> L1_SingleMu10er1p5 = {fReader, "L1_SingleMu10er1p5"};
   TTreeReaderValue<UChar_t> L1_SingleMu12er1p5 = {fReader, "L1_SingleMu12er1p5"};
   TTreeReaderValue<UChar_t> L1_SingleMu22 = {fReader, "L1_SingleMu22"};
   TTreeReaderValue<UInt_t> nTrigObj = {fReader, "nTrigObj"};
   TTreeReaderArray<Float_t> TrigObj_pt = {fReader, "TrigObj_pt"};
   TTreeReaderArray<Float_t> TrigObj_eta = {fReader, "TrigObj_eta"};
   TTreeReaderArray<Float_t> TrigObj_phi = {fReader, "TrigObj_phi"};
   TTreeReaderArray<Float_t> TrigObj_l1pt = {fReader, "TrigObj_l1pt"};
   TTreeReaderArray<Float_t> TrigObj_l1pt_2 = {fReader, "TrigObj_l1pt_2"};
   TTreeReaderArray<Float_t> TrigObj_l2pt = {fReader, "TrigObj_l2pt"};
   TTreeReaderArray<Int_t> TrigObj_id = {fReader, "TrigObj_id"};
   TTreeReaderArray<Int_t> TrigObj_l1iso = {fReader, "TrigObj_l1iso"};
   TTreeReaderArray<Int_t> TrigObj_l1charge = {fReader, "TrigObj_l1charge"};
   TTreeReaderArray<Int_t> TrigObj_filterBits = {fReader, "TrigObj_filterBits"};
   TTreeReaderValue<UInt_t> nOtherPV = {fReader, "nOtherPV"};
   TTreeReaderArray<Float_t> OtherPV_z = {fReader, "OtherPV_z"};
   TTreeReaderValue<Float_t> PV_ndof = {fReader, "PV_ndof"};
   TTreeReaderValue<Float_t> PV_x = {fReader, "PV_x"};
   TTreeReaderValue<Float_t> PV_y = {fReader, "PV_y"};
   TTreeReaderValue<Float_t> PV_z = {fReader, "PV_z"};
   TTreeReaderValue<Float_t> PV_chi2 = {fReader, "PV_chi2"};
   TTreeReaderValue<Float_t> PV_score = {fReader, "PV_score"};
   TTreeReaderValue<Int_t> PV_npvs = {fReader, "PV_npvs"};
   TTreeReaderValue<Int_t> PV_npvsGood = {fReader, "PV_npvsGood"};
   TTreeReaderValue<UInt_t> nSV = {fReader, "nSV"};
   TTreeReaderArray<Float_t> SV_dlen = {fReader, "SV_dlen"};
   TTreeReaderArray<Float_t> SV_dlenSig = {fReader, "SV_dlenSig"};
   TTreeReaderArray<Float_t> SV_pAngle = {fReader, "SV_pAngle"};
   TTreeReaderArray<Float_t> SV_chi2 = {fReader, "SV_chi2"};
   TTreeReaderArray<Float_t> SV_eta = {fReader, "SV_eta"};
   TTreeReaderArray<Float_t> SV_mass = {fReader, "SV_mass"};
   TTreeReaderArray<Float_t> SV_ndof = {fReader, "SV_ndof"};
   TTreeReaderArray<Float_t> SV_phi = {fReader, "SV_phi"};
   TTreeReaderArray<Float_t> SV_pt = {fReader, "SV_pt"};
   TTreeReaderArray<Float_t> SV_x = {fReader, "SV_x"};
   TTreeReaderArray<Float_t> SV_y = {fReader, "SV_y"};
   TTreeReaderArray<Float_t> SV_z = {fReader, "SV_z"};

   // branches intentionally pointed to wrong variable, to avoid crash when running on data
   TTreeReaderValue<UInt_t> nGenPart = {fReader, "nBsToPhiPhiTo4K"};
   TTreeReaderArray<Float_t> GenPart_eta = {fReader, "BsToPhiPhiTo4K_Bs_eta"};
   TTreeReaderArray<Float_t> GenPart_mass = {fReader, "BsToPhiPhiTo4K_Bs_eta"};
   TTreeReaderArray<Float_t> GenPart_phi = {fReader, "BsToPhiPhiTo4K_Bs_eta"};
   TTreeReaderArray<Float_t> GenPart_pt = {fReader, "BsToPhiPhiTo4K_Bs_eta"};
   TTreeReaderArray<Float_t> GenPart_vx = {fReader, "BsToPhiPhiTo4K_Bs_eta"};
   TTreeReaderArray<Float_t> GenPart_vy = {fReader, "BsToPhiPhiTo4K_Bs_eta"};
   TTreeReaderArray<Float_t> GenPart_vz = {fReader, "BsToPhiPhiTo4K_Bs_eta"};
   TTreeReaderArray<Int_t> GenPart_genPartIdxMother = {fReader, "BsToPhiPhiTo4K_k1_idx"};
   TTreeReaderArray<Int_t> GenPart_pdgId = {fReader, "BsToPhiPhiTo4K_k1_idx"};
   TTreeReaderArray<Int_t> GenPart_status = {fReader, "BsToPhiPhiTo4K_k1_idx"};
   TTreeReaderArray<Int_t> GenPart_statusFlags = {fReader, "BsToPhiPhiTo4K_k1_idx"};
   TTreeReaderArray<Int_t> Muon_genPartIdx = {fReader, "BsToPhiPhiTo4K_k1_idx"};
   // the correct value will be filled for MC only
   TTreeReaderValue<Int_t> Pileup_nPU = {fReader, "PV_npvs"};
   TTreeReaderValue<Float_t> Pileup_nTrueInt = {fReader, "PV_ndof"};

   BsToPhiPhiTo4KDumper(TTree * /*tree*/ =0) { }
   virtual ~BsToPhiPhiTo4KDumper() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();
   
   // output file
   TFile* my_file;  

   // some branches will be added only if sample is MC
   Bool_t isMC;
   Bool_t isSignalMC;

   // trees to fill
   TTree* signal_tree;
   
   // filling variables
   ULong64_t the_event = -99;
   Int_t the_run = -99;
   Int_t the_lumi = -99;
   Int_t the_pu_ntrueint = -99;
   Int_t the_pv_npvs = -99;

   Int_t the_ntriggermuon = -99;

   Float_t the_bs_beta = -99.;
   Float_t the_bs_gamma = -99.;
   Float_t the_bs_eta = -99.;
   Float_t the_bs_charge  = -99.;
   Float_t the_bs_cos2d = -99.;
   Float_t the_bs_cxx = -99.;
   Float_t the_bs_cyx = -99.;
   Float_t the_bs_cyy = -99.;
   Float_t the_bs_czx = -99.;
   Float_t the_bs_czy = -99.;
   Float_t the_bs_czz = -99.;
   Float_t the_bs_lxy = -99.;
   Float_t the_bs_lxyerr = -99.;
   Float_t the_bs_lxysig = -99.;
   Float_t the_bs_lxy_corr = -99.;
   Float_t the_bs_lxyerr_corr = -99.;
   Float_t the_bs_lxysig_corr = -99.;
   Float_t the_bs_lxyz = -99.;
   Float_t the_bs_lxyzerr = -99.;
   Float_t the_bs_lxyz_corr = -99.;
   Float_t the_bs_lxyzerr_corr = -99.;
   Float_t the_bs_lx = -99.;
   Float_t the_bs_ly = -99.;
   Float_t the_bs_lx_posbsz = -99.;
   Float_t the_bs_ly_posbsz = -99.;
   Float_t the_bs_lx_posbspv = -99.;
   Float_t the_bs_ly_posbspv = -99.;
   Float_t the_bs_ct_2d_cm = -99.;
   Float_t the_bs_ct_2d_cm_posbsz = -99.;
   Float_t the_bs_ct_2d_cm_posbspv = -99.;
   Float_t the_bs_ct_3d_cm = -99.;
   Float_t the_bs_mass = -99.;
   Float_t the_bs_mass_err = -99.;
   Float_t the_bs_mass_corr = -99.;
   Float_t the_bs_phi = -99.;
   Float_t the_bs_pt = -99.;
   Float_t the_bs_pterr = -99.;
   Float_t the_bs_px = -99.;
   Float_t the_bs_py = -99.;
   Float_t the_bs_sv_prob = -99.;
   Float_t the_bs_sv_chi2 = -99.;
   Float_t the_bs_sv_ndof = -99.;
   Float_t the_bs_vx = -99.;
   Float_t the_bs_vy = -99.;
   Float_t the_bs_vz = -99.;

   Float_t the_beamspot_x = -99.;
   Float_t the_beamspot_y = -99.;
   Float_t the_beamspot_z = -99.;

   Float_t the_cos_theta_star_phi1 = -99.;
   Float_t the_cos_theta_star_phi2 = -99.;

   Float_t the_deltar_k1k2 = -99.;
   Float_t the_deltar_k1k3 = -99.;
   Float_t the_deltar_k1k4 = -99.;
   Float_t the_deltar_k2k3 = -99.;
   Float_t the_deltar_k2k4 = -99.;
   Float_t the_deltar_k3k4 = -99.;
   Float_t the_deltar_min = -99.;
   Float_t the_deltar_max = -99.;
   Float_t the_deltar_phi1phi2 = -99.;

   Float_t the_k1k3_mass = -99.;
   Float_t the_k1k3_pt = -99.;
   Float_t the_k1k4_mass = -99.;
   Float_t the_k1k4_pt = -99.;
   Float_t the_k2k3_mass = -99.;
   Float_t the_k2k3_pt = -99.;
   Float_t the_k2k4_mass = -99.;
   Float_t the_k2k4_pt = -99.;

   Int_t phi1_idx = -99;
   Int_t phi2_idx = -99;
   Int_t k1_idx = -99;
   Int_t k2_idx = -99;
   Int_t k3_idx = -99;
   Int_t k4_idx = -99;

   Int_t ismatched = -99;

   Float_t the_phi1_eta = -99.;
   Float_t the_phi1_phi = -99.;
   Float_t the_phi1_pt = -99.;
   Float_t the_phi1_mass = -99.;
   Float_t the_phi1_cos_theta_star_k1 = -99.;
   Float_t the_phi1_cos_theta_star_k2 = -99.;
   Float_t the_phi1_deltar = -99.;
   Float_t the_phi1_cos2d = -99.;
   Float_t the_phi1_cxx = -99.;
   Float_t the_phi1_cyx = -99.;
   Float_t the_phi1_cyy = -99.;
   Float_t the_phi1_czx = -99.;
   Float_t the_phi1_czy = -99.;
   Float_t the_phi1_czz = -99.;
   Float_t the_phi1_sv_chi2 = -99.;
   Float_t the_phi1_sv_ndof = -99.;
   Float_t the_phi1_sv_prob = -99.;
   Float_t the_phi1_lxy = -99.;
   Float_t the_phi1_lxysig = -99.;
   Float_t the_phi1_mass_err = -99.;
   Float_t the_phi1_vx = -99.;
   Float_t the_phi1_vy = -99.;
   Float_t the_phi1_vz = -99.;
   Int_t the_phi1_is_matched = -99;

   Float_t the_phi2_eta = -99.;
   Float_t the_phi2_phi = -99.;
   Float_t the_phi2_pt = -99.;
   Float_t the_phi2_mass = -99.;
   Float_t the_phi2_cos_theta_star_k1 = -99.;
   Float_t the_phi2_cos_theta_star_k2 = -99.;
   Float_t the_phi2_deltar = -99.;
   Float_t the_phi2_cos2d = -99.;
   Float_t the_phi2_cxx = -99.;
   Float_t the_phi2_cyx = -99.;
   Float_t the_phi2_cyy = -99.;
   Float_t the_phi2_czx = -99.;
   Float_t the_phi2_czy = -99.;
   Float_t the_phi2_czz = -99.;
   Float_t the_phi2_sv_chi2 = -99.;
   Float_t the_phi2_sv_ndof = -99.;
   Float_t the_phi2_sv_prob = -99.;
   Float_t the_phi2_lxy = -99.;
   Float_t the_phi2_lxysig = -99.;
   Float_t the_phi2_mass_err = -99.;
   Float_t the_phi2_vx = -99.;
   Float_t the_phi2_vy = -99.;
   Float_t the_phi2_vz = -99.;
   Int_t the_phi2_is_matched = -99;

   Float_t the_phi1_pt_times_phi2_pt = -99.;

   Float_t the_k1_eta = -99.;
   Float_t the_k1_phi = -99.;
   Float_t the_k1_pt = -99.;
   Float_t the_k1_dcasig = -99.;
   Float_t the_k1_charge = -99.;
   Float_t the_k1_covlamlam = -99.;
   Float_t the_k1_covlamphi = -99.;
   Float_t the_k1_covphiphi = -99.;
   Float_t the_k1_covqoplam = -99.;
   Float_t the_k1_covqopphi = -99.;
   Float_t the_k1_covqopqop = -99.;
   Float_t the_k1_dxy = -99.;
   Float_t the_k1_dxysig = -99.;
   Float_t the_k1_dz = -99.;
   Float_t the_k1_dzsig = -99.;
   Float_t the_k1_normalisedchi2 = -99.;
   Float_t the_k1_pt_err = -99.;
   Float_t the_k1_validfraction = -99.;
   Float_t the_k1_vx = -99.;
   Float_t the_k1_vy = -99.;
   Float_t the_k1_vz = -99.;
   Int_t the_k1_highpurityflag = -99;
   Int_t the_k1_islost = -99;
   Int_t the_k1_ispacked = -99;
   Int_t the_k1_ndof = -99;
   Int_t the_k1_numberlosthits = -99;
   Int_t the_k1_numberpixelhits = -99;
   Int_t the_k1_numberpixellayers = -99;
   Int_t the_k1_numbertrackerlayers = -99;
   Int_t the_k1_numberofvalidhits = -99;
   Int_t the_k1_numberofvalidpixelhits = -99;
   Int_t the_k1_qualityindex = -99;

   Float_t the_k2_eta = -99.;
   Float_t the_k2_phi = -99.;
   Float_t the_k2_pt = -99.;
   Float_t the_k2_dcasig = -99.;
   Float_t the_k2_charge = -99.;
   Float_t the_k2_covlamlam = -99.;
   Float_t the_k2_covlamphi = -99.;
   Float_t the_k2_covphiphi = -99.;
   Float_t the_k2_covqoplam = -99.;
   Float_t the_k2_covqopphi = -99.;
   Float_t the_k2_covqopqop = -99.;
   Float_t the_k2_dxy = -99.;
   Float_t the_k2_dxysig = -99.;
   Float_t the_k2_dz = -99.;
   Float_t the_k2_dzsig = -99.;
   Float_t the_k2_normalisedchi2 = -99.;
   Float_t the_k2_pt_err = -99.;
   Float_t the_k2_validfraction = -99.;
   Float_t the_k2_vx = -99.;
   Float_t the_k2_vy = -99.;
   Float_t the_k2_vz = -99.;
   Int_t the_k2_highpurityflag = -99;
   Int_t the_k2_islost = -99;
   Int_t the_k2_ispacked = -99;
   Int_t the_k2_ndof = -99;
   Int_t the_k2_numberlosthits = -99;
   Int_t the_k2_numberpixelhits = -99;
   Int_t the_k2_numberpixellayers = -99;
   Int_t the_k2_numbertrackerlayers = -99;
   Int_t the_k2_numberofvalidhits = -99;
   Int_t the_k2_numberofvalidpixelhits = -99;
   Int_t the_k2_qualityindex = -99;

   Float_t the_k3_eta = -99.;
   Float_t the_k3_phi = -99.;
   Float_t the_k3_pt = -99.;
   Float_t the_k3_dcasig = -99.;
   Float_t the_k3_charge = -99.;
   Float_t the_k3_covlamlam = -99.;
   Float_t the_k3_covlamphi = -99.;
   Float_t the_k3_covphiphi = -99.;
   Float_t the_k3_covqoplam = -99.;
   Float_t the_k3_covqopphi = -99.;
   Float_t the_k3_covqopqop = -99.;
   Float_t the_k3_dxy = -99.;
   Float_t the_k3_dxysig = -99.;
   Float_t the_k3_dz = -99.;
   Float_t the_k3_dzsig = -99.;
   Float_t the_k3_normalisedchi2 = -99.;
   Float_t the_k3_pt_err = -99.;
   Float_t the_k3_validfraction = -99.;
   Float_t the_k3_vx = -99.;
   Float_t the_k3_vy = -99.;
   Float_t the_k3_vz = -99.;
   Int_t the_k3_highpurityflag = -99;
   Int_t the_k3_islost = -99;
   Int_t the_k3_ispacked = -99;
   Int_t the_k3_ndof = -99;
   Int_t the_k3_numberlosthits = -99;
   Int_t the_k3_numberpixelhits = -99;
   Int_t the_k3_numberpixellayers = -99;
   Int_t the_k3_numbertrackerlayers = -99;
   Int_t the_k3_numberofvalidhits = -99;
   Int_t the_k3_numberofvalidpixelhits = -99;
   Int_t the_k3_qualityindex = -99;

   Float_t the_k4_eta = -99.;
   Float_t the_k4_phi = -99.;
   Float_t the_k4_pt = -99.;
   Float_t the_k4_dcasig = -99.;
   Float_t the_k4_charge = -99.;
   Float_t the_k4_covlamlam = -99.;
   Float_t the_k4_covlamphi = -99.;
   Float_t the_k4_covphiphi = -99.;
   Float_t the_k4_covqoplam = -99.;
   Float_t the_k4_covqopphi = -99.;
   Float_t the_k4_covqopqop = -99.;
   Float_t the_k4_dxy = -99.;
   Float_t the_k4_dxysig = -99.;
   Float_t the_k4_dz = -99.;
   Float_t the_k4_dzsig = -99.;
   Float_t the_k4_normalisedchi2 = -99.;
   Float_t the_k4_pt_err = -99.;
   Float_t the_k4_validfraction = -99.;
   Float_t the_k4_vx = -99.;
   Float_t the_k4_vy = -99.;
   Float_t the_k4_vz = -99.;
   Int_t the_k4_highpurityflag = -99;
   Int_t the_k4_islost = -99;
   Int_t the_k4_ispacked = -99;
   Int_t the_k4_ndof = -99;
   Int_t the_k4_numberlosthits = -99;
   Int_t the_k4_numberpixelhits = -99;
   Int_t the_k4_numberpixellayers = -99;
   Int_t the_k4_numbertrackerlayers = -99;
   Int_t the_k4_numberofvalidhits = -99;
   Int_t the_k4_numberofvalidpixelhits = -99;
   Int_t the_k4_qualityindex = -99;

   ClassDef(BsToPhiPhiTo4KDumper,0);

};

#endif

#ifdef BsToPhiPhiTo4KDumper_cxx
void BsToPhiPhiTo4KDumper::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t BsToPhiPhiTo4KDumper::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef BsToPhiPhiTo4KDumper_cxx
