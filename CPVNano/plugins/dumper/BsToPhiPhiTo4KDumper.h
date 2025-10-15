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
   TTreeReaderValue<UInt_t> bunchCrossing = {fReader, "bunchCrossing"};
   TTreeReaderValue<UInt_t> orbitNumber = {fReader, "orbitNumber"};
   TTreeReaderValue<Int_t> nBsToPhiPhiTo4K = {fReader, "nBsToPhiPhiTo4K"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_genKaon1GrandGrandMother_genIdx = {fReader, "BsToPhiPhiTo4K_genKaon1GrandGrandMother_genIdx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_genKaon1GrandMother_genIdx = {fReader, "BsToPhiPhiTo4K_genKaon1GrandMother_genIdx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_genKaon1Mother_genIdx = {fReader, "BsToPhiPhiTo4K_genKaon1Mother_genIdx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_genKaon2GrandGrandMother_genIdx = {fReader, "BsToPhiPhiTo4K_genKaon2GrandGrandMother_genIdx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_genKaon2GrandMother_genIdx = {fReader, "BsToPhiPhiTo4K_genKaon2GrandMother_genIdx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_genKaon2Mother_genIdx = {fReader, "BsToPhiPhiTo4K_genKaon2Mother_genIdx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_genKaon3GrandGrandMother_genIdx = {fReader, "BsToPhiPhiTo4K_genKaon3GrandGrandMother_genIdx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_genKaon3GrandMother_genIdx = {fReader, "BsToPhiPhiTo4K_genKaon3GrandMother_genIdx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_genKaon3Mother_genIdx = {fReader, "BsToPhiPhiTo4K_genKaon3Mother_genIdx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_genKaon4GrandGrandMother_genIdx = {fReader, "BsToPhiPhiTo4K_genKaon4GrandGrandMother_genIdx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_genKaon4GrandMother_genIdx = {fReader, "BsToPhiPhiTo4K_genKaon4GrandMother_genIdx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_genKaon4Mother_genIdx = {fReader, "BsToPhiPhiTo4K_genKaon4Mother_genIdx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_isMatched = {fReader, "BsToPhiPhiTo4K_isMatched"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_k1_genIdx = {fReader, "BsToPhiPhiTo4K_k1_genIdx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_k1_idx = {fReader, "BsToPhiPhiTo4K_k1_idx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_k2_genIdx = {fReader, "BsToPhiPhiTo4K_k2_genIdx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_k2_idx = {fReader, "BsToPhiPhiTo4K_k2_idx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_k3_genIdx = {fReader, "BsToPhiPhiTo4K_k3_genIdx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_k3_idx = {fReader, "BsToPhiPhiTo4K_k3_idx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_k4_genIdx = {fReader, "BsToPhiPhiTo4K_k4_genIdx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_k4_idx = {fReader, "BsToPhiPhiTo4K_k4_idx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_phi1_idx = {fReader, "BsToPhiPhiTo4K_phi1_idx"};
   TTreeReaderArray<Int_t> BsToPhiPhiTo4K_phi2_idx = {fReader, "BsToPhiPhiTo4K_phi2_idx"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_beta = {fReader, "BsToPhiPhiTo4K_Bs_beta"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_charge = {fReader, "BsToPhiPhiTo4K_Bs_charge"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_cos2D = {fReader, "BsToPhiPhiTo4K_Bs_cos2D"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_ct_2D_cm = {fReader, "BsToPhiPhiTo4K_Bs_ct_2D_cm"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_ct_2D_cm_posbspv = {fReader, "BsToPhiPhiTo4K_Bs_ct_2D_cm_posbspv"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_ct_2D_cm_posbsz = {fReader, "BsToPhiPhiTo4K_Bs_ct_2D_cm_posbsz"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_ct_2D_cm_posthepv = {fReader, "BsToPhiPhiTo4K_Bs_ct_2D_cm_posthepv"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_cxx = {fReader, "BsToPhiPhiTo4K_Bs_cxx"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_cyx = {fReader, "BsToPhiPhiTo4K_Bs_cyx"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_cyy = {fReader, "BsToPhiPhiTo4K_Bs_cyy"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_czx = {fReader, "BsToPhiPhiTo4K_Bs_czx"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_czy = {fReader, "BsToPhiPhiTo4K_Bs_czy"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_czz = {fReader, "BsToPhiPhiTo4K_Bs_czz"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_eta = {fReader, "BsToPhiPhiTo4K_Bs_eta"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_gamma = {fReader, "BsToPhiPhiTo4K_Bs_gamma"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lx = {fReader, "BsToPhiPhiTo4K_Bs_lx"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lx_posbspv = {fReader, "BsToPhiPhiTo4K_Bs_lx_posbspv"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lx_posbsz = {fReader, "BsToPhiPhiTo4K_Bs_lx_posbsz"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lx_posthepv = {fReader, "BsToPhiPhiTo4K_Bs_lx_posthepv"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxy = {fReader, "BsToPhiPhiTo4K_Bs_lxy"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxyErr = {fReader, "BsToPhiPhiTo4K_Bs_lxyErr"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxyErr_corr = {fReader, "BsToPhiPhiTo4K_Bs_lxyErr_corr"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxy_corr = {fReader, "BsToPhiPhiTo4K_Bs_lxy_corr"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxy_sig = {fReader, "BsToPhiPhiTo4K_Bs_lxy_sig"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxy_sig_corr = {fReader, "BsToPhiPhiTo4K_Bs_lxy_sig_corr"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxyz = {fReader, "BsToPhiPhiTo4K_Bs_lxyz"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxyzErr = {fReader, "BsToPhiPhiTo4K_Bs_lxyzErr"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxyzErr_corr = {fReader, "BsToPhiPhiTo4K_Bs_lxyzErr_corr"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_lxyz_corr = {fReader, "BsToPhiPhiTo4K_Bs_lxyz_corr"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_ly = {fReader, "BsToPhiPhiTo4K_Bs_ly"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_ly_posbspv = {fReader, "BsToPhiPhiTo4K_Bs_ly_posbspv"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_ly_posbsz = {fReader, "BsToPhiPhiTo4K_Bs_ly_posbsz"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_Bs_ly_posthepv = {fReader, "BsToPhiPhiTo4K_Bs_ly_posthepv"};
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
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_cos_theta_k1 = {fReader, "BsToPhiPhiTo4K_cos_theta_k1"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_cos_theta_k3 = {fReader, "BsToPhiPhiTo4K_cos_theta_k3"};
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
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_gen_Bs_ct = {fReader, "BsToPhiPhiTo4K_gen_Bs_ct"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1_dcasig_pv = {fReader, "BsToPhiPhiTo4K_k1_dcasig_pv"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1_dxy = {fReader, "BsToPhiPhiTo4K_k1_dxy"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1_dxyS = {fReader, "BsToPhiPhiTo4K_k1_dxyS"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1_dz = {fReader, "BsToPhiPhiTo4K_k1_dz"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1_dzS = {fReader, "BsToPhiPhiTo4K_k1_dzS"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1_eta = {fReader, "BsToPhiPhiTo4K_k1_eta"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1_mass = {fReader, "BsToPhiPhiTo4K_k1_mass"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1_phi = {fReader, "BsToPhiPhiTo4K_k1_phi"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1_pt = {fReader, "BsToPhiPhiTo4K_k1_pt"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1k3_mass = {fReader, "BsToPhiPhiTo4K_k1k3_mass"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1k3_pt = {fReader, "BsToPhiPhiTo4K_k1k3_pt"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1k4_mass = {fReader, "BsToPhiPhiTo4K_k1k4_mass"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k1k4_pt = {fReader, "BsToPhiPhiTo4K_k1k4_pt"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2_dcasig_pv = {fReader, "BsToPhiPhiTo4K_k2_dcasig_pv"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2_dxy = {fReader, "BsToPhiPhiTo4K_k2_dxy"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2_dxyS = {fReader, "BsToPhiPhiTo4K_k2_dxyS"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2_dz = {fReader, "BsToPhiPhiTo4K_k2_dz"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2_dzS = {fReader, "BsToPhiPhiTo4K_k2_dzS"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2_eta = {fReader, "BsToPhiPhiTo4K_k2_eta"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2_mass = {fReader, "BsToPhiPhiTo4K_k2_mass"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2_phi = {fReader, "BsToPhiPhiTo4K_k2_phi"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2_pt = {fReader, "BsToPhiPhiTo4K_k2_pt"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2k3_mass = {fReader, "BsToPhiPhiTo4K_k2k3_mass"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2k3_pt = {fReader, "BsToPhiPhiTo4K_k2k3_pt"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2k4_mass = {fReader, "BsToPhiPhiTo4K_k2k4_mass"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k2k4_pt = {fReader, "BsToPhiPhiTo4K_k2k4_pt"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k3_dcasig_pv = {fReader, "BsToPhiPhiTo4K_k3_dcasig_pv"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k3_dxy = {fReader, "BsToPhiPhiTo4K_k3_dxy"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k3_dxyS = {fReader, "BsToPhiPhiTo4K_k3_dxyS"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k3_dz = {fReader, "BsToPhiPhiTo4K_k3_dz"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k3_dzS = {fReader, "BsToPhiPhiTo4K_k3_dzS"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k3_eta = {fReader, "BsToPhiPhiTo4K_k3_eta"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k3_mass = {fReader, "BsToPhiPhiTo4K_k3_mass"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k3_phi = {fReader, "BsToPhiPhiTo4K_k3_phi"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k3_pt = {fReader, "BsToPhiPhiTo4K_k3_pt"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k4_dcasig_pv = {fReader, "BsToPhiPhiTo4K_k4_dcasig_pv"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k4_dxy = {fReader, "BsToPhiPhiTo4K_k4_dxy"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k4_dxyS = {fReader, "BsToPhiPhiTo4K_k4_dxyS"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k4_dz = {fReader, "BsToPhiPhiTo4K_k4_dz"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_k4_dzS = {fReader, "BsToPhiPhiTo4K_k4_dzS"};
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
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_phi_star = {fReader, "BsToPhiPhiTo4K_phi_star"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_the_PV_chi2 = {fReader, "BsToPhiPhiTo4K_the_PV_chi2"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_the_PV_covXX = {fReader, "BsToPhiPhiTo4K_the_PV_covXX"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_the_PV_covXY = {fReader, "BsToPhiPhiTo4K_the_PV_covXY"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_the_PV_covXZ = {fReader, "BsToPhiPhiTo4K_the_PV_covXZ"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_the_PV_covYY = {fReader, "BsToPhiPhiTo4K_the_PV_covYY"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_the_PV_covYZ = {fReader, "BsToPhiPhiTo4K_the_PV_covYZ"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_the_PV_covZZ = {fReader, "BsToPhiPhiTo4K_the_PV_covZZ"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_the_PV_ndof = {fReader, "BsToPhiPhiTo4K_the_PV_ndof"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_the_PV_x = {fReader, "BsToPhiPhiTo4K_the_PV_x"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_the_PV_y = {fReader, "BsToPhiPhiTo4K_the_PV_y"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_the_PV_z = {fReader, "BsToPhiPhiTo4K_the_PV_z"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_trgmu_dxy = {fReader, "BsToPhiPhiTo4K_trgmu_dxy"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_trgmu_dxyS = {fReader, "BsToPhiPhiTo4K_trgmu_dxyS"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_trgmu_dz = {fReader, "BsToPhiPhiTo4K_trgmu_dz"};
   TTreeReaderArray<Float_t> BsToPhiPhiTo4K_trgmu_dzS = {fReader, "BsToPhiPhiTo4K_trgmu_dzS"};
   TTreeReaderValue<Int_t> nPhiToKK = {fReader, "nPhiToKK"};
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
   TTreeReaderArray<Float_t> PhiToKK_cos_theta_star_k1 = {fReader, "PhiToKK_cos_theta_star_k1"};
   TTreeReaderArray<Float_t> PhiToKK_cos_theta_star_k2 = {fReader, "PhiToKK_cos_theta_star_k2"};
   TTreeReaderArray<Float_t> PhiToKK_deltaR_postfit = {fReader, "PhiToKK_deltaR_postfit"};
   TTreeReaderArray<Float_t> PhiToKK_deltaR_prefit = {fReader, "PhiToKK_deltaR_prefit"};
   TTreeReaderArray<Float_t> PhiToKK_eta = {fReader, "PhiToKK_eta"};
   TTreeReaderArray<Float_t> PhiToKK_k1_pt = {fReader, "PhiToKK_k1_pt"};
   TTreeReaderArray<Float_t> PhiToKK_k2_pt = {fReader, "PhiToKK_k2_pt"};
   TTreeReaderArray<Float_t> PhiToKK_mass = {fReader, "PhiToKK_mass"};
   TTreeReaderArray<Float_t> PhiToKK_phi = {fReader, "PhiToKK_phi"};
   TTreeReaderArray<Float_t> PhiToKK_phi_charge = {fReader, "PhiToKK_phi_charge"};
   TTreeReaderArray<Float_t> PhiToKK_phi_cos2D = {fReader, "PhiToKK_phi_cos2D"};
   TTreeReaderArray<Float_t> PhiToKK_phi_ct_2D_cm = {fReader, "PhiToKK_phi_ct_2D_cm"};
   TTreeReaderArray<Float_t> PhiToKK_phi_ct_2D_cm_posbspv = {fReader, "PhiToKK_phi_ct_2D_cm_posbspv"};
   TTreeReaderArray<Float_t> PhiToKK_phi_ct_2D_cm_posbsz = {fReader, "PhiToKK_phi_ct_2D_cm_posbsz"};
   TTreeReaderArray<Float_t> PhiToKK_phi_cxx = {fReader, "PhiToKK_phi_cxx"};
   TTreeReaderArray<Float_t> PhiToKK_phi_cyx = {fReader, "PhiToKK_phi_cyx"};
   TTreeReaderArray<Float_t> PhiToKK_phi_cyy = {fReader, "PhiToKK_phi_cyy"};
   TTreeReaderArray<Float_t> PhiToKK_phi_czx = {fReader, "PhiToKK_phi_czx"};
   TTreeReaderArray<Float_t> PhiToKK_phi_czy = {fReader, "PhiToKK_phi_czy"};
   TTreeReaderArray<Float_t> PhiToKK_phi_czz = {fReader, "PhiToKK_phi_czz"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_DCASig = {fReader, "PhiToKK_phi_k1_DCASig"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_charge = {fReader, "PhiToKK_phi_k1_charge"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_chi2 = {fReader, "PhiToKK_phi_k1_chi2"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_covLamLam = {fReader, "PhiToKK_phi_k1_covLamLam"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_covLamPhi = {fReader, "PhiToKK_phi_k1_covLamPhi"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_covPhiPhi = {fReader, "PhiToKK_phi_k1_covPhiPhi"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_covQopLam = {fReader, "PhiToKK_phi_k1_covQopLam"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_covQopPhi = {fReader, "PhiToKK_phi_k1_covQopPhi"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_covQopQop = {fReader, "PhiToKK_phi_k1_covQopQop"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_drTrg = {fReader, "PhiToKK_phi_k1_drTrg"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_dxy = {fReader, "PhiToKK_phi_k1_dxy"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_dxyS = {fReader, "PhiToKK_phi_k1_dxyS"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_dz = {fReader, "PhiToKK_phi_k1_dz"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_dzS = {fReader, "PhiToKK_phi_k1_dzS"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_dzTrg = {fReader, "PhiToKK_phi_k1_dzTrg"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_eta = {fReader, "PhiToKK_phi_k1_eta"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_mass = {fReader, "PhiToKK_phi_k1_mass"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_normalisedChi2 = {fReader, "PhiToKK_phi_k1_normalisedChi2"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_pdgid = {fReader, "PhiToKK_phi_k1_pdgid"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k1_phi = {fReader, "PhiToKK_phi_k1_phi"};
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
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_drTrg = {fReader, "PhiToKK_phi_k2_drTrg"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_dxy = {fReader, "PhiToKK_phi_k2_dxy"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_dxyS = {fReader, "PhiToKK_phi_k2_dxyS"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_dz = {fReader, "PhiToKK_phi_k2_dz"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_dzS = {fReader, "PhiToKK_phi_k2_dzS"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_dzTrg = {fReader, "PhiToKK_phi_k2_dzTrg"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_eta = {fReader, "PhiToKK_phi_k2_eta"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_mass = {fReader, "PhiToKK_phi_k2_mass"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_normalisedChi2 = {fReader, "PhiToKK_phi_k2_normalisedChi2"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_pdgid = {fReader, "PhiToKK_phi_k2_pdgid"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_phi = {fReader, "PhiToKK_phi_k2_phi"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_ptErr = {fReader, "PhiToKK_phi_k2_ptErr"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_validFraction = {fReader, "PhiToKK_phi_k2_validFraction"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_vx = {fReader, "PhiToKK_phi_k2_vx"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_vy = {fReader, "PhiToKK_phi_k2_vy"};
   TTreeReaderArray<Float_t> PhiToKK_phi_k2_vz = {fReader, "PhiToKK_phi_k2_vz"};
   TTreeReaderArray<Float_t> PhiToKK_phi_lxy = {fReader, "PhiToKK_phi_lxy"};
   TTreeReaderArray<Float_t> PhiToKK_phi_lxy_sig = {fReader, "PhiToKK_phi_lxy_sig"};
   TTreeReaderArray<Float_t> PhiToKK_phi_masserr = {fReader, "PhiToKK_phi_masserr"};
   TTreeReaderArray<Float_t> PhiToKK_phi_vx = {fReader, "PhiToKK_phi_vx"};
   TTreeReaderArray<Float_t> PhiToKK_phi_vy = {fReader, "PhiToKK_phi_vy"};
   TTreeReaderArray<Float_t> PhiToKK_phi_vz = {fReader, "PhiToKK_phi_vz"};
   TTreeReaderArray<Float_t> PhiToKK_pt = {fReader, "PhiToKK_pt"};
   TTreeReaderArray<Float_t> PhiToKK_sv_chi2 = {fReader, "PhiToKK_sv_chi2"};
   TTreeReaderArray<Float_t> PhiToKK_sv_ndof = {fReader, "PhiToKK_sv_ndof"};
   TTreeReaderArray<Float_t> PhiToKK_sv_prob = {fReader, "PhiToKK_sv_prob"};
   TTreeReaderValue<Short_t> BeamSpot_type = {fReader, "BeamSpot_type"};
   TTreeReaderValue<Float_t> BeamSpot_sigmaZ = {fReader, "BeamSpot_sigmaZ"};
   TTreeReaderValue<Float_t> BeamSpot_sigmaZError = {fReader, "BeamSpot_sigmaZError"};
   TTreeReaderValue<Float_t> BeamSpot_z = {fReader, "BeamSpot_z"};
   TTreeReaderValue<Float_t> BeamSpot_zError = {fReader, "BeamSpot_zError"};
   TTreeReaderValue<Int_t> nMuon = {fReader, "nMuon"};
   TTreeReaderArray<Bool_t> Muon_inTimeMuon = {fReader, "Muon_inTimeMuon"};
   TTreeReaderArray<Bool_t> Muon_mediumID = {fReader, "Muon_mediumID"};
   TTreeReaderArray<Bool_t> Muon_softID = {fReader, "Muon_softID"};
   TTreeReaderArray<Bool_t> Muon_tightID = {fReader, "Muon_tightID"};
   TTreeReaderArray<Int_t> Muon_charge = {fReader, "Muon_charge"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu10_Barrel_L1HP11_IP6 = {fReader, "Muon_fired_HLT_Mu10_Barrel_L1HP11_IP6"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu10p5_IP3p5 = {fReader, "Muon_fired_HLT_Mu10p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu12_IP6 = {fReader, "Muon_fired_HLT_Mu12_IP6"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu6_Barrel_L1HP7_IP6 = {fReader, "Muon_fired_HLT_Mu6_Barrel_L1HP7_IP6"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu7_Barrel_L1HP8_IP6 = {fReader, "Muon_fired_HLT_Mu7_Barrel_L1HP8_IP6"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu7_IP4 = {fReader, "Muon_fired_HLT_Mu7_IP4"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_Barrel_L1HP9_IP6 = {fReader, "Muon_fired_HLT_Mu8_Barrel_L1HP9_IP6"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_IP3 = {fReader, "Muon_fired_HLT_Mu8_IP3"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_IP5 = {fReader, "Muon_fired_HLT_Mu8_IP5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8_IP6 = {fReader, "Muon_fired_HLT_Mu8_IP6"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu8p5_IP3p5 = {fReader, "Muon_fired_HLT_Mu8p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu9_Barrel_L1HP10_IP6 = {fReader, "Muon_fired_HLT_Mu9_Barrel_L1HP10_IP6"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu9_IP4 = {fReader, "Muon_fired_HLT_Mu9_IP4"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu9_IP5 = {fReader, "Muon_fired_HLT_Mu9_IP5"};
   TTreeReaderArray<Int_t> Muon_fired_HLT_Mu9_IP6 = {fReader, "Muon_fired_HLT_Mu9_IP6"};
   TTreeReaderArray<Int_t> Muon_isGlobalMuon = {fReader, "Muon_isGlobalMuon"};
   TTreeReaderArray<Int_t> Muon_isPF = {fReader, "Muon_isPF"};
   TTreeReaderArray<Int_t> Muon_isTrackerMuon = {fReader, "Muon_isTrackerMuon"};
   TTreeReaderArray<Int_t> Muon_isTriggering = {fReader, "Muon_isTriggering"};
   TTreeReaderArray<Int_t> Muon_looseID = {fReader, "Muon_looseID"};
   TTreeReaderArray<Int_t> Muon_numberOfPixelLayers = {fReader, "Muon_numberOfPixelLayers"};
   TTreeReaderArray<Int_t> Muon_numberOfStations = {fReader, "Muon_numberOfStations"};
   TTreeReaderArray<Int_t> Muon_numberOfTrackerLayers = {fReader, "Muon_numberOfTrackerLayers"};
   TTreeReaderArray<Int_t> Muon_numberOfValidMuonHits = {fReader, "Muon_numberOfValidMuonHits"};
   TTreeReaderArray<Int_t> Muon_numberOfValidPixelHits = {fReader, "Muon_numberOfValidPixelHits"};
   TTreeReaderArray<Int_t> Muon_pdgId = {fReader, "Muon_pdgId"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu10_Barrel_L1HP11_IP6 = {fReader, "Muon_prescale_HLT_Mu10_Barrel_L1HP11_IP6"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu10p5_IP3p5 = {fReader, "Muon_prescale_HLT_Mu10p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu12_IP6 = {fReader, "Muon_prescale_HLT_Mu12_IP6"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu6_Barrel_L1HP7_IP6 = {fReader, "Muon_prescale_HLT_Mu6_Barrel_L1HP7_IP6"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu7_Barrel_L1HP8_IP6 = {fReader, "Muon_prescale_HLT_Mu7_Barrel_L1HP8_IP6"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu7_IP4 = {fReader, "Muon_prescale_HLT_Mu7_IP4"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_Barrel_L1HP9_IP6 = {fReader, "Muon_prescale_HLT_Mu8_Barrel_L1HP9_IP6"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_IP3 = {fReader, "Muon_prescale_HLT_Mu8_IP3"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_IP5 = {fReader, "Muon_prescale_HLT_Mu8_IP5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8_IP6 = {fReader, "Muon_prescale_HLT_Mu8_IP6"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu8p5_IP3p5 = {fReader, "Muon_prescale_HLT_Mu8p5_IP3p5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu9_Barrel_L1HP10_IP6 = {fReader, "Muon_prescale_HLT_Mu9_Barrel_L1HP10_IP6"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu9_IP4 = {fReader, "Muon_prescale_HLT_Mu9_IP4"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu9_IP5 = {fReader, "Muon_prescale_HLT_Mu9_IP5"};
   TTreeReaderArray<Int_t> Muon_prescale_HLT_Mu9_IP6 = {fReader, "Muon_prescale_HLT_Mu9_IP6"};
   TTreeReaderArray<Int_t> Muon_trackerHighPurityFlag = {fReader, "Muon_trackerHighPurityFlag"};
   TTreeReaderArray<Float_t> Muon_caloCompatibility = {fReader, "Muon_caloCompatibility"};
   TTreeReaderArray<Float_t> Muon_covLamLam = {fReader, "Muon_covLamLam"};
   TTreeReaderArray<Float_t> Muon_covLamPhi = {fReader, "Muon_covLamPhi"};
   TTreeReaderArray<Float_t> Muon_covPhiPhi = {fReader, "Muon_covPhiPhi"};
   TTreeReaderArray<Float_t> Muon_covQopLam = {fReader, "Muon_covQopLam"};
   TTreeReaderArray<Float_t> Muon_covQopPhi = {fReader, "Muon_covQopPhi"};
   TTreeReaderArray<Float_t> Muon_covQopQop = {fReader, "Muon_covQopQop"};
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
   TTreeReaderArray<Float_t> Muon_ip3d = {fReader, "Muon_ip3d"};
   TTreeReaderArray<Float_t> Muon_kinkFinderChi2 = {fReader, "Muon_kinkFinderChi2"};
   TTreeReaderArray<Float_t> Muon_localPositionChi2 = {fReader, "Muon_localPositionChi2"};
   TTreeReaderArray<Float_t> Muon_mass = {fReader, "Muon_mass"};
   TTreeReaderArray<Float_t> Muon_matched_dpt = {fReader, "Muon_matched_dpt"};
   TTreeReaderArray<Float_t> Muon_matched_dr = {fReader, "Muon_matched_dr"};
   TTreeReaderArray<Float_t> Muon_pfiso03_all = {fReader, "Muon_pfiso03_all"};
   TTreeReaderArray<Float_t> Muon_pfiso03_sumChargedHadronPt = {fReader, "Muon_pfiso03_sumChargedHadronPt"};
   TTreeReaderArray<Float_t> Muon_pfiso03_sumChargedParticlePt = {fReader, "Muon_pfiso03_sumChargedParticlePt"};
   TTreeReaderArray<Float_t> Muon_pfiso03_sumNeutralHadronEt = {fReader, "Muon_pfiso03_sumNeutralHadronEt"};
   TTreeReaderArray<Float_t> Muon_pfiso03_sumNeutralHadronEtHighThreshold = {fReader, "Muon_pfiso03_sumNeutralHadronEtHighThreshold"};
   TTreeReaderArray<Float_t> Muon_pfiso03_sumPUPt = {fReader, "Muon_pfiso03_sumPUPt"};
   TTreeReaderArray<Float_t> Muon_pfiso03_sumPhotonEt = {fReader, "Muon_pfiso03_sumPhotonEt"};
   TTreeReaderArray<Float_t> Muon_pfiso03_sumPhotonEtHighThreshold = {fReader, "Muon_pfiso03_sumPhotonEtHighThreshold"};
   TTreeReaderArray<Float_t> Muon_pfiso04_all = {fReader, "Muon_pfiso04_all"};
   TTreeReaderArray<Float_t> Muon_pfiso04_sumChargedHadronPt = {fReader, "Muon_pfiso04_sumChargedHadronPt"};
   TTreeReaderArray<Float_t> Muon_pfiso04_sumChargedParticlePt = {fReader, "Muon_pfiso04_sumChargedParticlePt"};
   TTreeReaderArray<Float_t> Muon_pfiso04_sumNeutralHadronEt = {fReader, "Muon_pfiso04_sumNeutralHadronEt"};
   TTreeReaderArray<Float_t> Muon_pfiso04_sumNeutralHadronEtHighThreshold = {fReader, "Muon_pfiso04_sumNeutralHadronEtHighThreshold"};
   TTreeReaderArray<Float_t> Muon_pfiso04_sumPUPt = {fReader, "Muon_pfiso04_sumPUPt"};
   TTreeReaderArray<Float_t> Muon_pfiso04_sumPhotonEt = {fReader, "Muon_pfiso04_sumPhotonEt"};
   TTreeReaderArray<Float_t> Muon_pfiso04_sumPhotonEtHighThreshold = {fReader, "Muon_pfiso04_sumPhotonEtHighThreshold"};
   TTreeReaderArray<Float_t> Muon_phi = {fReader, "Muon_phi"};
   TTreeReaderArray<Float_t> Muon_pt = {fReader, "Muon_pt"};
   TTreeReaderArray<Float_t> Muon_ptErr = {fReader, "Muon_ptErr"};
   TTreeReaderArray<Float_t> Muon_segmentCompatibility = {fReader, "Muon_segmentCompatibility"};
   TTreeReaderArray<Float_t> Muon_sip3d = {fReader, "Muon_sip3d"};
   TTreeReaderArray<Float_t> Muon_validHitFraction = {fReader, "Muon_validHitFraction"};
   TTreeReaderArray<Float_t> Muon_vx = {fReader, "Muon_vx"};
   TTreeReaderArray<Float_t> Muon_vy = {fReader, "Muon_vy"};
   TTreeReaderArray<Float_t> Muon_vz = {fReader, "Muon_vz"};
   TTreeReaderValue<Int_t> nTriggerMuon = {fReader, "nTriggerMuon"};
   TTreeReaderArray<Int_t> TriggerMuon_charge = {fReader, "TriggerMuon_charge"};
   TTreeReaderArray<Int_t> TriggerMuon_pdgId = {fReader, "TriggerMuon_pdgId"};
   TTreeReaderArray<Float_t> TriggerMuon_eta = {fReader, "TriggerMuon_eta"};
   TTreeReaderArray<Float_t> TriggerMuon_mass = {fReader, "TriggerMuon_mass"};
   TTreeReaderArray<Float_t> TriggerMuon_phi = {fReader, "TriggerMuon_phi"};
   TTreeReaderArray<Float_t> TriggerMuon_pt = {fReader, "TriggerMuon_pt"};
   TTreeReaderArray<Float_t> TriggerMuon_vx = {fReader, "TriggerMuon_vx"};
   TTreeReaderArray<Float_t> TriggerMuon_vy = {fReader, "TriggerMuon_vy"};
   TTreeReaderArray<Float_t> TriggerMuon_vz = {fReader, "TriggerMuon_vz"};
   TTreeReaderValue<Float_t> Rho_fixedGridRhoAll = {fReader, "Rho_fixedGridRhoAll"};
   TTreeReaderValue<Float_t> Rho_fixedGridRhoFastjetAll = {fReader, "Rho_fixedGridRhoFastjetAll"};
   TTreeReaderValue<Float_t> Rho_fixedGridRhoFastjetCentral = {fReader, "Rho_fixedGridRhoFastjetCentral"};
   TTreeReaderValue<Float_t> Rho_fixedGridRhoFastjetCentralCalo = {fReader, "Rho_fixedGridRhoFastjetCentralCalo"};
   TTreeReaderValue<Float_t> Rho_fixedGridRhoFastjetCentralChargedPileUp = {fReader, "Rho_fixedGridRhoFastjetCentralChargedPileUp"};
   TTreeReaderValue<Float_t> Rho_fixedGridRhoFastjetCentralNeutral = {fReader, "Rho_fixedGridRhoFastjetCentralNeutral"};
   TTreeReaderValue<Int_t> nOtherPV = {fReader, "nOtherPV"};
   TTreeReaderArray<Float_t> OtherPV_z = {fReader, "OtherPV_z"};
   TTreeReaderArray<Float_t> OtherPV_score = {fReader, "OtherPV_score"};
   TTreeReaderValue<UChar_t> PV_npvs = {fReader, "PV_npvs"};
   TTreeReaderValue<UChar_t> PV_npvsGood = {fReader, "PV_npvsGood"};
   TTreeReaderValue<Float_t> PV_ndof = {fReader, "PV_ndof"};
   TTreeReaderValue<Float_t> PV_x = {fReader, "PV_x"};
   TTreeReaderValue<Float_t> PV_y = {fReader, "PV_y"};
   TTreeReaderValue<Float_t> PV_z = {fReader, "PV_z"};
   TTreeReaderValue<Float_t> PV_chi2 = {fReader, "PV_chi2"};
   TTreeReaderValue<Float_t> PV_score = {fReader, "PV_score"};
   TTreeReaderValue<Float_t> PV_sumpt2 = {fReader, "PV_sumpt2"};
   TTreeReaderValue<Float_t> PV_sumpx = {fReader, "PV_sumpx"};
   TTreeReaderValue<Float_t> PV_sumpy = {fReader, "PV_sumpy"};
   TTreeReaderValue<Int_t> nSV = {fReader, "nSV"};
   TTreeReaderArray<Short_t> SV_charge = {fReader, "SV_charge"};
   TTreeReaderArray<Float_t> SV_dlen = {fReader, "SV_dlen"};
   TTreeReaderArray<Float_t> SV_dlenSig = {fReader, "SV_dlenSig"};
   TTreeReaderArray<Float_t> SV_dxy = {fReader, "SV_dxy"};
   TTreeReaderArray<Float_t> SV_dxySig = {fReader, "SV_dxySig"};
   TTreeReaderArray<Float_t> SV_pAngle = {fReader, "SV_pAngle"};
   TTreeReaderArray<UChar_t> SV_ntracks = {fReader, "SV_ntracks"};
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
   TTreeReaderValue<Int_t> nGenPart = {fReader, "nBsToPhiPhiTo4K"};
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
   //TTreeReaderValue<UInt_t> Pileup_nPU = {fReader, "PV_npvs"};
   //TTreeReaderValue<Float_t> Pileup_nTrueInt = {fReader, "PV_ndof"};

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
   Float_t the_bs_lx_posthepv = -99.;
   Float_t the_bs_ly_posthepv = -99.;
   Float_t the_bs_ct_2d_cm = -99.;
   Float_t the_bs_ct_2d_cm_posbsz = -99.;
   Float_t the_bs_ct_2d_cm_posbspv = -99.;
   Float_t the_bs_ct_2d_cm_posthepv = -99.;
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

   Float_t the_cos_theta_k1 = -99.;
   Float_t the_cos_theta_k3 = -99.;
   Float_t the_phi_star = -99.;

   Float_t the_pv_chi2 = -99.;
   Float_t the_pv_covXX = -99.;
   Float_t the_pv_covXY = -99.;
   Float_t the_pv_covXZ = -99.;
   Float_t the_pv_covYY = -99.;
   Float_t the_pv_covYZ = -99.;
   Float_t the_pv_covZZ = -99.;
   Float_t the_pv_ndof = -99.;
   Float_t the_pv_x = -99.;
   Float_t the_pv_y = -99.;
   Float_t the_pv_z = -99.;

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
   Float_t the_phi1_ct_2d_cm = -99.;
   Float_t the_phi1_ct_2d_cm_posbsz = -99.;
   Float_t the_phi1_ct_2d_cm_posbspv = -99.;
   Float_t the_phi1_k1_drtrg = -99.;
   Float_t the_phi1_k2_drtrg = -99.;
   Float_t the_phi1_k1_dztrg = -99.;
   Float_t the_phi1_k2_dztrg = -99.;
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
   Float_t the_phi2_ct_2d_cm = -99.;
   Float_t the_phi2_ct_2d_cm_posbsz = -99.;
   Float_t the_phi2_ct_2d_cm_posbspv = -99.;
   Float_t the_phi2_k1_drtrg = -99.;
   Float_t the_phi2_k2_drtrg = -99.;
   Float_t the_phi2_k1_dztrg = -99.;
   Float_t the_phi2_k2_dztrg = -99.;
   Int_t the_phi2_is_matched = -99;

   Float_t the_phi1_pt_times_phi2_pt = -99.;

   Float_t the_k1_eta = -99.;
   Float_t the_k1_phi = -99.;
   Float_t the_k1_pt = -99.;
   Float_t the_k1_dcasig = -99.;
   Float_t the_k1_dcasig_bestpv = -99.;
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
   Float_t the_k1_dxy_bestpv = -99.;
   Float_t the_k1_dxysig_bestpv = -99.;
   Float_t the_k1_dz_bestpv = -99.;
   Float_t the_k1_dzsig_bestpv = -99.;
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
   Float_t the_k2_dcasig_bestpv = -99.;
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
   Float_t the_k2_dxy_bestpv = -99.;
   Float_t the_k2_dxysig_bestpv = -99.;
   Float_t the_k2_dz_bestpv = -99.;
   Float_t the_k2_dzsig_bestpv = -99.;
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
   Float_t the_k3_dcasig_bestpv = -99.;
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
   Float_t the_k3_dxy_bestpv = -99.;
   Float_t the_k3_dxysig_bestpv = -99.;
   Float_t the_k3_dz_bestpv = -99.;
   Float_t the_k3_dzsig_bestpv = -99.;
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
   Float_t the_k4_dcasig_bestpv = -99.;
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
   Float_t the_k4_dxy_bestpv = -99.;
   Float_t the_k4_dxysig_bestpv = -99.;
   Float_t the_k4_dz_bestpv = -99.;
   Float_t the_k4_dzsig_bestpv = -99.;
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

   Int_t ismatched = -99;
   Int_t the_k1_genidx = -99;
   Int_t the_k2_genidx = -99;
   Int_t the_k3_genidx = -99;
   Int_t the_k4_genidx = -99;
   Int_t the_k1mother_genidx = -99;
   Int_t the_k2mother_genidx = -99;
   Int_t the_k3mother_genidx = -99;
   Int_t the_k4mother_genidx = -99;
   Int_t the_k1grandmother_genidx = -99;
   Int_t the_k2grandmother_genidx = -99;
   Int_t the_k3grandmother_genidx = -99;
   Int_t the_k4grandmother_genidx = -99;
   Int_t the_k1grandgrandmother_genidx = -99;
   Int_t the_k2grandgrandmother_genidx = -99;
   Int_t the_k3grandgrandmother_genidx = -99;
   Int_t the_k4grandgrandmother_genidx = -99;

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
