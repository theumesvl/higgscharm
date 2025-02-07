import numpy as np
import awkward as ak
from coffea.nanoevents.methods import candidate


def delta_r_higher(first, second, threshold=0.4):
    # select objects from 'first' which are at least 'threshold' away from all objects in 'second'.
    mval = first.metric_table(second)
    return ak.all(mval > threshold, axis=-1)

def delta_r_lower(first, second, threshold=0.4):
    # select objects from 'first' which are at least 'threshold' within from all objects in 'second'.
    mval = first.metric_table(second)
    return ak.all(mval <= threshold, axis=-1)

# calculate transverse mass mT for a specific lepton or lepton pair and the total MET of the event
def transverse_mass(lep, met):
    return np.sqrt(2 * lep.pt * met.pt * (ak.ones_like(met.pt) - np.cos(lep.delta_phi(met))))
    
def select_leading_lepton (leptons):
    leading_lepton = leptons[ak.argmax(leptons.pt, axis=1, keepdims=True)]
    

def select_dileptons(objects, key):
    leptons = ak.zip(
        {
            "pt": objects[key].pt,
            "eta": objects[key].eta,
            "phi": objects[key].phi,
            "mass": objects[key].mass,
            "charge": objects[key].charge,
        },
        with_name="PtEtaPhiMCandidate",
        behavior=candidate.behavior,
    )
    # make sure they are sorted by transverse momentum
    leptons = leptons[ak.argsort(leptons.pt, axis=1)]
    # create pair combinations with all leptons
    dileptons = ak.combinations(leptons, 2, fields=["l1", "l2"])
    # add dimuon 4-momentum field
    dileptons["ll"] = dileptons.l1 + dileptons.l2
    dileptons["pt"] = dileptons.ll.pt
    return dileptons

# #NOT to be used
# def select_WW_candidates(events, muons, electrons, MET):
#     """
#     Selects WW candidates where one W decays to a muon and the other to an electron.
    
#     Parameters:
#     - events: NanoAOD event collection
#     - muons: Selected muon collection
#     - electrons: Selected electron collection
#     - MET: Missing transverse energy (neutrino proxy)
    
#     Returns:
#     - WW candidates as awkward arrays
#     """
#     leading_muon = select_leading_lepton(muons)
#     leading_electron = select_leading_lepton(electrons)
#     # W1 Candidate: Muon + MET
#     mt_w1 = transverse_mass(leading_muon, MET)

#     # W2 Candidate: Electron + MET
#     mt_w2 = transverse_mass(leading_electron, MET)
    
#     # Create Lorentz vectors for leptons
#     muon_p4 = ak.zip(
#         {
#             "pt": leading_muon.pt, 
#             "eta": leading_muon.eta, 
#             "phi": leading_muon.phi, 
#             "mass": leading_muon.mass,
#             "charge": leading_muon.charge
#         },
#         with_name="PtEtaPhiMCandidate",
#         behavior=candidate.behavior,
#     )
    
#     electron_p4 = ak.zip(
#         {
#             "pt": leading_electron.pt, 
#             "eta": leading_electron.eta, 
#             "phi": leading_electron.phi, 
#             "mass": leading_electron.mass,
#             "charge": leading_electron.charge
#         },
#         with_name="PtEtaPhiMCandidate",
#         behavior=candidate.behavior,
#     )

# def select_WWtoMuEle_WW_candidates(ll_pairs):
#     ww_pairs = ak.combinations(ll_pairs, 2, fields=["w1", "w2"])
#     ww_pairs = ak.zip(
#         {
#             "w1": ak.zip(
#                 {
#                     "l1": ww_pairs.w1.l1,
#                     "l2": ww_pairs.w1.l2,
#                     "p4": ww_pairs.w1.l1 + ww_pairs.w1.l2,
#                 }
#             ),
#             "w2": ak.zip(
#                 {
#                     "l1": ww_pairs.w2.l1,
#                     "l2": ww_pairs.w2.l2,
#                     "p4": ww_pairs.w2.l1 + ww_pairs.w2.l2,
#                 }
#             ),
#         }
#     )
#     # ghost removal: ∆R(η, φ) > 0.02 between each of the four leptons
#     ghost_removal_mask = (
#         (ww_pairs.w1.l1.delta_r(ww_pairs.w1.l2) > 0.02)
#         & (ww_pairs.w1.l1.delta_r(ww_pairs.w2.l1) > 0.02)
#         & (ww_pairs.w1.l1.delta_r(ww_pairs.w2.l2) > 0.02)
#         & (ww_pairs.w1.l2.delta_r(ww_pairs.w2.l1) > 0.02)
#         & (ww_pairs.w1.l2.delta_r(ww_pairs.w2.l2) > 0.02)
#         & (ww_pairs.w2.l1.delta_r(ww_pairs.w2.l2) > 0.02)
#     )
#     # Lepton pT: two of the four selected leptons should pass pT,i > 20 GeV and pT,j > 10
#     lepton_pt_mask = (
#         (
#             ak.any(ww_pairs.w1.l1.pt > 20, axis=-1)
#             & ak.any(ww_pairs.w1.l2.pt > 10, axis=-1)
#         )
#         | (
#             ak.any(ww_pairs.w1.l1.pt > 20, axis=-1)
#             & ak.any(ww_pairs.w2.l1.pt > 10, axis=-1)
#         )
#         | (
#             ak.any(ww_pairs.w1.l1.pt > 20, axis=-1)
#             & ak.any(ww_pairs.w2.l2.pt > 10, axis=-1)
#         )
#         | (
#             ak.any(ww_pairs.w1.l2.pt > 20, axis=-1)
#             & ak.any(ww_pairs.w1.l1.pt > 10, axis=-1)
#         )
#         | (
#             ak.any(ww_pairs.w1.l2.pt > 20, axis=-1)
#             & ak.any(ww_pairs.w2.l1.pt > 10, axis=-1)
#         )
#         | (
#             ak.any(ww_pairs.w1.l2.pt > 20, axis=-1)
#             & ak.any(ww_pairs.w2.l2.pt > 10, axis=-1)
#         )
#         | (
#             ak.any(ww_pairs.w2.l1.pt > 20, axis=-1)
#             & ak.any(ww_pairs.w1.l1.pt > 10, axis=-1)
#         )
#         | (
#             ak.any(ww_pairs.w2.l1.pt > 20, axis=-1)
#             & ak.any(ww_pairs.w1.l2.pt > 10, axis=-1)
#         )
#         | (
#             ak.any(ww_pairs.w2.l1.pt > 20, axis=-1)
#             & ak.any(ww_pairs.w2.l2.pt > 10, axis=-1)
#         )
#         | (
#             ak.any(ww_pairs.w2.l2.pt > 20, axis=-1)
#             & ak.any(ww_pairs.w1.l1.pt > 10, axis=-1)
#         )
#         | (
#             ak.any(ww_pairs.w2.l2.pt > 20, axis=-1)
#             & ak.any(ww_pairs.w1.l2.pt > 10, axis=-1)
#         )
#         | (
#             ak.any(ww_pairs.w2.l2.pt > 20, axis=-1)
#             & ak.any(ww_pairs.w2.l1.pt > 10, axis=-1)
#         )
#     )
#     # QCD suppression: all four opposite-sign pairs that can be built with the four leptons (regardless of lepton flavor) must satisfy m > 4 GeV
#     qcd_condition_mask = (
#         (ww_pairs.w1.p4.mass > 4)
#         & (ww_pairs.w2.p4.mass > 4)
#         & (
#             (ww_pairs.w1.l1.charge + ww_pairs.w2.l1.charge != 0)
#             | (
#                 (ww_pairs.w1.l1.charge + ww_pairs.w2.l1.charge == 0)
#                 & ((ww_pairs.w1.l1 + ww_pairs.w2.l1).mass > 4)
#             )
#         )
#         & (
#             (ww_pairs.w1.l1.charge + ww_pairs.w2.l2.charge != 0)
#             | (
#                 (ww_pairs.w1.l1.charge + ww_pairs.w2.l2.charge == 0)
#                 & ((ww_pairs.w1.l1 + ww_pairs.w2.l2).mass > 4)
#             )
#         )
#         & (
#             (ww_pairs.w1.l2.charge + ww_pairs.w2.l1.charge != 0)
#             | (
#                 (ww_pairs.w1.l2.charge + ww_pairs.w2.l1.charge == 0)
#                 & ((ww_pairs.w1.l2 + ww_pairs.w2.l1).mass > 4)
#             )
#         )
#         & (
#             (ww_pairs.w1.l2.charge + ww_pairs.w2.l2.charge != 0)
#             | (
#                 (ww_pairs.w1.l2.charge + ww_pairs.w2.l2.charge == 0)
#                 & ((ww_pairs.w1.l2 + ww_pairs.w2.l2).mass > 4)
#             )
#         )
#     )
#     # w1 mass > 40 GeV
#     mass_mask = ww_pairs.w1.p4.mass > 40

#     mask = ghost_removal_mask & lepton_pt_mask & qcd_condition_mask & mass_mask
#     return ww_pairs[ak.fill_none(mask, False)]
