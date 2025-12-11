import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import pandas as pd
from scipy.special import gamma

needed = [
    "lep_n", "lep_pt", "lep_eta", "lep_phi", "lep_E", "lep_isTightID", "lep_type", "lep_charge", "lep_ptcone30", "lep_z0", "lep_trackd0pvunbiased",
    "mcWeight", "XSection", "SumWeights",
    "scaleFactor_PILEUP",
    "scaleFactor_ELE",
    "scaleFactor_MUON",
    "scaleFactor_PHOTON",
    "scaleFactor_TAU",
    "scaleFactor_BTAG",
    "scaleFactor_LepTRIGGER",
    "scaleFactor_PhotonTRIGGER"
]

files1 = ["mc_363490.llll.4lep.root"]
files2 = [
    "mc_341947.ZH125_ZZ4lep.4lep.root",
    "mc_341964.WH125_ZZ4lep.4lep.root",
    "mc_344235.VBFH125_ZZ4lep.4lep.root",
    "mc_345060.ggH125_ZZ4lep.4lep.root"
]
files3 = ["data_A.4lep.root","data_B.4lep.root","data_C.4lep.root","data_D.4lep.root"]
cutflow_table = pd.DataFrame()

def count_events(label, data):
    if "Data" in label:
        return len(data)
    else:
        return float(ak.sum(data["weight"]))

def invariant_mass(E, px, py, pz):
        return np.sqrt(E**2 - px**2 - py**2 - pz**2) / 1000.0

# filtering motivated by https://www.sciencedirect.com/science/article/pii/S037026931200857X
def filter(files, label):
    data = ak.concatenate([uproot.open(f)["mini"].arrays(needed, library="ak") for f in files])
    # add invariant-mass and weight
    E = ak.sum(data["lep_E"], axis=1)
    px = ak.sum(data["lep_pt"]*np.cos(data["lep_phi"]), axis=1)
    py = ak.sum(data["lep_pt"]*np.sin(data["lep_phi"]), axis=1)
    pz = ak.sum(data["lep_pt"]*np.sinh(data["lep_eta"]), axis=1)
    m4l = invariant_mass(E, px, py, pz)
    data = ak.with_field(data, m4l, "m4l")
    data = ak.with_field(data, data["mcWeight"]*data["scaleFactor_PILEUP"]*data["scaleFactor_ELE"]*data["scaleFactor_MUON"]*data["scaleFactor_PHOTON"]*data["scaleFactor_TAU"]*data["scaleFactor_BTAG"]*data["scaleFactor_LepTRIGGER"]*data["scaleFactor_PhotonTRIGGER"]*data["XSection"]*10*1000/data["SumWeights"], "weight")
    cutflow_table.loc["All events (preselection)",label] = count_events(label, data)
    # Events with exactly 4 reconstructed leptons
    data = data[data["lep_n"]==4]
    cutflow_table.loc["Exactly 4 leptons",label] = count_events(label, data)
    # Pseudorapidity acceptance for leptons
    data = data[ak.all(ak.where(data["lep_type"] == 11, np.abs(data["lep_eta"]) < 2.47, True), axis=1)
        & ak.all(ak.where(data["lep_type"] == 13, np.abs(data["lep_eta"]) < 2.7, True), axis=1)]
    cutflow_table.loc["Lepton |eta| acceptance",label] = count_events(label, data)
    # Transverse momentum thresholds for the 4 leptons
    pT_sort_indices = ak.argsort(data["lep_pt"], ascending=False, axis=1)
    pt_sorted = data["lep_pt"][pT_sort_indices]
    lep1_pt = pt_sorted[:, 0]
    lep2_pt = pt_sorted[:, 1]
    lep3_pt = pt_sorted[:, 2]
    lep4_pt = pt_sorted[:, 3]
    data = data[(lep1_pt > 20000.0)&(lep2_pt > 15000.0)&(lep3_pt > 10000.0)&(lep4_pt > 7000.0)]
    cutflow_table.loc["Lepton pT thresholds cut",label] = count_events(label, data)
    # Minimum angular separation dR between leptons
    eta = data["lep_eta"]
    phi = data["lep_phi"]
    lep_type = data["lep_type"]
    def deltaR(eta1, phi1, eta2, phi2):
        dphi = np.abs(phi1 - phi2)
        dphi = ak.where(dphi > np.pi, 2*np.pi - dphi, dphi)
        return np.sqrt((eta1 - eta2)**2 + dphi**2)
    deltaR_mask = ak.ones_like(eta[:,0], dtype=bool)
    for i in range(4):
        for j in range(i+1, 4):
            dR = deltaR(eta[:,i], phi[:,i], eta[:,j], phi[:,j])
            cut = ak.where(lep_type[:,i] == lep_type[:,j], dR > 0.1, dR > 0.2)
            deltaR_mask = deltaR_mask & cut
    data = data[deltaR_mask]
    cutflow_table.loc["Lepton isolation (dR) cut",label] = count_events(label, data)
    # Leptons that pass the Tight identification criteria: https://cds.cern.ch/record/2842463/files/2211.16178.pdf
    data = data[ak.all(data["lep_isTightID"], axis=1)]
    cutflow_table.loc["Tight lepton ID cut",label] = count_events(label, data)
    # |z0| < 10 mm to ensure leptons originate near the primary vertex
    data = data[ak.all(np.abs(data["lep_z0"]) < 10.0, axis=1)]
    cutflow_table.loc["Lepton longitudinal impact (z0) cut",label] = count_events(label, data)
    # |d0| < 1 mm to reject cosmic-ray muons
    data = data[ak.all(ak.where(data["lep_type"] == 13, np.abs(data["lep_trackd0pvunbiased"]) < 1.0, True), axis=1)]
    cutflow_table.loc["Muon transverse impact parameter (d0) cut",label] = count_events(label, data)
    # Select isolated leptons pTcone/pT < 0.15 to reject jets
    data = data[ak.all(data["lep_ptcone30"]/data["lep_pt"] < 0.15, axis=1)]
    cutflow_table.loc["Lepton isolation (pTcone) cut",label] = count_events(label, data)
    # Select opposite-sign, same-flavor lepton pairs consistent with Z boson mass (m12 and m34 cuts)
    lep_idx = ak.local_index(data["lep_E"], axis=1)
    lep_p4 = ak.zip({
        "E": data["lep_E"], 
        "px": data["lep_pt"] * np.cos(data["lep_phi"]), 
        "py": data["lep_pt"] * np.sin(data["lep_phi"]),
        "pz": data["lep_pt"] * np.sinh(data["lep_eta"]), 
        "charge": data["lep_charge"], "type": data["lep_type"], "idx": lep_idx
    })
    pairs = ak.combinations(lep_p4, 2)
    pairs = pairs[(pairs["0"].charge * pairs["1"].charge == -1) & (pairs["0"].type == pairs["1"].type)]
    quads = ak.combinations(pairs, 2)
    quads = quads[(quads["0"]["0"].idx != quads["1"]["0"].idx) & \
                (quads["0"]["0"].idx != quads["1"]["1"].idx) & \
                (quads["0"]["1"].idx != quads["1"]["0"].idx) & \
                (quads["0"]["1"].idx != quads["1"]["1"].idx)]
    m_Z = 91.1880 #https://pdg.lbl.gov/
    m12 = invariant_mass(
        quads["0"]["0"].E + quads["0"]["1"].E,
        quads["0"]["0"].px + quads["0"]["1"].px,
        quads["0"]["0"].py + quads["0"]["1"].py,
        quads["0"]["0"].pz + quads["0"]["1"].pz
    )
    m34 = invariant_mass(
        quads["1"]["0"].E + quads["1"]["1"].E,
        quads["1"]["0"].px + quads["1"]["1"].px,
        quads["1"]["0"].py + quads["1"]["1"].py,
        quads["1"]["0"].pz + quads["1"]["1"].pz
    )
    swap = np.abs(m34 - m_Z) < np.abs(m12 - m_Z)
    m12, m34 = ak.where(swap, m34, m12), ak.where(swap, m12, m34)
    E_sum = quads["0"]["0"].E + quads["0"]["1"].E + quads["1"]["0"].E + quads["1"]["1"].E
    px_sum = quads["0"]["0"].px + quads["0"]["1"].px + quads["1"]["0"].px + quads["1"]["1"].px
    py_sum = quads["0"]["0"].py + quads["0"]["1"].py + quads["1"]["0"].py + quads["1"]["1"].py
    pz_sum = quads["0"]["0"].pz + quads["0"]["1"].pz + quads["1"]["0"].pz + quads["1"]["1"].pz
    m4l = invariant_mass(E_sum, px_sum, py_sum, pz_sum)
    m12_cut = (m12 > 50.0) & (m12 < 106.0)
    m_min = ak.where(
        m4l < 120, 17.5,
        ak.where(m4l > 190, 50.0, 17.5 + (50-17.5)/(190-120)*(m4l-120))
    )
    m34_cut = (m34 > m_min) & (m34 < 115.0)
    quad_cut = m12_cut & m34_cut
    data = data[ak.num(quads[quad_cut], axis=1) > 0]
    cutflow_table.loc["Z candidate selection",label] = count_events(label, data)
    #output
    return data


background = filter(files1, "Background (MC)")
signal = filter(files2, "Signal (MC)")
observed = filter(files3, "Data (observed)")
print(cutflow_table)
mc_count = len(ak.to_numpy(background["m4l"]))
obs_count = len(ak.to_numpy(observed["m4l"][(80<=observed["m4l"])&(observed["m4l"]<=250)]))

# plot normalized MC simulation data
bins = np.linspace(80, 250, 34)
B, edges = np.histogram(
    ak.to_numpy(background["m4l"]),
    bins=bins,
    weights=ak.to_numpy(background["weight"])
)
S, _ = np.histogram(
    ak.to_numpy(signal["m4l"]),
    bins=bins,
    weights=ak.to_numpy(signal["weight"])
)
plt.figure(figsize=(11, 6))
plt.bar((edges[:-1]+edges[1:])/2, B, width=np.diff(edges), color="#98c5e6", label=f"$4l$ Background (MC)")
plt.bar((edges[:-1]+edges[1:])/2, S, width=np.diff(edges), color="#f0bcbc", bottom=B, label="Higgs Signal (ZH, WH, VBF, ggH) (MC)")
# plot experimental data
counts, edges = np.histogram(ak.to_numpy(observed["m4l"]), bins=bins)
plt.errorbar((edges[:-1]+edges[1:])/2, counts, yerr=np.sqrt(counts), fmt='o', color='black', label="Observed Data")
# calculate expected significance and observed significance
def negative_log_likelihood(mu, k_i, s_i, b_i):
    lambda_i = mu * s_i + b_i
    return -np.sum(k_i * np.log(lambda_i) - lambda_i - np.log(gamma(k_i+1)))
mu_hat = minimize(negative_log_likelihood, x0=[1.0], args=(counts, S,B), bounds=[(0, None)],).x[0]
print("mu_hat_MLE =", mu_hat)
print("Expected Significance", np.sqrt(2*(negative_log_likelihood(0,S+B,S,B)-negative_log_likelihood(1,S+B,S,B))))
print("Observed Significance", np.sqrt(2*(negative_log_likelihood(0,counts,S,B)-negative_log_likelihood(mu_hat,counts,S,B))))
plt.axvline(125, color='red', linestyle='--')
plt.xlabel(f"$m_{{4l}}$ [GeV]")
plt.ylabel("Events / 5 GeV")
plt.title(f"Observed vs Expected $m_{{4l}}$ Distribution")
plt.legend()
plt.show()

