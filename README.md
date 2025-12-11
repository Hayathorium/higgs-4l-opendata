# Data Requirements

To run the analysis, **you must download the following ROOT files** from the official CERN Open Data portal and place them in the same directory as `code.py`:

**Download page:** [https://opendata.cern/record/15005](https://opendata.cern/record/15005)

## Required Files

### **Real Data Files**

These contain four-lepton candidate events collected during different ATLAS data-taking periods:

| File Name          | Description                                     |
| ------------------ | ----------------------------------------------- |
| `data_A.4lep.root` | ATLAS Run A data containing 4-lepton candidates |
| `data_B.4lep.root` | ATLAS Run B data containing 4-lepton candidates |
| `data_C.4lep.root` | ATLAS Run C data containing 4-lepton candidates |
| `data_D.4lep.root` | ATLAS Run D data containing 4-lepton candidates |

---

### **Signal Monte Carlo Samples**

Simulated Higgs boson production with ( H \to ZZ^* \to 4\ell ):

| File Name                            | Description                                |
| ------------------------------------ | ------------------------------------------ |
| `mc_341947.ZH125_ZZ4lep.4lep.root`   | ZH associated production ((m_H = 125) GeV) |
| `mc_341964.WH125_ZZ4lep.4lep.root`   | WH associated production                   |
| `mc_344235.VBFH125_ZZ4lep.4lep.root` | Vector Boson Fusion (VBF) Higgs            |
| `mc_345060.ggH125_ZZ4lep.4lep.root`  | Gluon–gluon fusion (ggF) Higgs             |

---

### **Background Monte Carlo Samples**

| File Name                  | Description                             |
| -------------------------- | --------------------------------------- |
| `mc_363490.llll.4lep.root` | Irreducible (ZZ^* \to 4\ell) background |

---

