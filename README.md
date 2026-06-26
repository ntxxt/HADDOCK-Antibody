# HADDOCK3 Antibody–Antigen Docking Tutorial Dataset

This repository contains all input data, restraint files, and configuration scripts required to reproduce a HADDOCK3 antibody–antigen docking tutorial. The dataset supports multiple docking scenarios using AlphaFold3 predictions, CDR-based restraints, and experimentally derived interface information.

---

## 📁 Contents

### 🧬 Structures
- `antigen.pdb` — AlphaFold3-prepared antigen structure (Chain A)
- `antibody_chainB.pdb` — processed antibody structure (Chain B)
- `complex_prepared.pdb` — antibody–antigen complex for restraint generation
- `4G6M_ref.pdb` — experimental reference structure (CAPRI evaluation)

### 📄 Sequences & annotations
- `antibody.fasta` — antibody heavy/light chain sequences
- `active.txt` — CDR-derived antibody active residues
- `passive.txt` — antigen surface residues
- `antigen.txt` — antigen residue list for restraint generation
- `access_ag.txt` — antigen solvent accessibility output

### 🔗 Restraints (HADDOCK3)
- `ti.tbl` — true interface restraints (reference structure)
- `ti_AF3.tbl` — AF3-derived interface restraints
- `ti_exp.tbl` — experimental + CDR-based restraints
- `unambig.tbl` — restraints keeping antibody chains together

### ⚙️ Configuration files
- `hd3.cfg` — Scenario 1 (AF3-guided docking)
- `hd3_exp.cfg` — Scenario 3 (experimental restraints)
- `hd3_af3.job` — job submission script

---

## 🧪 Docking scenarios

1. **AF3-guided docking**  
   Uses AlphaFold3-predicted interfaces to define restraints.

2. **CDR + surface-guided docking**  
   Uses antibody CDR residues and antigen surface accessibility.

3. **Experimental restraint docking**  
   Uses NMR/literature epitope information combined with CDR data.

---

## 🛠️ Tools used

- AlphaFold3 (antigen prediction)
- ImmuneBuilder / ABodyBuilder2 (antibody modeling)
- HADDOCK3 (docking and restraint generation)
- pdb-tools (structure processing)

---

## 🚀 Running HADDOCK3

Each scenario can be executed using:

```bash
haddock3 hd3.cfg
```

---

## 📊 Output analysis

Docking results are evaluated using CAPRI-style metrics (DockQ, iRMSD, fnat). Final models are ranked based on the HADDOCK score after refinement and clustering.
