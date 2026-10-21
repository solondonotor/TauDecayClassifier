## 1. Contents and requirements

The repository contains:

- `Analyzer_RTDM.h` — C++ implementation of the tau decay classification.
- `run_analysis.py` — runs the analysis over the selected sample.
- `show_table.py` — reads the output ROOT file and displays the decay-mode table.

The analysis requires an FCCAnalyses environment with ROOT, Python 3 and a C++17-compatible environment.

---

## 2. Running the analysis

The default sample is configured in `run_analysis.py`:

processList = {
    "p8_ee_Ztautau_ecm91": {},}

Run:
python run_analysis.py

The output is saved in:
result_RTDM/

For quick testing, the current script processes only the first 500,000 entries:
.Filter("0 < rdfentry_ && rdfentry_<= 500000")


To process the full sample, remove or comment out this filter.

To test another sample, replace its name in processList.

---

## 3. Tau decay classification

The classification is performed using:
FCCAnalyses::MCDecayMode::tauMinusPionicModes(Particle, Daughters)

and saved in the output branch:
tauMinusPionicMode

Hadronic decays use the following five-position code:
J Kc Kn Nc Nn

where:
- J = intermediate resonance/topology
- Kc = charged kaons
- Kn = neutral kaons
- Nc = charged pions
- Nn = neutral pions
The resonance/topology code J is:
J	Topology
0	direct / no identified resonance
1	via rho-
2	via rho0
3	via a1 -> rho0 pi
4	via a1 -> rho- pi0
5	via a1 (no identified resonance)
6	via omega
7	via K*
8	via eta
9	other


Special codes are:
Code	Meaning
1	e- nu_tau
2	mu- nu_tau
3	pi- nu_tau
4	K- nu_tau
145	e- nu_tau gamma
246	mu- nu_tau gamma
-1	no classification
-2	incorrect charge
-3	recursion truncated
-4	counter > 9


## 4. Inspecting the output

After running the analysis:

python show_table.py result_RTDM/SAMPLE_NAME.root

For example, for the default sample:
python show_table.py result_RTDM/p8_ee_Ztautau_ecm91.root

The script displays:
- decay-mode code
- number of occurrences
- relative fraction
- human-readable description
Since the codes are stored as integers, leading zeros are lost in the ROOT file. For example:
00110 -> 110

## 5. Typical workflow
FCC-ee MC sample
       |
       v
run_analysis.py
       |
       v
Analyzer_RTDM.h
       |
       v
tauMinusPionicMode
       |
       v
result_RTDM/*.root
       |
       v
show_table.py
       |
       v
Decay-mode classification table

## 6. Quick test
For the default sample, simply run:
python run_analysis.py

and then:
python show_table.py

To test another sample:
1. Change the sample name in processList.
2. Run python run_analysis.py.
3. Point show_table.py to the corresponding output ROOT file.
4. Run python show_table.py.
