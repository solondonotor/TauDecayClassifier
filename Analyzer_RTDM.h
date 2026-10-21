#ifndef FCCANA_ADDITIONAL_ANALYZERS_H
#define FCCANA_ADDITIONAL_ANALYZERS_H

// =====================================================================
// Tau decay classification from the MC record (EDM4hep)
//
// Three levels of classification:
//   1. classifyTauDecay    : mode 1..46 from the DIRECT daughters of the tau
//   2. classifyFullDecay   : resonance sub-mode (100..503), looks one level
//                            deeper (rho, K*, omega, a1)
//   3. classifyPionicCode  : numeric code ALWAYS 5 digits J Kc Kn Nc Nn
//                            (Kc=Kn=0 when no kaons, making their absence explicit)
//
// Conventions
//   * Photons are NOT counted in nTotal (FSR / PHOTOS must not change the mode).
//   * K0, anti-K0, K0S and K0L are all added to nK0.
//   * Everything is classified in the tau- frame. A tau+ is charge-conjugated
//     first (see conjugate / orient), so there is a single copy of the logic.
//   * Intermediate states (rho, omega, a1, K*, eta) are only assigned a J when
//     their decay is checked and allowed (charge, Bose symmetry). Otherwise J=9.
//
// Error codes returned by classifyPionicCode
//   -1  no hadrons found / not classified
//   -2  charge of the final state differs from the tau charge
//   -3  recursion depth exceeded (truncated tree)
//   -4  some counter is > 9 (code would be ambiguous)
//
// J (intermediate state)
//   0 pure (no resonance)          5 a1 -> 3 pions (non resonant)
//   1 via rho-                     6 via omega
//   2 via rho0 (verified rho0->pi+pi-; should be rare/absent)
//   3 via a1 -> rho0 pi            7 via K* (charged or neutral)
//   4 via a1 -> rho- pi0           8 via eta
//                                  9 other / unverified resonance
//
// Pionic code: ALWAYS J*10000 + Kc*1000 + Kn*100 + Nc*10 + Nn  (5 digits)
//   When no kaons: Kc=Kn=0, so 011 (old 3-digit) becomes 00011, 630 becomes 60030.
//  
//
// Special single-body codes: 1 e, 2 mu, 3 pi, 4 K, 145 e+gamma, 246 mu+gamma
//
// Mode table (classifyTauDecay), all written for tau-:
//   1 e nu nu | 2 mu nu nu | 3 pi nu | 4 K nu | 5 pi pi0 | 6 K pi0
//   7 pi 2pi0 | 8 pi pi pi+ | 10 K*- | 11 pi pi pi+ pi0 | 12 pi 3pi0
//   13 K K0 | 14 pi K0 K0bar | 15 K pi0 K0 | 16 pi K0 pi0 | 17 pi K0
//   18 omega pi | 19 eta pi | 20 eta K | 21 eta K*- | 22 eta K pi0
//   23 3pi- 2pi+ | 24 3pi- 2pi+ pi0 | 25 pi pi pi+ 2pi0 | 26 K pi+ pi-
//   27 K pi+ pi- pi0 | 28 K pi- K+ | 29 K K- K+ | 30 K K+ pi- pi0
//   31 e e e+ nu nu | 32 K1- | 33 K*0bar pi | 34 K*0 K | 35 phi pi
//   36 phi K | 37 pi K0 2pi0 | 38 pi K0 K0 pi0 | 40 f1 pi
//   41 eta pi pi0 | 42 eta K0 pi | 43 pi pi pi+ 3pi0 | 44 K0 rho-
//   45 e nu nu gamma | 46 mu nu nu gamma
//   (mode 9 removed: tau- -> K*0 nu violates charge; mode 39 removed:
//    identical to 30; 19/18 entries with extra pi0 removed)
//
// Full-decay codes (classifyFullDecay):
//   100 rho-  -> pi- pi0
//   110 K*-   -> K- pi0        111 K*- -> K0bar pi-
//   200 pi- rho0 -> pi+ pi-    (201 rho0 -> pi0 pi0 is forbidden: removed)
//   210 K*0 (K*0bar) -> K pi (charged)   211 -> K0 pi0   (modes 33, 34)
//   300 omega -> pi+ pi- pi0   301 omega -> pi0 gamma   302 omega -> pi+ pi-
//   400 a1 -> rho0 pi, rho0 -> pi+ pi-
//   401 a1 -> rho- pi0, rho- -> pi- pi0
//   402 a1 -> pi- pi0 pi0      403 a1 -> pi- pi- pi+  (non resonant)
//   500 pi- pi0 | 501 pi- pi- pi+ | 502 pi- 2pi0 | 503 pi- 3pi0   (pure)
// =====================================================================
 

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>

#include "ROOT/RVec.hxx"
#include "edm4hep/MCParticleData.h"

namespace FCCAnalyses::MCDecayMode
{

using PartVec = ROOT::VecOps::RVec<edm4hep::MCParticleData>;
using IntVec  = ROOT::VecOps::RVec<int>;

constexpr int KAON_CODE_OFFSET = 100000;

// =====================================================================
// Counters of direct daughters
// =====================================================================
struct DaughterCounts {
    int nEle = 0, nEleBar = 0, nMu = 0, nMuBar = 0;
    int nNuE = 0, nNuEBar = 0, nNuMu = 0, nNuMuBar = 0, nNuTau = 0, nNuTauBar = 0;
    int nGamma = 0;
    int nPiPlus = 0, nPiMinus = 0, nPi0 = 0;
    int nKPlus = 0, nKMinus = 0;
    int nK0 = 0;                 // K0, anti-K0, K0S, K0L
    int nK0L = 0, nK0S = 0;      // informational only (also included in nK0)
    int nRhoPlus = 0, nRhoMinus = 0, nRho0 = 0;
    int nKstarPlus = 0, nKstarMinus = 0, nKstar0 = 0, nKstar0Bar = 0;
    int nOmega = 0, nEta = 0, nPhi = 0;
    int nA1Plus = 0, nA1Minus = 0, nF1 = 0;
    int nK1minus = 0, nK1plus = 0;          // K1(1270)  10323
    int nK1_minus = 0, nK1_plus = 0;        // K1(1400)  20323
    int nK1_starminus = 0, nK1_starplus = 0;// K*(1410)  100323
    int nTotal = 0;                          // all daughters EXCEPT photons
};

// Charge conjugation of the counters (tau+ -> tau- frame)
inline DaughterCounts conjugate(DaughterCounts c)
{
    std::swap(c.nEle, c.nEleBar);
    std::swap(c.nMu, c.nMuBar);
    std::swap(c.nNuE, c.nNuEBar);
    std::swap(c.nNuMu, c.nNuMuBar);
    std::swap(c.nNuTau, c.nNuTauBar);
    std::swap(c.nPiPlus, c.nPiMinus);
    std::swap(c.nKPlus, c.nKMinus);
    std::swap(c.nRhoPlus, c.nRhoMinus);
    std::swap(c.nKstarPlus, c.nKstarMinus);
    std::swap(c.nKstar0, c.nKstar0Bar);
    std::swap(c.nA1Plus, c.nA1Minus);
    std::swap(c.nK1plus, c.nK1minus);
    std::swap(c.nK1_plus, c.nK1_minus);
    std::swap(c.nK1_starplus, c.nK1_starminus);
    return c;
}

inline DaughterCounts orient(const DaughterCounts& c, bool isTauMinus)
{
    return isTauMinus ? c : conjugate(c);
}

// Charge carried by the pions and charged kaons
inline int chargeOf(const DaughterCounts& c)
{
    return c.nPiPlus - c.nPiMinus + c.nKPlus - c.nKMinus;
}

inline DaughterCounts countDaughters(size_t i, const PartVec& particles, const IntVec& daughters)
{
    DaughterCounts c;
    const size_t nP = particles.size();
    const size_t nD = daughters.size();

    for (unsigned j = particles[i].daughters_begin; j < particles[i].daughters_end; ++j) {
        if (j >= nD) continue;
        int di = daughters[j];
        if (di < 0 || static_cast<size_t>(di) >= nP) continue;
        const int pdg = particles[di].PDG;

        if (pdg == 22) { c.nGamma++; continue; }   // photons never enter nTotal
        c.nTotal++;

        switch (pdg) {
            case   11: c.nEle++;       break;
            case  -11: c.nEleBar++;    break;
            case   13: c.nMu++;        break;
            case  -13: c.nMuBar++;     break;
            case   12: c.nNuE++;       break;
            case  -12: c.nNuEBar++;    break;
            case   14: c.nNuMu++;      break;
            case  -14: c.nNuMuBar++;   break;
            case   16: c.nNuTau++;     break;
            case  -16: c.nNuTauBar++;  break;
            case  211: c.nPiPlus++;    break;
            case -211: c.nPiMinus++;   break;
            case  111: c.nPi0++;       break;
            case  321: c.nKPlus++;     break;
            case -321: c.nKMinus++;    break;
            case  311: case -311: c.nK0++;                 break;
            case  130: c.nK0L++; c.nK0++;                  break;
            case  310: c.nK0S++; c.nK0++;                  break;
            case  213: c.nRhoPlus++;   break;
            case -213: c.nRhoMinus++;  break;
            case  113: c.nRho0++;      break;
            case  323: c.nKstarPlus++; break;
            case -323: c.nKstarMinus++;break;
            case  313: c.nKstar0++;    break;
            case -313: c.nKstar0Bar++; break;
            case  223: c.nOmega++;     break;
            case  221: c.nEta++;       break;
            case  333: c.nPhi++;       break;
            case  20213: c.nA1Plus++;  break;
            case -20213: c.nA1Minus++; break;
            case  20223: c.nF1++;      break;
            case  10323: c.nK1plus++;  break;
            case -10323: c.nK1minus++; break;
            case  20323: c.nK1_plus++; break;
            case -20323: c.nK1_minus++;break;
            case  100323: c.nK1_starplus++;  break;
            case -100323: c.nK1_starminus++; break;
            default: break;
        }
    }
    return c;
}

// =====================================================================
// Record navigation helpers
// =====================================================================

// Follow X -> X(copy) chains (same PDG, exactly one such daughter)
inline int lastCopy(int idx, const PartVec& p, const IntVec& d)
{
    for (int guard = 0; guard < 20; ++guard) {
        int next = -1, n = 0;
        for (unsigned j = p[idx].daughters_begin; j < p[idx].daughters_end; ++j) {
            if (j >= d.size()) continue;
            int di = d[j];
            if (di < 0 || static_cast<size_t>(di) >= p.size()) continue;
            if (p[di].PDG == p[idx].PDG) { next = di; ++n; }
        }
        if (n != 1) break;
        idx = next;
    }
    return idx;
}

// First daughter with that PDG, then follow its copy chain
inline int findDaughterWithPDG(size_t parent_idx, int target_pdg,
                               const PartVec& p, const IntVec& d)
{
    const size_t nP = p.size();
    const size_t nD = d.size();
    for (unsigned j = p[parent_idx].daughters_begin; j < p[parent_idx].daughters_end; ++j) {
        if (j >= nD) continue;
        int di = d[j];
        if (di < 0 || static_cast<size_t>(di) >= nP) continue;
        if (p[di].PDG == target_pdg) return lastCopy(di, p, d);
    }
    return -1;
}

// A tau is "last" if none of its daughters is a tau (skips tau -> tau copies)
inline bool isLastTau(size_t i, const PartVec& p, const IntVec& d)
{
    for (unsigned j = p[i].daughters_begin; j < p[i].daughters_end; ++j) {
        if (j >= d.size()) continue;
        int di = d[j];
        if (di < 0 || static_cast<size_t>(di) >= p.size()) continue;
        if (std::abs(p[di].PDG) == 15) return false;
    }
    return true;
}

// =====================================================================
// Checks on how a resonance decayed (all in the tau- frame)
// =====================================================================

// rho: charge -1 -> pi- pi0 ; charge 0 -> pi+ pi- only (pi0 pi0 is forbidden)
inline bool rhoDecayOK(const DaughterCounts& rc, int charge)
{
    if (rc.nTotal != 2) return false;
    if (charge == 0) return rc.nPiPlus == 1 && rc.nPiMinus == 1;
    return rc.nPi0 == 1 && (rc.nPiPlus + rc.nPiMinus) == 1 && chargeOf(rc) == charge;
}

// omega: returns 300 / 301 / 302, or -1 if not an allowed decay
inline int omegaDecayMode(const DaughterCounts& oc)
{
    if (oc.nTotal == 3 && oc.nPiPlus == 1 && oc.nPiMinus == 1 && oc.nPi0 == 1) return 300;
    if (oc.nTotal == 1 && oc.nPi0 == 1 && oc.nGamma >= 1)                     return 301;
    if (oc.nTotal == 2 && oc.nPiPlus == 1 && oc.nPiMinus == 1)                return 302;
    return -1;
}

// K*: one kaon + one pion, charge conserved, strangeness not violated.
// strangeSign: -1 for K*-, K*0bar ; +1 for K*0
inline bool kstarDecayOK(const DaughterCounts& kc, int charge, int strangeSign)
{
    const int nK  = kc.nKPlus + kc.nKMinus + kc.nK0;
    const int nPi = kc.nPiPlus + kc.nPiMinus + kc.nPi0;
    if (kc.nTotal != 2 || nK != 1 || nPi != 1) return false;
    if (chargeOf(kc) != charge) return false;
    if ((kc.nKPlus - kc.nKMinus) * strangeSign < 0) return false;
    return true;
}

// eta -> gamma gamma, 3pi0, pi+pi-pi0, pi+pi-gamma
inline bool etaDecayOK(const DaughterCounts& ec)
{
    if (ec.nTotal == 0 && ec.nGamma >= 2) return true;
    if (ec.nTotal == 3 && ec.nPi0 == 3) return true;
    if (ec.nTotal == 3 && ec.nPiPlus == 1 && ec.nPiMinus == 1 && ec.nPi0 == 1) return true;
    if (ec.nTotal == 2 && ec.nPiPlus == 1 && ec.nPiMinus == 1 && ec.nGamma >= 1) return true;
    return false;
}

// =====================================================================
// classifyTauDecay : mode 1..46 from the direct daughters
// =====================================================================
inline int classifyTauDecay(const DaughterCounts& cin, int tauPDG)
{
    const DaughterCounts c = orient(cin, tauPDG == 15);   // tau- frame
    if (c.nNuTau != 1) return -1;

    switch (c.nTotal) {
    case 2:
        if (c.nPiMinus == 1)    return 3;
        if (c.nKMinus == 1)     return 4;
        if (c.nKstarMinus == 1) return 10;
        if (c.nK1minus == 1 || c.nK1_minus == 1 || c.nK1_starminus == 1) return 32;
        break;
    case 3:
        if (c.nEle == 1 && c.nNuEBar  == 1) return c.nGamma == 0 ? 1 : 45;
        if (c.nMu  == 1 && c.nNuMuBar == 1) return c.nGamma == 0 ? 2 : 46;
        if (c.nPiMinus == 1 && c.nPi0 == 1)       return 5;
        if (c.nKMinus  == 1 && c.nPi0 == 1)       return 6;
        if (c.nKMinus  == 1 && c.nK0  == 1)       return 13;
        if (c.nPiMinus == 1 && c.nK0  == 1)       return 17;
        if (c.nK0 == 1 && c.nRhoMinus == 1)       return 44;
        if (c.nOmega == 1 && c.nPiMinus == 1)     return 18;
        if (c.nEta == 1 && c.nPiMinus == 1)       return 19;
        if (c.nEta == 1 && c.nKMinus == 1)        return 20;
        if (c.nEta == 1 && c.nKstarMinus == 1)    return 21;
        if (c.nKstar0Bar == 1 && c.nPiMinus == 1) return 33;
        if (c.nKstar0 == 1 && c.nKMinus == 1)     return 34;
        if (c.nPhi == 1 && c.nPiMinus == 1)       return 35;
        if (c.nPhi == 1 && c.nKMinus == 1)        return 36;
        if (c.nF1 == 1  && c.nPiMinus == 1)       return 40;
        break;
    case 4:
        if (c.nPiMinus == 2 && c.nPiPlus == 1)                   return 8;
        if (c.nPiMinus == 1 && c.nPi0 == 2)                      return 7;
        if (c.nPiMinus == 1 && c.nK0 == 2)                       return 14;
        if (c.nKMinus == 1 && c.nPi0 == 1 && c.nK0 == 1)         return 15;
        if (c.nPiMinus == 1 && c.nK0 == 1 && c.nPi0 == 1)        return 16;
        if (c.nKMinus == 1 && c.nPiPlus == 1 && c.nPiMinus == 1) return 26;
        if (c.nKMinus == 1 && c.nPiMinus == 1 && c.nKPlus == 1)  return 28;
        if (c.nKMinus == 2 && c.nKPlus == 1)                     return 29;
        if (c.nEta == 1 && c.nPiMinus == 1 && c.nPi0 == 1)       return 41;
        if (c.nEta == 1 && c.nKMinus == 1 && c.nPi0 == 1)        return 22;
        if (c.nEta == 1 && c.nK0 == 1 && c.nPiMinus == 1)        return 42;
        break;
    case 5:
        if (c.nPiMinus == 2 && c.nPiPlus == 1 && c.nPi0 == 1)               return 11;
        if (c.nPiMinus == 1 && c.nPi0 == 3)                                 return 12;
        if (c.nKMinus == 1 && c.nPiPlus == 1 && c.nPiMinus == 1 && c.nPi0 == 1) return 27;
        if (c.nKMinus == 1 && c.nKPlus == 1 && c.nPiMinus == 1 && c.nPi0 == 1)  return 30;
        if (c.nEle == 2 && c.nEleBar == 1 && c.nNuEBar == 1)                return 31;
        if (c.nPiMinus == 1 && c.nK0 == 1 && c.nPi0 == 2)                   return 37;
        if (c.nPiMinus == 1 && c.nK0 == 2 && c.nPi0 == 1)                   return 38;
        break;
    case 6:
        if (c.nPiMinus == 3 && c.nPiPlus == 2)                              return 23;
        if (c.nPiMinus == 2 && c.nPiPlus == 1 && c.nPi0 == 2)               return 25;
        break;
    case 7:
        if (c.nPiMinus == 3 && c.nPiPlus == 2 && c.nPi0 == 1)               return 24;
        if (c.nPiMinus == 2 && c.nPiPlus == 1 && c.nPi0 == 3)               return 43;
        break;
    default:
        break;
    }
    return -1;
}

// =====================================================================
// Resonance analysis (shared by classifyFullDecay and getIntermediateJ)
// =====================================================================
struct ResonanceInfo {
    int J        = 0;    // intermediate state index (see header)
    int fullMode = -1;   // 100..403, or -1
};

inline ResonanceInfo analyseResonance(size_t tau_idx, const PartVec& p, const IntVec& d,
                                      const DaughterCounts& raw, bool isTauMinus)
{
    ResonanceInfo r;
    const DaughterCounts q = orient(raw, isTauMinus);   // tau- frame
    const int sg = isTauMinus ? 1 : -1;                 // PDG sign for the real record

    auto decayOf = [&](int idx) {
        return orient(countDaughters(static_cast<size_t>(idx), p, d), isTauMinus);
    };

    // ---- J=1: rho- -------------------------------------------------
    if (q.nRhoMinus >= 1) {
        int ri = findDaughterWithPDG(tau_idx, -213 * sg, p, d);
        if (ri >= 0 && rhoDecayOK(decayOf(ri), -1)) {
            r.J = 1;
            if (q.nTotal == 2 && q.nRhoMinus == 1) r.fullMode = 100;
        } else r.J = 9;
        return r;
    }

    // ---- J=2: rho0 (only if it really decays to pi+ pi-) -------------
    if (q.nRho0 >= 1) {
        int ri = findDaughterWithPDG(tau_idx, 113, p, d);
        if (ri >= 0 && rhoDecayOK(decayOf(ri), 0)) {
            r.J = 2;
            if (q.nTotal == 3 && q.nRho0 == 1 && q.nPiMinus == 1) r.fullMode = 200;
        } else r.J = 9;
        return r;
    }

    // ---- J=3,4,5: a1- ----------------------------------------------
    if (q.nA1Minus >= 1) {
        int ai = findDaughterWithPDG(tau_idx, -20213 * sg, p, d);
        if (ai < 0) { r.J = 9; return r; }
        const DaughterCounts ac = decayOf(ai);
        int fm = -1;
        if (ac.nTotal == 2 && ac.nRho0 == 1 && ac.nPiMinus == 1) {
            int ri = findDaughterWithPDG(ai, 113, p, d);
            if (ri >= 0 && rhoDecayOK(decayOf(ri), 0)) { r.J = 3; fm = 400; }
            else r.J = 9;
        } else if (ac.nTotal == 2 && ac.nRhoMinus == 1 && ac.nPi0 == 1) {
            int ri = findDaughterWithPDG(ai, -213 * sg, p, d);
            if (ri >= 0 && rhoDecayOK(decayOf(ri), -1)) { r.J = 4; fm = 401; }
            else r.J = 9;
        } else if (ac.nTotal == 3 && (ac.nPiPlus + ac.nPiMinus + ac.nPi0) == 3 &&
                   chargeOf(ac) == -1) {
            r.J = 5;
            if (ac.nPi0 == 2) fm = 402;
            else if (ac.nPi0 == 0) fm = 403;
        } else r.J = 9;
        if (r.J != 9 && q.nTotal == 2 && q.nA1Minus == 1) r.fullMode = fm;
        return r;
    }

    // ---- J=6: omega -------------------------------------------------
    if (q.nOmega >= 1) {
        int oi = findDaughterWithPDG(tau_idx, 223, p, d);
        int om = (oi >= 0) ? omegaDecayMode(decayOf(oi)) : -1;
        if (om < 0) { r.J = 9; return r; }
        r.J = 6;
        if (q.nTotal == 3 && q.nOmega == 1 && q.nPiMinus == 1) r.fullMode = om;
        return r;
    }

    // ---- J=7: K* ----------------------------------------------------
    if (q.nKstarMinus >= 1 || q.nKstar0 >= 1 || q.nKstar0Bar >= 1) {
        bool ok = false;
        int fm = -1;
        if (q.nKstarMinus >= 1) {
            int ki = findDaughterWithPDG(tau_idx, -323 * sg, p, d);
            if (ki >= 0) {
                const DaughterCounts kc = decayOf(ki);
                if (kstarDecayOK(kc, -1, -1)) {
                    ok = true;
                    if (q.nTotal == 2 && q.nKstarMinus == 1) fm = (kc.nKMinus == 1) ? 110 : 111;
                }
            }
        } else {
            int ki, strange;
            if (q.nKstar0Bar >= 1) { ki = findDaughterWithPDG(tau_idx, -313 * sg, p, d); strange = -1; }
            else                   { ki = findDaughterWithPDG(tau_idx,  313 * sg, p, d); strange = +1; }
            if (ki >= 0) {
                const DaughterCounts kc = decayOf(ki);
                if (kstarDecayOK(kc, 0, strange)) {
                    ok = true;
                    // modes 33 (K*0bar pi-) and 34 (K*0 K-)
                    if (q.nTotal == 3 &&
                        ((q.nKstar0Bar == 1 && q.nPiMinus == 1) ||
                         (q.nKstar0 == 1 && q.nKMinus == 1)))
                        fm = (kc.nK0 == 1) ? 211 : 210;
                }
            }
        }
        r.J = ok ? 7 : 9;
        r.fullMode = ok ? fm : -1;
        return r;
    }

    // ---- J=8: eta ---------------------------------------------------
    if (q.nEta >= 1) {
        int ei = findDaughterWithPDG(tau_idx, 221, p, d);
        r.J = (ei >= 0 && etaDecayOK(decayOf(ei))) ? 8 : 9;
        return r;
    }

    // ---- J=9: other resonances or wrong-sign charged ones -----------
    if (q.nPhi || q.nF1 ||
        q.nK1minus || q.nK1plus || q.nK1_minus || q.nK1_plus ||
        q.nK1_starminus || q.nK1_starplus ||
        q.nRhoPlus || q.nKstarPlus || q.nA1Plus) {
        r.J = 9;
        return r;
    }

    r.J = 0;   // pure
    return r;
}

// =====================================================================
// classifyFullDecay : resonance sub-mode (100..403) or pure (500..503)
// =====================================================================
inline int classifyFullDecay(size_t tau_idx, const PartVec& p, const IntVec& d,
                             const DaughterCounts& tauC, bool isTauMinus)
{
    ResonanceInfo r = analyseResonance(tau_idx, p, d, tauC, isTauMinus);
    if (r.fullMode >= 0) return r.fullMode;

    const DaughterCounts q = orient(tauC, isTauMinus);
    if (q.nNuTau != 1) return -1;

    const bool noRes = (q.nRhoPlus == 0 && q.nRhoMinus == 0 && q.nRho0 == 0 &&
                        q.nKstarPlus == 0 && q.nKstarMinus == 0 &&
                        q.nKstar0 == 0 && q.nKstar0Bar == 0 &&
                        q.nOmega == 0 && q.nA1Plus == 0 && q.nA1Minus == 0 &&
                        q.nEta == 0 && q.nPhi == 0 && q.nF1 == 0);
    if (!noRes) return -1;

    if (q.nPiMinus == 1 && q.nPi0 == 1 && q.nTotal == 3)                    return 500;
    if (q.nPiMinus == 2 && q.nPiPlus == 1 && q.nTotal == 4)                 return 501;
    if (q.nPiMinus == 1 && q.nPi0 == 2 && q.nTotal == 4)                    return 502;
    if (q.nPiMinus == 1 && q.nPi0 == 3 && q.nTotal == 5)                    return 503;
    return -1;
}

// =====================================================================
// Final state (stable pions and kaons), recursive
// =====================================================================
struct FinalStateCounts {
    int  Nc = 0;      // charged pions
    int  Nn = 0;      // pi0
    int  Kc = 0;      // charged kaons
    int  Kn = 0;      // neutral kaons (K0, K0bar, KS, KL)
    int  nEta = 0;    // etas met while descending
    int  nGam = 0;    // photons met
    int  charge = 0;  // net charge of pi+- and K+-
    bool truncated = false;
};

// Count-and-stop at pions, kaons, photons; ignore leptons and neutrinos;
// keep descending through everything else (rho, omega, a1, K*, eta, ...).
inline void collectFinalState(size_t idx, const PartVec& p, const IntVec& d,
                              FinalStateCounts& fsc, int depth = 0)
{
    if (depth > 12) { fsc.truncated = true; return; }
    const size_t nP = p.size();
    const size_t nD = d.size();

    const int pdg  = p[idx].PDG;
    const int apdg = std::abs(pdg);

    switch (apdg) {
        case 211: fsc.Nc++; fsc.charge += (pdg > 0 ? +1 : -1); return;
        case 111: fsc.Nn++; return;
        case 321: fsc.Kc++; fsc.charge += (pdg > 0 ? +1 : -1); return;
        case 311: case 130: case 310: fsc.Kn++; return;
        case  22: fsc.nGam++; return;
        case  11: case 13: case 12: case 14: case 16: return;
        case 221: fsc.nEta++; break;     // eta: count it AND descend
        default: break;
    }

    if (p[idx].daughters_begin >= p[idx].daughters_end) return;

    for (unsigned j = p[idx].daughters_begin; j < p[idx].daughters_end; ++j) {
        if (j >= nD) continue;
        int di = d[j];
        if (di < 0 || static_cast<size_t>(di) >= nP) continue;
        const int dpdg = std::abs(p[di].PDG);
        if (dpdg == 12 || dpdg == 14 || dpdg == 16) continue;
        collectFinalState(static_cast<size_t>(di), p, d, fsc, depth + 1);
    }
}

inline int getIntermediateJ(size_t tau_idx, const PartVec& p, const IntVec& d,
                            const DaughterCounts& tauC, bool isTauMinus)
{
    return analyseResonance(tau_idx, p, d, tauC, isTauMinus).J;
}

// =====================================================================
// classifyPionicCode
// =====================================================================
inline int classifyPionicCode(size_t tau_idx, const PartVec& p, const IntVec& d,
                              const DaughterCounts& tauC, int tauMode, bool isTauMinus)
{
    // one-body and leptonic modes
    switch (tauMode) {
        case 1:  return 1;
        case 2:  return 2;
        case 3:  return 3;
        case 4:  return 4;
        case 45: return 145;
        case 46: return 246;
        default: break;
    }

    FinalStateCounts fsc;
    collectFinalState(tau_idx, p, d, fsc);

    const int Nc = fsc.Nc, Nn = fsc.Nn, Kc = fsc.Kc, Kn = fsc.Kn;

    if (Nc == 0 && Nn == 0 && Kc == 0 && Kn == 0) return -1;
    if (fsc.truncated)                            return -3;
    if (fsc.charge != (isTauMinus ? -1 : +1))     return -2;
    if (Nc > 9 || Nn > 9 || Kc > 9 || Kn > 9)     return -4;

    const int J = getIntermediateJ(tau_idx, p, d, tauC, isTauMinus);

  //  if (Kc > 0 || Kn > 0)
  //     return KAON_CODE_OFFSET + J * 10000 + Kc * 1000 + Kn * 100 + Nc * 10 + Nn;

    return J * 10000 + Kc * 1000 + Kn * 100 + Nc * 10 + Nn;
}

// =====================================================================
// One loop for all public functions
// =====================================================================
enum class ModeKind { Basic, Full, Pionic };

inline IntVec tauModesImpl(const PartVec& p, const IntVec& d, int pdgSel, ModeKind kind)
{
    IntVec modes;
    for (size_t i = 0; i < p.size(); ++i) {
        const int pdg = p[i].PDG;
        if (std::abs(pdg) != 15) continue;
        if (pdgSel != 0 && pdg != pdgSel) continue;
        if (!isLastTau(i, p, d)) continue;

        const bool isTauMinus = (pdg == 15);
        const DaughterCounts c = countDaughters(i, p, d);

        switch (kind) {
            case ModeKind::Basic:
                modes.push_back(classifyTauDecay(c, pdg));
                break;
            case ModeKind::Full:
                modes.push_back(classifyFullDecay(i, p, d, c, isTauMinus));
                break;
            case ModeKind::Pionic: {
                const int mode = classifyTauDecay(c, pdg);
                modes.push_back(classifyPionicCode(i, p, d, c, mode, isTauMinus));
                break;
            }
        }
    }
    return modes;
}

// =====================================================================
// Public functions (RDataFrame .Define)
//   .Define("tauMode",       "FCCAnalyses::MCDecayMode::tauDecayModes(MCParticles, MCParticles#1.index)")
//   .Define("tauFullMode",   "FCCAnalyses::MCDecayMode::tauFullModes(MCParticles, MCParticles#1.index)")
//   .Define("tauPionicMode", "FCCAnalyses::MCDecayMode::tauPionicModes(MCParticles, MCParticles#1.index)")
// =====================================================================
inline IntVec tauDecayModes(const PartVec& p, const IntVec& d)        { return tauModesImpl(p, d,   0, ModeKind::Basic);  }
inline IntVec tauMinusModes(const PartVec& p, const IntVec& d)        { return tauModesImpl(p, d,  15, ModeKind::Basic);  }
inline IntVec tauPlusModes(const PartVec& p, const IntVec& d)         { return tauModesImpl(p, d, -15, ModeKind::Basic);  }
inline IntVec tauFullModes(const PartVec& p, const IntVec& d)         { return tauModesImpl(p, d,   0, ModeKind::Full);   }
inline IntVec tauPionicModes(const PartVec& p, const IntVec& d)       { return tauModesImpl(p, d,   0, ModeKind::Pionic); }
inline IntVec tauMinusPionicModes(const PartVec& p, const IntVec& d)  { return tauModesImpl(p, d,  15, ModeKind::Pionic); }
inline IntVec tauPlusPionicModes(const PartVec& p, const IntVec& d)   { return tauModesImpl(p, d, -15, ModeKind::Pionic); }

// =====================================================================
// Debug: "daughter dump"
//
// Prints the decay tree (3 levels) of the last taus whose pionic code is
// targetCode, up to maxPrint prints in total. Use with a single thread.
//
//   .Define("dbg", "FCCAnalyses::MCDecayMode::dumpTauDecaysWithCode(
//                     MCParticles, MCParticles#1.index, 613, 20)")
//
// RDataFrame is lazy: the column must be used by some action
// (Snapshot, Histo, Max, ...) for the printing to happen.
// =====================================================================
inline std::string decayString(int idx, const PartVec& p, const IntVec& d, int depth = 0)
{
    std::string s = std::to_string(p[idx].PDG);
    const int a = std::abs(p[idx].PDG);
    if (depth >= 3 || a == 211 || a == 111 || a == 321 || a == 22 || a == 11 ||
        a == 13 || a == 12 || a == 14 || a == 16 || a == 130 || a == 310 || a == 311)
        return s;

    std::string kids;
    for (unsigned j = p[idx].daughters_begin; j < p[idx].daughters_end; ++j) {
        if (j >= d.size()) continue;
        int di = d[j];
        if (di < 0 || static_cast<size_t>(di) >= p.size()) continue;
        kids += " " + decayString(di, p, d, depth + 1);
    }
    if (!kids.empty()) s += " [" + kids + " ]";
    return s;
}

inline int dumpTauDecaysWithCode(const PartVec& p, const IntVec& d, int targetCode, int maxPrint)
{
    static int nPrinted = 0;
    for (size_t i = 0; i < p.size(); ++i) {
        const int pdg = p[i].PDG;
        if (std::abs(pdg) != 15 || !isLastTau(i, p, d)) continue;
        if (nPrinted >= maxPrint) return nPrinted;

        const DaughterCounts c = countDaughters(i, p, d);
        const int mode = classifyTauDecay(c, pdg);
        const int code = classifyPionicCode(i, p, d, c, mode, pdg == 15);
        if (code != targetCode) continue;

        std::cout << "[dump] code=" << code << "  " << decayString(static_cast<int>(i), p, d) << std::endl;
        ++nPrinted;
    }
    return nPrinted;
}

} // namespace FCCAnalyses::MCDecayMode
#endif
