import ROOT
import sys
from collections import Counter

if len(sys.argv) < 2:
    print("Usage: python show_table.py <input.root>")
    sys.exit(1)

input_file = sys.argv[1]

f = ROOT.TFile(input_file)
tree = f.Get("events")

pionic_counts = Counter()
for event in tree:
    for m in event.tauMinusPionicMode:
        pionic_counts[m] += 1

def code_str(code):
    if code in (1, 2, 3, 4, 145, 246, -1, -2, -3, -4):
        return str(code)
    return f"{code:05d}"   # siempre 5 dígitos

def describe(code):
    special = {
        1:   "e- nu_tau",
        2:   "mu- nu_tau",
        3:   "pi- nu_tau",
        4:   "K- nu_tau",
        145: "e- nu_tau gamma",
        246: "mu- nu_tau gamma",
        -1:  "no clasificado",
        -2:  "carga incorrecta",
        -3:  "recursion truncada",
        -4:  "contador > 9",
    }
    if code in special:
        return special[code]

    # Decodificar 5 dígitos: J Kc Kn Nc Nn
    s  = f"{code:05d}"
    J  = int(s[0])
    Kc = int(s[1])
    Kn = int(s[2])
    Nc = int(s[3])
    Nn = int(s[4])

    J_name = {
        0: "directo",
        1: "via rho-",
        2: "via rho0",
        3: "via a1->rho0 pi",
        4: "via a1->rho- pi0",
        5: "via a1 (no res.)",
        6: "via omega",
        7: "via K*",
        8: "via eta",
        9: "otro",
    }.get(J, f"J={J}")

    partes = []
    if Kc > 0: partes.append(f"{Kc}K±")
    if Kn > 0: partes.append(f"{Kn}K0")
    if Nc > 0: partes.append(f"{Nc}π±")
    if Nn > 0: partes.append(f"{Nn}π0")
    hadrones = " ".join(partes) if partes else "sin hadrones"

    return f"{J_name}, {hadrones}"

def sort_key(code):
    if code in (1, 2, 3, 4, 145, 246):
        return (0, code)
    if code < 0:
        return (99, code)
    return (1, code)

total = sum(pionic_counts.values())
print(f"\n{'Código':>8}  {'N':>8}  {'BR (%)':>8}   Descripción")
print("-" * 72)
for code, n in sorted(pionic_counts.items(), key=lambda x: sort_key(x[0])):
    br = 100.0 * n / total if total > 0 else 0
    print(f"{code_str(code):>8}  {n:>8}  {br:>7.3f}%   {describe(code)}")
print("-" * 72)
print(f"{'TOTAL':>8}  {total:>8}  100.000%")
