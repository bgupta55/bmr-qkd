"""
make_table2b.py -- regenerates Table 2b (reference values) of the BMR manuscript.

Usage:  python3 make_table2b.py            (prints Markdown, writes ../data/R11_table2b.{md,json})
        python3 make_table2b.py --check    (also asserts every cell against the corrected reference rows)

Definitions (all from the manuscript, Sec. 3.2.1 / 4.2.1):
  eta   = sqrt(32 ln(8/delta) / n)           Hoeffding radius, union bound over 4 correlators x 2 tails
  S_lb  = S_hat - eta                        evaluated at the expected S_hat (S0 honest, S_dag compromised)
  r_e   = h(S_lb) if S_lb >= 2 else 0        h(S) = 1 - log2(1 + sqrt(2 - S^2/4))
  gate  = (S_lb >= 2)
  floor = (r_e >= R),  R = 0.10  (S_lb threshold 2.2360)
  theta = S0 - z_alpha * sigma_W(S0, n),  sigma_W(S, n) = sqrt(4 - S^2/4)/sqrt(n)   (binomial / Werner noise)
  beta  = P(S_hat < theta | S_true) = Phi((theta - S_true)/sigma_W(S_true, n))        (detection power)

Corrections relative to the withdrawn v5 Table 2b:
  * the compromised link S_dag = 2.418 PASSES the gate at n = 10^4 (S_lb = 2.1925 >= 2); r_e = 0.0790, not 0.000
  * detection is by the alarm (and the request floor), not by the gate
  * power uses the Werner sigma only (no mixing with sigma = 0.02): beta(2.418) = 1.000
"""
import json, math, os, sys
from statistics import NormalDist

DELTA = 1e-6
S0, S_DAG = 2.687, 2.418
R_FLOOR = 0.10
ALPHA = 1e-6
PHI = NormalDist().cdf
Z_ALPHA = NormalDist().inv_cdf(1 - ALPHA)          # 4.7534
N_LIST = [1000, 1078, 2500, 5000, 10000, 20000, 40000]


def eta(n, delta=DELTA):
    return math.sqrt(32.0 * math.log(8.0 / delta) / n)


def h(S):
    if S < 2.0:
        return 0.0
    return 1.0 - math.log2(1.0 + math.sqrt(max(0.0, 2.0 - S * S / 4.0)))


def sigma_w(S, n):
    return math.sqrt(4.0 - S * S / 4.0) / math.sqrt(n)


def theta(n):
    return S0 - Z_ALPHA * sigma_w(S0, n)


def beta(S_true, n):
    return PHI((theta(n) - S_true) / sigma_w(S_true, n))


def min_n_honest(S, delta=DELTA):
    """smallest integer n with S - eta(n) >= 2."""
    n = math.ceil(32.0 * math.log(8.0 / delta) / (S - 2.0) ** 2)
    while S - eta(n, delta) < 2.0:
        n += 1
    return n


def row(n):
    e = eta(n)
    slh, slc = S0 - e, S_DAG - e
    rh, rc = h(slh), h(slc)
    return dict(n=n, four_n=4 * n, eta=e, S_lb_honest=slh, r_honest=rh, S_lb_comp=slc, r_comp=rc,
                gate_honest=slh >= 2.0, gate_comp=slc >= 2.0, floor_comp=rc >= R_FLOOR,
                theta_alpha=theta(n), beta_comp=beta(S_DAG, n), beta_2p62=beta(2.62, n))


# Corrected reference rows (BLUEPRINT v5 S21): n, eta, S_lb_h, r_h, S_lb_c, r_c, theta
REF = {
    1000:  (0.713,  1.974, 0.0,     1.705, 0.0,    None),
    1078:  (0.687,  2.000, 0.00004, 1.731, 0.0,    None),
    2500:  (0.451,  2.236, 0.09997, 1.967, 0.0,    2.5462),
    5000:  (0.319,  2.368, 0.1735,  2.099, 0.0381, 2.5874),
    10000: (0.2255, 2.4615, 0.2373, 2.1925, 0.0790, 2.6166),
    20000: (0.1595, 2.5275, 0.2909, 2.2585, 0.1114, 2.6372),
    40000: (0.1128, 2.5742, 0.3347, 2.3052, 0.1365, 2.6518),
}


def main(check=False):
    rows = [row(n) for n in N_LIST]
    out = ["| n (per setting) | 4n | eta | S_lb honest | r_e honest | S_lb comp. | r_e comp. | Gate (comp.) | r_e >= R (comp.) | theta_alpha | beta(2.418) |",
           "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        gate_h = "honest fails" if not r["gate_honest"] and r["n"] == 1000 else ("min n (honest)" if r["n"] == 1078 else "")
        gate_c = "passes" if r["gate_comp"] else "**fails**"
        if r["n"] == 1000:
            gate_c = "honest fails"
        if r["n"] == 1078:
            gate_c = "min n (honest)"
        fl = "-" if r["n"] in (1000, 1078) else ("**yes**" if r["floor_comp"] else "no")
        th = "-" if r["n"] in (1000, 1078) else f"{r['theta_alpha']:.4f}"
        be = "-" if r["n"] in (1000, 1078) else f"{r['beta_comp']:.3f}"
        out.append(f"| {r['n']:,} | {r['four_n']:,} | {r['eta']:.4f} | {r['S_lb_honest']:.4f} | {r['r_honest']:.5f} | "
                   f"{r['S_lb_comp']:.4f} | {r['r_comp']:.4f} | {gate_c} | {fl} | {th} | {be} |")
    md = "\n".join(out)
    extra = dict(
        z_alpha=Z_ALPHA,
        S_lb_threshold_R=None,
        min_n_honest=min_n_honest(S0),
        eta_n100=eta(100), S_lb_honest_n100=S0 - eta(100),
        r_e_compromised_n1e4=h(S_DAG - eta(10000)),
        sigma_W_S0_n1e4=sigma_w(S0, 10000),
        sigma_W_Sdag_n1e4=sigma_w(S_DAG, 10000),
        beta_2p418_n1e4=beta(S_DAG, 10000),
        beta_2p62_n1e4_sigmaW_at_S_true=beta(2.62, 10000),
        beta_2p62_n1e4_sigmaW_at_S0=PHI((theta(10000) - 2.62) / sigma_w(S0, 10000)),
        beta_2p62_n1e4_sigma_0p02=PHI((theta(10000) - 2.62) / 0.02),
    )
    # S_lb at which h(S_lb) = R (bisection)
    lo, hi = 2.0, 2.8
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if h(mid) < R_FLOOR else (lo, mid)
    extra["S_lb_threshold_R"] = 0.5 * (lo + hi)
    print(md)
    print()
    for k, v in extra.items():
        print(f"{k:34s} {v:.6f}" if isinstance(v, float) else f"{k:34s} {v}")

    here = os.path.dirname(os.path.abspath(__file__))
    ddir = os.path.join(here, "..", "data")
    os.makedirs(ddir, exist_ok=True)
    with open(os.path.join(ddir, "R11_table2b.md"), "w") as f:
        f.write(md + "\n")
    with open(os.path.join(ddir, "R11_table2b.json"), "w") as f:
        json.dump(dict(rows=rows, extra=extra), f, indent=1)

    if check:
        bad = 0
        for r in rows:
            e_ref, slh, rh, slc, rc, th = REF[r["n"]]
            tests = [("eta", r["eta"], e_ref, 6e-4), ("S_lb_h", r["S_lb_honest"], slh, 6e-4),
                     ("r_h", r["r_honest"], rh, 6e-5), ("S_lb_c", r["S_lb_comp"], slc, 6e-4),
                     ("r_c", r["r_comp"], rc, 6e-5)]
            if th is not None:
                tests.append(("theta", r["theta_alpha"], th, 6e-4))
            for name, got, ref, tol in tests:
                ok = abs(got - ref) <= tol
                bad += (not ok)
                if not ok:
                    print(f"  MISMATCH n={r['n']} {name}: got {got:.5f} ref {ref:.5f}")
        c = row(10000)
        corrected = [("S_dag passes the gate (S_lb=2.1925)", c["gate_comp"] and abs(c["S_lb_comp"] - 2.1925) < 6e-4),
                     ("r_e(S_dag) = 0.0790", abs(c["r_comp"] - 0.0790) < 6e-5),
                     ("r_e(S_dag) < R (request floor removes it)", c["r_comp"] < R_FLOOR),
                     ("beta(2.418) = 1.000 (Werner sigma)", abs(c["beta_comp"] - 1.0) < 5e-4),
                     ("n=100: eta=2.255, S_lb=0.432", abs(eta(100) - 2.255) < 1e-3 and abs(S0 - eta(100) - 0.432) < 1e-3),
                     ("min honest n = 1078", min_n_honest(S0) == 1078)]
        for name, ok in corrected:
            print(("  [PASS] " if ok else "  [FAIL] ") + name)
            bad += (not ok)
        print("TABLE 2b CHECK:", "ALL PASS" if bad == 0 else f"{bad} FAILURES")
        sys.exit(0 if bad == 0 else 1)


if __name__ == "__main__":
    main(check="--check" in sys.argv)
