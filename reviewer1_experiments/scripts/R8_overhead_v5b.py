"""
R8_overhead_v5b.py -- R8: overhead table with SIGNED certificates, both generators (ER p=0.25 and p=0.40), n = 10..150.
Certificate layout (L9): 6 x float64 + int32 + uint8 = 53 B payload; raw ECDSA-P256 signature r||s = 64 B  ->  117 B signed.
Sensitivity: if the freshness stamp tau_e is stored as int32 (4 B) the payload is 49 B and the signed certificate 113 B.
Edge counts: (i) the expected count p*n(n-1)/2 (rounded as in L9 for p=0.4), (ii) the realised mean over 30 seeded graphs.
Timings: re-measured here (2000 ECDSA-P256 sign / verify ops, cryptography library) and compared with the stored L9 values.
Epoch = 4*n_audit/pair_rate = 4*1e4/1e5 = 0.4 s.  MTU payload 1472 B.
"""
import hashlib, json, math, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import B_core as B
from cryptography.hazmat.primitives.asymmetric.ec import SECP256R1, generate_private_key, ECDSA
from cryptography.hazmat.primitives.asymmetric.utils import decode_dss_signature
from cryptography.hazmat.primitives.hashes import SHA256

PAYLOAD = 6 * 8 + 4 + 1; SIG = 64; SIGNED = PAYLOAD + SIG; MTU = 1472; EPOCH_S = 0.4
key = generate_private_key(SECP256R1()); pub = key.public_key(); msg = bytes(PAYLOAD)
st, vt = [], []
for _ in range(2000):
    t = time.perf_counter(); sig = key.sign(msg, ECDSA(SHA256())); st.append((time.perf_counter() - t) * 1e6)
for _ in range(2000):
    t = time.perf_counter(); pub.verify(sig, msg, ECDSA(SHA256())); vt.append((time.perf_counter() - t) * 1e6)
SIGN_US, VERIFY_US = float(np.median(st)), float(np.median(vt))
stored = json.load(open("../data/L9_overhead_results.json"))["timings"]
# the table uses the stored L9 values (18.0 / 40.17 us) so that it is the same machine/run as the submitted overhead numbers;
# the re-measurement is reported next to it.
SIGN_T, VERIFY_T = stored["ecdsa_sign_median_us"], stored["ecdsa_verify_median_us"]
print(f"re-measured: sign {SIGN_US:.1f} us, verify {VERIFY_US:.1f} us; stored: {SIGN_T} / {VERIFY_T}")
print("DER signature length (bytes) example:", len(sig), "(raw r||s = 64)")

rows = []
for fam, p in (("ER p=0.25", 0.25), ("ER p=0.40", 0.40)):
    for n in (10, 20, 30, 50, 75, 100, 150):
        E_exp = int(n * (n - 1) / 2 * p)
        real = []
        for g in range(30):
            rg = np.random.default_rng(np.random.SeedSequence([8000 + g, n, int(p * 100)]))
            iu, ju = np.triu_indices(n, 1); real.append(int((rg.random(len(iu)) < p).sum()))
        E = E_exp
        db = E * SIGNED
        rows.append(dict(family=fam, n=n, E_expected=E_exp, E_realised_mean=float(np.mean(real)), E_realised_sd=float(np.std(real, ddof=1)),
                         db_unsigned_B=E * PAYLOAD, db_signed_B=db, db_signed_int32stamp_B=E * (SIGNED - 4), datagrams=math.ceil(db / MTU),
                         sign_parallel_us=SIGN_T, sign_sequential_ms=E * SIGN_T / 1e3, pct_epoch_sequential=100 * E * SIGN_T / 1e6 / EPOCH_S,
                         pct_epoch_parallel=100 * SIGN_T / 1e6 / EPOCH_S))
        print(rows[-1])
# ledger / verification: hops from R1 (ER n=30 BMR mean hops 4.2) and hop range 3..5
led = {h: dict(ledger_B=h * SIGNED, verify_us=h * VERIFY_T) for h in (3, 4, 5)}
# BA m=2 for reference
ba = []
for n in (10, 20, 30, 50, 75, 100, 150):
    E = (3 * 2 // 2) + (n - 3) * 2
    ba.append(dict(n=n, E=E, db_signed_B=E * SIGNED, datagrams=math.ceil(E * SIGNED / MTU), sign_sequential_ms=E * SIGN_T / 1e3,
                   pct_epoch_sequential=100 * E * SIGN_T / 1e6 / EPOCH_S))
out = dict(payload_B=PAYLOAD, sig_B=SIG, signed_B=SIGNED, mtu=MTU, epoch_s=EPOCH_S, sign_us_stored=SIGN_T, verify_us_stored=VERIFY_T,
           sign_us_remeasured=SIGN_US, verify_us_remeasured=VERIFY_US, rows=rows, ledger=led, ba_reference=ba,
           expected_check={"n20_p0.4": 8892, "n50_p0.4": 57330, "n100_p0.4": 231660, "n150_p0.4": 522990})
json.dump(out, open("../data/R8_overhead.json", "w"), indent=1)
md = ["### T-R8-A  Signed-certificate overhead (117 B = 53 B payload + 64 B raw ECDSA-P256 signature); epoch = 0.4 s (n_audit = 10⁴, 10⁵ pairs/s)\n",
      f"Timings (cryptography library, Apple arm64, median of 2000 ops): stored L9 values sign {SIGN_T} µs / verify {VERIFY_T} µs; re-measured today {SIGN_US:.1f} / {VERIFY_US:.1f} µs.\n",
      "| generator | n | links E (expected) | E realised (mean ± SD, 30 graphs) | unsigned DB (B) | **signed DB (B)** | datagrams (1472 B) | sequential signing of all E (ms) | % of 0.4 s epoch (sequential) | % of epoch (per-link parallel, 18 µs) |",
      "|---|---|---|---|---|---|---|---|---|---|"]
for r in rows:
    md.append(f"| {r['family']} | {r['n']} | {r['E_expected']} | {r['E_realised_mean']:.1f} ± {r['E_realised_sd']:.1f} | {r['db_unsigned_B']:,} | **{r['db_signed_B']:,}** | {r['datagrams']} | {r['sign_sequential_ms']:.2f} | {r['pct_epoch_sequential']:.2f} | {r['pct_epoch_parallel']:.4f} |")
md += ["\n### T-R8-B  Per-route ledger and destination verification (signed 117 B per hop; verify 40.17 µs per certificate)\n",
       "| hops | ledger (B) | verification at destination (µs) |\n|---|---|---|"]
for h, v in led.items(): md.append(f"| {h} | {v['ledger_B']} | {v['verify_us']:.1f} |")
md.append(f"\nMean BMR route length in R1 (ER n=30, Δ=5×10⁻⁶): 4.18 hops → ledger ≈ {4.18*SIGNED:.0f} B, verification ≈ {4.18*VERIFY_T:.0f} µs.")
md += ["\n### T-R8-C  Sensitivity and reference\n",
       "| item | value |\n|---|---|",
       f"| signed certificate if τ_e is stored as int32 (49 B payload) | {SIGNED-4} B → ER p=0.4, n=150: {4470*(SIGNED-4):,} B |",
       f"| DER-encoded signature length (observed) | {len(sig)} B (raw r‖s = 64 B used above) |"]
md.append("| BA (m=2) signed DB, n = 50 / 100 / 150 | " + " / ".join(f"{b['db_signed_B']:,} B ({b['datagrams']} dg, {b['pct_epoch_sequential']:.2f}% epoch)" for b in ba if b['n'] in (50, 100, 150)) + " |")
open("../data/R8_tables.md", "w").write("\n".join(md) + "\n"); print("\n".join(md))
