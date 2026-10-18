"""
L5_flagship_sim_v5b.py -- R4: the flagship scenario as an actual simulation (replaces the analytic L5).

8-node ring, s=0, d=4, two disjoint 4-hop arcs:  A = 0-1-2-3-4 (links 0,1,2,3)   B = 0-7-6-5-4 (links 7,6,5,4)
Declared class assignment: ALL links class B (S0 = 2.687, eps = 1e-6, gamma = 8e-7).  Delta = 5e-6 (main) / 1e-6 (tight).
60 epochs; fresh audit of n = 10^4 pairs per link per epoch, exact Binomial(4n, (1+S/4)/2) noise.
From epoch 25 on, link 1 (1-2, on arc A) is compromised: true S = 2.418 (S_lb ~ 2.19 >= 2 -> the GATE passes; r_e = 0.079 < R = 0.10;
the ALARM fires).  500 trials.  Routers: the R1 router set (B_core.select_all), re-selected every epoch from that epoch's audit.
Tie-breaking: deterministic (arc A first) = a router with a fixed tie-break; --tiebreak random = random.
Outputs: ../data/R4_*.csv/json, ../figures/fig_R4_trace.png
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import B_core as B

N_AUDIT = 10_000; EPOCHS = 60; ATTACK_EPOCH = 25; ATT_LINK = 1; R = 0.10
S_B = 2.687; EPS_B = 1e-6; GAM_B = 8e-7
PATHS = np.array([[0, 1, 2, 3], [7, 6, 5, 4]], np.int32)       # padded link matrix (no -1 needed: both 4 hops)
HOPS = np.array([4, 4])
ROUTERS = ["BMR", "SPF", "SPF_valid", "SPF_budget", "SPF_valid_budget", "BMR_hop",
           "QPATHstyle", "QLEAPstyle", "QCASTstyle"]


class DetRng:
    """Deterministic tie-break: lowest index (arc A) first."""
    def random(self, k):
        return np.arange(k, dtype=float) / max(k, 1)


def epoch_links(rng, attacked):
    E = 8
    S_true = np.full(E, S_B)
    if attacked:
        S_true[ATT_LINK] = B.S_ATTACK
    p = (1.0 + S_true / 4.0) / 2.0
    K = rng.binomial(4 * N_AUDIT, p)
    S_hat = 2.0 * K / N_AUDIT - 4.0
    t = B.tau(N_AUDIT)
    S_lb = np.maximum(0.0, S_hat - t)
    r_e = B.h_vec(S_lb)
    gate = S_lb >= 2.0
    alarm = S_hat < S_B - B.Z_ALPHA * B.sigma_w(S_B, N_AUDIT)
    fresh = np.ones(E, bool)
    valid = gate & ~alarm & fresh
    v = np.clip(S_hat / (2 * B.SQ2), 0.0, 1.0)
    return dict(r=r_e, eps=np.full(E, EPS_B), gamma=np.full(E, GAM_B), S_hat=S_hat, S_lb=S_lb, gate=gate,
                alarm=alarm, fresh=fresh, valid=valid, comp=np.zeros(E, bool), v=v, pgen=np.maximum(0.01, v),
                S_true=S_true)


def run(a):
    ss = np.random.SeedSequence([a.seed, int(a.Delta * 1e9), 0 if a.tiebreak == "det" else 1])
    rng = np.random.default_rng(ss)
    tb = DetRng() if a.tiebreak == "det" else np.random.default_rng(ss.spawn(1)[0])
    rows = []; det = []; n_clean_dec = 0; fa_events = []; fa_pre = 0; fa_post_other = 0
    S_att_hat = []; Slb_att = []; r_att = []; gate_att = []
    for trial in range(a.trials):
        detect = None; clean_fa = 0
        for ep in range(1, EPOCHS + 1):
            attacked = ep >= ATTACK_EPOCH
            lk = epoch_links(rng, attacked)
            # clean link-epoch decisions: every link before the attack, all links except the attacked one afterwards
            clean = np.ones(8, bool)
            if attacked: clean[ATT_LINK] = False
            n_clean_dec += int(clean.sum())
            fa = lk["alarm"] & clean
            if fa.any():
                clean_fa += int(fa.sum())
                if attacked: fa_post_other += int(fa.sum())
                else: fa_pre += int(fa.sum())
                fa_events.append((trial, ep, np.flatnonzero(fa).tolist()))
            if attacked:
                S_att_hat.append(lk["S_hat"][ATT_LINK]); Slb_att.append(lk["S_lb"][ATT_LINK])
                r_att.append(lk["r"][ATT_LINK]); gate_att.append(bool(lk["gate"][ATT_LINK]))
                if detect is None and lk["alarm"][ATT_LINK]:
                    detect = ep
            M = B.path_metrics(PATHS, lk)
            sel, M, F, tput = B.select_all(PATHS, HOPS, lk, a.Delta, R, "sealed", tb, ROUTERS, M)
            for r_ in B.trial_rows(sel, M, F, tput, PATHS, HOPS, a.Delta, "sealed", 1):
                rows.append((trial, ep, r_["router"], r_["feasible"], int(r_["path"].startswith("0-1-2-3")),
                             r_["r_nom"], r_["r_sec"], r_["alarm_on_path"]))
        det.append(detect)
    df = pd.DataFrame(rows, columns=["trial", "epoch", "router", "feasible", "on_A", "r_nom", "r_sec", "alarm_on_path"])
    return df, det, dict(clean_decisions=n_clean_dec, fa_pre_attack=fa_pre, fa_post_attack_other_links=fa_post_other,
                         fa_events=fa_events[:20], S_att_hat_mean=float(np.mean(S_att_hat)), S_att_hat_sd=float(np.std(S_att_hat)),
                         S_lb_att_mean=float(np.mean(Slb_att)), S_lb_att_min=float(np.min(Slb_att)),
                         S_lb_att_max=float(np.max(Slb_att)), r_att_mean=float(np.mean(r_att)), r_att_max=float(np.max(r_att)),
                         gate_pass_rate=float(np.mean(gate_att)))


def summarise(df, det, info, a):
    out = dict(Delta=a.Delta, trials=a.trials, tiebreak=a.tiebreak, seed=a.seed, n_audit=N_AUDIT)
    out.update({k: v for k, v in info.items() if k != "fa_events"})
    out["fa_events_first20"] = info["fa_events"]
    out["expected_fa_at_alpha_1e-6"] = info["clean_decisions"] * 1e-6
    d_ = [x for x in det if x is not None]
    out["detected"] = len(d_); out["undetected"] = len(det) - len(d_)
    out["detect_epoch_hist"] = {str(k): int(v) for k, v in pd.Series(d_).value_counts().sort_index().items()}
    pre = df[df.epoch < ATTACK_EPOCH]; post = df[df.epoch > ATTACK_EPOCH]; at = df[df.epoch == ATTACK_EPOCH]
    tab = {}
    for rt in ROUTERS:
        rec = {}
        for nm, sub in (("pre_1_24", pre), ("attack_25", at), ("post_26_60", post)):
            s = sub[sub.router == rt]
            f = s.feasible.mean()
            rec[nm] = dict(feasible=float(f), r_nom_given_feasible=float(s.r_nom[s.feasible == 1].mean()) if f > 0 else 0.0,
                           r_sec=float(s.r_sec.mean()), on_attacked_arc=float(s.on_A[s.feasible == 1].mean()) if f > 0 else 0.0,
                           alarm_on_path=float(s.alarm_on_path.mean()))
        # per-trial means -> CI over trials for post-attack r_sec
        pt = post[post.router == rt].groupby("trial").r_sec.mean().values
        rec["post_r_sec_ci95"] = [float(pt.mean() - 1.96 * pt.std(ddof=1) / np.sqrt(len(pt))), float(pt.mean() + 1.96 * pt.std(ddof=1) / np.sqrt(len(pt)))]
        tab[rt] = rec
    out["table"] = tab
    # epochs-to-recover for each router: first epoch >= 25 where router's path avoids the alarm and r_sec>0 (per trial), then mean
    rec_ep = {}
    for rt in ROUTERS:
        s = df[(df.router == rt) & (df.epoch >= ATTACK_EPOCH)]
        ok = s[s.r_sec > 0].groupby("trial").epoch.min()
        rec_ep[rt] = dict(recovered_trials=int(len(ok)), of=int(a.trials), mean_first_secure_epoch=float(ok.mean()) if len(ok) else None)
    out["recovery"] = rec_ep
    return out


def make_figure(dfs, outpng):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True)
    for k, (title, df) in enumerate(dfs):
        g = df.groupby(["router", "epoch"]).r_sec.mean().unstack(0)
        for rt, st in [("BMR", "-"), ("SPF_valid_budget", "--"), ("SPF_valid", "--"), ("SPF", ":"), ("SPF_budget", ":"),
                       ("QPATHstyle", "-."), ("QLEAPstyle", "-."), ("QCASTstyle", "-.")]:
            ax[k].plot(g.index, g[rt], st, label=rt, lw=1.8 if rt == "BMR" else 1.2)
        ax[k].axvline(ATTACK_EPOCH, color="k", lw=0.8); ax[k].text(ATTACK_EPOCH + .3, ax[k].get_ylim()[1] * .92, "attack", fontsize=8)
        ax[k].set_title(title, fontsize=10); ax[k].set_xlabel("epoch"); ax[k].grid(alpha=.3)
    ax[0].set_ylabel("mean r_sec (500 trials)"); ax[0].legend(fontsize=7, ncol=2)
    plt.tight_layout(); plt.savefig(outpng, dpi=150)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=500)
    ap.add_argument("--Delta", type=float, default=5e-6)
    ap.add_argument("--tiebreak", choices=["det", "random"], default="det")
    ap.add_argument("--seed", type=int, default=20261012)
    ap.add_argument("--tag", default="")
    ap.add_argument("--outdir", default="../data")
    a = ap.parse_args()
    t0 = time.time()
    df, det, info = run(a)
    out = summarise(df, det, info, a)
    tag = a.tag or f"D{a.Delta:g}_{a.tiebreak}"
    out["runtime_s"] = time.time() - t0
    json.dump(out, open(f"{a.outdir}/R4_summary_{tag}.json", "w"), indent=1)
    df.groupby(["router", "epoch"]).agg(feasible=("feasible", "mean"), r_nom=("r_nom", "mean"), r_sec=("r_sec", "mean"),
                                        on_A=("on_A", "mean")).reset_index().to_csv(f"{a.outdir}/R4_trace_{tag}.csv", index=False)
    print(json.dumps({k: out[k] for k in ["Delta", "clean_decisions", "fa_pre_attack", "fa_post_attack_other_links", "expected_fa_at_alpha_1e-6",
                                          "detected", "detect_epoch_hist", "S_lb_att_mean", "r_att_mean", "gate_pass_rate"]}, indent=1))
    for rt, r in out["table"].items():
        print(f"{rt:18s} pre: feas={r['pre_1_24']['feasible']:.3f} rsec={r['pre_1_24']['r_sec']:.4f} | "
              f"post: feas={r['post_26_60']['feasible']:.3f} rnom|f={r['post_26_60']['r_nom_given_feasible']:.4f} rsec={r['post_26_60']['r_sec']:.4f} onA={r['post_26_60']['on_attacked_arc']:.2f}")
