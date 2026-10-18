"""Exact expectation of SPF R_sec over the uniform tie-break among minimum-hop paths (p_exp=0, Delta=5e-6, H=5), same seeds as R1/R5."""
import numpy as np, sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import B_core as B
SEED=20261012; res={"clean":[], "compromised":[]}
for g in range(50):
    rg=np.random.default_rng(np.random.SeedSequence([SEED,30,1,g,0])); edges=B.build_graph("er",30,rg); E=len(edges)
    adj=B.adj_matrix(30,edges); L,hops=B.enumerate_paths(adj,0,29,3,5)
    if len(hops)==0: 
        for k in res: res[k].extend([0.0]*5)
        continue
    for draw in range(5):
        for f,cond in ((0.0,"clean"),(0.10,"compromised")):
            rl=np.random.default_rng(np.random.SeedSequence([SEED,30,1,g,draw,1])); rc=np.random.default_rng(np.random.SeedSequence([SEED,30,1,g,draw,2,int(f*1000)]))
            lk=B.make_links(rl,E,10000,f=f,noise="binomial",p_exp=0.0,rng_comp=rc); M=B.path_metrics(L,lk)
            m=hops==hops.min(); ok=M["allvalid"]&(M["seps"]<=5e-6*(1+1e-9))
            res[cond].append(float(np.mean(np.where(ok[m],M["rmin"][m],0.0))))
out={k:dict(mean=float(np.mean(v)),n=len(v)) for k,v in res.items()}
json.dump(out,open("../data/R5_spf_exact.json","w"),indent=1); print(out)
