import sys,time,json; sys.path.insert(0,'/home/claude/b1')
from sympy import symbols, Poly
from sympy.polys.numberfields.galoisgroups import galois_group
from artin.run import run_pipeline
from artin.certificate import load_json
x=symbols('x')
cases={
 "q4_13sq61": [4,1,2,0,1],            # x^4+2x^2+x+4
 "q4_7_47sq": [-4,-3,-1,0,1],         # x^4-x^2-3x-4
 "q4_5_7sq_61": [4,1,1,0,1],          # x^4+x^2+x+4
 "s6_ord72_601": [1,0,-2,-2,0,0,1],   # x^6-2x^4-2x^3+1
 "s6_ord48_59sq": [1,0,0,0,2,0,1],    # x^6+2x^4+1
 "s6_ord18_13sq": [1,2,0,-2,-1,0,1],  # x^6-x^4-2x^3+2x+1
 "s6_ord24_23sq": [-1,0,0,0,1,0,1],   # x^6+x^4-1
 "d7_PSL27": [3,-7,0,0,0,0,0,1],      # x^7-7x+3, PSL(2,7)
}
which=sys.argv[1:]
for name,f in cases.items():
    if which and name not in which: continue
    P=Poly(list(reversed(f)),x); G,alt=galois_group(P)
    gens=[[int(g(i))+1 for i in range(len(f)-1)] for g in G.generators]
    d="/tmp/hard_"+name
    t=time.time()
    try:
        run_pipeline({"f":f,"generators":gens,"run_dir":d,"options":{"class_bound_X":150,"precision_check":False,"class_confirmation":False,"local_unramified_check":2,"quiet":True,"fe_test":False}})
        L=load_json(d+"/local.json"); CJ=load_json(d+"/conductors.json"); RN=load_json(d+"/rootnumbers.json"); T=load_json(d+"/timings.json")
        print(name,"|G|",G.order(),"ok",round(time.time()-t,1),"primes",{k:(v["status"][:4],v.get("e"),v.get("f")) for k,v in L["ramified"].items()},
              "conductors",[c["partial_conductor"] for c in CJ["conductors"]],"W",[c.get("W_exact") or (c.get("W_complex") and [round(c["W_complex"][0],3),round(c["W_complex"][1],3)]) for c in RN["characters"]][:6], "top stages",sorted(T.items(),key=lambda kv:-kv[1])[:2], flush=True)
    except Exception as ex:
        print(name,"|G|",G.order(),"FAILED",round(time.time()-t,1),str(ex)[:300], flush=True)
