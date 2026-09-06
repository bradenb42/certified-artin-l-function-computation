import sys,json,time; sys.path.insert(0,'/home/claude/b1')
from artin.check_newforms import check_dihedral
from artin.ramified import discriminant
from artin.chartable import is_prime
sample=[]
seen=set()
for a in range(-6,7):
    for b in range(-8,9):
        f=[b,a,0,1]; d=discriminant(f)
        if d<0 and is_prime(-d) and -d<=2000 and -d not in seen:
            seen.add(-d); sample.append((f,d))
sample.sort(key=lambda t:-t[1])
extra=[([-1,-1,0,0,1],-283)]   # S4 quartic: quadratic subfield Q(sqrt -283)
which=sys.argv[1:]
lo=int(which[0]) if which else 0; hi=int(which[1]) if len(which)>1 else 10**9
out=[]
for f,d in sample+extra:
    if not (lo<=-d<=hi): continue
    gens=[[2,3,1],[2,1,3]] if len(f)==4 else [[2,3,4,1],[2,1,3,4]]
    r=check_dihedral(f,gens,d,log=lambda m:None)
    out.append(r)
    ch=r.get("characters",[{}])[0]
    print(f, "D",d, r.get("all_ok"), "mism",ch.get("n_mismatches"), "level",ch.get("level_ok"), "neb",ch.get("nebentypus_ok"), "W",ch.get("W_pipeline"), r.get("seconds"), flush=True)
json.dump(out,open("/tmp/check3_%d_%d.json"%(lo,hi),"w"),indent=1)
