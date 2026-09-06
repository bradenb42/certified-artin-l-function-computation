import sys,json,time; sys.path.insert(0,'/home/claude/b1')
from artin.check_splitting_field import check_polynomial, group_generators
sample = {
 "x^3-x^2-2x-8 (S3)": ([-8,-2,-1,1], ("S",3)),
 "x^3-x-1 (S3, disc -23)": ([-1,-1,0,1], ("S",3)),
 "x^3+x^2-2x-1 (C3, 7)": ([-1,-2,1,1], ("C",3)),
 "x^3-7x+14 (S3, 2 wild)": ([14,-7,0,1], ("S",3)),
 "x^4-x-1 (S4)": ([-1,-1,0,0,1], ("S",4)),
 "x^4+x^3+x^2+x+1 (C4)": ([1,1,1,1,1], ("C",4)),
 "x^4-2 (D4, 2 wild)": ([-2,0,0,0,1], ("D",4)),
 "(x^2-5)(x^2-13) (V4)": ([65,0,-18,0,1], ("V4",4)),
 "x^5-x-1 (S5)": ([-1,-1,0,0,0,1], ("S",5)),
 "x^5+x^4-4x^3-3x^2+3x+1 (C5)": ([1,3,-3,-4,1,1], ("C",5)),
 "x^5-5x+12 (D5, wild)": ([12,-5,0,0,0,1], ("D",5)),
 "x^5+20x+16": ([16,20,0,0,0,1], ("A",5)),
 "x^5-4x-7 (S5b)": ([-7,-4,0,0,0,1], ("S",5)),
 "x^5-2x-7 (S5c, 3 tame)": ([-7,-2,0,0,0,1], ("S",5)),
 "x^4+x+1 (S4b)": ([1,1,0,0,1], ("S",4)),
 "x^5-4x-7 (S5b)": ([-7,-4,0,0,0,1], ("S",5)),
 "x^5-2x-7 (S5c, 3 tame)": ([-7,-2,0,0,0,1], ("S",5)),
 "x^4+x+1 (S4b)": ([1,1,0,0,1], ("S",4)),
}
which=sys.argv[1:]
out={}
for name,(f,(gname,n)) in sample.items():
    if which and not any(w in name for w in which): continue
    r=check_polynomial(f, group_generators(gname,n), log=lambda m:None)
    out[name]=r
    print(name, "|G|",r.get("|G|"), "ok",r.get("all_compared_ok"), r.get("seconds"))
    for ell,e in r.get("primes",{}).items():
        print("   ",ell,"v",e["v_disc_f"], "ref",e["reference"] if isinstance(e["reference"],str) else {k:e["reference"][k] for k in ("status","e","f","g","v_disc_N")}, "pipe",e.get("pipeline"), "restricted" if e.get("restricted") else "", "ok",e.get("ok"), "sum f",e.get("sum_chi1_f_ell"))
json.dump(out,open("/tmp/check2_%s.json"%("_".join(which) or "all"),"w"),indent=1)
