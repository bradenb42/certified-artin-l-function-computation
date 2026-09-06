import sys,time; sys.path.insert(0,'/home/claude/b1')
from artin.run import run_pipeline
from artin.certificate import load_json
# PSL(2,7) on 7 points: generators (1234567) and (2 3 5)(4 7 6)? we use the standard: a=(0..6), b=(1,2,4)(3,6,5) gives 21; use instead x->x+1 and x->2x on P^1(F_7)? that's 42. PSL(2,7) acting on the 7 points of the Fano plane: generators
a=[2,3,4,5,6,7,1]; b=[1,2,4,3,7,6,5]  # b = (2 3)(4 5)(6 7): with the 7-cycle generates PSL(2,7) (order 168)
from artin.perm import PermGroup, from_json
G=PermGroup([from_json(a),from_json(b)]); print("order",G.order())
t=time.time()
f=[3,-7,0,0,0,0,0,1]
run_pipeline({"f":f,"generators":[a,b],"run_dir":"/tmp/hard_psl27","options":{"class_bound_X":100,"precision_check":False,"class_confirmation":False,"local_unramified_check":1,"quiet":True,"fe_test":False,"compute_models":False}})
print("done",time.time()-t)
