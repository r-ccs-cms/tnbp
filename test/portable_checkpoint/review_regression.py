"""Gate representability and estimator-owned atomic publication, CPU or CUDA.
The launcher may be mpiexec or an srun adapter. No production fault switches.
"""
from pathlib import Path
import hashlib,json,os,re,shutil,subprocess,sys,tempfile,time
mpi,estimator,fixture,state_fixture,evidence=sys.argv[1:]
Path(evidence).mkdir(parents=True,exist_ok=True)
root=Path(tempfile.mkdtemp(prefix='accepted-',dir=evidence));cases=[]
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
def command(n,binary,args):return [mpi,'--oversubscribe','-np',str(n),binary,*map(str,args)]
def record(label,cmd,p,error=None):
 (root/(label+'.log')).write_text('COMMAND '+' '.join(cmd)+'\n'+p.stdout)
 cases.append({'label':label,'command':cmd,'exit':p.returncode,'expected_error':error})
 (root/'cases.json').write_text(json.dumps(cases,indent=2)+'\n')
 assert (p.returncode==0 if error is None else p.returncode!=0 and (error is True or error in p.stdout)),(label,p.stdout)
 return p.stdout
def run(label,n,binary,*args,error=None):
 cmd=command(n,binary,args)
 return record(label,cmd,subprocess.run(cmd,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=90),error)
def qasm(name,body):
 path=root/name;path.write_text('OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[4];\n'+body+'\nbarrier q;\n');return path
def observable(name,text):
 p=root/name;p.write_text('1.0 0.0 '+text+'\n');return p
def measure(text,value):
 values=re.findall(r'^ Expectation value of single site operator.* (\S+)$',text,re.M)
 assert len(values)==1 and abs(float(values[0])-value)<1e-8,(values,value)
def payloads(path):
 return {str(p.relative_to(path)):hashlib.sha256(p.read_bytes()).hexdigest() for sub in ('vertices','messages') for p in (path/sub).glob('*')}
prep=qasm('prepare.qasm','x q[0];\ncz q[0],q[1];\ncz q[2],q[3];')
z0=observable('z0.txt','IIIZ');z1=observable('z1.txt','IIZI');z2=observable('z2.txt','IZII')
saved=root/'disconnected'
run('prepare',1,estimator,'--circuit',prep,'--sparse_pauli',z0,'--step_end',1,'--portable_save',saved)
original=payloads(saved)
for gate,body in [('cx','cx q[0],q[2];'),('ccx','ccx q[0],q[1],q[2];')]:
 circuit=qasm(gate+'.qasm',body)
 dest=root/(gate+'-invalid-output')
 run(gate+'-unrepresentable',2,estimator,'--circuit',circuit,'--sparse_pauli',z2,
     '--portable_load',saved,'--step_end',1,'--portable_save',dest,error='cannot be represented')
 assert not dest.exists() and payloads(saved)==original
# Disconnected graphs themselves are valid, including barriers over all sites.
for label,body,obs,expect in [('within-component','cx q[0],q[1];',z1,-1),('single-site','x q[0];',z0,1)]:
 circuit=qasm(label+'.qasm',body)
 text=run(label,2,estimator,'--circuit',circuit,'--sparse_pauli',obs,
          '--portable_load',saved,'--step_end',1,'--measurement_barrier',0)
 measure(text,expect)
# Observe the actual wrapper while it is paused after the generic state save.
identity=qasm('identity.qasm','id q[0];\ncz q[0],q[1];\ncz q[1],q[2];\ncz q[2],q[3];')
for n in (1,2):
 target=root/f'published-{n}';staging=Path(str(target)+'.estimator.pending')
 cmd=command(n,fixture,['pause',target,identity]);proc=subprocess.Popen(cmd,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 try:
  deadline=time.monotonic()+60
  while not (staging/'ready').exists():
   if proc.poll() is not None:raise AssertionError(proc.communicate()[0])
   assert time.monotonic()<deadline,'waiting for publication pause'
   time.sleep(.05)
  assert (staging/'state/COMMITTED').exists()
  assert not (staging/'state/estimator.resume').exists() and not target.exists()
  (staging/'release').write_text('continue\n')
  stdout,_=proc.communicate(timeout=60)
  record(f'publish-{n}',cmd,subprocess.CompletedProcess(cmd,proc.returncode,stdout))
 finally:
  if proc.poll() is None:proc.kill();proc.communicate()
 assert target.is_dir() and (target/'estimator.resume').is_file() and not staging.exists()
 copy=root/f'published-copy-{n}'
 text=run(f'resume-{n}',2 if n==1 else 1,estimator,'--circuit',identity,'--sparse_pauli',z0,
          '--portable_load',target,'--step_end',3,'--portable_save',copy)
 assert 'PORTABLE_RESTART circuit=same step_start=3' in text
 assert payloads(target)==payloads(copy)
 # Both failure points must leave the final name unpublished.
 for mode in ('fail-after-state','fail-write'):
  target_bad=root/f'{mode}-{n}';pending=Path(str(target_bad)+'.estimator.pending')
  run(f'{mode}-{n}',n,fixture,mode,target_bad,identity,error=True)
  assert not target_bad.exists() and (pending/'state/COMMITTED').exists()
  assert not (pending/'state/estimator.resume').exists()
# A failure after metadata write must not overwrite an independently created destination.
blocked=root/'blocked';run('blocked-publication',2,fixture,'block-destination',blocked,identity,error='destination already exists')
assert (blocked/'sentinel').read_text()=='preserve existing destination\n'
assert not (blocked/'COMMITTED').exists()
assert (Path(str(blocked)+'.estimator.pending')/'state/estimator.resume').exists()
# Old partially published app files must not silently trigger a step-zero restart.
partial=root/'old-partial';shutil.copytree(root/'published-1',partial)
(partial/'estimator.resume').rename(partial/'estimator.resume.pending')
run('reject-old-pending',2,estimator,'--circuit',identity,'--sparse_pauli',z0,
    '--portable_load',partial,'--step_end',3,error='incomplete estimator resume file')
run('library-ignores-app-pending',2,state_fixture,'vertices-only',partial)
# A genuine library-only state remains valid and starts at zero in estimator.
bare=root/'bare';run('bare-library-copy',2,state_fixture,'load',root/'published-1',bare)
text=run('bare-estimator',2,estimator,'--circuit',identity,'--sparse_pauli',z0,'--portable_load',bare,'--step_end',1)
assert 'PORTABLE_RESTART circuit=unspecified step_start=0' in text
# Nonadjacent gates in a connected graph remain representable.
routed=qasm('routed.qasm','x q[0];\ncx q[0],q[2];')
text=run('connected-routing',2,estimator,'--circuit',routed,'--sparse_pauli',z2,
         '--portable_load',root/'published-1','--step_end',1,'--measurement_barrier',0)
measure(text,-1)
# Public string paths preserve spaces and accept trailing directory separators.
spaced=root/'state with spaces'
run('string-path-save',2,estimator,'--circuit',identity,'--sparse_pauli',z0,
    '--portable_save',str(spaced)+'/', '--step_end',1)
copy=root/'string-copy'
text=run('string-path-load',1,estimator,'--circuit',identity,'--sparse_pauli',z0,
         '--portable_load',str(spaced)+'/', '--step_end',1,'--portable_save',copy)
assert 'PORTABLE_RESTART circuit=same step_start=1' in text
assert payloads(spaced)==payloads(copy)
bare_spaced=root/'bare state with spaces'
run('library-string-save',2,state_fixture,'load',spaced,str(bare_spaced)+'/')
run('library-string-load',1,state_fixture,'load',str(bare_spaced)+'/')
(root/'PASS').write_text(f'{len(cases)} review regression cases passed\n')
print(f'PASS {len(cases)} review regression cases; EVIDENCE {root}')
