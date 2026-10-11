"""Check probe output structure, and detect the known uniform mixing failure.

Usage: python3 check_samples.py output.log [--uniform]
General success checks execution/bitstrings, not distributional accuracy.
The --uniform check fails on the Metropolis all-zero reproducer.
Use --accuracy for a conservative per-site error bound (six binomial standard
errors at p=1/2); this is meaningful for the independent heat-bath product case.
"""
import math
import re
import sys
from pathlib import Path

text=Path(sys.argv[1]).read_text()
requested=re.search(r'^mode=.* samples=(\d+)$',text,re.M)
assert requested, 'missing run header'
samples=re.findall(r'^sample \d+: ([01]+)$',text,re.M)
assert len(samples)==int(requested[1]), 'missing retained samples'
assert all(len(s)==51 for s in samples), 'sample length is not 51'
assert re.search(r'SAMPLING (?:COMPLETED|PASS)',text), 'sampling did not complete'
marginals=re.findall(r'^marginal \d+ expected=(\S+) observed=(\S+)$',text,re.M)
assert len(marginals)==51 and all(math.isfinite(float(x)) and 0<=float(x)<=1 for row in marginals for x in row)
if '--uniform' in sys.argv and len(set(samples))==1:
    print('FAIL: uniform product state produced one identical bitstring for every retained sample')
    sys.exit(1)
if '--accuracy' in sys.argv:
    bound=3/math.sqrt(len(samples))
    error=max(abs(float(expected)-float(observed)) for expected,observed in marginals)
    assert error<bound, f'max marginal error {error} exceeds {bound}'
    print(f'ACCURACY: max marginal error {error:.6g} < {bound:.6g}')
print(f'COMPLETE: {len(samples)} valid 51-bit samples; {len(set(samples))} distinct states (no general accuracy assertion)')
