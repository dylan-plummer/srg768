"""Section 7 of the paper: the original seed (BASE_SEED, slopes {0,1,t,t^3+t}) equals z_{0,0,1} o M0^{-1}
with M0 = ((t^2+1, t^2), (0, 1)), i.e. ((5, 4), (0, 1)) in the integer encoding."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from srg768 import BASE_SEED, norm_seed, transport, support_slopes
z = transport(norm_seed(0, 0, 1), ((5, 4), (0, 1)))
assert (z == BASE_SEED).all()
print('original seed = z_{0,0,1} o M0^{-1}, M0 = ((t^2+1, t^2), (0, 1)); support slopes', sorted(support_slopes(z)))
