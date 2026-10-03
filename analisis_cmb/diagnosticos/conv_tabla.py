import sys; from pathlib import Path; sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np, ou_cmb as M, planck_datos as P
from scipy.signal import lfilter
from scipy.ndimage import gaussian_filter1d
om=M.omega_l0(); tt=P.load_tt(); err=np.minimum(tt["err_lo"],tt["err_hi"]); lens=P.load_lensing()
NT=2000; L0=np.log(M.A_MIN); htab=-L0/(NT-1)
KMAX=64; h=htab/KMAX
npad_tab=int(np.ceil(M.PAD/htab)); nfine=(NT-1+npad_tab)*KMAX+1
lna_f=L0+h*np.arange(nfine)
def ou_uniform(n,hh,theta,s,z):
    r=np.exp(-theta*hh); x0=s*z[0]
    y,_=lfilter([s*np.sqrt(1-r*r)],[1,-r],z[1:],zi=[r*x0]); return np.concatenate([[x0],y])
def smooth(x,hh,d):
    sg=d/hh; return gaussian_filter1d(x,sg,mode="constant",truncate=6)/gaussian_filter1d(np.ones_like(x),sg,mode="constant",truncate=6)
res={}
for seed in range(10):
  z=np.random.default_rng(M.SEED+seed).standard_normal(nfine)
  x=ou_uniform(nfine,h,1.0,M.SIGMA_X,z)
  for d in (0.05,0.1):
    for k in (64,16):
      step=1; xf=x; lf=lna_f
      xs=smooth(xf,h*step,d)
      tab=np.arange(0,len(lf),k)  # table points (spacing htab), includes pad
      lt=lf[tab]; xt=xs[tab]; mask=lt<=1e-12
      rho,w=M.rho_and_w(lt,mask,xt[None,:],om)
      r=M.run(M.make_params(a=np.exp(lt[mask]),w=w[0]),keep_results=False)
      res[(seed,d,k)]=(r["DlTT"][tt["ell"]],P.bin_pp(r["Clpp"],lens),r["thetastar100"])
for d in (0.05,0.1):
  for k1,k2 in ((64,16),):
    m=np.array([[np.max(np.abs(res[(s,d,k1)][0]-res[(s,d,k2)][0])/err),np.max(np.abs(res[(s,d,k1)][1]-res[(s,d,k2)][1])/lens["err"]),abs(res[(s,d,k1)][2]-res[(s,d,k2)][2])/P.THETASTAR100_ERR] for s in range(10)])
    print(f"D={d} tabla N={(NT-1)*KMAX//k1+1} vs N={(NT-1)*KMAX//k2+1} (camino h_tab/64) (10 semillas, max): TT {m[:,0].max():.2e} pp {m[:,1].max():.2e} th* {m[:,2].max():.2e} (mediana th* {np.median(m[:,2]):.1e})")
