"""
Maxwell + device -> M : derivation verification (numerology-free).

Derives the optical transfer operator for a real device (phase screen SLM +
free propagation + optional thin lens / spectral shaping) starting from the
paraxial wave equation, and checks the three quantitative claims that decide
whether a "tesseract / B4 / 8-cube" channel can appear:

  (1) Fresnel propagator P_d  and  quadratic-phase lens L_f  compose to the
      standard ABCD (metaplectic / Collins integral). Verified on a test field.

  (2) Controllable degrees of freedom of the device = number of independently
      addressable phase pixels (x spectral bins if pulse-shaped). Hardware,
      not group dimension.

  (3) The only legitimate entry of a 4-cube symmetry is as a discrete subgroup
      of Sp(4,R) (the 4D optical phase space x,theta_x,y,theta_y). Computing
      |B4 ∩ Sp(4,R)| and the dimension of its commutant inside sp(4,R) shows
      the symmetry CONSTRAINS (reduces) the operator, it does not amplify it.
"""
import numpy as np
from itertools import permutations, product

# ============================================================
# (1) Fresnel + quadratic phase == ABCD  (Collins integral)
# ============================================================
lam = 633e-9
k = 2*np.pi/lam

# 1D grid (generalises to 2D; the algebra is identical)
N = 2048
L = 6e-3                      # 6 mm window
x = np.linspace(-L/2, L/2, N)
dx = x[1]-x[0]

def gaussian_in():
    w0 = 0.30e-3              # 0.30 mm waist at z=0
    return np.exp(-x**2/w0**2)            # waist at z=0, flat phase

def fresnel_propagate(u, d):
    """Free-space propagation over distance d via the angular-spectrum (FFT) method.
    This is the exact discretization of the same Fresnel operator whose closed form
    is the Collins kernel with ABCD=(1 d; 0 1)."""
    kx = 2*np.pi*np.fft.fftfreq(N, dx)
    H = np.exp(-1j*kx**2*d/(2*k))          # paraxial transfer function exp(i k_z d)
    return np.fft.ifft(np.fft.fft(u)*H)

# analytic check: a Gaussian waist w0 at z=0 stays Gaussian, w(d)=w0*sqrt(1+(d/z_R)^2)
def thin_lens(u, f):
    return u * np.exp(-1j*k*x**2/(2*f))

# Device:  free d1 -> lens f -> free d2.  ABCD = (1 d2;0 1)(1 0;-1/f 1)(1 d1;0 1)
d1, f, d2 = 0.12, 0.20, 0.30
A = 1 - d2/f;  B = d1 + d2 - d1*d2/f;  C = -1/f;  D = 1 - d1/f
# sanity: det = 1
assert abs(A*D - B*C - 1) < 1e-12

u0 = gaussian_in()
# apply device via the *integral operators* (Fresnel + phase)
u1 = fresnel_propagate(u0, d1)
u2 = thin_lens(u1, f)
u_dev = fresnel_propagate(u2, d2)

# analytic free-propagation sanity (waist spread)
z_R = np.pi*(0.30e-3)**2/lam
w_pred = 0.30e-3*np.sqrt(1+(d1/z_R)**2)
u1p = fresnel_propagate(u0, d1)
sig = np.sqrt(np.sum((np.abs(u1p)**2)*(x**2))/np.sum(np.abs(u1p)**2))
w_num = 2*sig
print(f"    [sanity] Gaussian waist after d={d1} m: analytic w={w_pred*1e6:.2f} um, numeric={w_num*1e6:.2f} um")

# Apply the SAME device via the analytic ABCD law for the Gaussian q-parameter:
#   q_out = (A q_in + B)/(C q_in + D),  q_in = i z_R  (waist, flat phase at input)
# Then extract q from the numerically-propagated field and compare.
def extract_q(u):
    I = np.abs(u)**2; I /= I.sum()
    sig = np.sqrt(np.sum(I*x**2)); w = 2*sig
    phi = np.unwrap(np.angle(u))
    m = I > 0.3*I.max()                      # fit curvature on the bright core
    c = np.polyfit(x[m], phi[m], 2)          # phi ~ c2 x^2 + c1 x + c0
    invR = 2*c[0]/k
    inv_q = invR - 1j*lam/(np.pi*w**2)
    return 1.0/inv_q, w, 1.0/invR

q_in = 1j*z_R
q_abcd = (A*q_in + B)/(C*q_in + D)
q_num, w_num, invR_num = extract_q(u_dev)

print("(1) Device operator == ABCD  (Gaussian q-parameter test):")
print(f"    q_in  = {q_in:.5f}")
print(f"    q_ABCD= (Aq+B)/(Cq+D) = {q_abcd.real:.5f} + {q_abcd.imag:.5f}i")
print(f"    q_num (extracted from field) = {q_num.real:.5f} + {q_num.imag:.5f}i")
print(f"    -> relative |q_num - q_ABCD|/|q_ABCD| = {abs(q_num-q_abcd)/abs(q_abcd):.2e}")
print(f"    (device = Fresnel(+lens)+Fresnel reproduces the ABCD action on the mode)")
print()

# ============================================================
# (2) Controllable DOF = N_pix  (x N_omega if spectral shaping)
# ============================================================
SLM_W, SLM_H = 1920, 1080
N_pix = SLM_W*SLM_H
N_omega = 128                       # typical pulse-shaper bins
print("(2) Controllable degrees of freedom of the device:")
print(f"    phase-only SLM {SLM_W}x{SLM_H}: N_pix = {N_pix} real phase params")
print(f"    + spectral shaper ({N_omega} bins): {N_pix*N_omega} params")
print(f"    -> set by HARDWARE pixel/bin count, not by any polytope combinatorics.")
print()

# ============================================================
# (3) 4-cube symmetry B4 as a subgroup of Sp(4,R); commutant in sp(4,R)
# ============================================================
# Phase-space ordering: (x, theta_x, y, theta_y).  J = blockdiag(J2,J2), J2=[[0,1],[-1,0]]
J2 = np.array([[0.,1.],[-1.,0.]])
J  = np.block([[J2, np.zeros((2,2))],[np.zeros((2,2)), J2]])
n = 4

def is_sp(M, tol=1e-9):
    return np.allclose(M.T@J@M, J, atol=tol)

# B4 = signed permutations of the 4 coordinates
def signed_perms():
    mats=[]
    for perm in permutations(range(4)):
        P = np.eye(4)[:, list(perm)]
        for signs in product([1,-1], repeat=4):
            mats.append(np.diag(signs).astype(float)@P)
    return mats

B4 = signed_perms()
print(f"(3) |B4| = {len(B4)} (signed permutations of 4 coords = 2^4 * 4! = 384)")

G = [M for M in B4 if is_sp(M)]
print(f"    |B4 ∩ Sp(4,R)| = {len(G)}  (only these preserve the optical symplectic form)")

# Basis of sp(4,R)  <=>  symmetric 4x4 matrices via X = -J S, S symmetric (since J^2=-I)
sym_basis=[]
for i in range(n):
    for j in range(i,n):
        S=np.zeros((n,n)); S[i,j]=S[j,i]=1.0
        sym_basis.append(S)
sp_basis=[-J@S for S in sym_basis]          # 10 matrices
assert len(sp_basis)==10
def in_sp(X, tol=1e-9):                      # Lie-algebra condition: X^T J + J X = 0
    return np.allclose(X.T@J + J@X, 0.0, atol=tol)
assert all(in_sp(X) for X in sp_basis)
# orthonormalise (Frobenius) for a clean projection
def flat(X): return X.reshape(-1)
Bm = np.stack([flat(X)/np.linalg.norm(flat(X)) for X in sp_basis], axis=1)  # 16x10

# Ad_g : X -> g X g^{-1}  = g X g^T  (g orthogonal). Projector onto fixed subspace:
#   P_fix = (1/|G|) sum_g Ad_g   restricted to sp_basis
P = np.zeros((10,10))
for g in G:
    for a,Xa in enumerate(sp_basis):
        ad = g@Xa@g.T
        P[:,a] += Bm.T@flat(ad)
P /= len(G)
commutant_dim = int(np.linalg.matrix_rank(P, tol=1e-9))
print(f"    dim sp(4,R) = {len(sp_basis)} ;  dim( commutant of G in sp(4,R) ) = {commutant_dim}")
print(f"    -> B4-style symmetry leaves a {commutant_dim}-dim subspace: it CONSTRAINS")
print(f"       the operator (10 -> {commutant_dim}), it does NOT create extra channels.")
print()
print("Comb (8 cubes / 16 vertices / 24 faces / 32 edges of the 4-cube) never appears")
print("as a rank or DOF in any of the three checks. To obtain a 'tesseract M' one must")
print("INSERT it as an ansatz; it is not produced by Maxwell + device geometry.")
