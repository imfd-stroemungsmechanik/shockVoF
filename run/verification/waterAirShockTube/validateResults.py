import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size": 6,
    "lines.linewidth": 0.8,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4,
    "grid.linewidth": 0.4,
})

# Stiffened-gas EoS
class StiffenedGasEoS:
    def __init__(self, gamma, pInf, Cv):
        self.gamma, self.pInf, self.Cv = gamma, pInf, Cv
        self.gp1, self.gm1 = gamma + 1.0, gamma - 1.0
        self.g_ratio = self.gm1 / self.gp1
        self.R = Cv * self.gm1

    def c(self, rho, p):
        return np.sqrt(self.gamma * (p + self.pInf) / rho)

    def T(self, rho, p):
        return (p + self.pInf) / (rho * self.R)


# Exact multi-material Riemann solver
class ExactRiemann:
    def __init__(self, eosL, eosR, rhoL, pL, uL, rhoR, pR, uR):
        self.eosL, self.eosR = eosL, eosR
        self.rhoL, self.pL, self.uL = rhoL, pL, uL
        self.rhoR, self.pR, self.uR = rhoR, pR, uR
        self.cL = eosL.c(rhoL, pL)
        self.cR = eosR.c(rhoR, pR)

        # Acoustic initial guess
        zL, zR = rhoL * self.cL, rhoR * self.cR
        ps = max((zR*pL + zL*pR + zL*zR*(uL - uR)) / (zL + zR), 1e-6)

        # Newton iteration for p*
        for _ in range(200):
            fL, dfL = self._f(eosL, ps, rhoL, pL, self.cL)
            fR, dfR = self._f(eosR, ps, rhoR, pR, self.cR)
            Phi = fL + fR + (uR - uL)
            dPhi = dfL + dfR
            if abs(dPhi) < 1e-30: break
            dp = -Phi / dPhi
            p_min = max(-eosL.pInf, -eosR.pInf)
            p_new = ps + dp
            while p_new <= p_min + 1e-6:
                dp *= 0.5
                p_new = ps + dp
            if abs(dp) / (0.5 * abs(p_new + ps) + 1e-30) < 1e-10:
                ps = p_new; break
            ps = p_new

        self.pS = ps
        fL, _ = self._f(eosL, ps, rhoL, pL, self.cL)
        fR, _ = self._f(eosR, ps, rhoR, pR, self.cR)
        self.uS = 0.5 * (uL + uR) + 0.5 * (fR - fL)
        self.rhoSL = self._rho_star(eosL, rhoL, pL, ps)
        self.rhoSR = self._rho_star(eosR, rhoR, pR, ps)

    @staticmethod
    def _f(eos, ps, rho, p, c):
        pt = p + eos.pInf
        if ps > p:  # shock
            A = 2.0 / (eos.gp1 * rho)
            B = eos.g_ratio * pt + eos.pInf
            s = np.sqrt(A / (ps + B))
            return (ps - p) * s, s * (1.0 - 0.5 * (ps - p) / (ps + B))
        else:  # rarefaction
            exp = eos.gm1 / (2.0 * eos.gamma)
            r = ((ps + eos.pInf) / pt)**exp
            return (2.0 * c / eos.gm1 * (r - 1.0),
                    c / (eos.gamma * pt) * ((ps + eos.pInf) / pt)**(exp - 1.0))

    @staticmethod
    def _rho_star(eos, rho, p, ps):
        r = (ps + eos.pInf) / (p + eos.pInf)
        if ps > p:
            g = eos.g_ratio
            return rho * (r + g) / (g * r + 1.0)
        return rho * r**(1.0 / eos.gamma)

    def sample(self, x, t, x0=0.5):
        S = (x - x0) / t
        rho, u, p, alpha = [np.empty_like(x) for _ in range(4)]
        for i, s in enumerate(S):
            if s < self.uS:
                alpha[i] = 1.0
                rho[i], u[i], p[i] = self._left(s)
            else:
                alpha[i] = 0.0
                rho[i], u[i], p[i] = self._right(s)
        return rho, u, p, alpha

    def _left(self, s):
        eos = self.eosL
        if self.pS > self.pL:  # shock
            A = 2.0 / (eos.gp1 * self.rhoL)
            B = eos.g_ratio * (self.pL + eos.pInf) + eos.pInf
            shock_s = self.uL - np.sqrt((self.pS + B) / A) / self.rhoL
            if s < shock_s: return self.rhoL, self.uL, self.pL
            return self.rhoSL, self.uS, self.pS
        else:  # rarefaction
            cS = eos.c(self.rhoSL, self.pS)
            if s < self.uL - self.cL: return self.rhoL, self.uL, self.pL
            if s > self.uS - cS: return self.rhoSL, self.uS, self.pS
            pt = self.pL + eos.pInf
            cc = (2.0 / eos.gp1) * (self.cL + 0.5 * eos.gm1 * (self.uL - s))
            uu = s + cc
            pp = pt * (cc / self.cL)**(2.0 * eos.gamma / eos.gm1) - eos.pInf
            rr = self.rhoL * ((pp + eos.pInf) / pt)**(1.0 / eos.gamma)
            return rr, uu, pp

    def _right(self, s):
        eos = self.eosR
        if self.pS > self.pR:  # shock
            A = 2.0 / (eos.gp1 * self.rhoR)
            B = eos.g_ratio * (self.pR + eos.pInf) + eos.pInf
            shock_s = self.uR + np.sqrt((self.pS + B) / A) / self.rhoR
            if s > shock_s: return self.rhoR, self.uR, self.pR
            return self.rhoSR, self.uS, self.pS
        else:  # rarefaction
            cS = eos.c(self.rhoSR, self.pS)
            if s > self.uR + self.cR: return self.rhoR, self.uR, self.pR
            if s < self.uS + cS: return self.rhoSR, self.uS, self.pS
            pt = self.pR + eos.pInf
            cc = (2.0 / eos.gp1) * (self.cR + 0.5 * eos.gm1 * (s - self.uR))
            uu = s - cc
            pp = pt * (cc / self.cR)**(2.0 * eos.gamma / eos.gm1) - eos.pInf
            rr = self.rhoR * ((pp + eos.pInf) / pt)**(1.0 / eos.gamma)
            return rr, uu, pp


# Read OpenFOAM data
t_end = 0.00015
x, p_of, rho_of, U_of, alpha_of = np.loadtxt("postProcessing/sample/" + str(t_end) + "/data.xy", unpack=True)


# Exact solution
eosL = StiffenedGasEoS(gamma=4.4, pInf=6.0e8, Cv=602.1)
eosR = StiffenedGasEoS(gamma=1.4, pInf=0.0,   Cv=717.5)

pL, TL, uL = 1.0e9, 781.58, 0.0
pR, TR, uR = 1.0e5, 6.9686, 0.0
rhoL = (pL + eosL.pInf) / (eosL.R * TL)
rhoR = (pR + eosR.pInf) / (eosR.R * TR)

solver = ExactRiemann(eosL, eosR, rhoL, pL, uL, rhoR, pR, uR)

x_ex = np.linspace(0, 1, 2001)
rho_ex, u_ex, p_ex, alpha_ex = solver.sample(x_ex, t_end)


# Plot

fig, axes = plt.subplots(2, 2, figsize=(6.3, 4.25), layout="constrained")

axes[0][0].plot(x_ex, p_ex/1e9, color="black", linewidth=0.75, linestyle=(0, (5, 5)))
axes[0][0].plot(x, p_of/1e9, color="#4292C6")
axes[0][0].set_ylabel(r"Pressure [GPa]")

# Zoomed inset: rarefaction tail and star-region pressure plateau
axins = axes[0][0].inset_axes([0.46, 0.38, 0.50, 0.55])
axins.plot(x_ex, p_ex/1e9, color="black", linewidth=0.75, linestyle=(0, (5, 5)))
axins.plot(x, p_of/1e9, color="#4292C6")
axins.set_xlim(0.25, 0.65)
axins.set_ylim(-0.003, 0.033)
axins.set_yticks([0.0, 0.01, 0.02, 0.03])
axins.tick_params(labelsize=5.5, pad=1.5)
axins.grid(alpha=0.5)
axes[0][0].indicate_inset_zoom(axins)

axes[0][1].plot(x_ex, u_ex, color="black", linewidth=0.75, linestyle=(0, (5, 5)))
axes[0][1].plot(x, U_of, color="#4292C6")
axes[0][1].set_ylabel(r"Velocity [m/s]")

axes[1][0].plot(x_ex, rho_ex, color="black", linewidth=0.75, linestyle=(0, (5, 5)))
axes[1][0].plot(x, rho_of, color="#4292C6")
axes[1][0].set_ylabel(r"Density [kg/m$^3$]")

axes[1][1].plot(x_ex, alpha_ex, label="Exact", color="black", linewidth=0.75, linestyle=(0, (5, 5)))
axes[1][1].plot(x, alpha_of, label="Simulation", color="#4292C6")
axes[1][1].set_ylabel(r"Phase fraction [-]")
axes[1][1].legend(loc="lower left")

for ax in axes.flat:

    ax.set_xlabel(r"$x$-coordinate [m]")
    ax.grid(alpha=0.5)


plt.savefig("results_waterAirShockTube.png", dpi=450, bbox_inches="tight")
