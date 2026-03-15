import numpy as np
import matplotlib.pyplot as plt
from scipy.special import zeta
from scipy.integrate import quad
from pathlib import Path

# =============================
# Parameters
# =============================
hbar = 1.0
m = 1.0
c = 1.0
k_B = 1.0

xi = 0.1
phi = 0.5
Phi_0 = 0.3
g_d = 2

eps_enh = -(5 * xi * phi) / 2.0   # = -0.125  -> enhancement
eps_sup = +(5 * xi * phi) / 2.0   # = +0.125  -> suppression

# =============================
# Coefficients
# =============================
C_7 = (189 * g_d) / (200 * np.pi**2) * zeta(7) / (hbar**3 * m**3 * c**9)
C_4 = (21 * g_d) / (40 * np.pi**2) * zeta(4) * (xi * Phi_0) / (hbar**2 * m**2 * c**6)
K_UR = (7 * np.pi**2 / 240) * g_d / (hbar**3 * c**3)
a = np.sqrt(5 * hbar * m * c**3)

T_c_phi = (C_4 / C_7) ** (1 / 3)
T_cross_UR = (K_UR / C_7) ** (1 / 3)

# =============================
# Pressure functions
# =============================
def P_free(T):
    return C_7 * (k_B * T) ** 7


def P_const_proj(T):
    return C_7 * (k_B * T) ** 7 + C_4 * (k_B * T) ** 4


def P_UR(T):
    return K_UR * (k_B * T) ** 4


def P_parallel_direct(T, eps_s):
    """
    Robust evaluation from the grand-potential integral.

    This avoids relying on a naive termwise Fermi-series implementation for
    eps_s < 0, where the formal series in exp(-n eps_s/T) is not the safest
    numerical representation. The integral is finite and directly matches the
    thermodynamic definition.
    """
    pref = g_d / (2 * np.pi**2)

    def integrand(k):
        return k**2 * np.log1p(np.exp(-(a * np.sqrt(k) + eps_s) / T))

    total = 0.0
    for lo, hi in [(0.0, 4.0), (4.0, 25.0), (25.0, 100.0), (100.0, np.inf)]:
        val, _ = quad(integrand, lo, hi, limit=300, epsabs=1e-10, epsrel=1e-10)
        total += val
    return T * pref * total


# =============================
# Temperature grids
# =============================
T_main = np.logspace(-0.5, 1.0, 240)   # [0.316..., 10]
T_wide = np.logspace(np.log10(0.08), 1.0, 260)

# =============================
# Compute curves
# =============================
P_free_main = P_free(T_main)
P_const_main = P_const_proj(T_main)
P_UR_main = P_UR(T_main)
P_par_enh_main = np.array([P_parallel_direct(t, eps_enh) for t in T_main])
P_par_sup_main = np.array([P_parallel_direct(t, eps_sup) for t in T_main])

ratio_enh = P_par_enh_main / P_free_main
ratio_sup = P_par_sup_main / P_free_main
ratio_const = P_const_main / P_free_main

P_free_wide = P_free(T_wide)
P_const_wide = P_const_proj(T_wide)
P_UR_wide = P_UR(T_wide)
term_T7_wide = C_7 * T_wide**7
term_T4_wide = C_4 * T_wide**4

outdir = Path(".")

# =============================
# Figure 1: Main text
# =============================
fig1 = plt.figure(figsize=(14, 6.5))
gs = fig1.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.24)
ax1 = fig1.add_subplot(gs[0, 0])
ax2 = fig1.add_subplot(gs[0, 1])

ax1.loglog(T_main, P_free_main, color="black", lw=3.0, label="PCP free")
ax1.loglog(T_main, P_par_enh_main, color="tab:blue", lw=2.4, ls="--",
           label=r"Parallel branch, $\epsilon_s=-0.125$")
ax1.loglog(T_main, P_par_sup_main, color="tab:blue", lw=2.4, ls=":",
           label=r"Parallel branch, $\epsilon_s=+0.125$")
ax1.loglog(T_main, P_const_main, color="tab:red", lw=2.6, ls="-",
           label=r"Constant-projection branch")
ax1.loglog(T_main, P_UR_main, color="tab:green", lw=2.4, ls="-.",
           label="Standard UR fermions")

T_ref = np.array([0.45, 1.4])
y7 = 0.04 * T_ref**7
y4 = 0.18 * T_ref**4
ax1.loglog(T_ref, y7, color="0.35", lw=1.2, ls="--")
ax1.loglog(T_ref, y4, color="0.55", lw=1.2, ls="--")
ax1.text(0.88, 0.04 * 0.88**7 * 1.5, r"$T^7$", fontsize=11, color="0.25")
ax1.text(0.88, 0.18 * 0.88**4 * 1.1, r"$T^4$", fontsize=11, color="0.4")

ax1.set_xlabel(r"$T$")
ax1.set_ylabel(r"$P(T)$")
ax1.set_title("(a) Absolute pressures")
ax1.grid(True, which="both", alpha=0.25)
ax1.set_xlim(0.316, 10.0)
y_all = np.concatenate([P_free_main, P_par_enh_main, P_par_sup_main, P_const_main, P_UR_main])
ax1.set_ylim(y_all.min() * 0.7, y_all.max() * 1.35)
ax1.legend(fontsize=9, frameon=True, loc="upper left")

ax2.semilogx(T_main, ratio_enh, color="tab:blue", lw=2.4, ls="--",
             label=r"$P_{\parallel}(\epsilon_s=-0.125)/P_{\mathrm{free}}$")
ax2.semilogx(T_main, ratio_sup, color="tab:blue", lw=2.4, ls=":",
             label=r"$P_{\parallel}(\epsilon_s=+0.125)/P_{\mathrm{free}}$")
ax2.semilogx(T_main, ratio_const, color="tab:red", lw=2.6, ls="-",
             label=r"$P_{\Phi_0}/P_{\mathrm{free}}$")
ax2.axhline(1.0, color="black", lw=1.0)
ax2.axhspan(0.95, 1.05, color="gold", alpha=0.18)
ax2.annotate("Low-T deviations",
             xy=(0.36, 1.50), xytext=(0.46, 1.58),
             arrowprops=dict(arrowstyle="->", lw=0.9),
             fontsize=10)
ax2.annotate(r"High-T convergence" "\n" r"for $T \gtrsim 2.5$--$5$",
             xy=(4.5, 1.01), xytext=(2.2, 1.18),
             arrowprops=dict(arrowstyle="->", lw=0.9),
             fontsize=10)
ax2.set_xlabel(r"$T$")
ax2.set_ylabel(r"$P/P_{\mathrm{free}}$")
ax2.set_title("(b) Ratios relative to PCP free")
ax2.set_xlim(0.316, 10.0)
ax2.set_ylim(0.6, 1.65)
ax2.grid(True, which="both", alpha=0.25)
ax2.legend(fontsize=9, frameon=True, loc="upper right")

fig1.tight_layout()
fig1.savefig(outdir / "pcp_pressure_comparison_main.png", dpi=300, bbox_inches="tight")
fig1.savefig(outdir / "pcp_pressure_comparison_main.pdf", bbox_inches="tight")
plt.close(fig1)

# =============================
# Figure 2: Appendix crossover
# =============================
fig2, ax = plt.subplots(figsize=(8.5, 6))
ax.loglog(T_main, P_free_main, color="black", lw=3.0, label=r"$P_{\mathrm{free}}=\mathcal{{C}}_7 T^7$")
ax.loglog(T_main, P_UR_main, color="tab:green", lw=2.6, ls="-.", label=r"$P_{\mathrm{{UR}}}=K_{\mathrm{{UR}}} T^4$")
ax.loglog(T_main, P_const_main, color="tab:red", lw=2.6, label=r"$P_{{\Phi_0}}=\mathcal{{C}}_7 T^7+\mathcal{{C}}_4 T^4$")
P_cross = P_free(T_cross_UR)
ax.plot(T_cross_UR, P_cross, marker="o", ms=8, color="black")
ax.annotate(rf"$T_\times \approx {T_cross_UR:.3f}$",
            xy=(T_cross_UR, P_cross),
            xytext=(2.0, P_cross * 2.2),
            arrowprops=dict(arrowstyle="->", lw=0.9),
            fontsize=10)
ax.set_xlabel(r"$T$")
ax.set_ylabel(r"$P(T)$")
ax.set_title("PCP versus standard ultra-relativistic scaling")
ax.grid(True, which="both", alpha=0.25)
ax.legend(frameon=True, fontsize=10, loc="upper left")
ax.set_xlim(0.316, 10.0)
ax.set_ylim(min(P_UR_main.min(), P_free_main.min(), P_const_main.min()) * 0.8,
            max(P_UR_main.max(), P_free_main.max(), P_const_main.max()) * 1.3)
fig2.tight_layout()
fig2.savefig(outdir / "pcp_pressure_crossover_appendix.png", dpi=300, bbox_inches="tight")
fig2.savefig(outdir / "pcp_pressure_crossover_appendix.pdf", bbox_inches="tight")
plt.close(fig2)

# =============================
# Figure 3: Appendix decomposition
# =============================
fig3, ax = plt.subplots(figsize=(8.5, 6))
ax.loglog(T_wide, P_const_wide, color="tab:red", lw=2.8, label=r"$P_{\Phi_0}$")
ax.loglog(T_wide, term_T7_wide, color="orange", lw=2.5, ls="--", label=r"$\mathcal{{C}}_7 T^7$")
ax.loglog(T_wide, term_T4_wide, color="purple", lw=2.5, ls=":", label=r"$\mathcal{{C}}_4 T^4$")
P_tc = C_7 * T_c_phi**7 + C_4 * T_c_phi**4
ax.plot(T_c_phi, P_tc, marker="o", ms=8, color="tab:red")
ax.annotate(rf"Formal equality point" "\n" rf"$T_c^{{(\Phi)}} \approx {T_c_phi:.4f}$",
            xy=(T_c_phi, P_tc),
            xytext=(0.42, P_tc * 2.0),
            arrowprops=dict(arrowstyle="->", lw=0.9),
            fontsize=10)
ax.set_xlabel(r"$T$")
ax.set_ylabel(r"$P(T)$")
ax.set_title("Internal decomposition of the constant-projection branch")
ax.grid(True, which="both", alpha=0.25)
ax.legend(frameon=True, fontsize=10, loc="upper left")
ax.set_xlim(0.08, 10.0)
ax.set_ylim(min(term_T4_wide.min(), term_T7_wide.min(), P_const_wide.min()) * 0.8,
            max(term_T4_wide.max(), term_T7_wide.max(), P_const_wide.max()) * 1.3)
fig3.tight_layout()
fig3.savefig(outdir / "pcp_pressure_decomposition_appendix.png", dpi=300, bbox_inches="tight")
fig3.savefig(outdir / "pcp_pressure_decomposition_appendix.pdf", bbox_inches="tight")
plt.close(fig3)

# =============================
# Console cross-checks
# =============================
for T_test in [10**(-0.5), 5.0]:
    r_enh = P_parallel_direct(T_test, eps_enh) / P_free(T_test)
    r_sup = P_parallel_direct(T_test, eps_sup) / P_free(T_test)
    r_const = P_const_proj(T_test) / P_free(T_test)
    print(f"T = {{T_test:.6f}}")
    print(f"  P_parallel(eps_enh)/P_free = {{r_enh:.9f}}")
    print(f"  P_parallel(eps_sup)/P_free = {{r_sup:.9f}}")
    print(f"  P_const_proj/P_free      = {{r_const:.9f}}")

print(f"C_7 = {{C_7:.9f}}")
print(f"C_4 = {{C_4:.9f}}")
print(f"C_7/C_4 = {{C_7/C_4:.9f}}")
print(f"eps_enh = {{eps_enh:.3f}}")
print(f"eps_sup = {{eps_sup:.3f}}")
print(f"T_c_phi = {{T_c_phi:.9f}}")
print(f"T_cross_UR = {{T_cross_UR:.9f}}")
