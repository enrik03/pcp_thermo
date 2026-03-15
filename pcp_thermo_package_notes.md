# pcp-thermodynamics: Technical Notes
## Post-Carroll-Pauli (PCP) Thermodynamic Pressure Analysis

---

### 1. Overview

This document collects the technical derivations, numerical values, figure captions,
and LaTeX fragments for the PCP thermodynamic analysis implemented in
`pcp_thermo_figures_corrected.py`.

---

### 2. Notation and Units

| Symbol | Meaning | Units |
|--------|---------|-------|
| $M$ | PCP heavy-field mass scale | GeV (set to 1 in dimensionless runs) |
| $T$ | Temperature | GeV |
| $N_c$ | Number of colours | 3 |
| $N_f$ | Number of light flavours | 2 |
| $g_{\rm eff}$ | Effective relativistic degrees of freedom | dimensionless |
| $C_4$ | Coefficient of $T^4$ in PCP pressure | dimensionless |
| $C_7$ | Coefficient of $T^7/M^3$ in PCP pressure | $M^{-3}$ absorbed |
| $T_c^{(\Phi)}$ | PCP condensate / characteristic temperature | GeV |
| $T_x$ | PCP–UR crossover temperature | GeV |

We work in natural units $\hbar = c = k_B = 1$.

---

### 3. Ultra-Relativistic Pressure

The ideal-gas Stefan–Boltzmann pressure is

$$
P_{\rm UR}(T) = C_{\rm UR}\,T^4, \qquad
C_{\rm UR} = \frac{\pi^2}{90}\,g_{\rm eff},
$$

with the effective degrees of freedom

$$
g_{\rm eff} = 2(N_c^2-1) + \tfrac{7}{8}\cdot 2\cdot 2\,N_c N_f.
$$

For $N_c = 3$, $N_f = 2$:

$$
g_{\rm eff} = 16 + \tfrac{7}{8}\cdot 24 = 16 + 21 = 37, \qquad
C_{\rm UR} = \frac{37\pi^2}{90}.
$$

---

### 4. PCP Pressure Expansion

The PCP sector modifies the dispersion relation of the quark and gluon quasiparticles.
At one loop the thermal pressure can be written as a power series in $T/M$:

$$
\boxed{
P_{\rm PCP}(T) = C_4\,T^4 + C_7\,\frac{T^7}{M^3} + \mathcal{O}\!\left(\frac{T^{10}}{M^6}\right).
}
$$

#### 4.1 Coefficient $C_4$

$C_4$ is determined by the one-loop bosonic and fermionic thermal integrals in the
massless limit,

$$
C_4 = -\frac{g_B}{2\pi^2}\int_0^\infty du\,u^2\ln(1-e^{-u})
      + \frac{g_F}{2\pi^2}\int_0^\infty du\,u^2\ln(1+e^{-u}),
$$

where $g_B = 2(N_c^2-1)$ (gluons) and $g_F = 2 N_c N_f$ (quarks, spin and
colour already absorbed).  Numerically:

$$
C_4 \approx 0.2560 \quad \text{(PCP leading term)}.
$$

#### 4.2 Coefficient $C_7$

$C_7$ arises from the leading PCP correction to the heavy-field propagator
at finite temperature:

$$
C_7 = \frac{g_{\rm PCP}}{2\pi^2 M^3}\int_0^\infty du\,u^5 \left(-\frac{e^u}{(e^u-1)^2}\right),
\qquad g_{\rm PCP} = \frac{1}{16\pi^2}.
$$

Numerically,

$$
C_7 \approx 1.25 \times 10^{-5}\,M^{-3}.
$$

---

### 5. Characteristic Temperatures

#### 5.1 PCP Condensate Temperature $T_c^{(\Phi)}$

The Phi field develops a thermal mass

$$
m_{\rm th}^2(T) = M^2 - \frac{g\,T^2}{12},
$$

where $g \approx 1/4$ is a representative PCP quartic coupling.  The condensate
temperature is defined as $m_{\rm th}(T_c^{(\Phi)}) = 0$:

$$
T_c^{(\Phi)} = \sqrt{\frac{12}{g}}\,M = \sqrt{48}\,M \approx 6.93\,M.
$$

#### 5.2 Crossover Temperature $T_x$

The crossover between the PCP and UR descriptions occurs when

$$
P_{\rm PCP}(T_x) = P_{\rm UR}(T_x)
\implies
C_7\,\frac{T_x^3}{M^3} = C_{\rm UR} - C_4,
$$

which gives

$$
T_x = M\left(\frac{C_{\rm UR}-C_4}{C_7}\right)^{1/3}.
$$

---

### 6. Figure Captions (LaTeX)

#### Figure 1 – Main text

```latex
\begin{figure}[htb]
  \centering
  \includegraphics[width=0.8\textwidth]{pcp_pressure_comparison_main}
  \caption{%
    Normalised pressure $P(T)/T^4$ as a function of temperature $T/M$
    for the PCP sector (solid blue) and the ultra-relativistic ideal gas
    (dashed red).  The vertical dotted line marks the crossover temperature
    $T_x$ at which $P_{\rm PCP}(T_x) = P_{\rm UR}(T_x)$.
    Dimensionless units $M=1$ are used throughout.
  }
  \label{fig:pcp_main}
\end{figure}
```

#### Figure 2 – Appendix (crossover)

```latex
\begin{figure}[htb]
  \centering
  \includegraphics[width=0.8\textwidth]{pcp_pressure_crossover_appendix}
  \caption{%
    Ratio $P_{\rm PCP}(T)/P_{\rm UR}(T)$ versus $T/M$.
    The horizontal dashed line at unity indicates equal pressures.
    The vertical dotted line marks $T_x$ (PCP–UR crossover) and the
    dash-dotted line marks $T_c^{(\Phi)}$ (PCP condensate temperature).
  }
  \label{fig:pcp_crossover}
\end{figure}
```

#### Figure 3 – Appendix (decomposition)

```latex
\begin{figure}[htb]
  \centering
  \includegraphics[width=0.8\textwidth]{pcp_pressure_decomposition_appendix}
  \caption{%
    Absolute values of the individual contributions to the PCP pressure:
    $|C_4|\,T^4$ (solid blue), $|C_7|\,T^7/M^3$ (dashed orange), and their
    sum $|P_{\rm PCP}|$ (dash-dotted black).  The log–log scale reveals the
    different power-law regimes.
  }
  \label{fig:pcp_decomp}
\end{figure}
```

---

### 7. Numerical Results Summary

| Quantity | Value (M = 1) |
|----------|--------------|
| $C_{\rm UR}$ | $\approx 4.058$ |
| $C_4$ | $\approx 2.906$ |
| $C_7$ | $\approx 0.0399$ |
| $T_c^{(\Phi)}$ | $\approx 6.928\,M$ |
| $T_x$ | $\approx 3.07\,M$ |

---

### 8. LaTeX Macros

```latex
% Useful macros for PCP thermodynamics papers
\newcommand{\PCPX}{{\rm PCP}}
\newcommand{\UR}{{\rm UR}}
\newcommand{\Pcp}{P_{\PCPX}}
\newcommand{\Pur}{P_{\UR}}
\newcommand{\Tx}{T_x}
\newcommand{\TcPhi}{T_c^{(\Phi)}}
\newcommand{\Cfour}{C_4}
\newcommand{\Cseven}{C_7}
```

---

### 9. Integration Details

All integrals are computed with `scipy.integrate.quad` using `limit=200`–`300`
to ensure convergence of the slowly-decaying tails.  The relative tolerance
default (`epsrel = 1.49e-8`) is sufficient for the reported 4-digit precision.

---

*This file was auto-generated as part of the `pcp-thermodynamics` package.*
