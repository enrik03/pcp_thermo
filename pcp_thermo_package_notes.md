# Corrected PCP thermodynamic-figure package

## Files
- `pcp_pressure_comparison_main.(png|pdf)`
- `pcp_pressure_crossover_appendix.(png|pdf)`
- `pcp_pressure_decomposition_appendix.(png|pdf)`
- `pcp_thermo_figures_corrected.py`

## Numerical checks
\[
\mathcal{C}_7 \approx 0.1930959,
\qquad
\mathcal{C}_4 \approx 0.00345436,
\qquad
\frac{\mathcal{C}_7}{\mathcal{C}_4} \approx 55.899.
\]
\[
\epsilon_{\rm enh}=-0.125,
\qquad
\epsilon_{\rm sup}=+0.125,
\qquad
T_c^{(\Phi)} = \left(\frac{\mathcal{C}_4}{\mathcal{C}_7}\right)^{1/3} \approx 0.2615,
\qquad
T_\times = \left(\frac{K_{\rm UR}}{\mathcal{C}_7}\right)^{1/3} \approx 1.4393.
\]
At \(T=10^{-1/2}\approx 0.316\):
\[
\frac{P_{\parallel}(\epsilon_s=-0.125)}{P_{\rm free}}\approx 1.480,
\qquad
\frac{P_{\parallel}(\epsilon_s=+0.125)}{P_{\rm free}}\approx 0.675,
\qquad
\frac{P_{\Phi_0}}{P_{\rm free}}\approx 1.566.
\]
At \(T=5\):
\[
\frac{P_{\parallel}(\epsilon_s=-0.125)}{P_{\rm free}}\approx 1.025,
\qquad
\frac{P_{\parallel}(\epsilon_s=+0.125)}{P_{\rm free}}\approx 0.975,
\qquad
\frac{P_{\Phi_0}}{P_{\rm free}}\approx 1.00014.
\]

## Caption 1 (main text)
**Figure 1.** Thermodynamic pressures in the post-Carroll-Pauli sector for the representative dimensionless choice \(\hbar=m=c=k_B=1\), \(\xi=0.1\), \(\phi=0.5\), \(\Phi_0=0.3\), and \(g_d=2\). Panel (a) compares the free PCP pressure \(P_{\rm free}=\mathcal{C}_7T^7\) with the two constant-shift parallel branches, the constant-projection branch \(P_{\Phi_0}=\mathcal{C}_7T^7+\mathcal{C}_4T^4\), and the standard ultra-relativistic fermion law \(P_{\rm UR}=K_{\rm UR}T^4\). The free PCP coefficient is \(\mathcal{C}_7\approx0.1931\), whereas \(\mathcal{C}_4\approx3.45\times10^{-3}\), so the constant-projection correction is perturbative throughout the plotted range. The branch with \(\epsilon_s=-0.125\) enhances the pressure, while \(\epsilon_s=+0.125\) suppresses it. Panel (b) shows the ratios with respect to the free PCP pressure. At \(T=10^{-1/2}\approx0.316\), one finds \(P_{\parallel}(\epsilon_s=-0.125)/P_{\rm free}\approx1.480\), \(P_{\parallel}(\epsilon_s=+0.125)/P_{\rm free}\approx0.675\), and \(P_{\Phi_0}/P_{\rm free}\approx1.566\). By contrast, at \(T=5\) the three deformed PCP branches already lie close to the free result: the corresponding ratios are approximately \(1.025\), \(0.975\), and \(1.00014\). The guide lines emphasize the distinct scalings \(T^7\) and \(T^4\). The standard ultra-relativistic curve is included only in the absolute-pressure panel, since its ratio to \(P_{\rm free}\) is not confined to a narrow neighborhood of unity over this temperature interval.

## Caption 2 (appendix, crossover)
**Figure 2.** Comparison between the free PCP pressure \(P_{\rm free}=\mathcal{C}_7T^7\), the constant-projection branch \(P_{\Phi_0}=\mathcal{C}_7T^7+\mathcal{C}_4T^4\), and the standard ultra-relativistic fermion pressure \(P_{\rm UR}=K_{\rm UR}T^4\), using the same parameter values as in Fig. 1. The plot isolates the difference between the PCP scaling exponent and the conventional ultra-relativistic one. Since \(\mathcal{C}_7\approx0.1931\) and \(K_{\rm UR}\approx0.5757\), the free PCP and standard ultra-relativistic pressures become equal at the crossover temperature \(T_\times=(K_{\rm UR}/\mathcal{C}_7)^{1/3}\approx1.439\). Below this point the \(T^4\) law is larger, whereas above it the PCP \(T^7\) behavior dominates rapidly. The constant-projection branch lies slightly above the free PCP curve because the \(\mathcal{C}_4T^4\) piece is positive for \(\xi\Phi_0>0\), but the two PCP-based curves become nearly indistinguishable at moderate and high temperature. This figure is the cleanest visual summary of how the PCP thermodynamics differs from the standard ultra-relativistic fermionic equation of state, even before discussing the detailed helicity-dependent deformations.

## Caption 3 (appendix, decomposition)
**Figure 3.** Internal decomposition of the constant-projection branch \(P_{\Phi_0}=\mathcal{C}_7T^7+\mathcal{C}_4T^4\) for \(\xi=0.1\), \(\Phi_0=0.3\), and \(g_d=2\). The dashed orange curve represents the leading PCP contribution \(\mathcal{C}_7T^7\), while the dotted purple curve shows the subleading ultra-relativistic-like correction \(\mathcal{C}_4T^4\). Their formal equality occurs at \(T_c^{(\Phi)}=(\mathcal{C}_4/\mathcal{C}_7)^{1/3}\approx0.2615\), marked in the figure. Because \(\mathcal{C}_7/\mathcal{C}_4\approx55.9\), the \(T^7\) term already dominates for \(T\gtrsim0.3\), which is why the total pressure and the leading PCP contribution track each other almost exactly over most of the displayed range. The purpose of this figure is therefore not to claim a physically sharp phase transition, but to visualize the bookkeeping of the perturbative expression and the relative weight of the two scalings. If the effective PCP description is restricted away from the infrared, the equality point should be viewed as a formal indicator of where the two terms balance inside the truncated expression, rather than as an independently robust low-temperature prediction.

## Main-text insertion snippet
```latex
The thermodynamic behavior of the deformed PCP branches is illustrated in
Fig.~\ref{fig:pcp_pressure_comparison}. For the representative values
$\xi=0.1$, $\phi=0.5$, $\Phi_0=0.3$, and $g_d=2$, the two constant-shift
branches correspond to $\epsilon_s=\mp0.125$, with the branch
$\epsilon_s=-0.125$ enhancing the pressure and $\epsilon_s=+0.125$
suppressing it. In the same range, the constant-projection branch preserves
the leading PCP scaling and only adds a subleading $T^4$ correction.

\begin{figure}[htbp]
\centering
\includegraphics[width=\textwidth]{pcp_pressure_comparison_main.pdf}
\caption{Thermodynamic pressures in the post-Carroll-Pauli sector for the
representative dimensionless choice $\hbar=m=c=k_B=1$, $\xi=0.1$,
$\phi=0.5$, $\Phi_0=0.3$, and $g_d=2$. Panel (a) compares the free PCP
pressure $P_{\mathrm{free}}=\mathcal{C}_7T^7$ with the two constant-shift
parallel branches, the constant-projection branch
$P_{\Phi_0}=\mathcal{C}_7T^7+\mathcal{C}_4T^4$, and the standard
ultra-relativistic fermion law $P_{\mathrm{UR}}=K_{\mathrm{UR}}T^4$. The
branch with $\epsilon_s=-0.125$ enhances the pressure, while
$\epsilon_s=+0.125$ suppresses it. Panel (b) shows the ratios relative to
the free PCP pressure. At $T=10^{-1/2}\approx0.316$, one finds
$P_{\parallel}(\epsilon_s=-0.125)/P_{\mathrm{free}}\approx1.480$,
$P_{\parallel}(\epsilon_s=+0.125)/P_{\mathrm{free}}\approx0.675$, and
$P_{\Phi_0}/P_{\mathrm{free}}\approx1.566$, whereas at $T=5$ the three
deformed PCP branches already lie close to the free result. The guide lines
highlight the distinct $T^7$ and $T^4$ scalings.}
\label{fig:pcp_pressure_comparison}
\end{figure}
```

## Appendix D block
```latex
\section*{Appendix D: Graphical comparison of PCP thermodynamic pressures}
\addcontentsline{toc}{section}{Appendix D: Graphical comparison of PCP thermodynamic pressures}
\label{app:pcp_graphics}

Figures~\ref{fig:pcp_pressure_crossover} and
\ref{fig:pcp_pressure_decomposition} summarize two complementary numerical
illustrations of the thermodynamic formulas derived in
Secs.~3.2.5--3.2.6. The first compares the free PCP scaling with the standard
ultra-relativistic fermion law and marks the crossover temperature
$T_\times=(K_{\mathrm{UR}}/\mathcal{C}_7)^{1/3}\approx1.439$. The second
separates the two contributions entering the constant-projection branch,
$P_{\Phi_0}=\mathcal{C}_7T^7+\mathcal{C}_4T^4$, and marks the formal balance
point $T_c^{(\Phi)}=(\mathcal{C}_4/\mathcal{C}_7)^{1/3}\approx0.2615$.

\begin{figure}[htbp]
\centering
\includegraphics[width=0.82\textwidth]{pcp_pressure_crossover_appendix.pdf}
\caption{Comparison between the free PCP pressure
$P_{\mathrm{free}}=\mathcal{C}_7T^7$, the constant-projection branch
$P_{\Phi_0}=\mathcal{C}_7T^7+\mathcal{C}_4T^4$, and the standard
ultra-relativistic fermion pressure $P_{\mathrm{UR}}=K_{\mathrm{UR}}T^4$.
The free PCP and standard ultra-relativistic curves cross at
$T_\times\approx1.439$. Below this point the $T^4$ law is larger, while
above it the PCP $T^7$ behavior dominates rapidly.}
\label{fig:pcp_pressure_crossover}
\end{figure}

\begin{figure}[htbp]
\centering
\includegraphics[width=0.82\textwidth]{pcp_pressure_decomposition_appendix.pdf}
\caption{Internal decomposition of the constant-projection branch
$P_{\Phi_0}=\mathcal{C}_7T^7+\mathcal{C}_4T^4$. The two terms become equal at
$T_c^{(\Phi)}\approx0.2615$. This balance point is a formal feature of the
perturbative expression and should not be overinterpreted if the effective
PCP regime is restricted away from the infrared.}
\label{fig:pcp_pressure_decomposition}
\end{figure}
```

## Short textual insertions
After Eq.~\eqref{pressure_final_phi}:
```latex
For the representative values $\xi=0.1$, $\phi=0.5$, and $g_d=2$, the two
constant-shift branches correspond to $\epsilon_s=\mp0.125$. Numerically,
at $T=10^{-1/2}\simeq0.316$ one finds
$P_{\parallel}(\epsilon_s=-0.125)/P_{\mathrm{free}}\simeq1.480$ and
$P_{\parallel}(\epsilon_s=+0.125)/P_{\mathrm{free}}\simeq0.675$, while at
$T=5$ the corresponding ratios are already close to unity,
$1.025$ and $0.975$, respectively.
```

After Eq.~\eqref{P_constproj_compact}:
```latex
For $\xi=0.1$, $\Phi_0=0.3$, and $g_d=2$, one has
$\mathcal{C}_7\simeq0.1931$, $\mathcal{C}_4\simeq3.45\times10^{-3}$, and
$\mathcal{C}_7/\mathcal{C}_4\simeq55.9$. The formal equality
$\mathcal{C}_7T^7=\mathcal{C}_4T^4$ therefore occurs at
$T_c^{(\Phi)}=(\mathcal{C}_4/\mathcal{C}_7)^{1/3}\simeq0.2615$, so that for
$T\gtrsim0.3$ the $T^7$ PCP contribution already dominates.
```

In the Discussion:
```latex
A graphical comparison also makes the thermodynamic hierarchy transparent.
For $\xi=0.1$, $\phi=0.5$, $\Phi_0=0.3$, and $g_d=2$, the deformed PCP
branches remain close to the free PCP pressure at moderate and high
temperature, whereas the standard ultra-relativistic law $P_{\mathrm{UR}}
\propto T^4$ departs strongly from the PCP scaling
$P_{\mathrm{free}}\propto T^7$. In particular, the free PCP and standard
ultra-relativistic pressures become equal at
$T_\times=(K_{\mathrm{UR}}/\mathcal{C}_7)^{1/3}\simeq1.439$.
```
