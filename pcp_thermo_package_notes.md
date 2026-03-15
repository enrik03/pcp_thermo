# PCP thermodynamic numerical details

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


```latex
For the representative values $\xi=0.1$, $\phi=0.5$, and $g_d=2$, the two
constant-shift branches correspond to $\epsilon_s=\mp0.125$. Numerically,
at $T=10^{-1/2}\simeq0.316$ one finds
$P_{\parallel}(\epsilon_s=-0.125)/P_{\mathrm{free}}\simeq1.480$ and
$P_{\parallel}(\epsilon_s=+0.125)/P_{\mathrm{free}}\simeq0.675$, while at
$T=5$ the corresponding ratios are already close to unity,
$1.025$ and $0.975$, respectively.
```


