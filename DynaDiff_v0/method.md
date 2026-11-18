# How DynaDiff simulates dynamical electron diffraction

DynaDiff propagates a plane electron wave through a perfect, parallel-sided crystal slab. The beam enters along the current structure view, perpendicular to the slab. At the entrance all intensity is in the direct beam, 000. Inside the crystal, any retained beam can scatter into any other. Their exit amplitudes interfere, so a reflection can brighten and then fade as thickness increases. The implementation is in `dynadiff/physics.py`. [Bloch-wave propagation, equations 1–8](https://journals.iucr.org/a/issues/2024/02/00/tw5006/index.html)

**Try it on silicon.** Open the bundled example, choose **View → [110]**, and keep the default 200 kV, 20 nm, 1.25 Å⁻¹ radius and 101-beam limit. With **Multiple scattering** off, the (0 0 2) reflection is extinct because the contributions from silicon sites cancel. Turn the effect on and about **2.8% of the incident intensity** reaches (0 0 2) through other beams. The reflection table shows this change on a linear scale.

## From structure to pattern

1. **Build the beam set.** For direct lattice matrix $A$ with rows in Å, a reflection $h=(h,k,l)$ has reciprocal vector $g=hA^{-T}$ in Å⁻¹. Its dot products with the screen-right and screen-up axes give its plot coordinates. The solver enumerates beams inside the requested reciprocal radius and checks the actual motif for conventional centering translations. Among accessible forward beams, it selects those with the smallest absolute excitation errors first, using reciprocal length to break ties. A motif extinction can still acquire intensity through indirect scattering.

2. **Compute the elastic structure factors.** For each beam difference $d$, sum the neutral-atom electron form factors over the symmetry-expanded unit cell. With fractional position $x_j$, occupancy $o_j$, and isotropic displacement parameter $B_j$ in Å², the sum is

   $$F_d=\sum_j o_j f_j(|g_d|/2)\,e^{-B_j|g_d|^2/4}\,e^{-2\pi i d\cdot x_j}.$$

   Gemmi supplies the five-Gaussian electron factors in Å. [Gemmi scattering tables](https://gemmi.readthedocs.io/en/stable/scattering.html#electron-form-factors)

3. **Propagate through the thickness.** Let $\lambda$ be the electron wavelength, $\gamma$ the relativistic mass multiplier, $\Omega$ the cell volume and $g_z$ the beam-parallel component. The propagation matrix has diagonal entries $2\pi s_g$, where $s_g=-g_z-\lambda |g|^2/2$, and off-diagonal entries $\gamma\lambda F_{h-h'}/\Omega$. The uniform mean-potential phase drops out of intensities. For thickness $t$ in Å, the exit wave is $a(t)=\exp(iMt)a(0)$. Each plotted intensity is $|a_h(t)|^2$ as a fraction of the incident intensity. The Hermitian matrix conserves total intensity within the chosen beam set. [Structure matrix and extinction distance](https://journals.iucr.org/a/issues/2024/02/00/tw5006/index.html)

The default basis contains up to 101 beams, and the user can request up to 9,999 including (000). Up to 401 retained beams, DynaDiff diagonalizes the propagation matrix directly. For larger bases it applies the same beam-to-beam couplings as a zero-padded FFT convolution over integer Miller indices, with the excitation term applied separately. A Chebyshev expansion with Bessel coefficients propagates the incident wave to the exit surface. This avoids storing a dense 9,999-by-9,999 matrix. The beam count can fall below the requested limit when the reciprocal radius contains fewer candidates or the cutoff would split a symmetry-equivalent shell. [Jacobi–Anger expansion](https://dlmf.nist.gov/10.12)

## What the switches change

| Effect | Enabled by default | When disabled |
| --- | --- | --- |
| Relativistic correction | Voltage determines the relativistic de Broglie wavelength and mass multiplier. | The same voltage sets the classical wavelength, with $\gamma=1$. |
| Debye–Waller factor | Each site's $B$ scales linearly with the entered temperature relative to the structure's measurement temperature. Its default reference is 293.15 K. Missing values use an estimated $B=0.5$ Å² at 293.15 K. | All thermal damping factors equal one. |
| Multiple scattering | All beam pairs couple throughout the slab. | A first-Born calculation uses the same beam set and the finite-thickness $\mathrm{sinc}^2(\pi t s_g)$ rocking curve, with $\mathrm{sinc}(x)=\sin(x)/x$. The direct beam stays undepleted. |
| Ewald-sphere curvature | The excitation error includes the wavelength-dependent curvature term. | A flat reciprocal-space section uses $s_g=-g_z$. |

Temperature scaling is a high-temperature approximation. At low temperature it omits zero-point motion. Thermal damping affects coherent amplitudes only. The model omits diffuse and inelastic scattering, absorption, bonding charge, beam convergence, and variation in thickness or orientation. The electron form-factor table covers $|g_d|/2$ through 2 Å⁻¹. For quantitative comparisons, increase both the reciprocal radius and the beam limit until the reflections of interest stabilize, keeping beam differences within this tabulated range. [Relativistic electron optics](https://onlinelibrary.wiley.com/iucr/itc/Ba/ch5o2v0001/) · [Temperature-dependent dynamical scattering](https://arxiv.org/abs/2209.09465)

## Kinematic reference patterns

Both kinematic modes compute one scattering event per reflection and apply the slab rocking curve. The X-ray mode uses neutral-atom IT92 form factors, the forward Ewald sphere, unpolarized Thomson scattering, and the exit-beam direction cosine for normal flux through the slab. A monochromatic X-ray pattern can be weak when the view is far from a Bragg orientation. Each table row reports an individual reflection, including when several reflections share a detector coordinate. [Gemmi X-ray factors](https://gemmi.readthedocs.io/en/stable/scattering.html#x-ray-form-factors) · [Kinematic scattering by nanocrystals](https://journals.iucr.org/j/issues/2023/01/00/xx5007/index.html)

## From reflection intensities to a detector image

**Simulated detector** spreads each reflection across pixels. The default is a circular 2D Gaussian, an empirical model used to fit individual spots in selected-area electron diffraction. Its full width at half maximum (FWHM) is adjustable in Å⁻¹. To approximate a particular experiment, measure the FWHM of an isolated spot in reciprocal-space units and enter that value. The alternative **Airy** profile models the central maximum and rings produced by a coherent circular aperture. [Gaussian fits to experimental diffraction spots](https://doi.org/10.1093/jmicro/dfx121) · [Electron diffraction through a circular selected-area aperture](https://pubmed.ncbi.nlm.nih.gov/23545433/)

Each reflection contributes its calculated integrated intensity to a normalized profile on a 513 × 513 detector. Contributions at the same pixel add, and tails that extend beyond the field of view are cropped. The scale bar gives reciprocal distance. The colorbar gives incident-beam intensity fraction per pixel. Instrument noise, diffuse background, sensor saturation and camera-specific calibration are outside this image model.

The outlier setting clips the brightest percentage of pixels for display. Set it to zero to use the full maximum. Log scaling displays four decades below the upper limit. Colormap, reversal, clipping and log scaling affect the rendered PNG. The reflection table and CSV retain linear, integrated intensities.
