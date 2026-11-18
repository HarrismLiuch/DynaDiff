# DynaDiff

DynaDiff is a Python desktop app for exploring single-crystal diffraction. Rotate a crystal, and the simulated pattern updates along the same beam direction. You can compare kinematic X-ray and electron diffraction, or see how multiple scattering redistributes electron intensity as a crystal gets thicker.

## Launch on macOS

1. Clone or download this repository and open the `DynaDiff_v0` folder in Finder.
2. With **Python 3.10 or later** installed, double-click `setup_venv.command`. The first setup needs an internet connection to download dependencies into a local `.venv`. Wait until Terminal says **DynaDiff is ready**.
3. Double-click `DynaDiff.command`. The app opens with a silicon crystal viewed along [001].

You can also launch from Terminal at the repository root. The launcher runs setup automatically on first use.

```sh
./DynaDiff_v0/DynaDiff.command
```

## Try it on silicon

1. In **Atomic structure → View**, choose **[110]**. Double-click the atomic view to fit the whole cell.
2. Keep **Dynamical electron** at the default **200 kV** and **20 nm**. In **Reflections**, select the row with indices **(0 0 2)**. Its spot is highlighted in the pattern. The **I / I₀** column reads about **0.028**, or **2.8% of the incident intensity**.
3. Uncheck **Multiple scattering**. The (0 0 2) reflection drops effectively to zero. Turn it back on and the spot returns as electrons can reach that reflection through other beams.

![DynaDiff showing silicon along [110] at 200 kV and 20 nm, with the (0 0 2) spot circled and its intensity selected in the reflection table.](docs/silicon-110.png)

## Work with your own crystal

| Task | How to do it |
| --- | --- |
| Change the interface theme | Click **Theme** in the top bar and choose **Vintage**, **White** (default), **Studio Light**, **Midnight**, or **High Contrast**. DynaDiff remembers your choice. |
| Load a structure | Click **Open structure…**. Supported formats include CIF, crystallographic PDB, VASP POSCAR/CONTCAR, XSF, periodic extended XYZ and Quantum ESPRESSO input (`.pwi`). |
| Set the beam direction | Drag the atomic view to tilt, Shift-drag to roll, or enter three zone-axis indices or Cartesian components and click **Set**. Scroll to zoom. |
| Choose the scattering model | In **Simulation**, choose **Kinematic X-ray** and enter a wavelength, or choose an electron mode and set the accelerating voltage. In **Dynamical electron**, vary thickness or toggle the scattering effects to compare patterns. |
| View detector peaks | Below **Diffraction**, switch **Indexed spots** to **Simulated detector**. Open **Settings…** to change the peak profile, width, colormap or intensity scale. Scroll over a spot to zoom and double-click to reset. |
| Save results | In **Reflections → Export**, save the displayed pattern as PNG or the linear reflection intensities as CSV. In **Atomic structure → View**, save the atomic view as PNG or the complete unit cell as CIF. |

For all controls, see the [application guide](DynaDiff_v0/README.md). For equations, assumptions and limits of the dynamical model, see the [simulation method](DynaDiff_v0/method.md).
