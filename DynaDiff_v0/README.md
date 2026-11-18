# DynaDiff

DynaDiff is a Python desktop application for exploring how a single crystal scatters a beam along the direction shown in the atomic preview.

## Launch on macOS

1. Install Python 3.10 or later if needed. Double-click `setup_venv.command` to create the project's local `.venv` and download its dependencies.
2. Double-click `DynaDiff.command`. The workspace opens with diamond-cubic silicon viewed along [001].
3. Use **Open structure…** to load a crystal. The file picker remembers the last successfully opened folder across launches. From Terminal you can also run `./DynaDiff.command "/path/to/crystal.cif"`.

The launcher runs setup automatically on first use. To distribute DynaDiff, share the sources, example, requirements and both executable command files. Git excludes the local environment. Each recipient runs setup on their own machine. Set `DYNADIFF_PYTHON` to a Python executable to choose the interpreter when creating a new environment.

Choose **Theme** in the top bar or **Appearance → Theme** in the menu bar to switch between **Vintage**, **White** (default), **Studio Light**, **Midnight**, and **High Contrast**. The choice applies immediately and is saved for the next launch. Atom colors and the detector colormap have their own controls.

## Simulate a pattern

- **Load a crystal.** Supported inputs include CIF, crystallographic PDB, VASP POSCAR/CONTCAR, XSF, periodic extended XYZ and Quantum ESPRESSO input (`.pwi`). Extended XYZ needs a three-dimensional periodic lattice. CIF symmetry operations expand the atomic sites into the unit cell.
- **Orient the crystal.** Drag the atomic structure to tilt it, Shift-drag to roll, or enter a zone axis [u v w] or Cartesian beam vector and click **Set**. Both views use the resulting orientation. Rotation preserves magnification by default. Scroll to zoom or double-click to fit the current view. Enable **View → Auto-fit on rotation** to keep the complete cell in the field of view.
- **Color the atoms.** Under **View → Atom colors**, choose Mineral, Jmol or CPK. Use **Customize element…** to pick an individual element’s color. The legend follows your choices.
- **Choose a model.** Select **Kinematic X-ray**, **Kinematic electron** or **Dynamical electron**. Set the X-ray wavelength or electron voltage and the crystal thickness along the beam. The pattern updates automatically.
- **Explore dynamical effects.** Multiple scattering, Debye–Waller attenuation, relativistic corrections and Ewald-sphere curvature start enabled. Change thickness to see intensity transfer between reflections. The default beam limit is 101, adjustable from 2 to 9,999 including the direct beam. Increase the beam limit and reciprocal radius to check whether a reflection's intensity has converged.
- **Render a detector image.** Beneath the diffraction pattern, switch **Indexed spots** to **Simulated detector**. Open **Settings…** to choose a Gaussian or Airy profile, set its FWHM, or browse colormaps by provider and category. Matplotlib follows its reference gallery, CMasher follows the type and lightness order in its introduction, and Crameri follows its scientific colour map gallery. Typing a name filters the list in the same order. Each preview follows **Reverse colormap**. Bright-pixel clipping and log scaling control the intensity range. [Matplotlib reference](https://matplotlib.org/stable/gallery/color/colormap_reference.html) · [CMasher overview](https://cmasher.readthedocs.io/user/introduction.html#colormap-overview) · [Crameri gallery](https://github.com/callumrollo/cmcrameri#extra-instructions)
- **Inspect a pattern.** Scroll over a reflection to zoom around the pointer in either display mode. Double-click to restore the full pattern, or use the zoom commands under **Display**. The scale bar follows the zoom. Click a spot to select its Miller indices in the table.
- **Export.** Under the atomic structure's **View** menu or **File**, export the current atomic view as PNG or the complete unit cell as CIF. The CIF lists all sites in P1 and preserves occupancies and displacement factors. In the reflections panel, **Export** saves the current diffraction view as PNG, including its zoom, or all integrated reflection intensities as CSV.

The dynamical calculation models coherent elastic scattering in a periodic crystal. The equations, temperature approximation and model limits are described in [method.md](method.md).
