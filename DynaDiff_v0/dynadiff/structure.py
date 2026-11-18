"""Periodic crystal structures with row-vector lattice and fractional coordinates."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import gemmi
import numpy as np
from ase.data import atomic_numbers
from ase.geometry import cell_to_cellpar
from ase.io import read as ase_read

ROOM_TEMPERATURE_K = 293.15
U_TO_B = 8 * np.pi**2


@dataclass(frozen=True)
class Crystal:
    name: str
    cell: np.ndarray
    frac: np.ndarray
    symbols: tuple[str, ...]
    occupancy: np.ndarray
    b_iso: np.ndarray
    reference_temperature_k: float = ROOM_TEMPERATURE_K
    spacegroup: str = ""
    numbers: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        cell = np.asarray(self.cell, dtype=float)
        if (cell.shape != (3, 3) or not np.isfinite(cell).all()
                or np.linalg.matrix_rank(cell) != 3):
            raise ValueError("A crystal needs three independent unit-cell vectors in ångströms.")
        count = len(self.symbols)
        if not count:
            raise ValueError("The structure has no atomic sites.")
        frac = np.asarray(self.frac, dtype=float)
        if frac.shape != (count, 3) or not np.isfinite(frac).all():
            raise ValueError("Each atomic site needs three finite coordinates.")
        occupancy = np.asarray(self.occupancy, dtype=float)
        if (occupancy.shape != (count,) or not np.isfinite(occupancy).all()
                or np.any((occupancy < 0) | (occupancy > 1))):
            raise ValueError("Site occupancies must be numbers between zero and one.")
        b_iso = np.asarray(self.b_iso, dtype=float)
        if b_iso.shape != (count,) or np.any(np.isinf(b_iso) | (b_iso < 0)):
            raise ValueError("A reported atomic displacement factor must be finite and nonnegative.")
        numbers = np.array([atomic_numbers.get(symbol, 0) for symbol in self.symbols])
        if np.any(numbers == 0):
            bad = self.symbols[np.flatnonzero(numbers == 0)[0]]
            raise ValueError(f"Unrecognized chemical element '{bad}' in an atomic site.")
        object.__setattr__(self, "cell", cell)
        object.__setattr__(self, "frac", frac % 1)
        object.__setattr__(self, "occupancy", occupancy)
        object.__setattr__(self, "b_iso", b_iso)
        object.__setattr__(self, "numbers", numbers)

    @property
    def cart(self):
        """Cartesian positions in ångströms, within one conventional unit cell."""
        return self.frac @ self.cell


def crystal_cif(crystal: Crystal) -> str:
    """Write one expanded unit cell in P1 with fractional coordinates."""
    cell, frac = crystal.cell, crystal.frac
    if np.linalg.det(cell) < 0:
        # CIF reconstructs a right-handed cell from lengths and angles. Swap
        # two axes and their coordinate columns to preserve the physical motif.
        order = [1, 0, 2]
        cell, frac = cell[order], frac[:, order]

    block = gemmi.cif.Block("DynaDiff")
    block.set_pair("_chemical_name_common", gemmi.cif.quote(crystal.name))
    parameters = ("length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma")
    for key, value in zip(parameters, cell_to_cellpar(cell)):
        block.set_pair("_cell_" + key, f"{value:.16g}")
    block.set_pair("_space_group_name_H-M_alt", "'P 1'")
    block.set_pair("_space_group_IT_number", "1")
    block.set_pair("_diffrn_ambient_temperature", f"{crystal.reference_temperature_k:.16g}")
    sites = block.init_loop("_atom_site_", [
        "label", "type_symbol", "fract_x", "fract_y", "fract_z", "occupancy", "B_iso_or_equiv",
    ])
    for index, (symbol, position, occupancy, b_iso) in enumerate(
            zip(crystal.symbols, frac, crystal.occupancy, crystal.b_iso), start=1):
        sites.add_row([f"{symbol}{index}", symbol, *(f"{v:.16g}" for v in position),
                       f"{occupancy:.16g}", "?" if np.isnan(b_iso) else f"{b_iso:.16g}"])
    return block.as_string()


def _from_gemmi(small, name, temperature):
    cell = np.asarray(small.cell.orth.mat).T
    sites = small.get_all_unit_cell_sites()
    symbols = tuple("H" if site.element.name == "D" else site.element.name for site in sites)
    return Crystal(
        name, cell, [list(site.fract) for site in sites], symbols,
        [site.occ for site in sites],
        [U_TO_B * site.u_iso for site in sites],
        temperature, small.spacegroup.xhm() if small.spacegroup else small.spacegroup_hm,
    )


def _numeric_column(block, tag, count):
    values = block.find_values(tag)
    if not values:
        return np.full(count, np.nan)
    if len(values) != count:
        raise ValueError(f"{tag} has a different length from the atomic site list.")
    return np.array([gemmi.cif.as_number(gemmi.cif.as_string(v)) for v in values])


def _set_cif_symmetry(small, block):
    # magCIF extends each spatial operation with a time-reversal sign.
    # Electrostatic/X-ray diffraction uses its spatial part and centering.
    def magnetic_ops(name):
        prefix = "_space_group_symop_magn_" + name
        column = block.find_values(prefix + ".xyz") or block.find_values(prefix + "_xyz")
        return [gemmi.Op(",".join(gemmi.cif.as_string(value).split(",")[:3])) for value in column]

    if not small.symops:
        operations = magnetic_ops("operation")
        if operations:
            centering = magnetic_ops("centering") or [gemmi.Op("x,y,z")]
            small.symops = sorted({(c * op).wrap().triplet() for c in centering for op in operations})
    # Prefer explicit XYZ operations, then Hall, H-M, then the IT number.
    small.determine_and_set_spacegroup("S.H2N")
    if small.symops or small.spacegroup_hall:
        if not small.spacegroup and not small.cell.images:
            raise ValueError("The CIF symmetry operations or Hall symbol could not be identified.")
        return
    if not small.spacegroup:
        if small.spacegroup_hm or small.spacegroup_number:
            raise ValueError("The CIF space-group symbol or number could not be identified.")
        return  # A coordinate list without symmetry information describes a full cell.
    tag = (block.find_value("_space_group_IT_coordinate_system_code")
           or block.find_value("_space_group.IT_coordinate_system_code"))
    code = gemmi.cif.as_string(tag or "?").lower()
    if code:
        for group in gemmi.spacegroup_table():
            ext = group.ext.strip("\x00").lower()
            qualifier = group.qualifier.lower()
            if group.number == small.spacegroup.number and code in (ext + qualifier, ext + (qualifier or "abc")):
                small.spacegroup_hm = group.xhm()
                small.determine_and_set_spacegroup("2")
                return
        raise ValueError(f"Unrecognized coordinate system '{code}' for {small.spacegroup.hm}.")


def _read_core_cif(block, path):
    for field_name in ("length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"):
        tag = "_cell_" + field_name
        value = gemmi.cif.as_number(gemmi.cif.as_string(block.find_value(tag) or "?"))
        if not np.isfinite(value) or value <= 0 or (field_name.startswith("angle") and value >= 180):
            raise ValueError(f"A valid {tag} is required.")

    # A few structure exporters provide element symbols but omit site labels.
    types = block.find_values("_atom_site_type_symbol")
    if types and not block.find_values("_atom_site_label"):
        labels = [f"site{i + 1}" for i in range(len(types))]
        if types.get_loop():
            types.get_loop().add_columns(["_atom_site_label"], ".")
            label_column = block.find_values("_atom_site_label")
            for i, label in enumerate(labels):
                label_column[i] = label
        else:
            block.set_pair("_atom_site_label", labels[0])
    small = gemmi.make_small_structure_from_block(block)
    _set_cif_symmetry(small, block)
    count = len(small.sites)
    if not count:
        raise ValueError("The CIF needs atomic site labels or element symbols.")
    fractional = all(block.find_values("_atom_site_fract_" + axis) for axis in "xyz")
    if not fractional:
        cart = np.column_stack([_numeric_column(block, "_atom_site_Cartn_" + ax, count) for ax in "xyz"])
        if not np.isfinite(cart).all():
            raise ValueError("The CIF needs three fractional or Cartesian coordinates for each site.")
        for site, xyz in zip(small.sites, cart):
            site.fract = small.cell.fractionalize(gemmi.Position(*xyz))

    # Gemmi substitutes zero for missing U. Keep missing values distinct from
    # explicitly measured zero, including when just one atom lacks a B factor.
    b = _numeric_column(block, "_atom_site_B_iso_or_equiv", count)
    u = _numeric_column(block, "_atom_site_U_iso_or_equiv", count)
    b = np.where(np.isfinite(b), b, u * U_TO_B)
    for site, value in zip(small.sites, b):
        if not np.isfinite(value) and site.aniso.nonzero():
            value = small.cell.calculate_u_eq(site.aniso) * U_TO_B
        site.u_iso = value / U_TO_B

    temperature = ROOM_TEMPERATURE_K
    for tag in ("_diffrn_ambient_temperature", "_cell_measurement_temperature"):
        value = gemmi.cif.as_number(gemmi.cif.as_string(block.find_value(tag) or "?"))
        if np.isfinite(value) and value > 0:
            temperature = value
            break
    name = gemmi.cif.as_string(block.find_value("_chemical_name_common") or "?")
    return _from_gemmi(small, name or block.name or path.stem, temperature)


def _read_macromolecular(model, name):
    if not model.cell.is_crystal():
        raise ValueError("The file needs crystallographic unit-cell parameters.")
    if not len(model):
        raise ValueError("The file has no atom coordinates.")
    small = gemmi.mx_to_sx_structure(model, 0)
    return _from_gemmi(small, name, ROOM_TEMPERATURE_K)


_ASE_FORMATS = {
    ".vasp": "vasp", ".poscar": "vasp", ".contcar": "vasp",
    ".xsf": "xsf", ".extxyz": "extxyz", ".xyz": "extxyz",
    ".pwi": "espresso-in",
}


def load_structure(path) -> Crystal:
    """Read the first structural CIF block or first periodic coordinate frame."""
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"Structure file not found: {path}")
    filename = Path(path.stem) if path.suffix.lower() == ".gz" else path
    suffix = filename.suffix.lower()
    try:
        if suffix in (".cif", ".mcif", ".mmcif"):
            for block in gemmi.cif.read(str(path)):
                if block.find_values("_atom_site.fract_x") or block.find_values("_atom_site.Cartn_x"):
                    return _read_macromolecular(gemmi.make_structure_from_block(block), block.name)
                if block.find_values("_atom_site_fract_x") or block.find_values("_atom_site_Cartn_x"):
                    return _read_core_cif(block, filename)
            raise ValueError("No data block contains atomic coordinates.")
        if suffix in (".pdb", ".ent"):
            return _read_macromolecular(gemmi.read_structure(str(path)), filename.stem)
        fmt = _ASE_FORMATS.get(suffix)
        if filename.name.upper().split(".")[0] in ("POSCAR", "CONTCAR"):
            fmt = "vasp"
        if not fmt:
            raise ValueError("Choose a CIF, POSCAR/CONTCAR, XSF, periodic extended XYZ, PDB or PWI file.")
        atoms = ase_read(str(path), format=fmt, index=0)
        if atoms.cell.rank != 3 or not np.all(atoms.pbc):
            raise ValueError("A three-dimensional periodic unit cell is required. For XYZ, include Lattice in the extended XYZ header.")
        count = len(atoms)
        return Crystal(
            filename.stem, atoms.cell.array, atoms.get_scaled_positions(wrap=False),
            tuple(atoms.get_chemical_symbols()), atoms.arrays.get("occupancy", np.ones(count)),
            atoms.arrays.get("b_iso", np.full(count, np.nan)),
        )
    except Exception as exc:
        raise ValueError(f"Could not read {path.name}. {str(exc) or type(exc).__name__}") from exc
