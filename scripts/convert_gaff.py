#!/usr/bin/env python3
"""Convert GAFF force field files (gaff.prm + gaff.dat) to molpy XML format.

Usage:
    python scripts/convert_gaff.py <gaff_prm> <gaff_dat> -o <output_xml>

The gaff.prm file contains SMARTS atom type definitions.
The gaff.dat file contains bond, angle, dihedral, improper, and VDW parameters.

Unit conversions (AMBER real units -> molpy XML storage):
- Bond k:     kcal/mol/Å² (stored as-is, AMBER real units)
- Bond r0:    Å (stored as-is)
- Angle k:    kcal/mol/rad² (stored as-is)
- Angle θ0:   degrees -> radians
- Dihedral:   kcal/mol, degrees -> radians for phase
- VDW R*:     Å (R_min/2) -> σ in Å via σ = R* * 2^(5/6) ... no, store as R* directly
- VDW ε:      kcal/mol (stored as-is)
"""

import argparse
import math
import re
import xml.etree.ElementTree as ET
from xml.dom import minidom


# ---- Element lookup from atomic number ----
ATOMIC_NUM_TO_ELEMENT = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O",
    9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",
    16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 21: "Sc", 22: "Ti",
    23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu",
    30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr",
    37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo", 43: "Tc",
    44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn",
    51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba", 57: "La",
    58: "Ce", 59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd",
    65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu",
    72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir", 78: "Pt",
    79: "Au", 80: "Hg", 81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At",
    86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th", 91: "Pa", 92: "U",
    93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk", 98: "Cf", 99: "Es",
    100: "Fm", 101: "Md", 102: "No", 103: "Lr",
}

# Mass to element mapping for GAFF atom types
MASS_TO_ELEMENT = {
    1.008: "H", 12.01: "C", 14.01: "N", 16.00: "O",
    19.00: "F", 30.97: "P", 32.06: "S", 35.45: "Cl",
    79.90: "Br", 126.9: "I",
}

# Map GAFF type prefix to element
TYPE_PREFIX_TO_ELEMENT = {
    "h": "H", "c": "C", "n": "N", "o": "O", "f": "F", "p": "P",
    "s": "S", "cl": "Cl", "br": "Br", "i": "I",
}

# All elements by symbol for direct matching
ELEMENTS_BY_SYMBOL = {
    "He": "He", "Li": "Li", "Be": "Be", "B": "B", "Ne": "Ne",
    "Na": "Na", "Mg": "Mg", "Al": "Al", "Si": "Si", "Ar": "Ar",
    "K": "K", "Ca": "Ca", "Sc": "Sc", "Ti": "Ti", "V": "V",
    "Cr": "Cr", "Mn": "Mn", "Fe": "Fe", "Co": "Co", "Ni": "Ni",
    "Cu": "Cu", "Zn": "Zn", "Ga": "Ga", "Ge": "Ge", "As": "As",
    "Se": "Se", "Kr": "Kr", "Rb": "Rb", "Sr": "Sr", "Y": "Y",
    "Zr": "Zr", "Nb": "Nb", "Mo": "Mo", "Tc": "Tc", "Ru": "Ru",
    "Rh": "Rh", "Pd": "Pd", "Ag": "Ag", "Cd": "Cd", "In": "In",
    "Sn": "Sn", "Sb": "Sb", "Te": "Te", "Xe": "Xe", "Cs": "Cs",
    "Ba": "Ba", "La": "La", "Ce": "Ce", "Pr": "Pr", "Nd": "Nd",
    "Pm": "Pm", "Sm": "Sm", "Eu": "Eu", "Gd": "Gd", "Tb": "Tb",
    "Dy": "Dy", "Ho": "Ho", "Er": "Er", "Tm": "Tm", "Yb": "Yb",
    "Lu": "Lu", "Hf": "Hf", "Ta": "Ta", "W": "W", "Re": "Re",
    "Os": "Os", "Ir": "Ir", "Pt": "Pt", "Au": "Au", "Hg": "Hg",
    "Tl": "Tl", "Pb": "Pb", "Bi": "Bi", "Po": "Po", "At": "At",
    "Rn": "Rn", "Fr": "Fr", "Ra": "Ra", "Ac": "Ac", "Th": "Th",
    "Pa": "Pa", "U": "U", "Np": "Np", "Pu": "Pu", "Am": "Am",
    "Cm": "Cm", "Bk": "Bk", "Cf": "Cf", "Es": "Es", "Fm": "Fm",
    "Md": "Md", "No": "No", "Lr": "Lr",
}


def gaff_type_to_element(type_name: str, mass: float) -> str:
    """Determine element from GAFF type name and mass."""
    # Check if type_name is a direct element symbol (e.g., He, Li, Fe)
    if type_name in ELEMENTS_BY_SYMBOL:
        return ELEMENTS_BY_SYMBOL[type_name]

    # Try mass lookup
    for m, elem in MASS_TO_ELEMENT.items():
        if abs(mass - m) < 0.1:
            return elem

    # Try type prefix
    for prefix, elem in sorted(TYPE_PREFIX_TO_ELEMENT.items(), key=lambda x: -len(x[0])):
        if type_name.startswith(prefix):
            return elem

    # Special cases
    if type_name == "X":
        return "*"
    if type_name == "DU":
        return "X"  # dummy

    return "X"  # unknown


def parse_gaff_prm(filepath: str) -> list[dict]:
    """Parse gaff.prm SMARTS atom type definitions.

    Returns list of dicts: {smarts, type, desc, context_smarts}
    where context_smarts is the full SMARTS pattern including environment.
    """
    rules = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("//"):
                continue

            if not line.startswith("atom "):
                continue

            # Format: atom <SMARTS_pattern>  <optional_env>  <type> "<description>"
            # The SMARTS pattern starts with [ and the type is a short string
            rest = line[5:].strip()

            # Find the SMARTS pattern - it starts with [ and we need to handle
            # the environment part after it (e.g., [#1X1]O means H bonded to O)
            # The type name comes after whitespace
            # Strategy: split from the right to find type and description

            # Extract description in quotes if present
            desc = ""
            desc_match = re.search(r'"([^"]*)"', rest)
            if desc_match:
                desc = desc_match.group(1)
                rest = rest[:desc_match.start()].strip()

            # Now rest is: <SMARTS_with_env>  <type>  <optional_parens>
            # Remove any trailing parenthesized content like (-*)
            rest = re.sub(r'\s+\(-\*\)\s*$', '', rest)
            rest = rest.strip()

            # Split by whitespace - the last token is the type
            parts = rest.split()
            if len(parts) < 2:
                continue

            type_name = parts[-1]
            full_smarts = " ".join(parts[:-1])

            rules.append({
                "smarts": full_smarts,
                "type": type_name,
                "desc": desc,
            })

    return rules


def _parse_amber_type_field(field: str) -> list[str]:
    """Parse AMBER fixed-width type fields separated by '-'.

    AMBER uses fixed-width fields where type names can have trailing spaces.
    E.g., 'X -c -c -X' contains types ['X', 'c', 'c', 'X'].
    """
    return [t.strip() for t in field.split("-")]


def parse_gaff_dat(filepath: str) -> dict:
    """Parse gaff.dat parameter file.

    Returns dict with keys: atom_types, bonds, angles, dihedrals, impropers, vdw
    """
    with open(filepath) as f:
        lines = f.readlines()

    result = {
        "atom_types": [],
        "bonds": [],
        "angles": [],
        "dihedrals": [],
        "impropers": [],
        "vdw": [],
    }

    # Skip first line (title)
    i = 1

    # Parse atom types section (until blank line)
    while i < len(lines) and lines[i].strip():
        line = lines[i].strip()
        # Format: <type>  <mass>  <polarizability>  <description>
        parts = line.split()
        if len(parts) >= 3:
            type_name = parts[0]
            try:
                mass = float(parts[1])
                polarizability = float(parts[2])
                desc = " ".join(parts[3:]) if len(parts) > 3 else ""
                result["atom_types"].append({
                    "type": type_name,
                    "mass": mass,
                    "polarizability": polarizability,
                    "desc": desc,
                })
            except ValueError:
                pass
        i += 1

    # Skip blank line and hydrophilic atom line
    i += 1
    while i < len(lines) and not re.match(r'\S+\s*-\s*\S+\s+\d', lines[i]):
        i += 1

    # Parse bonds section
    # Format: xx-yy  <k>  <r0>  <source>
    # Types are separated by '-', may have spaces around them
    while i < len(lines) and lines[i].strip():
        line = lines[i]
        # Split into type-field and numeric part
        # The type field ends where the first number starts after the last '-'
        match = re.match(r'(.+?)\s+(\d+\.?\d*)\s+(\d+\.?\d*)', line)
        if match:
            type_field = match.group(1)
            types = _parse_amber_type_field(type_field)
            if len(types) == 2:
                k = float(match.group(2))
                r0 = float(match.group(3))
                result["bonds"].append({
                    "type1": types[0], "type2": types[1],
                    "k": k, "r0": r0,
                })
        i += 1

    # Skip blank line
    i += 1

    # Parse angles section
    # Format: xx-yy-zz  <k>  <theta0>  <source>
    while i < len(lines) and lines[i].strip():
        line = lines[i]
        match = re.match(r'(.+?)\s+(\d+\.?\d*)\s+(\d+\.?\d*)', line)
        if match:
            type_field = match.group(1)
            types = _parse_amber_type_field(type_field)
            if len(types) == 3:
                k = float(match.group(2))
                theta0_deg = float(match.group(3))
                result["angles"].append({
                    "type1": types[0], "type2": types[1], "type3": types[2],
                    "k": k, "theta0": theta0_deg,
                })
        i += 1

    # Skip blank line
    i += 1

    # Parse dihedrals section
    # Format: X -a1-a2-X  <divider>  <barrier>  <phase>  <periodicity>
    # Note: negative periodicity means more terms follow for same dihedral
    while i < len(lines) and lines[i].strip():
        line = lines[i]
        # Match type field (4 types separated by -), then numeric values
        match = re.match(
            r'(.+?)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+(-?[\d.]+)',
            line
        )
        if match:
            type_field = match.group(1)
            types = _parse_amber_type_field(type_field)
            if len(types) == 4:
                divider = int(match.group(2))
                barrier = float(match.group(3))
                phase = float(match.group(4))
                periodicity = float(match.group(5))

                result["dihedrals"].append({
                    "type1": types[0], "type2": types[1],
                    "type3": types[2], "type4": types[3],
                    "divider": divider,
                    "barrier": barrier,
                    "phase": phase,
                    "periodicity": periodicity,
                })
        i += 1

    # Skip blank line
    i += 1

    # Parse impropers section
    # Format: X -a1-a2-a3  <barrier>  <phase>  <periodicity>
    while i < len(lines) and lines[i].strip():
        line = lines[i]
        match = re.match(
            r'(.+?)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)',
            line
        )
        if match:
            type_field = match.group(1)
            types = _parse_amber_type_field(type_field)
            if len(types) == 4:
                barrier = float(match.group(2))
                phase = float(match.group(3))
                periodicity = float(match.group(4))

                result["impropers"].append({
                    "type1": types[0], "type2": types[1],
                    "type3": types[2], "type4": types[3],
                    "barrier": barrier,
                    "phase": phase,
                    "periodicity": periodicity,
                })
        i += 1

    # Skip to VDW section (look for "MOD4      RE")
    while i < len(lines):
        if "MOD4" in lines[i] and "RE" in lines[i]:
            i += 1
            break
        i += 1

    # Parse VDW section
    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith("END"):
            break
        # Format: <type>  <R*>  <epsilon>  <description>
        parts = line.split()
        if len(parts) >= 3:
            try:
                type_name = parts[0]
                rstar = float(parts[1])
                epsilon = float(parts[2])
                result["vdw"].append({
                    "type": type_name,
                    "rstar": rstar,
                    "epsilon": epsilon,
                })
            except ValueError:
                pass
        i += 1

    return result


def consolidate_dihedrals(raw_dihedrals: list[dict]) -> list[dict]:
    """Consolidate multi-term dihedrals (negative periodicity indicates more terms).

    AMBER format uses negative periodicity to indicate that more terms follow
    for the same dihedral. We consolidate these into single entries with
    multiple terms.
    """
    consolidated = []
    current = None

    for d in raw_dihedrals:
        periodicity = d["periodicity"]
        is_continuation = periodicity < 0
        abs_period = abs(periodicity)

        term = {
            "periodicity": int(abs_period),
            "barrier": d["barrier"] / d["divider"],  # V/n
            "phase": d["phase"],
        }

        if current is not None and (
            current["type1"] == d["type1"] and
            current["type2"] == d["type2"] and
            current["type3"] == d["type3"] and
            current["type4"] == d["type4"] and
            current["_has_more"]
        ):
            # This is a continuation of the previous dihedral
            current["terms"].append(term)
            current["_has_more"] = is_continuation
        else:
            # Save previous if exists
            if current is not None:
                del current["_has_more"]
                consolidated.append(current)
            current = {
                "type1": d["type1"],
                "type2": d["type2"],
                "type3": d["type3"],
                "type4": d["type4"],
                "terms": [term],
                "_has_more": is_continuation,
            }

    if current is not None:
        del current["_has_more"]
        consolidated.append(current)

    return consolidated


def generate_xml(prm_rules, dat_params, output_path):
    """Generate GAFF XML file."""

    root = ET.Element("ForceField")
    root.set("name", "GAFF")
    root.set("version", "1.4")
    root.set("combining_rule", "lorentz-berthelot")

    # Build type->mass and type->element maps from dat
    type_mass = {}
    type_element = {}
    for at in dat_params["atom_types"]:
        type_mass[at["type"]] = at["mass"]
        type_element[at["type"]] = gaff_type_to_element(at["type"], at["mass"])

    # Build VDW lookup
    vdw_map = {}
    for v in dat_params["vdw"]:
        vdw_map[v["type"]] = v

    # ---- AtomTypes ----
    atomtypes_elem = ET.SubElement(root, "AtomTypes")
    for rule in prm_rules:
        type_name = rule["type"]
        mass = type_mass.get(type_name, 0.0)
        element = type_element.get(type_name, "X")

        type_elem = ET.SubElement(atomtypes_elem, "Type")
        type_elem.set("name", type_name)
        type_elem.set("class", type_name)
        if element != "X" and element != "*":
            type_elem.set("element", element)
        type_elem.set("mass", f"{mass:.4f}")
        type_elem.set("def", rule["smarts"])
        if rule["desc"]:
            type_elem.set("desc", rule["desc"])

    # ---- HarmonicBondForce ----
    bonds_elem = ET.SubElement(root, "HarmonicBondForce")
    for bond in dat_params["bonds"]:
        bond_elem = ET.SubElement(bonds_elem, "Bond")
        bond_elem.set("type1", bond["type1"])
        bond_elem.set("type2", bond["type2"])
        bond_elem.set("class1", bond["type1"])
        bond_elem.set("class2", bond["type2"])
        # Store in AMBER real units (kcal/mol/Å², Å)
        bond_elem.set("k", f"{bond['k']:.4f}")
        bond_elem.set("length", f"{bond['r0']:.4f}")

    # ---- HarmonicAngleForce ----
    angles_elem = ET.SubElement(root, "HarmonicAngleForce")
    for angle in dat_params["angles"]:
        angle_elem = ET.SubElement(angles_elem, "Angle")
        angle_elem.set("type1", angle["type1"])
        angle_elem.set("type2", angle["type2"])
        angle_elem.set("type3", angle["type3"])
        angle_elem.set("class1", angle["type1"])
        angle_elem.set("class2", angle["type2"])
        angle_elem.set("class3", angle["type3"])
        angle_elem.set("k", f"{angle['k']:.4f}")
        # Convert theta0 from degrees to radians
        theta0_rad = angle["theta0"] * math.pi / 180.0
        angle_elem.set("angle", f"{theta0_rad:.6f}")

    # ---- PeriodicTorsionForce ----
    dihedrals_consolidated = consolidate_dihedrals(dat_params["dihedrals"])
    torsion_elem = ET.SubElement(root, "PeriodicTorsionForce")
    for dih in dihedrals_consolidated:
        proper_elem = ET.SubElement(torsion_elem, "Proper")
        proper_elem.set("type1", dih["type1"])
        proper_elem.set("type2", dih["type2"])
        proper_elem.set("type3", dih["type3"])
        proper_elem.set("type4", dih["type4"])
        proper_elem.set("class1", dih["type1"])
        proper_elem.set("class2", dih["type2"])
        proper_elem.set("class3", dih["type3"])
        proper_elem.set("class4", dih["type4"])

        for j, term in enumerate(dih["terms"], 1):
            # Store barrier in kcal/mol, phase in radians
            phase_rad = term["phase"] * math.pi / 180.0
            proper_elem.set(f"periodicity{j}", str(term["periodicity"]))
            proper_elem.set(f"k{j}", f"{term['barrier']:.6f}")
            proper_elem.set(f"phase{j}", f"{phase_rad:.6f}")

    # ---- PeriodicImproperForce ----
    improper_elem = ET.SubElement(root, "PeriodicImproperForce")
    for imp in dat_params["impropers"]:
        imp_elem = ET.SubElement(improper_elem, "Improper")
        imp_elem.set("type1", imp["type1"])
        imp_elem.set("type2", imp["type2"])
        imp_elem.set("type3", imp["type3"])
        imp_elem.set("type4", imp["type4"])
        imp_elem.set("class1", imp["type1"])
        imp_elem.set("class2", imp["type2"])
        imp_elem.set("class3", imp["type3"])
        imp_elem.set("class4", imp["type4"])
        phase_rad = imp["phase"] * math.pi / 180.0
        imp_elem.set("periodicity1", str(int(imp["periodicity"])))
        imp_elem.set("k1", f"{imp['barrier']:.6f}")
        imp_elem.set("phase1", f"{phase_rad:.6f}")

    # ---- NonbondedForce ----
    nb_elem = ET.SubElement(root, "NonbondedForce")
    nb_elem.set("coulomb14scale", "0.8333")  # AMBER default 1/1.2
    nb_elem.set("lj14scale", "0.5")  # AMBER default

    for vdw in dat_params["vdw"]:
        atom_elem = ET.SubElement(nb_elem, "Atom")
        atom_elem.set("type", vdw["type"])
        atom_elem.set("charge", "0.0")  # charges come from AM1-BCC, not FF
        # Store R* (R_min/2) and epsilon directly in AMBER real units
        # sigma = R* * 2^(5/6) for conversion, but we store R* for clarity
        atom_elem.set("sigma", f"{vdw['rstar']:.4f}")
        atom_elem.set("epsilon", f"{vdw['epsilon']:.4f}")

    # Pretty print
    xml_str = ET.tostring(root, encoding="unicode")
    dom = minidom.parseString(xml_str)
    pretty_xml = dom.toprettyxml(indent="  ", encoding=None)

    # Remove extra blank lines from minidom
    lines = pretty_xml.split("\n")
    cleaned = [l for l in lines if l.strip()]
    # Re-add XML declaration at top
    output = "\n".join(cleaned) + "\n"

    with open(output_path, "w") as f:
        f.write(output)

    print(f"Generated {output_path}")
    print(f"  Atom types: {len(prm_rules)}")
    print(f"  Bonds: {len(dat_params['bonds'])}")
    print(f"  Angles: {len(dat_params['angles'])}")
    print(f"  Dihedrals: {len(dihedrals_consolidated)}")
    print(f"  Impropers: {len(dat_params['impropers'])}")
    print(f"  VDW params: {len(dat_params['vdw'])}")


def main():
    parser = argparse.ArgumentParser(description="Convert GAFF files to molpy XML")
    parser.add_argument("prm", help="Path to gaff.prm (SMARTS definitions)")
    parser.add_argument("dat", help="Path to gaff.dat (parameters)")
    parser.add_argument("-o", "--output", default="gaff.xml", help="Output XML path")
    args = parser.parse_args()

    prm_rules = parse_gaff_prm(args.prm)
    dat_params = parse_gaff_dat(args.dat)
    generate_xml(prm_rules, dat_params, args.output)


if __name__ == "__main__":
    main()
