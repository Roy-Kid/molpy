from pathlib import Path
import molpy as mp
from itertools import islice
from typing import Iterator


class LAMMPSForceFieldReader:

    def __init__(self, scripts: Path | list[Path], data: Path):

        self.scripts = scripts if isinstance(scripts, list) else [scripts]
        self.data = data

    def read(self, system) -> mp.ForceField:

        self.forcefield: mp.ForceField = system.forcefield
        lines = []
        for script in self.scripts:
            with open(script, "r") as f:
                lines.extend(f.readlines())
        with open(self.data, "r") as f:
            lines.extend(f.readlines())
        lines = filter(lambda line: line, map(LAMMPSForceFieldReader.sanitizer, lines))
        n_pairtypes = 0
        for line in lines:
            kw = line[0]

            if kw == "units":
                self.forcefield.unit = line[1]

            elif kw == "bond_style":
                self.read_bondstyle(line[1:])

            elif kw == "pair_style":
                self.read_pairstyle(line[1:])

            elif kw == "angle_style":
                self.read_anglestyle(line[1:])

            elif kw == "dihedral_style":
                self.read_dihedralstyle(line[1:])

            elif kw == "improper_style":
                self.read_improperstyle(line[1:])

            elif kw == "mass":
                self.mass_per_atomtype = self.read_mass_line(line[1:])

            elif kw == "bond_coeff":
                self.read_bondcoeff(self.bondstyle, line[1:])

            elif kw == "angle_coeff":
                self.read_angle_coeff(self.anglestyle, line[1:])

            elif kw == "dihedral_coeff":
                self.read_dihedral_coeff(self.dihedralstyle, line[1:])

            elif kw == "pair_coeff":
                self.read_pair_coeff(self.pairstyle, line[1:])

            elif kw == "pair_modify":
                self.read_pair_modify(line[1:])

            elif kw == "atom_style":
                self.read_atomstyle(line[1:])

            # define in data
            elif kw == "Masses":
                self.read_mass_section(islice(lines, n_atomtypes))

            elif "Coeffs" in line:

                if kw == "Bond":
                    if "#" in line:
                        bondstyle_name = line[line.index("#") + 1]
                    else:
                        bondstyle_name = self.bondstyle.name
                    self.read_bondcoeff_section(
                        bondstyle_name, islice(lines, n_bondtypes)
                    )

                elif kw == "Angle":
                    if "#" in line:
                        anglestyle_name = line[line.index("#") + 1]
                    else:
                        anglestyle_name = self.anglestyle.name
                    self.read_angle_coeff_section(
                        anglestyle_name, islice(lines, n_angletypes)
                    )

                elif kw == "Dihedral":
                    if "#" in line:
                        dihedralstyle_name = line[line.index("#") + 1]
                    else:
                        if isinstance(self.dihedralstyle, list):
                            dihedralstyle_name = None
                        else:
                            dihedralstyle_name = self.dihedralstyle.name
                    self.read_dihedral_coeff_section(
                        dihedralstyle_name, islice(lines, n_dihedraltypes)
                    )

                elif kw == "Improper":
                    if "#" in line:
                        improperstyle_name = line[line.index("#") + 1]
                    else:
                        improperstyle_name = self.improperstyle.name
                    self.read_improper_coeff_section(
                        improperstyle_name, islice(lines, n_impropertypes)
                    )

                elif kw == "Pair":
                    if "#" in line:
                        pairstyle_name = line[line.index("#") + 1]
                    else:
                        pairstyle_name = self.pairstyle.name
                    self.read_pair_coeff_section(
                        pairstyle_name, islice(lines, n_pairtypes)
                    )

            if line[-1] == "types":

                if line[-2] == "atom":
                    n_atomtypes = int(line[0])

                elif line[-2] == "bond":
                    n_bondtypes = int(line[0])

                elif line[-2] == "angle":
                    n_angletypes = int(line[0])

                elif line[-2] == "dihedral":
                    n_dihedraltypes = int(line[0])

                elif line[-2] == "improper":
                    n_impropertypes = int(line[0])

                elif line[-2] == "pair":
                    n_pairtypes = int(line[0])

        # assert self.forcefield.n_atomtypes == n_atomtypes, ValueError(
        #     f"Number of atom types mismatch: {self.forcefield.n_atomtypes} != {n_atomtypes}"
        # )
        # assert self.forcefield.n_bondtypes == n_bondtypes, ValueError(
        #     f"Number of bond types mismatch: {self.forcefield.n_bondtypes} != {n_bondtypes}"
        # )
        # assert self.forcefield.n_angletypes == n_angletypes, ValueError(
        #     f"Number of angle types mismatch: {self.forcefield.n_angletypes} != {n_angletypes}"
        # )
        # assert self.forcefield.n_dihedraltypes == n_dihedraltypes, ValueError(
        #     f"Number of dihedral types mismatch: {self.forcefield.n_dihedraltypes} != {n_dihedraltypes}"
        # )
        # assert self.forcefield.n_impropertypes == n_impropertypes, ValueError(
        #     f"Number of improper types mismatch: {self.forcefield.n_impropertypes} != {n_impropertypes}"
        # )
        # assert self.forcefield.n_pairtypes == n_atomtypes * n_atomtypes, ValueError(
        #     f"Number of pair types mismatch: {self.forcefield.n_pairtypes} != {n_atomtypes * n_atomtypes}"
        # )

        return system

    @staticmethod
    def sanitizer(line: str) -> str:
        return line.split()

    def read_atomstyle(self, line):
        self.atomstyle = self.forcefield.def_atomstyle(line[0])

    def read_bondstyle(self, line):

        if line[0] == "hybrid":
            self.read_bondstyle(line[1:])

        else:
            self.bondstyle = self.forcefield.def_bondstyle(line[0])

    def read_anglestyle(self, line):

        if line[0] == "hybrid":
            self.read_anglestyle(line[1:])

        else:
            self.anglestyle = self.forcefield.def_anglestyle(line[0])

    def read_dihedralstyle(self, line):

        if line[0] == "hybrid":
            results = {}
            style_ = ""
            i = 1
            while i < len(line):
                if not line[i].isdigit():
                    style_ = line[i]
                    results[style_] = []
                else:
                    results[style_].append(line[i])
                i += 1
            for style, coeffs in results.items():
                self.forcefield.def_dihedralstyle(style)
            self.dihedralstyle = self.forcefield.dihedralstyles

        else:
            self.dihedralstyle = self.forcefield.def_dihedralstyle(line[0])

    def read_improperstyle(self, line):

        if line[0] == "hybrid":
            self.read_improperstyle(line[1:])

        else:
            self.improperstyle = self.forcefield.def_improperstyle(line[0])

    def read_pairstyle(self, line):

        if line[0] == "hybrid":
            self.read_pairstyle(line[1:])

        else:
            self.pairstyle = self.forcefield.def_pairstyle(line[0], *line[1:])

    def read_mass_section(self, lines):

        for line in lines:
            type_, m = self.read_mass_line(line)
            self.forcefield.atomstyles[0].get_by(lambda atom: atom.name == str(type_))[
                "mass"
            ] = m

    def read_mass_line(self, line: list[str]):
        return line[0], float(line[1])

    def read_bondcoeff_section(self, stylename: str, lines: Iterator[str]):
        bondstyle = self.forcefield.get_bondstyle(stylename)
        if bondstyle is None:
            bondstyle = self.forcefield.def_bondstyle(stylename)
        for line in lines:
            self.read_bondcoeff(bondstyle, line)

    def read_angle_coeff_section(self, stylename: str, lines: Iterator[str]):
        anglestyle = self.forcefield.get_anglestyle(stylename)
        if anglestyle is None:
            anglestyle = self.forcefield.def_anglestyle(stylename)
        for line in lines:
            self.read_angle_coeff(anglestyle, line)

    def read_dihedral_coeff_section(self, stylename: str, lines: Iterator[str]):
        if stylename is not None:
            dihedralstyle = self.forcefield.get_dihedralstyle(stylename)
            if dihedralstyle is None:
                dihedralstyle = self.forcefield.def_dihedralstyle(stylename)
        else:
            dihedralstyle = None
        for line in lines:
            self.read_dihedral_coeff(dihedralstyle, line)

    def read_improper_coeff_section(self, stylename: str, lines: Iterator[str]):
        improperstyle = self.forcefield.get_improperstyle(stylename)
        if improperstyle is None:
            improperstyle = self.forcefield.def_improperstyle(stylename)
        for line in lines:
            type_id = line[0]
            if type_id.isalpha():
                break
            self.read_improper_coeff(improperstyle, line)

    def read_pair_coeff_section(self, stylename: str, lines: Iterator[str]):
        pairstyle = self.forcefield.get_pairstyle(stylename)
        if pairstyle is None:
            pairstyle = self.forcefield.def_pairstyle(stylename)
        for line in lines:
            # if line[0].isalpha():
            #     break
            # line.insert(0, line[0])  # pair_coeff i j ...
            self.read_pair_coeff(pairstyle, line)

    def read_bondcoeff(self, style, line):

        bond_type_id = line[0]

        if line[1].isalpha():  # hybrid
            bondstyle_name = line[1]
            style = self.forcefield.get_bondstyle(bondstyle_name)
            if style is None:
                style = self.forcefield.def_bondstyle(bondstyle_name)
            coeffs = line[2:]
        else:
            coeffs = line[1:]

        style.def_type(
            bond_type_id,
            None,
            None,
            *coeffs,
        )

    def read_angle_coeff(self, style, line):

        angle_type_id = line[0]

        if line[1].isalpha():  # hybrid
            anglestyle_name = line[1]
            style = self.forcefield.get_anglestyle(anglestyle_name)
            if style is None:
                style = self.forcefield.def_anglestyle(anglestyle_name)
            coeffs = line[2:]
        else:
            coeffs = line[1:]

        style.def_type(
            angle_type_id,
            None,
            None,
            None,
            *coeffs,
        )

    def read_dihedral_coeff(self, style, line):

        dihedral_type_id = line[0]

        if not line[1].isdigit():  # hybrid
            dihedralsyle_name = line[1]
            style = self.forcefield.get_dihedralstyle(dihedralsyle_name)
            if style is None:
                style = self.forcefield.def_dihedralstyle(dihedralsyle_name)
            coeffs = line[2:]
        else:
            coeffs = line[1:]
        style.def_type(
            dihedral_type_id,
            None,
            None,
            None,
            None,
            *coeffs,
        )

    def read_improper_coeff(self, style, line):

        improper_type_id = line[0]

        if line[1].isalpha():  # hybrid
            improperstyle_name = line[1]
            style = self.forcefield.get_improperstyle(improperstyle_name)
            if style is None:
                style = self.forcefield.def_improperstyle(improperstyle_name)
            coeffs = line[2:]
        else:
            coeffs = line[1:]

        style.def_type(
            improper_type_id,
            None,
            None,
            None,
            None,
            *coeffs,
        )

    def read_pair_coeff(self, style, line):

        i, j = line[0], line[1]
        # TODO: unfold * expression
        assert len(self.forcefield.atomstyles) > 0, ValueError("No atom style defined")
        atomstyle = self.forcefield.atomstyles[0]
        atomtype_i = atomstyle.get_by(lambda atom: atom.name == i)
        if atomtype_i is None:
            atomtype_i = atomstyle.def_type(i)
        atomtype_j = atomstyle.get_by(lambda atom: atom.name == j)
        if atomtype_j is None:
            atomtype_j = atomstyle.def_type(j)

        if line[2].isalpha():  # hybrid
            pairstyle_name = line[2]
            style = self.forcefield.get_pairstyle(pairstyle_name)
            if style is None:
                style = self.forcefield.def_pairstyle(pairstyle_name)
            coeffs = line[3:]
        else:
            coeffs = line[2:]

        style.def_type(
            f"{atomtype_i.name}-{atomtype_j.name}",
            atomtype_i,
            atomtype_j,
            *coeffs,
        )

    def read_pair_modify(self, line):

        if line[0] == "pair":
            raise NotImplementedError("pair_modify hybrid not implemented")
        else:
            assert self.forcefield.n_pairstyles == 1, ValueError(
                "pair_modify command requires one pair style"
            )
            pairstyle = self.forcefield.pairstyles[0]

            if "modified" in pairstyle:
                for l in line:
                    if l not in pairstyle["modified"]:
                        pairstyle["modified"].append(l)
            else:
                pairstyle["modified"] = line


class LAMMPSForceFieldWriter:
    def __init__(self, fpath: str | Path):
        self.fpath = fpath

    @staticmethod
    def _write_styles(lines: list[str], styles, style_type: str):
        """
        写出 bond/angle/dihedral/improper style 及其参数
        兼容 parms 命名和 TypeContainer/list 类型
        """
        if not styles:
            return
        if len(styles) == 1:
            style = styles[0]
            if not list(style.types):
                return
            parms = getattr(style, "parms", [])
            lines.append(
                f"{style_type}_style {style.name} {' '.join(f"{x:.3f}" for x in parms)}\n"
            )
            if "modified" in style:
                parms = " ".join(style["modified"])
                lines.append(f"{style_type}_modify {parms}\n")
            for typ in list(style.types):
                parms = getattr(typ, "parms", [])
                lines.append(
                    f"{style_type}_coeff {typ.name} {' '.join(f"{x:.3f}" for x in parms)}\n"
                )
        else:
            style_keywords = " ".join([style.name for style in styles])
            lines.append(f"{style_type}_style hybrid {style_keywords}\n")
            for style in styles:
                for typ in list(style.types):
                    parms = getattr(typ, "parms", [])
    
                    lines.append(
                        f"{style_type}_coeff {typ.name} {style.name} {' '.join(f"{x:.3f}" for x in parms)}\n"
                    )
        lines.append("\n")

    @staticmethod
    def _write_pair_styles(lines: list[str], styles, style_type: str):
        """
        写出 pair_style 及其参数，兼容 parms 命名和 TypeContainer/list 类型
        """
        if not styles:
            return
        if len(styles) == 1:
            style = styles[0]
            parms = getattr(style, "parms", [])
            lines.append(
                f"{style_type}_style {style.name} {' '.join(f"{x:.3f}" for x in parms)}\n"
            )
            if "modified" in style:
                parms = " ".join(style["modified"])
                lines.append(f"{style_type}_modify {parms}\n")
            for typ in list(style.types):
                parms = getattr(typ, "parms", [])

                ats = typ.atomtypes
                if ats[0] == ats[1]:
                    atomnames = ats[0].name
                else:
                    atomnames = f"{ats[0].name} {ats[1].name}"
                lines.append(
                    f"{style_type}_coeff {atomnames} {' '.join(f"{x:.3f}" for x in parms)}\n"
                )
        else:
            style_keywords = " ".join([style.name for style in styles])
            lines.append(f"{style_type}_style hybrid {style_keywords}\n")
            for style in styles:
                for typ in list(style.types):
                    parms = getattr(typ, "parms", [])
    
                    ats = typ.atomtypes
                    if ats[0] == ats[1]:
                        atomnames = ats[0].name
                    else:
                        atomnames = f"{ats[0].name} {ats[1].name}"
                    lines.append(
                        f"{style_type}_coeff {atomnames} {style.name} {' '.join(f"{x:.3f}" for x in parms)}\n"
                    )
        lines.append("\n")

    def write(self, forcefield: mp.ForceField):
        ff = forcefield
        lines = []
        # 写 atom_style
        if ff.atomstyles:
            if len(ff.atomstyles) == 1:
                lines.append(f"# atom_style {ff.atomstyles[0].name}\n")
            else:
                atomstyles = " ".join([atomstyle.name for atomstyle in ff.atomstyles])
                lines.append(f"# atom_style hybrid {atomstyles}\n")
        else:
            raise ValueError("No atom style defined")
        # 可选 mass 输出（如需）
        # for atomstyle in ff.atomstyles:
        #     for atomtype in atomstyle.types:
        #         if 'mass' in atomtype:
        #             lines.append(f"mass {atomtype.name} {atomtype['mass']}\n")
        self._write_styles(lines, ff.bondstyles, "bond")
        self._write_styles(lines, ff.anglestyles, "angle")
        self._write_styles(lines, ff.dihedralstyles, "dihedral")
        self._write_styles(lines, ff.improperstyles, "improper")
        self._write_pair_styles(lines, ff.pairstyles, "pair")
        lines.append("\n")
        with open(self.fpath, "w") as f:
            f.writelines(lines)
