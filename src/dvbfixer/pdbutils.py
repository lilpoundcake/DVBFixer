"""Shared PDB utilities for dvbfixer.

CONECT record remapping and other PDB format helpers used across modules.
"""


def build_serial_map(pdb_path):
    """Read a PDB file and build (chain, resid, atomname) -> serial map."""
    serial_map = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                serial = int(line[6:11])
                chain = line[21]
                resid = line[22:26].strip()
                atomname = line[12:16].strip()
                serial_map[(chain, resid, atomname)] = serial
    return serial_map


def remap_conect_records(input_path, new_serial_map):
    """Read CONECT from input PDB, remap serials to match a new output.

    Builds old_serial -> (chain, resid, atomname) from input, then maps
    to new serials via new_serial_map. Returns list of remapped CONECT lines.
    """
    with open(input_path) as f:
        lines = f.readlines()

    old_serial_to_key = {}
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            serial = int(line[6:11])
            chain = line[21]
            resid = line[22:26].strip()
            atomname = line[12:16].strip()
            old_serial_to_key[serial] = (chain, resid, atomname)

    result = []
    for line in lines:
        if not line.startswith('CONECT'):
            continue
        # Parse CONECT serials (5-char wide fields after "CONECT")
        serials = []
        s = line[6:]
        while len(s) >= 5:
            chunk = s[:5].strip()
            if chunk:
                serials.append(int(chunk))
            s = s[5:]
        if len(serials) < 2:
            continue

        new_serials = []
        all_mapped = True
        for old_s in serials:
            key = old_serial_to_key.get(old_s)
            if key and key in new_serial_map:
                new_serials.append(new_serial_map[key])
            else:
                all_mapped = False
                break

        if all_mapped and len(new_serials) >= 2:
            conect = f"CONECT{new_serials[0]:5d}"
            for ns in new_serials[1:]:
                conect += f"{ns:5d}"
            result.append(conect.ljust(80) + "\n")

    return result


def append_before_end(output_path, extra_lines):
    """Insert lines before the END record in a PDB file."""
    if not extra_lines:
        return
    with open(output_path) as f:
        lines = f.readlines()
    with open(output_path, 'w') as f:
        for line in lines:
            if line.startswith("END"):
                f.writelines(extra_lines)
            f.write(line)
