"""sbtools — a small toolbox for Homework 2.

Everything here uses only numpy, matplotlib and Biopython. Nothing needs DSSP,
PyMOL or an internet connection beyond downloading your structure once.

You are not expected to read every line. You are expected to know what each
function returns, because the homework script calls them by name.
"""
import os
import urllib.request

import numpy as np
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.SASA import ShrakeRupley

STRUCT_DIR = "structures"

# Theoretical maximum solvent-accessible surface area per residue, A^2.
# Tien M.Z. et al. PLoS ONE. 2013;8:e80635. doi:10.1371/journal.pone.0080635
MAX_ASA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167,
           "GLU": 223, "GLN": 225, "GLY": 104, "HIS": 224, "ILE": 197,
           "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159,
           "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174}

AA3 = set(MAX_ASA)


# --------------------------------------------------------------- getting files
def fetch_pdb(pdb_id, fmt="pdb"):
    """Download one PDB entry. Returns the local path."""
    os.makedirs(STRUCT_DIR, exist_ok=True)
    pdb_id = pdb_id.lower()
    ext = "cif" if fmt == "cif" else "pdb"
    path = os.path.join(STRUCT_DIR, f"{pdb_id}.{ext}")
    if not os.path.exists(path):
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.{ext}"
        urllib.request.urlretrieve(url, path)
    return path


def fetch_alphafold(uniprot):
    """Download the AlphaFold DB model for one UniProt accession.

    The model URL is read from the AlphaFold API rather than guessed, because the
    version number in the file name changes with every database release.
    """
    import json
    os.makedirs(STRUCT_DIR, exist_ok=True)
    path = os.path.join(STRUCT_DIR, f"AF_{uniprot}.pdb")
    if not os.path.exists(path):
        api = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"
        with urllib.request.urlopen(api) as r:
            meta = json.loads(r.read().decode())
        urllib.request.urlretrieve(meta[0]["pdbUrl"], path)
    return path


def load(path, model=0):
    """Parse a .pdb or .cif file and return one model."""
    parser = MMCIFParser(QUIET=True) if path.endswith(".cif") else PDBParser(QUIET=True)
    return list(parser.get_structure("s", path))[model]


def n_models(path):
    parser = MMCIFParser(QUIET=True) if path.endswith(".cif") else PDBParser(QUIET=True)
    return len(list(parser.get_structure("s", path)))


def header_facts(path):
    """A few things worth reading out of the file itself, as a dict of strings."""
    out = {}
    if path.endswith(".cif"):
        return out
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith("EXPDTA"):
            out["method"] = line[10:].strip()
        elif line.startswith("REMARK   2 RESOLUTION."):
            out["resolution"] = line[22:].strip()
        elif line.startswith("TITLE"):
            out["title"] = (out.get("title", "") + " " + line[10:].strip()).strip()
        elif line.startswith("SSBOND"):
            out["ssbond"] = str(int(out.get("ssbond", "0")) + 1)
        elif line.startswith("ATOM"):
            break
    return out


# ------------------------------------------------------------------- residues
def residues(model, chain_id=None, standard_only=True):
    """Amino-acid residues with a complete backbone, as a flat list."""
    out = []
    for ch in model:
        if chain_id is not None and ch.id != chain_id:
            continue
        for r in ch:
            if r.id[0] != " ":
                continue
            if standard_only and r.get_resname() not in AA3:
                continue
            if not all(a in r for a in ("N", "CA", "C")):
                continue
            out.append(r)
    return out


def chain_ids(model, min_len=20):
    return [ch.id for ch in model if len(residues(model, ch.id)) >= min_len]


def ca_coords(model, chain_id):
    return np.array([r["CA"].coord for r in residues(model, chain_id)], float)


# --------------------------------------------------------------- torsion angles
def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / (np.linalg.norm(b1) + 1e-9)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def phi_psi(model, chain_id):
    """(phi, psi) for every residue that has both, plus bookkeeping.

    Returns a dict of arrays: phi, psi, resname, resid, next_resname.
    """
    res = residues(model, chain_id)
    phi, psi, nm, rid, nxt = [], [], [], [], []
    for i in range(1, len(res) - 1):
        a, b, c = res[i - 1], res[i], res[i + 1]
        if b.id[1] != a.id[1] + 1 or c.id[1] != b.id[1] + 1:
            continue          # a gap in the chain: the angle is undefined
        phi.append(dihedral(a["C"].coord, b["N"].coord, b["CA"].coord, b["C"].coord))
        psi.append(dihedral(b["N"].coord, b["CA"].coord, b["C"].coord, c["N"].coord))
        nm.append(b.get_resname()); rid.append(b.id[1]); nxt.append(c.get_resname())
    return dict(phi=np.array(phi), psi=np.array(psi), resname=np.array(nm),
                resid=np.array(rid), next_resname=np.array(nxt))


def omega(model, chain_id):
    """Omega dihedral per peptide bond; values near 0 are cis bonds."""
    res = residues(model, chain_id)
    out = []
    for i in range(len(res) - 1):
        a, b = res[i], res[i + 1]
        if b.id[1] != a.id[1] + 1:
            continue
        out.append((b.id[1], b.get_resname(),
                    dihedral(a["CA"].coord, a["C"].coord, b["N"].coord, b["CA"].coord)))
    return out


# ------------------------------------------ Ramachandran allowed-region lookup
class Rama:
    """The sterically allowed region of the (phi, psi) plane.

    Loaded from `rama_allowed.npz`, computed for Lecture 2 by building an
    Ala-Ala-Ala tripeptide at every (phi, psi) on a 5 degree grid and measuring the
    worst hard-sphere overlap over the atom pairs the dihedrals control, using the
    contact criteria of Ramachandran, Ramakrishnan & Sasisekharan (1963).

    Overlap 0 means 'atoms exactly touch'; positive means they interpenetrate.

    Two tolerances are used, and the difference between them is itself worth
    understanding:

    tol = 0.30 A       the band shown in the lecture. It admits 42% of the (phi, psi)
                       square and contains 96% of residues in the mixed set of
                       structures used there.
    outlier_tol = 0.50 A   the threshold this homework uses to CALL a residue an
                       outlier. It has to be looser, because hard spheres are an
                       approximation: real atoms are slightly soft, and the strict
                       band clips the lower edge of the alpha-helix basin. In a
                       helical protein the 0.30 band alone flags 20-25% of perfectly
                       ordinary helix residues. At 0.50 A, >= 99% of residues in
                       well-refined structures are inside, which is the same order
                       as MolProbity's allowed-region statistics.
    """

    def __init__(self, path="rama_allowed.npz", tol=0.30, outlier_tol=0.50):
        z = np.load(path, allow_pickle=True)
        self.viol = z["viol"]
        self.step = 360 // self.viol.shape[0]
        self.grid = np.arange(-180, 180, self.step)
        self.tol = tol
        self.outlier_tol = outlier_tol

    def overlap(self, phi, psi):
        """Worst hard-sphere overlap (A) at these angles."""
        i = ((np.asarray(phi) + 180) // self.step).astype(int) % len(self.grid)
        j = ((np.asarray(psi) + 180) // self.step).astype(int) % len(self.grid)
        return self.viol[j, i]

    def is_allowed(self, phi, psi):
        """Inside the strict hard-sphere band used in the lecture."""
        return self.overlap(phi, psi) <= self.tol

    def is_outlier(self, phi, psi):
        """Outside even the loose band — worth looking at one by one."""
        return self.overlap(phi, psi) > self.outlier_tol

    def background(self, ax, color="#e7ebf0", loose="#f2f4f7"):
        """Shade both bands on an existing matplotlib axis."""
        ax.contourf(self.grid, self.grid, self.viol, levels=[-99, self.outlier_tol],
                    colors=[loose], zorder=0)
        ax.contourf(self.grid, self.grid, self.viol, levels=[-99, self.tol],
                    colors=[color], zorder=0)
        ax.set_xlim(-180, 180); ax.set_ylim(-180, 180); ax.set_aspect("equal")
        ax.set_xticks([-180, -90, 0, 90, 180]); ax.set_yticks([-180, -90, 0, 90, 180])
        ax.set_xlabel("φ  (degrees)"); ax.set_ylabel("ψ  (degrees)")
        return ax


# ------------------------------------------------------- secondary structure
def _amide_h(res):
    """Estimated backbone amide H positions, DSSP convention."""
    out = {}
    for i in range(1, len(res)):
        prev, cur = res[i - 1], res[i]
        if "C" not in prev or "O" not in prev or "N" not in cur:
            continue
        d = prev["C"].coord - prev["O"].coord
        n = np.linalg.norm(d)
        if n > 1e-6:
            out[i] = cur["N"].coord + d / n
    return out


def _hbond_energy(res, i, j, H):
    """DSSP electrostatic hydrogen-bond energy, kcal/mol, donor i to acceptor j."""
    if i not in H or res[i].get_resname() == "PRO":
        return 0.0
    N, h = res[i]["N"].coord, H[i]
    C, O = res[j]["C"].coord, res[j]["O"].coord
    r = lambda a, b: max(float(np.linalg.norm(a - b)), 1e-3)
    return 0.084 * (1 / r(O, N) + 1 / r(C, h) - 1 / r(O, h) - 1 / r(C, N)) * 332.0


def secondary_structure(model, chain_id, cutoff=-0.5):
    """DSSP-style H/E/L string for one chain, plus its backbone H-bonds.

    H = alpha-helix, E = beta-strand (bridge), L = everything else. Returns
    (ss_string, resids, hbonds) where hbonds is a list of (donor, acceptor, energy).
    """
    res = [r for r in residues(model, chain_id) if "O" in r]
    n = len(res)
    H = _amide_h(res)
    ca = np.array([r["CA"].coord for r in res])
    bonds = {}
    for i in range(n):
        for j in range(n):
            if abs(i - j) < 2 or np.linalg.norm(ca[i] - ca[j]) > 9.0:
                continue
            e = _hbond_energy(res, i, j, H)
            if e < cutoff:
                bonds[(i, j)] = e
    has = lambda i, j: (i, j) in bonds
    ss = ["L"] * n
    for i in range(n - 5):
        if has(i + 4, i) and has(i + 5, i + 1):
            for k in range(i + 1, i + 5):
                ss[k] = "H"
    bridge = [False] * n
    for i in range(1, n - 1):
        for j in range(1, n - 1):
            if abs(i - j) < 3:
                continue
            anti = (has(i, j) and has(j, i)) or (has(j + 1, i - 1) and has(i + 1, j - 1))
            para = (has(i, j - 1) and has(j + 1, i)) or (has(j, i - 1) and has(i + 1, j))
            if anti or para:
                bridge[i] = bridge[j] = True
    for i in range(n):
        if bridge[i] and ss[i] != "H":
            ss[i] = "E"
    out = list(ss)
    i = 0
    while i < n:
        j = i
        while j < n and out[j] == out[i]:
            j += 1
        if (out[i] == "H" and j - i < 4) or (out[i] == "E" and j - i < 2):
            out[i:j] = ["L"] * (j - i)
        i = j
    resids = np.array([r.id[1] for r in res])
    hb = [(res[i].id[1], res[j].id[1], round(e, 2)) for (i, j), e in bonds.items()]
    return "".join(out), resids, hb


# --------------------------------------------------------------- burial / SASA
def relative_sasa(model, chain_id=None, isolate_chain=False):
    """Relative solvent accessibility per residue: 0 buried, ~1 fully exposed.

    Waters, ions and ligands are removed first — otherwise ordered waters sit on the
    surface and make everything look buried.

    isolate_chain=False  the residue is surrounded by every other protein chain in
                         the file, so interface residues come out buried.
    isolate_chain=True   only `chain_id` is kept, so the same residues come out
                         exposed. The difference between the two is exactly the
                         'an interface is a core turned inside out' slide.

    Keys are (chain, resid, resname). Reference maximum areas: Tien et al. 2013.
    """
    import copy
    m = copy.deepcopy(model)
    for ch in list(m):
        if isolate_chain and chain_id is not None and ch.id != chain_id:
            m.detach_child(ch.id)
            continue
        for r in list(ch):
            if r.id[0] != " " or r.get_resname() not in AA3:
                ch.detach_child(r.id)
    ShrakeRupley().compute(m, level="R")
    out = {}
    for r in residues(m, chain_id):
        mx = MAX_ASA.get(r.get_resname())
        if mx:
            out[(r.get_parent().id, r.id[1], r.get_resname())] = r.sasa / mx
    return out


def interface_residues(model, chain_a, chain_b, cutoff=4.5):
    """Residues of chain_a with any atom within cutoff of chain_b."""
    atoms_b = np.array([a.coord for r in residues(model, chain_b) for a in r])
    hits = []
    for r in residues(model, chain_a):
        xyz = np.array([a.coord for a in r])
        d = np.linalg.norm(xyz[:, None, :] - atoms_b[None, :, :], axis=-1)
        if d.min() <= cutoff:
            hits.append((r.id[1], r.get_resname(), round(float(d.min()), 2)))
    return hits


# ------------------------------------------------------------ superposition
def kabsch(P, Q):
    """Optimal superposition of P onto Q. Returns (rmsd, P_fitted).

    Kabsch W. Acta Crystallogr A. 1976;32:922-923. doi:10.1107/S0567739476001873
    """
    pc, qc = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    Pf = (P - pc) @ R.T + qc
    return float(np.sqrt(((Pf - Q) ** 2).sum(1).mean())), Pf


def matched_ca(model1, chain1, model2, chain2):
    """Cα coordinates of residues present in both structures, matched by number."""
    d1 = {r.id[1]: r["CA"].coord for r in residues(model1, chain1)}
    d2 = {r.id[1]: r["CA"].coord for r in residues(model2, chain2)}
    keys = sorted(set(d1) & set(d2))
    return (np.array([d1[k] for k in keys]), np.array([d2[k] for k in keys]),
            np.array(keys))


def rmsd(model1, chain1, model2, chain2, resid_subset=None):
    """RMSD after optimal superposition, optionally over a subset of residues."""
    P, Q, keys = matched_ca(model1, chain1, model2, chain2)
    if resid_subset is not None:
        m = np.isin(keys, list(resid_subset))
        P, Q, keys = P[m], Q[m], keys[m]
    if len(keys) < 3:
        raise ValueError("fewer than 3 residues in common — check chain IDs and numbering")
    return kabsch(P, Q)[0], len(keys)
