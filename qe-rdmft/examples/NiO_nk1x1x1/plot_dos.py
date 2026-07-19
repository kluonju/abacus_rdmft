#!/usr/bin/env python3
"""Plot RDMFT total and atom-projected DOS from *.rdmft.dos / *.rdmft.pdos_at* files."""

import argparse
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

RYTOEV = 13.605693009
L_NAMES = {0: "s", 1: "p", 2: "d", 3: "f"}
SECTION_RE = re.compile(
    r"^#\s*l=(\d+)\s+spin=(\d+)\s+atom=(\d+)\s+\(([^)]+)\)\s*$"
)


def ry_to_ev_plot(energy: np.ndarray, dos: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Convert QE RDMFT output (Ry, states/Ry/cell) to plot units (eV)."""
    return energy * RYTOEV, dos / RYTOEV


def convert_sections_to_ev(
    sections: Dict[Tuple[int, int], Tuple[np.ndarray, np.ndarray]],
) -> Dict[Tuple[int, int], Tuple[np.ndarray, np.ndarray]]:
    return {
        key: ry_to_ev_plot(energy, pdos)
        for key, (energy, pdos) in sections.items()
    }


def read_total_dos(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return energy (Ry), spin-up DOS, spin-down DOS (negative, as written)."""
    energy: List[float] = []
    up: List[float] = []
    down: List[float] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        energy.append(float(parts[0]))
        up.append(float(parts[1]))
        down.append(float(parts[2]))
    return np.asarray(energy), np.asarray(up), np.asarray(down)


def read_pdos(path: Path) -> Tuple[str, Dict[Tuple[int, int], Tuple[np.ndarray, np.ndarray]]]:
    """Parse one *.pdos_at* file into per-(l, spin) curves."""
    species = ""
    sections: Dict[Tuple[int, int], Tuple[np.ndarray, np.ndarray]] = {}
    current_key = None  # type: Optional[Tuple[int, int]]
    energy: List[float] = []
    pdos: List[float] = []

    def flush() -> None:
        nonlocal current_key, energy, pdos
        if current_key is not None and energy:
            sections[current_key] = (np.asarray(energy), np.asarray(pdos))
        energy = []
        pdos = []

    for line in path.read_text().splitlines():
        m = SECTION_RE.match(line.strip())
        if m:
            flush()
            l_val, spin, _atom, sp = int(m.group(1)), int(m.group(2)), m.group(3), m.group(4)
            species = sp
            current_key = (l_val, spin)
            continue
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        energy.append(float(parts[0]))
        pdos.append(float(parts[1]))
    flush()
    return species, sections


def discover_prefix(calc_dir: Path) -> str:
    matches = sorted(calc_dir.glob("*.rdmft.dos"))
    if not matches:
        raise FileNotFoundError(f"No *.rdmft.dos file in {calc_dir}")
    if len(matches) > 1:
        names = ", ".join(p.name for p in matches)
        raise ValueError(f"Multiple DOS files found ({names}); use --prefix")
    return matches[0].name[: -len(".dos")]


def sum_channels(
    sections: Dict[Tuple[int, int], Tuple[np.ndarray, np.ndarray]],
    l_values=None,  # type: Optional[Iterable[int]]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sum selected l channels into total spin-resolved PDOS."""
    keys = sorted(sections)
    if not keys:
        raise ValueError("empty PDOS sections")

    energy = sections[keys[0]][0]
    up = np.zeros_like(energy)
    down = np.zeros_like(energy)
    for (l_val, spin), (e, d) in sections.items():
        if l_values is not None and l_val not in l_values:
            continue
        if spin == 1:
            up += d
        else:
            down += d
    return energy, up, down


def add_fermi_line(ax, ylim: Tuple[float, float]) -> None:
    ax.axvline(0.0, color="k", linewidth=0.8, linestyle="--", alpha=0.7, zorder=0)
    ax.text(
        0.02,
        0.97,
        r"$E_F$",
        transform=ax.get_xaxis_transform(),
        ha="left",
        va="top",
        fontsize=9,
    )


def plot_spin_mirror(
    ax,
    energy: np.ndarray,
    up: np.ndarray,
    down: np.ndarray,
    label: str,
    color: str,
    alpha: float = 0.85,
    up_color: Optional[str] = None,
    down_color: Optional[str] = None,
) -> None:
    c_up = up_color or color
    c_dn = down_color or color
    ax.fill_between(energy, 0.0, up, color=c_up, alpha=0.35, linewidth=0)
    ax.fill_between(energy, down, 0.0, color=c_dn, alpha=0.35, linewidth=0)
    ax.plot(energy, up, color=c_up, linewidth=1.2, label=f"{label} ↑", alpha=alpha)
    ax.plot(
        energy, down, color=c_dn, linewidth=1.2, label=f"{label} ↓", alpha=alpha, linestyle="--"
    )


def spin_unpolarized(up: np.ndarray, down: np.ndarray) -> np.ndarray:
    """Sum spin channels into a single positive spectral density."""
    return up + np.abs(down)


def aggregate_pdos(
    atom_data: List[Tuple[str, int, Dict[Tuple[int, int], Tuple[np.ndarray, np.ndarray]]]],
    species_match,
    l_values,
) -> np.ndarray:
    """Sum PDOS over selected atoms and l channels (both spins, positive)."""
    energy = atom_data[0][2][sorted(atom_data[0][2])[0]][0]
    total = np.zeros_like(energy)
    for species, _atom_idx, sections in atom_data:
        if not species_match(species):
            continue
        for (l_val, _spin), (_e, pdos) in sections.items():
            if l_val in l_values:
                total += np.abs(pdos)
    return total


def split_ni_d(ni_d: np.ndarray, energy: np.ndarray, split_e: float) -> Tuple[np.ndarray, np.ndarray]:
    """Approximate t2g / eg from l-resolved Ni 3d using a crystal-field energy split."""
    t2g = np.where(energy < split_e, ni_d, 0.0)
    eg = np.where(energy >= split_e, ni_d, 0.0)
    return t2g, eg


def peak_in_window(
    energy: np.ndarray, curve: np.ndarray, emin: float, emax: float
) -> Tuple[float, float]:
    mask = (energy >= emin) & (energy <= emax)
    if not np.any(mask):
        return float(energy[0]), float(curve[0])
    idx = int(np.argmax(curve[mask]))
    e_win = energy[mask]
    return float(e_win[idx]), float(curve[mask][idx])


def plot_nio_style(
    energy: np.ndarray,
    total: np.ndarray,
    o_p: np.ndarray,
    t2g: np.ndarray,
    eg: np.ndarray,
    emin: float,
    emax: float,
    label: str,
    output: Path,
    scale: Optional[float] = None,
) -> None:
    """Publication-style NiO DOS: total outline plus O-p, t2g, and eg projections."""
    mask = (energy >= emin) & (energy <= emax)
    curves = [total, o_p, t2g, eg]
    peak = max(np.max(curve[mask]) for curve in curves)
    factor = scale if scale is not None else (50.0 / peak if peak > 0 else 1.0)

    total = total * factor
    o_p = o_p * factor
    t2g = t2g * factor
    eg = eg * factor

    fig, ax = plt.subplots(figsize=(4.5, 4.0))
    ax.fill_between(energy, 0.0, t2g, color="#1f4e79", alpha=0.85, linewidth=0, zorder=1)
    ax.fill_between(energy, 0.0, eg, color="#7ec8e3", alpha=0.85, linewidth=0, zorder=2)
    ax.fill_between(energy, 0.0, o_p, color="#ffb6c1", alpha=0.45, linewidth=0, zorder=3)
    ax.plot(energy, o_p, color="#c71585", linewidth=1.4, zorder=4)
    ax.plot(energy, total, color="black", linewidth=2.0, zorder=5)

    ax.axvline(0.0, color="#8b4513", linewidth=1.0, linestyle="--", zorder=0)
    ax.set_xlim(emin, emax)
    ax.set_ylim(0.0, max(total[mask].max() * 1.08, 1.0))
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("Spectral density (arb. units)")
    ax.text(0.03, 0.97, label, transform=ax.transAxes, ha="left", va="top", fontsize=12)

    t2g_e, t2g_y = peak_in_window(energy, t2g, emin, -1.5)
    eg_e, eg_y = peak_in_window(energy, eg, -1.0, emax)
    op_e, op_y = peak_in_window(energy, o_p, -3.0, 2.0)

    ax.annotate(
        "O-p",
        xy=(op_e, op_y),
        xytext=(op_e - 4.5, op_y + 0.18 * ax.get_ylim()[1]),
        fontsize=9,
        color="#c71585",
        arrowprops=dict(arrowstyle="-|>", color="#c71585", lw=0.8, shrinkA=0, shrinkB=2),
    )
    ax.annotate(
        r"$t_{2g}$",
        xy=(t2g_e, t2g_y),
        xytext=(t2g_e - 4.0, t2g_y + 0.12 * ax.get_ylim()[1]),
        fontsize=9,
        color="#1f4e79",
        arrowprops=dict(arrowstyle="-|>", color="#1f4e79", lw=0.8, shrinkA=0, shrinkB=2),
    )
    ax.annotate(
        r"$e_g$",
        xy=(eg_e, eg_y),
        xytext=(eg_e + 3.5, eg_y + 0.10 * ax.get_ylim()[1]),
        fontsize=9,
        color="#007aa3",
        arrowprops=dict(arrowstyle="-|>", color="#007aa3", lw=0.8, shrinkA=0, shrinkB=2),
    )

    ax.tick_params(direction="in", top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)

    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)


def style_dos_axis(ax, emin: float, emax: float, title: str) -> None:
    ax.axhline(0.0, color="0.5", linewidth=0.6)
    add_fermi_line(ax, ax.get_ylim())
    ax.set_xlim(emin, emax)
    ax.set_ylabel("DOS (states/eV/cell)")
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8, loc="upper right")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Calculation directory (default: script location)",
    )
    parser.add_argument(
        "--prefix",
        default=None,
        help="RDMFT DOS prefix, e.g. nio.rdmft (default: auto-detect)",
    )
    parser.add_argument("--emin", type=float, default=None, help="Minimum energy (eV rel. to mu)")
    parser.add_argument("--emax", type=float, default=None, help="Maximum energy (eV rel. to mu)")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output PNG (default: <prefix>_dos.png in calc dir)",
    )
    parser.add_argument(
        "--style",
        choices=("default", "nio"),
        default="default",
        help="Plot style: default spin-resolved panels or NiO publication style",
    )
    parser.add_argument(
        "--label",
        default="NiO",
        help="Material label for --style nio (default: NiO)",
    )
    parser.add_argument(
        "--d-split",
        type=float,
        default=-1.0,
        help="Energy (eV) split between t2g and eg Ni 3d for --style nio",
    )
    parser.add_argument(
        "--total-only",
        action="store_true",
        help="Plot total spin-resolved DOS only (no atom projections)",
    )
    parser.add_argument(
        "--by-orbital",
        action="store_true",
        help="Plot s/p/d channels instead of l-summed atom totals",
    )
    args = parser.parse_args()

    calc_dir = args.dir.resolve()
    prefix = args.prefix or discover_prefix(calc_dir)
    dos_path = calc_dir / f"{prefix}.dos"
    if not dos_path.is_file():
        raise SystemExit(f"Missing total DOS file: {dos_path}")

    if args.emin is None:
        args.emin = -15.0
    if args.emax is None:
        args.emax = 15.0

    e_tot_ry, dos_up_ry, dos_dn_ry = read_total_dos(dos_path)
    e_tot, dos_up = ry_to_ev_plot(e_tot_ry, dos_up_ry)
    _, dos_dn = ry_to_ev_plot(e_tot_ry, dos_dn_ry)

    if args.style == "nio":
        pdos_paths = sorted(calc_dir.glob(f"{prefix}.pdos_at*_*"))
        if not pdos_paths:
            raise SystemExit(f"No PDOS files matching {prefix}.pdos_at*_* in {calc_dir}")

        atom_data = []
        for path in pdos_paths:
            species, sections = read_pdos(path)
            atom_idx = int(path.name.rsplit("_", 1)[-1])
            atom_data.append((species, atom_idx, convert_sections_to_ev(sections)))

        total = spin_unpolarized(dos_up, dos_dn)
        o_p = aggregate_pdos(atom_data, lambda sp: sp == "O", {1})
        ni_d = aggregate_pdos(atom_data, lambda sp: sp.startswith("Ni"), {2})
        t2g, eg = split_ni_d(ni_d, e_tot, args.d_split)

        out = args.output or (calc_dir / f"{prefix}_nio_style.png")
        plot_nio_style(
            e_tot, total, o_p, t2g, eg, args.emin, args.emax, args.label, out
        )
        print(f"Wrote {out}")
        return

    if args.total_only:
        out = args.output or (calc_dir / f"{prefix}_tdos.png")
        fig, ax = plt.subplots(figsize=(8, 5))
        plot_spin_mirror(
            ax, e_tot, dos_up, dos_dn, "total", "black", up_color="#2166ac", down_color="#b2182b"
        )
        style_dos_axis(ax, args.emin, args.emax, f"{prefix}: total DOS")
        ax.set_xlabel("Energy (eV, relative to RDMFT $\\mu$)")
        fig.tight_layout()
        fig.savefig(out, dpi=150)
        print(f"Wrote {out}")
        return

    pdos_paths = sorted(calc_dir.glob(f"{prefix}.pdos_at*_*"))
    if not pdos_paths:
        raise SystemExit(f"No PDOS files matching {prefix}.pdos_at*_* in {calc_dir}")

    atom_data: List[Tuple[str, int, Dict[Tuple[int, int], Tuple[np.ndarray, np.ndarray]]]] = []
    for path in pdos_paths:
        species, sections = read_pdos(path)
        atom_idx = int(path.name.rsplit("_", 1)[-1])
        atom_data.append((species, atom_idx, convert_sections_to_ev(sections)))

    out = args.output or (calc_dir / f"{prefix}_dos.png")

    if args.by_orbital:
        nrows = 1 + len(atom_data)
        fig, axes = plt.subplots(nrows, 1, figsize=(8, 2.8 * nrows), sharex=True)
        if nrows == 1:
            axes = [axes]
    else:
        fig, axes = plt.subplots(
            2, 1, figsize=(8, 8), sharex=True, gridspec_kw={"height_ratios": [1.2, 2.0]}
        )

    colors = plt.cm.tab10(np.linspace(0, 1, max(len(atom_data), 3)))

    # Total DOS
    ax_tot = axes[0]
    plot_spin_mirror(ax_tot, e_tot, dos_up, dos_dn, "total", "black")
    style_dos_axis(ax_tot, args.emin, args.emax, f"{prefix}: total DOS")

    if args.by_orbital:
        for ax, (species, atom_idx, sections), color in zip(axes[1:], atom_data, colors):
            for (l_val, spin), (energy, pdos) in sorted(sections.items()):
                ls = "-" if spin == 1 else "--"
                ax.plot(
                    energy,
                    pdos,
                    color=color,
                    linewidth=1.0,
                    linestyle=ls,
                    label=f"{L_NAMES.get(l_val, str(l_val))} {'↑' if spin == 1 else '↓'}",
                )
            style_dos_axis(ax, args.emin, args.emax, f"atom {atom_idx} ({species})")
    else:
        ax_pdos = axes[1]
        for (species, atom_idx, sections), color in zip(atom_data, colors):
            energy, up, down = sum_channels(sections)
            plot_spin_mirror(ax_pdos, energy, up, down, f"{species} (at.{atom_idx})", color)
        style_dos_axis(ax_pdos, args.emin, args.emax, f"{prefix}: atom-projected DOS (l-summed)")
        ax_pdos.set_xlabel("Energy (eV, relative to RDMFT $\\mu$)")

    if args.by_orbital:
        axes[-1].set_xlabel("Energy (eV, relative to RDMFT $\\mu$)")

    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
