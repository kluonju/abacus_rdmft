#!/usr/bin/env python3
"""Compare ELK TDOS.OUT (Ha) with QE nio.rdmft.dos (Ry)."""

import argparse
import sys
from pathlib import Path

RY_TO_HA = 0.5
HA_TO_RY = 2.0


def read_elk_tdos(path):
    """Return energies (Ha) and spin columns from TDOS.OUT."""
    blocks = []
    current = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                if current:
                    blocks.append(current)
                    current = []
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            e, d = float(parts[0]), float(parts[1])
            current.append((e, d))
    if current:
        blocks.append(current)
    if not blocks:
        raise ValueError(f"no data in {path}")
    energies = [e for e, _ in blocks[0]]
    cols = [[d for _, d in block] for block in blocks]
    return energies, cols


def read_qe_dos(path):
    """Return energies (Ry) and spin columns from *.rdmft.dos."""
    energies = []
    cols = []
    nspin = 1
    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                if "DOS_up" in line and "DOS_down" in line:
                    nspin = 2
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            energies.append(float(parts[0]))
            if nspin == 2 and len(parts) >= 3:
                if not cols:
                    cols = [[], []]
                cols[0].append(float(parts[1]))
                cols[1].append(-float(parts[2]))
            else:
                if not cols:
                    cols = [[]]
                cols[0].append(float(parts[1]))
    if not energies:
        raise ValueError(f"no data in {path}")
    return energies, cols


def interp_linear(x, xp, fp):
    if x <= xp[0]:
        return fp[0]
    if x >= xp[-1]:
        return fp[-1]
    for i in range(len(xp) - 1):
        if xp[i] <= x <= xp[i + 1]:
            t = (x - xp[i]) / (xp[i + 1] - xp[i])
            return fp[i] * (1.0 - t) + fp[i + 1] * t
    return fp[-1]


def compare(elk_e_ha, elk_cols, qe_e_ry, qe_cols):
    qe_e_ha = [e * RY_TO_HA for e in qe_e_ry]
    nspin = min(len(elk_cols), len(qe_cols))
    max_abs = 0.0
    l2 = 0.0
    npts = 0
    for isp in range(nspin):
        for e_ha, qe_val in zip(qe_e_ha, qe_cols[isp]):
            elk_val = interp_linear(e_ha, elk_e_ha, elk_cols[isp])
            # QE: states/Ry/cell; ELK: states/Ha/cell
            qe_val_ha = qe_val / HA_TO_RY
            diff = abs(elk_val - qe_val_ha)
            max_abs = max(max_abs, diff)
            l2 += diff * diff
            npts += 1
    rmse = (l2 / max(npts, 1)) ** 0.5
    return {"max_abs": max_abs, "rmse": rmse, "npts": float(npts)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--elk", type=Path, required=True, help="ELK TDOS.OUT path")
    parser.add_argument("--qe", type=Path, required=True, help="QE *.rdmft.dos path")
    args = parser.parse_args()

    elk_e, elk_cols = read_elk_tdos(args.elk)
    qe_e, qe_cols = read_qe_dos(args.qe)
    stats = compare(elk_e, elk_cols, qe_e, qe_cols)

    print(f"Compared {int(stats['npts'])} points on {len(qe_cols)} spin channel(s)")
    print(f"Max |ΔDOS| (states/Ha/cell): {stats['max_abs']:.6e}")
    print(f"RMSE (states/Ha/cell):      {stats['rmse']:.6e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
