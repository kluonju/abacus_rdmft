#!/usr/bin/env bash
# Scan H–H distance (Å) for RDMFT (Müller) total energy and write dissoc_curve.csv.
# Uses the same INPUT/KPT/orbitals/PP as this directory; generates STRU with
# Cartesian_angstrom positions so R is the true bond length in Å.
#
# Usage:
#   ./dissoc_curve.sh              # default distance list
#   R_LIST="0.5 0.74 1.0 2.0" ./dissoc_curve.sh
#
# Requires: abacus built with RDMFT (e.g. build/abacus_3p on PATH).

set -euo pipefail

EXAMPLE_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${EXAMPLE_DIR}"

export PATH="${PATH}:/home/kluo/Documents/repo/abacus-develop/build"
ABACUS="${ABACUS:-abacus_3p}"
NP="${NP:-4}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

export OMPI_ALLOW_RUN_AS_ROOT="${OMPI_ALLOW_RUN_AS_ROOT:-1}"
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="${OMPI_ALLOW_RUN_AS_ROOT_CONFIRM:-1}"
export OMPI_MCA_btl_vader_single_copy_mechanism="${OMPI_MCA_btl_vader_single_copy_mechanism:-none}"

if [[ -f /opt/intel/oneapi/setvars.sh ]]; then
	set +u
	# shellcheck source=/dev/null
	source /opt/intel/oneapi/setvars.sh
	set -u
fi

# Default: equilibrium region + dissociation (Å)
R_LIST="${R_LIST:-0.5 0.6 0.7 0.74 0.8 0.9 1.0 1.1 1.2 1.4 1.6 1.8 2.0 2.5 3.0 3.5 4.0 4.5 5.0}"

OUTCSV="${OUTCSV:-dissoc_curve.csv}"
RUN_ROOT="${RUN_ROOT:-dissoc_runs}"
mkdir -p "${RUN_ROOT}"

echo "R_Angstrom E_total_eV" >"${OUTCSV}"

write_stru() {
	local r="$1"
	local f="$2"
	cat >"${f}" <<EOF
ATOMIC_SPECIES
H 1.000 H.LDA.UPF

NUMERICAL_ORBITAL
H_gga_8au_60Ry_2s1p.orb

LATTICE_CONSTANT
1.0

LATTICE_VECTORS
10.0 0.0 0.0
0.0 10.0 0.0
0.0 0.0 10.0

ATOMIC_POSITIONS
Cartesian_angstrom

H
0.0
2
0.0 0.0 0.0
0.0 0.0 ${r}
EOF
}

for R in ${R_LIST}; do
	tag=$(echo "${R}" | tr '.' '_')
	rdir="${RUN_ROOT}/d_${tag}"
	mkdir -p "${rdir}"
	write_stru "${R}" "${rdir}/STRU"
	cp INPUT KPT "${rdir}/"
	ln -sf "${EXAMPLE_DIR}/H.LDA.UPF" "${rdir}/H.LDA.UPF"
	ln -sf "${EXAMPLE_DIR}/H_gga_8au_60Ry_2s1p.orb" "${rdir}/H_gga_8au_60Ry_2s1p.orb"

	echo "=== R = ${R} Å  (dir: ${rdir}) ==="
	(
		cd "${rdir}"
		mpirun -np "${NP}" "${ABACUS}" 2>&1 | tee log
	)
	# Final total energy (Ry / eV block)
	energy=""
	if [[ -f "${rdir}/OUT.ABACUS/running_scf.log" ]]; then
		energy="$(grep '!FINAL_ETOT_IS' "${rdir}/OUT.ABACUS/running_scf.log" | tail -1 | awk '{print $2}')"
	fi
	if [[ -z "${energy}" ]]; then
		echo "WARN: missing energy for R=${R} Å" >&2
		echo "${R} nan" >>"${OUTCSV}"
	else
		echo "${R} ${energy}" >>"${OUTCSV}"
	fi
done

echo "Wrote ${OUTCSV}"
