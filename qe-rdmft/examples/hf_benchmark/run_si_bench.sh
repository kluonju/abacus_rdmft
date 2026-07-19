#!/bin/bash
# Benchmark sd / cg / lbfgs for both alternating and joint strategies
# on Si bulk (gamma-only) with rdmft_functional = 'hf'.
#
# Pre-requisites in the working directory:
#   - pw.x  (or update PWX below)
#   - Si.upf (e.g. copy from QEHeat/examples/pseudo/Si_ONCV_PBE-1.1.upf)
set -eu
PWX="${PWX:-./bin/pw.x}"

# First produce the KS-HF reference once.
cat > si_ref.in <<EOF
&CONTROL
  calculation = 'scf'
  prefix      = 'si'
  outdir      = './tmp_si_ref'
  pseudo_dir  = './'
  verbosity   = 'low'
/
&SYSTEM
  ibrav     = 2
  celldm(1) = 10.20
  nat       = 2
  ntyp      = 1
  ecutwfc   = 18.0
  nbnd      = 8
  input_dft = 'hf'
  ace = .true.
  exxdiv_treatment = 'gygi-baldereschi'
  x_gamma_extrapolation = .false.
  nqx1 = 1, nqx2 = 1, nqx3 = 1
/
&ELECTRONS
  conv_thr    = 1.0d-7
  mixing_beta = 0.5d0
/
ATOMIC_SPECIES
Si 28.086 Si.upf
ATOMIC_POSITIONS {alat}
Si 0.00 0.00 0.00
Si 0.25 0.25 0.25
K_POINTS Gamma
EOF
rm -rf tmp_si_ref; mkdir tmp_si_ref
$PWX < si_ref.in > si_ref.out 2>&1
KS_HF_REF=$(grep "^!!" si_ref.out | awk '{print $5}' | tail -1)
if [ -z "$KS_HF_REF" ]; then
  KS_HF_REF=$(grep "^!  " si_ref.out | grep "total energy" | awk '{print $5}' | tail -1)
fi
echo "KS-HF reference: $KS_HF_REF Ry"
echo
printf "%-16s %-7s %-9s %-22s %-10s %-12s\n" "strategy" "opt" "nouter" "RDMFT_ETOT" "dE_KS" "wall(s)"
printf "%-16s %-7s %-9s %-22s %-10s %-12s\n" "----------------" "-------" "---------" "----------------------" "----------" "------------"

for strategy in alternating joint; do
  for opt in sd cg lbfgs; do
    cat > si_bench.in <<EOF
&CONTROL
  calculation = 'scf'
  prefix      = 'si'
  outdir      = './tmp'
  pseudo_dir  = './'
  verbosity   = 'low'
/
&SYSTEM
  ibrav     = 2
  celldm(1) = 10.20
  nat       = 2
  ntyp      = 1
  ecutwfc   = 18.0
  nbnd      = 8
  input_dft = 'hf'
  ace = .true.
  exxdiv_treatment = 'gygi-baldereschi'
  x_gamma_extrapolation = .false.
  nqx1 = 1, nqx2 = 1, nqx3 = 1
/
&ELECTRONS
  conv_thr    = 1.0d-7
  mixing_beta = 0.5d0
/
ATOMIC_SPECIES
Si 28.086 Si.upf
ATOMIC_POSITIONS {alat}
Si 0.00 0.00 0.00
Si 0.25 0.25 0.25
K_POINTS Gamma

&rdmft
  do_rdmft               = .true.
  rdmft_functional       = 'hf'
  rdmft_solver_strategy  = '$strategy'
  rdmft_occ_optimizer    = '$opt'
  rdmft_orb_optimizer    = '$opt'
  rdmft_joint_optimizer  = '$opt'
  rdmft_constraint       = 'projected_gradient'
  rdmft_outer_maxiter    = 30
  rdmft_occ_maxiter      = 20
  rdmft_orb_maxiter      = 10
  rdmft_energy_tol       = 1.0d-9
  rdmft_occ_grad_tol     = 1.0d-7
  rdmft_orb_grad_tol     = 1.0d-7
  rdmft_occ_tol          = 1.0d-6
  rdmft_lbfgs_memory     = 8
  rdmft_occ_init_mode    = 'perturbed'
  rdmft_occ_init_perturb = 0.05d0
  rdmft_occ_ls_stepsize  = 0.5d0
  rdmft_orb_ls_stepsize  = 0.3d0
  rdmft_verbose          = 1
/
EOF
    rm -rf tmp; mkdir tmp
    t0=$(date +%s.%N)
    timeout 90 $PWX < si_bench.in > si_bench.out 2>&1 || true
    t1=$(date +%s.%N)
    dt=$(python3 -c "print(f'{$t1-$t0:.2f}')")
    if [ "$strategy" = "alternating" ]; then
      nouter=$(grep -c 'RDMFT outer' si_bench.out || echo 0)
    else
      nouter=$(grep -c 'RDMFT joint outer' si_bench.out || echo 0)
    fi
    etot=$(grep 'RDMFT_ETOTAL' si_bench.out | awk '{print $3}' | tail -1)
    if [ -n "$etot" ]; then
      dE=$(python3 -c "print(f'{$etot - ($KS_HF_REF):+.2e}')")
    else
      etot="FAILED"
      dE="-"
    fi
    printf "%-16s %-7s %-9s %-22s %-10s %-12s\n" "$strategy" "$opt" "$nouter" "$etot" "$dE" "$dt"
  done
done
