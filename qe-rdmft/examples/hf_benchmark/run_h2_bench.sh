#!/bin/bash
# Benchmark sd / cg / lbfgs for both alternating and joint strategies
# on H2 with rdmft_functional = 'hf'.
#
# Pre-requisites in the working directory:
#   - pw.x  (or update PWX below)
#   - H_HSCV_PBE-1.0.UPF  (copy from CPV/examples/EXX-wf-example/)
set -eu
PWX="${PWX:-./bin/pw.x}"

# First produce the KS-HF reference once.
cat > h2_ref.in <<EOF
&CONTROL
  calculation = 'scf'
  prefix      = 'h2_ref'
  outdir      = './tmp_h2_ref'
  pseudo_dir  = './'
  verbosity   = 'low'
/
&SYSTEM
  ibrav     = 1
  celldm(1) = 12.0
  nat       = 2
  ntyp      = 1
  ecutwfc   = 25.0
  ecutrho   = 100.0
  nbnd      = 4
  assume_isolated = 'martyna-tuckerman'
  input_dft = 'hf'
  ace = .true.
  exxdiv_treatment = 'gygi-baldereschi'
  x_gamma_extrapolation = .false.
/
&ELECTRONS
  conv_thr    = 1.0d-8
  mixing_beta = 0.5d0
/
ATOMIC_SPECIES
H 1.00 H_HSCV_PBE-1.0.UPF
ATOMIC_POSITIONS {angstrom}
H 0.0 0.0 0.0
H 0.0 0.0 0.74
K_POINTS Gamma
EOF
rm -rf tmp_h2_ref && mkdir tmp_h2_ref
$PWX < h2_ref.in > h2_ref.out 2>&1
KS_HF_REF=$(grep "^!!" h2_ref.out | awk '{print $5}' | tail -1)
echo "H2 KS-HF reference: $KS_HF_REF Ry"
echo
printf "%-16s %-7s %-9s %-22s %-10s %-12s\n" "strategy" "opt" "nouter" "RDMFT_ETOT" "dE_KS" "wall(s)"
printf "%-16s %-7s %-9s %-22s %-10s %-12s\n" "----------------" "-------" "---------" "----------------------" "----------" "------------"

for strategy in alternating joint; do
  for opt in sd cg lbfgs; do
    # Build INPUT file.
    cat > h2_bench.in <<EOF
&CONTROL
  calculation = 'scf'
  prefix      = 'h2'
  outdir      = './tmp'
  pseudo_dir  = './'
  verbosity   = 'low'
/
&SYSTEM
  ibrav     = 1
  celldm(1) = 12.0
  nat       = 2
  ntyp      = 1
  ecutwfc   = 25.0
  ecutrho   = 100.0
  nbnd      = 4
  assume_isolated = 'martyna-tuckerman'
  input_dft = 'hf'
  ace = .true.
  exxdiv_treatment = 'gygi-baldereschi'
  x_gamma_extrapolation = .false.
/
&ELECTRONS
  conv_thr    = 1.0d-8
  mixing_beta = 0.5d0
/
ATOMIC_SPECIES
H 1.00 H_HSCV_PBE-1.0.UPF
ATOMIC_POSITIONS {angstrom}
H 0.0 0.0 0.0
H 0.0 0.0 0.74
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
  rdmft_occ_init_perturb = 0.10d0
  rdmft_occ_ls_stepsize  = 1.0d0
  rdmft_orb_ls_stepsize  = 0.5d0
  rdmft_verbose          = 1
/
EOF
    rm -rf tmp; mkdir tmp
    t0=$(date +%s.%N)
    timeout 60 $PWX < h2_bench.in > h2_bench.out 2>&1 || true
    t1=$(date +%s.%N)
    dt=$(python3 -c "print(f'{$t1-$t0:.2f}')")
    if [ "$strategy" = "alternating" ]; then
      nouter=$(grep -c 'RDMFT outer' h2_bench.out || echo 0)
    else
      nouter=$(grep -c 'RDMFT joint outer' h2_bench.out || echo 0)
    fi
    etot=$(grep 'RDMFT_ETOTAL' h2_bench.out | awk '{print $3}' | tail -1)
    if [ -n "$etot" ]; then
      dE=$(python3 -c "print(f'{$etot - ($KS_HF_REF):+.2e}')")
    else
      etot="FAILED"
      dE="-"
    fi
    printf "%-16s %-7s %-9s %-22s %-10s %-12s\n" "$strategy" "$opt" "$nouter" "$etot" "$dE" "$dt"
  done
done
