#!/usr/bin/env bash
# ----------------------------------------------------------------
# RDMFT-HF on spin-polarised (nspin = 2) inputs.
#
# Verifies that the RDMFT-HF post-SCF stationary point coincides
# with the converged KS-HF total energy on:
#
#   * H2 closed shell, nspin=2, tot_magnetization=0
#   * Si 2x2x2 multi-k closed shell, nspin=2, tot_magnetization=0
#   * H2^- (3 electrons), nspin=2, tot_magnetization=1
#   * H2 (closed shell) starting from a rough PBE SCF, nspin=2
#   * Si 2x2x2 starting from a rough PBE SCF, nspin=2
#
# Usage:
#   ./run.sh                               # uses ./pw.x
#   PWX=/path/to/pw.x ./run.sh             # custom binary
# ----------------------------------------------------------------
set -eu
PWX="${PWX:-./pw.x}"

run_case () {
  local name="$1"
  local infile="$2"
  local label="$3"
  rm -rf "tmp_$name" && mkdir "tmp_$name"
  timeout 600 "$PWX" < "$infile" > "$name.out" 2>&1 || true
  local ks=$(grep "^!!" "$name.out" | awk '{print $5}' | tail -1)
  if [ -z "$ks" ]; then
     # Some inputs (rough PBE) never reach EXX self-consistency;
     # fall back to the last printed total energy.
     ks=$(grep "^!" "$name.out" | awk '{print $5}' | tail -1)
  fi
  local rd=$(grep "RDMFT_ETOTAL" "$name.out" | awk '{print $3}' | tail -1)
  local mag=$(grep "total magnetization" "$name.out" | awk '{print $4}' | tail -1)
  if [ -z "$ks" ];  then ks="-";        fi
  if [ -z "$rd" ];  then rd="FAILED";   fi
  if [ -z "$mag" ]; then mag="-";       fi
  if [ "$rd" != "FAILED" ] && [ "$ks" != "-" ]; then
    local de=$(python3 -c "print(f'{$rd - ($ks):+.2e}')")
  else
    local de="-"
  fi
  printf "%-44s %16s %16s %8s %14s\n" "$label" "$ks" "$rd" "$mag" "$de"
}

printf "%-44s %16s %16s %8s %14s\n" "case" "KS-HF / KS-PBE" "RDMFT-HF" "Mz" "dE_KS"
printf "%-44s %16s %16s %8s %14s\n" "--------------------------------------------" \
       "----------------" "----------------" "--------" "--------------"

# ---- 1) H2 closed shell, nspin=2, HF SCF -> RDMFT-HF ---------
cat > h2_n2.in <<'EOF'
&CONTROL
  calculation='scf', prefix='h2n2', outdir='./tmp_h2_n2', pseudo_dir='./'
/
&SYSTEM
  ibrav=1, celldm(1)=12.0, nat=2, ntyp=1, ecutwfc=25.0, nbnd=4
  assume_isolated='martyna-tuckerman'
  input_dft='hf', ace=.true.
  exxdiv_treatment='gygi-baldereschi', x_gamma_extrapolation=.false.
  nspin=2, tot_magnetization=0
  occupations='smearing', smearing='gaussian', degauss=0.0001
/
&ELECTRONS
  conv_thr=1.0d-8, mixing_beta=0.5d0
/
ATOMIC_SPECIES
H 1.00 H_HSCV_PBE-1.0.UPF
ATOMIC_POSITIONS {angstrom}
H 0.0 0.0 0.0
H 0.0 0.0 0.74
K_POINTS Gamma
&rdmft
  do_rdmft=.true., rdmft_functional='hf'
  rdmft_solver_strategy='alternating'
  rdmft_orb_optimizer='cg', rdmft_orb_strategy='joint'
  rdmft_outer_maxiter=10, rdmft_occ_maxiter=10, rdmft_orb_maxiter=15
  rdmft_energy_tol=1.0d-9, rdmft_occ_grad_tol=1.0d-7
  rdmft_orb_grad_tol=1.0d-7, rdmft_occ_tol=1.0d-6
  rdmft_occ_init_mode='ks'
/
EOF
run_case h2_n2_hf h2_n2.in "H2 nspin=2 Mz=0, HF -> RDMFT-HF"

# ---- 2) Si multi-k closed shell, nspin=2, HF SCF -> RDMFT-HF
cat > si_n2.in <<'EOF'
&CONTROL
  calculation='scf', prefix='sin2', outdir='./tmp_si_n2', pseudo_dir='./'
/
&SYSTEM
  ibrav=2, celldm(1)=10.20, nat=2, ntyp=1, ecutwfc=18.0, nbnd=8
  input_dft='hf', ace=.true.
  exxdiv_treatment='gygi-baldereschi', x_gamma_extrapolation=.false.
  nqx1=2, nqx2=2, nqx3=2
  nspin=2, tot_magnetization=0
  occupations='smearing', smearing='gaussian', degauss=0.001
/
&ELECTRONS
  conv_thr=1.0d-7, mixing_beta=0.5d0
/
ATOMIC_SPECIES
Si 28.086 Si.upf
ATOMIC_POSITIONS {alat}
Si 0.00 0.00 0.00
Si 0.25 0.25 0.25
K_POINTS automatic
2 2 2 0 0 0
&rdmft
  do_rdmft=.true., rdmft_functional='hf'
  rdmft_solver_strategy='alternating'
  rdmft_orb_optimizer='cg', rdmft_orb_strategy='joint'
  rdmft_outer_maxiter=10, rdmft_occ_maxiter=10, rdmft_orb_maxiter=15
  rdmft_energy_tol=1.0d-9, rdmft_occ_grad_tol=1.0d-7
  rdmft_orb_grad_tol=1.0d-7, rdmft_occ_tol=1.0d-6
  rdmft_occ_init_mode='ks'
/
EOF
run_case si_n2_hf si_n2.in "Si 2x2x2 nspin=2 Mz=0, HF -> RDMFT-HF"

# ---- 3) H2^- (already in multik_charged but re-run here) -----
cat > h2m_n2.in <<'EOF'
&CONTROL
  calculation='scf', prefix='h2m', outdir='./tmp_h2m_n2', pseudo_dir='./'
/
&SYSTEM
  ibrav=1, celldm(1)=14.0, nat=2, ntyp=1, ecutwfc=25.0, nbnd=4
  assume_isolated='martyna-tuckerman'
  tot_charge=-1, nspin=2, tot_magnetization=1
  input_dft='hf', ace=.true.
  exxdiv_treatment='gygi-baldereschi', x_gamma_extrapolation=.false.
  occupations='smearing', smearing='gaussian', degauss=0.001
/
&ELECTRONS
  conv_thr=1.0d-7, mixing_beta=0.4d0
/
ATOMIC_SPECIES
H 1.00 H_HSCV_PBE-1.0.UPF
ATOMIC_POSITIONS {angstrom}
H 0.0 0.0 0.0
H 0.0 0.0 0.74
K_POINTS Gamma
&rdmft
  do_rdmft=.true., rdmft_functional='hf'
  rdmft_solver_strategy='alternating'
  rdmft_orb_optimizer='cg', rdmft_orb_strategy='joint'
  rdmft_outer_maxiter=10, rdmft_occ_maxiter=10, rdmft_orb_maxiter=15
  rdmft_energy_tol=1.0d-9, rdmft_occ_grad_tol=1.0d-7
  rdmft_orb_grad_tol=1.0d-7, rdmft_occ_tol=1.0d-6
  rdmft_occ_init_mode='ks'
/
EOF
run_case h2m_n2_hf h2m_n2.in "H2^- nspin=2 Mz=1, HF -> RDMFT-HF"

# ---- 4) H2 closed shell, nspin=2, rough PBE -> RDMFT-HF ------
cat > h2_pbe_n2.in <<'EOF'
&CONTROL
  calculation='scf', prefix='h2pn2', outdir='./tmp_h2_pbe_n2', pseudo_dir='./'
/
&SYSTEM
  ibrav=1, celldm(1)=12.0, nat=2, ntyp=1, ecutwfc=25.0, nbnd=4
  assume_isolated='martyna-tuckerman'
  input_dft='pbe', ace=.true.
  exxdiv_treatment='gygi-baldereschi', x_gamma_extrapolation=.false.
  nspin=2, tot_magnetization=0
  occupations='smearing', smearing='gaussian', degauss=0.0001
/
&ELECTRONS
  conv_thr=1.0d-2, electron_maxstep=2, mixing_beta=0.7d0
/
ATOMIC_SPECIES
H 1.00 H_HSCV_PBE-1.0.UPF
ATOMIC_POSITIONS {angstrom}
H 0.0 0.0 0.0
H 0.0 0.0 0.74
K_POINTS Gamma
&rdmft
  do_rdmft=.true., rdmft_functional='hf'
  rdmft_solver_strategy='alternating'
  rdmft_orb_optimizer='cg', rdmft_orb_strategy='joint'
  rdmft_outer_maxiter=20, rdmft_occ_maxiter=10, rdmft_orb_maxiter=20
  rdmft_energy_tol=1.0d-9, rdmft_occ_grad_tol=1.0d-7
  rdmft_orb_grad_tol=1.0d-7, rdmft_occ_tol=1.0d-6
  rdmft_occ_init_mode='ks'
/
EOF
# For the rough-PBE rows we need a *separate* HF reference: the
# table's "KS-HF / KS-PBE" column would otherwise show the
# converged-PBE energy (which is what the script picks off the
# "!" line of the SCF), and dE would then be measured against PBE
# rather than HF.  We re-use case 1 / case 2's converged HF
# energies as the proper reference and pass them as a manual
# override via the env var KS_REF_OVERRIDE.
H2_KS_HF=$(grep "^!!" h2_n2_hf.out | awk '{print $5}' | tail -1)
SI_KS_HF=$(grep "^!!" si_n2_hf.out | awk '{print $5}' | tail -1)

run_case_vs_ref () {
  local name="$1"; local infile="$2"; local label="$3"; local ref="$4"
  rm -rf "tmp_$name" && mkdir "tmp_$name"
  timeout 600 "$PWX" < "$infile" > "$name.out" 2>&1 || true
  local rd=$(grep "RDMFT_ETOTAL" "$name.out" | awk '{print $3}' | tail -1)
  local mag=$(grep "total magnetization" "$name.out" | awk '{print $4}' | tail -1)
  if [ -z "$rd" ];  then rd="FAILED";   fi
  if [ -z "$mag" ]; then mag="-";       fi
  if [ "$rd" != "FAILED" ]; then
    local de=$(python3 -c "print(f'{$rd - ($ref):+.2e}')")
  else
    local de="-"
  fi
  printf "%-44s %16s %16s %8s %14s\n" "$label" "$ref" "$rd" "$mag" "$de"
}
run_case_vs_ref h2_pbe_n2_to_hf h2_pbe_n2.in \
                "H2 nspin=2 Mz=0, PBE(2it) -> RDMFT-HF" "$H2_KS_HF"

# ---- 5) Si 2x2x2, nspin=2, rough PBE -> RDMFT-HF -------------
cat > si_pbe_n2.in <<'EOF'
&CONTROL
  calculation='scf', prefix='sipn2', outdir='./tmp_si_pbe_n2', pseudo_dir='./'
/
&SYSTEM
  ibrav=2, celldm(1)=10.20, nat=2, ntyp=1, ecutwfc=18.0, nbnd=8
  input_dft='pbe', ace=.true.
  exxdiv_treatment='gygi-baldereschi', x_gamma_extrapolation=.false.
  nqx1=2, nqx2=2, nqx3=2
  nspin=2, tot_magnetization=0
  occupations='smearing', smearing='gaussian', degauss=0.001
/
&ELECTRONS
  conv_thr=1.0d-2, electron_maxstep=2, mixing_beta=0.5d0
/
ATOMIC_SPECIES
Si 28.086 Si.upf
ATOMIC_POSITIONS {alat}
Si 0.00 0.00 0.00
Si 0.25 0.25 0.25
K_POINTS automatic
2 2 2 0 0 0
&rdmft
  do_rdmft=.true., rdmft_functional='hf'
  rdmft_solver_strategy='alternating'
  rdmft_orb_optimizer='cg', rdmft_orb_strategy='joint'
  rdmft_outer_maxiter=20, rdmft_occ_maxiter=10, rdmft_orb_maxiter=20
  rdmft_energy_tol=1.0d-9, rdmft_occ_grad_tol=1.0d-7
  rdmft_orb_grad_tol=1.0d-7, rdmft_occ_tol=1.0d-6
  rdmft_occ_init_mode='ks'
/
EOF
run_case_vs_ref si_pbe_n2_to_hf si_pbe_n2.in \
                "Si 2x2x2 nspin=2 Mz=0, PBE(2it) -> RDMFT-HF" "$SI_KS_HF"
