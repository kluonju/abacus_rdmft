#!/bin/bash
# Verify RDMFT-HF on multi-k Si and on charged molecules (tot_charge != 0).
set -eu
PWX="${PWX:-./bin/pw.x}"

run_case () {
  local name="$1"
  local infile="$2"
  local label="$3"
  rm -rf "tmp_$name" && mkdir "tmp_$name"
  timeout 240 $PWX < "$infile" > "$name.out" 2>&1 || true
  local ks=$(grep "^!!" "$name.out" | awk '{print $5}' | tail -1)
  local rd=$(grep "RDMFT_ETOTAL" "$name.out" | awk '{print $3}' | tail -1)
  local ne=$(grep -E "RDMFT outer" "$name.out" | tail -1 | awk -F'Ne =' '{print $2}' | awk '{print $1}')
  if [ -z "$ks" ]; then ks="-"; fi
  if [ -z "$rd" ]; then rd="FAILED"; ne="-"; fi
  if [ "$rd" != "FAILED" ] && [ -n "$ks" ] && [ "$ks" != "-" ]; then
    local de=$(python3 -c "print(f'{$rd - ($ks):+.2e}')")
  else
    local de="-"
  fi
  printf "%-40s %16s %16s %10s %14s\n" "$label" "$ks" "$rd" "$ne" "$de"
}

printf "%-40s %16s %16s %10s %14s\n" "case" "KS-HF" "RDMFT-HF" "Ne" "dE_KS"
printf "%-40s %16s %16s %10s %14s\n" "----------------------------------------" "----------------" "----------------" "----------" "--------------"

# ---- Multi-k Si: occ-only alternating ----
cat > si_alt_occ.in <<'EOF'
&CONTROL
  calculation='scf', prefix='si', outdir='./tmp_si_occ', pseudo_dir='./'
/
&SYSTEM
  ibrav=2, celldm(1)=10.20, nat=2, ntyp=1, ecutwfc=18.0, nbnd=8
  input_dft='hf', ace=.true., exxdiv_treatment='gygi-baldereschi'
  x_gamma_extrapolation=.false., nqx1=2, nqx2=2, nqx3=2
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
  rdmft_constraint='projected_gradient'
  rdmft_outer_maxiter=4, rdmft_occ_maxiter=10, rdmft_orb_maxiter=0
  rdmft_energy_tol=1.0d-8, rdmft_occ_grad_tol=1.0d-6
  rdmft_occ_tol=1.0d-6, rdmft_occ_init_mode='perturbed'
  rdmft_occ_init_perturb=0.05d0
/
EOF
run_case si_alt_occ si_alt_occ.in "Multi-k Si (2x2x2), alt, occ-only"

# ---- Closed-shell H4^2+ (tot_charge=+2) ----
run_case h4_2plus h4_charged.in "H4^2+ (tot_charge=+2), nspin=1, alt"

# ---- Open-shell H2^- (tot_charge=-1, nspin=2) ----
sed -i "s/tot_charge = +1.0/tot_charge = -1.0/; s/tot_charge = 1.0/tot_charge = -1.0/" h2_charged.in 2>/dev/null || true
cat > h2_anion.in <<'EOF'
&CONTROL
  calculation='scf', prefix='h2a', outdir='./tmp_h2a', pseudo_dir='./'
/
&SYSTEM
  ibrav=1, celldm(1)=12.0, nat=2, ntyp=1, ecutwfc=25.0, nbnd=4
  tot_charge=-1.0, nspin=2, tot_magnetization=1.0
  assume_isolated='martyna-tuckerman'
  input_dft='hf', ace=.true., exxdiv_treatment='gygi-baldereschi'
  x_gamma_extrapolation=.false.
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
  rdmft_constraint='projected_gradient'
  rdmft_outer_maxiter=4, rdmft_occ_maxiter=10, rdmft_orb_maxiter=0
  rdmft_energy_tol=1.0d-9, rdmft_occ_grad_tol=1.0d-7
  rdmft_occ_tol=1.0d-6
/
EOF
run_case h2_anion h2_anion.in "H2^- (tot_charge=-1), nspin=2, alt"

# ---- Multi-k Si: full alternating with orb_maxiter=3 (joint multi-k Stiefel descent) ----
sed 's/rdmft_orb_maxiter=0/rdmft_orb_maxiter=3/' si_alt_occ.in > si_alt_orb_joint.in
sed -i "s|outdir='./tmp_si_occ'|outdir='./tmp_si_orb_joint'|" si_alt_orb_joint.in
run_case si_alt_orb_joint si_alt_orb_joint.in "Multi-k Si (2x2x2), alt, orb_strategy=joint"

# ---- Multi-k Si: block-coordinate per-k Stiefel descent ----
sed "s/rdmft_orb_maxiter=0/rdmft_orb_maxiter=3, rdmft_orb_strategy='block_k'/" si_alt_occ.in > si_alt_orb_blockk.in
sed -i "s|outdir='./tmp_si_occ'|outdir='./tmp_si_orb_blockk'|" si_alt_orb_blockk.in
run_case si_alt_orb_blockk si_alt_orb_blockk.in "Multi-k Si (2x2x2), alt, orb_strategy=block_k"
