#!/usr/bin/env bash
# ----------------------------------------------------------------
# Verify that RDMFT-HF correctly recovers the converged Hartree-
# Fock total energy from "PBE-flavoured" starting orbitals: a
# pure-PBE SCF (no exact exchange) ran for only 1-2 iterations is
# used as the orbital seed, then the RDMFT-HF solver is asked to
# converge those orbitals to the HF stationary point.
#
# Usage:
#   ./run.sh                               # uses ./pw.x
#   PWX=/path/to/pw.x ./run.sh             # custom binary
# ----------------------------------------------------------------
set -eu
PWX="${PWX:-./pw.x}"

run_one () {
    local label="$1"; local in="$2"; local ref="$3"
    rm -rf "tmp_$label" && mkdir "tmp_$label"
    timeout 600 "$PWX" < "$in" > "$label.out" 2>&1 || true
    local scf=$(grep "total energy" "$label.out" | awk '{print $5}' | tail -1)
    local rd=$(grep "RDMFT_ETOTAL" "$label.out" | awk '{print $3}' | tail -1)
    if [ -z "$rd" ]; then rd="FAIL"; fi
    if [ "$rd" != "FAIL" ]; then
        local de=$(python3 -c "print(f'{$rd - ($ref):+.2e}')")
    else
        local de="-"
    fi
    printf "%-46s %16s %16s %14s\n" "$label" "$scf" "$rd" "$de"
}

# --- KS-HF reference for H2 (12-Bohr cubic, MT) ---------------
cat > h2_hf_ref.in <<'EOF'
&CONTROL
  calculation='scf', prefix='h2r', outdir='./tmp_h2_hf_ref', pseudo_dir='./'
/
&SYSTEM
  ibrav=1, celldm(1)=12.0, nat=2, ntyp=1, ecutwfc=25.0, nbnd=4
  assume_isolated='martyna-tuckerman'
  input_dft='hf', ace=.true.
  exxdiv_treatment='gygi-baldereschi', x_gamma_extrapolation=.false.
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
EOF
rm -rf tmp_h2_hf_ref && mkdir tmp_h2_hf_ref
$PWX < h2_hf_ref.in > h2_hf_ref.out 2>&1
H2_REF=$(grep "^!!" h2_hf_ref.out | awk '{print $5}' | tail -1)

# --- KS-HF reference for Si (2x2x2 multi-k) -------------------
cat > si_hf_ref.in <<'EOF'
&CONTROL
  calculation='scf', prefix='sir', outdir='./tmp_si_hf_ref', pseudo_dir='./'
/
&SYSTEM
  ibrav=2, celldm(1)=10.20, nat=2, ntyp=1, ecutwfc=18.0, nbnd=8
  input_dft='hf', ace=.true.
  exxdiv_treatment='gygi-baldereschi', x_gamma_extrapolation=.false.
  nqx1=2, nqx2=2, nqx3=2
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
EOF
rm -rf tmp_si_hf_ref && mkdir tmp_si_hf_ref
$PWX < si_hf_ref.in > si_hf_ref.out 2>&1
SI_REF=$(grep "^!!" si_hf_ref.out | awk '{print $5}' | tail -1)

echo "KS-HF reference:  H2 = $H2_REF Ry,  Si = $SI_REF Ry"
echo
printf "%-46s %16s %16s %14s\n" \
       "case (PBE 2 SCF iters -> RDMFT-HF)" "PBE(stopped)" "RDMFT-HF" "dE_KSHF"
printf "%-46s %16s %16s %14s\n" \
       "----------------------------------------------" "----------------" \
       "----------------" "--------------"

# --- H2 pure PBE (2 iters) -> RDMFT-HF ------------------------
cat > h2_pbe_to_hf.in <<'EOF'
&CONTROL
  calculation='scf', prefix='h2p', outdir='./tmp_h2_pbe', pseudo_dir='./'
/
&SYSTEM
  ibrav=1, celldm(1)=12.0, nat=2, ntyp=1, ecutwfc=25.0, nbnd=4
  assume_isolated='martyna-tuckerman'
  input_dft='pbe', ace=.true.
  exxdiv_treatment='gygi-baldereschi', x_gamma_extrapolation=.false.
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
  rdmft_constraint='projected_gradient'
  rdmft_outer_maxiter=20, rdmft_occ_maxiter=10, rdmft_orb_maxiter=20
  rdmft_energy_tol=1.0d-9, rdmft_orb_grad_tol=1.0d-7
  rdmft_occ_grad_tol=1.0d-7, rdmft_occ_tol=1.0d-6
  rdmft_occ_init_mode='ks'
/
EOF
run_one h2_pbe_to_hf h2_pbe_to_hf.in "$H2_REF"

# --- Si pure PBE (2 iters) -> RDMFT-HF (multi-k) --------------
cat > si_pbe_to_hf.in <<'EOF'
&CONTROL
  calculation='scf', prefix='sip', outdir='./tmp_si_pbe', pseudo_dir='./'
/
&SYSTEM
  ibrav=2, celldm(1)=10.20, nat=2, ntyp=1, ecutwfc=18.0, nbnd=8
  input_dft='pbe', ace=.true.
  exxdiv_treatment='gygi-baldereschi', x_gamma_extrapolation=.false.
  nqx1=2, nqx2=2, nqx3=2
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
  rdmft_constraint='projected_gradient'
  rdmft_outer_maxiter=20, rdmft_occ_maxiter=10, rdmft_orb_maxiter=20
  rdmft_energy_tol=1.0d-9, rdmft_orb_grad_tol=1.0d-7
  rdmft_occ_grad_tol=1.0d-7, rdmft_occ_tol=1.0d-6
  rdmft_occ_init_mode='ks'
/
EOF
run_one si_pbe_to_hf si_pbe_to_hf.in "$SI_REF"
