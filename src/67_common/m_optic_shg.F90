!!****m* ABINIT/m_optic_shg
!! NAME
!! m_optic_shg
!!
!! FUNCTION
!!  Independent SHG formulas for the optic executable.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_optic_shg

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi

 use defs_datatypes, only : ebands_t
 use m_crystal,      only : crystal_t
 use m_io_tools,     only : open_file

 implicit none

 private

 public :: optic_shg

contains

subroutine optic_shg(icomp, itemp, nband_sum, cryst, ks_ebands, pmat, v1, v2, v3, nmesh, de, sc, brod, tol, &
                     fnam, shg_formula, shg_use_wrong_scissor, comm)

 integer, intent(in) :: icomp, itemp, nband_sum, v1, v2, v3, nmesh, shg_formula, shg_use_wrong_scissor, comm
 type(crystal_t), intent(in) :: cryst
 type(ebands_t), intent(in) :: ks_ebands
 complex(dpc), intent(in) :: pmat(ks_ebands%mband, ks_ebands%mband, ks_ebands%nkpt, 3, ks_ebands%nsppol)
 real(dp), intent(in) :: de, sc, brod, tol
 character(len=*), intent(in) :: fnam

 integer, parameter :: master = 0
 integer :: fout, ifreq, ik, ierr, isp, mband, my_k1, my_k2, my_rank, nfreq
 real(dp) :: unit_factor
 real(dp), allocatable :: energy(:), occ(:), weight(:,:,:)
 complex(dpc), allocatable :: chi(:), delta(:,:,:), p(:,:,:), r(:,:,:), r_d(:,:,:,:)
 character(len=fnlen) :: outfile
 character(len=500) :: msg

 my_rank = xmpi_comm_rank(comm)
 mband = ks_ebands%mband
 ABI_CHECK(nband_sum <= mband, "nband_sum <= mband")

 if (v1 <= 0 .or. v1 > 3 .or. v2 <= 0 .or. v2 > 3 .or. v3 <= 0 .or. v3 > 3) then
   ABI_ERROR("optic_shg tensor components must be between 1 and 3")
 end if
 if (shg_formula /= 1 .and. shg_formula /= 2) then
   ABI_ERROR("shg_formula must be 1 (Rashkeev1998 dynamic) or 2 (Lin1999 static)")
 end if

 if (shg_formula == 2) then
   nfreq = 1
 else
   nfreq = nmesh
 end if
 ABI_CALLOC(chi, (nfreq))
 ABI_MALLOC(energy, (nband_sum))
 ABI_MALLOC(occ, (nband_sum))
 ABI_MALLOC(weight, (3, 3, 3))
 ABI_MALLOC(p, (3, nband_sum, nband_sum))
 ABI_MALLOC(delta, (3, nband_sum, nband_sum))
 ABI_MALLOC(r, (3, nband_sum, nband_sum))
 ABI_MALLOC(r_d, (3, 3, nband_sum, nband_sum))

 call component_projection_weights(cryst, [v1, v2, v3], weight)
 unit_factor = shg_unit_pmv(cryst%ucvol)

 call xmpi_split_work(ks_ebands%nkpt, comm, my_k1, my_k2)

 do ik = my_k1, my_k2
   do isp = 1, ks_ebands%nsppol
     energy(:) = ks_ebands%eig(1:nband_sum, ik, isp)
     where (ks_ebands%eig(1:nband_sum, ik, isp) <= ks_ebands%fermie)
       occ(:) = one
     elsewhere
       occ(:) = zero
     end where

     call copy_pmat_kpoint(pmat, ik, isp, nband_sum, p)
     if (shg_formula == 2 .and. shg_use_wrong_scissor /= 0) call apply_wrong_scissor_to_p(p, occ, energy, sc, tol)

     if (shg_formula == 1) then
       call calc_delta(p, delta)
       call calc_r(p, energy, tol, r)
       call calc_r_derivative(r, delta, energy, occ, sc, shg_use_wrong_scissor /= 0, tol, r_d)
     end if

     do ifreq = 1, nfreq
       chi(ifreq) = chi(ifreq) + formula_spin_weight(shg_formula, ks_ebands%nsppol)*ks_ebands%wtk(ik) &
                    *calc_projected_component(shg_formula, p, r, r_d, delta, occ, energy, weight, &
                    cmplx(real(ifreq - 1, dp)*de, brod, kind=dp), sc, shg_use_wrong_scissor /= 0, tol)
     end do
   end do
 end do

 call xmpi_sum(chi, comm, ierr)

 if (my_rank == master) then
   if (shg_formula == 1) then
     outfile = trim(fnam)//"-RashkeevSHG.out"
   else
     outfile = trim(fnam)//"-Lin1999StaticSHG.out"
   end if
   if (open_file(outfile, msg, newunit=fout, action='WRITE', form='FORMATTED') /= 0) then
     ABI_ERROR(msg)
   end if
   call write_shg_output(fout, chi*unit_factor, nfreq, de, sc, brod, tol, v1, v2, v3, shg_formula, &
                         shg_use_wrong_scissor, nband_sum, itemp, icomp)
   close(fout)
 end if

 ABI_FREE(chi)
 ABI_FREE(energy)
 ABI_FREE(occ)
 ABI_FREE(weight)
 ABI_FREE(p)
 ABI_FREE(delta)
 ABI_FREE(r)
 ABI_FREE(r_d)

end subroutine optic_shg

subroutine copy_pmat_kpoint(pmat, ik, isp, nband, p)
 integer, intent(in) :: ik, isp, nband
 complex(dpc), intent(in) :: pmat(:, :, :, :, :)
 complex(dpc), intent(out) :: p(3, nband, nband)
 integer :: a

 do a = 1, 3
   p(a, :, :) = pmat(1:nband, 1:nband, ik, a, isp)
 end do
end subroutine copy_pmat_kpoint

function formula_spin_weight(formula, nsppol) result(weight)
 integer, intent(in) :: formula, nsppol
 real(dp) :: weight

 weight = two
 if (formula == 2) weight = one
 if (nsppol == 2) weight = one
end function formula_spin_weight

subroutine apply_wrong_scissor_to_p(p, occ, energy, sc, etol)
 complex(dpc), intent(inout) :: p(:, :, :)
 real(dp), intent(in) :: occ(:), energy(:), sc, etol
 integer :: m, n, nb
 real(dp) :: correction, e_mn

 if (abs(sc) < etol) return
 nb = size(energy)
 do n = 1, nb
   if (occ(n) <= half) cycle
   do m = 1, nb
     if (occ(m) >= half) cycle
     e_mn = energy(m) - energy(n)
     if (abs(e_mn) <= etol) cycle
     correction = (e_mn + sc)/e_mn
     p(:, n, m) = p(:, n, m)*correction
     p(:, m, n) = p(:, m, n)*correction
   end do
 end do
end subroutine apply_wrong_scissor_to_p

function calc_projected_component(formula, p, r, r_d, delta, occ, energy, weights, omega, sc, use_wrong_scissor, etol) result(chi)
 integer, intent(in) :: formula
 complex(dpc), intent(in) :: p(:, :, :), r(:, :, :), r_d(:, :, :, :), delta(:, :, :)
 real(dp), intent(in) :: occ(:), energy(:), weights(3, 3, 3), sc, etol
 complex(dpc), intent(in) :: omega
 logical, intent(in) :: use_wrong_scissor
 complex(dpc) :: chi
 integer :: a, b, c

 chi = czero
 do a = 1, 3
   do b = 1, 3
     do c = 1, 3
       if (abs(weights(a, b, c)) <= tol12) cycle
       if (formula == 1) then
         chi = chi + weights(a, b, c)*calc_chi_rashkeev_component(a, b, c, r, r_d, delta, occ, energy, omega, sc, etol)
       else
         chi = chi + weights(a, b, c)*calc_chi_lin1999_component(a, b, c, p, occ, energy, sc, use_wrong_scissor, etol)
       end if
     end do
   end do
 end do
end function calc_projected_component

subroutine calc_delta(p, delta)
 complex(dpc), intent(in) :: p(:, :, :)
 complex(dpc), intent(out) :: delta(:, :, :)
 integer :: a, m, n, nb

 nb = size(p, 2)
 do a = 1, 3
   do n = 1, nb
     do m = 1, nb
       delta(a, n, m) = p(a, n, n) - p(a, m, m)
     end do
   end do
 end do
end subroutine calc_delta

subroutine calc_r(p, energy, etol, r)
 complex(dpc), intent(in) :: p(:, :, :)
 real(dp), intent(in) :: energy(:), etol
 complex(dpc), intent(out) :: r(:, :, :)
 integer :: a, m, n, nb
 real(dp) :: de

 nb = size(energy)
 r = czero
 do n = 1, nb
   do m = 1, nb
     de = energy(n) - energy(m)
     if (abs(de) > etol) then
       do a = 1, 3
         r(a, n, m) = p(a, n, m)/(j_dpc*de)
       end do
     end if
   end do
 end do
end subroutine calc_r

subroutine calc_r_derivative(r, delta, energy, occ, sc, use_wrong_scissor, etol, r_d)
 complex(dpc), intent(in) :: r(:, :, :), delta(:, :, :)
 real(dp), intent(in) :: energy(:), occ(:), sc, etol
 logical, intent(in) :: use_wrong_scissor
 complex(dpc), intent(out) :: r_d(:, :, :, :)
 integer :: a, b, l, m, n, nb
 real(dp) :: de, e_lm, e_nl
 real(dp), allocatable :: energy_deriv(:)
 complex(dpc) :: term1, term2

 nb = size(energy)
 ABI_MALLOC(energy_deriv, (nb))
 energy_deriv = energy
 if (use_wrong_scissor .and. abs(sc) > etol) then
   where (occ < half)
     energy_deriv = energy_deriv + sc
   end where
 end if

 r_d = czero
 do n = 1, nb
   do m = 1, nb
     de = energy_deriv(n) - energy_deriv(m)
     if (abs(de) <= etol) cycle
     do a = 1, 3
       do b = 1, 3
         term1 = (r(a, n, m)*delta(b, m, n) + r(b, n, m)*delta(a, m, n))/de
         term2 = czero
         do l = 1, nb
           e_lm = energy_deriv(l) - energy_deriv(m)
           e_nl = energy_deriv(n) - energy_deriv(l)
           term2 = term2 + j_dpc*(e_lm*r(a, n, l)*r(b, l, m) - e_nl*r(b, n, l)*r(a, l, m))/de
         end do
         r_d(a, b, n, m) = term1 + term2
       end do
     end do
   end do
 end do
 ABI_FREE(energy_deriv)
end subroutine calc_r_derivative

function calc_chi_rashkeev_component(a, b, c, r, r_d, delta, occ, energy, omega, sc, etol) result(chi)
 integer, intent(in) :: a, b, c
 complex(dpc), intent(in) :: r(:, :, :), r_d(:, :, :, :), delta(:, :, :), omega
 real(dp), intent(in) :: occ(:), energy(:), sc, etol
 complex(dpc) :: chi

 chi = calc_chi_e(a, b, c, r, occ, energy, omega, sc, etol) &
       + calc_chi_i(a, b, c, r, r_d, delta, occ, energy, omega, sc, etol)
end function calc_chi_rashkeev_component

function calc_chi_e(a, b, c, r, occ, energy, omega, sc, etol) result(chi_e)
 integer, intent(in) :: a, b, c
 complex(dpc), intent(in) :: r(:, :, :), omega
 real(dp), intent(in) :: occ(:), energy(:), sc, etol
 complex(dpc) :: chi_e, tmp
 integer :: l, m, n, nb
 real(dp) :: e_lm, e_ln, e_ml, e_mn, f_ln, f_ml, f_nm

 nb = size(energy)
 chi_e = czero
 do n = 1, nb
   do m = 1, nb
     do l = 1, nb
       f_nm = occ(n) - occ(m)
       f_ln = occ(l) - occ(n)
       f_ml = occ(m) - occ(l)
       e_ln = energy(l) - energy(n) - f_ln*sc
       e_ml = energy(m) - energy(l) - f_ml*sc
       e_mn = energy(m) - energy(n) + f_nm*sc
       if (abs(e_ln - e_ml) < etol) cycle
       tmp = czero
       if (abs(f_nm) > half .and. safe_complex_denom(e_mn - two*omega, etol)) tmp = tmp + two*f_nm/(e_mn - two*omega)
       if (abs(f_ln) > half .and. safe_complex_denom(e_ln - omega, etol)) tmp = tmp + f_ln/(e_ln - omega)
       if (abs(f_ml) > half .and. safe_complex_denom(e_ml - omega, etol)) tmp = tmp + f_ml/(e_ml - omega)
       chi_e = chi_e + sym_r3(a, b, c, n, m, l, r)*tmp/(e_ln - e_ml)
     end do
   end do
 end do
end function calc_chi_e

function calc_chi_i(a, b, c, r, r_d, delta, occ, energy, omega, sc, etol) result(chi_i)
 integer, intent(in) :: a, b, c
 complex(dpc), intent(in) :: r(:, :, :), r_d(:, :, :, :), delta(:, :, :), omega
 real(dp), intent(in) :: occ(:), energy(:), sc, etol
 complex(dpc) :: chi_i, term1, term2, term3, term4
 integer :: m, n, nb
 real(dp) :: e_mn, f_nm

 nb = size(energy)
 chi_i = czero
 do n = 1, nb
   do m = 1, nb
     f_nm = occ(n) - occ(m)
     if (abs(f_nm) <= half) cycle
     e_mn = energy(m) - energy(n) + f_nm*sc
     if (abs(e_mn) <= etol) cycle
     if (.not. safe_complex_denom(e_mn - omega, etol)) cycle
     if (.not. safe_complex_denom(e_mn - two*omega, etol)) cycle
     term1 = two*r(a, n, m)*(r_d(c, b, m, n) + r_d(b, c, m, n))/(e_mn*(e_mn - two*omega))
     term2 = (r_d(c, a, n, m)*r(b, m, n) + r_d(b, a, n, m)*r(c, m, n))/(e_mn*(e_mn - omega))
     term3 = r(a, n, m)*(r(b, m, n)*delta(c, m, n) + r(c, m, n)*delta(b, m, n)) &
             *(one/(e_mn - omega) - four/(e_mn - two*omega))/(e_mn*e_mn)
     term4 = (r_d(a, b, n, m)*r(c, m, n) + r_d(a, c, n, m)*r(b, m, n))/(two*e_mn*(e_mn - omega))
     chi_i = chi_i + (term1 + term2 + term3 - term4)*f_nm*j_dpc/two
   end do
 end do
end function calc_chi_i

function calc_chi_lin1999_component(a, b, c, p, occ, energy, sc, use_wrong_scissor, etol) result(chi)
 integer, intent(in) :: a, b, c
 complex(dpc), intent(in) :: p(:, :, :)
 real(dp), intent(in) :: occ(:), energy(:), sc, etol
 logical, intent(in) :: use_wrong_scissor
 complex(dpc) :: chi

 chi = cmplx(calc_chi_lin1999_vh(a, b, c, p, occ, energy, sc, use_wrong_scissor, etol) &
             + calc_chi_lin1999_ve(a, b, c, p, occ, energy, sc, use_wrong_scissor, etol), zero, kind=dp)
end function calc_chi_lin1999_component

function calc_chi_lin1999_vh(a, b, c, p, occ, energy, sc, use_wrong_scissor, etol) result(res)
 integer, intent(in) :: a, b, c
 complex(dpc), intent(in) :: p(:, :, :)
 real(dp), intent(in) :: occ(:), energy(:), sc, etol
 logical, intent(in) :: use_wrong_scissor
 real(dp) :: e_cv, e_cvp, e_vc, e_vpc, res, s_cv, s_vpc
 integer :: cb, nb, vb, vbp

 nb = size(energy)
 res = zero
 do vb = 1, nb
   if (occ(vb) <= half) cycle
   do vbp = 1, nb
     if (occ(vbp) <= half) cycle
     do cb = 1, nb
       if (occ(cb) >= half) cycle
       if (use_wrong_scissor) then
         e_cv = energy(cb) - energy(vb) + sc
         e_vpc = energy(vbp) - energy(cb) - sc
         e_vc = -e_cv
         e_cvp = -e_vpc
         if (abs(e_cv) < etol .or. abs(e_vpc) < etol) cycle
         res = res + apply_p_imag(a, b, c, p, vb, vbp, cb)*(one/(e_cv**3*e_vpc**2) + two/(e_vc**4*e_cvp))
       else
         e_cv = energy(cb) - energy(vb)
         e_vpc = energy(vbp) - energy(cb)
         e_vc = -e_cv
         e_cvp = -e_vpc
         s_cv = energy(cb) - energy(vb) + sc
         s_vpc = energy(vbp) - energy(cb) - sc
         if (abs(e_cv) < etol .or. abs(e_vpc) < etol .or. abs(s_cv) < etol) cycle
         res = res + apply_p_imag(a, b, c, p, vb, vbp, cb)*(one/(s_cv**2*e_cv*s_vpc*e_vpc) + two/(s_cv**2*e_vc**2*e_cvp))
       end if
     end do
   end do
 end do
end function calc_chi_lin1999_vh

function calc_chi_lin1999_ve(a, b, c, p, occ, energy, sc, use_wrong_scissor, etol) result(res)
 integer, intent(in) :: a, b, c
 complex(dpc), intent(in) :: p(:, :, :)
 real(dp), intent(in) :: occ(:), energy(:), sc, etol
 logical, intent(in) :: use_wrong_scissor
 real(dp) :: e_cp_v, e_cv, e_vc, e_vcp, res, s_cv, s_vcp
 integer :: cb, cbp, nb, vb

 nb = size(energy)
 res = zero
 do vb = 1, nb
   if (occ(vb) <= half) cycle
   do cb = 1, nb
     if (occ(cb) >= half) cycle
     do cbp = 1, nb
       if (occ(cbp) >= half) cycle
       if (use_wrong_scissor) then
         e_cv = energy(cb) - energy(vb) + sc
         e_vc = -e_cv
         e_vcp = energy(vb) - energy(cbp) - sc
         e_cp_v = -e_vcp
         if (abs(e_cv) < etol .or. abs(e_vcp) < etol) cycle
         res = res + apply_p_imag(a, b, c, p, vb, cb, cbp)*(one/(e_cv**3*e_vcp**2) + two/(e_vc**4*e_cp_v))
       else
         s_cv = energy(cb) - energy(vb) + sc
         e_cv = energy(cb) - energy(vb)
         e_vc = -e_cv
         e_vcp = energy(vb) - energy(cbp)
         e_cp_v = -e_vcp
         s_vcp = energy(vb) - energy(cbp) - sc
         if (abs(e_cv) < etol .or. abs(e_vcp) < etol .or. abs(s_cv) < etol) cycle
         res = res + apply_p_imag(a, b, c, p, vb, cb, cbp)*(one/(s_cv**2*e_cv*s_vcp*e_vcp) + two/(s_cv**2*e_vc**2*e_cp_v))
       end if
     end do
   end do
 end do
end function calc_chi_lin1999_ve

function apply_p_imag(a, b, c, p, n, m, l) result(res)
 integer, intent(in) :: a, b, c, n, m, l
 complex(dpc), intent(in) :: p(:, :, :)
 real(dp) :: res

 res = aimag(p(a, n, m)*p(b, m, l)*p(c, l, n)) &
       + aimag(p(a, n, m)*p(c, m, l)*p(b, l, n)) &
       + aimag(p(b, n, m)*p(a, m, l)*p(c, l, n)) &
       + aimag(p(b, n, m)*p(c, m, l)*p(a, l, n)) &
       + aimag(p(c, n, m)*p(a, m, l)*p(b, l, n)) &
       + aimag(p(c, n, m)*p(b, m, l)*p(a, l, n))
end function apply_p_imag

logical function safe_complex_denom(value, etol)
 complex(dpc), intent(in) :: value
 real(dp), intent(in) :: etol

 safe_complex_denom = abs(value) > etol
end function safe_complex_denom

function sym_r3(a, b, c, n, m, l, r) result(value)
 integer, intent(in) :: a, b, c, n, m, l
 complex(dpc), intent(in) :: r(:, :, :)
 complex(dpc) :: value

 value = half*r(a, n, m)*(r(b, m, l)*r(c, l, n) + r(c, m, l)*r(b, l, n))
end function sym_r3

function shg_unit_pmv(volume_bohr3) result(unit_factor)
 real(dp), intent(in) :: volume_bohr3
 real(dp) :: au2esu, unit_factor

 au2esu = (5.29177d-11*2.99792458d0*1.0d4)/(13.60569172d0*two)
 unit_factor = (one/volume_bohr3)*au2esu*four*pi/(2.99792458d0*1.0d4)*1.0d12
end function shg_unit_pmv

subroutine component_projection_weights(cryst, component, weights)
 type(crystal_t), intent(in) :: cryst
 integer, intent(in) :: component(3)
 real(dp), intent(out) :: weights(3, 3, 3)
 integer :: aa, bb, cc, isym

 weights = zero
 do isym = 1, cryst%nsym
   do aa = 1, 3
     do bb = 1, 3
       do cc = 1, 3
         weights(aa, bb, cc) = weights(aa, bb, cc) + cryst%symrel_cart(component(1), aa, isym) &
                               *cryst%symrel_cart(component(2), bb, isym)*cryst%symrel_cart(component(3), cc, isym)
       end do
     end do
   end do
 end do
 weights = weights/real(cryst%nsym, dp)
end subroutine component_projection_weights

subroutine write_shg_output(unit, chi, nfreq, de, sc, brod, etol, v1, v2, v3, formula, use_wrong_scissor, nband_sum, itemp, icomp)
 integer, intent(in) :: unit, nfreq, v1, v2, v3, formula, use_wrong_scissor, nband_sum, itemp, icomp
 complex(dpc), intent(in) :: chi(:)
 real(dp), intent(in) :: de, sc, brod, etol
 integer :: ifreq

 if (formula == 1) then
   write(unit, '(a)') '# ABINIT optic Rashkeev1998 dynamic SHG coefficients'
 else
   write(unit, '(a)') '# ABINIT optic Lin1999 static SHG coefficients with explicit Kleinman symmetry'
 end if
 write(unit, '(a)') '# Columns: energy_eV a b c Re_chi_pm_per_V Im_chi_pm_per_V'
 write(unit, '(a,1x,i0)') '# Formula:', formula
 write(unit, '(a,1x,i0)') '# Component index:', icomp
 write(unit, '(a,1x,i0)') '# Temperature index:', itemp
 write(unit, '(a,1x,es16.8)') '# Scissor operator (Ha):', sc
 write(unit, '(a,1x,es16.8)') '# Frequency broadening eta (Ha):', brod
 write(unit, '(a,1x,es16.8)') '# Energy tolerance (Ha):', etol
 write(unit, '(a,1x,i0)') '# Wrong scissor scheme:', use_wrong_scissor
 write(unit, '(a,1x,i0)') '# Number of bands included:', nband_sum
 do ifreq = 1, nfreq
   if (formula == 2) then
     write(unit, '(1x,f16.8,3(1x,i1),2(1x,es24.16))') zero, v1, v2, v3, real(chi(ifreq), dp), aimag(chi(ifreq))
   else
     write(unit, '(1x,f16.8,3(1x,i1),2(1x,es24.16))') real(ifreq - 1, dp)*de*Ha_eV, v1, v2, v3, &
       real(chi(ifreq), dp), aimag(chi(ifreq))
   end if
 end do
end subroutine write_shg_output

end module m_optic_shg
