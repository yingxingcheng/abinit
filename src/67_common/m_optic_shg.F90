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
 use nlokit_shg,     only : nlokit_shg_apply_wrong_scissor_to_p, nlokit_shg_calc_delta, &
                            nlokit_shg_calc_projected_component, nlokit_shg_calc_r, &
                            nlokit_shg_calc_r_derivative, nlokit_shg_formula_spin_weight

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
 real(dp), allocatable :: energy(:), energy_deriv(:), occ(:), weight(:,:,:)
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
 ABI_MALLOC(energy_deriv, (nband_sum))
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
     if (shg_formula == 2 .and. shg_use_wrong_scissor /= 0) call nlokit_shg_apply_wrong_scissor_to_p(p, occ, energy, sc, tol, tol)

     if (shg_formula == 1) then
       energy_deriv(:) = energy(:)
       if (shg_use_wrong_scissor /= 0 .and. abs(sc) > tol) then
         where (occ(:) < tol)
           energy_deriv(:) = energy_deriv(:) + sc
         end where
       end if
       call nlokit_shg_calc_delta(p, delta)
       call nlokit_shg_calc_r(p, energy, tol, r)
       call nlokit_shg_calc_r_derivative(r, delta, energy_deriv, tol, r_d)
     end if

     do ifreq = 1, nfreq
       chi(ifreq) = chi(ifreq) + nlokit_shg_formula_spin_weight(shg_formula, ks_ebands%nsppol)*ks_ebands%wtk(ik) &
                    *nlokit_shg_calc_projected_component(shg_formula, p, r, r_d, delta, occ, energy, weight, &
                    cmplx(real(ifreq - 1, dp)*de, brod, kind=dp), sc, shg_use_wrong_scissor /= 0, tol, tol)
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
 ABI_FREE(energy_deriv)
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
