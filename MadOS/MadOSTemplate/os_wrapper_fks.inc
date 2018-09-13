subroutine smatrix_%(suffix)s_wrapper(p, wgt_os)
C this function computes the "on-shell subtraction term", which is subtracted 
C by the full real-emission matrix-element when some on-shell aprticles appear.
implicit none
include 'nexternal.inc'
C Arguments
double precision p(0:3, nexternal)
double precision wgt_os
C internal variables
double precision p_os(0:3, nexternal)
double precision p_reord(0:3, nexternal)
double precision wgt
double precision mom_mass, mom_wdth, dau1_mass, dau2_mass, spect_mass
C for the OS subtraction
logical os_include_pdf, os_include_flux
integer iossubtr
common /to_os_reshuf/ os_include_pdf, os_include_flux, iossubtr
integer idau1, idau2, ispect
parameter (idau1 = %(idau1)d)
parameter (idau2 = %(idau2)d)
parameter (ispect = %(ispect)d)
integer mom_perm(nexternal), mom_perm_pass(nexternal), ii, i, j
data mom_perm / %(mom_perm)s /
include 'coupl.inc'
double precision ZERO
parameter (ZERO = 0d0)
logical firsttime
data firsttime /.true./
double precision pdfratio, bwratio, fluxratio
integer ibw
integer stat
integer os_countall, os_countbad
common /to_os_count/ os_countall, os_countbad
C Functions 
double precision dot

os_countall = os_countall + 1

C do nothing for diagram removal without interference (iossubtr==1) or when no subtraction is performed (iossubtr==0), amplitudes are set to zero directly in the matrix_X.f
wgt_os = 0d0
if (iossubtr.eq.1.or.iossubtr.eq.0) return

mom_wdth = %(mom_wdth)s_keep

mom_mass = %(mom_mass)s
dau1_mass = %(dau1_mass)s
dau2_mass = %(dau2_mass)s
spect_mass = %(spect_mass)s

C if daughters are heavier than the mother, nothing has to be done
if ((dau1_mass+dau2_mass).gt.(mom_mass)) return


C first reorder the momenta
do j = 1, nexternal
do i = 0, 3
p_reord(i, j) = p(i, mom_perm(j))
enddo
enddo

C consistency check to assure the reshuffling was consistent
if (firsttime) then
if (dabs(dsqrt(dot(p_reord(0, idau1), p_reord(0, idau1))) - dau1_mass) / mom_mass .gt. 1d-3) then
write(*,*) 'DAUGHTER 1 NOT ON SHELL', dsqrt(dot(p_reord(0, idau1), p_reord(0, idau1))), dau1_mass
stop
endif

if (dabs(dsqrt(dot(p_reord(0, idau2), p_reord(0, idau2))) - dau2_mass) / mom_mass .gt. 1d-3) then
write(*,*) 'DAUGHTER 2 NOT ON SHELL', dsqrt(dot(p_reord(0, idau2), p_reord(0, idau2))), dau1_mass
stop
endif

firsttime = .false.
endif

stat=0 !(not all the reshuffling techniques provide stat, so set it to 0 from scratch)

C now do the reshuffling of the momenta according to the different values of iossubtr
C 0 and 1 have already been treated above
if (iossubtr.eq.2) then
C  iossubtr==2 -> DR with interferece
  call transform_os_ident(p_reord, p_os)
elseif (iossubtr.eq.3.or.iossubtr.eq.4) then
C  iossubtr = 3 -> DS with reshuffling on initial state, standard BW
C  iossubtr = 4 -> DS with reshuffling on initial state, running BW
  call transform_os_init(p_reord, p_os, idau1, idau2, dau1_mass, dau2_mass, mom_mass)
else if (iossubtr.eq.5.or.iossubtr.eq.6) then
C  iossubtr = 5 -> DS with reshuffling on all FS particles, standard BW
C  iossubtr = 6 -> DS with reshuffling on all FS particles, running BW
call transform_os_final(p_reord, p_os, idau1, idau2, dau1_mass, dau2_mass, mom_mass, stat)
else if (iossubtr.eq.7.or.iossubtr.eq.8) then
C  iossubtr = 7 -> DS with reshuffling on spectator, standard BW
C  iossubtr = 8-> DS with reshuffling on spectator, running BW
  call transform_os_spect(p_reord, p_os, idau1, idau2, ispect, dau1_mass, dau2_mass, spect_mass, mom_mass, stat)
else
  write(*,*) 'ERROR, iossubtr not implemented', iossubtr
  stop 1
endif

C if stat != 0, the reshuffling was not possible. just exit
if (stat.ne.0) then
  os_countbad = os_countbad + 1
  return
endif

C now set iBW for the BW ratio
C ibw==0, no BW
C ibw==1, standard BW
C ibw==2, running BW

if (iossubtr.eq.2) then
  ibw=0
else if (iossubtr.eq.3.or.iossubtr.eq.5.or.iossubtr.eq.7) then
  ibw=1
else if (iossubtr.eq.4.or.iossubtr.eq.6.or.iossubtr.eq.8) then
  ibw=2
endif

C compute BW and flux/PDF ratio, and include them if needed
call get_bw_ratio(p_reord, mom_mass, mom_wdth, idau1, idau2, ibw, bwratio) 
call get_pdf_flux_ratio(p_reord, p_os, pdfratio, fluxratio)
if (.not.os_include_pdf) pdfratio = 1d0
if (.not.os_include_flux) fluxratio = 1d0

C finally call the resonant matrix element with the reshuffled momenta
call smatrix_%(suffix)s(p_os, wgt)
C and add the weigth, compensating for the reshuffling
wgt_os = wgt * pdfratio * bwratio * fluxratio

999 continue
return
end
