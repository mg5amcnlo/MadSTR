#####################################################
#                                                   #
#  Source file of the MadSTR plugin                  #
#  Use only with consent of its authors.            #
#                                                   #
#              authors: M. Zaro,                    #
#                                                   #
#                                                   #
#####################################################
"""Compatibility layer between MadSTR and the code base it is plugged into:
MG5_aMC@NLO (2.9 to 3.6) or MadGraph7.

MadGraph7 has its own version numbering (0.x), which can not be compared to
the MG5_aMC@NLO one. Its NLO code descends from MG5_aMC@NLO 3.8.2, so the rest
of the plugin sees MadGraph7 as that version (mg5amcnlo_version) and asks
is_mg7() for the MadGraph7 specific parts.
"""

import os
import re

from madgraph import MG5DIR
import madgraph.various.misc as misc

pjoin = os.path.join

# the MG5_aMC@NLO version the NLO code of MadGraph7 descends from
MG7_AS_MG5AMCNLO_VERSION = ['3', '8', '2']


class MadSTRCompatError(Exception):
    pass


def is_mg7():
    """True when running inside MadGraph7: its executable is bin/madgraph and
    it has no bin/mg5_aMC."""
    return os.path.exists(pjoin(MG5DIR, 'bin', 'madgraph')) and \
        not os.path.exists(pjoin(MG5DIR, 'bin', 'mg5_aMC'))


def mg5amcnlo_version():
    """The version of the host code in the MG5_aMC@NLO numbering, as a list of
    strings (as misc.get_pkg_info()['version'].split('.') in MG5_aMC@NLO)."""
    if is_mg7():
        return list(MG7_AS_MG5AMCNLO_VERSION)
    return misc.get_pkg_info()['version'].split('.')


def _replace_once(text, old, new, what):
    if text.count(old) != 1:
        raise MadSTRCompatError(
            "MadSTR can not patch the real-emission template of this MadGraph7 "
            "version (%s: found %d times '%s'). Please report it to the MadSTR "
            "authors." % (what, text.count(old), old.strip()))
    return text.replace(old, new)


def mg7_real_matrix_template(os_subtraction=False):
    """Build the template of the real-emission matrix elements (matrix_N.f)
    and of the on-shell subtraction terms (matrix_N_os_M.f, when
    os_subtraction is True) from the template of the MadGraph7 installation,
    so that it follows the MadGraph7 version in use (wavefunction type, flavour
    handling, real-diagram integration channels, ...).

    The MadSTR changes are the ones of realmatrix_splitorders_madstr362.inc
    with respect to the MG5_aMC@NLO template:
     - the split-order contributions go to amp_split_os (common to_amp_split_os)
       instead of amp_split;
     - MATRIX declares AMP_KEEP_WIDTH and the common block to_os_reshuf, used by
       the on-shell diagram lines that the exporter writes.
    In the on-shell subtraction terms, the weights of the real-diagram
    integration channels are not touched: they belong to the real-emission
    matrix element the term is subtracted from, which smatrix_real evaluates
    first.
    """
    path = pjoin(MG5DIR, 'madgraph', 'iolibs', 'template_files',
                 'realmatrix_splitorders_fks.inc')
    text = open(path).read()

    text = _replace_once(text,
        "      include 'orders.inc'\n",
        "      include 'orders.inc'\n"
        "      double precision amp_split_os(amp_split_size)\n"
        "      common /to_amp_split_os/amp_split_os\n",
        'amp_split_os declaration')
    text = _replace_once(text,
        "amp_split(1:amp_split_size) = 0d0",
        "amp_split_os(1:amp_split_size) = 0d0",
        'amp_split_os reset')
    text = _replace_once(text,
        "amp_split(orders_to_amp_split_pos(amp_orders)) = ans(i)",
        "amp_split_os(orders_to_amp_split_pos(amp_orders)) = ans(i)",
        'amp_split_os filling')

    # the wavefunction declaration of MATRIX
    wavefunctions = re.findall(r'^ +.*\bW\(NWAVEFUNCS\)\n', text, re.M | re.I)
    if len(wavefunctions) != 1:
        raise MadSTRCompatError(
            "MadSTR can not patch the real-emission template of this MadGraph7 "
            "version (wavefunction declaration found %d times). Please report "
            "it to the MadSTR authors." % len(wavefunctions))
    text = text.replace(wavefunctions[0], wavefunctions[0] +
        "C for the OS subtraction\n"
        "      complex*16 AMP_KEEP_WIDTH(NGRAPHS)\n"
        "      logical str_include_pdf, str_include_flux\n"
        "      integer istr\n"
        "      common /to_os_reshuf/ str_include_pdf, str_include_flux, istr\n")

    if os_subtraction:
        text = text.replace(
            "      if (get_real_diagram_weights) real_amp2(:)=0d0\n", "")
        text = re.sub(r"      IF \(NTRY\.LT\.2 \.AND\. get_real_diagram_weights\) THEN\n"
                      r".*?\n      ENDIF\n", "", text, flags=re.S)
        text = text.replace("      if (get_real_diagram_weights) then\n"
                            "%(real_amp2_lines)s\n      endif\n", "")
        left = [l for l in text.split('\n') if 'real_amp2' in l.lower()
                and 'real_matrix_weights.inc' not in l.lower()]
        if left:
            raise MadSTRCompatError(
                "MadSTR can not remove the real-diagram channel weights from the "
                "on-shell subtraction terms for this MadGraph7 version (left: %s). "
                "Please report it to the MadSTR authors." % left)
    return text
