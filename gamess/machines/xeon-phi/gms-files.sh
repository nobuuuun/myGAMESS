#!/bin/sh
#  file assignments.
#
#  All binary files should be put on a node's local disk ($SCR directory),
#  for the highest speed access possible.  These .Fxx files are typically
#  not saved for the next run, but they may be big and/or I/O intensive.
#
#  It is convenient to write ASCII output files (PUNCH, RESTART, TRAJECT,
#  and MAKEFP) to the user's permanent disk, on your file server.  They
#  are small, written only by the master process, and are useful outputs
#  for further runs.
#
#                        ASCII input files
#             You must edit a, but will probably skip b+c.
#  Some data files may be read by a run, each is read only once, so
#  that storage of one (1) copy on your file server is appropriate.
#  a) AUXDATA is directory of data #-- sets containing
#     Note: change only AUXDATA, not ERICFMT,MCPPATH,BASPATH,QUANPOL!
#        1. a file of Fm(t) data for ERI computations
#        2. a BASES subdirectory, containing files of some basis #-- sets
#        3. a MCP subdirectory, containing files of MCP bases and potentials
#        4. data #-- sets for the Quantum chemistry Polarizable force field
#        5. a EFP subdirectory, containing standard EFP2 potentials
#  b) The EXTBAS file contains any user-supplied basis #-- sets.
#  c) The NUCBAS or POSBAS files are nuclear or positron basis #-- sets,
#     used by the NEO method.  See NEO's documentation for more details.
#  d) there are 3 places where you might want to uncomment a '#-- set echo',
#     to see all the file definitions (one is just below).
#
if [ "$VERBOSE" = true ] ; then
  set -x
fi
export AUXDATA=$GMSPATH/auxdata
export  EXTBAS=/dev/null
export  NUCBAS=/dev/null
export  POSBAS=/dev/null
#
export ERICFMT=$AUXDATA/ericfmt.dat
export MCPPATH=$AUXDATA/MCP
export BASPATH=$AUXDATA/BASES
export QUANPOL=$AUXDATA/QUANPOL
export  MAKEFP=$USERSCR/$JOB.efp
export   GAMMA=$USERSCR/$JOB.gamma
export TRAJECT=$USERSCR/$JOB.trj
export RESTART=$USERSCR/$JOB.rst
export  QMWAVE $USERSCR/$JOB.qmw
export   MDDIP $USERSCR/$JOB.dip
export OPTHES1 $USERSCR/$JOB.hs1
export OPTHES2 $USERSCR/$JOB.hs2
export     XYZ $SCR/$JOB.xyz
export   INPUT=$SCR/$JOB.F05
export   PUNCH=$USERSCR/$JOB.dat
export  AOINTS=$SCR/$JOB.F08
export  MOINTS=$SCR/$JOB.F09
export DICTNRY=$SCR/$JOB.F10
export DRTFILE=$SCR/$JOB.F11
export CIVECTR=$SCR/$JOB.F12
export CASINTS=$SCR/$JOB.F13
export  CIINTS=$SCR/$JOB.F14
export  WORK15=$SCR/$JOB.F15
export  WORK16=$SCR/$JOB.F16
export CSFSAVE=$SCR/$JOB.F17
export FOCKDER=$SCR/$JOB.F18
export  WORK19=$SCR/$JOB.F19
export  DASORT=$SCR/$JOB.F20
export DIABDAT=$SCR/$JOB.F21
export DFTINTS=$SCR/$JOB.F21
export DFTGRID=$SCR/$JOB.F22
export  JKFILE=$SCR/$JOB.F23
export  ORDINT=$SCR/$JOB.F24
export  EFPIND=$SCR/$JOB.F25
export PCMDATA=$SCR/$JOB.F26
export PCMINTS=$SCR/$JOB.F27
export SVPWRK1=$SCR/$JOB.F26
export SVPWRK2=$SCR/$JOB.F27
export COSCAV=$SCR/$JOB.F26
export COSDATA=$USERSCR/$JOB.cosmo
export COSPOT=$USERSCR/$JOB.pot
export  MLTPL=$SCR/$JOB.F28
export  MLTPLT=$SCR/$JOB.F29
export  DAFL30=$SCR/$JOB.F30
export  SOINTX=$SCR/$JOB.F31
export  SOINTY=$SCR/$JOB.F32
export  SOINTZ=$SCR/$JOB.F33
export  SORESC=$SCR/$JOB.F34
#   35 is used by RESTART, see above
export GCILIST=$SCR/$JOB.F37
export HESSIAN=$SCR/$JOB.F38
export QMMMTEI=$SCR/$JOB.F39
export SOCCDAT=$SCR/$JOB.F40
export  AABB41=$SCR/$JOB.F41
export  BBAA42=$SCR/$JOB.F42
export  BBBB43=$SCR/$JOB.F43
export  REMD=$SCR/$JOB.F44
export  MCQD50=$SCR/$JOB.F50
export  MCQD51=$SCR/$JOB.F51
export  MCQD52=$SCR/$JOB.F52
export  MCQD53=$SCR/$JOB.F53
export  MCQD54=$SCR/$JOB.F54
export  MCQD55=$SCR/$JOB.F55
export  MCQD56=$SCR/$JOB.F56
export  MCQD57=$SCR/$JOB.F57
export  MCQD58=$SCR/$JOB.F58
export  MCQD59=$SCR/$JOB.F59
export  MCQD60=$SCR/$JOB.F60
export  MCQD61=$SCR/$JOB.F61
export  MCQD62=$SCR/$JOB.F62
export  MCQD63=$SCR/$JOB.F63
export  MCQD64=$SCR/$JOB.F64
export NMRINT1=$SCR/$JOB.F61
export NMRINT2=$SCR/$JOB.F62
export NMRINT3=$SCR/$JOB.F63
export NMRINT4=$SCR/$JOB.F64
export NMRINT5=$SCR/$JOB.F65
export NMRINT6=$SCR/$JOB.F66
export DCPHFH2=$SCR/$JOB.F67
export DCPHF21=$SCR/$JOB.F68
export ELNUINT=$SCR/$JOB.F67
export NUNUINT=$SCR/$JOB.F68
export   GVVPT=$SCR/$JOB.F69
export NUMOIN=$SCR/$JOB.F69
export NUMOCAS=$SCR/$JOB.F70
export NUELMO=$SCR/$JOB.F71
export NUELCAS=$SCR/$JOB.F72

#    next files are for RI-MP2
export RIVMAT=$SCR/$JOB.F51
export RIT2A=$SCR/$JOB.F52
export RIT3A=$SCR/$JOB.F53
export RIT2B=$SCR/$JOB.F54
export RIT3B=$SCR/$JOB.F55

#    Next files are for GMCQDPT
export GMCREF=$SCR/$JOB.F70
export GMCO2R=$SCR/$JOB.F71
export GMCROC=$SCR/$JOB.F72
export GMCOOC=$SCR/$JOB.F73
export GMCCC0=$SCR/$JOB.F74
export GMCHMA=$SCR/$JOB.F75
export GMCEI1=$SCR/$JOB.F76
export GMCEI2=$SCR/$JOB.F77
export GMCEOB=$SCR/$JOB.F78
export GMCEDT=$SCR/$JOB.F79
export GMCERF=$SCR/$JOB.F80
export GMCHCR=$SCR/$JOB.F81
export GMCGJK=$SCR/$JOB.F82
export GMCGAI=$SCR/$JOB.F83
export GMCGEO=$SCR/$JOB.F84
export GMCTE1=$SCR/$JOB.F85
export GMCTE2=$SCR/$JOB.F86
export GMCHEF=$SCR/$JOB.F87
export GMCMOL=$SCR/$JOB.F88
export GMCMOS=$SCR/$JOB.F89
export GMCWGT=$SCR/$JOB.F90
export GMCRM2=$SCR/$JOB.F91
export GMCRM1=$SCR/$JOB.F92
export GMCR00=$SCR/$JOB.F93
export GMCRP1=$SCR/$JOB.F94
export GMCRP2=$SCR/$JOB.F95
export GMCVEF=$SCR/$JOB.F96
export GMCDIN=$SCR/$JOB.F97
export GMC2SZ=$SCR/$JOB.F98
export GMCCCS=$SCR/$JOB.F99
#
#    Next files are used only during closed shell coupled cluster runs.
#    Display the numerous definitions iff they are going to be used.
#
export  CCREST=$SCR/$JOB.F70
export  CCDIIS=$SCR/$JOB.F71
export  CCINTS=$SCR/$JOB.F72
export CCT1AMP=$SCR/$JOB.F73
export CCT2AMP=$SCR/$JOB.F74
export CCT3AMP=$SCR/$JOB.F75
export    CCVM=$SCR/$JOB.F76
export    CCVE=$SCR/$JOB.F77
export CCQUADS=$SCR/$JOB.F78
export QUADSVO=$SCR/$JOB.F79
export EOMSTAR=$SCR/$JOB.F80
export EOMVEC1=$SCR/$JOB.F81
export EOMVEC2=$SCR/$JOB.F82
export  EOMHC1=$SCR/$JOB.F83
export  EOMHC2=$SCR/$JOB.F84
export EOMHHHH=$SCR/$JOB.F85
export EOMPPPP=$SCR/$JOB.F86
export EOMRAMP=$SCR/$JOB.F87
export EOMRTMP=$SCR/$JOB.F88
export EOMDG12=$SCR/$JOB.F89
export    MMPP=$SCR/$JOB.F90
export   MMHPP=$SCR/$JOB.F91
export MMCIVEC=$SCR/$JOB.F92
export MMCIVC1=$SCR/$JOB.F93
export MMCIITR=$SCR/$JOB.F94
export  EOMVL1=$SCR/$JOB.F95
export  EOMVL2=$SCR/$JOB.F96
export EOMLVEC=$SCR/$JOB.F97
export  EOMHL1=$SCR/$JOB.F98
export  EOMHL2=$SCR/$JOB.F99
export  CCVVVV=$SCR/$JOB.F80
#
#    Next files are used only during open shell coupled cluster runs.
#
export AMPROCC=$SCR/$JOB.F70
export ITOPNCC=$SCR/$JOB.F71
export FOCKMTX=$SCR/$JOB.F72
export  LAMB23=$SCR/$JOB.F73
export   VHHAA=$SCR/$JOB.F74
export   VHHBB=$SCR/$JOB.F75
export   VHHAB=$SCR/$JOB.F76
export    VMAA=$SCR/$JOB.F77
export    VMBB=$SCR/$JOB.F78
export    VMAB=$SCR/$JOB.F79
export    VMBA=$SCR/$JOB.F80
export  VHPRAA=$SCR/$JOB.F81
export  VHPRBB=$SCR/$JOB.F82
export  VHPRAB=$SCR/$JOB.F83
export  VHPLAA=$SCR/$JOB.F84
export  VHPLBB=$SCR/$JOB.F85
export  VHPLAB=$SCR/$JOB.F86
export  VHPLBA=$SCR/$JOB.F87
export    VEAA=$SCR/$JOB.F88
export    VEBB=$SCR/$JOB.F89
export    VEAB=$SCR/$JOB.F90
export    VEBA=$SCR/$JOB.F91
export   VPPPP=$SCR/$JOB.F92
export INTERM1=$SCR/$JOB.F93
export INTERM2=$SCR/$JOB.F94
export INTERM3=$SCR/$JOB.F95
export ITSPACE=$SCR/$JOB.F96
export INSTART=$SCR/$JOB.F97
export  ITSPC3=$SCR/$JOB.F98
#
#    Next files are used only during elongation method runs.
#    Display the numerous definitions iff they are going to be used.
set +x
if [ -e $SCR/$JOB.F05 ]; then
  elgtyp=`grep -i NELONG= $SCR/$JOB.F05 | wc -l`
  if [ ! "$elgtyp" -gt "0" ]; then
      ELGNAME=$4
      if [ ! -z "$ELGNAME" ]; then
        ELGNAME=ELGFILE
      fi
      if [ "$VERBOSE" = true ] ; then
        set -x
      fi
      export AOINTS=$SCR/$ELGNAME.F08
      export ELGDOS=$USERSCR/$JOB.ldos
      export ELGDAT=$SCR/$ELGNAME.F71
      export ELGPAR=$SCR/$ELGNAME.F72
      export ELGCUT=$SCR/$ELGNAME.F74
      export ELGVEC=$SCR/$ELGNAME.F75
      export EGINTA=$SCR/$ELGNAME.F77
      export EGINTB=$SCR/$ELGNAME.F78
      export EGTDHF=$SCR/$ELGNAME.F79
      export EGTEST=$SCR/$ELGNAME.F80
      set +x
  fi
fi
#
#    Next files are used only during extended TDHF package runs.
#    Display the numerous definitions iff they are going to be used.
if [ "$VERBOSE" = true ] ; then
  set -x
fi
export  OLI201=$SCR/$JOB.F201
export  OLI202=$SCR/$JOB.F202
export  OLI203=$SCR/$JOB.F203
export  OLI204=$SCR/$JOB.F204
export  OLI205=$SCR/$JOB.F205
export  OLI206=$SCR/$JOB.F206
export  OLI207=$SCR/$JOB.F207
export  OLI208=$SCR/$JOB.F208
export  OLI209=$SCR/$JOB.F209
export  OLI210=$SCR/$JOB.F210
export  OLI211=$SCR/$JOB.F211
export  OLI212=$SCR/$JOB.F212
export  OLI213=$SCR/$JOB.F213
export  OLI214=$SCR/$JOB.F214
export  OLI215=$SCR/$JOB.F215
export  OLI216=$SCR/$JOB.F216
export  OLI217=$SCR/$JOB.F217
export  OLI218=$SCR/$JOB.F218
export  OLI219=$SCR/$JOB.F219
export  OLI220=$SCR/$JOB.F220
export  OLI221=$SCR/$JOB.F221
export  OLI222=$SCR/$JOB.F222
export  OLI223=$SCR/$JOB.F223
export  OLI224=$SCR/$JOB.F224
export  OLI225=$SCR/$JOB.F225
export  OLI226=$SCR/$JOB.F226
export  OLI227=$SCR/$JOB.F227
export  OLI228=$SCR/$JOB.F228
export  OLI229=$SCR/$JOB.F229
export  OLI230=$SCR/$JOB.F230
export  OLI231=$SCR/$JOB.F231
export  OLI232=$SCR/$JOB.F232
export  OLI233=$SCR/$JOB.F233
export  OLI234=$SCR/$JOB.F234
export  OLI235=$SCR/$JOB.F235
export  OLI236=$SCR/$JOB.F236
export  OLI237=$SCR/$JOB.F237
export  OLI238=$SCR/$JOB.F238
export  OLI239=$SCR/$JOB.F239
#
#    Next files are used only during divide-and-conquer runs
#
export   DCSUB=$SCR/$JOB.F250
export   DCVEC=$SCR/$JOB.F251
export   DCEIG=$SCR/$JOB.F252
export   DCDM=$SCR/$JOB.F253
export   DCDMO=$SCR/$JOB.F254
export   DCQ=$SCR/$JOB.F255
export   DCW=$SCR/$JOB.F256
export   DCEDM=$SCR/$JOB.F257
#
#    Next files are used only during LMO hyperpolarizability analysis
#
export LHYPWRK=$SCR/$JOB.F297
export LHYPWK2=$SCR/$JOB.F298
export BONDDPF=$SCR/$JOB.F299
#
#    Next two values are used only within the VB2000 add-on code
#
export VB2000PATH=$GMSPATH/vb2000
export GMSJOBNAME=$JOB
set +x
#
#    Next files are used during EFMO runs
#
if [ -e $SCR/$JOB.F05 ]; then
  efmo=`grep -i IEFMO= $SCR/$JOB.F05 | wc -l`
  if [ ! -z "$efmo" ]; then
    if [ "$VERBOSE" = true ] ; then
      set -x
    fi
    export EFMOI=$SCR/$JOB.F102
    export EFMOF=$SCR/$JOB.F103
    set +x
  fi
fi
#
#    Next files are used for explicitly correlated methods
#
if [ -e $SCR/$JOB.F05 ]; then
  pt2r12=`egrep -i '(PT212=.TRUE.|PT2R12=.T.)' $SCR/$JOB.F05 | wc -l`
  if [ ! -z "$pt2r12" ]; then
   interface=`egrep -i  '(RUNR12=.T.|RUNTYP=.TRUE.|SINGLS=.T.|SINGLES=.TRUE)' $SCR/$JOB.F05 | wc -l`
  if [ ! -z "$interface" ]; then
    if [ "$VERBOSE" = true ] ; then
      set -x
    fi
    export R12INP=$SCR/$JOB
    export PT2BAS=$SCR/$JOB.F300
  # export PT2R12_LOC /path/to/mpqc
    export FILE_LOC=$SCR
    set +x
   else 
    if [ "$VERBOSE" = true ] ; then
      set -x
    fi
    export R12INP=$USERSCR/$JOB
    export PT2BAS=$SCR/$JOB.F300
    set +x
   fi
  fi
fi
#
#    Next files are used only during CIM runs
#
if [ "$VERBOSE" = true ] ; then
  set -x
fi
export CIMFILE=$USERSCR/$JOB.cim
export CIMDMN=$USERSCR/$JOB.dmn
export CIMDCT=$SCR/$JOB.Fcdt
export CIMAOI=$SCR/$JOB.Fcao
export CIMMOI=$SCR/$JOB.Fcmo
#
#    CASINO trial wavefunction file
#
export  CASINO=$USERSCR/$JOB.casino
#    Next values for for spin-space 2nd order density.
export  DEN2P1 $SCR/$JOB.F70
export  DEN2P2 $SCR/$JOB.F71
export  DEN2P3 $SCR/$JOB.F72
export  DEN2P4 $SCR/$JOB.F73
set +x
