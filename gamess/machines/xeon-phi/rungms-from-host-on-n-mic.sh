#!/bin/sh

##############################################################
#  Usage: This script will allow you to run GAMESS natively  #
#         on the Intel Xeon Phi by launching this script     #
#         from the host machine.  You must set the number of #
#         number of MIC in NUMMICS and populate the miccards #
#         array.                                             #
##############################################################

#  Control verbosity of the output
VERBOSE=true

#  Set your directories
#  SCR is your scratch directory - use fast storage medium
SCR=/work2/quantum/sarom-ssok1/scratch

#  USERSCR is your directory to store restart files
USERSCR=/work2/quantum/sarom-ssok1/restart

#  Some trickery below:
#  grep | trim spaces | space deliminate the results and extract the 3rd field

#  Get GMSPATH from install.info
#  i.e. GMSPATH=/work1/quantum/sarom-ssok1/GAMESS-simple-build-intel-mic
GMSPATH=`grep 'setenv GMS_PATH' install.info | tr -s ' ' | cut -d ' ' -f 3`

#  Get GMS_MPI_PATH from install.info
#  i.e. GMS_MPI_PATH=/shared/intel/impi_latest
GMS_MPI_PATH=`grep 'setenv GMS_MPI_PATH' install.info | tr -s ' ' | cut -d ' ' -f 3`

#  Get Intel Compiler MIC library path from install.info
#  Only works if MKL path includes specific version information: 
#  i.e. composer_xe_2015.1.133 or compilers_and_libraries_2016.1.150
#  i.e. INTEL_PATH=/shared/intel/composer_xe_2015.1.133
INTEL_PATH=`grep 'setenv GMS_MATHLIB_PATH' install.info | tr -s ' ' | cut -d ' ' -f 3 | sed 's/mkl\/lib\/intel64//g'`

#  Source compilervars.sh to setup Intel copiler variables
source $INTEL_PATH/bin/compilervars.sh intel64

#  Source mpivars.sh to setup IMPI variables
source $GMS_MPI_PATH/intel64/bin/mpivars.sh

#  Set LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$INTEL_PATH/compiler/lib/mic:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$INTEL_PATH/mkl/lib/mic:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$GMS_MPI_PATH/mic/lib:$LD_LIBRARY_PATH

#  Enable Intel Xeon Phi coprocessor recognition
export I_MPI_MIC=1

#  Set MIC Proxy Path
#  export I_MPI_MIC_PROXY_PATH=$INTEL_PATH/mpirt/bin/mic
export I_MPI_MIC_PROXY_PATH=$GMS_MPI_PATH/mic/bin

#  Set network fabric for the Intel Xeon Phi
#  I_MPI_FABRICS=<intra-node fabric>:<inter-nodes fabric>
#  shm, tcp, ofa, dapl
#  
#  Note: Working on Cyence: I_MPI_FABRICS=shm:tcp
#  Note: Working on Bolt  : I_MPI_FABRICS=
export I_MPI_FABRICS=shm:tcp

#  Set MKL parameter
export MKL_NUM_THREADS=4

#  Set the number of MIC cards we are using
NUMMICS=N

#  Specify available mic nodes by name
miccards[0]=mic0
miccards[1]=mic1
#...
miccards[N]=mic(N-1)

#  Or specify via IPs
#miccards[0]=172.16.31.1
#miccards[1]=172.16.32.1
#...
#miccards[N]=172.16.32.(N-1)

##############################################################
#  You should not need to modify anything below              #
##############################################################
#  Last update = 26 Feb 2016
#
#  This is a shell script to execute GAMESS, by typing
#       rungms-mic.sh JOB VERNO NCPUS >& JOB.log &
#  JOB    is the name of the 'JOB.inp' file to be executed,
#  VERNO  is the number of the executable you chose at 'lked' time,
#  NCPUS  is the number of processors to be used.

TARGET=mpi  # hard-coded to mpi since we are using mpi on Intel Xeon Phi
JOB=$1      # name of the input file xxx.inp, give only the xxx part
VERNO=$2    # revision number of the executable created by 'lked' step
NCPUS=$3    # number of compute processes to be run

#  Provide defaults if last two arguments are not given to this script
if [ -z "$VERNO" ]; then
   VERNO=00
fi

if [ -z "$NCPUS" ]; then
   NCPUS=1
fi

#  ---- the top third of the script is input and other file assignments ----
#
echo "----- GAMESS execution script 'rungms' -----"
master=`hostname`
echo This job is running on host $master
echo under operating system `uname` at `date`

echo "Available scratch disk space (Kbyte units) at beginning of the job is"
df -k $SCR
echo "GAMESS temporary binary files will be written to $SCR"
echo "GAMESS supplementary output files will be written to $USERSCR"

#        this added as experiment, February 2007, as 8 MBytes
#        increased to 32 MB in October 2013 for the VB2000 code.
#        its intent is to detect large arrays allocated off the stack
#limit stacksize 32768

#  Grab a copy of the input file.
#  In the case of examNN jobs, file is in tests/standard subdirectory.
#  In the case of exam-vbNN jobs, file is in vb2000's tests subdirectory.
filenameNOPATH=$(basename "$JOB")
filename="${filenameNOPATH%.*}"
extension="${filenameNOPATH##*.}"

#  Strip off possible .inp
if [ "$extension" == "inp" ]; then
   JOB=$filename
fi

echo "Copying input file $JOB.inp to your run's scratch directory..."
if [ -e $JOB.inp ]; then
   cp  $JOB.inp  $SCR/$JOB.F05
else
   if [ -e tests/standard/$JOB.inp ]; then
      cp  tests/standard/$JOB.inp  $SCR/$JOB.F05
   else
      if [ -e tests/$JOB.inp ]; then
         cp  tests/$JOB.inp  $SCR/$JOB.F05
      else
         echo "Input file $JOB.inp does not exist."
         echo "This job expected the input file to be in directory `pwd`"
         echo "Please fix your file name problem, and resubmit."
         exit 4
    fi
  fi
fi

#    define many environment variables setting up file names.
#    anything can be overridden by a user's own choice, read 2nd.
#
source $GMSPATH/gms-files.sh
if [ -e $HOME/.gmsrc ]; then
   echo "reading your own $HOME/.gmsrc"
   source $HOME/.gmsrc
fi

#    Choose remote shell execution program.
#    Parallel run do initial launch of GAMESS on remote nodes by the
#    following program.  Note that the authentication keys for ssh
#    must have been set up correctly.
#    If you wish, choose 'rsh/rcp' using .rhosts authentication instead.
export DDI_RSH=ssh
export DDI_RCP=scp

#    Remove data left over from a previous run.
#
if [ -e $PUNCH ]; then
   echo Removing $PUNCH
   rm -rf $PUNCH
fi
if [ -e $MAKEFP ]; then
   echo Removing $MAKEFP
   rm -rf $MAKEFP
fi
if [ -e $TRAJECT ]; then
   echo Removing $TRAJECT
   rm -rf $TRAJECT
fi
if [ -e $RESTART ]; then
   echo Removing $RESTART
   rm -rf $RESTART
fi

if [ "$TARGET" == "mpi" ]; then

#  Processors Per Node
   PPN=$4

   if [ -z "$PPN" ]; then
      PPN=$NCPUS
   fi

#  NPROCS will always be equal to double the number of compute processes requested
#  We multiple this number by the $NUMMICS
   NPROCS=$((($NCPUS + $NCPUS) * $NUMMICS))

#  Set MPI path
   DDI_MPI_CHOICE=impi
   DDI_MPI_ROOT=$GMS_MPI_PATH/mic
   path=$DDI_MPI_ROOT/bin
   hash -r

#  Select MPI kick-off style. Both should work.
   MPI_KICKOFF_STYLE=hydra

   export HOSTFILE=$SCR/$JOB.nodes.mpd
   if [ -e $HOSTFILE ]; then
     rm $HOSTFILE
   fi
   touch $HOSTFILE

   if [ "$NUMMICS" = 1 ] ; then
     echo ${miccards[0]} >> $HOSTFILE
   else
     for mic in "${miccards[@]}"
     do
       echo $mic >> $HOSTFILE
     done
   fi

#  Host configuration.
   if [ "$VERBOSE" = true ] ; then
     echo '-----debug----'
     echo HOSTFILE $HOSTFILE contains
     cat $HOSTFILE
     echo '--------------'
   fi

#  B. the next file forces explicit "which process on what node" rules.
#     The contents depend on the kickoff style.  This file is how
#     we tell MPI to double-book the cores with two processes,
#     thus accounting for both compute processes and data servers.
#
   export PROCFILE=$SCR/$JOB.processes.mpd
   if [ -e $PROCFILE ]; then
     rm $PROCFILE
   fi
   touch $PROCFILE

   if [ "$MPI_KICKOFF_STYLE" == "hydra" ]; then
      if [ "$NUMMICS" -eq 1 ]; then
         PPN2=$(expr $PPN + $PPN)
         echo ${miccards[0]}":$PPN2" > $PROCFILE
      else
         PPN2=$(expr $PPN + $PPN)
         for mic in "${miccards[@]}"
         do
           echo $mic":$PPN2" >> $PROCFILE
         done
      fi
   fi

#  Host configuration.
   if [	"$VERBOSE" = true ] ; then
     echo '-----debug----'
     echo PROCFILE $PROCFILE contains
     cat $PROCFILE
     echo '--------------'
   fi

#     ==== values that influence the MPI operation ====
#
#     tunings below are specific to Intel MPI 3.2 and/or 4.0:
#        a very important option avoids polling for incoming messages
#           which allows us to compile DDI in pure "mpi" mode,
#           and get sleeping data servers if the run is SCF level.
#        trial and error showed process pinning slows down GAMESS runs,
#        set debug option to 5 to see messages while kicking off,
#        set debug option to 200 to see even more messages than that,
#        set statistics option to 1 or 2 to collect messaging info,
#        iMPI 4.0 on up defaults fabric to shm,dapl: dapl only is faster.
#
   if [ $DDI_MPI_CHOICE == impi ]; then

#     Turns out that you need to disable I_MPI_WAIT_MODE to prevent
#     GAMESS from hanging after a job has completed (NORMALLY or ABNORMALLY)
#
#     export I_MPI_WAIT_MODE=enable
      export I_MPI_PIN=disable
      export I_MPI_DEBUG=0
      export I_MPI_STATS=0

#              next two select highest speed mode of an Infiniband
#     export I_MPI_FABRICS=dapl
#     export I_MPI_DAT_LIBRARY=libdat2.so

#              next two select TCP/IP, a slower way to use Infiniband.
#              The device could be eth0 if IP over IB is not enabled.
#     setenv I_MPI_FABRICS tcp
#     setenv I_MPI_TCP_NETMASK ib0

#              in case someone wants to try the "tag matching interface",
#              an option which unfortunately ignores the WAIT_MODE in 4.0.2!
#     setenv I_MPI_FABRICS tmi
#     setenv I_MPI_TMI_LIBRARY libtmi.so
#     setenv I_MPI_TMI_PROVIDER psm
#     setenv TMI_CONFIG $DDI_MPI_ROOT/etc/tmi.conf
   fi

#   =========== runtime path/library setup is now finished! ===========
#     any issues with paths and libraries can be debugged just below:
#
   if [	"$VERBOSE" = true ] ; then
     echo "-----debug----"
     echo the execution path is
     echo $path
     echo " "
     echo the library path is
     echo $LD_LIBRARY_PATH
     echo " "
     echo The dynamically linked libraries for this binary are
     ldd $GMSPATH/gamess.$VERNO.x
     echo "--------------"
   fi

#  Now, at last, we can actually kick-off the MPI processes...
#
   echo "MPI kickoff will run GAMESS on $NCPUS cores in $NNODES nodes."
   echo "The binary to be executed is $GMSPATH/gamess.$VERNO.x"
   echo "MPI will run $NCPUS compute processes and $NCPUS data servers,"
   echo "    placing $PPN of each process type onto each node."
   echo "The scratch disk space on each node is $SCR, with free space"
   df -k $SCR

   cd $SCR

#  Hydra
   if [ "$MPI_KICKOFF_STYLE" == "hydra" ]; then
      if [ "$DDI_MPI_CHOICE" == "impi" ]; then
         export I_MPI_HYDRA_ENV=all
         export I_MPI_PERHOST=$PPN2
      fi

      set -x

      mpiexec.hydra -f $PROCFILE -n $NPROCS \
            $GMSPATH/gamess.$VERNO.x < /dev/null

      mpiexec.hydra -f $PROCFILE -n $NPROCS \
            /home/ssok1/bin/kill-stray-mic < /dev/null

      unset

   fi
   rm -f $PROCFILE
fi

#  ---- the bottom third of the script is to clean up all disk files ----
#  It is quite useful to display to users how big the disk files got to be.
#
echo ----- accounting info -----

#   Clean up the master's scratch directory.
#
echo Files used on the master node $master were:
ls -lF $SCR/$JOB.*
rm -f  $SCR/$JOB.F*

#   Clean up scratch directory of remote nodes.
#
#    This particular example is for the combination iMPI, w/SGE or PBS.
#    We have inherited a file of unique node names from above.
#    There is an option to rescue the output files from group DDI runs,
#    such as FMO, in case you need to see the other group's outputs.
#          clean off the last file on the master's scratch disk.
rm -f $HOSTFILE

#    And this is the end
#
date
exit
