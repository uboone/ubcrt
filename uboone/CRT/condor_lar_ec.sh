#! /bin/bash
#------------------------------------------------------------------
#
# Purpose: A general purpose larsoft batch worker script.
#
# Adapted from condor_lBdetMC.sh by E. Church.
#
# Usage:
#
# condor_lar.sh [options]
#
# Lar options:
#
# -c, --config <arg>      - Configuration (fcl) file (required).
# -s, --source <arg>      - Input file (full path).
# -S, --source-list <arg> - Input file list (full path, one per line).
# -o, --output <arg>      - Output file name.
# -T, --TFileName  <arg>  - TFile output file name
# -n, --nevts <arg>       - Number of events to process.
# --nskip <arg>           - Number of events to skip.
# --nfile <arg>           - Number of files to process per worker.
# --nfile_skip <arg>      - Number of files to skip (use with option -S).
# --inputmode <arg>       - Input mode ('textfile' or '', default '')
# --args <args...>        - Arguments for lar command line (place at end).
#
# Sam and parallel project options.
#
# --sam_user <arg>        - Specify sam user (default $GRID_USER).
# --sam_group <arg>       - Specify sam group (default --group option).
# --sam_station <arg>     - Specify sam station (default --group option).
# --sam_defname <arg>     - Sam dataset definition name.
# --sam_project <arg>     - Sam project name.
# --sam_start             - Specify that this worker should be responsible for
#                           starting and stopping the sam project.
# --sam_schema <arg>      - Use this option with argument "root" to stream files using
#                           xrootd.  Leave this option out for standard file copy.
# --njobs <arg>           - Parallel project with specified number of jobs (default one).
# --single                - Specify that the output and log directories will be emptied
#                           by the batch worker, and therefore the output and log
#                           directories will only ever contain output from a single
#                           worker.
#
# Mix input options (second input stream).
#
# --mix_defname <arg>     - Specify mix input sam dataset definition.
# --mix_project <arg>     - Specify mix input sam project.
#
# Larsoft options.
#
# --ups <arg>             - Comma-separated list of top level run-time ups products.
# -r, --release <arg>     - Release tag.
# -q, -b, --build <arg>   - Release build qualifier (default "debug", or "prof").
# --localdir <arg>        - Larsoft local test release directory (default none).
# --localtar <arg>        - Tarball of local test release.
# --mrb                   - Ignored (for compatibility).
# --srt                   - Exit with error status (SRT run time no longer supported).
#
# Other options.
#
# -h, --help              - Print help.
# -i, --interactive       - For interactive use.
# -g, --grid              - Be grid-friendly.
# --group <arg>           - Group or experiment (required).
# --workdir <arg>         - Work directory (required).
# --outdir <arg>          - Output directory (required).
# --logdir <arg>          - Log directory (required).
# --scratch <arg>         - Scratch directory (only for interactive).
# --cluster <arg>         - Job cluster (override $CLUSTER)
# --process <arg>         - Process within cluster (override $PROCESS).
# --procmap <arg>         - Name of process map file (override $PROCESS).
# --init-script <arg>     - User initialization script execute.
# --init-source <arg>     - User initialization script to source (bash).
# --end-script <arg>      - User end-of-job script to execute.
# --exe <arg>             - Specify art-like executable (default "lar").
#
# End options.
#
# Run time environment setup.
#
# MRB run-time environmental setup is controlled by four options:
#  --release (-r), --build (-b, -q), --localdir, and --localtar.  
#
# a) Use option --release or -r to specify version of top-level product(s).  
# b) Use option --build or -b to specify build full qualifiers (e.g. 
#    "debug:e5" or "e5:prof").
# c) Options --localdir or --localtar are used to specify your local
#    test release.  Use one or the other (not both).
#
#    Use --localdir to specify the location of your local install
#    directory ($MRB_INSTALL).
#
#    Use --localtar to specify thye location of a tarball of your
#    install directory (made relative to $MRB_INSTALL).
#
#    Note that --localdir is not grid-friendly.
#
# Notes.
#
# 1.  Each batch worker is uniquely identified by two numbers stored
#     in environment variables $CLUSTER and $PROCESS (the latter is 
#     a small integer that starts from zero and varies for different
#     jobs in a parallel job group).  These environment variables are
#     normally set by the batch system, but can be overridden by options 
#     --cluster, --process, and --procmap (e.g. to rerun failed jobs).
#
# 2.  The work directory must be set to an existing directory owned
#     by the submitter and readable by the batch worker.  Files from the 
#     work directory are copied to the batch worker scratch directory at
#     the start of the job.
#
# 3.  The job configuration file (-c option), initialization and end-of-job
#     scripts (optins --init-script, --init-source, --end-script) may
#     be stored in the work directory specified by option --workdir, or they
#     may be specified as absolute paths visible on the worker node.
#
# 4.  A local test release may be specified as an absolute path using
#     --localdir, or a tarball using --localtar.  The location of the tarball
#     may be specified as an absolute path visible on the worker, or a 
#     relative path relative to the work directory.
#
# 5.  The output directory must exist and be writable by the batch
#     worker (i.e. be group-writable for grid jobs).  The worker
#     makes a new subdirectory called ${CLUSTER}_${PROCESS} in the output
#     directory and copies all files in the batch scratch directory there 
#     at the end of the job.  If the output directory is not specified, the
#     default is /grid/data/<group>/outstage/<user> (user is defined as 
#     owner of work directory).
#
# 6.  Parallel projects are specified whenever --njobs is specified to
#     be greater than one.  Parallel projects are supported for single file,
#     file list, and sam project input.
#
#     In all cases, each worker processes some number of complete files.
#     If the number of jobs is greater than the number of input files, some
#     workers will not have any input files to process.
#
#     In any case, options --nfile and --nevts can be used to limit the
#     number of files or events that are processed by a single worker, 
#     regardless of the way files are divided among the workers.
#
#     Option --njobs is incompatible with options --nskip, and --nfile_skip.
#
#     a) Non-sam (single file or file list) input.
#
#     In this case, input files are preassigned to workers such that all input
#     files are approximately evenly divided among the workers.  All files
#     preassigned to this worker are copied to the scratch directory at the 
#     start of the job.
#
#     b) Sam project input.
#
#     In this case, files are assigned to workers in a non-deterministic
#     manner by the sam system.  The sam system fetches input files to the
#     scratch directory and deletes processed input files during job execution.
#
#
# 7.  Using option -n or --nevts to limit number of events processed: 
#
#     a) If no input files are specified (e.g. mc generation), --nevts
#        specifies total number of events among all workers.
#
#     b) If input files are specified, --nevts specifies total number of 
#        events processed by each worker or from each input file, whichever
#        is less.
#
# 8.  The interactive option (-i or --interactive) allows this script
#     to be run interactively by overriding some settings that are normally
#     obtained from the batch system, including $CLUSTER, $PROCESS, and
#     the scratch directory.  Interactive jobs always set PROCESS=0 (unless
#     overridden by --process).
#
# 9.  The grid option (-g or --grid) instructs this script to use grid-
#     friendly tools.  This means that there must be no direct access to
#     bluearc disks.  File transfers are done using gridftp or other
#     grid-friendly protocol.  Local test releases are not allowed to 
#     be specified as directories (--localdir), but may be specified as
#     tarballs (--localtar).
#
# 10. Mix options (--mix_defname, --mix_project) are only partially handled
#     in this script.  These options are parsed and their values are stored
#     in shell variables.  It is assumed that the sam project specified
#     by --mix_project has been started externally, unless --sam_start is
#     also specified, in which case this script will start the project.  
#     This script does not include any provision for joining the project.
#     Further processing of these options (joining sam project, generating
#     command line options or fcl wrappers) should be handled by user
#     provided initialization scripts (--init-script, --init-source).
#
#
# Created: H. Greenlee, 29-Aug-2012
#
#------------------------------------------------------------------

cd

# Parse arguments.

FCL=""
INFILE=""
INLIST=""
INMODE=""
OUTFILE=""
TFILE=""
NEVT=0
NSKIP=0
FIRST_EVENT=0
SUBRUN=1
NFILE=0
NFILE_SKIP=0
NJOBS=1
SINGLE=0
ARGS=""
UPS_PRDS=""
REL=""
QUAL=""
LOCALDIR=""
LOCALTAR=""
INTERACTIVE=0
GRP=""
WORKDIR=""
OUTDIR=""
LOGDIR=""
SCRATCH=""
CLUS=""
PROC=""
PROCMAP=""
INITSCRIPT=""
INITSOURCE=""
ENDSCRIPT=""
SAM_USER=$GRID_USER
SAM_GROUP=""
SAM_STATION=""
SAM_DEFNAME=""
SAM_PROJECT=""
SAM_START=0
SAM_SCHEMA=""
USE_SAM=0
MIX_DEFNAME=""
MIX_PROJECT=""
MIX_SAM=0
GRID=0
IFDH_OPT=""
DECLARE_IN_JOB=0
VALIDATE_IN_JOB=0
COPY_TO_FTS=0
MAINTAIN_PARENTAGE=0
EXE="lar"

while [ $# -gt 0 ]; do
  case "$1" in

    # Help.
    -h|--help )
      awk '/^# Usage:/,/^# End options/{print $0}' $0 | cut -c3- | head -n -2
      exit
      ;;

    # Config file.
    -c|--config )
      if [ $# -gt 1 ]; then
        FCL=$2
        shift
      fi
      ;;

    # Input file.
    -s|--source )
      if [ $# -gt 1 ]; then
        INFILE=$2
        shift
      fi
      ;;

    # Input file list. 
    -S|--source-list )
      if [ $# -gt 1 ]; then
        INLIST=$2
        shift
      fi
      ;;

    # Input file mode.
    --inputmode )
      if [ $# -gt 1 ]; then
        INMODE=$2
        shift
      fi
      ;;

    # Output file.
    -o|--output )
      if [ $# -gt 1 ]; then
        OUTFILE=$2
        shift
      fi
      ;;
      
    # Output TFile.
    -T|--TFileName )
      if [ $# -gt 1 ]; then
        TFILE=$2
        shift
      fi
      ;;    

    # Number of events.
    -n|--nevts )
      if [ $# -gt 1 ]; then
        NEVT=$2
        shift
      fi
      ;;

    # Number of events to skip.
    --nskip )
      if [ $# -gt 1 ]; then
        NSKIP=$2
        shift
      fi
      ;;

    # Number of files to process.
    --nfile )
      if [ $# -gt 1 ]; then
        NFILE=$2
        shift
      fi
      ;;

    # Number of files to skip.
    --nfile_skip )
      if [ $# -gt 1 ]; then
        NFILE_SKIP=$2
        shift
      fi
      ;;

    # Number of parallel jobs.
    --njobs )
      if [ $# -gt 1 ]; then
        NJOBS=$2
        shift
      fi
      ;;

    # Single worker mode.
    --single )
      SINGLE=1
      ;;

    # Sam user.
    --sam_user )
      if [ $# -gt 1 ]; then
        SAM_USER=$2
        shift
      fi
      ;;

    # Sam group.
    --sam_group )
      if [ $# -gt 1 ]; then
        SAM_GROUP=$2
        shift
      fi
      ;;

    # Sam station.
    --sam_station )
      if [ $# -gt 1 ]; then
        SAM_STATION=$2
        shift
      fi
      ;;

    # Sam dataset definition name.
    --sam_defname )
      if [ $# -gt 1 ]; then
        SAM_DEFNAME=$2
        USE_SAM=1
        shift
      fi
      ;;

    # Sam project name.
    --sam_project )
      if [ $# -gt 1 ]; then
        SAM_PROJECT=$2
        USE_SAM=1
        shift
      fi
      ;;

    # Sam start/stop project flag.
    --sam_start )
      SAM_START=1
      ;;

    # Sam schema.
    --sam_schema )
      if [ $# -gt 1 ]; then
        SAM_SCHEMA=$2
        shift
      fi
      ;;

    # General arguments for lar command line.
    --args )
      if [ $# -gt 1 ]; then
        shift
        ARGS=$@
        break
      fi
      ;;

    # Top level ups products (comma-separated list).
    --ups )
      if [ $# -gt 1 ]; then
        UPS_PRDS=$2
        shift
      fi
      ;;

    # Release tag.
    -r|--release )
      if [ $# -gt 1 ]; then
        REL=$2
        shift
      fi
      ;;

    # Release build qualifier.
    -q|-b|--build )
      if [ $# -gt 1 ]; then
        QUAL=$2
        shift
      fi
      ;;

    # Local test release directory.
    --localdir )
      if [ $# -gt 1 ]; then
        LOCALDIR=$2
        shift
      fi
      ;;

    # Local test release tarball.
    --localtar )
      if [ $# -gt 1 ]; then
        LOCALTAR=$2
        shift
      fi
      ;;

    # MRB flag.
    --mrb )
      ;;

    # SRT flag.
    --srt )
      echo "SRT run time environment is no longer supported."
      exit 1
      ;;

    # Interactive flag.
    -i|--interactive )
      INTERACTIVE=1
      ;;

    # Grid flag.
    -g|--grid )
      GRID=1
      ;;

    # Group.
    --group )
      if [ $# -gt 1 ]; then
        GRP=$2
        shift
      fi
      ;;

    # Work directory.
    --workdir )
      if [ $# -gt 1 ]; then
        WORKDIR=$2
        shift
      fi
      ;;

    # Output directory.
    --outdir )
      if [ $# -gt 1 ]; then
        OUTDIR=$2
        shift
      fi
      ;;

    # Log directory.
    --logdir )
      if [ $# -gt 1 ]; then
        LOGDIR=$2
        shift
      fi
      ;;

    # Scratch directory.
    --scratch )
      if [ $# -gt 1 ]; then
        SCRATCH=$2
        shift
      fi
      ;;

    # Job cluster.
    --cluster )
      if [ $# -gt 1 ]; then
        CLUS=$2
        shift
      fi
      ;;

    # Process within cluster.
    --process )
      if [ $# -gt 1 ]; then
        PROC=$2
        shift
      fi
      ;;

    # Process map.
    --procmap )
      if [ $# -gt 1 ]; then
        PROCMAP=$2
        shift
      fi
      ;;

    # User initialization script.
    --init-script )
      if [ $# -gt 1 ]; then
        INITSCRIPT=$2
        shift
      fi
      ;;

    # User source initialization script.
    --init-source )
      if [ $# -gt 1 ]; then
        INITSOURCE=$2
        shift
      fi
      ;;

    # User end-of-job script.
    --end-script )
      if [ $# -gt 1 ]; then
        ENDSCRIPT=$2
        shift
      fi
      ;;
    
    # Declare good output root files to SAM.
    --declare )
      DECLARE_IN_JOB=1
      ;;
      
    # Run validation steps in project.py on root outputs directly in the job.
    --validate )
      VALIDATE_IN_JOB=1
      ;;

   # Copy Output to FTS.
    --copy )
      COPY_TO_FTS=1
      ;;

    # Mix input sam dataset.
    --mix_defname )
      if [ $# -gt 1 ]; then
        MIX_DEFNAME=$2
        MIX_SAM=1
        shift
      fi
      ;;

    # Mix input sam project.
    --mix_project )
      if [ $# -gt 1 ]; then
        MIX_PROJECT=$2
        MIX_SAM=1
        shift
      fi
      ;;

    # Alter the output file's parentage such that it's parent(s) are from the input list OR sam process
    --maintain_parentage )
      MAINTAIN_PARENTAGE=1
      ;;  
    
    # Specify alternate art-like executable.
    --exe )
      if [ $# -gt 1 ]; then
        EXE=$2
        shift
      fi
      ;;

    # Other.
    * )
      echo "Unknown option $1"
      exit 1
  esac
  shift
done

#echo "FCL=$FCL"
#echo "INFILE=$INFILE"
#echo "INLIST=$INLIST"
#echo "OUTFILE=$OUTFILE"
#echo "TFILE=$TFILE"
#echo "NEVT=$NEVT"
#echo "NSKIP=$NSKIP"
#echo "NFILE=$NFILE"
#echo "NFILE_SKIP=$NFILE_SKIP"
#echo "NJOBS=$NJOBS"
#echo "ARGS=$ARGS"
#echo "REL=$REL"
#echo "QUAL=$QUAL"
#echo "LOCALDIR=$LOCALDIR"
#echo "LOCALTAR=$LOCALTAR"
#echo "INTERACTIVE=$INTERACTIVE"
#echo "GRP=$GRP"
#echo "WORKDIR=$WORKDIR"
#echo "OUTDIR=$OUTDIR"
#echo "LOGDIR=$LOGDIR"
#echo "SCRATCH=$SCRATCH"
#echo "CLUS=$CLUS"
#echo "PROC=$PROC"
#echo "INITSCRIPT=$INITSCRIPT"
#echo "INITSOURCE=$INITSOURCE"
#echo "ENDSCRIPT=$ENDSCRIPT"
#echo "VALIDATE_IN_JOB=$VALIDATE_IN_JOB"

# Done with arguments.

echo "Nodename: `hostname -f`"
id
echo "Load average:"
cat /proc/loadavg

# Set defaults.

if [ x$QUAL = x ]; then
  QUAL="prof:e9"
fi

if [ x$SAM_GROUP = x ]; then
  SAM_GROUP=$GRP
fi

if [ x$SAM_STATION = x ]; then
  SAM_STATION=$GRP
fi

# Fix for sites with newer linux kernels:

case `uname -r` in
    3.*) export UPS_OVERRIDE="-H Linux64bit+2.6-2.12";;
    4.*) export UPS_OVERRIDE="-H Linux64bit+2.6-2.12";;
esac
echo "uname -r: `uname -r`"
echo "UPS_OVERRIDE: $UPS_OVERRIDE"

# Make sure work directory is defined and exists.

if [ x$WORKDIR = x ]; then
  echo "Work directory not specified."
  exit 1
fi
if [ $GRID -eq 0 -a ! -d $WORKDIR ]; then
  echo "Work directory $WORKDIR does not exist."
  exit 1
fi
echo "Work directory: $WORKDIR"



echo "Condor dir input: $CONDOR_DIR_INPUT"

# Initialize experiment ups products and mrb.

echo "Initializing ups and mrb."

echo "Sourcing setup_experiment.sh"
source ${CONDOR_DIR_INPUT}/setup_experiment.sh

echo PRODUCTS=$PRODUCTS

# Ifdh may already be setup by jobsub wrapper.
# If not, set it up here.

echo "IFDHC_DIR=$IFDHC_DIR"
if [ x$IFDHC_DIR = x ]; then
  echo "Setting up ifdhc, because jobsub did not set it up."
  setup ifdhc
fi
echo "IFDHC_DIR=$IFDHC_DIR"

# Set GROUP environment variable.

unset GROUP
if [ x$GRP != x ]; then
  GROUP=$GRP
else
  echo "GROUP not specified."
  exit 1  
fi
export GROUP
echo "Group: $GROUP"

# Set options for ifdh.

if [ $GRID -ne 0 ]; then

  # Figure out if this is a production job.
  # This option is only used when copying back output.
  # It affects the ownership of copied back files.

  echo "X509_USER_PROXY = $X509_USER_PROXY"
  if ! echo $X509_USER_PROXY | grep -q Production; then
    FORCE=expgridftp
    IFDH_OPT="--force=$FORCE"
  else
    FORCE=gridftp
    IFDH_OPT="--force=$FORCE"
  fi
fi
echo "IFDH_OPT=$IFDH_OPT"

# Make sure fcl file argument was specified.

if [ x$FCL = x ]; then
  echo "No configuration option (-c|--config) was specified."
  exit 1
fi

# Make sure output directory exists and is writable.

if [ x$OUTDIR = x ]; then
  echo "Output directory not specified."
  exit 1
fi
if [ $GRID -eq 0 -a \( ! -d $OUTDIR -o ! -w $OUTDIR \) ]; then
  echo "Output directory $OUTDIR does not exist or is not writable."
  exit 1
fi
echo "Output directory: $OUTDIR"

# Make sure log directory exists and is writable.

if [ x$LOGDIR = x ]; then
  echo "Log directory not specified."
  exit 1
fi
if [ $GRID -eq 0 -a \( ! -d $LOGDIR -o ! -w $LOGDIR \) ]; then
  echo "Log directory $LOGDIR does not exist or is not writable."
  exit 1
fi
echo "Log directory: $LOGDIR"

# See if we need to set umask for group write.

if [ $GRID -eq 0 ]; then
  OUTUSER=`stat -c %U $OUTDIR`
  LOGUSER=`stat -c %U $LOGDIR`
  CURUSER=`whoami`
  if [ $OUTUSER != $CURUSER -o $LOGUSER != $CURUSER ]; then
    echo "Setting umask for group write."
    umask 002
  fi
fi

# Make sure scratch directory is defined.
# For batch, the scratch directory is always $_CONDOR_SCRATCH_DIR
# For interactive, the scratch directory is specified by option 
# --scratch or --outdir.

if [ $INTERACTIVE -eq 0 ]; then
  SCRATCH=$_CONDOR_SCRATCH_DIR
else
  if [ x$SCRATCH = x ]; then
    SCRATCH=$OUTDIR
  fi
fi
if [ x$SCRATCH = x -o ! -d "$SCRATCH" -o ! -w "$SCRATCH" ]; then
  echo "Local scratch directory not defined or not writable."
  exit 1
fi

# Create the scratch directory in the condor scratch diretory.
# Copied from condor_lBdetMC.sh.
# Scratch directory path is stored in $TMP.
# Scratch directory is automatically deleted when shell exits.

# Do not change this section.
# It creates a temporary working directory that automatically cleans up all
# leftover files at the end.
TMP=`mktemp -d ${SCRATCH}/working_dir.XXXXXXXXXX`
TMP=${TMP:-${SCRATCH}/working_dir.$$}

{ [[ -n "$TMP" ]] && mkdir -p "$TMP"; } || \
  { echo "ERROR: unable to create temporary directory!" 1>&2; exit 1; }
trap "[[ -n \"$TMP\" ]] && { cd ; rm -rf \"$TMP\"; }" 0
chmod 755 $TMP
cd $TMP
# End of the section you should not change.

echo "Scratch directory: $TMP"

# Copy files from work directory to scratch directory.

echo "No longer fetching files from work directory."
echo "that's now done with using jobsub -f commands"
mkdir work
cp ${CONDOR_DIR_INPUT}/* ./work/
cd work
find . -name \*.tar -exec tar xf {} \;
find . -name \*.py -exec chmod +x {} \;
find . -name \*.sh -exec chmod +x {} \;
echo "Local working directoroy:"
pwd
ls
echo

# Save the hostname and condor job id.

hostname > hostname.txt
echo ${CLUSTER}.${PROCESS} > jobid.txt

# Set default CLUSTER and PROCESS environment variables for interactive jobs.

if [ $INTERACTIVE -ne 0 ]; then
  CLUSTER=`date +%s`   # From time stamp.
  PROCESS=0            # Default zero for interactive.
fi

# Override CLUSTER and PROCESS from command line options.

if [ x$CLUS != x ]; then
  CLUSTER=$CLUS
fi
if [ x$PROC != x ]; then
  PROCESS=$PROC
fi
if [ x$PROCMAP != x ]; then
  if [ -f $PROCMAP ]; then
    PROCESS=`sed -n $(( $PROCESS + 1 ))p $PROCMAP`
  else
    echo "Process map file $PROCMAP not found."
    exit 1
  fi
fi
if [ x$CLUSTER = x ]; then
  echo "CLUSTER not specified."
  exit 1
fi
if [ x$PROCESS = x ]; then
  echo "PROCESS not specified."
  exit 1
fi
echo "Procmap: $PROCMAP"
echo "Cluster: $CLUSTER"
echo "Process: $PROCESS"

# Construct name of output subdirectory.

OUTPUT_SUBDIR=${CLUSTER}_${PROCESS}
echo "Output subdirectory: $OUTPUT_SUBDIR"

# Make sure fcl file exists.

if [ ! -f $FCL ]; then
  echo "Configuration file $FCL does not exist."
  exit 1
fi

# Make sure init script exists and is executable (if specified).

if [ x$INITSCRIPT != x ]; then
  if [ -f "$INITSCRIPT" ]; then
    chmod +x $INITSCRIPT
  else
    echo "Initialization script $INITSCRIPT does not exist."
    exit 1
  fi
fi

# Make sure init source script exists (if specified).

if [ x$INITSOURCE != x -a ! -f "$INITSOURCE" ]; then
  echo "Initialization source script $INITSOURCE does not exist."
  exit 1
fi

# Make sure end-of-job script exists and is executable (if specified).

if [ x$ENDSCRIPT != x ]; then
  if [ -f "$ENDSCRIPT" ]; then
    chmod +x $ENDSCRIPT
  else
    echo "Initialization script $ENDSCRIPT does not exist."
    exit 1
  fi
fi

# MRB run time environment setup goes here.

# Setup local test release, if any.

if [ x$LOCALDIR != x ]; then
  mkdir $TMP/local
  cd $TMP/local

  # Copy test release directory recursively.

  echo "Copying local test release from directory ${LOCALDIR}."

  # Make sure ifdhc is setup.

  if [ x$IFDHC_DIR = x ]; then
    echo "Setting up ifdhc before fetching local directory."
    setup ifdhc
  fi
  echo "IFDHC_DIR=$IFDHC_DIR"
  ifdh cp -r $IFDH_OPT $LOCALDIR .
  stat=$?
  if [ $stat -ne 0 ]; then
    echo "ifdh cp failed with status ${stat}."
    exit $stat
  fi
  find . -name \*.py -exec chmod +x {} \;
  find . -name \*.sh -exec chmod +x {} \;

  # Setup the environment.

  cd $TMP/work
  echo "Initializing localProducts from ${LOCALDIR}."
  if [ ! -f $TMP/local/setup ]; then
    echo "Local test release directory $LOCALDIR does not contain a setup script."
    exit 1
  fi
  sed "s@setenv MRB_INSTALL.*@setenv MRB_INSTALL ${TMP}/local@" $TMP/local/setup | \
  sed "s@setenv MRB_TOP.*@setenv MRB_TOP ${TMP}@" > $TMP/local/setup.local
  . $TMP/local/setup.local
  #echo "MRB_INSTALL=${MRB_INSTALL}."
  #echo "MRB_QUALS=${MRB_QUALS}."
  echo "Setting up all localProducts."
  if [ x$IFDHC_DIR != x ]; then
    unsetup ifdhc
  fi
  mrbslp
fi
cd $TMP/work

# Setup local larsoft test release from tarball.

if [ x$LOCALTAR != x ]; then
  mkdir $TMP/local
  cd $TMP/local

  # Fetch the tarball.

  echo "Fetching test release tarball ${LOCALTAR}."

  # Make sure ifdhc is setup.

  if [ x$IFDHC_DIR = x ]; then
    echo "Setting up ifdhc before fetching tarball."
    setup ifdhc
  fi
  echo "IFDHC_DIR=$IFDHC_DIR"
  ifdh cp $LOCALTAR local.tar
  stat=$?
  if [ $stat -ne 0 ]; then
    echo "ifdh cp failed with status ${stat}."
    exit $stat
  fi 

  # Extract the tarball.

  tar -xf local.tar

  # Setup the environment.

  cd $TMP/work
  echo "Initializing localProducts from tarball ${LOCALTAR}."
  sed "s@setenv MRB_INSTALL.*@setenv MRB_INSTALL ${TMP}/local@" $TMP/local/setup | \
  sed "s@setenv MRB_TOP.*@setenv MRB_TOP ${TMP}@" > $TMP/local/setup.local
  . $TMP/local/setup.local
  #echo "MRB_INSTALL=${MRB_INSTALL}."
  #echo "MRB_QUALS=${MRB_QUALS}."
  echo "Setting up all localProducts."
  if [ x$IFDHC_DIR != x ]; then
    unsetup ifdhc
  fi
  mrbslp
fi

# Setup specified version of top level run time products
# (if specified, and if local test release did not set them up).

for prd in `echo $UPS_PRDS | tr , ' '`
do
  if ! ups active | grep -q $prd; then
    echo "Setting up $prd $REL -q ${QUAL}."
    if [ x$IFDHC_DIR != x ]; then
      unsetup ifdhc
    fi
    setup $prd $REL -q $QUAL
  fi
done

ups active

cd $TMP/work

# In case mrb setup didn't setup a version of ifdhc, set up ifdhc again.

if [ x$IFDHC_DIR = x ]; then
  echo "Setting up ifdhc again, because larsoft did not set it up."
  setup ifdhc
fi
echo "IFDH_ART_DIR=$IFDH_ART_DIR"
echo "IFDHC_DIR=$IFDHC_DIR"

# Get input files to process, either single file, file list, or sam.
#
# For non-sam input, copy all files local using ifdh cp, and make a 
# local file list called condor_lar_input.list.  Save the remote file names (uri's)
# in another file called transferred_uris.list

rm -f condor_lar_input.list
rm -f transferred_uris.list
NFILE_TOTAL=0
parent_files=()
aunt_files=() #for data overaly, the data files being brought in are the output's aunts. 

if [ $USE_SAM -eq 0 -a x$INFILE != x ]; then

  # Single file case.

  # Don't allow any list-related options in single file case:
  # -S, --source-list, --nfile, --nfile_skip

  if [ x$INLIST != x -o $NFILE -ne 0 -o $NFILE_SKIP -ne 0 ]; then
    echo "File list options specified with single input file."
    exit 1
  fi

  #set the parent file to be the input file
  parent_files=("${parent_files[@]}" $INFILE)

  # Copy input file to scratch directoroy.

  NFILE_TOTAL=1
  LOCAL_INFILE_FULL=$INFILE
  LOCAL_INFILE=`basename $INFILE`
  echo "Copying $INFILE"
  ifdh cp $INFILE $LOCAL_INFILE
  stat=$?
  if [ $stat -ne 0 ]; then
    echo "ifdh cp failed with status ${stat}."
    exit $stat
  fi 
  if [ -f $LOCAL_INFILE -a $stat -eq 0 ]; then
    echo $INFILE > transferred_uris.list
    echo $LOCAL_INFILE_FULL > condor_lar_input.list
  else
    echo "Error fetching input file ${INFILE}."
    exit 1
  fi

elif [ $USE_SAM -eq 0 -a x$INLIST != x ]; then

  # Input list case.

  # Make sure input file list exists.

  if [ ! -f $INLIST ]; then
    echo "Input file list $INLIST does not exist."
    exit 1
  fi

  # Remember how many files are in the input file list.

  NFILE_TOTAL=`cat $INLIST | wc -l`
  echo "Input file list contains $NFILE_TOTAL total files."

  # Clamp the total number of files to be a maximum of NFILE * NJOBS, where
  # NFILE and NJOBS are specified via command line options.  In project.py
  # terms, NFILE is <maxfilesperjob> and NOJBS is <numjobs>.

  MAX_TOTAL=$(( $NFILE * $NJOBS ))
  if [ $MAX_TOTAL -gt 0 -a $NFILE_TOTAL -gt $MAX_TOTAL ]; then
    NFILE_TOTAL=$MAX_TOTAL
    echo "Number of files to be processed will be limited to ${NFILE_TOTAL}."
  fi

  # If --njobs was specified, calculate how many files
  # to skip and process in this worker.

  if [ $NJOBS -ne 0 ]; then

    # Don't allow option --nfile_skip in this case.

    if [ $NFILE_SKIP -ne 0 ]; then
      echo "Illegal options specified with --njobs."
      exit 1
    fi

    # Clamp NJOBS to be a maximum of $NFILE_TOTAL.
    # This means that workers with $PROCESS >= $NFILE_TOTAL will not have 
    # any input files to process.

    MYNJOBS=$NJOBS
    if [ $MYNJOBS -gt $NFILE_TOTAL ]; then
      MYNJOBS=$NFILE_TOTAL
    fi

    # Calculate number of files to skip and number of files to process.

    NFILE_SKIP=$(( $PROCESS * $NFILE_TOTAL / $MYNJOBS ))
    MYNFILE=$(( ( $PROCESS + 1 ) * $NFILE_TOTAL / $MYNJOBS - $NFILE_SKIP ))
    if [ $MYNFILE -eq 0 -o $NFILE_SKIP -ge $NFILE_TOTAL ]; then
      echo "This worker did not get any input files."
      exit 1
    fi
    if [ $MYNFILE -lt $NFILE -o $NFILE -eq 0 ]; then
      NFILE=$MYNFILE
    fi
  fi

  # Report number of files to skip and process.

  echo "Skipping $NFILE_SKIP files."
  if [ $NFILE -eq 0 ]; then
    echo "Processing all remaining files."
  else
    echo "Processing $NFILE files."
  fi

  # Copy input files and construct local input file list.

  nfile=0
  nfskip=$NFILE_SKIP
  nmax=$NFILE
  while read infile; do
    if [ $nfskip -gt 0 ]; then
      nfskip=$(( $nfskip - 1 ))
    else

      # Retain the original file name as the local file name, if possible.
      # Otherwise, generate a new (hopefully) unique name.

      if [ ! -f condor_lar_input.list ]; then
        touch condor_lar_input.list
      fi
      LOCAL_INFILE=`basename $infile`
      if grep -q $LOCAL_INFILE condor_lar_input.list; then
        LOCAL_INFILE=input${nfile}.root
	if [ "$INMODE" = 'textfile' ]; then
	  LOCAL_INFILE=input${nfile}.txt
	fi
      fi
      echo "Copying $infile"
      ifdh cp $infile $LOCAL_INFILE
      stat=$?
      if [ $stat -ne 0 ]; then
        echo "ifdh cp failed with status ${stat}."
        exit $stat
      fi 
      if [ -f $LOCAL_INFILE -a $stat -eq 0 ]; then
        echo $infile >> transferred_uris.list
        echo $LOCAL_INFILE >> condor_lar_input.list
	parent_files=("${parent_files[@]}" $LOCAL_INFILE)
      else
        echo "Error fetching input file ${infile}."
        exit 1
      fi
      nmax=$(( $nmax - 1 ))
      if [ $nmax -eq 0 ]; then
        break
      fi
    fi
    nfile=$(( $nfile + 1 ))
  done < $INLIST
fi

NFILE_LOCAL=0
if [ $USE_SAM -eq 0 ]; then
  if [ -f condor_lar_input.list ]; then

    # Sort input list by decreasing size so we don't get a file with
    # zero events as the first file.

    #ls -S1 `cat condor_lar_input.list` > condor_lar_input.list
    xargs ls -s1 < condor_lar_input.list | sort -nr | awk '{print $2}' > newcondor_lar_input.list
    mv -f newcondor_lar_input.list condor_lar_input.list
    echo "Local input file list:"
    cat condor_lar_input.list
    NFILE_LOCAL=`cat condor_lar_input.list | wc -l`
  else
    echo "No local input files."
  fi
  echo "Local input list has $NFILE_LOCAL files."
fi

#Break the master wrapper fcl into each stage
nfcls=0

while read -r line
do

  if [ "$(echo $line | awk '{print $1}')" = "#---STAGE" ]; then
    stage="$(echo $line | awk '{print $2}')"
    stage_fcl="Stage$stage.fcl"
    nfcls=$(( $nfcls + 1 )) 
    continue
  fi

  if [ "$line" = "#---END_STAGE" ]; then
     #cat EOF >> $fcl
    continue
  fi
  echo $line >> $stage_fcl
done < $FCL

#We now have nStage fcl files, each which need to be run serially 
stage=0

echo "Start loop over stages"
while [ $stage -lt $nfcls ]; do
  FCL="Stage$stage.fcl"
 
  # In case no input files were specified, and we are not getting input
  # from sam (i.e. mc generation), recalculate the first event number,
  # the subrun number, and the number of events to generate in this worker.
  # This also applies to the textfile inputmode.
  # Note this only applies to the first stage by definition
 
  if [ $stage -eq 0 -a $USE_SAM -eq 0 ] && [ $NFILE_TOTAL -eq 0 -o "$INMODE" = 'textfile' ]; then #need to ask what is going on here


    # Don't allow --nskip.

    if [ $NSKIP -gt 0 ]; then
      echo "Illegal option --nskip specified with no input."
      exit 1
    fi

    # Do calculation.

    NSKIP=$(( $PROCESS * $NEVT / $NJOBS ))
    NEV=$(( ( $PROCESS + 1 ) * $NEVT / $NJOBS - $NSKIP ))
    FIRST_EVENT=$(( $NSKIP + 1 ))
    NSKIP=0
    NEVT=$NEV

    # Set subrun=$PROCESS+1 in a wrapper fcl file.

    SUBRUN=$(( $PROCESS + 1))
    cat <<EOF > subrun_wrapper.fcl
#include "$FCL"  

source.firstSubRun: $SUBRUN

EOF
    if [ "$INMODE" = 'textfile' ]; then
    
      if [ $NFILE_LOCAL -ne 1 ]; then
        echo "Text file input mode specified with wrong number of input files."
        exit 1
      fi
      echo "physics.producers.generator.InputFileName: \"`cat condor_lar_input.list`\"" >> subrun_wrapper.fcl
    fi
  
    FCL=subrun_wrapper.fcl
  
    echo "First MC event: $FIRST_EVENT"
    echo "MC subrun: $SUBRUN"
    echo "Number of MC events: $NEVT"
 
  fi

  # Sam stuff for main input.

  PURL=''
  CPID=''
  if [ $USE_SAM -ne 0 -a $stage -eq 0 ]; then
    echo "In SAM if" 

    # Make sure a project name has been specified.

    if [ x$SAM_PROJECT = x ]; then
      echo "No sam project was specified."
      exit 1
    fi
    echo "Sam project: $SAM_PROJECT"

    # Start project (if requested).

    if [ $SAM_START -ne 0 ]; then
      if [ x$SAM_DEFNAME != x ]; then

        echo "Starting project $SAM_PROJECT using sam dataset definition $SAM_DEFNAME"
        ifdh startProject $SAM_PROJECT $SAM_STATION $SAM_DEFNAME $SAM_USER $SAM_GROUP
        if [ $? -eq 0 ]; then
          echo "Start project succeeded."
        else
          echo "Start projet failed."
          exit 1
        fi
      fi

      if [ x$SAM_DEFNAME = x ]; then

        echo "Start project requested, but no definition was specified."
        exit 1
      fi

    fi


    # Get the project url of a running project (maybe the one we just started,
    # or maybe started externally).  This command has to succeed, or we can't
    # continue.

    PURL=`ifdh findProject $SAM_PROJECT $SAM_STATION`
    if [ x$PURL = x ]; then
      echo "Unable to find url for project ${SAM_PROJECT}."
      exit 1
    else
      echo "Project url: $PURL"
    fi

    # Start the consumer process.  This command also has to succeed.

    NODE=`hostname`
    APPFAMILY=art

    # Parse fcl file to extract process_name, and use that
    # as the application name for starting the consumer process.

    APPNAME=`fhicl-dump $FCL | grep process_name: | tr -d '"' | awk '{print $2}'`
    if [ $? -ne 0 ]; then
      echo "fhicl-dump $FCL failed to run. May be missing a ups product, library, or fcl file."
      exit 1
    fi
    if [ x$APPNAME = x ]; then
      echo "Trouble determining application name."
      echo "cat $FCL"
      cat $FCL
      exit 1
    fi

    echo "Starting consumer process."
    echo "ifdh establishProcess $PURL $APPNAME $REL $NODE $SAM_USER $APPFAMILY $FCL $NFILE $SAM_SCHEMA"
    CPID=`ifdh establishProcess $PURL $APPNAME $REL $NODE $SAM_USER $APPFAMILY $FCL $NFILE $SAM_SCHEMA`
    if [ x$CPID = x ]; then
      echo "Unable to start consumer process for project url ${PURL}."
      exit 1
    else
      echo "Consumer process id $CPID"
    fi

    # Stash away the project name and consumer process id in case we need them
    # later for bookkeeping.

    echo $SAM_PROJECT > sam_project.txt
    echo $CPID > cpid.txt

  fi
 
  # Sam stuff for secondary input.

  if [ $MIX_SAM -ne 0 ]; then
    echo "In Mix SAM if" 

    # Make sure a project name has been specified.

    if [ x$MIX_PROJECT = x ]; then
      echo "No mix sam project was specified."
      exit 1
    fi
    echo "Mix project: $MIX_PROJECT"

    # Start mix project (if requested).

    if [ $SAM_START -ne 0 ]; then
      if [ x$MIX_DEFNAME != x ]; then

        echo "Starting project $MIX_PROJECT using sam dataset definition $MIX_DEFNAME"
        ifdh startProject $MIX_PROJECT $SAM_STATION $MIX_DEFNAME $SAM_USER $SAM_GROUP
        if [ $? -eq 0 ]; then
          echo "Start project succeeded."
        else
          echo "Start projet failed."
          exit 1
        fi
      fi

      if [ x$MIX_DEFNAME = x ]; then

        echo "Start project requested, but no mix definition was specified."
        exit 1
      fi
    fi
  fi

  #Figure out output file names.
  #If outfile is not defined and we are inputing a single file or file list, follow our 
  #convention that the output file should be %inputfilename_%systemtime_stage.root

  # Construct options for lar command line.

  LAROPT="-c $FCL --rethrow-default"
  echo "Laropt: $LAROPT"
  if [ -f condor_lar_input.list -a $stage -eq 0 ]; then
    if [ "$INMODE" != 'textfile' ]; then
      # EC: Blithely make this a -s instead of -S, cuz I happen to know I am submitting a list, and a Source Module (swizzler-esque)
      # can parse this list to actually read the contained binary files. Whereas, default condor_lar tries to force usage of -S when
      # it sees that it's a text list file, which will cause Source to choke. Source inputs can not deal with the -S without serious hoop-jumping.
      # See the unfathomable instructions (unfathomable by Eric) at https://cdcvs.fnal.gov/redmine/projects/art/wiki/Writing_your_own_input_sources.
        LAROPT="$LAROPT -s condor_lar_input.list" #artroot files to read in
      #AOUTFILE=`cat condor_lar_input.list`


    fi
  fi
 
  if [ x$OUTFILE != x ]; then
    LAROPT="$LAROPT -o `basename $OUTFILE .root`$stage.root"
    outstem=`basename $OUTFILE .root` 
  fi

  if [ x$TFILE != x ]; then
    LAROPT="$LAROPT -T $TFILE"
  fi

  if [ $NEVT -ne 0 ]; then
    LAROPT="$LAROPT -n $NEVT"  
  fi

  if [ $NSKIP -ne 0 ]; then
    LAROPT="$LAROPT --nskip $NSKIP"
  fi

  if [ $FIRST_EVENT -ne 0 ]; then
    LAROPT="$LAROPT -e $FIRST_EVENT"
  fi

  if [ x$PURL != x -a $stage -eq 0 ]; then
    LAROPT="$LAROPT --sam-web-uri $PURL"
  fi

  if [ x$CPID != x -a $stage -eq 0 ]; then
    LAROPT="$LAROPT --sam-process-id $CPID"
  fi

  if [ -n "$ARGS" ]; then
    LAROPT="$LAROPT $ARGS"  
  fi
 
  if [ $stage -ne 0 ]; then
    LAROPT="$LAROPT -s $next_stage_input"
  fi 

  # Run/source optional initialization scripts.

  if [ x$INITSCRIPT != x ]; then
    echo "Running initialization script ${INITSCRIPT}."
    if ! ./${INITSCRIPT}; then
      exit $?
    fi
  fi

  if [ x$INITSOURCE != x -a $stage -eq 0 ]; then
    echo "Sourcing initialization source script ${INITSOURCE}."
    . $INITSOURCE
    status=$?
    if [ $status -ne 0 ]; then
      exit $status
    fi
  fi

  # Save a copy of the environment, which can be helpful for debugging.

  env > env.txt

  # Save a canonicalized version of the fcl configuration.

  fhicl-dump $FCL > cfgStage$stage.fcl

  # Run lar.
  pwd
  echo "$EXE $LAROPT"
  echo "$EXE $LAROPT" > commandStage$stage.txt
  $EXE $LAROPT > larStage$stage.out 2> larStage$stage.err
  stat=$?
  echo $stat > larStage$stage.stat
  echo "$EXE completed with exit status ${stat}."
  #If lar returns a status other than 0, do not move on to other stages
  if [ $stat -ne 0 ]; then
    break
  fi
 
  if [ $stage -ne 0 ]; then
   rm -rf $next_stage_input
  fi
 
  #echo `ls -t1 *.root | head -n2 | grep -v 'hist*'`
 
  #echo "Outfile is $OUTFILE"
   

  next_stage_input=`ls -t1 *.root | head -n2 | grep -v 'hist*'`

  mixed_file=`sam_metadata_dumper $next_stage_input | grep mixparent | awk -F ":" '{gsub("\"" ,""); gsub(",",""); gsub(" ",""); print $2}'`
 
  if [ x$mixed_file != x ]; then
    aunt_files=("${aunt_files[@]}" $mixed_file)
  fi

  stage=$[$stage +1]
  FIRST_EVENT=0 #I don't think this does anything
 
  #rename the mem and time profile DBs by stage
  
  if [ -f time.db ]; then
    mv time.db time$stage.db
  fi
  if [ -f mem.db ]; then
    mv mem.db mem$stage.db
  fi

done



# Setup up current version of ifdhc (may be different than version setup by larsoft).

#echo "Setting up current version of ifdhc."
#if [ x$IFDHC_DIR != x ]; then
#  unsetup ifdhc
#fi
#setup ifdhc v1_3_2
echo "IFDHC_DIR=$IFDHC_DIR"

# Sam cleanups.

if [ $USE_SAM -ne 0 ]; then

  # Get list of consumed files.

  if [ x$CPID = x -a -f cpid.txt ]; then
    CPID=`cat cpid.txt`
  fi
  ifdh translateConstraints "consumer_process_id $CPID and consumed_status consumed" > consumed_files.list

  # End consumer process.

  ifdh endProcess $PURL $CPID

  # Stop project (if appropriate).

  if [ $SAM_START -ne 0 ]; then
    echo "Stopping project."
    ifdh endProject $PURL
  fi
fi

# Secondary sam cleanups.

if [ $MIX_SAM -ne 0 ]; then

  # Stop project (if appropriate).

  if [ $SAM_START -ne 0 ]; then
    echo "Stopping project."
    MURL=`ifdh findProject $MIX_PROJECT $SAM_STATION`
    ifdh endProject $MURL
  fi
fi

# Delete input files.

if [ $USE_SAM -eq 0 -a -f condor_lar_input.list ]; then
  while read file; do
    rm -f $file
  done < condor_lar_input.list
fi

# Run optional end-of-job script.

if [ x$ENDSCRIPT != x ]; then
  echo "Running end-of-job script ${ENDSCRIPT}."
  if ! ./${ENDSCRIPT}; then
    exit $?
  fi
fi

# Do root file checks.

# Randomize names of root files that have a corresponding json file.
# These are normally histogram files.  Art files do not have external
# json metadata at this point.

# Also randomize the names of root files if there is no input specified
# for this job (i.e. generator jobs).

# Also randomize and shorten names of root files that are longer than
# 200 characters.

ran=0
if [ $USE_SAM -eq 0 -a x$INFILE = x -a x$INLIST = x ]; then
  ran=1
fi

for root in *.root; do
  nc=`echo $root | wc -c`
  if [ -f ${root}.json -o $ran != 0 -o $nc -ge 200 ]; then
    echo "Move file 1 $root"
    base=`basename $root .root | cut -c1-150`_`uuidgen`
    echo "Move file 2 $base"
    mv $root ${base}.root
    if [ -f ${root}.json ]; then
      mv ${root}.json ${base}.root.json
    fi
  fi
done

# Calculate root metadata for all root files and save as json file.
# If json metadata already exists, merge with newly geneated root metadata.
# Extract a subrun number, if one exists.  Make remote (not necessarily unique) 
# and local directories for root files with identifiable subrun numbers.

subrun=''
declare -a outdirs
declare -a logdirs
declare -a subruns
for root in *.root; do
  json=${root}.json
  if [ -f $json ]; then
    root_metadata.py $root > ${json}2
    merge_json.py $json ${json}2 > ${json}3
    mv -f ${json}3 $json
    rm ${json}2
  else
    root_metadata.py $root > $json
  fi
  subrun=`subruns.py $root | awk 'NR==1{print $2}'`
  if [ x$subrun != x ]; then
    subruns[$subrun]=$subrun
    outdirs[$subrun]=`echo $OUTDIR | sed "s/@s/$subrun/"`
    echo "Output directory for subrun $subrun is ${outdirs[$subrun]}"
    mkdir out$subrun
    logdirs[$subrun]=`echo $LOGDIR | sed "s/@s/$subrun/"`    
    echo "Log directory for subrun $subrun is ${logdirs[$subrun]}"
    mkdir log$subrun
  fi
  
  
  
done

#create a master lar.stat file which contains the overall exit code of all stages
stageStat=0
overallStat=0
while [ $stageStat -lt $nfcls ]; do
  stat=`cat larStage$stageStat.stat`
  if [[ "$stat" = 65 && $ART_VERSION < v2_01 ]]; then
   # Workaround TimeTracker crash bug for input files with zero events. 
    for json in *.json; do  
        if grep -q '"events": *"0"' $json; then
         stat=0
        fi
    done
  fi    
  overallStat=$[$stat+$overallStat]
  
  #do some cleanup of intermediate files
  #rm Stage$stageStat.fcl  
  stageStat=$[$stageStat +1] 
done
echo $overallStat > lar.stat
valstat=0

# Make local output directories for files that don't have a subrun.

mkdir out
mkdir log

# Make local files group write, if appropriate.

if [ $GRID -eq 0 -a $OUTUSER != $CURUSER ]; then
  chmod -R g+rw .
fi



# Stash all of the files we want to save in a local
# directories with a unique name.  Then copy these directories
# and their contents recursively.

# First move .root and corresponding .json files into one subdirectory.
# Note that .root files never get replicated.

for root in *.root; do
  subrun=`subruns.py $root | awk 'NR==1{print $2}'`
  mv $root out$subrun
  mv ${root}.json log$subrun
done

# Copy any remaining files into all log subdirectories.
# These small files get replicated.

for outfile in *; do
  if [ -f $outfile ]; then
    cp $outfile log
    for subrun in ${subruns[*]}
    do
      cp $outfile log$subrun
      rm -f log/$outfile
    done
  fi
done

# Do validation (if requested).

if [ $VALIDATE_IN_JOB -eq 1 ]; then
    #If SAM was used, get the parent files based on the cpid
    if [ $USE_SAM -ne 0 ]; then
      id=`cat cpid.txt`
      parent_files=($(samweb list-files consumer_process_id=$id and consumed_status consumed))
    fi
    
    echo "The file's parents are: "
    
    for elt in ${parent_files[*]};
    do
      echo $elt
    done 
    
    echo "The file's aunts are: "
    for elt in ${aunt_files[*]};
    do
      echo $elt
    done  
    
    #if we are maintain the output's parentage, combine the file's parents and aunts into a flat string
    #this string will be interpretted by validate_in_job.py. If these are left empty, then validate_in_job will not change the file's parentage
    if [ $MAINTAIN_PARENTAGE -eq 1 ]; then
       export JOBS_PARENTS=`echo ${parent_files[*]}`
       export JOBS_AUNTS=`echo ${aunt_files[*]}`
    
    
    fi
    
    # Do validation function for the whole job.

    valstat=0
    curdir=`pwd`
    #cd $curdir/log
    #validate_in_job.py --dir $curdir/out --logfiledir $curdir/log --outdir $OUTDIR/$OUTPUT_SUBDIR --declare $DECLARE_IN_JOB --copy $COPY_TO_FTS
    #valstat=$?
    #cd $curdir

    # Do validation for each subrun.

    for subrun in ${subruns[*]}
    do
      cd $curdir/log$subrun
      
      validate_in_job.py --dir $curdir/out$subrun --logfiledir $curdir/log$subrun --outdir ${outdirs[$subrun]}/$OUTPUT_SUBDIR --declare $DECLARE_IN_JOB --copy $COPY_TO_FTS
      subvalstat=$?
      valstat=$(( $valstat + $subvalstat ))
    done
    cd $curdir
    
fi

# Make a tarball of each log directory, and save the tarball in its own log directory.

rm -f log0.tar
tar -cjf log0.tar -C log .
mv log0.tar log
for subrun in ${subruns[*]}
do
  rm -f log.tar
  tar -cjf log.tar -C log$subrun .
  mv log.tar log$subrun
done

# Clean remote output and log directories.

if [  1 -eq 0 ]; then
  export IFDH_FORCE=$FORCE #this isn't set when running interactive, causing problems...
fi

for dir in ${LOGDIR} ${OUTDIR}
do
  #echo $dir
  echo "Make sure directory0 ${dir}/$OUTPUT_SUBDIR exists."
  
  #mkdir ${dir}/$OUTPUT_SUBDIR 
  mkdir.py -v ${dir}/$OUTPUT_SUBDIR
  echo "Make sure directory0 ${dir}/$OUTPUT_SUBDIR is empty."
  emptydir.py -v ${dir}/$OUTPUT_SUBDIR
  mkdir.py -v ${dir}/$OUTPUT_SUBDIR
  echo "Directory0 ${dir}/$OUTPUT_SUBDIR clean ok."
done

if [ $SINGLE != 0 ]; then
  for dir in ${logdirs[*]} ${outdirs[*]}
  do
    echo "Make sure directory1 $dir exists."
    mkdir.py -v $dir
    echo "Make sure directory1 $dir is empty."
    emptydir.py -v $dir
    mkdir.py -v $dir/$OUTPUT_SUBDIR
    echo "Directory1 $dir/$OUTPUT_SUBDIR clean ok."
  done
else
  for dir in ${logdirs[*]} ${outdirs[*]}
  do
    echo "Make sure directory2 ${dir}/$OUTPUT_SUBDIR exists."
    mkdir.py -v ${dir}/$OUTPUT_SUBDIR
    echo "Make sure directory2 ${dir}/$OUTPUT_SUBDIR is empty."
    emptydir.py -v ${dir}/$OUTPUT_SUBDIR
    mkdir.py -v ${dir}/$OUTPUT_SUBDIR
    echo "Directory2 ${dir}/$OUTPUT_SUBDIR clean ok."
  done
fi

statout=0
export IFDH_CP_MAXRETRIES=0
echo "ls log"
ls log
echo "ifdh cp -D $IFDH_OPT log/* ${LOGDIR}/$OUTPUT_SUBDIR"
if [ "$( ls -A log )" ]; then
  ifdh cp -D $IFDH_OPT log/* ${LOGDIR}/$OUTPUT_SUBDIR
  stat=$?
  if [ $stat -ne 0 ]; then
    echo "ifdh cp failed with status ${stat}."
  fi
fi

for subrun in ${subruns[*]}
do
  echo "ls log$subrun"
  ls log$subrun
  echo "ifdh cp -D $IFDH_OPT log${subrun}/* ${logdirs[$subrun]}/$OUTPUT_SUBDIR"
  ifdh cp -D $IFDH_OPT log${subrun}/* ${logdirs[$subrun]}/$OUTPUT_SUBDIR
  stat=$?
  if [ $stat -ne 0 ]; then
    echo "ifdh cp failed with status ${stat}."
    statout=$stat
  fi
done

if [ "$( ls -A out )" ]; then
  echo "ifdh cp -D $IFDH_OPT out/* ${OUTDIR}/$OUTPUT_SUBDIR"
  ifdh cp -D $IFDH_OPT out/* ${OUTDIR}/$OUTPUT_SUBDIR
  stat=$?
  if [ $stat -ne 0 ]; then
    echo "ifdh cp failed with status ${stat}."
  fi
fi

for subrun in ${subruns[*]}
do
  echo "ifdh cp -D $IFDH_OPT out${subrun}/* ${outdirs[$subrun]}/$OUTPUT_SUBDIR"
  ifdh cp -D $IFDH_OPT out${subrun}/* ${outdirs[$subrun]}/$OUTPUT_SUBDIR
  stat=$?
  if [ $stat -ne 0 ]; then
    echo "ifdh cp failed with status ${stat}."
  fi
    statout=$stat 

done   

if [ $statout -eq 0 ]; then
  statout=`cat lar.stat`
fi

if [ $statout -eq 0 ]; then
  statout=$valstat
fi  

exit $statout
