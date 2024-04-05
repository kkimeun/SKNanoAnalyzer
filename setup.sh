echo -e "\033[31m##################### WARNING ########################\033[0m"
echo -e "\033[31m####         THIS IS DEVELOPMENT VERSION          ####\033[0m"
echo -e "\033[31m######################################################\033[0m"
echo ""

# Set up environment
echo "@@@@ Working Directory: `pwd`"
export SKNANO_HOME=`pwd`
export SKNANO_RUNLOG="/home/$USER/workspace/SKNanoRunlog"
export SKNANO_OUTPUT="/home/$USER/workspace/SKNanoOutput"

# root configuration
RELEASE="`cat /etc/redhat-release`"
echo "@@@@ Configuring ROOT for $RELEASE"
if [[ $RELEASE == *"7."* ]]; then
    echo "@@@@ Running on $RELEASE"
    source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-centos7-gcc12-opt/setup.sh
elif [[ $RELEASE == *"8."* ]]; then
    echo "@@@@ Running on $RELEASE"
    source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-centos8-gcc12-opt/setup.sh
elif [[ $RELEASE == *"9."* ]]; then
    echo "@@@@ Running on $RELEASE"
    source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh
else
    echo "@@@@ Not running on redhat 7, 8, or 9"
    echo "@@@@ Assuming root is installed in conda environment 'nano'"
    source ~/.conda-activate
    conda activate nano
fi
echo "@@@@ ROOT path: $ROOTSYS"


export SKNANO_LIB=$SKNANO_HOME/lib
export SKNANO_VERSION="Run3UltraLegacy_v1"
export SKNANO_DATA=$SKNANO_HOME/data/$SKNANO_VERSION
mkdir -p $SKNANO_DATA

export SKNANO_BIN=$SKNANO_HOME/bin
export SKNANO_PYTHON=$SKNANO_HOME/python
export PATH=$SKNANO_PYTHON:$PATH
export PYTHONPATH=$PYTHONPATH:$SKNANO_PYTHON

export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$SKNANO_HOME/DataFormats/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SKNANO_LIB
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$SKNANO_LIB

# setting LHAPDFs
if [ -d "external/lhapdf" ]; then
    export PATH=$PATH:$SKNANO_HOME/external/lhapdf/bin
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SKNANO_HOME/external/lhapdf/lib
    export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$SKNANO_HOME/external/lhapdf/lib
    export LHAPDF_DATA_PATH=$SKNANO_HOME/external/lhapdf/data
elif [ -d "/cvmfs" ]; then
    echo ""@@@@ configuring LHAPDF from cvmfs
else
    echo "@@@@ LHAPDF not found"
    echo "@@@@ consider to install LHAPDF using ./scripts/install_lhapdf.sh"
    exit(1)
fi
export LHAPDF_INCLUDE_DIR=`lhapdf-config --incdir`
export LHAPDF_LIB_DIR=`lhapdf-config --libdir`

echo "@@@@ LHAPDF include: $LHAPDF_INCLUDE_DIR"
echo "@@@@ LHAPDF lib: $LHAPDF_LIB_DIR"
echo "@@@@ reading data from $LHAPDF_DATA_PATH"

# env for correctionlibs
export CORRECTION_INCLUDE_DIR=`correction config --incdir`
export CORRECTION_LIB_DIR=`correction config --libdir`
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CORRECTION_LIB_DIR
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$CORRECTION_LIB_DIR
export JSONPOG_INTEGRATION_PATH=$SKNANO_HOME/external/jsonpog-integration

echo "@@@@ Correction include: $CORRECTION_INCLUDE_DIR"
echo "@@@@ Correction lib: $CORRECTION_LIBS"
