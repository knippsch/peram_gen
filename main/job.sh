# @ job_name         = A40_test
# @ error            = $(job_name).$(jobid).$(stepid).out
# @ output           = $(job_name).$(jobid).$(stepid).out
# @ environment      = COPY_ALL;
# @ wall_clock_limit = 01:00:00
# @ notification     = always
# @ notify_user      = b.knippschild@gmx.de
# @ job_type         = bluegene
# @ bg_connectivity  = TORUS
# @ bg_size          = 128
# @ queue

export NP=$LOADL_BG_SIZE
export NT=64
export NPN=1
export OMP_NUM_THREADS=${NT}

export NAME=${LOADL_JOB_NAME}

# the executable
export EFILE=${HOME}/peram_gen/main/main
# the input file to copy to the working directory
#export IFILE=/peram_gen/main/invert.input

#export RUNDIR=runs/nf211/${NAME}

# working directory
export WDIR=${WORK}
# output directory
export ODIR=${WORK}

# the output file
#export OFILE=${LOADL_STEP_OUT}
#echo "std output file: ${OFILE}" 

#if [[ ! -d ${WDIR} ]]
#then
#  mkdir -p ${WDIR}
#fi

#cp invert.input ${WDIR}/invert.input

# figure out the correct mapping
echo LOADL_BG_SHAPE $LOADL_BG_SHAPE

# default mapping for BG_SIZE=512 and 32
MP=EABCDT

case ${LOADL_BG_SHAPE} in
  # mappings for BG_SIZE=1024

  2x1x1x1 )
    MP=EABCDT 
  ;;

  1x2x1x1 )
    MP=EBACDT
  ;;    

  1x1x2x1 )
    MP=ECABDT
  ;;    

  1x1x1x2 )
    MP=EDABCT
  ;;

# mappings for BG_SIZE=64 (8x2x2x2)
  2x2x2x4x2 )
    MP=EDABCT
  ;;
  2x2x4x2x2 )
    MP=ECABDT
  ;;
  2x4x2x2x2 )
    MP=EBACDT
  ;;
  4x2x2x2x2 )
    MP=EABCDT
  ;;
  # mappings for bg_size=128 (16x2x2x2)
  2x2x4x4x2 )
    MP=CDABET
  ;;
  2x4x2x4x2 )
    MP=BDACET
  ;;
  4x2x2x4x2 )
    MP=ADBCET
  ;;
  2x4x4x2x2 )
    MP=BCADET
  ;;
  4x2x4x2x2 )
    MP=ACBDET
  ;;
  4x4x2x2x2 )
    MP=ABCDET
  ;;
  # mappings for bg_size=256 (16x4x2x2)
  2x4x4x4x2 )
    MP=BCDAET
  ;;
  4x4x2x4x2 )
    MP=ABDCET
  ;;
  4x2x4x4x2 )
    MP=ACDBET
  ;;
  4x4x4x2x2 )
    MP=ABCDET
  ;;
esac

cd ${HOME}/peram_gen/main

date

runjob --mapping ${MP} --envs "MUSPI_NUMINJFIFOS=8" --envs "MUSPI_NUMRECFIFOS=8" --envs "MUSPI_NUMBATIDS=2" --np ${NP} --ranks-per-node ${NPN} --cwd ${WDIR} --exe ./main --args ""
RETVAL=$?

date

