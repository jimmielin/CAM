#! /bin/bash

##############################################################################
###
### A stub data assimilation script that prints out information but makes
### no modifications to model data.
### Script checks for proper pre and post data assimilation output
### Tests using this script should be BFB with a non-data assimilation run
###
##############################################################################

errcode=0
if [ $# -ne 2 ]; then
  echo "ERROR: Wrong number of arguments, $# (should be 2)"
  errcode=$(( errcode + 1 ))
else
  caseroot=$1
  cycle=$2
  echo "caseroot: ${caseroot}"
  echo "cycle: ${cycle}"
  cd ${caseroot}
  res=$?
  if [ $res -ne 0 ]; then
    echo "ERROR: Unable to cd to caseroot, ${caseroot}"
    errcode=$(( errcode + 1 ))
  else
    ./xmlchange DATA_ASSIMILATION_ATM=TRUE
    res=$?
    if [ $res -ne 0 ]; then
      echo "ERROR: Unable to change DATA_ASSIMILATION_ATM to TRUE"
      errcode=$(( errcode + 1 ))
    fi
    rundir="`./xmlquery --value RUNDIR`"
    ninst=`./xmlquery --value NINST_ATM`
    if [ -n "${rundir}" -a -d "${rundir}" ]; then
      cd ${rundir}
      res=$?
      if [ $res -ne 0 ]; then
        echo "ERROR: Unable to cd to rundir, ${rundir}"
        errcode=$(( errcode + 1 ))
      else
        # Check the latest log file for a resume signal
        if [ $ninst -eq 1 ]; then
          lfiles="`ls -t atm.log.* 2> /dev/null | head -1`"
        else
          # Multi-instance, look for wav_nnnn.log*
          for inst in `seq 1 $ninst`; do
            ifilename="`printf "atm_%04d.log.*" $inst`"
            ifile="`ls -t ${ifilename} 2> /dev/null | head -1`"
            if [ -z "${ifile}" ]; then
              echo "No log files for instance $ninst found"
              errcode=$(( errcode + 1 ))
            elif [ -z "${lfiles}" ]; then
              lfiles="${ifile}"
            else
              lfiles="${lfiles} ${ifile}"
            fi
          done
        fi
        if [ -z "${lfiles}" ]; then
          echo "ERROR: Unable to find atm log file in `pwd -P`"
          errcode=$(( errcode + 1 ))
        else
          for atmfile in ${lfiles}; do
            dasig="`zgrep '^[ ]*DART run using CAM initial mode$' ${atmfile} 2> /dev/null`"
            initsig="`zgrep '^[ ]*Initial run$' ${atmfile} 2> /dev/null`"
            if [ $cycle -gt 0 ]; then
              if [ -n "${dasig}" ]; then
                echo "Post-DA resume signal found for cycle ${cycle}"
              else
                echo "No post-DA resume signal for cycle ${cycle}"
              fi
            elif [ -n "${dasig}" ]; then
              echo "Bad Post-DA resume signal found for cycle ${cycle}"
            fi
            if [ $cycle -eq 0 ]; then
              if [ -n "${initsig}" ]; then
                echo "Initial run signal found for cycle ${cycle}"
              else
                echo "No initial run signal found for cycle ${cycle}"
              fi
            elif [ -n "${initsig}" ]; then
              echo "Bad initial run signal found for cycle ${cycle}"
            fi
          done
        fi
      fi
    else
      echo "ERROR: RUNDIR (${rundir}) is not a valid directory"
      errcode=$(( errcode + 1 ))
    fi
  fi
fi

exit $errcode
