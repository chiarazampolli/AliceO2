#!/bin/bash

SEVERITY="detail"
ENABLE_METRICS=1
source $O2DPG_ROOT/DATA/common/setenv.sh
source $O2_ROOT/prodtests/full-system-test/workflow-setup.sh

# ---------------------------------------------------------------------------------------------------------------------
# Set general arguments
source $O2DPG_ROOT/DATA/common/getCommonArgs.sh

# Set up calibrations
source $O2DPG_ROOT/DATA/common/setenv_calib.sh

# check that WORKFLOW_DETECTORS is needed, otherwise the wrong calib wf will be built
if [[ -z $WORKFLOW_DETECTORS ]]; then echo "WORKFLOW_DETECTORS must be defined" 1>&2; exit 1; fi

# CCDB destination for uploads
[[ -z ${CCDBPATH+x} ]] && CCDBPATH="http://o2-ccdb.internal"
echo "CCDBPATH = $CCDBPATH" 1>&2

# Adding calibrations
EXTRA_WORKFLOW_CALIB=

echo "CALIB_PRIMVTX_MEANVTX = $CALIB_PRIMVTX_MEANVTX" 1>&2
echo "CALIB_TOF_LHCPHASE = $CALIB_TOF_LHCPHASE" 1>&2
echo "CALIB_TOF_CHANNELOFFSETS = $CALIB_TOF_CHANNELOFFSETS" 1>&2
echo "CALIB_TOF_DIAGNOSTICS = $CALIB_TOF_DIAGNOSTICS" 1>&2
echo "CALIB_EMC_CHANNELCALIB = $CALIB_EMC_CHANNELCALIB" 1>&2
echo "CALIB_PHS_ENERGYCALIB = $CALIB_PHS_ENERGYCALIB" 1>&2
echo "CALIB_PHS_BADMAPCALIB = $CALIB_PHS_BADMAPCALIB" 1>&2
echo "CALIB_PHS_TURNONCALIB = $CALIB_PHS_TURNONCALIB" 1>&2
echo "CALIB_PHS_RUNBYRUNCALIB = $CALIB_PHS_RUNBYRUNCALIB" 1>&2
echo "CALIB_TRD_VDRIFTEXB = $CALIB_TRD_VDRIFTEXB" 1>&2

# PrimVertex
if [[ $CALIB_PRIMVTX_MEANVTX == 1 ]]; then
    EXTRA_WORKFLOW_CALIB+="o2-calibration-mean-vertex-calibration-workflow $ARGS_ALL | "
fi

# TOF
if [[ $CALIB_TOF_LHCPHASE == 1 ]] || [[ $CALIB_TOF_CHANNELOFFSETS == 1 ]]; then
    if [[ $CALIB_TOF_LHCPHASE == 1 ]]; then
  EXTRA_WORKFLOW_CALIB+="o2-calibration-tof-calib-workflow $ARGS_ALL --do-lhc-phase --tf-per-slot 10 | "
    fi
    if [[ $CALIB_TOF_CHANNELOFFSETS == 1 ]]; then
  EXTRA_WORKFLOW_CALIB+="o2-calibration-tof-calib-workflow $ARGS_ALL --do-channel-offset --update-at-end-of-run-only --min-entries 8 --range 100000 | "
    fi
fi
if [[ $CALIB_TOF_DIAGNOSTICS == 1 ]]; then
    EXTRA_WORKFLOW_CALIB+="o2-calibration-tof-diagnostic-workflow $ARGS_ALL --tf-per-slot 26400 | "
fi

# TRD
if [[ $CALIB_TRD_VDRIFTEXB == 1 ]]; then
    EXTRA_WORKFLOW_CALIB+="o2-calibration-trd-vdrift-exb $ARGS_ALL | "
fi

# Calo cal
# EMC
if [[ $CALIB_EMC_CHANNELCALIB == 1 ]]; then
    EXTRA_WORKFLOW_CALIB+="o2-calibration-emcal-channel-calib-workflow --calibMode timeCalib $ARGS_ALL | "
fi

# PHS
if [[ $CALIB_PHS_ENERGYCALIB == 1 ]]; then
    EXTRA_WORKFLOW_CALIB+="o2-phos-calib-workflow --energy $ARGS_ALL | "
fi
if [[ $CALIB_PHS_BADMAPCALIB == 1 ]]; then
    EXTRA_WORKFLOW_CALIB+="o2-phos-calib-workflow --badmap --mode 0 $ARGS_ALL | "
fi
if [[ $CALIB_PHS_TURNONCALIB == 1 ]]; then
    EXTRA_WORKFLOW_CALIB+="o2-phos-calib-workflow --turnon $ARGS_ALL | "
fi
if [[ $CALIB_PHS_RUNBYRUNCALIB == 1 ]]; then
    EXTRA_WORKFLOW_CALIB+="o2-phos-calib-workflow --runbyrun $ARGS_ALL | "
fi

# starting with empty workflow
if workflow_has_parameters CALIB_PROXIES; then
    WORKFLOW=
    if [[ ! -z $CALIBDATASPEC_BARREL ]]; then
  WORKFLOW+="o2-dpl-raw-proxy ${ARGS_ALL} --dataspec \"$CALIBDATASPEC_BARREL\" $(get_proxy_connection barrel input) | "
    fi
    if [[ ! -z $CALIBDATASPEC_CALO ]]; then
  WORKFLOW+="o2-dpl-raw-proxy ${ARGS_ALL} --dataspec \"$CALIBDATASPEC_CALO\" $(get_proxy_connection calo input) | "
    fi
fi

WORKFLOW+=$EXTRA_WORKFLOW_CALIB

if [[ $CCDBPATH != "none" ]]; then WORKFLOW+="o2-calibration-ccdb-populator-workflow --ccdb-path $CCDBPATH $ARGS_ALL | "; fi

if workflow_has_parameters CALIB_PROXIES; then
    WORKFLOW+="o2-dpl-run $ARGS_ALL $GLOBALDPLOPT -b"
    if [ $WORKFLOWMODE == "print" ]; then
  echo Workflow command adding aggregator:
  echo $WORKFLOW | sed "s/| */|\n/g"
    else
  # Execute the command we have assembled
  WORKFLOW+=" --$WORKFLOWMODE"
  eval $WORKFLOW
    fi
fi
