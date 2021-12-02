
#path to reco file
ONIAFILEPATH=${1:-"../rootfiles/datasets/trigger/OniaTree_Run3_2021_PG_3_100_merged.root"}
#path to hltobj triggers file
TRIGGERFILEPATH=${2:-"../rootfiles/datasets/trigger/HLT/HLT_MC_ParticleGun_Muon_Pt_3_100_1130pre5_v1.root"}
#name of the trigger to process
TRIGGERNAME=${3:-"HLT_HIL1SingleMu0_v"}
#path to directory to place output
OUTPUTPATH=${4:-"../rootfiles/analysis/triggerStudy/HLT_MC_ParticleGun_Muon_Pt_3_100_1130pre5_v1/HLT_HIL1DoubleMuOpen_v1"}
#reco file is low pt or high pt : "lowpt" or "highpt"
PTRANGE=${5:-"highpt"}
#type of onia for dimuons (for acceptance selection) : "JPsi" or "Y"
ONIATYPE=${6:-"JPsi"}

mkdir -p ${OUTPUTPATH}

echo "Trigger efficiency study"

echo "reading reco file '${ONIAFILEPATH}'"
echo "reading hltobj file '${TRIGGERFILEPATH}'"

echo "output to:"
echo "  ${OUTPUTPATH}/$( basename $ONIAFILEPATH)"

./TrigEff/trigeff ${ONIAFILEPATH} ${TRIGGERFILEPATH} ${TRIGGERNAME} $OUTPUTPATH ${ONIATYPE}
echo "generating plots for ${TRIGGERNAME}"
./PlotEff/ploteff "${OUTPUTPATH}/output.root" "${OUTPUTPATH}" ${PTRANGE}