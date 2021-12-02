
#path to reco file
ONIAFILEPATH="../rootfiles/datasets/trigger/OniaTree_Run3_2021_PG_0p5_3_merged.root"
#path to hltobj triggers file
TRIGGERFILEPATH="../rootfiles/datasets/trigger/HLT/HLT_MC_ParticleGun_Muon_Pt_3_100_1130pre5_v1.root"
#path to directory to place output
OUTPUTPATH="../rootfiles/analysis/triggerStudy"
#reco file is low pt or high pt : "lowpt" or "highpt"
PTRANGE="highpt"
#type of onia for dimuons (for acceptance selection) : "JPsi" or "Y"
ONIATYPE="JPsi"
#maximum number of cores/threads to use simultaneusly, one core/thread process one trigger
MAXJOBS=7

#list of triggers to process, listed in file triggersLx.txt
TRIGGERSL1=$( cat triggersL1.txt )
TRIGGERSL2=$( cat triggersL2.txt )
TRIGGERSL3=$( cat triggersL3.txt )
TRIGGERS=( ${TRIGGERSL1[@]} ${TRIGGERSL2[@]} ${TRIGGERSL3[@]} )

OUTPUTDIR="${OUTPUTPATH}/$( basename ${ONIAFILEPATH%.*})"

mkdir -p ${OUTPUTPATH}

echo "Full HLT Trigger efficiency study"

#do the processing

echo "processing..."
echo "reading reco file '${ONIAFILEPATH}'"
echo "reading hltobj file '${TRIGGERFILEPATH}'"
echo "output to:"
echo "  ${OUTPUTDIR}"

for TRIGGERNAME in ${TRIGGERS[@]}
do
    OUTPATH="${OUTPUTDIR}/${TRIGGERNAME}"
    mkdir -p $OUTPATH
    echo "processing ${TRIGGERNAME}"
    ./trigeff.sh ${ONIAFILEPATH} ${TRIGGERFILEPATH} ${TRIGGERNAME} $OUTPATH $PTRANGE $ONIATYPE &> "${OUTPATH}/output.log"  &
    JOBS=( $(jobs -p) )
    JOBNUM="${#JOBS[@]}"
    if [ $MAXJOBS = $JOBNUM ]
    then
        wait -n
    fi
done
wait

./PlotEff/ploteff "-multi" "${OUTPUTDIR}" $TRIGGERSL1 "${OUTPUTDIR}/triggersL1"
./PlotEff/ploteff "-multi" "${OUTPUTDIR}" $TRIGGERSL2 "${OUTPUTDIR}/triggersL2"
./PlotEff/ploteff "-multi" "${OUTPUTDIR}" $TRIGGERSL3 "${OUTPUTDIR}/triggersL3"

mkdir "${OUTPUTDIR}/report"
cp "HelperScripts/exampleReport.tex" "${OUTPUTDIR}/report/presentation.tex"
cd "${OUTPUTDIR}/report"
pdflatex -interaction batchmode "presentation.tex"
echo "all done"
