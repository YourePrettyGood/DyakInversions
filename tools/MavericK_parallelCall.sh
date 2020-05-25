#!/bin/bash

JOBDIR=$1
PARAMS=$2
PREFIX=$3
TASKID=$4

TASKPREFIX="${PREFIX}_window${TASKID}"

if [[ ! "${JOBDIR}" =~ /\/$/ ]]; then
   JOBDIR="${JOBDIR}/"
   echo "Added terminal slash to job directory"
fi

if [[ ! -e "${JOBDIR}${TASKPREFIX}_forSTRUCTURE.tsv" ]]; then
   echo "Unable to find input data ${JOBDIR}${TASKPREFIX}_forSTRUCTURE.tsv"
   exit 2
fi

if [[ ! -e "${JOBDIR}${PARAMS}" ]]; then
   echo "Unable to find MavericK parameters file ${JOBDIR}${PARAMS}"
   exit 3
fi

mkdir -p logs
mkdir -p Evidence
mkdir -p Likelihoods
mkdir -p QMatrices

/usr/bin/time -v MavericK -parameters ${PARAMS} -masterRoot ${JOBDIR} -data ${TASKPREFIX}_forSTRUCTURE.tsv -EMalgorithm_on 1 -EMiterations 100 -EMrepeats 10 -outputEvidence_on 1 -outputEvidence Evidence/${TASKPREFIX}_Evidence.csv -outputEvidenceNormalised_on 1 -outputEvidenceNormalised Evidence/${TASKPREFIX}_NormalizedEvidence.csv -outputLikelihood_on 1 -outputLikelihood Likelihoods/${TASKPREFIX}_Likelihood.csv -outputLog_on 1 -outputLog logs/${TASKPREFIX}_Log.txt -outputQmatrix_ind_on 1 -outputQmatrix_ind QMatrices/${TASKPREFIX}_Qmatrix_perInd.csv
