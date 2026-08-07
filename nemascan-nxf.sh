#!/bin/bash
#
# This script acts as a wrapper around the execution of Nextflow passing environment variables as arguments
#
###################################################################################################################

DEFAULT_GOOGLE_PROJECT="caendr"
DEFAULT_QUEUE_REGION="us-east1"
DEFAULT_GOOGLE_SERVICE_ACCOUNT_EMAIL="caendr-pipeline-user@caendr.iam.gserviceaccount.com"
DEFAULT_SPECIES="c_elegans"
DEFAULT_VCF_VERSION="20250625"
DEFAULT_AWS_BUCKET="https://caendr-open-access-data-bucket.s3.us-east-2.amazonaws.com"


# Environment variables with default values:

if [[ -z "${GOOGLE_PROJECT}" ]]; then
  GOOGLE_PROJECT=${DEFAULT_GOOGLE_PROJECT}
  echo "GOOGLE_PROJECT environment variable is not set - defaulting to ${GOOGLE_PROJECT}"
fi

if [[ -z "${QUEUE_REGION}" ]]; then
  QUEUE_REGION=${DEFAULT_QUEUE_REGION}
  echo "QUEUE_REGION environment variable is not set - defaulting to ${QUEUE_REGION}"
fi

if [[ -z "${GOOGLE_SERVICE_ACCOUNT_EMAIL}" ]]; then
  GOOGLE_SERVICE_ACCOUNT_EMAIL=${DEFAULT_GOOGLE_SERVICE_ACCOUNT_EMAIL}
  echo "GOOGLE_SERVICE_ACCOUNT_EMAIL environment variable is not set - defaulting to ${GOOGLE_SERVICE_ACCOUNT_EMAIL}"
fi

if [[ -z "${SPECIES}" ]]; then
  SPECIES=${DEFAULT_SPECIES}
  echo "SPECIES environment variable is not set - defaulting to ${SPECIES}"
fi

if [[ -z "${VCF_VERSION}" ]]; then
  VCF_VERSION=${DEFAULT_VCF_VERSION}
  echo "VCF_VERSION environment variable is not set - defaulting to ${VCF_VERSION}"
fi

if [[ -z "${AWS_BUCKET}" ]]; then
  AWS_BUCKET=${DEFAULT_AWS_BUCKET}
  echo "AWS_BUCKET environment variable is not set - defaulting to ${AWS_BUCKET}"
fi


# Environment variables that MUST be set

if [[ -z "${TRAIT_FILE}" ]]; then
  echo "TRAIT_FILE environment variable must be set to the Google Storage path of the data"
  exit 1
fi

if [[ -z "${OUTPUT_DIR}" ]]; then
  echo "OUTPUT_DIR environment variable must be set to the Google Storage path of the output directory"
  exit 1
fi

if [[ -z "${WORK_DIR}" ]]; then
  echo "WORK_DIR environment variable must be set to the Google Storage path of the working directory"
  exit 1
fi

gcloud storage cp ${TRAIT_FILE} data_raw.tsv
dos2unix -O data_raw.tsv > traits.tsv
gcloud storage cp traits.tsv ${OUTPUT_DIR}/data_unix.tsv

nextflow run main.nf \
  -profile gcp \
  --google_project "${GOOGLE_PROJECT}" \
  --google_zone "${QUEUE_REGION}" \
  --google_service_account_email "${GOOGLE_SERVICE_ACCOUNT_EMAIL}" \
  --traits traits.tsv \
  --species "${SPECIES}" \
  --vcf "${AWS_BUCKET}/dataset_release/${SPECIES}/${VCF_VERSION}/variation/WI.${VCF_VERSION}.small.hard-filter.isotype.vcf.gz" \
  --imputed "${AWS_BUCKET}/dataset_release/${SPECIES}/${VCF_VERSION}/variation/WI.${VCF_VERSION}.impute.isotype.vcf.gz" \
  --annotation "${AWS_BUCKET}/dataset_release/${SPECIES}/${VCF_VERSION}/annotation/WI.${VCF_VERSION}.csq.strain-annotation.csv.gz" \
  --haplotypes "${AWS_BUCKET}/dataset_release/${SPECIES}/${VCF_VERSION}/haplotype/haplotype_df_isotype.bed" \
  --genes "${AWS_BUCKET}/dataset_release/${SPECIES}/${VCF_VERSION}/browser_tracks/${VCF_VERSION}_${SPECIES}_transcripts.bed.gz" \
  --isogroups "https://caendr.org/request-strains/isotype_list/download" \
  --significance_threshold EIGEN \
  --git_info /nemascan/assets/git_info.tsv \
  --work_dir "${WORK_DIR}" \
  -output-dir "${OUTPUT_DIR}"

EXITCODE=$?

gcloud storage cp .nextflow.log ${OUTPUT_DIR}/

exit ${EXITCODE}