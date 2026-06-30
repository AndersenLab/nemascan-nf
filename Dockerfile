##############################################################################################################################
#
# This container includes all necessary components for initializing the NemaScan-nf Nextflow pipeline in Google Cloud
# Additional configuration options can be passed in via environment variables
#
##############################################################################################################################

# Base image includes Google Cloud SDK tools
FROM google/cloud-sdk:slim

# Install OpenJDK JRE for Nextflow
RUN apt-get update && apt-get install -y --no-install-recommends openjdk-21-jre wget procps

LABEL Name="NemaScan-NF"

# Specify Nextflow version and mode 
ENV NXF_VER=26.04.1 \
  NXF_MODE=google \
  NXF_EDGE=0

# RUN git clone https://github.com/AndersenLab/NemaScan.git /nemascan && \
#     cd /nemascan && \
#     git checkout 2f19f5f80dcc397d73698fd5cd3cb571c53299b6

WORKDIR /nemascan

COPY assets/ /nemascan/assets/
COPY modules/ /nemascan/modules/
COPY subworkflows/ /nemascan/subworkflows/
COPY workflows/ /nemascan/workflows/
COPY conf/ /nemascan/conf/
COPY main.nf /nemascan/
COPY nextflow.config /nemascan/
COPY nemascan-nxf.sh /nemascan/


# Run the Nextflow install script (version and mode must be piped in to bash during install 
# or nextflow will initially download the latest version and only download and switch to NXF_VER when the container runs)
RUN wget -qO- https://get.nextflow.io | bash

# add nextflow and nemarun directory to te system path and make them executable
ENV PATH="/nemascan:${PATH}"
RUN chmod +x /nemascan/nemascan-nxf.sh /nemascan/nextflow

# cache nf-google plugin
RUN mkdir /nemascan/plugins
ENV NXF_PLUGINS_DIR=/nemascan/plugins
RUN nextflow plugin install nf-google@1.27.2

# Link executable scripts to bin directory
RUN mkdir /nemascan/bin && \
  for D in $(ls /nemascan/modules/local/); do \
    if [ -e /nemascan/modules/local/${D}/resources/usr/bin ]; then \
      for SCRIPT in $(ls /nemascan/modules/local/${D}/resources/usr/bin/); do \
        if [ $(echo "$SCRIPT" | grep -E "\.(py|R|r|sh)$") ]; then \
          ln -sf "../modules/local/${D}/resources/usr/bin/${SCRIPT}" /nemascan/bin/; \
        fi \
      done \
    fi \
  done && \
  if [ "$(ls -A /nemascan/bin)" ]; then chmod 777 /nemascan/bin/*; fi

