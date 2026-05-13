################## BASE IMAGE ######################
FROM davidsongroup/jaffa:2.4@sha256:df76cb343b340320c2e6a853e403bb33d253de746300f09eaedbbe2d53774898

################## METADATA ######################
LABEL base_image="davidsongroup/jaffa:2.4"
LABEL version="2.4"
LABEL extra.binaries="ps"

################## INSTALLATIONS ######################
RUN apt-get update && \
    apt-get install -y procps && \
    rm -rf /var/lib/apt/lists/*
