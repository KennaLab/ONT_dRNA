################## BASE IMAGE ######################
FROM quay.io/biocontainers/nanopolish:0.14.0--h5ca1c30_6

################## METADATA ######################
LABEL base_image="quay.io/biocontainers/nanopolish:0.14.0--h5ca1c30_6"
LABEL version="0.14.0--h5ca1c30_6"

################## INSTALLATIONS ######################
ENV HDF5_PLUGIN_PATH=/usr/local/hdf5/lib/plugin/
ENV PATH=$PATH:$HDF5_PLUGIN_PATH
