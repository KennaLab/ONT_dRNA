################## BASE IMAGE ######################
FROM python:3.13

################## METADATA ######################
LABEL base_image="python:3.13"
LABEL version="0.5.0"
LABEL extra.binaries="blue-crab"

################## INSTALLATIONS ######################
# Note: if not using the bin folder, you might want to add dir to environment variable PATH.
WORKDIR /usr/local/bin

# Install the application dependencies
RUN pip install --no-cache-dir blue-crab==0.5.0

ENTRYPOINT ["python"]
