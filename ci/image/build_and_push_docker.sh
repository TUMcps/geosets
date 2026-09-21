#! /bin/bash

# This shell script is used to automatically build the containers for the python image
# ------------------------------------------------------------------------------------

# Information
echo "Starting to build and push the containers for geosets. This may take some time."
echo "You have to be in the docker group to run this script successfully."
echo -e "\n  -- lrz gitlab login --"

# Docker gitlab login (You need to use access token for password)
docker login gitlab.lrz.de:5005

# Build and push cpp image (gcc + all C++ deps)
# docker build --no-cache -f Dockerfile -t gitlab.lrz.de:5005/cps/geosets/gcc-all-dependencies .
# docker push gitlab.lrz.de:5005/cps/geosets/gcc-all-dependencies

# Build and push python image (uv + C++ deps for scikit-build)
docker build -f Dockerfile.python -t gitlab.lrz.de:5005/cps/geosets/python-all-dependencies .
docker push gitlab.lrz.de:5005/cps/geosets/python-all-dependencies

# docker gitlab logout
docker logout gitlab.lrz.de:5005
