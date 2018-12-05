#!/usr/bin/env bash
set -e

REPO=scottx611x/sigma

export STAMP=`date +"%Y-%m-%d_%H-%M-%S"`
export DOCKER_IMAGE_TAG="image-$STAMP"
echo "$DOCKER_IMAGE_TAG"

docker pull $REPO

docker build -f docker-context/Dockerfile \
             --cache-from $REPO \
             --tag image-$STAMP \

# Run a container in CI so that we can grab the tag from the docker cli
if [[ -z "${CONTINUOUS_INTEGRATION}" ]]; then
    docker run -d -p $DOCKER_IMAGE_TAG 
fi