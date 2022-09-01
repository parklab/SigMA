#!/usr/bin/env bash
set -e

REPO=parklab/sigma

export STAMP=`date +"%Y-%m-%d_%H-%M-%S"`
export DOCKER_IMAGE_TAG="image-$STAMP"
echo "$DOCKER_IMAGE_TAG"

docker pull $REPO

if [[ -v CONTINUOUS_INTEGRATION ]]; then 
    # Speed up CI build times with docker caching workaround. See: https://github.com/moby/moby/issues/20316#issuecomment-358260810
    docker build \
        -f docker-context/Dockerfile \
        --cache-from rocker/shiny:3.5.0,$REPO:latest \
        --tag image-$STAMP \
        .;

    # Run a container in CI so that we can grab its tag from the docker cli
    docker run -d -p 3242:3242 $DOCKER_IMAGE_TAG; 
else
    # If building locally, just use cache relative to the local machine
    docker build -f docker-context/Dockerfile .;
fi
