#!/usr/bin/env bash
set -e

REPO=parklab/sigma

export STAMP=`date +"%Y-%m-%d_%H-%M-%S"`
export DOCKER_IMAGE_TAG="image-$STAMP"
echo "$DOCKER_IMAGE_TAG"

docker pull $REPO

docker build -f docker-context/Dockerfile_travis .;
