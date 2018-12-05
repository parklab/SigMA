#!/usr/bin/env bash
set -e

REPO=scottx611x/sigma

export STAMP=`date +"%Y-%m-%d_%H-%M-%S"`
export DOCKER_IMAGE="image-$STAMP"
echo "$DOCKER_IMAGE"

docker pull $REPO

docker build -f docker-context/Dockerfile \
             --cache-from $REPO \
             --tag image-$STAMP \
             .