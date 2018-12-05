#!/usr/bin/env bash
set -e

REPO=scottx611x/sigma

export STAMP=`date +"%Y-%m-%d_%H-%M-%S"`
export SUFFIX=-standalone
echo "image-$STAMP"

CONTAINER_NAME="container-$STAMP$SUFFIX"
mkdir "/tmp/$CONTAINER_NAME"

docker pull $REPO

docker build -f docker-context/Dockerfile \
             --cache-from $REPO:latest \
             --tag image-$STAMP \
             .
# docker run --name $CONTAINER_NAME \
#            --detach \
#            -p 3242:3242 \
#            image-$STAMP
rm -r "/tmp/$CONTAINER_NAME"