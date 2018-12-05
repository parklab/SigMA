#!/usr/bin/env bash
set -e

REPO=scottx611x/sigma

export STAMP=`date +"%Y-%m-%d_%H-%M-%S"`
echo "image-$STAMP"

docker pull $REPO

docker build -f docker-context/Dockerfile \
             --cache-from $REPO:latest \
             --tag image-$STAMP \
             .