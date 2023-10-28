#! /bin/bash
docker build -t $REGISTRY/$REPOSITORY:$IMAGE_TAG -f infrastructure/DockerfileRuntime .
docker push $REGISTRY/$REPOSITORY:$IMAGE_TAG
