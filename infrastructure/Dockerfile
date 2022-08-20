FROM ubuntu:18.04
SHELL ["/bin/bash", "-c"]
RUN pwd
# Setup dependency installers
ARG IMAGE_BUILD_PATH=/tmp/linux
ARG DEPENDENCIES_PATH=${IMAGE_BUILD_PATH}/dependencies
RUN mkdir -p ${IMAGE_BUILD_PATH}

COPY ./common/linux ${IMAGE_BUILD_PATH}
WORKDIR ${IMAGE_BUILD_PATH}
RUN ./apt.sh
RUN ./cmake.sh
RUN ./git.sh
RUN ${IMAGE_BUILD_PATH}/dependencies/install.sh

WORKDIR /FEBio
CMD ["sleep", "infinity"]
