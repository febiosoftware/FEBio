FROM febiosoftware/febiobase:ubuntu-22.04
SHELL ["/bin/bash", "-c"]
RUN pwd
# Setup dependency installers
ARG IMAGE_BUILD_PATH=/tmp/linux
ARG DEPENDENCIES_PATH=${IMAGE_BUILD_PATH}/dependencies
RUN mkdir -p ${IMAGE_BUILD_PATH}

COPY ./common/linux ${IMAGE_BUILD_PATH}
ARG IMAGE_BUILD_PATH=/tmp/linux
ARG DEPENDENCIES_PATH=${IMAGE_BUILD_PATH}/dependencies
WORKDIR ${IMAGE_BUILD_PATH}
RUN ${IMAGE_BUILD_PATH}/dependencies/install.sh
USER febio
WORKDIR /FEBio
CMD ["sleep", "infinity"]
