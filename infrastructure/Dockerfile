FROM febiosoftware/febiobase:ubuntu-22.04
SHELL ["/bin/bash", "-c"]
RUN pwd
# Setup dependency installers
ARG IMAGE_BUILD_PATH=/tmp/linux
ARG DEPENDENCIES_PATH=${IMAGE_BUILD_PATH}/dependencies
RUN mkdir -p ${IMAGE_BUILD_PATH}

COPY ./common/linux ${IMAGE_BUILD_PATH}
ARG IMAGE_BUILD_PATH=/tmp/linux
ARG BUILD_PATH=/tmp/src
RUN mkdir -p ${BUILD_PATH}
ARG DEPENDENCIES_PATH=${IMAGE_BUILD_PATH}/dependencies
WORKDIR ${DEPENDENCIES_PATH}
RUN ${IMAGE_BUILD_PATH}/dependencies/install.sh
RUN ${DEPENDENCIES_PATH}/hypre.sh
RUN ${DEPENDENCIES_PATH}/levmar.sh
RUN ${DEPENDENCIES_PATH}/mmg.sh
RUN ${DEPENDENCIES_PATH}/tetgen.sh
RUN ${DEPENDENCIES_PATH}/itk.sh
RUN ${DEPENDENCIES_PATH}/sitk.sh
RUN ${DEPENDENCIES_PATH}/occt.sh
RUN ${DEPENDENCIES_PATH}/netgen.sh
WORKDIR /FEBio
CMD ["sleep", "infinity"]
