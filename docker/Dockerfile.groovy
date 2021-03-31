# DockerFile for a Fludity development container

# Use a Groovy base image
FROM ubuntu:groovy

# This DockerFile is looked after by
MAINTAINER Tim Greaves

# Installs shouldn't expect input
ENV DEBIAN_FRONTEND=noninteractive

# Package updates and installs
RUN apt-get update && \
     apt-get -y dist-upgrade && \
     apt-get -y install gnupg dirmngr && \
     echo "deb http://ppa.launchpad.net/fluidity-core/ppa/ubuntu groovy main" > /etc/apt/sources.list.d/fluidity-core-ppa-groovy.list && \
     gpg --keyserver keyserver.ubuntu.com --recv 0D45605A33BAC3BE && \
     gpg --export --armor 33BAC3BE | apt-key add - && \
     apt-get update && \
     echo "Europe/London" > /etc/timezone && \
     apt-get -y install fluidity-dev texlive-pstricks texlive texlive-latex-extra texlive-science && \
     rm -rf /var/cache/apt && \
     rm -rf /var/lib/apt/lists

ENV PETSC_DIR /usr/lib/petscdir/3.13
ENV LDFLAGS -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi
ENV CPPFLAGS -I/usr/include/hdf5/openmpi
ENV OMPI_MCA_btl_vader_single_copy_mechanism none

# Add a Fluidity user who will be the default user for this container
# Make sure the user has a userid matching the host system
# -- pass this as an argument at build time
ARG userid=1000
RUN adduser --disabled-password --gecos "" -u $userid fluidity

USER fluidity
WORKDIR /home/fluidity
