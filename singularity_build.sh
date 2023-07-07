Bootstrap: docker
From: ubuntu:latest

%files
../calibrate

%post
echo "Updating apt repositories"
apt update && apt upgrade -y
DEBIAN_FRONTEND=noninteractive apt install -y \
            wget \
            rsync \
            zip \
            git \
	    build-essential \
	    cmake \
	    libcfitsio-dev \
	    libgsl-dev \
	    casacore-dev \
	    libboost-dev  \
            libboost-filesystem-dev \
            libboost-date-time-dev \
            libboost-system-dev \
            libboost-thread-dev 

cd / \
    && ln -s /usr/include 

cd /calibrate && \
	. ./build_container.sh

