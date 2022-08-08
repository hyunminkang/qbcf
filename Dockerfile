# Source Image
 FROM ubuntu:latest

# Set noninterative mode
 ENV DEBIAN_FRONTEND noninteractive

# apt update and install global requirements
 RUN apt-get clean all && \
     apt-get update && \
     apt-get upgrade -y && \
     apt-get install -y  \
         autoconf \
         build-essential \
         cmake \
         git \
         libbz2-dev \
         libcurl4-openssl-dev \
         libssl-dev \
         zlib1g-dev \
         liblzma-dev

# apt clean and remove cached source lists
 RUN apt-get clean && \
     rm -rf /var/lib/apt/lists/*

 RUN git clone https://github.com/samtools/htslib.git
 RUN cd htslib && \
     git submodule update --init --recursive && \
     autoreconf -i && \
     ./configure --prefix=/usr/local/ && \
     make && \
     make install
  
# Install qgenlib
 RUN git clone https://github.com/hyunminkang/qgenlib.git
 RUN cd qgenlib && \
     mkdir build && \
     cd build && \
     cmake .. && \
     make
 RUN cp /qgenlib/lib/libqgen.a /usr/local/lib/

# Install qbcf
 RUN git clone https://github.com/hyunminkang/qbcf.git
 RUN cd qbcf && \
     mkdir build && \
     cd build && \
     cmake .. && \
     make
 RUN cp /qbcf/bin/qbcf /usr/local/bin

# Define default command
# COPY ./entrypoint.sh /
# RUN chmod 755 /entrypoint.sh
# ENTRYPOINT ["/entrypoint.sh"]
 CMD ["qbcf"]
