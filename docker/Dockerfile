# Build Docker Image with `sudo docker build --tag patchwork .` (from inside this directory).
# This will create an image named "patchwork" and tag it "latest". 
# You can also tag your image after running `sudo docker build .`: 
# `sudo docker tag IMAGEID REPONAME/IMAGENAME:TAG`
# Run `sudo docker images` to see all your images. 
# Run `sudo docker run -it --name CONTAINER --user $(id -u) \
# --mount type=bind,src=HOSTPATH,dst=CONTPATH IMAGENAME:TAG bash`
# (e.g. `sudo docker run -it --name patchwork_container --user $(id -u) \
# --mount type=bind,src=/home/clara/Patchwork/test,dst=/home/patchwork/data patchwork:latest bash`)
# to create and start a container in interactive mode with an open bash terminal. 
# HOSTPATH is the absolute path to the host's directory containing the data needed to execute patchwork 
# (e.g. contig and reference fasta files). CONTPATH is the absolute path to an empty directory in the 
# container where the contents of HOSTPATH is to be bind-mounted. The host's data is the accessible 
# through CONTPATH. 
# The --user option allows you to e.g. create an output directory inside the bind-mounted host's directory.
# Without this option, you might get an EACCESS: permission denied error when trying to modify the host's 
# directory because the host's user ID and your user ID inside the container could differ. 
# Not sure how to fix this in a more elegant way. 
# To run Patchwork inside the container, type something like 
# `patchwork --contigs CONTPATH/path/to/contigs --reference CONTPATH/path/to/reference --output-dir CONTPATH/path/to/outdir`
# Type `exit` to exit the docker container. 
# Run `sudo docker start CONTAINER` to restart the container and `sudo docker exec CONTAINER CMD` to
# enter the container/ execute CMD. 

FROM ubuntu:20.04 AS build

RUN apt-get update --fix-missing && apt-get install -y sudo apt-utils --fix-missing
RUN sudo apt-get update --fix-missing && sudo apt-get install -y wget gcc --fix-missing
# For debugging using the git repo instead of a release, use the following line: 
#RUN sudo apt-get update --fix-missing && sudo apt-get install -y wget git gcc --fix-missing

RUN groupadd --gid 1000 patchwork \
    && useradd --uid 1000 --gid patchwork --shell /bin/bash --create-home patchwork \
    && chmod -R a+rwx /home/patchwork
WORKDIR /home/patchwork
USER patchwork
ENV VERSION 0.1.2-pre-alpha9

RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.3-linux-x86_64.tar.gz \
    && tar xzvf julia-1.6.3-linux-x86_64.tar.gz 
RUN wget http://github.com/bbuchfink/diamond/releases/download/v2.0.11/diamond-linux64.tar.gz \
    && tar xvzf diamond-linux64.tar.gz
ENV PATH /home/patchwork/julia-1.6.3/bin:$PATH
RUN wget https://github.com/fethalen/Patchwork/archive/v$VERSION.tar.gz \
    && tar xvzf v$VERSION.tar.gz && cd Patchwork-$VERSION \ 
    && julia src/compile.jl . src/precompiled.jl ../patchwork-$VERSION
# For debugging using the git repo instead of a finished release:
#RUN git clone -b sequence-output --single-branch https://github.com/fethalen/Patchwork \
#    && cd Patchwork && julia src/compile.jl . src/precompiled.jl ../patchwork-$VERSION


# To reduce the final image size, use a multi-stage build:
# The final image doesn't contain a Julia installation, 
# just Ubuntu as base image, the Patchwork binary and DIAMOND. 
FROM ubuntu:20.04
RUN groupadd --gid 1000 patchwork \
    && useradd --uid 1000 --gid patchwork --shell /bin/bash --create-home patchwork \
    && chmod -R a+rwx /home/patchwork
WORKDIR /home/patchwork
USER patchwork
ENV VERSION 0.1.2-pre-alpha9

COPY --from=build /home/patchwork/patchwork-$VERSION patchwork-$VERSION
COPY --from=build /home/patchwork/diamond diamond
RUN ln -s patchwork-$VERSION/bin/patchwork patchwork
ENV PATH /home/patchwork:$PATH

CMD patchwork --help