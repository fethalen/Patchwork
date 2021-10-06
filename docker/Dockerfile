# Build Docker Image with `sudo docker build --tag patchwork .` (from inside this directory).
# This will create an image named "patchwork" and tag it "latest". 
# You can also tag your image after running `sudo docker build .`: 
# `sudo docker tag IMAGEID REPONAME/IMAGENAME:TAG`
# Run `sudo docker images` to see all your images. 
# Run `sudo docker run patchwork` to view Patchwork's help menu (default behaviour).
# Else run like `sudo docker run patchwork ./patchwork ARGS`.
FROM ubuntu:20.04 AS build

RUN apt-get update --fix-missing && apt-get install -y sudo apt-utils --fix-missing
RUN sudo apt-get update --fix-missing && sudo apt-get install -y wget gcc --fix-missing

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
ENV PATH /home/patchwork/julia-1.6.3/bin:/home/patchwork/diamond:$PATH
RUN wget https://github.com/fethalen/Patchwork/archive/v$VERSION.tar.gz \
    && tar xvzf v$VERSION.tar.gz && cd Patchwork-$VERSION \ 
    && julia src/compile.jl . src/precompiled.jl ../patchwork-$VERSION

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
ENV PATH /home/patchwork/patchwork:/home/patchwork/diamond:$PATH

CMD ./patchwork --help