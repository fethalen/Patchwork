# FROM ubuntu:latest
FROM ubuntu:20.04 AS build

RUN apt-get update --fix-missing && apt-get install -y sudo apt-utils --fix-missing
RUN sudo apt-get update --fix-missing && sudo apt-get install -y wget gcc --fix-missing

RUN groupadd --gid 1000 patchwork \
    && useradd --uid 1000 --gid patchwork --shell /bin/bash --create-home patchwork \
    && chmod -R a+rwx /home/patchwork
WORKDIR /home/patchwork
USER patchwork
ENV VERSION 0.1.2-pre-alpha8

# Can we do this without conda? 
RUN wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
RUN bash Anaconda3-2020.11-Linux-x86_64.sh -b
RUN rm Anaconda3-2020.11-Linux-x86_64.sh 
ENV PATH /home/patchwork/anaconda3/bin:$PATH
RUN conda update conda
#RUN conda update anaconda
RUN conda update --all
RUN conda install -c conda-forge julia
RUN conda install -c bioconda diamond

# Try without conda: 
#RUN wget https://julialang.org/downloads/julia-1.6.3-linux-x86_64.tar.gz \
#    && tar xvzf julia-1.6.3-linux-x86_64.tar.gz \
#RUN wget http://github.com/bbuchfink/diamond/releases/download/v2.0.11/diamond-linux64.tar.gz \
#    && tar xvzf diamond-linux64.tar.gz
#ENV PATH /home/patchwork/julia-1.6.3/bin:/home/patchwork/diamond:$PATH
RUN wget https://github.com/fethalen/Patchwork/archive/v$VERSION.tar.gz \
    && tar xvzf v$VERSION.tar.gz && cd Patchwork-$VERSION \ 
    && julia src/compile.jl . src/precompiled.jl ../patchwork-$VERSION
RUN ls 


# To get the final image size down, use a multi-stage build:
# This contains neither a Conda/Anaconda nor a Julia installation, 
# just Ubuntu, the Patchwork binary and DIAMOND. 
FROM ubuntu:20.04
RUN groupadd --gid 1000 patchwork \
    && useradd --uid 1000 --gid patchwork --shell /bin/bash --create-home patchwork \
    && chmod -R a+w /home/patchwork
WORKDIR /home/patchwork
USER patchwork
ENV VERSION 0.1.2-pre-alpha8

COPY --from=build /home/patchwork/patchwork-$VERSION patchwork-VERSION
COPY --from=build /home/patchwork/anaconda3/bin/diamond diamond
#COPY --from=build /home/patchwork/diamond diamond
RUN ln -s patchwork-VERSION/bin/patchwork patchwork
ENV PATH patchwork:$PATH

CMD patchwork --help