FROM nvidia/cuda:11.1.1-base-ubuntu20.04
LABEL authors="Adyota Gupta"
ENV DEBIAN_FRONTEND=noninteractive

RUN rm /etc/apt/sources.list.d/cuda.list

## TODO: Clean up dependencies!
RUN apt-get update && apt-get install -y \
    curl \
    ca-certificates \
    sudo \
    bzip2 \
    libx11-6 \
    bc

## install anaconda
ENV CONDA_PATH=/AnacondaInstall
WORKDIR /tmp
RUN curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
RUN chmod +x /tmp/Miniforge3-Linux-x86_64.sh
RUN mkdir /root/.conda
RUN bash -c "/tmp/Miniforge3-Linux-x86_64.sh -b -p $CONDA_PATH"
RUN $CONDA_PATH/bin/mamba init bash
RUN $CONDA_PATH/bin/mamba update -n base -c defaults mamba -y
ENV PATH="/AnacondaInstall/bin/:${PATH}"

## install matplotlib, sympy, and cloudpickle
RUN bash -c "mamba install matplotlib sympy cloudpickle scipy -y"

## set normal and set code to env
WORKDIR /
ENV CODE_PATH="/code"

## no errors print out  when running mpirun
ENV OMPI_MCA_btl_vader_single_copy_mechanism=none
