FROM --platform=linux/amd64 ubuntu:20.04

RUN apt-get update && apt-get install -y vim wget curl git

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
RUN bash ~/miniconda.sh -b -p /miniconda
ENV PATH="/miniconda/bin:$PATH"

RUN conda install -y -c conda-forge jupyter jupyter_contrib_nbextensions python=3.11
RUN pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
RUN conda install -y -c conda-forge scanpy seaborn tifffile ome-types scikit-image scikit-learn pyyaml einops

# get env file seperately so it doesn't reinstall every time
# COPY ./env.yaml /pecgs-cnv/env.yaml
# RUN conda env create --file /pecgs-cnv/env.yaml 

# RUN conda create -n sc_analysis -y
COPY . /src

# ENV PATH="/miniconda/envs/sc_analysis/bin:$PATH"


CMD /bin/bash
