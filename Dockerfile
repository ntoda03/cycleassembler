FROM ubuntu:20.04
LABEL authors="Nicholas Toda" \
      description="Docker image containing all software requirements for the cyclerassembler pipeline"

RUN apt-get update 
RUN apt-get install -y wget 

RUN wget https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh \
    && bash Anaconda3-2019.03-Linux-x86_64.sh -b -p /usr/local/bin/anaconda \
    && rm Anaconda3-2019.03-Linux-x86_64.sh
ENV PATH=/usr/local/bin/anaconda/bin:$PATH

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/cyclerassembler-1.0/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name cyclerassembler-1.0 > cyclerassembler-1.0.yml

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
