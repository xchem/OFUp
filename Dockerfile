FROM continuumio/miniconda3

ARG name
ARG institution
ARG email
ARG city
ARG country

ENV name=$name
ENV institution=$institution
ENV email=$email
ENV city=$city
ENV country=$country

SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y libxrender-dev && apt-get install -y autodock-vina

# add conda to bashrc
RUN echo 'export PATH="/opt/conda/bin:$PATH"' >> ~/.bashrc
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc

# source bashrc
RUN source ~/.bashrc


ADD htmd_register.py /tmp
RUN /bin/bash -c "source ~/.bashrc; conda create --name ofup; conda activate ofup; \
conda install -y -c acellera -c psi4 htmd; conda upgrade nglview --force; conda upgrade -y jupyter --force; \
python3 -m pip install papermill; pip install nteract-scrapbook; \
jupyter-nbextension enable nglview --py --sys-prefix; python /tmp/htmd_register.py"

WORKDIR /
#RUN rm -rf /tmp

RUN mkdir /code
ADD . /code

RUN chmod 775 /code


