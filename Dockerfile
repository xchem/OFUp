FROM continuumio/miniconda3

#ARG name
#ARG institution
#ARG email
#ARG city
#ARG country
#
#ENV name=$name
#ENV institution=$institution
#ENV email=$email
#ENV city=$city
#ENV country=$country

SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y libxrender-dev && apt-get install -y autodock-vina

# add conda to bashrc
RUN echo 'export PATH="/opt/conda/bin:$PATH"' >> ~/.bashrc
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc

# source bashrc
RUN source ~/.bashrc

# create conda env from environment.yml
ADD environment.yml /tmp/environment.yml
WORKDIR /tmp
RUN conda env create


#ADD htmd_register.py /tmp
RUN /bin/bash -c "source ~/.bashrc; conda activate ofup; \
jupyter-nbextension enable nglview --py --sys-prefix; \
jupyter-nbextension enable --py --sys-prefix widgetsnbextension"
#python /tmp/htmd_register.py"

WORKDIR /
RUN rm -rf /tmp

RUN mkdir /code
ADD . /code

RUN chmod 775 /code


