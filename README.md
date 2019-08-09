# Open Follow-Up (OFUp)

## Introduction

This repository contains open-source tools for follow-up chemistry from fragment screening campaigns at XChem (https://www.diamond.ac.uk/Instruments/Mx/Fragment-Screening.html), but the tools can be applied anywhere else that is relevant.

There are example notebooks (see: https://github.com/xchem/OFUp/tree/master/example_notebooks) that are designed to guide the user through the process of selecting compounds for follow-up

These notebooks can be run from:
1. The containerised environment provided here (via. Docker container)
2. Manual installation (recommended: via. a conda environment)  


## Installation

This has only been  tested on Mac, but the beauty of docker containers is that they should run identically on any OS  


### Container install

Requirements: 
- Docker (https://www.docker.com/products/docker-desktop)
- git
- basic knowledge of the terminal  


**1. Clone this repository**

```
git clone https://github.com/xchem/OFUp.git
```  


**2. Build the docker container**

The HTMD package (https://www.acellera.com/products/high-throughput-molecular-dynamics/, https://github.com/Acellera/htmd) is used in this repository, and requires licensing (please refer to the HTMD documentation for licensing terms). Don't worry, we've looked after that for you, but you will have to provide some information when you are building your container in order to register.  


To build the container:
```
docker build -t ofup . --build-arg name='Some Person' --build-arg institution='Institution of Fun' --build-arg email='some.person@iofun.com' --build-arg city='Funland' --build-arg country='Funplace'
```

Replacing the name, institution, email, city and country values with your own details.   


That's it! We'll show you how to run the container in the next section  


### Anaconda install

Coming soon...  


## Usage

### Running example notebooks from the docker container

**1. Run the container as an interactive session**

```
docker run -it --rm -p 8888:8888 ofup /bin/bash
```

This does:
- ```docker run```: runs the container
- ```-it```: interactivley
- ```--rm```: removes the container automatically when you exit
- ```-p 8888:8888```: forwards port 8888 in the container to 8888 on your computer (default port for running jupyter notebooks)
- ```ofup```: the name of the container
- ```/bin/bash```: uses bash as the shell to interact with the container

