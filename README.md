# :exclamation: Under construction and untested :exclamation:

# Open Follow-Up (OFUp)

## Introduction

This repository contains open-source tools for follow-up chemistry from fragment screening campaigns at XChem (https://www.diamond.ac.uk/Instruments/Mx/Fragment-Screening.html), but the tools can be applied anywhere else that is relevant.

There are example notebooks (see: https://github.com/xchem/OFUp/tree/master/example_notebooks) that are designed to guide the user through the process of selecting compounds for follow-up

These notebooks can be run from:
1. The containerised environment provided here (via. Docker container)
2. Manual installation (recommended: via. a conda environment)  
<br/><br/>

## Installation

This has only been  tested on Mac, but the beauty of docker containers is that they should run identically on any OS  


### Container install

Requirements: 
- Docker (https://www.docker.com/products/docker-desktop)
- git
- basic knowledge of the terminal  
<br/>

**1. Clone this repository**

```
git clone https://github.com/xchem/OFUp.git
```  
<br/>

**2. Build the docker container**

The HTMD package (https://www.acellera.com/products/high-throughput-molecular-dynamics/, https://github.com/Acellera/htmd) is used in this repository, and requires licensing (please refer to the HTMD documentation for licensing terms). Don't worry, we've looked after that for you, but you will have to provide some information when you are building your container in order to register.  


To build the container:
```
docker build -t ofup . \
--build-arg name='Some Person' \
--build-arg institution='Institution of Fun' \
--build-arg email='some.person@iofun.com' \
--build-arg city='Funland' \
--build-arg country='Funplace'
```

Replacing the name, institution, email, city and country values with your own details.   


That's it! We'll show you how to run the container in the next section 

<br/>

---
### Anaconda install (Coming soon...)
---

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

---
  
**You may also want to mount a folder from your machine into the container, so that work you do in it is not lost:**

```
docker run -it --rm -p 8888:8888 \
--mount type=bind,source=/Users/res3/michellab/XChem-examples/KALRNA/,target=/Rachael/KALRNA \
ofup /bin/bash
```

The mount command here binds ```/Users/res3/michellab/XChem-examples/KALRNA/``` on the local machine to ```/Rachael/KALRNA``` in the container. 

That means that any files I change in the container under ```/Rachael/KALRNA``` will also change on my computer in ```/Users/res3/michellab/XChem-examples/KALRNA/```.

---

<br/>

**2. Run the notebook server from the container**

Theres a script to do this! From insude of the container:

```
/code/start_jupyter.sh
```

If there's a permissions error:

```
chmod 775 /code/start_jupyter.sh
/code/start_jupyter.sh
```

## Advanced usage

The example notebooks can also be executed and analysed in new jupyter notebooks, or via the command line. This uses an awesome package called papermill: https://github.com/nteract/papermill/

Interacting with the output from notebooks run and produced with papermill can be done with another awesome package, scrapbook: https://github.com/nteract/scrapbook


### Example:

to run a notebook from an interactive python session or from a new jupyter notebook:

```python
import papermill as pm

pm.execute_notebook('dock_multi_against_ref.ipynb', 
                    'test_interactive_docking.ipynb',
                    parameters=dict(pdb_file = '/data/XX02KALRNA-x1389_1/XX02KALRNA-x1389_1.pdb',
                                   ligands_file = '/data/XX02KALRNA-x1389_1/250_surprise_sucos.sdf',
                                   ref_file = '/data/XX02KALRNA-x1389_1/XX02KALRNA-x1389_1.sdf',
                                   conf_file = 'conf.txt',
                                   docking_directory = '/data/docking/XX02KALRNA-x1389_1/surprise_set'))
```



