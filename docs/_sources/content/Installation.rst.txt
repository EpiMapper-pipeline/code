=================
Installation
=================


.. contents::
    :local:


Download
============

epimapper is written in Python. It can be installed and accessed from the command line and is available for both Linux and macOS operating systems. The package can be downloaded from the  `EpiMapper GitHub page <https://github.com/EpiMapper-pipeline/code/tree/main/epimapper>`_. Alternatively, you can use the following command:

.. code-block:: bash

   $ wget 

Install
============

It is highly recommended to create a separate virtual environment for the package to avoid any library conflicts. You can create a virtual environment using the following commands. We recommend installing and using Miniconda/Anaconda (`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_). You can find a tutorial on creating and updating virtual environments here at `the Conda website <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_.

If Miniconda is already installed, you can proceed with the following step-by-step installation. To set-up the EpiMapper Python package you can follow the step-by-step details provided below:

1. Create a virtual environment:

.. code-block:: bash

   $ conda create --name epimapper python==3.8.15

2. Activate the virtual environment:

.. code-block:: bash

   $ conda activate epimapper

3. Install pip if not already installed:

.. code-block:: bash

   $ conda install pip


Please allow any other installations when prompted.

Requirements
============

Before installing the package, the dependencies must be fulfilled. It is advised to install the dependencies using Miniconda. The list of dependencies is as follows:

- matplotlib=3.7.1
- numpy=1.23.0
- deeptools=3.5.1 
- picard=3.0.0
- openjdk=17.0.8 
- seaborn=0.12.2
- plotnine=0.10.1 
- pandas=1.2.5 
- scipy=1.10.1
- fastqc=0.12.1
- MACS2=2.2.9.1
- bedtools=2.31.0
- samtools=1.18 
- bowtie2=2.5.1



These dependencies can be installed one by one using the Conda package manager. For example:

.. code-block:: bash

   $ conda install numpy==1.23.0

A requirements.txt file is provided with the package. You can automatically install all the requirements using the following command:

.. code-block:: bash

   $ conda install --file requirements.txt

Alternatively, you can install the requirements using pip:

.. code-block:: bash

   $ pip install -r requirements.txt


Package Installation
======================

To install the package, navigate to the ``epimapper`` directory (the folder containing setup.py and pyproject.toml) and run the following command:

.. code-block:: bash

   $ python setup.py install

For more details, refer to the readme file in the package.

Package Contents
==================

The package folder will contain the following:

- ``demo``: Contains function scripts.
- ``epimapper``: Contains the Python source code of the pipeline.
- ``readme.txt``: Instructions about the usage of the package.
- ``requirements.txt``: List of requirements that can be used for automatic installation using Miniconda or pip.
- ``setup.py``: Setup file for the package.
- ``seacr_tool``: Contains the SEACR shell scripts. 