============================
Detailed Installation Guide
============================



.. contents::
    :local:


Bioinformatics can be a daunting task at first, therefore we decided to create a super simple introduction to how you can get started analyzing your own data!


Downloading the EpiMapper Python Package
========================================
To begin, please visit the `EpiMapper GitHub page <https://github.com/EpiMapper-pipeline/code>`_, click the green "<> Code" button, and click "download zip" to download the package to your computer.
Then you may unzip the EpiMapper folder and move it wherever you wish to have keep your Python Packages. I would recommend to create a new folder on your desktop or in you documents folder and keep all downloaded tools there to easily locate them later. 
Now, before installing the Epimapper Python package, you most download the software that are dependencies of EpiMapper. 


Package Manager
=================


A quick little disclaimer: 
The EpiMapper Python package works for MAC OS and Linux environments. If you are using a Windows computer, more steps may be necessary, like setting up Ubuntu (https://ubuntu.com/download/desktop), an open-source operation system on Linux. 

Firstly, you have to open your terminal window, either locally on your computer or sign in to a SSH. Remember, when copying and pasting the following code bellow, you do not need include the "$" sign. This just symbolizes a new line in the terminal.
Here is a quick simple guide to get you started in the terminal: `Command Line Basics and Usefull Tricks With the Terminal <https://code.tutsplus.com/command-line-basics-and-useful-tricks-with-the-terminal--cms-29356t>`_

Next, we strongly advise you download and use a package manager to create virtual environments, if you do not already have it. Think of a package manager and a virtual environment as a separation box between software that might not like each other. Imagine you are cooking your favorite dish, let’s say pancakes: you’ll need eggs, milk, flour etc, but we can agree that you don’t need broccoli. Would you take the broccoli out of the fridge then? No, you would not, you don’t need it, so why risk taking out the broccoli, when there is a chance some of it will get into your mixing bowl and ruin the pancakes? This is the same with software, some software don´t mix very well, so it´s smart to just collect the ones we need for a specific task. 

Another example is if your vegan friend is coming over for pancakes, then you will need vegan milk, not cow’s milk. The result is pretty much the same: pancakes, but two different versions: vegan or ordinary. For software, one version may have some desirable features while the next version does not, therefore it is important to collect and use the version we need to produce the desired product. 

We recommend installing and using Miniconda/Anaconda (`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_). You can find a tutorial on creating and updating virtual environments `here <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_.

Further, to make installation easier, I would recommend to install `Libmamba <https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community>`_, this makes it faster to find and solve issues when downloading a lot of packages. It may be easily applyed in condas ``--solver`` parameter:

.. code-block:: bash

   $ conda install numpy==1.23.0 --solver=libmamba



Installing packages - Conda
===========================

Like the example with numpy above, we can donwload and install packages with conda by simply saying which package conda should install.
However, for some packages, there are a special parameter that tells conda where to look for the specific package you want to install.
In the example bellow we are installing the commenly used tool `BEDtools <https://bedtools.readthedocs.io/en/latest/>`_ , we have in the ``-c`` paramter specified that conda should look in bioconda to find the BEDtools software. 

.. code-block:: bash

    $ conda install -c bioconda bedtools


As a one time initial set up of bioconda, we need to write this into the terminal:

.. code-block:: bash

    $ conda config --add channels defaults
    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge
    $ conda config --set channel_priority strict 

Remember, this is just a one time set up, so it does not need to be done before donwloading every bioconda package. 



Installing packages - Pip
=========================

`Pip <https://pypi.org/project/pip/>`_ is a package installer for Python. It is easy to use, and may be used in collaberation with conda. Once you have created and activated a virutal enviourment by using conda, you may download packages just like with conda:

.. code-block:: bash

    $ pip install plotnine==0.12.1 

Why use pip when you already have conda? Well, sometimes conda does not want to co-operate, so if that happens, pip might save your day.  Other times we just need to give the computer some time, and it will fix the problem, remember to have some patience, it might take some time to download and install the packages.


Requirements for Epimapper
==========================

 The Epimapper Python package uses a lot of packages, and it is very importat that the version of each package is correct. Therefore, we recommend you to use the ``requirements.txt`` file, avalible in the package, to download the necessary packages.
 
 In short:

1. Make your desired virtual environment:

.. code-block:: bash

    $ conda create --name epimapper python==3.8.15


2. ``cd`` into the EpiMapper directory you downloaded from Github. 

Now you may install all the required software with conda:


3. 

.. code-block:: bash

    $ conda install --file requirements.txt

Or install the software with pip:


4.

.. code-block:: bash

    $ pip install -r requirements.txt


Now this step may take some time, be patient. If the ``requirements.txt`` file failes, please try to install the packages one by one as described above.


When you think you have downloaded all the software, you can check by running these command and see if all the software in the ``requirements.txt`` comes up:

.. code-block:: bash

    $ conda list


.. code-block:: bash

   $ pip freeze 

If everything seems okay you may move on to installing the EpiMapper Python package, but remember it may take some time and if you run into any problems try to search on the internet to find solutions. 



Installing the EpiMapper Python package
=======================================

The last step is to install the EpiMapper, which is easily done. ``cd`` into the EpiMapper folder you downloaded from GitHub and run this simple command 

.. code-block:: bash

    $ python setup.py install 


After this you should be ready to use EpiMapper Python package. You may try to see if everything is ready by running the following command:

.. code-block:: bash

    $ epimapper -h


The output should be something like this: 

.. code-block:: bash

    usage: epimapper <task> [<args>]

    Tasks available for using:

    fastqc                      Quality control for raw reads in fastq files

    bowtie2_alignment           Alignment of fastq files to a reference genome

    remove_duplicates           Removes duplicated reads

    fragment_length             Finds and plots fragment lengths in sam files

    filtering                   Filters and converts files from bam to bed

    spike_in_calibration        Spike in normalizes input files

    peak_calling                Calls for enriched regiones in bed files

    heatmap                     Plots heatmaps of enriched regions

    differential_analysis       Preforms differential analysis on input bed files

    Python Package for analysing epigenetic data such as ChIP-seq, CUT&RUN, ATAC-seq and CUT&Tag

    positional arguments:
    task        Pipeline task to run

    optional arguments:
    -h, --help  show this help message and exit   


