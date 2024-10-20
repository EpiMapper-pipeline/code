#linux shall use mamba
#mac with conda install
#mamba=2.0.2
#conda=23.3.1
########################
## Add channel in configure file .condarc such as
########################
# channels:
#  - defaults
#  - bioconda
#  - anaconda
#  - conda-forge


######################## 
#0. remove old env before installing
########################
#mamba remove -n epimapper_mc --all --yes
#or
#conda remove -n epimapper_mc --all --yes

########################
#clean up old files
########################
#mamba clean --all --yes
#and
#conda clean --all --yes

########################
#1.1 Type 1. make env works in both mac and linux
########################
mamba env create -f environment.yml --name epimapper_mc
#or 
#conda env create -f environment.yml --name epimapper_mc 

#######################
#1.2 Type 2. make env works in mac but mamba may not work in linux
#######################
#first
#conda create --name epimapper_mc --yes

#then
##conda activate epimapper_mc
#or
#source activate epimapper_mc

#finally
#mamba install --insecure --yes --file requirements.txt
#or
##conda install --yes --file requirements.txt

#######################
#1.3 Type 3. if both type 1 and typ2 do not work then please install all packages manually
#######################
#For example
#mamba remove -n epimapper_mc --all
#mamba create -n epimapper_mc python=3.9.19 --yes
#conda activate epimapper_mc
#mamba install samtools=1.6 deeptools=3.5.1 bedtools=2.31.1 --yes
#mamba install seaborn=0.13.0 setuptools=69.5.1 argparse pandas=1.4.4 numpy=1.23.0 matplotlib=3.8.4 --yes
#mamba install plotnine=0.10.1 scipy=1.9.1 pathlib=1.0.1 scikit-learn=1.5.0 openjdk=22.0.1 --yes
#mamba install fastqc=0.12.1 MACS2=2.2.7.1 cutadapt picard=3.1.1 bowtie2=2.2.5 --yes

#######################
#2. install epimapper
#######################
echo "install epimapper"
conda activate epimapper_mc

# uninstall 
#python -m pip uninstall epimapper

# install
python -m pip install .
#or
#python setup.py install

#######################
#3. run demo
#######################
echo "run demo"
#go to 
cd ../demos/demo_histone/ 

#then run 
./demo_histone.sh 


