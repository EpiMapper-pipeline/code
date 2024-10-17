#!/bin/bash
echo "remove epimapper_lx"
conda remove -n epimapper_lx --all --yes

echo "create epimapper_lx"
conda create --name epimapper_lx python==3.9.19 --yes

echo "activate epimapper_lx"
conda activate epimapper_lx 

echo "install required packages in epimapper_lx"
conda install --yes --file requirements.txt

#alternatively, if conda install for requirements not work then try to install the requirements using conda environment :
#conda env create -f environment.yml_conda_lx

echo "install epimapper"
conda activate epimapper_lx
python setup.py install

#or  python -m pip uninstall epimapper
python -m pip install . 
#python -m build 

echo "run demo"
../demo/demo_histone.sh



