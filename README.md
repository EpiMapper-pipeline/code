
# EpiMapper: A New Tool for Analyzing High Throughput Sequencing from CUT&Tag
## epimapper Documentation

EpiMapper is a complete pipeline for data analysis of epigenomic sequencing.

[Home](README.md) | [Quality Control of Raw Reads](docs/fastqc.md) | [Bowtie2 Alignment](docs/bowtie2_alignment.md) | [Removal of Duplicated Reads](docs/remove_duplicates.md) | [Evaluation of Fragment Length Distribution](docs/fragment_length.md) | [Quality Filtering and File Conversion](docs/filtering.md) | [Spike-in Normalization](docs/spike_in_calibration.md) | [Peak Calling](docs/peak_calling.md) | [Enrichment Visualization in Heatmaps](docs/heatmaps.md) | [Differential Analysis and Annotaion](docs/differential_analysis.md)

## Download:

epimapper is written in Python. It can be installed and accessed from command line and is avalible for both linux and MAC operating systems. 
	

## Installation:
<p>It is highly recommended to create a separate virtual environment for the package to avoid any library conflicts problem. You you create virtual environment using the following commands. We recommend to use install and use miniconda/anaconda (https://docs.conda.io/en/latest/miniconda.html). Tutorial of creating and updating virtual commands can be found at (https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) </p> 




<pre>conda create --name epimapper python==3.8
conda activate epimapper </pre>

<p>Install pip if not already installed: </p>
<pre> conda install pip</pre>

Please allow any other installations when prompted

<p>Prior to installing the package, dependencies must be fulfilled. It is advised to install dependencies using miniconda. List of dependencies is as follows: </p>
<ul>
  <li>matplotlib==3.6.2</li>
  <li>numpy==1.23.0</li>
  <li>pandas==1.5.2</li>
  <li>seaborn==0.12.2
  <li>plotnine==0.10.1</li>
  <li>setuptools==65.6.3</li>
  <li>bedtools==2.27.0</li>
  <li>fastqc==0.11.9</li>
  <li>bowtie2==2.2.5</li>
  <li>samtools==1.6</li>
  <li>deeptools==3.5.1</li>


</ul>

These dependencies can be installed one by one using conda manager. For example:

<pre>conda install numpy==1.23.0</pre>
	
A requirements.txt file is given with the package. All requiremnts can be automatically installed using one command:
<pre>conda install --file requirements.txt</pre>

Or can be installed using pip.

<pre>pip install numpy==1.23.0</pre>

A requirements.txt file is given with the package. All requiremnts can be automatically installed using one command:
<pre>pip install -r requirements.txt</pre>

You can install the package using following command, go to the dmr_analysis directory (folder containing setup.py and pyproject.toml) and type the following command
<pre>pip install .</pre>



		
## Contents of the package:
		
<p>The package folder will contain following:
	</p>
<ul>
	
	<li><code>epimapper</code> : Contains python soruce code of pipeline.</li>
	<li><code>readme.txt</code> : Instructions about usage of package.</li>
	<li><code>requirments.txt</code> :  List of requirements. Can be used for automatic installation from miniconda or pip.</li>
	<li><code>setup.py</code>: Setup file for package.</li>
	


</ul>	
	

	
## Pipeline Tasks:
	
<p>The pipeline consists of follwoing tasks. To run a task, type epimapper task args. To see what are the options for each task of the pipeline, please run: epimapper -h </p>

<ul>
 <li><code>fastqc </code> : Performs quality control on raw reads fastq files from high-thoughput sequencing.</li>
	<li><code>bowtie2_alignment </code> :Mapping reads to a reference genome with Bowtie2 alignment of fastq sequencing reads files from high-thoughput sequecing, and visulizing results. </li>
	<li><code>remove_duplicates</code> : Remove duplicated reads mapped to the same place in of a reference genome during alignment, and visulizing results.</li>
	<li><code>fragment_length</code> :Evaluation of mapped fragment length distribution of input SAM files exported from high-thoughput sequencing alignment, and visulizing results.</li>
	<li><code>filtering</code> : Performs data filtering for mapped reads based on their alignment quality, and file format conversion before high-level data analysis,before  visulizing reproducibility among biological replicates.</li>
	<li><code>spike_in_calibration</code> : Removes experimental bias by normalizing fragment counts based on sequencing depth to a spike-in genome and visulizes results.</li>
	<li><code>peak_calling</code> :Finds enriched regions/calls for peaks from chromatin profiling data with SEACR or MACS2, then visulizes results.</li>
	<li><code>heatmap</code> : Visualizes the enrichment of target protein in predefined genomic regions and peaks by creating heatmaps.</li>
	<li><code>differential_analysis</code>:Preforms differntial analysis on enriched reagions/peaks before annotating the stastitically significant changes to spesific genomic reagions and visulizing the results. </li>
	
</ul>
