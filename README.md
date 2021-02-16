
# BART3D v1.1

BART3D is a computational method for inferring transcriptional regulators (TRs) associated with genome-wide differential chromatin interactions (DCIs) detected by comparing two Hi-C maps.

# Introduction

BART3D leverages over 7,000 human TR binding profiles and over 5,000 mouse TR binding profiles from the public domain (collected in <a href="http://cistrome.org/db/">Cistrome Data Browser</a>).

BART3D is implemented in Python and distributed as an open-source package along with necessary data libraries.

**BART3D is developed and maintained by the <a href="https://faculty.virginia.edu/zanglab/">Chongzhi Zang Lab</a> at the University of Virginia.**


# Installation
### Prerequisites

BART3D uses Python's distutils tools for source installation. Before installing BART3D, please make sure Python3 and the following python packages are installed. We highly recommend the <a href="https://docs.anaconda.com/anaconda/install/">Anaconda environment</a>, which is easy to install all the required python packages.

- setuptools
- numpy
- pandas
- scipy
- cooler
- requests


### Download the source package and setup the configuration file

You have to download the Human or Mouse Data Library under your own directory before install BART3D. The unpacked libraries occupy 14GB hard drive storage in the download directory. 

```shell
wget https://virginia.box.com/shared/static/byufe8n6ft47hs4q7l6xxhyqnifjohm8.gz -O hg38_library.tar.gz
tar zxf hg38_library.tar.gz
wget https://virginia.box.com/shared/static/bxli2cc6zfj7h1llt9vx9km9u0d4vwhc.gz -O mm10_library.tar.gz
tar zxf mm10_library.tar.gz
```

To install a source distribution of BART3D, clone or download the package and go to the directory where you unpacked BART3D.

```shell
git clone https://github.com/zanglab/bart3d.git
cd bart3d
```

Modify the configure file (bart3d/bart3d.conf). For example, if you have the hg38_library (or mm10_library) downloaded in this directory: /abc/def/hg38_library (or /abc/def/mm10_library), then the bart.conf file should read:

```shell
[path]
hg38_library_dir = /abc/def/
mm10_library_dir = /abc/def/
```

### Global installation 
Install with root/administrator permission, or you have the <a href="https://docs.anaconda.com/anaconda/install/">Anaconda environment</a> prepared. By default, the script will install python library and executable codes globally.

```shell
python setup.py install
```

### Local installation 
If you want to install everything under a specific directory, for example, a directory as /path/to/bart3d/, use the following commands.

```shell
mkdir -p /path/to/bart3d/lib/pythonX.Y/site-packages 
export PYTHONPATH=/path/to/bart3d/lib/pythonX.Y/site-packages/:$PYTHONPATH 
python setup.py install --prefix /path/to/bart3d 
export PATH=/path/to/bart3d/bin/:$PATH
```

In this value, X.Y stands for the major–minor version of Python you are using (such as 3.5 ; you can find this with sys.version[:3] from a Python command line).

You’ll need to modify the environment variables and add those lines in your bash file (varies on each platform, usually is ~/.bashrc or ~/.bash_profile).

```shell
export PYTHONPATH= "/path/to/bart3d/lib/pythonX.Y/site-packages/:$PYTHONPATH"
export PATH="/path/to/bart3d/bin/:$PATH"
```


### Download the test data 
The test data can be downloaded here:
<!---
```shell
wget https://www.dropbox.com/s/3ghxmk9b6o3ak89/bart3d_test_data.tar.gz?dl=0 -O bart3d_test_data.tar.gz
tar zxf bart3d_test_data.tar.gz 
```
--->

- For hicpro input
```shell
wget https://www.dropbox.com/s/wfrjwm5pq2k6qo3/bart3d_test_data_hicpro_input.tar.gz?dl=0 -O bart3d_test_data_hicpro_input.tar.gz
tar zxf bart3d_test_data_hicpro_input.tar.gz
```

- For .hic input
```shell
wget https://www.dropbox.com/s/727gtfih44zcrxy/bart3d_test_data_juicer_input.tar.gz?dl=0 -O bart3d_test_data_juicer_input.tar.gz
tar zxf bart3d_test_data_juicer_input.tar.gz
```

- For .cool input
```shell
wget https://www.dropbox.com/s/gvtwug9qp52sps8/bart3d_test_data_cool_input.tar.gz?dl=0 -O bart3d_test_data_cool_input.tar.gz
tar zxf bart3d_test_data_cool_input.tar.gz
```


# Tutorial

### Essential options

**`-t TREATMENT, --treatment TREATMENT`**  

Treatment contact matrix file. The input file could be
.matrix file from HiC-Pro, .hic file from Jurcer, or
.cool file. Multiple samples can be provided as comma
separated list (no spaces around commas, e.g., -t
A,B,C).

**`-c CONTROL, --control CONTROL`**    

Control contact matrix file. Must be in the same
format/resolution as the treatment matrix. Multiple
samples can be provided as comma separated list (no
spaces around commas, e.g., -c A,B,C). 

**`-f {hicpro,hic,cool}, --fileFormat {hicpro,hic,cool}`**  

Format of the input matrix files. The following
options are available: "hicpro", "hic" or "cool". If
use "hicpro", additional abs/ord .bed file need to be
provided through --bedFileHicpro.

**`-s {hg38,mm10}, --species {hg38,mm10}`**  

Species, please choose from "hg38" or "mm10".

**`--genomicDistance GENOMICDISTANCE`**  

Genomic flanking regions for detecting differential
chromatin interactions. Default: 200000 (in bp).

**`--bedFileHicpro BEDFILEHICPRO`**  

For .hicpro format, please provide the abs/ord .bed
file from HiC-Pro. NOTES: for one species and one
resolution, the HiC-Pro output abs/ord .bed files are
identical.

**`--resolution RESOLUTION`**  

For .hic format, please specify a resolution. Default:
5000 (in bp).

**`--coverageNormalization`**                

Whether to normalize the DCI profile by the coverage
changes of Hi-C data. Default: FALSE.

**`--outdir OUTDIR`**         

If specified, all output files will be written to that
directory. Default: bart3d_output under the current
working directory

**`--outFileName OUTFILENAME`**  

Name string of the output files. Default: joint
basename of the first treatment/control files.


### Examples
              
BART3D accepts three formats of Hi-C contact maps as input:

1. Raw count maps from <a href="https://github.com/nservant/HiC-Pro">HiC-Pro</a>. In this case, additional coordinate index file (*_abs.bed or *_ord.bed file from HiC-Pro) need to be provided through **--bedFileHicpro**.
                        
```shell
bart3d -t t.matrix -c c.matrix -f hicpro --bedFileHicpro ord.bed -s hg38 --outdir hicpro_out --outFileName prename
```

2. .hic format files from <a href="https://github.com/aidenlab/juicer">Juicer</a>.

```shell
bart3d -t t1.hic,t2.hic -c c1.hic,c2.hic -f hic -s hg38 --outdir hic_out --resolution 5000 --outFileName prename
```

3. Hi-C contact maps in .cool format.

```shell
bart3d -t t1.cool,t2.cool,t3.cool -c c1.cool,c2.cool -f cool -s hg38 --outdir cool_out --outFileName prename
```
                      

### Output files

1. **\*_differential_score.bed**
is the genome-wide DCI profile which provides the DCI score at each genomic bin.

2. **\*_differential_score.bed.sorted.bw** 
is the bigwig format DCI profile for visualization.

3. **\*_Interaction_Decreased_auc.txt** 
provides the association score of each TR ChIP-seq dataset with the cis-regulatory profile that are related to DECREASED chromatin interactions in treat Hi-C map compared to control.

4. **\*_Interaction_Decreased_bart_results.txt** 
is a rank of all TRs with multiple quantification scores. Top ranked TRs are associated with DECREASED chromatin interactions in treat Hi-C map compared to control.

5. **\*_Interaction_Increased_auc.txt** 
provides the association score of each TR ChIP-seq dataset with the cis-regulatory profile that are related to INCREASED chromatin interactions in treat Hi-C map compared to control.

6. **\*_Interaction_Increased_bart_results.txt** 
is a rank of all TRs with multiple quantification scores. Top ranked TRs are associated with INCREASED chromatin interactions in treat Hi-C map compared to control.



<!-- 
# Citation

If you use BART in your data analysis, please cite: 

<a href="https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty194/4956015" target="_blank">BART: a transcription factor prediction tool with query gene sets or epigenomic profiles</a> <br>
Zhenjia Wang, Mete Civelek, Clint Miller, Nathan Sheffield, Michael J. Guertin, Chongzhi Zang. <i><b>Bioinformatics</b></i> 34, 2867–2869 (2018)
 -->

<!---
If you use "geneset" mode, please also cite:
<a href="http://genome.cshlp.org/content/26/10/1417" target="_blank">Modeling cis-regulation with a compendium of genome-wide histone H3K27ac profiles</a> <br>
Su Wang, Chongzhi Zang, Tengfei Xiao, Jingyu Fan, Shenglin Mei, Qian Qin, Qiu Wu, Xujuan Li, Kexin Xu, Housheng Hansen He, Myles Brown, Clifford A. Meyer, X. Shirley Liu. <i><b>Genome Research</b></i> 26, 1417–1429 (2016)
-->



