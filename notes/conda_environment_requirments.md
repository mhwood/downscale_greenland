## Python Packages
To use the python scripts provided in this directory, several packages are required. 

It is recommended that a fresh conda environment is established to include these packages. To create a new environment, use
```
conda create --name mitgcm
conda activate mitgcm
```
Then, activate the package and download the following packages:
```
conda install numpy
conda install matplotlib
conda install -c conda-forge xesmf
pip install ecco-v4-py
```
The required [simplegrid](https://github.com/nasa/simplegrid) package is not available by `pip` or `conda install` so it must be cloned and then installed: 
```
git clone https://github.com/nasa/simplegrid.git
cd simplegrid
pip install .
```
Then, add the MITgcm utils to your environment:
```
cd /path/to/MITgcm/utils/python/MITgcmutils
python setup.py install
```
