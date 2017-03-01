# ctamacros

Project description
==============

ROOT-based macros used by CTA PHYS group to compute sensitivity 
 

Dependencies
--------------

We should try to stick with just 3 dependencied: ROOT, MARS and CFITSIO.

Probably these packages would be useful:

```shell
git dpkg-dev make g++ gcc binutils libx11-dev libxpm-dev libxft-dev libxext-dev cmake libcfitsio3-dev
```


Installation
--------------

Currently, there is a test macro working: mfits/testMEventList.C

It only uses CFITSIO and MARS/ROOT as dependencies. 

To install CFITSIO, download the package, untar and install following these instructions:

```shell
wget http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio_latest.tar.gz
tar -xzvf cfitsio_latest.tar.gz
cd cfitsio
mkdir build
./configure --prefix=$PWD/build
make
make install
echo "export CFITSIO=$PWD/build" >> ~/.bashrc
source ~/.bashrc
```

Once this is done, then you may build this repository using (within the magic_dl3/ path):

```shell
./setup.sh
cd build && make all && cd ..
```

This creates an executable in the "bin/" folder, so to execute it:

```shell
cd bin 
./testMEventList
cd ..
```

The output should be similar to:

```shell
!../examples/testFile.fits
../data/Output_flute.root
../data/20150313_05041823_Q_CrabNebula-W0.40+035.root
A total amount of 11940 events are classified as Gammas.
Are you serious? It works!!! :O
```
The generated output file is located at "examples/testFile.fits".

The file already contains the 5 required columns within the EVENTS HDU extracted from the "data/20150313_05041823_Q_CrabNebula-W0.40+035.root" CrabNebula melibea file Giovanna provided, with 11940 awesome "gamma-like" events.


Documentation
--------------

MAGIC DL3 Working Group Wikipage:
 
* http://wiki.magic.pic.es/index.php/MAGIC_DL3_WG

Git tutorials (most of them are written for Github, but the workflow is almost identical to the GitLab's one):

* The [GitLab Documentation] (http://goo.gl/mAoqA2), especially the parts "Collaborate" and "Prioritize"
* Interactive Git tutorial: https://try.github.io 
* The infamous 'Hello World!' tutorial for creating your first git project: https://guides.github.com/activities/hello-world/
* Official beginners guide to the Github workflow: https://guides.github.com/introduction/flow/
