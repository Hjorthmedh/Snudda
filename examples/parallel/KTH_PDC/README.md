# Setting up Snudda on Dardel / PDC

First clone Snudda to your home directory:
```
git clone https://github.com/Hjorthmedh/Snudda.git
```

Run the ```Dardel_install_Snudda.sh``` script in ```Snudda/examples/parallel/KTH_PDC```

```
./Dardel_install_Snudda.sh
```

To install ```NEURON``` you need to download Miniconda and the NEURON source code.

In the ```Snudda/examples/parallel/KTH_PDC``` folder run:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod u+x Miniconda3-latest-Linux-x86_64.sh
```

```
pushd ~/
mkdir local/dardel
cd local/dardel
git clone https://github.com/neuronsimulator/nrn -b 8.0.0
popd
```

Now run the installation script:

```
sbatch Dardel_NEURON_install.job
```

This will create files in your home directory ```local/dardel``` subdirectory.


# Setting up Snudda on Beskow

First create your miniconda environment on Beskow:

Since git does not work on the compute nodes, first manually download Neuron:

```
mkdir /cfs/klemming/nobackup/${USER:0:1}/$USER/local/beskow
pushd /cfs/klemming/nobackup/${USER:0:1}/$USER/local/beskow
git clone https://github.com/neuronsimulator/nrn -b 8.0.0
popd
```

Then to install Snudda and to compile NEURON run:

```
sbatch Beskow_install.job
```

Setup a symbolic link to your favourite mechanisms from the ```Snudda/examples/parallel```directory:

```
ln -s /cfs/klemming/nobackup/"${USER:0:1}"/$USER/Snudda/snudda/data/neurons/mechanisms
```

Then to run a Snudda simulation:

```
sbatch Beskow_simulate.job
```

If you need to change the installation path, then update the line: 
```
L=/cfs/klemming/nobackup/${USER:0:1}/$USER/local/$SNIC_RESOURCE
```
which exists in ```Beskow_install.sh```, ```Miniconda_install.sh``` and ```activate_miniconda.sh```.


# Setting up Snudda on Tegner

Log on to Tegner:
```
kinit --forwardable hjorth@NADA.KTH.SE
ssh -X -o GSSAPIDelegateCredentials=yes -o GSSAPIKeyExchange=yes -o GSSAPIAuthentication=yes hjorth@tegner.pdc.kth.se
```

Get the latest version of Snudda:
```
git clone https://github.com/Hjorthmedh/Snudda.git
```

or

```
git clone git@github.com:Hjorthmedh/Snudda.git
```

In the ```Snudda/examples/parallel``` folder run:
```
./Tegner_install.sh
```

If you need to change the installation path, then you update the line: 
```
L=/cfs/klemming/nobackup/${USER:0:1}/$USER/local/$SNIC_RESOURCE
```
which exists in ```Tegner_install.sh```, ```Miniconda_install.sh```, ```activate_miniconda.sh```, ```Tegner_runSnudda.job```.

## Test your Snudda installation on Tegner

To generate a test network, in the ```Snudda/examples/parallel``` folder run:

```
sbatch Tegner_runSnudda.job
```

You can check the queueing of your job using 

```squeue -u <yourusername>```

When the job starts, Snudda log files will appear after a few minutes here:

```ls -lrt ../../networks/TegnerNetwork/log```

You can also find slurm log files in the ```Snudda/examples/parallel/save/``` folder


If you for some reason want to start the conda environment on the local login node:
```
source activate_miniconda.sh
conda activate
```
