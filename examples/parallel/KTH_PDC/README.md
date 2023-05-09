# Setting up Snudda on Dardel / PDC

First clone Snudda to your home directory. You might also need to get BasalGangliaData repository (private link):
```
git clone https://github.com/Hjorthmedh/Snudda.git
```

Run the ```Dardel_install_Snudda.sh``` script in ```Snudda/examples/parallel/KTH_PDC```

```
./Dardel_install_Snudda.sh
```

Finally to test that it worked, you can run:

```
sbatch Dardel_runSnudda.job
```

After the network is created, you can start a simulation:

```
sbatch Dardel_simulate.job
```

You can find the generated network files and simulation in ```Snudda/examples/parallel/KTH_PDC/networks/test_10k```.


If you for some reason want to start the conda environment on the local login node:
```
source activate_miniconda.sh
conda activate
```
