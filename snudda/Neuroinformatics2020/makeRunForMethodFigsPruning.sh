export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=Snudda_LOCAL

ipcluster start -n 6 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

snudda place Neuroinformatics2020/Net10062-var-1/
snudda place Neuroinformatics2020/Net10062-var-2/
snudda place Neuroinformatics2020/Net10062-var-3/
snudda place Neuroinformatics2020/Net10062-var-4/
snudda place Neuroinformatics2020/Net10062-var-5/

snudda detect Neuroinformatics2020/Net10062-var-1/ --volumeID Striatum
snudda prune Neuroinformatics2020/Net10062-var-1/

snudda detect Neuroinformatics2020/Net10062-var-2/ --volumeID Striatum
snudda prune Neuroinformatics2020/Net10062-var-2/

snudda detect Neuroinformatics2020/Net10062-var-3/ --volumeID Striatum
snudda prune Neuroinformatics2020/Net10062-var-3/

snudda detect Neuroinformatics2020/Net10062-var-4/ --volumeID Striatum
snudda prune Neuroinformatics2020/Net10062-var-4/

snudda detect Neuroinformatics2020/Net10062-var-5/ --volumeID Striatum
snudda prune Neuroinformatics2020/Net10062-var-5/


ipcluster stop


python3 analyse_striatum.py Neuroinformatics2020/Net10062-var-1/network-pruned-synapses.hdf5
python3 analyse_striatum.py Neuroinformatics2020/Net10062-var-2/network-pruned-synapses.hdf5
python3 analyse_striatum.py Neuroinformatics2020/Net10062-var-3/network-pruned-synapses.hdf5
python3 analyse_striatum.py Neuroinformatics2020/Net10062-var-4/network-pruned-synapses.hdf5
python3 analyse_striatum.py Neuroinformatics2020/Net10062-var-5/network-pruned-synapses.hdf5

python3 methodsPaperFigure2.py

