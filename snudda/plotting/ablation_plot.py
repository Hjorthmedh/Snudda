from snudda.plotting.Compare_network_activity import CompareNetwork

fileNames = ['/afs/pdc.kth.se/home/j/jofn/Snudda/networks/TegnerRun.895338/simulation/network-output-spikes-LTS-ablation.csv',\
            '/afs/pdc.kth.se/home/j/jofn/Snudda/networks/TegnerRun.895338/simulation/network-output-spikes-FS-ablation.csv',\
            '/afs/pdc.kth.se/home/j/jofn/Snudda/networks/TegnerRun.895338/simulation/network-output-spikes-SPN-ablation.csv',\
            '/afs/pdc.kth.se/home/j/jofn/Snudda/networks/TegnerRun.895338/simulation/network-output-spikes-fullnetwork.csv']
networkFiles = '/afs/pdc.kth.se/home/j/jofn/Snudda/networks/TegnerRun.895338/network-pruned-synapses.hdf5'
endTime = 0.1

experiment = ['LTS ablation', 'FS ablation', 'SPN ablation', 'Full network']
CompareNetworkactivity = CompareNetwork(fileNames=fileNames,networkFiles = networkFiles,CompareNeuronType = ["dSPN","iSPN"],skipTime=0.0,
                                  endTime=endTime,
                                  typeOrder=["FSN","dSPN","LTS","iSPN","ChIN"])
