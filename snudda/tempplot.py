from snudda.plotting.Compare_network_activity import CompareNetwork

fileNames = ['/home/jofrony/Desktop/TegnerRun.895338/simulation/network-output-spikes-LTS-ablation.csv',\
            '/home/jofrony/Desktop/TegnerRun.895338/simulation/network-output-spikes-FS-ablation.csv',\
            '/home/jofrony/Desktop/TegnerRun.895338/simulation/network-output-spikes-SPN-ablation.csv',\
            '/home/jofrony/Desktop/TegnerRun.895338/simulation/network-output-spikes-fullnetwork.csv']
networkFiles = '/home/jofrony/Desktop/TegnerRun.895338/network-pruned-synapses.hdf5'
endTime = 0.1

experiment = ['LTS ablation', 'FS ablation', 'SPN ablation', 'Full network']
CompareNetworkactivity = CompareNetwork(fileNames=fileNames,networkFiles = networkFiles,CompareNeuronType = ["dSPN","iSPN"],skipTime=0.0,
                                  endTime=endTime,
                                  typeOrder=["FSN","dSPN","LTS","iSPN","ChIN"])
