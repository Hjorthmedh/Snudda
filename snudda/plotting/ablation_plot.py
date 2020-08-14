from snudda.plotting.compare_frequency_activity import ComparePlotTraces

fileNames = ['/afs/pdc.kth.se/home/j/jofn/Snudda/networks/TegnerRun.895338/simulation/volt-LTS-ablation.csv','/afs/pdc.kth.se/home/j/jofn/Snudda/networks/TegnerRun.895338/simulation/volt-FS-ablation.csv','/afs/pdc.kth.se/home/j/jofn/Snudda/networks/TegnerRun.895338/simulation/volt-SPN-ablation.csv','/afs/pdc.kth.se/home/j/jofn/Snudda/networks/TegnerRun.895338/simulation/volt-fullnetwork.csv']

networkFile = '/afs/pdc.kth.se/home/j/jofn/Snudda/networks/TegnerRun.895338/network-pruned-synapses.hdf5'


fileNames = ['/home/jofrony/Desktop/TegnerRun.895338/simulation/volt-LTS-ablation.csv',\
            '/home/jofrony/Desktop/TegnerRun.895338/simulation/volt-FS-ablation.csv',\
            '/home/jofrony/Desktop/TegnerRun.895338/simulation/volt-SPN-ablation.csv',\
            '/home/jofrony/Desktop/TegnerRun.895338/simulation/volt-fullnetwork.csv']
networkFile = '/home/jofrony/Desktop/TegnerRun.895338/network-pruned-synapses.hdf5'

plotLabels = ['LTS ablation', 'FS ablation', 'SPN ablation', 'Full network']

colour = ['purple','red','blue','black']

Network = ComparePlotTraces(fileNames=fileNames,networkFiles=networkFile,labels=plotLabels,colours=colour)

plotOffset = 0 
skipTime = 0 
nTracesMax = 100

Network.plotTraceNeuronType(neuronType="dSPN",nTraces=nTracesMax,offset=plotOffset,skipTime=skipTime)
Network.plotTraceNeuronType(neuronType="iSPN",nTraces=nTracesMax,offset=plotOffset,skipTime=skipTime)

