from snudda.plotting.Network_plot_spike_raster import NetworkPlotSpikeRaster

fileName = '/afs/pdc.kth.se/home/j/jofn/Snudda/networks/TegnerRun.895338/simulation/network-output-spikes-LTS-ablation.csv'
networkFile = '/afs/pdc.kth.se/home/j/jofn/Snudda/networks/TegnerRun.895338/network-pruned-synapses.hdf5'
endTime = 0.1
npsr = NetworkPlotSpikeRaster(fileName,networkFile,skipTime=0.0,
                                  endTime=endTime,
                                  typeOrder=["FSN","dSPN","LTS","iSPN","ChIN"])
