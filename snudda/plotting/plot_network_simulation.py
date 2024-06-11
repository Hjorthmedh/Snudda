#Plot control file. Currently used by batch scripts in examples/parallel.
#Place this in snudda/plotting at a later stage
import sys
from snudda.plotting import PlotSpikeRaster
from snudda.plotting import SnuddaPlotSpikeRaster2

from snudda.plotting import PlotTraces
import matplotlib.pyplot as plt
import os

def plot_network_simulation(nametag,doHistRast=True,doTraces=True,doInpRast=True,inputFileName="input-spikes.hdf5",input_type='',neuronType='dSPN'):
    print('Starting plotting network simulation')
    print(nametag)
    print("..............")
    networkName=nametag
    network_path = os.path.join("networks", networkName)
    
    ##############################################################
    if doHistRast:#tes2
        print("Starting SnuddaPlotSpikeRaster2")
        type_order = ["chin", "dspn", "lts", "ispn", "fs", "fsn"]
        type_division = [[ "dspn", "ispn"],["chin", "fs", "fsn", "lts"]]
        sp = SnuddaPlotSpikeRaster2(network_path=network_path)
        #sp.plot_hist_raster(type_order=type_order, fig_size=(5, 10))
        #sp.plot_hist_raster(type_order=type_order,type_division=type_division, fig_size=(5, 10))
        sp.plot_hist_raster(type_division=type_division)   
    
    ##############################################################
    if doTraces:
        print("Starting plot_traces")
        network_file = os.path.join(network_path,"network-synapses.hdf5")
        network_output = os.path.join(network_path,"simulation", "output.hdf5")
        pt=PlotTraces(network_output, network_file=network_file)
        pt.plot_traces_sep(folder_name=input_type)
        pt.plot_traces(fig_name="traces.pdf")
    
    ##############################################################
    if doInpRast:
        print("Starting plot_input")
        from snudda.plotting import PlotInput
        input_file = os.path.join(network_path, inputFileName)
        print("input_file:")
        print(input_file)
        spi = PlotInput(input_file)
        spi.plot_input(neuronType, 10, fig_size=(15, 5))
        plt.xlim([0,10])
        plt.savefig(os.path.join(network_path,"figures",f"inputspikes_{input_type}"), dpi=300)
        print("asd")
    
    ##############################################################
if __name__ == '__main__':
    from argparse import ArgumentParser, RawTextHelpFormatter
    parser = ArgumentParser("Plot network", formatter_class=RawTextHelpFormatter)
    parser.add_argument("--nametag")
    parser.add_argument("--doHistRast",type=int, default=1)
    parser.add_argument("--doTraces",type=int, default=1)
    parser.add_argument("--doInpRast",type=int, default=1)
    parser.add_argument("--inputFileName",type=str, default="input-spikes.hdf5")
    parser.add_argument("--input_type",type=str, default='')
    parser.add_argument("--neuronType",type=str, default='dSPN')
    args = parser.parse_args()
    plot_network_simulation(args.nametag,doHistRast=args.doHistRast,doTraces=args.doTraces,doInpRast=args.doInpRast,inputFileName=args.inputFileName,input_type=args.input_type,neuronType=args.neuronType)

