# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import sys
from snudda.plotting import PlotSpikeRaster
from snudda.plotting import SnuddaPlotSpikeRaster2

from snudda.plotting import PlotTraces
import matplotlib.pyplot as plt
import os
from snudda.plotting.plot_network_simulation import SnuddaPlotInputTest

def plotNetworkSimulation(nametag,doHistRast=True,doinpFreqVsSpikerate=True,doTraces=True,doInpRast=True,inputFileName="input-spikes.hdf5",inputType='',neuronType='dSPN'):
    print('Starting plotting network simulation')
    print(nametag)
    print("..............")
    #networkName="nosyn_1k_test1"
    #networkName="noinp_1k_test1"
    networkName=nametag#+"_1k_spn01"
    network_path = os.path.join("networks", networkName)
    
    ##############################################################
    if doHistRast:
        print("Starting SnuddaPlotSpikeRaster2")
        type_order = ["chin", "dspn", "lts", "ispn", "fs", "fsn"]
        type_division = [[ "dspn", "ispn"],["chin", "fs", "fsn", "lts"]]
        sp = SnuddaPlotSpikeRaster2(network_path=network_path)
        #sp = PlotSpikeRaster(network_path=network_path)
    
        #sp.plot_hist_raster(type_order=type_order, fig_size=(5, 10))
        #sp.plot_hist_raster(type_order=type_order,type_division=type_division, fig_size=(5, 10))
        sp.plot_hist_raster(type_division=type_division)
        #sp.plot_colour_raster()

    ##############################################################
    #Plot input vs output
    if doinpFreqVsSpikerate:
        print("Starting plot_input_output")
        input_config = os.path.join(network_path,'..','..','input',"external-input-dSTR-scaled-v4-bobek11.json")
        it=SnuddaPlotInputTest(network_path=network_path)
        #it.plot_input_output2(input_config=input_config)
        it.plot_input_output(input_config=input_config)
        
    
    ##############################################################
    if doTraces:
        print("Starting plot_traces")
        network_file = os.path.join(network_path,"network-synapses.hdf5")
        network_output = os.path.join(network_path,"simulation", "network-output.hdf5")
        pt=PlotTraces(network_output, network_file=network_file)
        pt.plot_traces_sep(folderName=inputType)
        pt.plot_traces(fig_name="traces.pdf",numNeurons=1000)
    
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
        plt.savefig(os.path.join(network_path,"figures",f"inputspikes_{inputType}"), dpi=300)
        print("asd")
    
    ##############################################################
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    from argparse import ArgumentParser, RawTextHelpFormatter
    parser = ArgumentParser("Plot network", formatter_class=RawTextHelpFormatter)
    parser.add_argument("--nametag")
    parser.add_argument("--doHistRast",type=int, default=1)
    parser.add_argument("--doinpFreqVsSpikerate",type=int, default=1)
    parser.add_argument("--doTraces",type=int, default=1)
    parser.add_argument("--doInpRast",type=int, default=1)
    parser.add_argument("--inputFileName",type=str, default="input-spikes.hdf5")
    parser.add_argument("--inputType",type=str, default='')
    parser.add_argument("--neuronType",type=str, default='dSPN')
    args = parser.parse_args()
    plotNetworkSimulation(args.nametag,doHistRast=args.doHistRast,doinpFreqVsSpikerate=args.doinpFreqVsSpikerate,doTraces=args.doTraces,doInpRast=args.doInpRast,inputFileName=args.inputFileName,inputType=args.inputType,neuronType=args.neuronType)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
