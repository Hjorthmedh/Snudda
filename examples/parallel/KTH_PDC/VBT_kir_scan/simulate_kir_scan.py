from snudda import Snudda
from argparse import ArgumentParser, RawTextHelpFormatter


def simulate(network_path, output_path, kir_factor=1.0, time=1.0):
    snd = Snudda(network_path=network_path)
    sim = snd.simulate(time=0, output_file=output_path)

    # Next we need to reduce kir by a factor
    print("Reduce the KIR channel conductance in SPN")

    print(f"Before: {sim.neurons[0].icell.soma[0](0.5).kir_ms.gbar = }")
    
    for sec in sim.sim.neuron.h.allsec():
        for seg in sec:
            if sim.sim.neuron.h.ismembrane("kir_ms", sec=sec):
                # Scale the conductance parameter for each segment
                seg.kir_ms.gbar *= kir_factor


    print(f"After: {sim.neurons[0].icell.soma[0](0.5).kir_ms.gbar = }")
    

    sim.run(t=time*1e3)
    sim.write_output()

    

if __name__ == "__main__":

    import sys

    if '-python' in sys.argv:
        print("Network_simulate.py called through nrniv, fixing arguments")
        pythonidx = sys.argv.index('-python')
        if len(sys.argv) > pythonidx:
            sys.argv = sys.argv[pythonidx + 1:]

    parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument("network_path", type=str)
    parser.add_argument("--kir_factor", type=float, default=1.0)
    parser.add_argument("--output", type=str, default=None)
    parser.add_argument("--time", type=float, default=1.0)
    args = parser.parse_args()

    simulate(network_path=args.network_path,
             output_path=args.output,
             kir_factor=args.kir_factor,
             time=args.time)

    
