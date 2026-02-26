import h5py
import numpy as np

with h5py.File("../topology100/network-synapses.hdf5", "r+") as f:
    neurons = f["network/neurons"]

    if "reaction_diffusion_file" in neurons:
        print("Already exists, skipping.")
    else:
        # determine number of neurons
        n = neurons["neuron_id"].shape[0]

        dt = h5py.string_dtype(encoding="utf-8")

        neurons.create_dataset(
            "reaction_diffusion_file",
            shape=(n,),
            dtype=dt
        )

        neurons["reaction_diffusion_file"][:] = ""

        print(f"Created reaction_diffusion_file for {n} neurons.")            
