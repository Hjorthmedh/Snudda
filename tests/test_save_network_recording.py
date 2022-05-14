import unittest
import os
import neuron


from snudda.simulate.save_network_recording import *

class TestSaveNetworkRecording(unittest.TestCase):

    def testSavingTimestep(self):

        pass

    def testNeuronRecordings(self):

        test_vector = neuron.h.Vector([1, 2, 3, 4])
        test_recordings = NeuronRecordings(neuron_id=12345)
        test_recordings.register_compartment_data(data=test_vector,
                                                  data_type="voltage",
                                                  sec_id=0,
                                                  sec_x=0.5)
        test_recordings.register_synapse_data(data=test_vector,
                                              data_type="synaptic_current",
                                              synapse_type="tmGlut",
                                              presynaptic_id=0,
                                              sec_id=1,
                                              sec_x=0.1,
                                              cond=1e-3)
        test_recordings.register_spike_data(data=test_vector,
                                            sec_id=0,
                                            sec_x=0.5)

        self.assertTrue("voltage" in test_recordings.data)
        self.assertTrue("synaptic_current" in test_recordings.data)
        self.assertTrue("spikes" in test_recordings.data)

        self.assertTrue(isinstance(test_recordings.data["voltage"], CompartmentData))
        self.assertTrue(isinstance(test_recordings.data["synaptic_current"], SynapseData))
        self.assertTrue(isinstance(test_recordings.data["spikes"], SpikeData))





