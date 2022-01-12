import neuron
import bluepyopt.ephys as ephys
from bluepyopt.ephys.simulators import NrnSimulatorException


class NrnSimulatorParallel(ephys.simulators.NrnSimulator):

    def __init__(self, dt=None, cvode_active=False, cvode_minstep=None,
                 random123_globalindex=None, mechanisms_directory=None):
        """Constructor
        Args:
            dt (float): the integration time step used by neuron.
            cvode_active (bool): should neuron use the variable time step
                integration method (Default: False for Snudda)
            cvode_minstep (float): the minimum time step allowed for a cvode
                step. Default is 0.0.
            random123_globalindex (int): used to set the global index used by
                all instances of the Random123 instances of Random
            mechanisms_directory (str): path to the parent directory of the
                directory containing the mod files. If the mod files are in
                "./data/mechanisms", then mechanisms_directory should be
                "./data/".
        """

        super().__init__(dt=dt, cvode_active=cvode_active, cvode_minstep=cvode_minstep,
                         random123_globalindex=random123_globalindex, mechanisms_directory=mechanisms_directory)

        self.pc = self.neuron.h.ParallelContext()
        self.pc.set_maxstep(10)



    def run(
            self,
            tstop=None,
            dt=None,
            cvode_active=None,
            random123_globalindex=None):
        """Run protocol"""

        self.neuron.h.tstop = tstop

        if cvode_active and dt is not None:
            raise ValueError(
                'NrnSimulator: Impossible to combine the dt argument when '
                'cvode_active is True in the NrnSimulator run method')

        if cvode_active is None:
            cvode_active = self.cvode_active

        if not cvode_active and dt is None:  # use dt of simulator
            if self.neuron.h.dt != self.dt:
                raise Exception(
                    'NrnSimulator: Some process has changed the '
                    'time step dt of Neuron since the creation of this '
                    'NrnSimulator object. Not sure this is intended:\n'
                    'current dt: %.6g\n'
                    'init dt: %.6g' % (self.neuron.h.dt, self.dt))
            dt = self.dt

        self.neuron.h.cvode_active(1 if cvode_active else 0)
        if self.cvode_minstep_value is not None:
            save_minstep = self.cvode_minstep
            self.cvode_minstep = self.cvode_minstep_value

        if cvode_active:
            print(f'Running Neuron simulator {tstop:.6g} ms, with cvode')
        else:
            self.neuron.h.dt = dt
            self.neuron.h.steps_per_ms = 1.0 / dt
            print(f'Running Neuron simulator {tstop:.6g} ms, with dt={dt}%r')

        if random123_globalindex is None:
            random123_globalindex = self.random123_globalindex

        if random123_globalindex is not None:
            rng = self.neuron.h.Random()
            rng.Random123_globalindex(random123_globalindex)

        try:
            # self.neuron.h.run()
            self.neuron.h.finitialize()

            self.pc.psolve(tstop)
            self.pc.barrier()

            # Using pc.barrier instead of runworker and done (which doesn't terminate)
            # https://www.neuron.yale.edu/phpBB/viewtopic.php?f=6&t=2781
            # self.pc.runworker()
            # self.pc.done()

        except Exception as e:
            raise NrnSimulatorException('Neuron simulator error', e)

        if self.cvode_minstep_value is not None:
            self.cvode_minstep = save_minstep

        print('Neuron simulation finished')
