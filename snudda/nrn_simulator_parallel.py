import neuron
import bluepyopt.ephys as ephys
from bluepyopt.ephys.simulators import NrnSimulatorException


class NrnSimulatorParallel(ephys.simulators.NrnSimulator):

    def __init__(self, dt=None, cvode_active=False):

        self.disable_banner = True
        self.banner_disabled = False

        self.neuron.h.load_file('stdrun.hoc')

        self.dt = dt if dt is not None else self.neuron.h.dt

        self.cvode_active = cvode_active
        if cvode_active:
            self.neuron.h.cvode_active(1)
            # assert False, "Make sure Cvode works"
        else:
            self.neuron.h.cvode_active(0)

        self.pc = self.neuron.h.ParallelContext()
        self.pc.set_maxstep(10)

    def run(self, tstop=None, dt=None):

        self.neuron.h.tstop = tstop

        if dt is None:
            dt = self.dt

        self.neuron.h.dt = dt
        self.neuron.h.steps_per_ms = 1.0 / dt

        try:
            # self.neuron.h.stdinit()
            self.neuron.h.finitialize()

            self.pc.psolve(tstop)
            self.pc.barrier()

            # Using pc.barrier instead of runworker and done (which doesn't terminate)
            # https://www.neuron.yale.edu/phpBB/viewtopic.php?f=6&t=2781
            # self.pc.runworker()
            # self.pc.done()

        except Exception as e:
            raise NrnSimulatorException('Neuron simulator error', e)
