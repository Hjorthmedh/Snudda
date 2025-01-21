import bluepyopt.ephys as ephys


class NrnFileMorphology_axon_fix(ephys.morphologies.NrnFileMorphology):

    def __init__(self,
                 morphology_path,
                 axon_length,
                 axon_diameter,
                 do_set_nseg=True,
                 comment='',
                 replace_axon_hoc=None,
                 axon_nseg_frequency=40,
                 nseg_frequency=40,
                 replace_axon_myelin_length=None,
                 replace_axon_myelin_diameter=None,
                 morph_modifiers=None,
                 morph_modifiers_hoc=None,
                 morph_modifiers_kwargs=None):

        super().__init__(morphology_path=morphology_path,
                         do_replace_axon=False,
                         do_set_nseg=do_set_nseg,
                         comment=comment,
                         replace_axon_hoc=replace_axon_hoc,
                         axon_nseg_frequency=axon_nseg_frequency*1e6,
                         nseg_frequency=nseg_frequency,
                         morph_modifiers=morph_modifiers,
                         morph_modifiers_hoc=morph_modifiers_hoc,
                         morph_modifiers_kwargs=morph_modifiers_kwargs)

        self.axon_diameter = axon_diameter
        self.axon_length = axon_length
        self.replace_axon_myelin_length = replace_axon_myelin_length
        self.replace_axon_myelin_diameter = replace_axon_myelin_diameter

    def instantiate(self, sim=None, icell=None):

        super().instantiate(sim=sim, icell=icell)

        self.replace_axon2(sim=sim, icell=icell)
        # Call the new replace_axon

    def replace_axon2(self, sim=None, icell=None):
        """Replace axon (length and diameter)"""

        for section in icell.axonal:
            sim.neuron.h.delete_section(sec=section)

        # Create new axon array
        sim.neuron.h.execute(f"create axon[{len(self.axon_length)}]", icell)

        assert len(self.axon_length) == len(self.axon_diameter)

        for index, (section, ax_len, ax_dia) in enumerate(zip(icell.axon, self.axon_length, self.axon_diameter)):
            section.L = ax_len * 1e6  # Convert to micrometers
            section.nseg = 1 + 2 * int(section.L / self.axon_nseg_frequency)
            section.diam = ax_dia * 1e6  # Convert to micrometers
            icell.axonal.append(sec=section)
            icell.all.append(sec=section)

        icell.axon[0].connect(icell.soma[0], 1.0, 0.0)
        for index in range(1, len(self.axon_length)):
            icell.axon[index].connect(icell.axon[index-1], 1.0, 0.0)

        print(f"replace_axon2: {self.axon_length =}, {self.axon_diameter =}")

        if self.replace_axon_myelin_length is not None:
            print(f"Adding myelinated section to axon: {self.replace_axon_myelin_length =}, "
                  f"{self.replace_axon_myelin_diameter =}")

            for section in icell.myelinated:
                sim.neuron.h.delete_section(sec=section)

            sim.neuron.h.execute(f"create myelin[{len(self.replace_axon_myelin_length)}]", icell)

            assert len(self.replace_axon_myelin_length) == len(self.replace_axon_myelin_diameter)

            for index, (section, ax_len, ax_dia) in enumerate(zip(icell.myelin,
                                                                  self.replace_axon_myelin_length,
                                                                  self.replace_axon_myelin_diameter)):
                section.L = ax_len * 1e6  # Convert to micrometers
                section.nseg = 1 + 2 * int(section.L / self.axon_nseg_frequency)
                section.diam = ax_dia * 1e6  # Convert to micrometers
                icell.myelinated.append(sec=section)
                icell.all.append(sec=section)

            icell.myelin[0].connect(icell.axon[len(self.axon_length)-1], 1.0, 0.0)
            for index in range(1, len(self.replace_axon_myelin_length)):
                icell.myelin[index].connect(icell.myelin[index - 1], 1.0, 0.0)

