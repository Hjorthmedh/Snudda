from neuron import h

class Spine:

    def __init__(self,
                 parent_section,
                 parent_x,
                 name,
                 head_length,
                 head_diameter,
                 head_axial_resistance,
                 neck_length,
                 neck_diameter,
                 neck_axial_resistance,
                 membrane_capacitance=1.0,
                 parameter_list_neck=tuple(),
                 parameter_list_head=tuple(),
                 neck_mechanism_list = tuple(),
                 head_mechanism_list = tuple()
                 ):

        # Example:
        # mechanism_list = ["pas", "cal13_ms.mod", "cadyn_spine"]
        # parameter_list = [("g_pas", 1.25e-5), ("e_pas", -85)]

        self.parent_section = parent_section
        self.name = name

        self.neck = self.create_spine_neck(length=neck_length,
                                           diameter=neck_diameter,
                                           axial_resistance=neck_axial_resistance,
                                           membrane_capacitance=membrane_capacitance,
                                           mechanism_list=neck_mechanism_list,
                                           parameter_list=parameter_list_neck)

        self.head = self.create_head(spine_neck=self.neck,
                                     length=head_length,
                                     diameter=head_diameter,
                                     axial_resistance=head_axial_resistance,
                                     membrane_capacitance=membrane_capacitance,
                                     mechanism_list=head_mechanism_list,
                                     parameter_list=parameter_list_head)

        # Connecting the spine to the parent compartment
        # Parameters: sec, parent_x, child_x
        self.neck.connect(self.parent_section(parent_x), 0)


    def create_spine_neck(self,
                          length,
                          diameter,
                          axial_resistance,
                          membrane_capacitance,
                          mechanism_list,
                          parameter_list):

        name_sec = self.name + "_neck"
        neck = h.Section(name=name_sec)
        neck.nseg = 1
        neck.L = length * 1e6  # SI -> Natural units in NEURON
        neck.diam = diameter * 1e6
        neck.Ra = axial_resistance
        neck.cm = membrane_capacitance

        for mech in mechanism_list:
            neck.insert(mech)

        for par_name, par_value in parameter_list:
            for seg in neck:
                setattr(seg, par_name, par_value)

        return neck

    def create_head(self,
                    spine_neck,
                    length,
                    diameter,
                    axial_resistance,
                    membrane_capacitance,
                    mechanism_list,
                    parameter_list):

        name_sec = self.name + "_head"
        head = h.Section(name=name_sec)

        head.nseg = 1
        head.L = length * 1e6
        head.diam = diameter * 1e6
        head.Ra = axial_resistance
        head.cm = membrane_capacitance

        for mech in mechanism_list:
            head.insert(mech)

        for par_name, par_value in parameter_list:
            for seg in head:
                setattr(seg, par_name, par_value)

        # Connect spine head to the neck
        head.connect(spine_neck(1), 0)

        return head


