# 1. Load master_config into master dictionary
# 2. Identify slave configs (configs referred in master_config), and load them into master dictionary
# 3. Enrich master dictionary, e.g with additional neuron information
# 4. Additional modifications?
# 5. Export network-config.json

class ConfigParser:

    def __init__(self, path=None, snudda_data=None, parse=True):

        self.path = path
        self.network_info = dict()  # master dictionary
        self.snudda_data = snudda_data

        if path is not None and parse:
            self.parse(path)

    def parse_config(self, path=None):

        if path is not None:
            self.path = path

        pass

    def write_config(self, output_path=None):

        pass


