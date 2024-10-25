from snudda import Snudda
from argparse import ArgumentParser, RawTextHelpFormatter


if __name__ == "__main__":

    parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument("network_path", type=str)
    args = parser.parse_args()

    snd = Snudda(network_path=args.network_path)
    snd.import_config("config/network.json", overwrite=True)
    snd.create_network(honor_morphology_stay_inside=False)

    snd.setup_input(input_config="input.json")
    
