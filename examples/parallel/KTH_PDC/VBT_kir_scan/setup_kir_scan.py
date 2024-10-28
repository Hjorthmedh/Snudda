from snudda import Snudda
from argparse import ArgumentParser, RawTextHelpFormatter


if __name__ == "__main__":

    parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument("network_path", type=str)
    parser.add_argument("--ipython_profile", default=None)
    args = parser.parse_args()

    if args.ipython_profile:
        parallel_flag=True
    else:
        parallel_flag=False
    
    snd = Snudda(network_path=args.network_path, parallel=parallel_flag,
                 ipython_profile=args.ipython_profile)
    
    snd.import_config("config/network.json", overwrite=True)
    snd.create_network(honor_morphology_stay_inside=False)

    snd.setup_input(input_config="input.json")
    
