import parameter_transfer as pt
import morphology_transfer as mt
import mechanisms_transfer as mecht
import selected_models as sm
import creation_of_meta as wm
from argparse import ArgumentParser, RawTextHelpFormatter
import os
from math import * 


#########################################################################################
parser = ArgumentParser("Create meta.json and parameters.json", formatter_class=RawTextHelpFormatter)   
parser.add_argument("--source", help="Parameters source")
parser.add_argument("--morphSource", help="Source neurons path")
parser.add_argument("--destination", help="Destination neurons path")     
args = parser.parse_args()
BG_data = os.environ["SNUDDA_DATA"]
#########################################################################################

print("Starting bgmod2bgdata")
pt.combine_hall_of_fame_with_optimisation_parameters(args.source,args.destination,selected=False)
mt.transfer_morphologies(source=args.morphSource, destination=args.destination, selected=False)
mecht.transfer_mechanisms(source=args.source, destination=args.destination)
sm.select_validated_models(source=args.source, destination=args.destination)
wm.write_meta(args.destination, selected=False)