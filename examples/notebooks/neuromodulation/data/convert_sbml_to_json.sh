# python ../../../../snudda/utils/sbml_to_snudda.py SBML/MODEL_speedy_reduced2.xml JSON/robert_reaction_diffusion.json

# nM --> mM (in Snudda, derived from M/m3 = mM/litre)
# python ../../../../snudda/utils/sbml_to_snudda.py SBML/Robert-MODEL_speedy_reduced2_UPDATED.xml JSON/reaction_diffusion_D1.json --conc_scale_factor 1e-6

python ../../../../snudda/utils/sbml_to_snudda.py SBML/Robert-MODEL_speedy_reduced2_UPDATED.xml JSON/reaction_diffusion_D1_with_target.json --conc_scale_factor 1e-6

python ../../../../snudda/utils/sbml_to_snudda.py SBML/Robert-MODEL_speedy_reduced2_UPDATED_notarget.xml JSON/reaction_diffusion_D1.json --conc_scale_factor 1e-6

# python ../../../../snudda/utils/sbml_to_snudda.py SBML/Nair_2016_optimized_UPDATED.xml JSON/reaction_diffusion_D1.json

# python ../../../../snudda/utils/sbml_to_snudda.py SBML/Nair2015-D1-BIOMD0000000635_url.xml JSON/reaction_diffusion_D1.json

python ../../../../snudda/utils/sbml_to_snudda.py SBML/Nair2015-D2-BIOMD0000000636_url_UPDATED.xml JSON/reaction_diffusion_D2.json --conc_scale_factor 1e-6

#python ../../../../snudda/utils/sbml_to_snudda.py SBML/Yapo2017-D1-MODEL1701170000_url.xml JSON/reaction_diffusion_D1.json

# python ../../../../snudda/utils/sbml_to_snudda.py SBML/Yapo2017-D2-MODEL1701170001_url.xml JSON/reaction_diffusion_D2.json

# python ../../../../snudda/utils/sbml_to_snudda.py SBML/Robert-MODEL_speedy_reduced2.xml JSON/robert_reaction_diffusion.json



