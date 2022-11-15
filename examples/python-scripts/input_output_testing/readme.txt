First clone and install snudda branch input_output_analysis.
git clone -b input_output_analysis. https://github.com/Hjorthmedh/Snudda.git
cd Snudda
pip install -r requirements.txt
pip install -e .[dev]

In the branch I've included Snudda/examples/20220930_3, which is same morphologies and parameters as BasalGangliaData/Parkinson/Bo2022
Dont forget to compile mechanisms, for example:
cd C:\Users\bo.bekkouche\PycharmProjects\input_output01\Snudda\examples\python-scripts\input_output
nrnivmodl C:\Users\bo.bekkouche\PycharmProjects\input_output01\Snudda\examples\20220930_3\mechanisms


for current injection vs output experiment including normalized current injection run:
currinj_swap_norm.py
curinj_simulate_pd0.py
curinj_simulate_pd2.py
curinj_simulate_pd2_norm.py
plot_curinj_norm.py

for current injection vs output experiment including normalized current injection run with parallel simulation (currently not working):
currinj_swap_norm.py
sh .\run_currinj_sim.sh
plot_curinj_norm.py