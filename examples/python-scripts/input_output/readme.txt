First clone and install snudda branch input_output_analysis.
git clone -b input_output_analysis. https://github.com/Hjorthmedh/Snudda.git
cd Snudda
pip install -r requirements.txt
pip install -e .[dev]

In the branch I've included Snudda/examples/20220930_3, which is same morphologies and parameters as BasalGangliaData/Parkinson/Bo2022
Dont forget to compile mechanisms, for example:
cd C:\Users\bo.bekkouche\PycharmProjects\input_output01\Snudda\examples\python-scripts\input_output
nrnivmodl C:\Users\bo.bekkouche\PycharmProjects\input_output01\Snudda\examples\20220930_3\mechanisms


for current injection vs output experiment run:
currinj_swap.py
curinj_simulate_pd0.py
curinj_simulate_pd2.py
plot_curinj.py

for spike input frequency vs output frequency run:
input_spikes_vs_output.py


