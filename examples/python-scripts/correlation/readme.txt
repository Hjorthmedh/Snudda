First clone and install snudda branch input_output_analysis.
git clone -b input_output_analysis. https://github.com/Hjorthmedh/Snudda.git
cd Snudda
pip install -r requirements.txt
pip install -e .[dev]

In the branch I've included Snudda/examples/20220930_3, which is same morphologies and parameters as BasalGangliaData/Parkinson/Bo2022
Dont forget to compile mechanisms, for example:
cd C:\Users\bo.bekkouche\PycharmProjects\input_output01\Snudda\examples\python-scripts\correlation
nrnivmodl C:\Users\bo.bekkouche\PycharmProjects\input_output01\Snudda\examples\20220930_3\mechanisms


Correlation experiment.
Make sure to match network folder names highlighted in each file (network_name_pd0 and network_name_pd2).
Then run:
swap_correlation.py
run_corr_sim.sh
plotSwapActivity.py

