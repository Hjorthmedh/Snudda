import json
import numpy as np
import matplotlib.pyplot as plt


param_key = "abc"

with open("data/modulation-bevan2020.json", "rt") as f:
    par_data = json.load(f)[param_key]

param = {}
for d in par_data:
    param[d["param_name"]] = d["value"]
    

channels = np.unique(["_".join(x.split("_")[-2:]) for x in param.keys()])

plt.figure()

for chan in channels:

    if f"mod_pka_g_min_{chan}" in param:
        mod_min = param[f"mod_pka_g_min_{chan}"]
        mod_max = param[f"mod_pka_g_max_{chan}"]
        half_activation = param[f"mod_pka_g_half_{chan}"]
        hill_coefficient = param[f"mod_pka_g_hill_{chan}"]
    elif f"mod_pka_p_min_{chan}" in param:
        mod_min = param[f"mod_pka_p_min_{chan}"]
        mod_max = param[f"mod_pka_p_max_{chan}"]
        half_activation = param[f"mod_pka_p_half_{chan}"]
        hill_coefficient = param[f"mod_pka_p_hill_{chan}"]
    else:
        print("Skipping {chan}")
        continue
        
    conc = np.arange(0, half_activation*3, half_activation/100)
    modulation = mod_min + np.divide((mod_max-mod_min) * np.power(conc, hill_coefficient),
                                     (np.power(conc, hill_coefficient) + np.power(half_activation, hill_coefficient)))


    plt.plot(conc*1e6, modulation, label=chan)

plt.xlabel("Concentration (nM)")
plt.ylabel("Modulation")
plt.title("Modulation")
plt.legend(loc="upper right")
plt.ion()
plt.show()
plt.savefig(f"ion_channel_modulation_curves.png")


input("Press a key")
