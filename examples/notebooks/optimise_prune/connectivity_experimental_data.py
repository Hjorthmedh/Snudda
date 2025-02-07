all_experimental_data = {}

all_experimental_data["dSPN", "iSPN"] = [(0, 50e-6, 3/47.0), (0, 100e-6, 3/66.0)]
all_experimental_data["dSPN", "dSPN"] = [(0, 50e-6, 5/19.0), (0, 100e-6, 3/43.0)]
all_experimental_data["iSPN", "iSPN"] = [(0, 50e-6, 14/39.0), (0, 100e-6, 7/31.0)]
all_experimental_data["iSPN", "dSPN"] = [(0, 50e-6, 13/47.0), (0, 100e-6, 10/80.0)]

all_experimental_data["FS", "FS"] = [(0, 250e-6, 7/12.0)]
all_experimental_data["FS", "iSPN"] = [(0, 100e-6, 6/9.0), (0, 150e-6, 21/54.0), (0, 250e-6, 27/77.0)]
all_experimental_data["FS", "dSPN"] = [(0, 100e-6, 8/9.0), (0, 150e-6, 29/48.0), (0, 250e-6, 48/90.0)]
all_experimental_data["FS", "LTS"] = [(0, 250e-6, 2/12.0)]

# NPY population: 79% LTS, 21% NGF according to Ibanes-Sandoval 2011
# Optogenetical stimulation of FS population, 62.5% of LTS population react
# Old ref (which one?) FS -> LTS 2/12 pairs connected
# Optogenetical stimulation of FS, yields 9/9 NGF activated
# -- but we dont know how many FS on average connect
#
# (1 - (1 - 2/12) ** N_FS) = 62.5% of connection to LTS
# --> 5 eller 6 FS celler
#
# (1 - (1 - 0.4) ** 5) = 92% (or higher)
# Ie, we need 40% or higher connectivity FS -> NGF to get 9/9 NGF connected
#
# We implicitly assume LTS and NGF are approximately the same dendritic size
#

all_experimental_data["FS", "NGF"] = [(0, 150e-6, 0.4)]   # Se comment above

all_experimental_data["LTS", "dSPN"] = [(0, 250e-6, 2/60.0)]
all_experimental_data["LTS", "iSPN"] = [(0, 250e-6, 2/60.0)]

all_experimental_data["ChIN", "dSPN"] = [(0, 250e-6, 0.05)]
all_experimental_data["ChIN", "iSPN"] = [(0, 250e-6, 0.05)]
all_experimental_data["ChIN", "LTS"] = [(0, 250e-6, 53/72.0)]  # Lou ea 2013
all_experimental_data["ChIN", "FS"] = [(0, 250e-6, 19/26.0)]   # Lou ea 2013
all_experimental_data["ChIN", "NGF"] = [(0, 250e-6, 8/14)]  # English 2012


all_experimental_data["NGF", "dSPN"] = [(0, 100e-6, 25/29.0),  # Ibanez-Sandoval 2011
                                        (0, 250e-6, 11/14.0),  # English 2012
                                        (0, 250e-6, 30/50.0)]  # Luo 2013
all_experimental_data["NGF", "iSPN"] = [(0, 100e-6, 25/29.0),  # Ibanez-Sandoval 2011
                                        (0, 250e-6, 11/14.0),  # English 2012
                                        (0, 250e-6, 30/50.0)]  # Luo 2013
all_experimental_data["NGF", "ChIN"] = [(0, 250e-6, 3/14.0)]  # English 2012

all_experimental_data["TH", "LTS"] = [(0, 250e-6, 0.3)]   # No data, GUESSING 30%
all_experimental_data["TH", "dSPN"] = [(0, 250e-6, 13/44)]   # Luo 2013 (dist?)
all_experimental_data["TH", "iSPN"] = [(0, 250e-6, 13/44)]   # Luo 2013 (dist?)

all_experimental_data["ChIN", "TH"] = [(0, 250e-6, 0.3)]   # No data, GUESSING 30%



