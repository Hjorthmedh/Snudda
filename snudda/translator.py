
translation ={
    "dendrite" : "basal",
    "soma" : "soma",
    "axon" : "axon"
}

re_translation ={
    "dend" : "dendrite",
    "soma" : "soma",
    "axon" : "axon"
}


conductance = {
	"naf_ms": "gbar",
	"kir_ms": "gbar",
	"kas_ms": "gbar",
	"kaf_ms": "gbar",
	"cal12_ms": "pbar",
	"cal13_ms": "pbar",
	"can_ms": "pbar",
	"car_ms": "pbar"
}


def translate(word):

    return translation[word]

    
def re_translate(word):

    return re_translation[word]

def conductance_translate(word):

    return conductance[word]
