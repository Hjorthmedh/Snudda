
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


def translate(word):

    return translation[word]

    
def re_translate(word):

    return re_translation[word]

