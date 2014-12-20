# Required information for generating peptides that can be cleaved by a given protease
enzymes = {
    "trypsin": {
        "cleavage_start": [""],
        "cleavage_end": ["K", "R"],
        "name": "trypsin"
    }
}


def get_enzyme(name):
    return enzymes[name.lower()]
