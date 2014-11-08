

# Compute the peptide sequence lengths, which is harder than taking the len()
# of the sequence string. Must exclude all modifications in "()" and the glycan
# composition at the end with "[]".
def get_sequence_length(data):
    data["peptideLens"] = data.Peptide.str.split(
        r'\(.*?\)').str.join('').str.len()
    return data


def extract_modifications(sequence):
    state = "aa"
    mods = []
    glycan = ""
    current_aa = ""
    current_mod = ""
    paren_level = 0
    i = 0
    while i < len(sequence):
        next_char = sequence[i]
        if next_char is "(":
            if state == "aa":
                state = "mod"
                assert paren_level == 0
                paren_level += 1
            elif state == "mod":
                paren_level += 1
                current_mod += next_char
        elif next_char == ")":
            if state == "aa":
                raise Exception("Invalid Sequence. ) found outside of modification.")
            elif state == "mod":
                paren_level -= 1
                if paren_level == 0:
                    state = 'aa'
                    mods.append(current_mod)
                    current_mod = ""
                    current_aa = ""
                else:
                    current_mod += next_char
        elif next_char == "[":
            glycan = sequence[i:]
            break
        elif state == "aa":
            if(current_aa != ""):
                current_mod = ""
                current_aa = ""
            current_aa += next_char
        elif state == "mod":
            current_mod += next_char
        else:
            raise Exception("Unknown Tokenizer State", current_aa, current_mod, i, next_char)
        i += 1
    if current_mod != "":
        mods.append(current_mod)
    mods = sorted(mods)

    return ','.join(mods) + glycan


def get_modification_signature(data_struct):
    mod_data = data_struct.Glycopeptide_identifier.apply(extract_modifications)
    data_struct["modificationSignature"] = mod_data
    return data_struct
