from ..utils.memoize import memoize


# Compute the peptide sequence lengths, which is harder than taking the len()
# of the sequence string. Must exclude all modifications in "()" and the glycan
# composition at the end with "[]".
def get_sequence_length(data):
    data["peptideLens"] = data.Peptide.str.split(
        r'\(.*?\)').str.join('').str.len()
    return data


@memoize(10000)
def sequence_tokenizer(sequence):
    state = "aa"
    mods = []
    chunks = []
    glycan = ""
    current_aa = ""
    current_mod = ""
    paren_level = 0
    i = 0
    while i < len(sequence):
        next_char = sequence[i]
        if next_char == "(":
            if state == "aa":
                state = "mod"
                assert paren_level == 0
                paren_level += 1
            elif state == "mod":
                paren_level += 1
                current_mod += next_char
        elif next_char == ")":
            if state == "aa":
                raise Exception(
                    "Invalid Sequence. ) found outside of modification, Position {0}. {1}".format(i, sequence))
            elif state == "mod":
                paren_level -= 1
                if paren_level == 0:
                    state = 'aa'
                    mods.append(current_mod)
                    chunks.append([current_aa, current_mod])
                    current_mod = ""
                    current_aa = ""
                else:
                    current_mod += next_char
        elif next_char == "[" and state == 'aa':
            glycan = sequence[i:]
            break
        elif state == "aa":
            if(current_aa != ""):
                chunks.append([current_aa, current_mod])
                current_mod = ""
                current_aa = ""
            current_aa += next_char
        elif state == "mod":
            current_mod += next_char
        else:
            raise Exception(
                "Unknown Tokenizer State", current_aa, current_mod, i, next_char)
        i += 1
    if current_aa != "":
        chunks.append([current_aa, current_mod])
    if current_mod != "":
        mods.append(current_mod)

    return chunks, mods, glycan


def extract_modifications(sequence):
    sequence_tokens, mods, glycan = sequence_tokenizer(sequence)
    mods = sorted(mods)
    return ','.join(mods) + glycan


def modification_signature(data_struct):
    mod_data = data_struct.Glycopeptide_identifier.apply(extract_modifications)
    data_struct["modificationSignature"] = mod_data
    return data_struct


def glycosites(sequence):
    tokens, mods, glycan = sequence_tokenizer(sequence)
    glycosites = [i for i in range(len(tokens)) if "HexNAc" in tokens[i][1]]
    return glycosites


def total_expected_ions_with_hexnac(row, kind):
    if kind not in ["b", "y"]:
        raise KeyError("Can calculate total_expected_ions_with_hexnac\
         for b or y ions. {kind} requested".format(kind=kind))
    if kind == "b":
        # exclude b1 ion
        return row.peptideLens - min(row.glycosylation_sites) - 1
    else:
        # Distance from the end to last glycosylation site
        return row.peptideLens - (row.peptideLens - max(row.glycosylation_sites))


def percent_expected_ions_with_hexnac_observed(row):
    expected_b_hexnac = total_expected_ions_with_hexnac(row, "b")
    expected_y_hexnac = total_expected_ions_with_hexnac(row, "y")
    b_ions_hexnac_observed = len(row.b_ions_with_HexNAc)
    y_ions_hexnac_observed = len(row.y_ions_with_HexNAc)

    # if(expected_b_hexnac < 2):
    #     print("{seq} does not generate good b ions with hexnac".format(seq=row.Glycopeptide_identifier))
    # if(expected_y_hexnac < 2):
    #     print("{seq} does not generate good y ions with hexnac".format(seq=row.Glycopeptide_identifier))

    percent_expected_observed = (b_ions_hexnac_observed + y_ions_hexnac_observed) \
        / float(expected_b_hexnac + expected_y_hexnac)

    return percent_expected_observed
