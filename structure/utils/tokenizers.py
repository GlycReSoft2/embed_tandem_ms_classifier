from memoize import memoize_partial_sequence


@memoize_partial_sequence(15000, [1, 4])
def sequence_tokenizer(sequence):
    state = "aa"
    chunks = []
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
                    new_chunk = [current_aa, current_mod]
                    chunks.append(new_chunk)
                    current_mod = ""
                    current_aa = ""
                else:
                    current_mod += next_char
        elif state == "aa":
            if(current_aa != ""):
                new_chunk = [current_aa, current_mod]
                chunks.append(new_chunk)
                current_mod = ""
                current_aa = ""
            current_aa += next_char
        elif state == "mod":
            current_mod += next_char
        else:
            raise Exception("Unknown Tokenizer State", current_aa, current_mod, i, next_char)
        i += 1
    if current_aa != "":
        new_chunk = [current_aa, ""]
        chunks.append(new_chunk)
    return chunks
