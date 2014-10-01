
cdef int AA_STATE = 0
cdef int MOD_STATE = 0

cdef list sequence_tokenizer(str sequence):
    cdef int state
    state = AA_STATE
    cdef list chunks = list()
    cdef str current_aa = ""
    cdef str current_mod = ""
    cdef int paren_level = 0
    cdef int i = 0
    cdef int length = len(sequence)
    cdef str next_char = None
    while i < length:
        next_char = sequence[i]
        if next_char is "(":
            if state == AA_STATE:
                state = MOD_STATE
                assert paren_level == 0
                paren_level += 1
            elif state == MOD_STATE:
                paren_level += 1
                current_mod += next_char
        elif next_char == ")":
            if state == AA_STATE:
                raise Exception("Invalid Sequence. ) found outside of modification.")
            elif state == MOD_STATE:
                paren_level -= 1
                if paren_level == 0:
                    state = AA_STATE
                    new_chunk = tuple(current_aa, current_mod)
                    chunks.append(new_chunk)
                    current_mod = ""
                    current_aa = ""
                else:
                    current_mod += next_char
        elif state == AA_STATE:
            if(current_aa != ""):
                new_chunk = tuple(current_aa, current_mod)
                chunks.append(new_chunk)
                current_mod = ""
                current_aa = ""
            current_aa += next_char
        elif state == MOD_STATE:
            current_mod += next_char
        else:
            raise Exception("Unknown Tokenizer State", current_aa, current_mod, i, next_char)
        i += 1
    if current_aa != "":
        new_chunk = tuple(current_aa, "")
        chunks.append(new_chunk)
    return chunks
