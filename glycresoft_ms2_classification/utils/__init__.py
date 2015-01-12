import json
import os


splitext = os.path.splitext

__all__ = ["memoize", "try_deserialize", "try_get_outfile", "collectiontools"]


def try_deserialize(input_parameter):
    if isinstance(input_parameter, basestring):
        results = json.load(open(input_parameter))
    elif isinstance(input_parameter, file):
        results = json.load(input_parameter)
    else:
        results = input_parameter
    return results


def try_get_outfile(input_parameter, ext):
    if isinstance(input_parameter, basestring):
        return splitext(splitext(input_parameter)[0])[0] + "." + ext
    elif isinstance(input_parameter, file):
        return try_get_outfile(input_parameter.name, ext)
    elif "metadata" in input_parameter:
        origin_name = input_parameter["metadata"]["ms1_output_file"]
        origin_name = splitext(origin_name)[0]
        tag_name = input_parameter["metadata"].get("tag", "")
        if tag_name != "":
            origin_name += '.' + tag_name
        return origin_name + "." + ext
    elif hasattr(input_parameter, "metadata"):
        origin_name = input_parameter.metadata["ms1_output_file"]
        origin_name = splitext(origin_name)[0]
        tag_name = input_parameter.metadata.get("tag", "")
        if tag_name != "":
            origin_name += '-' + tag_name
        return origin_name + "." + ext
    else:
        raise ValueError("Could not interpolate an output file name for {0}".format(input_parameter))
