import re

from collections import defaultdict

try:
    from lxml import etree as ET
except:
    try:
        from xml.etree import cElementTree as ET
    except:
        from xml.etree import ElementTree as ET


def to_string(xml_tree):
    return ET.tostring(xml_tree)


def recursive_dict(element, strip_ns=True):
    data_dict = {k: convert_str(v) for k, v in element.attrib.items()}
    children = map(recursive_dict, element)
    children_nodes = defaultdict(list)
    clean_nodes = {}
    for node, data in children:
        children_nodes[node].append(data)
    for node, data_list in children_nodes.items():
        if len(data_list) == 1:
            clean_nodes[node] = data_list[0]
        else:
            clean_nodes[node] = data_list

    param_to_key_value(clean_nodes)

    # If there are nodes mapped from children, use them
    if clean_nodes:
        data_dict.update(clean_nodes)
    elif element.text is not None and not element.text.isspace():
        data_dict["text"] = element.text

    # Simple entries shouldn't be dictionaries
    if len(data_dict) == 1 and "text" in data_dict:
        data_dict = convert_str(data_dict["text"])
    elif len(data_dict) == 1 and "name" in data_dict:
        data_dict = data_dict["name"]
    tag = element.tag
    if strip_ns:
        tag = re.sub(r"\{.*\}", "", tag)
    return tag, data_dict


def remove_namespace(doc, namespace):
    """Remove namespace in the passed document in place."""
    ns = u'{%s}' % namespace
    nsl = len(ns)
    for elem in doc.getiterator():
        if elem.tag.startswith(ns):
            elem.tag = elem.tag[nsl:]


def convert_str(content_string):
    content_string_lower = content_string.lower()
    if content_string_lower in falsey_strings:
        return False
    elif content_string_lower in truey_strings:
        return True
    else:
        try:
            return int(content_string_lower)
        except ValueError:
            try:
                return float(content_string_lower)
            except ValueError:
                return content_string


falsey_strings = set(("false",))
truey_strings = set(("true",))


def cv_term_update(data_dict, cv_param):
    if 'value' in cv_param:
        value = cv_param["value"]
        try:
            data_dict[cv_param["name"]] = float(value)
        except ValueError:
            data_dict[cv_param["name"]] = value
    else:
        data_dict["name"] = cv_param["name"]


def param_to_key_value(data_dict):
    try:
        cvs = data_dict.pop("cvParam")
        if isinstance(cvs, list):
            for cv_param in cvs:
                cv_term_update(data_dict, cv_param)
        else:
            cv_term_update(data_dict, cvs)
    except:
        pass
    try:
        cvs = data_dict.pop("userParam")
        if isinstance(cvs, list):
            for cv_param in cvs:
                cv_term_update(data_dict, cv_param)
        else:
            cv_term_update(data_dict, cvs)
    except:
        pass
