'''Centralized Definition of Exception Definitions'''
import os
import json
import traceback
import warnings


# Forward Declaration so that the child class can be
# used in the parent's definition.
ErrorCodeErrorException = Exception


class ErrorCodingMeta(type):

    '''Defines a metaclass that assigns a unique negative integer
    for labeling the type of error that occurred, for easy communication
    with another runtime.
    '''
    errcode = -1
    ERRCODE_CAP = -255
    errmap = {}

    def __new__(cls, name, parents, attrs):
        attrs['errcode'] = cls.errcode
        cls.errmap[cls.errcode] = name
        cls.errcode -= 1

        if cls.errcode < cls.ERRCODE_CAP:
            raise ErrorCodeErrorException(
                "Out of POSIX compliant error codes!")
        if "error_code_interface" not in __name__:
            warnings.warn(
                "Error codes should be defined in the error_code_interface.py file.\
                 This guarantees deterministic code assignment. " + __name__)
        return type.__new__(cls, name, parents, attrs)

    @classmethod
    def build_error_code_map(cls, output_file=None):
        if output_file is None:
            output_file = os.path.join(os.path.dirname(__file__), "exit_code_map.json")
        json.dump(cls.errmap, open(output_file, 'wb'), indent=4)
        return output_file


class GlycReSoftInterprocessCommunicationException(Exception):
    __metaclass__ = ErrorCodingMeta

    def __init__(self, msg):
        super(GlycReSoftInterprocessCommunicationException, self).__init__(msg)
        traceback.print_exc()


class ErrorCodeErrorException(GlycReSoftInterprocessCommunicationException):

    '''An error that is raised when there are no more
    error codes to assign. Figures, huh?'''
    __metaclass__ = ErrorCodingMeta
    pass


class RscriptException(GlycReSoftInterprocessCommunicationException):

    '''An error that is raised when a task communicating
    with R fails.'''
    __metaclass__ = ErrorCodingMeta
    pass


class NoIonsMatchedException(GlycReSoftInterprocessCommunicationException):
    __metaclass__ = ErrorCodingMeta
    pass


class NoSitesFoundWrapperException(GlycReSoftInterprocessCommunicationException):
    __metaclass__ = ErrorCodingMeta
    pass


class UnqualifiedModificationWrapperException(GlycReSoftInterprocessCommunicationException):
    __metaclass__ = ErrorCodingMeta
    pass


class ModelFitException(GlycReSoftInterprocessCommunicationException):
    __metaclass__ = ErrorCodingMeta
    pass


class ClassificationException(GlycReSoftInterprocessCommunicationException):
    __metaclass__ = ErrorCodingMeta
    pass


class ImportErrorWrapperException(GlycReSoftInterprocessCommunicationException):
    __metaclass__ = ErrorCodingMeta
    pass
