'''Centralized Definition of Exception Definitions'''
import os
import json
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
        if __name__ != "error_code_interface":
            warnings.warn(
                "Error codes should be defined in the error_code_interface.py file.\
                 This guarantees deterministic code assignment.")
        return type.__new__(cls, name, parents, attrs)

    @classmethod
    def build_error_code_map(cls, path=os.path.join(os.path.dirname(__file__), "exit_code_map.json")):
        json.dump(cls.errmap, open(path, 'wb'), indent=4)
        return path


class GlycReSoftInterprocessCommunicationException(Exception):
    __metaclass__ = ErrorCodingMeta
    pass


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
