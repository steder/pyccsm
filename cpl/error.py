"""
error.py

This module defines a series of exceptions to make it easier to deal with
error handling in the coupler.
"""
class CPLException(Exception):
    def __init__(self, value="Nobody should throw me!"):
        Exception.__init__(self)
        self.value = value
    def __str__(self):
        return repr(self.value)

class AttributeVectorError(CPLException):   
    def __init__(self, value="No error message provided"):
        Exception.__init__(self)
        self.value = value  
    def __str__(self):
        return repr(self.value)

class BundleError(CPLException):   
    def __init__(self, value="No error message provided"):
        Exception.__init__(self)
        self.value = value  
    def __str__(self):
        return repr(self.value)

class CommError(CPLException):   
    def __init__(self, value="No error message provided"):
        Exception.__init__(self)
        self.value = value  
    def __str__(self):
        return repr(self.value)

class ContractError(CPLException):   
    def __init__(self, value="No error message provided"):
        Exception.__init__(self)
        self.value = value  
    def __str__(self):
        return repr(self.value)
    
class ControlError(CPLException):
    def __init__(self, value="No error message provided"):
        Exception.__init__(self)
        self.value = value  
    def __str__(self):
        return repr(self.value)

class DomainError(CPLException):   
    def __init__(self, value="No error message provided"):
        Exception.__init__(self)
        self.value = value  
    def __str__(self):
        return repr(self.value)

class FieldsError(CPLException):   
    def __init__(self, value="No error message provided"):
        Exception.__init__(self)
        self.value = value  
    def __str__(self):
        return repr(self.value)

class InfobufferError(CPLException):   
    def __init__(self, value="No error message provided"):
        Exception.__init__(self)
        self.value = value  
    def __str__(self):
        return repr(self.value)

