"""
Collection of custom exception classes for specific cases
"""


class Error(Exception):
    pass


class FrameshiftException(Error):
    pass


class InternalStopCodonException(Error):
    pass


class InvalidEntryException(Error):
    pass


class Cys104NotFoundException(Error):
    pass
