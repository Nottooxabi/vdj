class Error(Exception):
    pass


class FrameshiftException(Error):
    pass


class InternalStopCodonException(Error):
    pass

class InvalidEntryException(Error):
    pass