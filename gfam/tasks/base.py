import logging


class LoggedTask:
    def __init__(self, logger: logging.Logger = None):
        if logger is None:
            self.log = logging.getLogger(__file__)
            self.log.setLevel(logging.DEBUG)
            logging.basicConfig(level=logging.DEBUG, format="%(message)s")
        else:
            self.log = logger
