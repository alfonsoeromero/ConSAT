from functools import wraps


def needs_config(func):
    """Decorator for methods in `ConSATMasterScript` that require a valid
    configuration file.
    """
    @wraps(func)
    def wrapper(self, *args, **kwds):
        if self.config is None:
            self.error("Configuration file not found: %s. Try gfam init "
                       "to generate a new one." % self.options.config_file)
        return func(self, *args, **kwds)
    return wrapper
