import os
import sys
import textwrap
from abc import abstractmethod
from configparser import ConfigParser
from io import StringIO

from gfam.modula.hash import sha1
from gfam.scripts import CommandLineApp


class BaseMasterScript(CommandLineApp):
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.modula = None
        self.config = None

    def create_parser(self):
        """Creates the command line parser for the ConSAT master script"""
        parser = super().create_parser()
        parser.add_option("-f", "--force", dest="force", action="store_true",
                          help="force recalculation of results even "
                               "when consat thinks everything is up-to-date")
        parser.add_option("-s", "--silent", dest="silent", action="store_true",
                          help="does not print any output to the terminal")
        return parser

    @abstractmethod
    def get_modula_config_as_string(self) -> str:
        """Return modula config as a string"""
        pass

    @abstractmethod
    def get_default_config_filename(self) -> str:
        """Return default config filename for the script"""
        pass

    def get_modula_config(self, config):
        """Based on the given `ConfigParser` instance in `config`, constructs
        another `ConfigParser` that tells Modula which tasks to execute and
        where each of the input files are to be found."""
        modula_config_str = textwrap.dedent(
            self._get_modula_config_as_string())

        modula_config = ConfigParser()
        modula_config.readfp(StringIO(modula_config_str))

        # Store the name of the config file
        modula_config.set("@global", "config_file", self.options.config_file)

        # Store the hash of the configuration as a default parameter for
        # all the algorithms
        config_str = StringIO()
        config.write(config_str)
        modula_config.set("DEFAULT", "config_file_hash",
                          sha1(config_str.getvalue()).hexdigest())

        # Set up the module and storage path
        modula_config.set("@paths", "modules",
                          os.path.dirname(sys.modules[__name__].__file__))
        modula_config.set("@paths", "storage",
                          config.get("DEFAULT", "folder.work"))

        # Add the input files
        for name, value in config.items("generated"):
            if name.startswith("file."):
                modula_config.set("@inputs", name, value)

        return modula_config

    def read_config(self):
        """Reads the configuration from the given file and returns an
        appropriate `ConfigParser` instance."""
        self.options.config_file = self.options.config_file or "consat.cfg"

        config_file = self.options.config_file
        if not os.path.exists(config_file):
            return None

        config = ConfigParser()
        config.read([config_file])
        return config
