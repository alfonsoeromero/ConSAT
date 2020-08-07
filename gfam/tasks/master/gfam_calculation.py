import os

import gfam.modula as modula
from gfam.modula.module import CalculationModule
from gfam.scripts import CommandLineApp
from gfam.utilities.file_utils import redirected


class GFamCalculation(CalculationModule):
    """Class representing a GFam calculation step. This is a subclass of
    `modula.CalcuationModule`_ and it assumes that the name of the module
    refers to a valid Python module in `gfam.scripts`_."""

    def __init__(self, *args, **kwds):
        super(GFamCalculation, self).__init__(*args, **kwds)
        self.modula = None
        self.extra_args = None

    def add_extra_args(self, args):
        """ Adds the possibility of "hotplugging" extra arguments after
            the configuration file has been read
        """
        self.extra_args = args

    def run(self):
        """Runs the calculation"""
        self.logger.info("Starting module %s", self.name)

        self.prepare()

        # Search for the CommandLineApp object in the module
        app = []
        for value in self.module.__dict__.values():
            if isinstance(value, type) and value != CommandLineApp \
                    and issubclass(value, CommandLineApp):
                app.append(value)

        if len(app) != 1:
            raise ValueError("more than one CommandLineApp in %s" % self.name)

        # Create the application
        app = app[0](logger=self.logger)
        args = ["-c", self.config.get("@global.config_file")]
        # add some extra args, if any
        try:
            self.extra_args
        except NameError:
            self.extra_args = []
        if self.extra_args is not None:
            args.extend(self.extra_args)

        for param, value in self.parameters.items():
            if not param.startswith("switch."):
                continue
            switch, value = value.split(" ", 1)
            value = modula.STORAGE_ENGINE.get_filename(value.strip())
            args.extend([switch, value])

        if "infile" in self.parameters:
            infiles = self.parameters["infile"].split(",")
            for infile in infiles:
                infile = modula.STORAGE_ENGINE.get_filename(infile.strip())
                args.append(infile)

        if "stdin" in self.parameters:
            stdin = modula.STORAGE_ENGINE.get_source(self.parameters["stdin"])
        else:
            stdin = None

        out_fname = modula.STORAGE_ENGINE.get_filename(self.name)
        stdout = modula.STORAGE_ENGINE.get_result_stream(self, mode="wb")
        try:
            with redirected(stdin=stdin, stdout=stdout):
                retcode = app.run(args)
            stdout.close()
            if retcode:
                raise RuntimeError("non-zero return code from child module")
        except RuntimeError:
            # If an error happens, remove the output file and re-raise
            # the exception
            stdout.close()
            os.unlink(out_fname)
            raise

        self.logger.info("Finished module %s", self.name)
