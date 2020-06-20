"""
Modula -- a modular calculation framework for Python

Especially for scientific calculations and stuff
"""
from __future__ import print_function

import os
import sys
from collections import deque

from gfam.modula.configuration import Configuration
from gfam.modula.log import init_master_logger
from gfam.modula.module import DefaultModuleManager
from gfam.modula.storage import DiskStorageEngine

__version__ = "0.1"

CONFIG = None
LOGGER = None
MODULE_MANAGER = None
STORAGE_ENGINE = None
if sys.version_info[0] >= 3:
    unicode = str


def init(configuration=".",
         module_manager_factory=DefaultModuleManager,
         storage_engine_factory=DiskStorageEngine,
         debug=False):
    """Initializes the Modula engine

    `configuration` is either a `Configuration` instance or the name of the
    directory containing the module configuration file.

    `module_manager_factory` is a factory routine that constructs
    concrete `ModuleManager` instances. It is safe to leave it at its
    default value.

    `debug` is ``True`` if debug messages should be printed; otherwise
    it is ``False``.
    """
    global CONFIG, MODULE_MANAGER, STORAGE_ENGINE, LOGGER

    if isinstance(configuration, (str, unicode)):
        CONFIG = Configuration(rootdir=configuration)
    else:
        CONFIG = Configuration(cfg=configuration)

    LOGGER = init_master_logger(debug=debug)
    MODULE_MANAGER = module_manager_factory(CONFIG)
    STORAGE_ENGINE = storage_engine_factory(CONFIG, MODULE_MANAGER)


def init_project(rootdir):
    """Initializes a Modula project in a directory"""

    if not os.path.isdir(rootdir):
        os.mkdir(rootdir)

    for directory in ["lib", "figures", "modules", "storage"]:
        full_path = os.path.join(rootdir, directory)
        if not os.path.isdir(full_path):
            os.mkdir(full_path)

    modules_fh = os.path.join(rootdir, "modules.cfg")
    if not os.path.isfile(modules_fh):
        modules_fh = open(modules_fh, "w")
        print("""[@global]

[@inputs]
# Enter input files here in the following format:
# id1=path
# id2=path
# ...
""", file=modules_fh)
        modules_fh.close()


def run(module_name, force=False, extra_args=None):
    """Runs the given module in the Modula framework"""
    global CONFIG, MODULE_MANAGER, STORAGE_ENGINE, LOGGER

    # Check dependencies, run them if needed
    to_check, to_run = deque([module_name]), []
    while to_check:
        name = to_check.popleft()

        module = MODULE_MANAGER.get(name)
        last_updated_at = module.get_last_updated_at()
        depends = module.get_dependencies()

        should_run = False
        for dependency in depends:
            dep = MODULE_MANAGER.get(dependency)
            if dep.get_last_updated_at() >= last_updated_at:
                should_run = True
                if dependency not in to_check:
                    to_check.append(dependency)
            else:
                LOGGER.debug("%s is newer than %s, not running",
                             dependency, name)

        if should_run:
            to_run.append(name)

    to_run.reverse()

    if not to_run and force:
        to_run = [module_name]

    if not to_run:
        LOGGER.info("Nothing to do")
        return

    # Run the modules that we collected, store the results
    for name in to_run:
        module = MODULE_MANAGER.get(name)
        if name == module_name and extra_args is not None:
            module.add_extra_args(extra_args)
        result = module.run()
        if result is not None:
            STORAGE_ENGINE.store(module, result)
