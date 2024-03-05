# -*- coding: utf-8 -*-
# Created by Felipe Lopes de Oliveira
# Distributed under the terms of the MIT License.

"""
The logger creation function.
"""

import logging.config
import logging.handlers


def create_logger(level: str = "DEBUG",
                  format: str = "detailed",
                  save_to_file: bool = False,
                  log_filename: str = 'pycofbuilder.log'):
    """
    Build a logger with the given level and format.

    Parameters
    ----------
    level : str
        The logging level. Default is "DEBUG".
        Can be one of "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL".
    format : str
        The logging format. Default is "detailed".
        Can be one of "simple", "detailed".
    save_to_file : bool
        Whether to save the logs to a file. Default is False.
    log_filename : str
        The file to save the logs to. Default is "pycofbuilder.log".

    Returns
    -------
    logger : logging.Logger
        The logger object.
    """

    # Check if the parameters are valid
    allowed_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
    assert level.upper() in allowed_levels, "Invalid level, must be one of {}".format(allowed_levels)

    allowed_formats = ["simple", "detailed"]
    assert format.lower() in allowed_formats, "Invalid format, must be one of {}".format(allowed_formats)

    # Set up logging
    config_log = {
      "version": 1,
      "disable_existing_loggers": False,
      "formatters": {
        "simple": {
          "format": "%(message)s"
        },
        "detailed": {
          "format": "[%(levelname)s|%(module)s|L%(lineno)d] %(asctime)s: %(message)s",
          "datefmt": "%Y-%m-%dT%H:%M:%S%z"
        }
      },
      "handlers": {
        "stderr": {
          "class": "logging.StreamHandler",
          "level": level.upper(),
          "formatter": format.lower(),
          "stream": "ext://sys.stderr"
        },
      },
      "loggers": {
          "root": {"level": level.upper(), "handlers": ["stderr"]}
          }
    }

    if save_to_file:
        config_log["handlers"]["file"] = {
          "class": "logging.FileHandler",
          "level": level.upper(),
          "formatter": format.lower(),
          "filename": log_filename,
          "mode": "a",
          "encoding": "utf-8"
        }
        config_log["loggers"]["root"]["handlers"].append("file")

    logger = logging.getLogger("pycofbuilder")

    logging.config.dictConfig(config_log)

    return logger
