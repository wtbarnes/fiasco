"""
Logging utilities for fiasco.
"""
import logging

from astropy.logger import AstropyLogger

__all__ = ['_init_log', '_config_to_loggerConf']

def _init_log(config=None):
    """
    Initializes the fiasco log.

    In most circumstances this is called automatically when importing
    fiasco. This code is based on that provided by Astropy see
    "licenses/ASTROPY.rst".
    """
    orig_logger_cls = logging.getLoggerClass()
    logging.setLoggerClass(AstropyLogger)
    try:
        log = logging.getLogger('fiasco')
        if config is not None:
            _config_to_loggerConf(config)
        log._set_defaults()
    finally:
        logging.setLoggerClass(orig_logger_cls)

    return log


def _config_to_loggerConf(config):
    """
    Translates a user-provided config to ~`astropy.logger.LoggerConf`.
    """

    if config.has_section('logger'):
        from astropy.logger import Conf as LoggerConf
        conf = LoggerConf()
        loggerconf_option_list = ['log_level', 'use_color', 'log_warnings', 'log_exceptions', 'log_to_file',
                                  'log_file_path', 'log_file_level', 'log_file_format']
        for this_option in loggerconf_option_list:
            if config.has_option('logger', this_option):
                setattr(conf, this_option, config.get('logger', this_option))
    return conf
