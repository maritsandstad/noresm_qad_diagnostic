# Optional helper functions for applications that want to configure package logging with a recommended format.
import logging
from noresm_qad_diagnostic import PACKAGE_LOGGER_NAME

def get_logger(
    name:   str | None = None
) -> logging.Logger:
    """Return a package-scoped logger.
    Example: get_logger('module') -> 'noresm_qad_diagnostic.module'

    Parameters
    ----------
    name : str | None, optional
        name for logger object.

    Returns
    -------
    logging.Logger : logging.Logger
        Logger object with the name.
    """
    return logging.getLogger(PACKAGE_LOGGER_NAME + (f".{name}" if name else ""))

def configure_default_logging(
    level:  int | None = None,
    fmt:    str | None = None
) -> None:
    """Optionally configure a package-level StreamHandler with a recommended format.
    Usecase is if you want a different logging level for noresm-qad-diagnostic functions.

    Parameters
    ----------
    level : int | None, optional
        If provided, sets the package logger level.
        None will leave it NOTSET to defer to root), by default None.
    fmt : str | None, optional
        Optional custom format string.
        None will include time and a literal 'noresm-qad-diagnostic' tag, by default None.
    """
    logger = logging.getLogger(PACKAGE_LOGGER_NAME)

    if level is not None:
        # setLevel is optional — setting it to NOTSET is the safest default to defer decisions to root logger
        logger.setLevel(level)

    # If a handler of the type StreamHandler is already present, don't add another (avoids duplicate messages)
    if any(isinstance(h, logging.StreamHandler) for h in logger.handlers):
        return

    handler = logging.StreamHandler()
    fmt = fmt if fmt is not None else "%(asctime)s - noresm-qad-diagnostic - %(name)s - %(levelname)s - %(message)s"
    formatter = logging.Formatter(fmt)
    handler.setFormatter(formatter)

    logger.addHandler(handler)