from .noresm_qad_diagnostic import NorESMQADD

# Package-level logger setup:
# we just provide a logger and a NullHandler
import logging
from .constants import PACKAGE_LOGGER_NAME
# Use a constant package name so logs are consistently tagged "noresm_qad_diagnostic"
logger = logging.getLogger(PACKAGE_LOGGER_NAME)
# Prevent 'No handler found' warnings if the application didn't configure logging.
# This does NOT prevent messages propagating to the root logger if the app configured logging.
logger.addHandler(logging.NullHandler())

# Export the package logger for convenience
__all__ = ["logger", "PACKAGE_LOGGER_NAME", "NorESMQADD"]
