import os
import sys
import logging
from typing import Final


class LoggerConfig:
    DEFAULT_LOG_LEVEL: Final = "INFO"
    LOG_LEVELS: Final = ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"]

    def __init__(self, config={}, env=os.environ):
        log_level = env.get("LOG_LEVEL", LoggerConfig.DEFAULT_LOG_LEVEL)
        self.log_level = config.get("log_level", log_level)
        self.__assert_values()
        logger = logging.getLogger()
        handler = logging.StreamHandler(sys.stdout)
        logger.setLevel(self.log_level)
        logger.addHandler(handler)
        self.logger = logger

    def __assert_values(self):
        if self.log_level not in LoggerConfig.LOG_LEVELS:
            self.log_level = LoggerConfig.DEFAULT_LOG_LEVEL
