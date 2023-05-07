from logger import LoggerConfig
from handler import Handler

logger_config = LoggerConfig()
logger = logger_config.logger


def handle(event, context):
    try:
        handler = Handler(logger)
        handler.run(event, context)
    except Exception as e:
        logger.error(e)
