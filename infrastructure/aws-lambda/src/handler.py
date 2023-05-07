import os
import subprocess  # noqa: F401
from pathlib import Path
from typing import Final

from lib.aws.s3_event import S3Event
from lib.aws.s3 import S3


class Targets:
    TARGET: Final = "/tmp"

    def __init__(self, path):
        self.path = Path(path)
        # Get the input file name sans its extension
        self.stem = self.path.stem
        self.target_dir = os.path.join(Targets.TARGET, self.path.parent)
        self.target_file = os.path.join(self.target_dir, self.path.name)
        self.target_log = os.path.join(self.target_dir, F"{self.stem}.log")
        self.target_plot = os.path.join(self.target_dir, F"{self.stem}.xplt")


class Handler:
    def __init__(self, logger, s3=S3()):
        self.logger = logger
        self.s3 = s3

    def run(self, event, context):
        self.logger.info("Handling event: %s" % event)
        self.logger.info("With context: %s" % context)
        self.logger.info("Parsing s3 event")
        s3_event = S3Event(event)

        self.logger.info("Parsing Targets")
        targets = Targets(s3_event.key)

        self.logger.info("Downloading febio input: s3://%s:%s" % s3_event.bucket, s3_event.key)  # noqa: E501
        self.__download_input(s3_event, targets)

        self.logger.info("Building command")
        cmd = self.__build_command(targets)

        self.__run(cmd)

    def __download_input(self, s3_event, targets):
        if not os.path.exists(targets.target_dir):
            os.makedirs(targets.target_dir)

        with open(targets.target_file) as target:
            self.s3.download_file(s3_event.bucket, s3_event.key, target)

    def __build_command(self, targets):
        cmd = [
            "febio4",
            "-i",
            targets.target_file,
            "-silent",
            "-o",
            targets.target_log,
            "-p",
            targets.target_plot,
        ]

        self.logger.info("Command built: %s" % cmd)
        return cmd

    def __run(self, cmd):
        self.logger.info("Building FEBio plot")
        # result = subprocess.run(cmd, capture_output=True, text=True)
        # self.logger.debug(result.stdout)
