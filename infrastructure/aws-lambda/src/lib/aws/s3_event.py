import urllib.parse


class S3Event:
    def __init__(self, event):
        self.event = event
        self.__parse_record()

    def __parse_record(self):
        self.bucket = self.event["Records"][0]["s3"]["bucket"]["name"]
        self.key = urllib.parse.unquote_plus(
            self.event["Records"][0]["s3"]["object"]["key"],
            encoding="utf8"
        )
