import boto3


class S3:
    def __init__(self):
        self.s3 = boto3.resource("s3")

    def download_file(self, bucket, key, path):
        object = self.s3.Object(bucket, key)
        object.download_fileobj(path)

    def write_object(self, bucket, key, data):
        object = self.s3.Object(bucket, key)
        object.upload_fileobj(data)
