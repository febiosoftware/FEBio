resource "aws_s3_bucket" "packages" {
  bucket = "febio-packages"

  lifecycle {
    prevent_destroy = false
  }
}

resource "aws_s3_bucket_acl" "packages" {
  bucket = aws_s3_bucket.packages.id
  acl    = "private"
}

resource "aws_ecr_repository" "febio-runtime" {
  name                 = "febio-runtime"
  image_tag_mutability = "MUTABLE"
}

output "repository_url" {
  value = aws_ecr_repository.febio-runtime.repository_url
}
