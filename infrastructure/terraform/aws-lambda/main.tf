locals {
  name    = "febio-lambda"
  project = "${local.name}-${terraform.workspace}"

  filter_prefix = var.filter_prefix == "" ? null : var.filter_prefix

  filter_suffix = var.filter_suffix == "" ? null : var.filter_suffix
}

data "aws_ecr_repository" "febio-runtime-lambda" {
  name = "febio-runtime-lambda"
}


resource "aws_s3_bucket" "lambda" {
  bucket = local.project

  lifecycle {
    prevent_destroy = false
  }
}

resource "aws_lambda_permission" "s3" {
  statement_id  = "AllowBucketInvocation"
  action        = "lambda:InvokeFunction"
  function_name = aws_lambda_function.lambda.arn
  principal     = "s3.amazonaws.com"
  source_arn    = aws_s3_bucket.lambda.arn
}

resource "aws_lambda_function" "lambda" {
  image_uri     = "${data.aws_ecr_repository.febio-runtime-lambda.repository_url}:latest"
  package_type  = "Image"
  function_name = local.project
  role = aws_iam_role.default.arn
  source_code_hash = trimprefix(data.aws_ecr_repository.febio-runtime-lambda.id, "sha256:")
  timeout = 60
}

resource "aws_s3_bucket_notification" "lambda" {
  bucket = aws_s3_bucket.lambda.id

  lambda_function {
    lambda_function_arn = aws_lambda_function.lambda.arn
    events              = ["s3:ObjectCreated:*"]
    filter_prefix       = local.filter_prefix
    filter_suffix       = local.filter_suffix
  }

  depends_on = [aws_lambda_permission.s3]
}
