data "aws_iam_policy_document" "assume_role" {
  statement {
    sid    = ""
    effect = "Allow"

    principals {
      identifiers = ["lambda.amazonaws.com"]
      type        = "Service"
    }

    actions = ["sts:AssumeRole"]
  }
}

resource "aws_iam_role" "default" {
  name               = local.project
  assume_role_policy = data.aws_iam_policy_document.assume_role.json
}

data "aws_iam_policy_document" "grant_bucket_access" {
  statement {
    sid = ""
    actions = [
      "s3:Get*",
      "s3:List*",
      "s3:Put*",
    ]

    resources = [
      aws_s3_bucket.lambda.arn,
      "${aws_s3_bucket.lambda.arn}/*"
    ]

  }
}

resource "aws_iam_role_policy" "s3" {
  name   = "${local.project}-s3"
  role   = aws_iam_role.default.id
  policy = data.aws_iam_policy_document.grant_bucket_access.json
}

data "aws_iam_policy_document" "lambda" {
  statement {
    effect = "Allow"
    actions = [
      "logs:CreateLogGroup",
      "logs:CreateLogStream",
      "logs:PutLogEvents",
    ]
    resources = ["*"]
  }

  statement {
    effect    = "Allow"
    actions   = ["lambda:InvokeFunction"]
    resources = [aws_lambda_function.lambda.arn]
  }
}

resource "aws_iam_role_policy" "lambda" {
  name   = "${local.project}-lambda"
  role   = aws_iam_role.default.id
  policy = data.aws_iam_policy_document.lambda.json
}
