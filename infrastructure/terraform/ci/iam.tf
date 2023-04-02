resource "aws_iam_user" "gh" {
  name = "tf-gh-actions-febio"
}

resource "aws_iam_access_key" "gh" {
  user = aws_iam_user.gh.name
}

resource "aws_iam_user_group_membership" "gh" {
  user = aws_iam_user.gh.name

  groups = [
    aws_iam_group.gh.name,
  ]
}

resource "aws_iam_group" "gh" {
  name = "github-ci-group"
  path = "/users/"
}

resource "aws_iam_group_policy" "gh" {
  name   = "github-ci-policy"
  group  = aws_iam_group.gh.name
  policy = data.aws_iam_policy_document.gh.json
}


data "aws_iam_policy_document" "gh" {
  statement {
    actions = [
      "ecr:BatchCheckLayerAvailability",
      "ecr:GetDownloadUrlForLayer",
      "ecr:GetRepositoryPolicy",
      "ecr:DescribeRepositories",
      "ecr:ListImages",
      "ecr:DescribeImages",
      "ecr:BatchGetImage",
      "ecr:GetLifecyclePolicy",
      "ecr:GetLifecyclePolicyPreview",
      "ecr:ListTagsForResource",
      "ecr:DescribeImageScanFindings",
      "ecr:InitiateLayerUpload",
      "ecr:UploadLayerPart",
      "ecr:CompleteLayerUpload",
      "ecr:PutImage"
    ]

    resources = [
      "${aws_ecr_repository.febio-runtime.arn}"
    ]

  }

  statement {
    actions = [
      "ecr:GetAuthorizationToken",
    ]
    resources = ["*"]
  }

  statement {
    actions = [
      "ec2:DescribeInstances",
      "ec2:TerminateInstances",
      "ec2:RunInstances",
      "ec2:DescribeInstanceStatus"
    ]

    resources = ["*"]
  }

  statement {
    actions = [
      "ec2:CreateTags"
    ]

    condition {
      test     = "StringEquals"
      variable = "ec2:CreateAction"
      values   = ["RunInstances"]
    }

    resources = ["*"]
  }

  statement {
    actions = [
      "s3:ListBucket",
      "s3:PutObject",
      "s3:GetObject",
      "s3:PutObjectAcl",
      "s3:GetObjectAcl"
    ]
    resources = ["${aws_s3_bucket.packages.arn}"]
  }

}
