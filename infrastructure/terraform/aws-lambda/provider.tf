terraform {
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 4.16"
    }

  }

  backend "s3" {
    bucket               = "febio-tf-state"
    key                  = "febio-aws-lambda.tfstate"
    region               = "us-east-1"
    workspace_key_prefix = "febio-aws-lambda"
    profile              = "febio"
    encrypt              = true
  }
}

provider "aws" {
  default_tags {
    tags = {
      provisioner = "Terraform"
      workspace   = terraform.workspace
    }
  }
}
