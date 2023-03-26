terraform {
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 4.16"
    }

    github = {
      source  = "integrations/github"
      version = "~> 5.0"
    }
  }

  backend "s3" {
    bucket               = "febio-tf-state"
    key                  = "febio.tfstate"
    region               = "us-east-1"
    workspace_key_prefix = "febio"
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
