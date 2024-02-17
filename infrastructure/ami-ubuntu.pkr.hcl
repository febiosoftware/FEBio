packer {
  required_plugins {
    amazon = {
      version = ">= 0.0.2"
      source  = "github.com/hashicorp/amazon"
    }
  }
}

locals {
  buildtime = formatdate("YYYYMMDDhhmm", timestamp())
}

data "amazon-ami" "ubuntu" {
  filters = {
    name             = "ubuntu/images/hvm-ssd/ubuntu-jammy-22.04-amd64-server*"
    root-device-type = "ebs"
  }

  most_recent = true
  owners      = ["099720109477"]
  region      = "us-east-1"
}

variable "skip_create_ami" {
  type    = bool
  default = false
}

source "amazon-ebs" "ubuntu" {
  ami_name      = "packer-provisioned-ubuntu-22.04-intel-oneapi-${local.buildtime}"
  instance_type = "c5a.4xlarge"
  source_ami    = data.amazon-ami.ubuntu.id
  ssh_username  = "ubuntu"

  skip_create_ami = var.skip_create_ami

  iam_instance_profile = "s3-read-access"

  aws_polling {
    delay_seconds = 60
    max_attempts  = 90
  }

  launch_block_device_mappings {
    device_name           = "/dev/sda1"
    volume_size           = 50
    volume_type           = "gp2"
    delete_on_termination = true
  }
}
variable "image_build_path" {
  type    = string
  default = "/tmp"
}

variable "dependencies_folder" {
  type    = string
  default = "/tmp/linux/dependencies/"
}

build {
  name = "febio"
  sources = [
    "source.amazon-ebs.ubuntu"
  ]

  provisioner "file" {
    source      = "./common/linux"
    destination = var.image_build_path
  }

  provisioner "shell" {
    remote_folder = "${var.image_build_path}/linux"
    script        = "./common/linux/apt.sh"
  }

  # awscli
  provisioner "shell" {
    script = "./common/linux/aws.sh"
  }

  provisioner "shell" {
    remote_folder = "${var.image_build_path}/linux"
    script        = "./common/linux/install-builder.sh"
  }

  provisioner "shell" {
    script = "./common/linux/qt.sh"
  }

  provisioner "shell" {
    script = "./common/linux/openapi.sh"
  }

  # Latest version of cmake (v3.23.2)
  provisioner "shell" {
    script = "./common/linux/cmake.sh"
  }

  # Latest version of git (v2.38.1)
  provisioner "shell" {
    script = "./common/linux/git.sh"
  }

  # awscli
  provisioner "shell" {
    script = "./common/linux/aws.sh"
  }

  provisioner "shell" {
    environment_vars = [
      "DEPENDENCIES_PATH=${var.dependencies_folder}"
    ]
    script = "./common/linux/dependencies/install.sh"
  }
}
