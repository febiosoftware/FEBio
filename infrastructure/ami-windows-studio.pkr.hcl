packer {
  required_plugins {
    amazon = {
      version = ">= 0.0.2"
      source  = "github.com/hashicorp/amazon"
    }
  }
}

locals {
  buildtime         = formatdate("YYYYMMDDhhmm", timestamp())
  installation_path = var.installation_path
  python            = "python-3.11.1-amd64.exe"
  python_uri        = "https://www.python.org/ftp/python/3.11.1/python-3.11.1-amd64.exe"
  security_group_id = "sg-0ff74da6968687b23"
}

variable "installation_path" {
  default = "c:\\local"
}

data "amazon-parameterstore" "winrm_password" {
  name            = "/packer/winrm_password"
  with_decryption = true
}

data "amazon-parameterstore" "winrm_username" {
  name            = "/packer/winrm_username"
  with_decryption = true
}

variable "skip_create_ami" {
  type    = bool
  default = true
}

source "amazon-ebs" "windows" {
  ami_name      = "packer-provisioned-windows-2019-febio-studio-${local.buildtime}"
  instance_type = "c4.2xlarge"
  source_ami_filter {
    filters = {
      name             = "packer-provisioned-windows-2019-intel-oneapi-*"
      root-device-type = "ebs"
    }

    most_recent = true
    owners      = ["353328028284"]
  }

  security_group_ids = [local.security_group_id]

  communicator = "winrm"

  skip_create_ami = var.skip_create_ami

  aws_polling {
    delay_seconds = 60
    max_attempts  = 90
  }


  user_data_file = "./common/windows/user-data.txt"
  winrm_password = data.amazon-parameterstore.winrm_password.value
  winrm_username = data.amazon-parameterstore.winrm_username.value

  launch_block_device_mappings {
    device_name           = "/dev/sda1"
    volume_size           = 100
    volume_type           = "gp2"
    delete_on_termination = true
  }
}

build {
  name = "febio-studio"
  sources = [
    "source.amazon-ebs.windows"
  ]


  # qt
  provisioner "windows-shell" {
    script = "./common/windows/qt.bat"
  }


  provisioner "windows-shell" {

    inline = [<<EOF
timeout /T 600 /NOBREAK >NUL
EOF
    ]
  }


  # sysprep for next launch
  #provisioner "powershell" {
  #  inline = [
  #    "C:\\ProgramData\\Amazon\\EC2-Windows\\Launch\\Scripts\\InitializeInstance.ps1 -Schedule",
  #  ]
  #}
}
