packer {
  required_plugins {
    amazon = {
      version = ">= 0.0.2"
      source  = "github.com/hashicorp/amazon"
    }
  }
}

locals {
  buildtime          = formatdate("YYYYMMDDhhmm", timestamp())
  installation_path  = var.installation_path
  src_path           = "${local.installation_path}\\src"
  vcpkg_package_path = "${local.installation_path}\\febio"
  environment = {
    "INSTALLATION_PATH"  = local.installation_path
    "SOURCE_PATH"        = local.src_path
    "VCPKG_PACKAGE_PATH" = local.vcpkg_package_path
  }
}

variable "installation_path" {
  default = "c:\\usr\\local"
}

data "amazon-ami" "windows" {
  filters = {
    name             = "packer-provisioned-windows-2019-intel-oneapi-*"
    root-device-type = "ebs"
  }

  most_recent = true
  owners      = ["353328028284"]
  region      = "us-east-1"
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
  instance_type = "c5a.8xlarge"
  source_ami    = data.amazon-ami.windows.id

  communicator = "winrm"

  skip_create_ami = var.skip_create_ami

  iam_instance_profile = "s3-read-access"

  aws_polling {
    delay_seconds = 60
    max_attempts  = 90
  }


  user_data_file = "./common/windows/user-data.txt"
  winrm_password = data.amazon-parameterstore.winrm_password.value
  winrm_username = data.amazon-parameterstore.winrm_username.value

  launch_block_device_mappings {
    device_name           = "/dev/sda1"
    volume_size           = 150
    volume_type           = "gp2"
    delete_on_termination = true
  }
}

build {
  name = "febiostudio"
  sources = [
    "source.amazon-ebs.windows"
  ]

  # Conditionally maximize partition
  provisioner "powershell" {
    script = "./common/windows/resize-partition.ps1"
    env    = local.environment
  }

  # paths
  provisioner "powershell" {
    script = "./common/windows/pathprep.ps1"
    env    = local.environment
  }

  provisioner "powershell" {
    script = "./common/windows/jq.ps1"
    env    = local.environment
  }

  provisioner "powershell" {
    script = "./common/windows/msmpi.ps1"
    env    = local.environment
  }

  provisioner "powershell" {
    script = "./common/windows/aws.ps1"
    env    = local.environment
  }

  provisioner "powershell" {
    script = "./common/windows/install-builder.ps1"
    env    = local.environment
  }

  # Lua 5.3
  provisioner "powershell" {
    script = "./common/windows/lua.ps1"
    env    = local.environment
  }

  provisioner "powershell" {
    script = "./common/windows/ffmpeg.ps1"
    env    = local.environment
  }

  # vcpkg
  provisioner "powershell" {
    script = "./common/windows/vcpkg-installer.ps1"
    env    = local.environment
  }

  #vcpkg packages
  provisioner "file" {
    source      = "./common/windows/vcpkg.json"
    destination = "${local.vcpkg_package_path}\\"
  }

  provisioner "powershell" {
    script = "./common/windows/vcpkg-package-install.ps1"
    env    = local.environment
  }

  # LEVMAR
  provisioner "windows-shell" {
    script = "./common/windows/levmar.bat"
    env    = local.environment
  }

  # HYPRE
  provisioner "windows-shell" {
    script = "./common/windows/hypre.bat"
    env    = local.environment
  }

  # # mmg
  provisioner "windows-shell" {
    script = "./common/windows/mmg.bat"
    env    = local.environment
  }

  # # tetgen
  provisioner "windows-shell" {
    script = "./common/windows/tetgen.bat"
    env    = local.environment
  }

  # itk
  provisioner "windows-shell" {
    script = "./common/windows/itk.bat"
    env    = local.environment
  }

  # sitk
  provisioner "windows-shell" {
    script = "./common/windows/sitk.bat"
    env    = local.environment
  }

  # occt
  provisioner "windows-shell" {
    script = "./common/windows/occt.bat"
    env    = local.environment
  }

  # netgen
  provisioner "windows-shell" {
    script = "./common/windows/netgen.bat"
    env    = local.environment
  }

  # libzip
  provisioner "windows-shell" {
    script = "./common/windows/libzip.bat"
    env    = local.environment
  }

  # sysprep for next launch
  provisioner "powershell" {
    inline = [
      "C:\\ProgramData\\Amazon\\EC2-Windows\\Launch\\Scripts\\InitializeInstance.ps1 -Schedule",
    ]
  }
}
