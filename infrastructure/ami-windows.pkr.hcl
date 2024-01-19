packer {
  required_plugins {
    amazon = {
      version = ">= 0.0.2"
      source  = "github.com/hashicorp/amazon"
    }
  }
}

locals {
  buildtime              = formatdate("YYYYMMDDhhmm", timestamp())
  intel_basekit          = "w_BaseKit_p_2022.2.0.252_offline.exe"
  intel_basekit_uri      = "https://registrationcenter-download.intel.com/akdlm/IRC_NAS/18674/${local.intel_basekit}"
  intel_install_windows  = "https://raw.githubusercontent.com/oneapi-src/oneapi-ci/master/scripts/install_windows.bat"
  vs_2019_buildtools_bin = "vs_buildtools.exe"
  vs_2019_buildtools_uri = "https://aka.ms/vs/16/release/${local.vs_2019_buildtools_bin}"
  installation_path      = var.installation_path
}

variable "installation_path" {
  default = "c:\\usr\local"
}

data "amazon-parameterstore" "winrm_password" {
  name            = "/packer/winrm_password"
  with_decryption = true
}

data "amazon-parameterstore" "winrm_username" {
  name            = "/packer/winrm_username"
  with_decryption = true
}


data "amazon-ami" "windows" {
  filters = {
    name             = "Windows_Server-2019-English-Full-Base-*"
    root-device-type = "ebs"
  }

  most_recent = true
  owners      = ["801119661308"]
  region      = "us-east-1"
}

variable "skip_create_ami" {
  type    = bool
  default = false
}

source "amazon-ebs" "windows" {
  ami_name      = "packer-provisioned-windows-2019-intel-oneapi-${local.buildtime}"
  instance_type = "c4.2xlarge"
  source_ami    = data.amazon-ami.windows.id
  communicator  = "winrm"

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
  name = "febio"
  sources = [
    "source.amazon-ebs.windows"
  ]

  # Install Choco
  provisioner "powershell" {
    script = "./common/windows/choco.ps1"
  }

  # VS Build tools
  provisioner "windows-shell" {
    inline = [
      "curl -L -O ${local.vs_2019_buildtools_uri}",
      "start /wait ${local.vs_2019_buildtools_bin} --add Microsoft.VisualStudio.Workload.VCTools --includeOptional --includeRecommended --quiet --nocache --wait",
      "del ${local.vs_2019_buildtools_bin}",
    ]
  }

  # Intel basekit
  provisioner "windows-shell" {
    inline = [
      "curl -L -O ${local.intel_install_windows}",
      "install_windows.bat ${local.intel_basekit_uri}"
    ]
  }

  # paths
  #provisioner "powershell" {
  #  inline = [<<EOF
$u#serpath=[Environment]::GetEnvironmentVariable("Path", "User")
se#tx PATH "$userpath;${local.installation_path}"
EO#F
  #  ]
  #}

  #vcpkg
  #provisioner "powershell" {
  #  script = "./common/windows/vcpkg.ps1"
  #}

  # ZLIB
  # provisioner "windows-shell" {
  #   script = "./common/windows/zlib.bat"
  # }

  # # HYPRE
  # provisioner "windows-shell" {
  #   script = "./common/windows/hypre.bat"
  # }

  # # levmar
  # provisioner "windows-shell" {
  #   script = "./common/windows/levmar.bat"
  # }

  # # levmar header path
  # provisioner "powershell" {
  #   inline = [<<EOF
#$userpath=[Environment]::GetEnvironmentVariable("Path", "User")
#setx PATH "$userpath;${local.installation_path}\include\levmar"
#EOF
  #   ]
  # }

  # # mmg
  # provisioner "windows-shell" {
  #   script = "./common/windows/mmg.bat"
  # }

  # # sysprep for next launch
  provisioner "powershell" {
    inline = [
      "C:\\ProgramData\\Amazon\\EC2-Windows\\Launch\\Scripts\\InitializeInstance.ps1 -Schedule",
    ]
  }
}
