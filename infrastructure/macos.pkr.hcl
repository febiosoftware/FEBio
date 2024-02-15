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
  sudo              = "echo 'packer' | sudo -S /bin/zsh -l -c '{{ .Vars }} {{ .Path }}'"
  zsh_command       = "echo 'packer' | /bin/zsh -c -l '{{ .Vars }} {{ .Path }}'"
  installation_path = var.installation_path
  homebrew_prefix   = "${local.installation_path}/homebrew"
  homebrew_bin      = "${local.homebrew_prefix}/bin/brew"
  src_path          = "${local.installation_path}/src"
  ssh_user          = var.ssh_user
  environment = {
    "HOMEBREW_PREFIX"   = local.homebrew_prefix
    "HOMEBREW_BIN"      = local.homebrew_bin
    "SSH_USER"          = local.ssh_user
    "INSTALLATION_PATH" = local.installation_path
    "SOURCE_PATH"       = local.src_path
  }
}

source "null" "macos" {
  ssh_pty              = true
  ssh_host             = var.ssh_host
  ssh_username         = var.ssh_user
  ssh_private_key_file = var.ssh_private_key
  ssh_port             = var.ssh_port
  communicator         = "ssh"
}

variable "installation_path" {
  default = "/usr/local/x86_64"
}

variable "ssh_host" {
  type    = string
  default = "midnight.local"
}

variable "ssh_user" {
  type = string
}

variable "ssh_port" {
  type    = string
  default = "22"
}

variable "ssh_private_key" {
  type = string
}

build {
  name = "febio"
  sources = [
    "source.null.macos"
  ]

  provisioner "shell" {
    inline = [
      "sudo xcode-select --install || echo \"xcode tools  are already installed\""
    ]
  }

  provisioner "shell" {
    execute_command = local.sudo
    script          = "./common/macos/rosetta.sh"
  }

  provisioner "shell" {
    execute_command = local.sudo
    script          = "./common/macos/openapi.sh"
  }

  provisioner "shell" {
    # execute_command = local.zsh_command
    script = "./common/macos/installation_prep.sh"
    env    = local.environment
  }

  provisioner "shell" {
    script = "./common/macos/homebrew-x86.sh"
    env    = local.environment
  }

  provisioner "shell" {
    script = "./common/macos/homebrew-packages.sh"
    env    = local.environment
  }

  provisioner "shell" {
    execute_command = local.zsh_command
    script = "./common/macos/ffmpeg.sh"
    env    = local.environment
  }

  provisioner "shell" {
    execute_command = local.zsh_command
    script          = "./common/macos/qt.sh"
    env             = local.environment
  }

  provisioner "shell" {
    script = "./common/macos/lua.sh"
    env    = local.environment
  }

  # LEVMAR
  provisioner "shell" {
    execute_command = local.zsh_command
    script          = "./common/macos/levmar.sh"
    env             = local.environment
  }

  # HYPRE
  provisioner "shell" {
    execute_command = local.zsh_command
    script          = "./common/macos/hypre.sh"
    env             = local.environment
  }

  # mmg
  provisioner "shell" {
    execute_command = local.zsh_command
    script          = "./common/macos/mmg.sh"
    env             = local.environment
  }

  # tetgen
  provisioner "shell" {
    execute_command = local.zsh_command
    script          = "./common/macos/tetgen.sh"
    env             = local.environment
  }

  # itk
  provisioner "shell" {
    execute_command = local.zsh_command
    script          = "./common/macos/itk.sh"
    env             = local.environment
  }

  # sitk
  provisioner "shell" {
    execute_command = local.zsh_command
    script          = "./common/macos/sitk.sh"
    env             = local.environment
  }

  # occt
  provisioner "shell" {
    execute_command = local.zsh_command
    script          = "./common/macos/occt.sh"
    env             = local.environment
  }

  # netgen
  provisioner "shell" {
    execute_command = local.zsh_command
    script          = "./common/macos/netgen.sh"
    env             = local.environment
  }

  # libzip
  provisioner "shell" {
    execute_command = local.zsh_command
    script          = "./common/macos/libzip.sh"
    env             = local.environment
  }
}
