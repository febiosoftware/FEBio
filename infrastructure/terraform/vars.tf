variable "repo_name" {
  type    = string
  default = "FEBio"
}

variable "aws_access_key_id" {
  type      = string
  sensitive = true
}

variable "aws_secret_access_key" {
  type      = string
  sensitive = true
}

variable "gh_token" {
  type      = string
  sensitive = true
}

variable "repo_host" {
  type      = string
  sensitive = false
}

variable "repo_user" {
  type      = string
  sensitive = false
}

variable "repo_key" {
  type      = string
  sensitive = true
}

variable "repo_action_permissions" {
  type    = list(string)
  default = ["aws/*", "tundranerd/*"]
}

