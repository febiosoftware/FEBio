data "aws_ami_ids" "github-ci" {
  owners = ["self"]

  filter {
    name   = "name"
    values = ["*intel-oneapi*", "*febio*"]
  }
}

resource "aws_ec2_tag" "github-ci" {
  for_each    = toset(data.aws_ami_ids.github-ci.ids)
  resource_id = each.value
  key         = "for"
  value       = "Github CI"
}
