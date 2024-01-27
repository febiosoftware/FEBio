#! /bin/bash
parent_path=$(dirname $0)
sh --login -i -c "${parent_path}/create-sdk-wrapped.sh"
