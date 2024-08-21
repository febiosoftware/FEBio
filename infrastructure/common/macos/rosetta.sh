#!/bin/bash
set -e

if arch -x86_64 /usr/bin/true 2> /dev/null; then
	echo "Rosetta previously installed"
else
	/usr/sbin/softwareupdate --install-rosetta --agree-to-license
fi
