#! /bin/bash
mkdir -p ~/.ssh/
echo "$SSH_KEY" > ~/.ssh/repo.key
chmod 600 ~/.ssh/repo.key
cat >>~/.ssh/config <<END
Host repo
	HostName $SSH_HOST
	User $SSH_USER
	IdentityFile ~/.ssh/repo.key
	StrictHostKeyChecking no
END
