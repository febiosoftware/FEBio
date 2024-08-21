#! /bin/bash

download_source() {
	local source=$1

	curl -L -O "$source"
}

extract_source() {
	local archive=$1
	unzip -o "$archive"
}

