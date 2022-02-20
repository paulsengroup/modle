#!/usr/bin/env sh

shasum -a512 *.xz | tee checksums.sha512
