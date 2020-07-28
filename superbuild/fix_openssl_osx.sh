#!/bin/bash
#
# Prerequisites to compile superbuild in a 
# MAC OS X environment.
#
# Arnau Miro, OGS 2020

OPENSSL_VERS=OpenSSL_1_1_1

git clone -b $OPENSSL_VERS https://github.com/openssl/openssl.git openssl
cd openssl

./config
make
sudo make install

cd ..
rm -rf openssl