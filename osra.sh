#!/bin/sh
export OMP_NUM_THREADS=1
GS_PATH=/usr/local/bin
export PATH=$PATH:$GS_PATH
export MAGICK_CONFIGURE_PATH=/opt/local/osra/1.3.0/
/opt/local/osra/1.3.0/osra.bin $*
