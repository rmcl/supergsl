#!/bin/bash

SRC_DIR=./
DST_DIR=./stubs

#protoc -I=$SRC_DIR --python_out=$DST_DIR $SRC_DIR/sgsl.proto
mkdir -p $DST_DIR
touch $DST_DIR/__init__.py
python -m grpc_tools.protoc \
    -I$SRC_DIR \
    --python_out=$DST_DIR \
    --mypy_out=$DST_DIR \
    --grpc_python_out=$DST_DIR \
    sgsl.proto
