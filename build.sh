#!/bin/bash
#----------

BUILD_DIR=build
USE_OPENMP="ON"

# Parse the command line argument if provided
for arg in "$@"; do
  case $arg in
    USE_OPENMP=*)
      USE_OPENMP="${arg#*=}"
      shift # Remove argument name from processing
      ;;
    *)
      echo "Unknown argument: $arg"
      exit 1
      ;;
  esac
done

echo "USE_OPENMP=${USE_OPENMP}"

mkdir -p ${BUILD_DIR}
( cd ${BUILD_DIR}                          \
  && cmake "-DUSE_OPENMP=${USE_OPENMP}" .. \
  && make                                  \
  && make install                          )
