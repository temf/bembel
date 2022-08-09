#!/bin/bash

cd Bembel;
cpplint --quiet --extensions=hpp,cpp --root=.. --recursive ./src 2>&1 | tee ../scripts/cpplint_output.txt
