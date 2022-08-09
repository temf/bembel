#!/bin/bash

cd Bembel;
cpplint --quiet --extensions=hpp,cpp --root=.. --recursive ./src 2>> ../scripts/cpplint_output.txt;
cat ../scripts/cpplint_output.txt
