#!/bin/bash

cd Bembel;
cpplint --quiet --extensions=hpp,cpp --root=.. --recursive ./src
