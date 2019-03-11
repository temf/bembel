#!/bin/bash
wget http://bitbucket.org/eigen/eigen/get/3.3.7.zip
rm -r src/include/eigen3
unzip 3.3.7.zip -d eigen3/
mv eigen3 src/include
rm 3.3.7.zip*
cd src/include/eigen3/eigen*
mv * ..
