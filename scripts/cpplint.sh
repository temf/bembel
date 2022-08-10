#!/bin/bash

cpplint --quiet --recursive --exclude=Bembel/src/util/Doxygen.hpp Bembel tests examples
