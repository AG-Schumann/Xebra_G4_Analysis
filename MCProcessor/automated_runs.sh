#!/bin/bash

make clean; cmake .; make

./XebraMCProcessor -d $1 -i $2



