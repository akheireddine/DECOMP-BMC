#! /bin/bash

echo "`grep -e "^[0-9][0-9]*$"  $1 | paste -sd+ | bc` nanoseconds"
