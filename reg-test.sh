#!/bin/bash

declare -A svg2nc_args
svg2nc_args[parquetry]="-t -.01 -x .15 -c 008000:.25 -c 008080:.2 -d .125 -m .25"
svg2nc_args[sample]="-t -.01 -x .1 -c 5a131f:.25 -c 000080:.2 -d .125 -m .25 -e 000000 -a"

fail_count=0
for name in ${!svg2nc_args[*]}; do
	base=test/${name}
	args=${svg2nc_args[${name}]}
	set -xe
	./svg2nc -v -f 12 -l .2 ${args} -n ${base}.ngc~ -p ${base}.ps~ ${base}.svg
	set +xe
	for type in ngc ps; do
		diff ${base}.${type} ${base}.${type}~
		if [ $? -ne 0 ]; then
			echo "FAILURE: Found differences in ${base}.${type}"
			fail_count=$((fail_count+1))
		else
			echo "PASS: ${base}.${type}"
		fi
	done
done

if [ ${fail_count} -ne 0 ]; then
	echo "TESTS FAILED"
else
	echo "ALL TESTS PASSED"
fi
