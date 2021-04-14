#!/bin/bash

if [ -r /opt/crc/Modules/current/init/bash ]; then
        source /opt/crc/Modules/current/init/bash
fi

./main.out ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10}

bzip2 ${6}
bzip2 ${7}
bzip2 ${8}
bzip2 ${9}
