#!/bin/bash

if [ -r /opt/crc/Modules/current/init/bash ]; then
        source /opt/crc/Modules/current/init/bash
fi

./main.out ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9}

bzip2 participants/${6}
bzip2 recurrent/${7}
bzip2 trial/${8}

