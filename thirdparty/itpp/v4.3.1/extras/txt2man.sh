#!/bin/sh
if [ -f itpp-config.t2t ]; then
	txt2tags -i itpp-config.t2t -o ../itpp-config.1.in
fi
