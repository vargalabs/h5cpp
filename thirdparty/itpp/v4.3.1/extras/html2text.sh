#!/bin/sh
w3m -dump -cols 78 $1 \
	| sed '1,/IT++ Compilation and Installation using Microsoft Visual/d' \
	| sed '/How To Set Up a Local/,$d' \
  | recode -f utf8..ascii

