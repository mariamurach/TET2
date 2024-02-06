#!/bin/bash
awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( I in a ) print I a[i]}' counts/*counts.txt 