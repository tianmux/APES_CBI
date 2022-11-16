#!/bin/bash
find . -name '*.h' | xargs wc -l > lines_h.txt
find . -name '*.cpp' | xargs wc -l > lines_cpp.txt
find . -name '*.cu' | xargs wc -l > lines_cu.txt


