cat All2880_64-QAM.txt| parallel --pipe -N1 --results 'best-64/{#}' --line-buffer --progress ../src/best-parity 4 12 16 64-QAM.txt -
