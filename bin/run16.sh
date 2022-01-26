cat All24_16-QAM.txt| parallel --pipe -N1 --results 'best-16/{#}' --line-buffer --progress ../src/crible 4 4 16 16-QAM.txt - H4-16.txt
#cat All24_16-QAM.txt| parallel --pipe -N1 --results 'best-16/{#}' --line-buffer --progress ../src/best-parity 4 4 16 16-QAM.txt -
