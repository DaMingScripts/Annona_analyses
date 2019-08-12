#!/bin/sh


cat $1 | while read line
do
	scaffold=$(echo $line | cut -d" " -f1 | cut -d"." -f3)
	nb_genes=$(grep $scaffold $1 | wc -l)
	echo $scaffold" "$nb_genes >> $1".tmp"
done

sort $1".tmp" | uniq > $1"_genes-per-scaffold-count.txt"
rm $1".tmp"


