#!/bin/bash

common_dir="/home/anna/Dropbox/AltSp/analysis/common_files"
cwd=$(pwd)
cd $common_dir

awk '$10 > 1{print $4}' Gencode24.bed | awk 'FNR==NR{a[$0]++; next}{split($4, exonname, "_"); if (a[exonname[1]]) print}' - gencode.exons.bed > gencode.multiexons.bed

awk ' BEGIN{OFS="\t"}
$1~/^chr(X|Y|[0-9]*)$/ && $3 - $2 > 50  {
    if ($6 == "+") {
	left="start"
	right="end"
    } else {
	left="end"
	right="start"
    }
    name1=$1":"$2-400":"$2+20":"$6
    print $1, $2-400, $2+20, name1":"left, 0, $6

    name2=$1":"$3-20":"$3+400":"$6
    print $1, $3-20, $3+400, name2":"right, 0, $6 
} ' gencode.multiexons.bed | sort -k1,1 -k2,2n -k4,4 | uniq > _gencode.exons.p20m400.bed


awk ' BEGIN{OFS="\t"}
$1~/^chr(X|Y|[0-9]*)$/ && $3 - $2 > 50  {
    if ($6 == "+") {
	left="start"
	right="end"
    } else {
	left="end"
	right="start"
    }
    name1=$1":"$2-400":"$2+20":"$6

    name2=$1":"$3-20":"$3+400":"$6

    $5=name1"|"name2
    print $0

} ' gencode.exons.bed > gencode.exon.end.mapping.bed



awk 'BEGIN{OFS="\t"}$3 - $2 < 400 {print $1, $2, $3}' gencode.introns.bed | sort | uniq > _gencode.introns.short400.bed 

bedtools intersect -v -a _gencode.exons.p20m400.bed -b _gencode.introns.short400.bed  > _gencode.eoi.p20m400.bed  ###
######## excludes exon ends spanning short introns


awk 'BEGIN{OFS="\t"}
{ 
for (i = 0; i<21; i++) {
  if ($6 == "+")
      print $1, $2 + i * 20, $2 + (i+1) * 20, $4":"i, $5, $6
  else
      print $1, $3 - (i+1) * 20, $3 - i * 20, $4":"i, $5, $6
  }
}
' _gencode.eoi.p20m400.bed > gencode.eoi.20nt.bins.bed

 sort -k1,1 -k2,2n  gencode.eoi.20nt.bins.bed > gencode.eoi.20nt.bins.sorted
rm gencode.eoi.20nt.bins.bed _gencode*


awk '
BEGIN{OFS="\t"}
FNR==NR {
    split($4, exname, "_")
    if ( a[exname[1]] < exname[3] ) {
	a[exname[1]] = exname[3]
    }
    next
}
{
    split($4,exname, "_")
    if (exname[3] == 0 || exname[3] == a[exname[1]]) 
	next
    else
	print 
}
' gencode.exons.bed gencode.exons.bed > gencode.middle.exons.bed

awk ' BEGIN{OFS="\t"}
$1~/^chr(X|Y|[0-9]*)$/ && $3 - $2 > 50  {

    if ($6 == "+") {
	left="start"
	right="end"
    } else {
	left="end"
	right="start"
    }
    name1=$1":"$2":"$3":"$6
    print $1, $2, $2+1, name1":"left, 0, $6

    name2=$1":"$2":"$3":"$6
    print $1, $3-1, $3, name2":"right, 0, $6 
} ' gencode.middle.exons.bed | sort -k1,1 -k2,2n -k4,4 | uniq > gencode.exon.ends.bed
    
# awk ' BEGIN{OFS="\t"}
# $1~/^chr(X|Y|[0-9]*)$/ && $3 - $2 > 50  {

#     if ($6 == "+") {
# 	left="start"
# 	right="end"
#     } else {
# 	left="end"
# 	right="start"
#     }
#     name1=$1":"$2":"$3":"$6
#     print $1, $2, $2+1, name1":"left, 0, $6

#     name2=$1":"$2":"$3":"$6
#     print $1, $3-1, $3, name2":"right, 0, $6 
# } ' gencode.exons.bed | sort -k1,1 -k2,2n -k4,4 | uniq > gencode.exon.ends.bed

cd $cwd
