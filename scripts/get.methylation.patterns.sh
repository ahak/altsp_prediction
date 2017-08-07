#!/bin/bash

dropbox="/home/anna/Dropbox/AltSp"

analysis="$dropbox/analysis"
apps_dir="$dropbox/apps"
common_dir="$analysis/common_files"
FSOM=$analysis/../apps/FSOM/FSOM_linux_X86_64_V0

# clname="A427"
# covfile="/home/anna/Dropbox/AltSp/data/A427_[BSseq]/DRR016652_pe.bedGraph.gz.bismark.zero.cov"
# transcript="/home/anna/Dropbox/AltSp/data/A427_RNAseq/transcripts.gtf"
for i in "$@"
do
    case $i in
	-c=*|--clname=*)
	    clname="${i#*=}"
	    shift # past argument=value
	    ;;
	*)
	    otherargs="${i}"
	    echo -e "#### Arguments not recognized:\n$otherargs"
	    echo "Proceeding nevertheless..."
            # unknown option
	    ;;
    esac
done


covfile=$(ls $analysis/../data/$clname/bsseq/*.bedGraph.gz.bismark.zero.cov)
transcript="$analysis/../data/$clname/rnaseq/transcripts.gtf"

mkdir $analysis/$clname
cd $analysis/$clname

awk '$4 >= 90' $covfile  > _$clname.bs.positive.bed
awk '$4 <= 10' $covfile  > _$clname.bs.negative.bed

ln -s $analysis/common_files/gencode.eoi.20nt.bins.sorted ./
ln -s $analysis/common_files/bins_wo_CpG.bed ./

cut -f1-6 gencode.eoi.20nt.bins.sorted | sort | uniq | sort -k1,1 -k2,2n > gencode.bins

sort -k1,1 -k2,2n _$clname.bs.positive.bed > _$clname.bs.positive.sorted
bedtools intersect -c -sorted -a gencode.bins   -b _$clname.bs.positive.sorted | awk '$7 > 0'  > _meth.positive.bed

sort -k1,1 -k2,2n _$clname.bs.negative.bed > _$clname.bs.negative.sorted
bedtools intersect -c -sorted  -a  gencode.bins  -b _$clname.bs.negative.sorted | awk '$7 > 0' > _meth.negative.bed

awk 'BEGIN{OFS="\t"}FNR == NR {a[$4]++; print; next} {if (! ($4 in a)) {$7 = 0; print} }' _meth.positive.bed _meth.negative.bed > _meth.condensed_red.bed
cut -f1-6,8 bins_wo_CpG.bed | sort | uniq >> _meth.condensed_red.bed

sort _meth.condensed_red.bed | uniq > tmp && mv tmp _meth.condensed_red.bed


awk 'BEGIN{OFS="\t"}{a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6] += $7;} END{for (i in a) print i, a[i] }' _meth.condensed_red.bed | sort -k4,4 > _meth.condensed.bed


awk 'BEGIN{OFS="\t"}{sub(/:[0-9]*$/, "", $4); if ($4 in a) a[$4] += 1; else a[$4] = 1}END{for (var in a) print var, a[var];}' _meth.condensed.bed | sort -nrk 2,2  | awk '$2 >= 10 ' > _ends.selected.dat

awk 'BEGIN{OFS = "\t"}FNR == NR{a[$4] = $0; next} {output = ""; for (i = 0; i <21; i++) {id = $1":"i; if (id in a) {split(a[id], bed); output=output"\t"bed[7]} else output=output"\t""NA" }  print $1, output  }' _meth.condensed.bed _ends.selected.dat   > FSOM.dat

# awk '
# BEGIN{OFS="\t"}
# FNR==NR { 
#     split($5, ends, "|")
#     fullname=$1":"$2":"$3":"$NF
#     map[ends[1]] = fullname
#     map[ends[2]] = fullname
#     next
# }
# {
#     split($1, name, ":")
#     endkey = name[1]":"name[2]":"name[3]":"name[4]
#     if (map[endkey]) {
# 	$1=map[endkey]":"name[5];
# 	print;
#     } else {
# 	print endkey "not found" > "/dev/stderr" 
#     }
    
# }' ../common_files/gencode.exon.end.mapping.bed FSOM.dat > FSOM_names.adjusted.dat

# cut -f1 FSOM_names.adjusted.dat | sort | uniq -c | awk 'FNR == NR && $1 > 1 { a[$2]++; next }FNR != NR {if (! a[$1]) print}' -  FSOM_names.adjusted.dat > FSOM.mapped.dat
# rm FSOM_names.adjusted.dat

## awk 'BEGIN{OFS = "\t"}FNR == NR{a[$4] = $0; next} {output = ""; for (i = 0; i <21; i++) {id = $1":"i; if (id in a) {split(a[id], bed); output=output"\t"bed[7]} else output=output"\t""NA" } print ">" $1; print "DSQ", output  }' _meth.condensed.bed _ends.selected.dat   > FSOM.input


Rscript -<<EOF
library(impute)
fsom = read.table('FSOM.dat')
rownames(fsom) = fsom\$V1
fsom\$V1 = NULL
fsom.pos = fsom[ rowSums(fsom, na.rm = T) != 0, ]
fsom.empty = fsom[ rowSums(fsom, na.rm = T) == 0, ]
fsom.empty[is.na(fsom.empty)] = 0

nr = nrow(fsom.pos)
n = 5
splitted = split(fsom.pos, sample(rep(1:n, each=ceiling(nr/n), length.out=nr)))

imputed = do.call(rbind, lapply(splitted, function(x) impute.knn(as.matrix(x))\$data))

# imputed = impute.knn(as.matrix(fsom.pos))\$data


write.table(rbind(imputed, fsom.empty), file = "fsom.imputed.dat", quote = F, row.names = T, col.names=F)
EOF


#awk 'BEGIN{OFS = "\t"}{print "> "$1; $1 = "DSQ"; print }' fsom.imputed.dat > fsom.imputed.input

awk 'BEGIN{OFS="\t"} $3 == "exon"{ gsub(/[";]/, "",  $10); print $1, $4-1, $5, $10, 0, $7}' $transcript | sort -k1,1 -k2,2n  > _$clname.exons.bed

awk 'BEGIN{OFS="\t"} $3 == "transcript" {gsub(/[";]/, "",  $10); print $1, $4-1, $5, $10}' $transcript > _$clname.transcripts.bed

awk 'BEGIN{OFS="\t"}{split($1, end,":"); print end[1], end[2], end[3], $1, 0,  end[4]}' fsom.imputed.dat > fsom.ends.bed



get420 _$clname.exons.bed > $clname.exon.ends.bed #### produces 420bp regions similar to fsom 


bedtools intersect  -u -a fsom.ends.bed -b _$clname.transcripts.bed > _covered.fsom.ends.bed


awk 'BEGIN{OFS="\t"}FNR==NR{a[$1":"$2":"$3":"$6]++; next}{if (a[$1":"$2":"$3":"$6]) print $1, $2, $3, $6, 1; else print $1, $2, $3, $6, 0}' $clname.exon.ends.bed  _covered.fsom.ends.bed > end.presence.dat

