

#find most frequent gene. 
cd ../../analysis/pairedtag_rna/age_diff/
cat *.txt |cut -f1 -d ' '|sort|uniq -c |sort -k1,1nr > diff_count.tsv

