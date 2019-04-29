head -n 384 "barcodes_nextSeq_comp.txt" > barcodes_hiSeq2500.txt
tail -n +385 barcodes_nextSeq.txt >> barcodes_hiSeq2500.txt
