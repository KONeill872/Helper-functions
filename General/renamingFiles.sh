# To edit multiple file names using unix command 

for file in *;
do 
mv "file" "${file/old/new}"
done

# Example: Changing BAM file name. Here we add the prefix SeqBatch_A to all listed files. The #* signifies everything before '_ID'.
# Remember to escape special characters e.g. \.
for f in *.sort.bam; do mv "$f" "${f/#*_ID/SeqBatch_A_ID}";done
