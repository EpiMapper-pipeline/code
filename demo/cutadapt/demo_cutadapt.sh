# This demo shows how you can use cutadapt to clean your FASTQ reads

# -f:  input FASTQ folder, not named fastq

# -fa: forward adapter

# -ba: backwards adaptor

# -o: output folder



python cutadapt.py -f raw -fa CTGTCTCTTATACACATCT -ba AGATGTGTATAAGAGACAG -o out 


