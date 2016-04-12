Balaur
=====
Privacy-preserving read alignment for hybrid clouds using MinHash and kmer voting.

####Commands:  
1.```index``` construct the MinHash reference genome index   
```balaur index [options] <seq_fasta>```   

2.```align``` align reads  
```balaur align [options] <seq_fasta> <reads_fastq>```  

#####MinHash options:  
```-h <arg>``` length of the MinHash fingerprint (default: 128)  
```-T <arg> ``` number of hash tables (default: 78)  
```-k <arg> ``` length of the sequence kmers (default: 16)  
```-b <arg> ``` length of the fingerprint projections (default: 2)  
```-H <arg> ``` [index-only] upper bound on kmer occurrence in the reference (default: 800) 
```-w <arg> ``` [index-only] length of the reference windows to hash (should be set to the expected read length for optimal results) (default: 150) 

#####Alignment options:  
```-m <arg>```  candidate contig filtering: min required number of index buckets shared with the read  (default: 1)  
```-N <arg> ``` candidate contig filtering: max distance from the best found number of buckets shared with a contig (default: 20)  
```-v <arg> ``` length the voting kmers (default: 20)  
```-d <arg> ``` votes array convolution radius (default: 10)  
```-x <arg> ``` min required distance between separate mapping hits (default: 40)  
```-c <arg> ``` cutoff on the min number of votes (default: 0)  
```-f <arg> ``` mapq scaling factor (default: 50)  
```-I <arg> ``` voting contig kmer sampling rate (default: 3)  

#####Privacy-related options:
```-V ```  enable vanilla mode (non-cryptographic hashing, no repeat filtering)  
```-B <arg> ```  voting kmer discretized position range (default: 20)  
```-S <arg> ```  voting task batching: number of contigs per read encrypted with same keys (default: 1)  
```-M ```  [recommended] enable masking repeat kmer neighbors

#####Other options:  
```-t <arg> ``` number of threads (default: 1) 
