# Run FastQC for .fastq files to view their quality metrics pre-trimming

# command:
fastqc 136x2_R1.fastq 136x2_R2.fastq -o ../analysis/FastQC_out >
../analysis/analysis_log.txt

# output:
null
null
Analysis complete for 136x2_R1.fastq
Analysis complete for 136x2_R2.fastq



# Run Trimmomatic on the .fastq files to remove low-quality bases

# command:
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
136x2_R1.fastq 136x2_R2.fastq \
136x2_R1-pe.fastq /dev/null \
136x2_R2-pe.fastq /dev/null \
ILLUMINACLIP:/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 \
SLIDINGWINDOW:4:15 MINLEN:50 |& tee -a ../analysis/analysis_log.txt

# output:
Picked up JAVA_TOOL_OPTIONS: -Xmx2g
TrimmomaticPE: Started with arguments:
 136x2_R1.fastq 136x2_R2.fastq 136x2_R1-pe.fastq /dev/null 136x2_R2-pe.fastq /dev/null ILLUMINACLIP:/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
java.io.FileNotFoundException: /TruSeq3-PE-2.fa (No such file or directory)
	at java.base/java.io.FileInputStream.open0(Native Method)
	at java.base/java.io.FileInputStream.open(FileInputStream.java:216)
	at java.base/java.io.FileInputStream.<init>(FileInputStream.java:157)
	at org.usadellab.trimmomatic.fasta.FastaParser.parse(FastaParser.java:54)
	at org.usadellab.trimmomatic.trim.IlluminaClippingTrimmer.loadSequences(IlluminaClippingTrimmer.java:110)
	at org.usadellab.trimmomatic.trim.IlluminaClippingTrimmer.makeIlluminaClippingTrimmer(IlluminaClippingTrimmer.java:71)
	at org.usadellab.trimmomatic.trim.TrimmerFactory.makeTrimmer(TrimmerFactory.java:32)
	at org.usadellab.trimmomatic.Trimmomatic.createTrimmers(Trimmomatic.java:59)
	at org.usadellab.trimmomatic.TrimmomaticPE.run(TrimmomaticPE.java:552)
	at org.usadellab.trimmomatic.Trimmomatic.main(Trimmomatic.java:80)
Quality encoding detected as phred33
Input Read Pairs: 1140172 Both Surviving: 1090924 (95.68%) Forward Only Surviving: 32620 (2.86%) Reverse Only Surviving: 10959 (0.96%) Dropped: 5669 (0.50%)
TrimmomaticPE: Completed successfully



# Run FastQC for trimmed .fastq files to vew their quality metrics

# command: 
fastqc 136x2_R1-pe.fastq 136x2_R2-pe.fastq -o ../analysis/FastQC_out >> ../analysis/analysis_log.txt

# output:
null
null
Analysis complete for 136x2_R1-pe.fastq
Analysis complete for 136x2_R2-pe.fastq



# Run spades assembly on the trimmed fastq files in order to assemble the genome via short-read assembly

# command:
module load spades

spades.py -1 ./136x2_R1_pe.fastq \
          -2 ./136x2_R2_pe.fastq \
          -o ../analysis/spades_genome_assembly >> ../analysis/analysis_log.txt

# output:
== Warning ==  No assembly mode was specified! If you intend to assemble high-coverage multi-cell/isolate data, use '--isolate' option.


Command line: /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/spades/4.0.0/bin/spades.py	-1	/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/data/136x2_R1-pe.fastq	-2	/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/data/136x2_R2-pe.fastq	-o	/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly	

System information:
  SPAdes version: 4.0.0
  Python version: 3.11.5
  OS: Linux-3.10.0-1160.80.1.el7.x86_64-x86_64-Intel-R-_Xeon-R-_Gold_6238_CPU_@_2.10GHz-with-glibc2.37

Output dir: /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly
Mode: read error correction and assembling
Debug mode is turned OFF

Dataset parameters:
  Standard mode
  For multi-cell/isolate data we recommend to use '--isolate' option; for single-cell MDA data use '--sc'; for metagenomic data use '--meta'; for RNA-Seq use '--rna'.
  Reads:
    Library number: 1, library type: paired-end
      orientation: fr
      left reads: ['/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/data/136x2_R1-pe.fastq']
      right reads: ['/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/data/136x2_R2-pe.fastq']
      interlaced reads: not specified
      single reads: not specified
      merged reads: not specified
Read error correction parameters:
  Iterations: 1
  PHRED offset will be auto-detected
  Corrected reads will be compressed
Assembly parameters:
  k: automatic selection based on read length
  Repeat resolution is enabled
  Mismatch careful mode is turned OFF
  MismatchCorrector will be SKIPPED
  Coverage cutoff is turned OFF
  Assembly graph output will use GFA v1.2 format
Other parameters:
  Dir for temp files: /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/tmp
  Threads: 16
  Memory limit (in Gb): 187


======= SPAdes pipeline started. Log can be found here: /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/spades.log

/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/data/136x2_R1-pe.fastq: max reads length: 151
/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/data/136x2_R2-pe.fastq: max reads length: 151

Reads length: 151

Default k-mer sizes were set to [21, 33, 55, 77] because estimated read length (151) is equal to or greater than 150

===== Before start started. 


===== Read error correction started. 


===== Read error correction started. 


== Running: /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/spades/4.0.0/bin/spades-hammer /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/corrected/configs/config.info

  0:00:00.000     1M / 16M   INFO    General                 (main.cpp                  :  76)   Starting BayesHammer, built from N/A, git revision N/A
  0:00:00.025     1M / 16M   INFO    General                 (main.cpp                  :  77)   Loading config from "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/corrected/configs/config.info"
  0:00:00.036     1M / 16M   INFO    General                 (main.cpp                  :  79)   Maximum # of threads to use (adjusted due to OMP capabilities): 16
  0:00:00.042     1M / 16M   INFO    General                 (memory_limit.cpp          :  55)   Memory limit set to 187 Gb
  0:00:00.051     1M / 16M   INFO    General                 (main.cpp                  :  87)   Trying to determine PHRED offset
  0:00:00.060     1M / 16M   INFO    General                 (main.cpp                  :  93)   Determined value is 33
  0:00:00.067     1M / 16M   INFO    General                 (hammer_tools.cpp          :  40)   Hamming graph threshold tau=1, k=21, subkmer positions = [ 0 10 ]
  0:00:00.076     1M / 16M   INFO    General                 (main.cpp                  : 114)   Size of aux. kmer data 24 bytes
     === ITERATION 0 begins ===
  0:00:00.089     1M / 16M   INFO    General                 (kmer_index_builder.hpp    : 258)   Splitting kmer instances into 16 files using 16 threads. This might take a while.
  0:00:00.096     1M / 16M   INFO    General                 (file_limit.hpp            :  43)   Open file limit set to 51200
  0:00:00.110     1M / 16M   INFO    General                 (kmer_splitter.hpp         :  94)   Memory available for splitting buffers: 3.89583 Gb
  0:00:00.119     1M / 16M   INFO    General                 (kmer_splitter.hpp         : 102)   Using cell size of 4194304
  0:00:00.130  9217M / 9217M INFO   K-mer Splitting          (kmer_data.cpp             :  98)   Processing "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/data/136x2_R1-pe.fastq"
  0:00:04.530  9217M / 9217M INFO   K-mer Splitting          (kmer_data.cpp             : 108)   Processed 1090924 reads
  0:00:04.549  9217M / 9217M INFO   K-mer Splitting          (kmer_data.cpp             :  98)   Processing "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/data/136x2_R2-pe.fastq"
  0:00:09.283  9217M / 9217M INFO   K-mer Splitting          (kmer_data.cpp             : 108)   Processed 2181848 reads
  0:00:09.292  9217M / 9217M INFO   K-mer Splitting          (kmer_data.cpp             : 113)   Total 2181848 reads processed
  0:00:09.585     1M / 4243M INFO    General                 (kmer_index_builder.hpp    : 264)   Starting k-mer counting.
  0:00:10.001     1M / 4243M INFO    General                 (kmer_index_builder.hpp    : 275)   K-mer counting done. There are 32246620 kmers in total.
  0:00:10.020     1M / 4243M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 410)   Building perfect hash indices
  0:00:10.668    24M / 4243M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 446)   Index built. Total 32246620 kmers, 23302168 bytes occupied (5.78099 bits per kmer).
  0:00:10.668    24M / 4243M INFO   K-mer Counting           (kmer_data.cpp             : 355)   Arranging kmers in hash map order
  0:00:11.200   520M / 4243M INFO    General                 (main.cpp                  : 149)   Clustering Hamming graph.
  0:00:24.963   520M / 4243M INFO    General                 (main.cpp                  : 156)   Extracting clusters:
  0:00:24.963   520M / 4243M INFO    General                 (concurrent_dsu.cpp        :  19)   Connecting to root
  0:00:25.027   520M / 4243M INFO    General                 (concurrent_dsu.cpp        :  35)   Calculating counts
  0:00:31.423   792M / 4243M INFO    General                 (concurrent_dsu.cpp        :  64)   Writing down entries
  0:00:40.402   520M / 4243M INFO    General                 (main.cpp                  : 168)   Clustering done. Total clusters: 13438029
  0:00:40.438   272M / 4243M INFO   K-mer Counting           (kmer_data.cpp             : 372)   Collecting K-mer information, this takes a while.
  0:00:40.837  1012M / 4243M INFO   K-mer Counting           (kmer_data.cpp             : 378)   Processing "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/data/136x2_R1-pe.fastq"
  0:00:49.777  1012M / 4243M INFO   K-mer Counting           (kmer_data.cpp             : 378)   Processing "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/data/136x2_R2-pe.fastq"
  0:00:58.414  1012M / 4243M INFO   K-mer Counting           (kmer_data.cpp             : 385)   Collection done, postprocessing.
  0:00:58.556  1012M / 4243M INFO   K-mer Counting           (kmer_data.cpp             : 398)   There are 32246620 kmers in total. Among them 20771302 (64.4139%) are singletons.
  0:00:58.556  1012M / 4243M INFO    General                 (main.cpp                  : 174)   Subclustering Hamming graph
  0:01:02.796  1012M / 4243M INFO   Hamming Subclustering    (kmer_cluster.cpp          : 637)   Subclustering done. Total 6 non-read kmers were generated.
  0:01:02.802  1012M / 4243M INFO   Hamming Subclustering    (kmer_cluster.cpp          : 638)   Subclustering statistics:
  0:01:02.810  1012M / 4243M INFO   Hamming Subclustering    (kmer_cluster.cpp          : 639)     Total singleton hamming clusters: 5833290. Among them 2321756 (39.8018%) are good
  0:01:02.814  1012M / 4243M INFO   Hamming Subclustering    (kmer_cluster.cpp          : 640)     Total singleton subclusters: 25929. Among them 25752 (99.3174%) are good
  0:01:02.823  1012M / 4243M INFO   Hamming Subclustering    (kmer_cluster.cpp          : 641)     Total non-singleton subcluster centers: 7688675. Among them 5661220 (73.6306%) are good
  0:01:02.835  1012M / 4243M INFO   Hamming Subclustering    (kmer_cluster.cpp          : 642)     Average size of non-trivial subcluster: 3.43536 kmers
  0:01:02.839  1012M / 4243M INFO   Hamming Subclustering    (kmer_cluster.cpp          : 643)     Average number of sub-clusters per non-singleton cluster: 1.01445
  0:01:02.847  1012M / 4243M INFO   Hamming Subclustering    (kmer_cluster.cpp          : 644)     Total solid k-mers: 8008728
  0:01:02.856  1012M / 4243M INFO   Hamming Subclustering    (kmer_cluster.cpp          : 645)     Substitution probabilities: (     0.940017    0.0309646   0.00544922    0.0235687 )
(    0.0244424     0.959171    0.0126699   0.00371663 )
(    0.0037001    0.0126554     0.959181    0.0244638 )
(    0.0235466    0.0054303    0.0309437     0.940079 )
  0:01:02.883  1012M / 4243M INFO    General                 (main.cpp                  : 179)   Finished clustering.
  0:01:02.891  1012M / 4243M INFO    General                 (main.cpp                  : 198)   Starting solid k-mers expansion in 16 threads.
  0:01:07.523  1012M / 4243M INFO    General                 (main.cpp                  : 219)   Solid k-mers iteration 0 produced 86740 new k-mers.
  0:01:11.783  1012M / 4243M INFO    General                 (main.cpp                  : 219)   Solid k-mers iteration 1 produced 5825 new k-mers.
  0:01:15.968  1012M / 4243M INFO    General                 (main.cpp                  : 219)   Solid k-mers iteration 2 produced 191 new k-mers.
  0:01:20.198  1012M / 4243M INFO    General                 (main.cpp                  : 219)   Solid k-mers iteration 3 produced 8 new k-mers.
  0:01:20.198  1012M / 4243M INFO    General                 (main.cpp                  : 223)   Solid k-mers finalized
  0:01:20.198  1012M / 4243M INFO    General                 (hammer_tools.cpp          : 226)   Starting read correction in 16 threads.
  0:01:20.198  1012M / 4243M INFO    General                 (hammer_tools.cpp          : 239)   Correcting pair of reads: "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/data/136x2_R1-pe.fastq" and "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/data/136x2_R2-pe.fastq"
  0:01:23.463  2108M / 4243M INFO    General                 (hammer_tools.cpp          : 173)   Prepared batch 0 of 1090924 reads.
  0:01:29.143  2202M / 4243M INFO    General                 (hammer_tools.cpp          : 180)   Processed batch 0
  0:01:30.585  2202M / 4243M INFO    General                 (hammer_tools.cpp          : 190)   Written batch 0
  0:01:31.451  1012M / 4243M INFO    General                 (hammer_tools.cpp          : 280)   Correction done. Changed 832176 bases in 455939 reads.
  0:01:31.481  1012M / 4243M INFO    General                 (hammer_tools.cpp          : 281)   Failed to correct 59 bases out of 319521289.
  0:01:31.603     1M / 4243M INFO    General                 (main.cpp                  : 256)   Saving corrected dataset description to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/corrected/corrected.yaml
  0:01:31.633     1M / 4243M INFO    General                 (main.cpp                  : 263)   All done. Exiting.

===== Read error correction finished. 


===== corrected reads compression started. 


== Running: /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/python/3.11.5/bin/python3 /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/spades/4.0.0/share/spades/spades_pipeline/scripts/compress_all.py --input_file /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/corrected/corrected.yaml --ext_python_modules_home /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/spades/4.0.0/share/spades --max_threads 16 --output_dir /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/corrected --gzip_output

== Compressing corrected reads (with pigz)
== Files to compress: ['/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/corrected/136x2_R1-pe00.0_0.cor.fastq', '/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/corrected/136x2_R2-pe00.0_0.cor.fastq', '/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/corrected/136x2_R_unpaired00.0_0.cor.fastq']
== Files compression is finished
== Dataset yaml file is updated

===== corrected reads compression finished. 


===== Read error correction finished. 


===== Assembling started. 


===== K21 started. 


== Running: /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/spades/4.0.0/bin/spades-core /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K21/configs/config.info

  0:00:00.000     1M / 16M   INFO    General                 (main.cpp                  :  94)   Loaded config from "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K21/configs/config.info"
  0:00:00.028     1M / 16M   INFO    General                 (memory_limit.cpp          :  55)   Memory limit set to 187 Gb
  0:00:00.037     1M / 16M   INFO    General                 (main.cpp                  : 102)   Starting SPAdes, built from N/A, git revision N/A
  0:00:00.049     1M / 16M   INFO    General                 (main.cpp                  : 103)   Maximum k-mer length: 128
  0:00:00.057     1M / 16M   INFO    General                 (main.cpp                  : 104)   Assembling dataset ("/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/dataset.info") with K=21
  0:00:00.066     1M / 16M   INFO    General                 (main.cpp                  : 105)   Maximum # of threads to use (adjusted due to OMP capabilities): 16
  0:00:00.074     1M / 16M   INFO    General                 (pipeline.cpp              : 212)   SPAdes started
  0:00:00.082     1M / 16M   INFO    General                 (pipeline.cpp              : 225)   Starting from stage: read_conversion
  0:00:00.091     1M / 16M   INFO    General                 (pipeline.cpp              : 234)   Two-step repeat resolution disabled
  0:00:00.100     1M / 16M   INFO   GraphCore                (graph_core.hpp            : 689)   Graph created, vertex min_id: 3, edge min_id: 3
  0:00:00.109     1M / 16M   INFO   GraphCore                (graph_core.hpp            : 690)   Vertex size: 48, edge size: 40
  0:00:00.118     1M / 16M   INFO    General                 (edge_index.hpp            : 132)   Size of edge index entries: 12/8
  0:00:00.127     1M / 16M   INFO   StageManager             (stage.cpp                 : 189)   STAGE == Binary Read Conversion (id: read_conversion)
  0:00:00.141     1M / 16M   INFO    General                 (read_converter.cpp        :  78)   Converting reads to binary format for library #0 (takes a while)
  0:00:00.149     1M / 16M   INFO    General                 (read_converter.cpp        :  99)   Converting paired reads
  0:00:00.567    82M / 121M  INFO    General                 (binary_converter.cpp      : 127)   16384 reads processed
  0:00:00.598    82M / 125M  INFO    General                 (binary_converter.cpp      : 127)   32768 reads processed
  0:00:00.644    82M / 130M  INFO    General                 (binary_converter.cpp      : 127)   65536 reads processed
  0:00:00.814    82M / 162M  INFO    General                 (binary_converter.cpp      : 127)   131072 reads processed
  0:00:01.181    52M / 204M  INFO    General                 (binary_converter.cpp      : 127)   262144 reads processed
  0:00:01.790     0M / 259M  INFO    General                 (binary_converter.cpp      : 127)   524288 reads processed
  0:00:02.992     0M / 338M  INFO    General                 (binary_converter.cpp      : 127)   1048576 reads processed
  0:00:03.281     0M / 338M  INFO    General                 (binary_converter.cpp      : 143)   1089342 reads written
  0:00:03.308     0M / 338M  INFO    General                 (read_converter.cpp        : 113)   Converting single reads
  0:00:03.349     0M / 338M  INFO    General                 (binary_converter.cpp      : 143)   1519 reads written
  0:00:03.362     0M / 338M  INFO    General                 (read_converter.cpp        : 119)   Converting merged reads
  0:00:03.377     0M / 338M  INFO    General                 (binary_converter.cpp      : 143)   0 reads written
  0:00:03.538     1M / 338M  INFO   StageManager             (stage.cpp                 : 189)   STAGE == de Bruijn graph construction (id: construction)
  0:00:03.612     1M / 338M  INFO    General                 (construction.cpp          : 150)   Max read length 151
  0:00:03.623     1M / 338M  INFO    General                 (construction.cpp          : 156)   Average read length 146.498
  0:00:03.635     1M / 338M  INFO    General                 (stage.cpp                 : 121)   PROCEDURE == k+1-mer counting (id: construction:kpomer_counting)
  0:00:03.644     1M / 338M  INFO    General                 (kmer_index_builder.hpp    : 258)   Splitting kmer instances into 160 files using 16 threads. This might take a while.
  0:00:03.657     1M / 338M  INFO    General                 (file_limit.hpp            :  43)   Open file limit set to 51200
  0:00:03.670     1M / 338M  INFO    General                 (kmer_splitter.hpp         :  94)   Memory available for splitting buffers: 3.89582 Gb
  0:00:03.678     1M / 338M  INFO    General                 (kmer_splitter.hpp         : 102)   Using cell size of 419430
  0:00:06.371  9601M / 9601M INFO    General                 (kmer_splitters.hpp        : 128)   Processed 4360406 reads
  0:00:06.383     1M / 2466M INFO    General                 (kmer_splitters.hpp        : 134)   Used 4360406 reads
  0:00:06.463     1M / 2466M INFO    General                 (kmer_index_builder.hpp    : 264)   Starting k-mer counting.
  0:00:07.943     1M / 2466M INFO    General                 (kmer_index_builder.hpp    : 275)   K-mer counting done. There are 5446398 kmers in total.
  0:00:07.955     1M / 2466M INFO    General                 (stage.cpp                 : 121)   PROCEDURE == Extension index construction (id: construction:extension_index_construction)
  0:00:08.042     1M / 2466M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 453)   Building kmer index
  0:00:08.050     1M / 2466M INFO    General                 (kmer_index_builder.hpp    : 258)   Splitting kmer instances into 160 files using 16 threads. This might take a while.
  0:00:08.064     2M / 2466M INFO    General                 (file_limit.hpp            :  43)   Open file limit set to 51200
  0:00:08.077     2M / 2466M INFO    General                 (kmer_splitter.hpp         :  94)   Memory available for splitting buffers: 3.89581 Gb
  0:00:08.085     2M / 2466M INFO    General                 (kmer_splitter.hpp         : 102)   Using cell size of 419430
  0:00:09.777  9602M / 9602M INFO    General                 (kmer_splitters.hpp        : 197)   Processed 5446398 kmers
  0:00:09.796  9602M / 9602M INFO    General                 (kmer_splitters.hpp        : 202)   Used 5446398 kmers.
  0:00:09.806     2M / 2466M INFO    General                 (kmer_index_builder.hpp    : 264)   Starting k-mer counting.
  0:00:11.490     2M / 2466M INFO    General                 (kmer_index_builder.hpp    : 275)   K-mer counting done. There are 5433771 kmers in total.
  0:00:11.502     2M / 2466M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 410)   Building perfect hash indices
  0:00:11.611     7M / 2466M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 446)   Index built. Total 5433771 kmers, 4046504 bytes occupied (5.95756 bits per kmer).
  0:00:11.611     7M / 2466M INFO    General                 (kmer_index_builder.hpp    : 168)   Merging final buckets.
  0:00:12.822    12M / 2466M INFO   DeBruijnExtensionIndexBu (kmer_extension_index_build: 101)   Building k-mer extensions from k+1-mers
  0:00:13.033    12M / 2466M INFO   DeBruijnExtensionIndexBu (kmer_extension_index_build: 106)   Building k-mer extensions from k+1-mers finished.
  0:00:13.045    11M / 2466M INFO    General                 (stage.cpp                 : 121)   PROCEDURE == Early tip clipping (id: construction:early_tip_clipper)
  0:00:13.050    11M / 2466M INFO    General                 (construction.cpp          : 298)   Early tip clipper length bound set as (RL - K)
  0:00:13.058    11M / 2466M INFO   Early tip clipping       (early_simplification.hpp  :  48)   Early tip clipping
  0:00:13.211    14M / 2466M INFO   Early tip clipping       (early_simplification.hpp  :  83)   #tipped junctions: 31542
  0:00:13.229    14M / 2466M INFO   Early tip clipping       (early_simplification.hpp  :  94)   Clipped tips: 31791
  0:00:13.242    11M / 2466M INFO   Early tip clipping       (early_simplification.hpp  :  50)   217289 22-mers were removed by early tip clipper
  0:00:13.247    11M / 2466M INFO    General                 (stage.cpp                 : 121)   PROCEDURE == Condensing graph (id: construction:graph_condensing)
  0:00:13.314    12M / 2466M INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 381)   Extracting unbranching paths
  0:00:13.467    14M / 2466M INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 400)   Extracting unbranching paths finished. 37413 sequences extracted
  0:00:13.606    14M / 2466M INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 336)   Collecting perfect loops
  0:00:13.721    14M / 2466M INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 369)   Collecting perfect loops finished. 0 loops collected
  0:00:13.736    14M / 2466M INFO   DeBruijnGraphConstructor (debruijn_graph_constructor: 586)   Sorting edges...
  0:00:13.767    14M / 2466M INFO   DeBruijnGraphConstructor (debruijn_graph_constructor: 588)   Edges sorted
  0:00:13.775    14M / 2466M INFO    General                 (debruijn_graph_constructor: 516)   Total 74826 edges to create
  0:00:13.789    17M / 2466M INFO    General                 (debruijn_graph_constructor: 519)   Collecting link records
  0:00:13.862    18M / 2466M INFO    General                 (debruijn_graph_constructor: 521)   Ordering link records
  0:00:13.896    18M / 2466M INFO    General                 (debruijn_graph_constructor: 524)   Sorting done
  0:00:13.905    18M / 2466M INFO    General                 (debruijn_graph_constructor: 537)   Sorting LinkRecords...
  0:00:13.938    18M / 2466M INFO    General                 (debruijn_graph_constructor: 540)   LinkRecords sorted
  0:00:13.947    18M / 2466M INFO    General                 (debruijn_graph_constructor: 542)   Total 24786 vertices to create
  0:00:13.955    21M / 2466M INFO    General                 (debruijn_graph_constructor: 545)   Connecting the graph
  0:00:14.008    19M / 2466M INFO    General                 (stage.cpp                 : 121)   PROCEDURE == Filling coverage indices (PHM) (id: construction:coverage_filling_phm)
  0:00:14.017    19M / 2466M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 410)   Building perfect hash indices
  0:00:14.115    24M / 2466M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 446)   Index built. Total 5446398 kmers, 4055960 bytes occupied (5.95764 bits per kmer).
  0:00:14.129    45M / 2466M INFO    General                 (coverage_hash_map_builder.:  49)   Collecting k-mer coverage information from reads, this takes a while.
  0:00:17.965    45M / 2466M INFO    General                 (construction.cpp          : 427)   Filling coverage and flanking coverage from PHM
  0:00:18.086    45M / 2466M INFO    General                 (coverage_filling.hpp      :  83)   Processed 74776 edges
  0:00:18.575     9M / 2466M INFO   StageManager             (stage.cpp                 : 189)   STAGE == EC Threshold Finding (id: ec_threshold_finder)
  0:00:18.577     9M / 2466M INFO    General                 (kmer_coverage_model.cpp   : 182)   Kmer coverage valley at: 4
  0:00:18.578     9M / 2466M INFO    General                 (kmer_coverage_model.cpp   : 202)   K-mer histogram maximum: 45
  0:00:18.578     9M / 2466M INFO    General                 (kmer_coverage_model.cpp   : 238)   Estimated median coverage: 45. Coverage mad: 16.3086
  0:00:18.578     9M / 2466M INFO    General                 (kmer_coverage_model.cpp   : 260)   Fitting coverage model
  0:00:18.625     9M / 2466M INFO    General                 (kmer_coverage_model.cpp   : 296)   ... iteration 2
  0:00:18.755     9M / 2466M INFO    General                 (kmer_coverage_model.cpp   : 296)   ... iteration 4
  0:00:19.229     9M / 2466M INFO    General                 (kmer_coverage_model.cpp   : 296)   ... iteration 8
  0:00:19.734     9M / 2466M INFO    General                 (kmer_coverage_model.cpp   : 296)   ... iteration 16
  0:00:20.076     9M / 2466M INFO    General                 (kmer_coverage_model.cpp   : 310)   Fitted mean coverage: 45.211. Fitted coverage std. dev: 15.627
  0:00:20.078     9M / 2466M INFO    General                 (kmer_coverage_model.cpp   : 335)   Probability of erroneous kmer at valley: 0.692352
  0:00:20.078     9M / 2466M INFO    General                 (kmer_coverage_model.cpp   : 359)   Preliminary threshold calculated as: 18
  0:00:20.078     9M / 2466M INFO    General                 (kmer_coverage_model.cpp   : 363)   Threshold adjusted to: 18
  0:00:20.078     9M / 2466M INFO    General                 (kmer_coverage_model.cpp   : 376)   Estimated genome size (ignoring repeats): 4676848
  0:00:20.078     9M / 2466M INFO    General                 (genomic_info_filler.cpp   :  56)   Mean coverage was calculated as 45.211
  0:00:20.078     9M / 2466M INFO    General                 (genomic_info_filler.cpp   :  71)   EC coverage threshold value was calculated as 18
  0:00:20.078     9M / 2466M INFO    General                 (genomic_info_filler.cpp   :  72)   Trusted kmer low bound: 0
  0:00:20.078     9M / 2466M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Raw Simplification (id: raw_simplification)
  0:00:20.079     9M / 2466M INFO    General                 (simplification.cpp        : 129)   PROCEDURE == Initial cleaning
  0:00:20.079     9M / 2466M INFO    General                 (graph_simplification.hpp  : 674)   Flanking coverage based disconnection disabled
  0:00:20.079     9M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Self conjugate edge remover
  0:00:20.084     9M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Self conjugate edge remover triggered 0 times
  0:00:20.084     9M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Initial tip clipper
  0:00:20.088     9M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Initial tip clipper triggered 300 times
  0:00:20.088     9M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Initial ec remover
  0:00:20.225     9M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Initial ec remover triggered 8226 times
  0:00:20.225     9M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Initial isolated edge remover
  0:00:20.226     9M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Initial isolated edge remover triggered 68 times
  0:00:20.226     8M / 2466M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Simplification (id: simplification)
  0:00:20.226     8M / 2466M INFO    General                 (simplification.cpp        : 397)   Graph simplification started
  0:00:20.226     8M / 2466M INFO    General                 (graph_simplification.hpp  : 646)   Creating parallel br instance
  0:00:20.226     8M / 2466M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 1
  0:00:20.226     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.227     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 62 times
  0:00:20.228     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.283     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 2022 times
  0:00:20.284     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.285     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 23 times
  0:00:20.286     8M / 2466M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 2
  0:00:20.286     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.287     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 1 times
  0:00:20.287     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.287     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 1 times
  0:00:20.287     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.292     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 309 times
  0:00:20.292     8M / 2466M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 3
  0:00:20.292     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.292     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 5 times
  0:00:20.292     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.296     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 108 times
  0:00:20.296     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.298     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 222 times
  0:00:20.298     8M / 2466M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 4
  0:00:20.299     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.299     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 3 times
  0:00:20.299     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.302     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 101 times
  0:00:20.302     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.303     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 139 times
  0:00:20.303     8M / 2466M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 5
  0:00:20.303     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.304     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 2 times
  0:00:20.304     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.306     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 64 times
  0:00:20.306     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.306     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 87 times
  0:00:20.307     8M / 2466M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 6
  0:00:20.307     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.307     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 1 times
  0:00:20.307     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.308     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 55 times
  0:00:20.308     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.309     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 67 times
  0:00:20.309     8M / 2466M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 7
  0:00:20.309     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.309     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 1 times
  0:00:20.309     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.310     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 26 times
  0:00:20.310     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.311     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 55 times
  0:00:20.311     8M / 2466M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 8
  0:00:20.311     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.311     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 1 times
  0:00:20.311     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.312     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 32 times
  0:00:20.312     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.313     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 44 times
  0:00:20.313     8M / 2466M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 9
  0:00:20.313     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.313     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 1 times
  0:00:20.313     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.313     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 21 times
  0:00:20.313     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.314     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 38 times
  0:00:20.314     8M / 2466M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 10
  0:00:20.314     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.314     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 1 times
  0:00:20.314     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.315     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 19 times
  0:00:20.315     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.315     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 39 times
  0:00:20.315     8M / 2466M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 11
  0:00:20.315     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.316     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 2 times
  0:00:20.316     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.319     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 14 times
  0:00:20.319     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.319     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 0 times
  0:00:20.319     8M / 2466M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 12
  0:00:20.319     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.319     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 0 times
  0:00:20.319     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.319     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 0 times
  0:00:20.320     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.320     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 0 times
  0:00:20.320     8M / 2466M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Simplification Cleanup (id: simplification_cleanup)
  0:00:20.320     8M / 2466M INFO    General                 (simplification.cpp        : 189)   PROCEDURE == Post simplification
  0:00:20.320     8M / 2466M INFO    General                 (graph_simplification.hpp  : 455)   Disconnection of relatively low covered edges disabled
  0:00:20.320     8M / 2466M INFO    General                 (graph_simplification.hpp  : 494)   Complex tip clipping disabled
  0:00:20.320     8M / 2466M INFO    General                 (graph_simplification.hpp  : 646)   Creating parallel br instance
  0:00:20.321     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.321     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 0 times
  0:00:20.321     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.323     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 0 times
  0:00:20.323     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.324     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 0 times
  0:00:20.324     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.326     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 0 times
  0:00:20.326     8M / 2466M INFO    General                 (simplification.cpp        : 348)   Disrupting self-conjugate edges
  0:00:20.327     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Removing isolated edges
  0:00:20.328     8M / 2466M INFO   Simplification           (parallel_processing.hpp   : 171)   Removing isolated edges triggered 0 times
  0:00:20.328     8M / 2466M INFO    General                 (simplification.cpp        : 508)   After simplification:
  0:00:20.328     8M / 2466M INFO    General                 (simplification.cpp        : 509)     Average coverage = 54.8474
  0:00:20.328     8M / 2466M INFO    General                 (simplification.cpp        : 510)     Total length = 4981171
  0:00:20.329     8M / 2466M INFO    General                 (simplification.cpp        : 511)     Median edge length: 10192
  0:00:20.329     8M / 2466M INFO    General                 (simplification.cpp        : 512)     Edges: 6255
  0:00:20.329     8M / 2466M INFO    General                 (simplification.cpp        : 513)     Vertices: 4308
  0:00:20.329     8M / 2466M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Contig Output (id: contig_output)
  0:00:20.332     8M / 2466M INFO    General                 (read_converter.cpp        : 135)   Outputting contigs to "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K21/simplified_contigs"
  0:00:20.353    12M / 2466M INFO    General                 (binary_converter.cpp      : 143)   3135 reads written
  0:00:20.389     8M / 2466M INFO    General                 (pipeline.cpp              : 292)   SPAdes finished
  0:00:20.393     1M / 2466M INFO    General                 (main.cpp                  : 131)   Assembling time: 0 hours 0 minutes 20 seconds

===== K21 finished. 


===== K33 started. 


== Running: /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/spades/4.0.0/bin/spades-core /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K33/configs/config.info

  0:00:00.000     1M / 16M   INFO    General                 (main.cpp                  :  94)   Loaded config from "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K33/configs/config.info"
  0:00:00.049     1M / 16M   INFO    General                 (memory_limit.cpp          :  55)   Memory limit set to 187 Gb
  0:00:00.078     1M / 16M   INFO    General                 (main.cpp                  : 102)   Starting SPAdes, built from N/A, git revision N/A
  0:00:00.094     1M / 16M   INFO    General                 (main.cpp                  : 103)   Maximum k-mer length: 128
  0:00:00.127     1M / 16M   INFO    General                 (main.cpp                  : 104)   Assembling dataset ("/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/dataset.info") with K=33
  0:00:00.144     1M / 16M   INFO    General                 (main.cpp                  : 105)   Maximum # of threads to use (adjusted due to OMP capabilities): 16
  0:00:00.177     1M / 16M   INFO    General                 (pipeline.cpp              : 212)   SPAdes started
  0:00:00.194     1M / 16M   INFO    General                 (pipeline.cpp              : 225)   Starting from stage: read_conversion
  0:00:00.227     1M / 16M   INFO    General                 (pipeline.cpp              : 234)   Two-step repeat resolution disabled
  0:00:00.243     1M / 16M   INFO   GraphCore                (graph_core.hpp            : 689)   Graph created, vertex min_id: 3, edge min_id: 3
  0:00:00.277     1M / 16M   INFO   GraphCore                (graph_core.hpp            : 690)   Vertex size: 48, edge size: 40
  0:00:00.293     1M / 16M   INFO    General                 (edge_index.hpp            : 132)   Size of edge index entries: 12/8
  0:00:00.326     1M / 16M   INFO   StageManager             (stage.cpp                 : 189)   STAGE == Binary Read Conversion (id: read_conversion)
  0:00:00.351     1M / 16M   INFO    General                 (read_converter.cpp        :  57)   Binary reads detected
  0:00:00.373     1M / 16M   INFO   StageManager             (stage.cpp                 : 189)   STAGE == de Bruijn graph construction (id: construction)
  0:00:00.377     1M / 16M   INFO    General                 (construction.cpp          : 115)   Contigs from previous K will be used: /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K21/simplified_contigs
  0:00:00.428     1M / 16M   INFO    General                 (construction.cpp          : 150)   Max read length 151
  0:00:00.436     1M / 16M   INFO    General                 (construction.cpp          : 156)   Average read length 146.498
  0:00:00.446     1M / 16M   INFO    General                 (stage.cpp                 : 121)   PROCEDURE == k+1-mer counting (id: construction:kpomer_counting)
  0:00:00.454     1M / 16M   INFO    General                 (kmer_index_builder.hpp    : 258)   Splitting kmer instances into 160 files using 16 threads. This might take a while.
  0:00:00.467     1M / 16M   INFO    General                 (file_limit.hpp            :  43)   Open file limit set to 51200
  0:00:00.480     1M / 16M   INFO    General                 (kmer_splitter.hpp         :  94)   Memory available for splitting buffers: 3.89582 Gb
  0:00:00.488     1M / 16M   INFO    General                 (kmer_splitter.hpp         : 102)   Using cell size of 209715
  0:00:04.734  9602M / 9602M INFO    General                 (kmer_splitters.hpp        : 128)   Processed 4366676 reads
  0:00:04.755     2M / 4347M INFO    General                 (kmer_splitters.hpp        : 134)   Used 4366676 reads
  0:00:04.902     2M / 4347M INFO    General                 (kmer_index_builder.hpp    : 264)   Starting k-mer counting.
  0:00:06.218     2M / 4347M INFO    General                 (kmer_index_builder.hpp    : 275)   K-mer counting done. There are 5574759 kmers in total.
  0:00:06.229     1M / 4347M INFO    General                 (stage.cpp                 : 121)   PROCEDURE == Extension index construction (id: construction:extension_index_construction)
  0:00:06.318     2M / 4347M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 453)   Building kmer index
  0:00:06.324     2M / 4347M INFO    General                 (kmer_index_builder.hpp    : 258)   Splitting kmer instances into 160 files using 16 threads. This might take a while.
  0:00:06.336     2M / 4347M INFO    General                 (file_limit.hpp            :  43)   Open file limit set to 51200
  0:00:06.349     2M / 4347M INFO    General                 (kmer_splitter.hpp         :  94)   Memory available for splitting buffers: 3.89581 Gb
  0:00:06.358     2M / 4347M INFO    General                 (kmer_splitter.hpp         : 102)   Using cell size of 209715
  0:00:08.291  9602M / 9602M INFO    General                 (kmer_splitters.hpp        : 197)   Processed 5574759 kmers
  0:00:08.309  9602M / 9602M INFO    General                 (kmer_splitters.hpp        : 202)   Used 5574759 kmers.
  0:00:08.320     2M / 4347M INFO    General                 (kmer_index_builder.hpp    : 264)   Starting k-mer counting.
  0:00:09.641     2M / 4347M INFO    General                 (kmer_index_builder.hpp    : 275)   K-mer counting done. There are 5565321 kmers in total.
  0:00:09.656     2M / 4347M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 410)   Building perfect hash indices
  0:00:09.780     7M / 4347M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 446)   Index built. Total 5565321 kmers, 4141592 bytes occupied (5.95343 bits per kmer).
  0:00:09.780     7M / 4347M INFO    General                 (kmer_index_builder.hpp    : 168)   Merging final buckets.
  0:00:11.351    12M / 4347M INFO   DeBruijnExtensionIndexBu (kmer_extension_index_build: 101)   Building k-mer extensions from k+1-mers
  0:00:11.543    12M / 4347M INFO   DeBruijnExtensionIndexBu (kmer_extension_index_build: 106)   Building k-mer extensions from k+1-mers finished.
  0:00:11.550    12M / 4347M INFO    General                 (stage.cpp                 : 121)   PROCEDURE == Early tip clipping (id: construction:early_tip_clipper)
  0:00:11.550    12M / 4347M INFO    General                 (construction.cpp          : 298)   Early tip clipper length bound set as (RL - K)
  0:00:11.550    12M / 4347M INFO   Early tip clipping       (early_simplification.hpp  :  48)   Early tip clipping
  0:00:11.694    14M / 4347M INFO   Early tip clipping       (early_simplification.hpp  :  83)   #tipped junctions: 34149
  0:00:11.696    14M / 4347M INFO   Early tip clipping       (early_simplification.hpp  :  94)   Clipped tips: 34420
  0:00:11.703    12M / 4347M INFO   Early tip clipping       (early_simplification.hpp  :  50)   284763 34-mers were removed by early tip clipper
  0:00:11.703    12M / 4347M INFO    General                 (stage.cpp                 : 121)   PROCEDURE == Condensing graph (id: construction:graph_condensing)
  0:00:11.752    12M / 4347M INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 381)   Extracting unbranching paths
  0:00:11.902    15M / 4347M INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 400)   Extracting unbranching paths finished. 29564 sequences extracted
  0:00:12.035    14M / 4347M INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 336)   Collecting perfect loops
  0:00:12.160    15M / 4347M INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 369)   Collecting perfect loops finished. 0 loops collected
  0:00:12.166    14M / 4347M INFO   DeBruijnGraphConstructor (debruijn_graph_constructor: 586)   Sorting edges...
  0:00:12.168    14M / 4347M INFO   DeBruijnGraphConstructor (debruijn_graph_constructor: 588)   Edges sorted
  0:00:12.168    14M / 4347M INFO    General                 (debruijn_graph_constructor: 516)   Total 59128 edges to create
  0:00:12.168    17M / 4347M INFO    General                 (debruijn_graph_constructor: 519)   Collecting link records
  0:00:12.218    18M / 4347M INFO    General                 (debruijn_graph_constructor: 521)   Ordering link records
  0:00:12.219    18M / 4347M INFO    General                 (debruijn_graph_constructor: 524)   Sorting done
  0:00:12.219    18M / 4347M INFO    General                 (debruijn_graph_constructor: 537)   Sorting LinkRecords...
  0:00:12.220    18M / 4347M INFO    General                 (debruijn_graph_constructor: 540)   LinkRecords sorted
  0:00:12.220    18M / 4347M INFO    General                 (debruijn_graph_constructor: 542)   Total 20126 vertices to create
  0:00:12.220    20M / 4347M INFO    General                 (debruijn_graph_constructor: 545)   Connecting the graph
  0:00:12.255    18M / 4347M INFO    General                 (stage.cpp                 : 121)   PROCEDURE == Filling coverage indices (PHM) (id: construction:coverage_filling_phm)
  0:00:12.255    18M / 4347M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 410)   Building perfect hash indices
  0:00:12.363    23M / 4347M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 446)   Index built. Total 5574759 kmers, 4148440 bytes occupied (5.95318 bits per kmer).
  0:00:12.375    45M / 4347M INFO    General                 (coverage_hash_map_builder.:  49)   Collecting k-mer coverage information from reads, this takes a while.
  0:00:15.985    45M / 4347M INFO    General                 (construction.cpp          : 427)   Filling coverage and flanking coverage from PHM
  0:00:16.113    45M / 4347M INFO    General                 (coverage_filling.hpp      :  83)   Processed 59110 edges
  0:00:16.649     7M / 4347M INFO   StageManager             (stage.cpp                 : 189)   STAGE == EC Threshold Finding (id: ec_threshold_finder)
  0:00:16.652     8M / 4347M INFO    General                 (kmer_coverage_model.cpp   : 182)   Kmer coverage valley at: 6
  0:00:16.652     7M / 4347M INFO    General                 (kmer_coverage_model.cpp   : 202)   K-mer histogram maximum: 39
  0:00:16.652     7M / 4347M INFO    General                 (kmer_coverage_model.cpp   : 238)   Estimated median coverage: 41. Coverage mad: 14.826
  0:00:16.652     7M / 4347M INFO    General                 (kmer_coverage_model.cpp   : 260)   Fitting coverage model
  0:00:16.708     8M / 4347M INFO    General                 (kmer_coverage_model.cpp   : 296)   ... iteration 2
  0:00:16.848     8M / 4347M INFO    General                 (kmer_coverage_model.cpp   : 296)   ... iteration 4
  0:00:17.364     8M / 4347M INFO    General                 (kmer_coverage_model.cpp   : 296)   ... iteration 8
  0:00:18.502     8M / 4347M INFO    General                 (kmer_coverage_model.cpp   : 296)   ... iteration 16
  0:00:19.522     8M / 4347M INFO    General                 (kmer_coverage_model.cpp   : 310)   Fitted mean coverage: 40.8547. Fitted coverage std. dev: 14.4534
  0:00:19.532     8M / 4347M INFO    General                 (kmer_coverage_model.cpp   : 335)   Probability of erroneous kmer at valley: 0.451448
  0:00:19.540     8M / 4347M INFO    General                 (kmer_coverage_model.cpp   : 359)   Preliminary threshold calculated as: 17
  0:00:19.548     8M / 4347M INFO    General                 (kmer_coverage_model.cpp   : 363)   Threshold adjusted to: 17
  0:00:19.552     8M / 4347M INFO    General                 (kmer_coverage_model.cpp   : 376)   Estimated genome size (ignoring repeats): 4659887
  0:00:19.561     7M / 4347M INFO    General                 (genomic_info_filler.cpp   :  56)   Mean coverage was calculated as 40.8547
  0:00:19.569     7M / 4347M INFO    General                 (genomic_info_filler.cpp   :  71)   EC coverage threshold value was calculated as 17
  0:00:19.577     7M / 4347M INFO    General                 (genomic_info_filler.cpp   :  72)   Trusted kmer low bound: 0
  0:00:19.586     7M / 4347M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Raw Simplification (id: raw_simplification)
  0:00:19.594     7M / 4347M INFO    General                 (simplification.cpp        : 129)   PROCEDURE == Initial cleaning
  0:00:19.602     7M / 4347M INFO    General                 (graph_simplification.hpp  : 674)   Flanking coverage based disconnection disabled
  0:00:19.611     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Self conjugate edge remover
  0:00:19.623     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Self conjugate edge remover triggered 1 times
  0:00:19.632     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Initial tip clipper
  0:00:19.648     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Initial tip clipper triggered 415 times
  0:00:19.654     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Initial ec remover
  0:00:19.803     8M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Initial ec remover triggered 6404 times
  0:00:19.816     8M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Initial isolated edge remover
  0:00:19.836     8M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Initial isolated edge remover triggered 254 times
  0:00:19.844     7M / 4347M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Simplification (id: simplification)
  0:00:19.853     7M / 4347M INFO    General                 (simplification.cpp        : 397)   Graph simplification started
  0:00:19.861     7M / 4347M INFO    General                 (graph_simplification.hpp  : 646)   Creating parallel br instance
  0:00:19.869     7M / 4347M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 1
  0:00:19.878     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:19.895     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 79 times
  0:00:19.903     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:19.972     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 1572 times
  0:00:19.985     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.002     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 6 times
  0:00:20.009     7M / 4347M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 2
  0:00:20.017     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.026     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 2 times
  0:00:20.034     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.042     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 1 times
  0:00:20.050     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.062     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 346 times
  0:00:20.067     7M / 4347M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 3
  0:00:20.071     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.080     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 12 times
  0:00:20.088     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.102     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 159 times
  0:00:20.106     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.116     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 240 times
  0:00:20.123     7M / 4347M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 4
  0:00:20.131     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.140     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 7 times
  0:00:20.148     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.161     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 135 times
  0:00:20.165     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.174     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 86 times
  0:00:20.182     7M / 4347M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 5
  0:00:20.190     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.198     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 3 times
  0:00:20.207     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.216     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 50 times
  0:00:20.223     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.232     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 66 times
  0:00:20.240     7M / 4347M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 6
  0:00:20.248     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.256     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 3 times
  0:00:20.265     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.274     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 30 times
  0:00:20.281     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.290     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 44 times
  0:00:20.298     7M / 4347M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 7
  0:00:20.306     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.315     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 0 times
  0:00:20.323     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.332     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 20 times
  0:00:20.339     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.348     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 41 times
  0:00:20.356     7M / 4347M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 8
  0:00:20.364     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.373     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 0 times
  0:00:20.381     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.390     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 24 times
  0:00:20.398     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.406     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 32 times
  0:00:20.414     7M / 4347M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 9
  0:00:20.423     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.431     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 3 times
  0:00:20.439     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.448     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 12 times
  0:00:20.456     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.464     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 22 times
  0:00:20.472     7M / 4347M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 10
  0:00:20.481     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.485     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 3 times
  0:00:20.493     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.502     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 9 times
  0:00:20.510     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.515     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 15 times
  0:00:20.522     7M / 4347M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 11
  0:00:20.531     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.547     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 1 times
  0:00:20.556     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.573     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 6 times
  0:00:20.581     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.597     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 0 times
  0:00:20.605     7M / 4347M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 12
  0:00:20.610     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.618     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 0 times
  0:00:20.626     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.635     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 0 times
  0:00:20.643     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.651     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 0 times
  0:00:20.668     7M / 4347M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Simplification Cleanup (id: simplification_cleanup)
  0:00:20.676     7M / 4347M INFO    General                 (simplification.cpp        : 189)   PROCEDURE == Post simplification
  0:00:20.684     7M / 4347M INFO    General                 (graph_simplification.hpp  : 455)   Disconnection of relatively low covered edges disabled
  0:00:20.693     7M / 4347M INFO    General                 (graph_simplification.hpp  : 494)   Complex tip clipping disabled
  0:00:20.701     7M / 4347M INFO    General                 (graph_simplification.hpp  : 646)   Creating parallel br instance
  0:00:20.709     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.726     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 0 times
  0:00:20.734     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.752     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 0 times
  0:00:20.756     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.772     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 0 times
  0:00:20.781     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.798     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 0 times
  0:00:20.806     7M / 4347M INFO    General                 (simplification.cpp        : 348)   Disrupting self-conjugate edges
  0:00:20.814     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Removing isolated edges
  0:00:20.831     7M / 4347M INFO   Simplification           (parallel_processing.hpp   : 171)   Removing isolated edges triggered 0 times
  0:00:20.839     7M / 4347M INFO    General                 (simplification.cpp        : 508)   After simplification:
  0:00:20.847     7M / 4347M INFO    General                 (simplification.cpp        : 509)     Average coverage = 49.3804
  0:00:20.856     7M / 4347M INFO    General                 (simplification.cpp        : 510)     Total length = 4997914
  0:00:20.864     7M / 4347M INFO    General                 (simplification.cpp        : 511)     Median edge length: 21564
  0:00:20.872     7M / 4347M INFO    General                 (simplification.cpp        : 512)     Edges: 3101
  0:00:20.881     7M / 4347M INFO    General                 (simplification.cpp        : 513)     Vertices: 2362
  0:00:20.889     7M / 4347M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Contig Output (id: contig_output)
  0:00:20.900     7M / 4347M INFO    General                 (read_converter.cpp        : 135)   Outputting contigs to "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K33/simplified_contigs"
  0:00:20.913    10M / 4347M INFO    General                 (binary_converter.cpp      : 143)   1551 reads written
  0:00:20.928     7M / 4347M INFO    General                 (pipeline.cpp              : 292)   SPAdes finished
  0:00:20.934     1M / 4347M INFO    General                 (main.cpp                  : 131)   Assembling time: 0 hours 0 minutes 20 seconds

===== K33 finished. 


===== K55 started. 


== Running: /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/spades/4.0.0/bin/spades-core /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K55/configs/config.info

  0:00:00.000     1M / 16M   INFO    General                 (main.cpp                  :  94)   Loaded config from "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K55/configs/config.info"
  0:00:00.034     1M / 16M   INFO    General                 (memory_limit.cpp          :  55)   Memory limit set to 187 Gb
  0:00:00.059     1M / 16M   INFO    General                 (main.cpp                  : 102)   Starting SPAdes, built from N/A, git revision N/A
  0:00:00.075     1M / 16M   INFO    General                 (main.cpp                  : 103)   Maximum k-mer length: 128
  0:00:00.088     1M / 16M   INFO    General                 (main.cpp                  : 104)   Assembling dataset ("/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/dataset.info") with K=55
  0:00:00.105     1M / 16M   INFO    General                 (main.cpp                  : 105)   Maximum # of threads to use (adjusted due to OMP capabilities): 16
  0:00:00.122     1M / 16M   INFO    General                 (pipeline.cpp              : 212)   SPAdes started
  0:00:00.138     1M / 16M   INFO    General                 (pipeline.cpp              : 225)   Starting from stage: read_conversion
  0:00:00.156     1M / 16M   INFO    General                 (pipeline.cpp              : 234)   Two-step repeat resolution disabled
  0:00:00.172     1M / 16M   INFO   GraphCore                (graph_core.hpp            : 689)   Graph created, vertex min_id: 3, edge min_id: 3
  0:00:00.189     1M / 16M   INFO   GraphCore                (graph_core.hpp            : 690)   Vertex size: 48, edge size: 40
  0:00:00.202     1M / 16M   INFO    General                 (edge_index.hpp            : 132)   Size of edge index entries: 12/8
  0:00:00.210     1M / 16M   INFO    General                 (pipeline.cpp              : 245)   Will need read mapping, kmer mapper will be attached
  0:00:00.219     1M / 16M   INFO   StageManager             (stage.cpp                 : 189)   STAGE == Binary Read Conversion (id: read_conversion)
  0:00:00.228     1M / 16M   INFO    General                 (read_converter.cpp        :  57)   Binary reads detected
  0:00:00.244     1M / 16M   INFO   StageManager             (stage.cpp                 : 189)   STAGE == de Bruijn graph construction (id: construction)
  0:00:00.248     1M / 16M   INFO    General                 (construction.cpp          : 115)   Contigs from previous K will be used: /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K33/simplified_contigs
  0:00:00.316     1M / 16M   INFO    General                 (construction.cpp          : 150)   Max read length 151
  0:00:00.328     1M / 16M   INFO    General                 (construction.cpp          : 156)   Average read length 146.498
  0:00:00.340     1M / 16M   INFO    General                 (stage.cpp                 : 121)   PROCEDURE == k+1-mer counting (id: construction:kpomer_counting)
  0:00:00.344     1M / 16M   INFO    General                 (kmer_index_builder.hpp    : 258)   Splitting kmer instances into 160 files using 16 threads. This might take a while.
  0:00:00.357     1M / 16M   INFO    General                 (file_limit.hpp            :  43)   Open file limit set to 51200
  0:00:00.373     1M / 16M   INFO    General                 (kmer_splitter.hpp         :  94)   Memory available for splitting buffers: 3.89582 Gb
  0:00:00.381     1M / 16M   INFO    General                 (kmer_splitter.hpp         : 102)   Using cell size of 209715
  0:00:03.523  9602M / 9602M INFO    General                 (kmer_splitters.hpp        : 128)   Processed 4363508 reads
  0:00:03.546     2M / 3525M INFO    General                 (kmer_splitters.hpp        : 134)   Used 4363508 reads
  0:00:03.689     2M / 3525M INFO    General                 (kmer_index_builder.hpp    : 264)   Starting k-mer counting.
  0:00:04.865     2M / 3525M INFO    General                 (kmer_index_builder.hpp    : 275)   K-mer counting done. There are 5726146 kmers in total.
  0:00:04.877     1M / 3525M INFO    General                 (stage.cpp                 : 121)   PROCEDURE == Extension index construction (id: construction:extension_index_construction)
  0:00:04.970     2M / 3525M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 453)   Building kmer index
  0:00:04.986     2M / 3525M INFO    General                 (kmer_index_builder.hpp    : 258)   Splitting kmer instances into 160 files using 16 threads. This might take a while.
  0:00:04.999     2M / 3525M INFO    General                 (file_limit.hpp            :  43)   Open file limit set to 51200
  0:00:05.016     2M / 3525M INFO    General                 (kmer_splitter.hpp         :  94)   Memory available for splitting buffers: 3.89581 Gb
  0:00:05.025     2M / 3525M INFO    General                 (kmer_splitter.hpp         : 102)   Using cell size of 209715
  0:00:06.908  9602M / 9602M INFO    General                 (kmer_splitters.hpp        : 197)   Processed 5726146 kmers
  0:00:06.918  9602M / 9602M INFO    General                 (kmer_splitters.hpp        : 202)   Used 5726146 kmers.
  0:00:06.929     2M / 3525M INFO    General                 (kmer_index_builder.hpp    : 264)   Starting k-mer counting.
  0:00:08.118     2M / 3525M INFO    General                 (kmer_index_builder.hpp    : 275)   K-mer counting done. There are 5721339 kmers in total.
  0:00:08.140     2M / 3525M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 410)   Building perfect hash indices
  0:00:08.279     7M / 3525M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 446)   Index built. Total 5721339 kmers, 4253880 bytes occupied (5.94809 bits per kmer).
  0:00:08.279     7M / 3525M INFO    General                 (kmer_index_builder.hpp    : 168)   Merging final buckets.
  0:00:09.365    12M / 3525M INFO   DeBruijnExtensionIndexBu (kmer_extension_index_build: 101)   Building k-mer extensions from k+1-mers
  0:00:09.585    12M / 3525M INFO   DeBruijnExtensionIndexBu (kmer_extension_index_build: 106)   Building k-mer extensions from k+1-mers finished.
  0:00:09.596    12M / 3525M INFO    General                 (stage.cpp                 : 121)   PROCEDURE == Condensing graph (id: construction:graph_condensing)
  0:00:09.694    12M / 3525M INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 381)   Extracting unbranching paths
  0:00:09.851    17M / 3525M INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 400)   Extracting unbranching paths finished. 92973 sequences extracted
  0:00:09.985    17M / 3525M INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 336)   Collecting perfect loops
  0:00:10.083    17M / 3525M INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 369)   Collecting perfect loops finished. 0 loops collected
  0:00:10.094    17M / 3525M INFO   DeBruijnGraphConstructor (debruijn_graph_constructor: 586)   Sorting edges...
  0:00:10.098    17M / 3525M INFO   DeBruijnGraphConstructor (debruijn_graph_constructor: 588)   Edges sorted
  0:00:10.098    17M / 3525M INFO    General                 (debruijn_graph_constructor: 516)   Total 185946 edges to create
  0:00:10.098    24M / 3525M INFO    General                 (debruijn_graph_constructor: 519)   Collecting link records
  0:00:10.263    27M / 3525M INFO    General                 (debruijn_graph_constructor: 521)   Ordering link records
  0:00:10.265    27M / 3525M INFO    General                 (debruijn_graph_constructor: 524)   Sorting done
  0:00:10.267    28M / 3525M INFO    General                 (debruijn_graph_constructor: 537)   Sorting LinkRecords...
  0:00:10.268    28M / 3525M INFO    General                 (debruijn_graph_constructor: 540)   LinkRecords sorted
  0:00:10.268    28M / 3525M INFO    General                 (debruijn_graph_constructor: 542)   Total 88166 vertices to create
  0:00:10.269    37M / 3525M INFO    General                 (debruijn_graph_constructor: 545)   Connecting the graph
  0:00:10.429    32M / 3525M INFO    General                 (stage.cpp                 : 121)   PROCEDURE == Filling coverage indices (PHM) (id: construction:coverage_filling_phm)
  0:00:10.430    32M / 3525M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 410)   Building perfect hash indices
  0:00:10.534    37M / 3525M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 446)   Index built. Total 5726146 kmers, 4257296 bytes occupied (5.94787 bits per kmer).
  0:00:10.548    59M / 3525M INFO    General                 (coverage_hash_map_builder.:  49)   Collecting k-mer coverage information from reads, this takes a while.
  0:00:13.463    59M / 3525M INFO    General                 (construction.cpp          : 427)   Filling coverage and flanking coverage from PHM
  0:00:13.595    59M / 3525M INFO    General                 (coverage_filling.hpp      :  83)   Processed 185943 edges
  0:00:14.032    21M / 3525M INFO   StageManager             (stage.cpp                 : 189)   STAGE == EC Threshold Finding (id: ec_threshold_finder)
  0:00:14.035    21M / 3525M INFO    General                 (kmer_coverage_model.cpp   : 182)   Kmer coverage valley at: 5
  0:00:14.035    21M / 3525M INFO    General                 (kmer_coverage_model.cpp   : 202)   K-mer histogram maximum: 32
  0:00:14.036    21M / 3525M INFO    General                 (kmer_coverage_model.cpp   : 238)   Estimated median coverage: 33. Coverage mad: 13.3434
  0:00:14.036    21M / 3525M INFO    General                 (kmer_coverage_model.cpp   : 260)   Fitting coverage model
  0:00:14.083    21M / 3525M INFO    General                 (kmer_coverage_model.cpp   : 296)   ... iteration 2
  0:00:14.202    21M / 3525M INFO    General                 (kmer_coverage_model.cpp   : 296)   ... iteration 4
  0:00:14.751    21M / 3525M INFO    General                 (kmer_coverage_model.cpp   : 296)   ... iteration 8
  0:00:15.786    21M / 3525M INFO    General                 (kmer_coverage_model.cpp   : 296)   ... iteration 16
  0:00:17.039    21M / 3525M INFO    General                 (kmer_coverage_model.cpp   : 310)   Fitted mean coverage: 33.0349. Fitted coverage std. dev: 12.2823
  0:00:17.049    21M / 3525M INFO    General                 (kmer_coverage_model.cpp   : 335)   Probability of erroneous kmer at valley: 0.476596
  0:00:17.058    21M / 3525M INFO    General                 (kmer_coverage_model.cpp   : 359)   Preliminary threshold calculated as: 15
  0:00:17.062    21M / 3525M INFO    General                 (kmer_coverage_model.cpp   : 363)   Threshold adjusted to: 15
  0:00:17.071    21M / 3525M INFO    General                 (kmer_coverage_model.cpp   : 376)   Estimated genome size (ignoring repeats): 4588115
  0:00:17.079    21M / 3525M INFO    General                 (genomic_info_filler.cpp   :  56)   Mean coverage was calculated as 33.0349
  0:00:17.088    21M / 3525M INFO    General                 (genomic_info_filler.cpp   :  71)   EC coverage threshold value was calculated as 15
  0:00:17.097    21M / 3525M INFO    General                 (genomic_info_filler.cpp   :  72)   Trusted kmer low bound: 0
  0:00:17.105    21M / 3525M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Gap Closer (id: early_gapcloser)
  0:00:17.119    34M / 3525M INFO    General                 (edge_index.hpp            : 132)   Size of edge index entries: 12/8
  0:00:17.137    36M / 3525M INFO    General                 (gap_closer.cpp            : 102)   Total edges in tip neighborhood: 41086 out of 185943, length: 615870
  0:00:17.148    36M / 3525M INFO    General                 (edge_index.hpp            : 196)   Using small index (max_id = 186878)
  0:00:17.169    36M / 3525M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 410)   Building perfect hash indices
  0:00:17.413    37M / 3525M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 446)   Index built. Total 578693 kmers, 418720 bytes occupied (5.78849 bits per kmer).
  0:00:17.424    41M / 3525M INFO    General                 (edge_index_builders.hpp   : 266)   Collecting edge information from graph, this takes a while.
  0:00:17.462    41M / 3525M INFO    General                 (sequence_mapper_notifier.h:  64)   Starting sequence mapping
  0:00:19.731    41M / 3525M INFO    General                 (sequence_mapper_notifier.h: 103)   Total 1089342 reads processed
  0:00:19.736    36M / 3525M INFO    General                 (gap_closer.cpp            : 491)   Initializing gap closer
  0:00:19.736    36M / 3525M INFO   GapCloser                (gap_closer.cpp            : 406)   Collecting gap candidates
  0:00:19.741    36M / 3525M INFO   GapCloser                (gap_closer.cpp            : 410)   Total 168 tips collected, total 182 connection candidates
  0:00:19.742    36M / 3525M INFO   GapCloser                (gap_closer.cpp            : 431)   Closing short gaps complete: filled 11 gaps after checking 170 candidates
  0:00:19.747    36M / 3525M INFO    General                 (gap_closer.cpp            : 495)   Gap closer done
  0:00:19.748    21M / 3525M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Raw Simplification (id: raw_simplification)
  0:00:19.748    21M / 3525M INFO    General                 (simplification.cpp        : 129)   PROCEDURE == Initial cleaning
  0:00:19.748    21M / 3525M INFO    General                 (graph_simplification.hpp  : 674)   Flanking coverage based disconnection disabled
  0:00:19.748    21M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Self conjugate edge remover
  0:00:19.752    21M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Self conjugate edge remover triggered 0 times
  0:00:19.752    21M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Initial tip clipper
  0:00:19.989    20M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Initial tip clipper triggered 37801 times
  0:00:19.989    20M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Initial ec remover
  0:00:20.115    20M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Initial ec remover triggered 3772 times
  0:00:20.115    20M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Initial isolated edge remover
  0:00:20.117    20M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Initial isolated edge remover triggered 1064 times
  0:00:20.117    19M / 3525M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Simplification (id: simplification)
  0:00:20.118    19M / 3525M INFO    General                 (simplification.cpp        : 397)   Graph simplification started
  0:00:20.118    19M / 3525M INFO    General                 (graph_simplification.hpp  : 646)   Creating parallel br instance
  0:00:20.118    19M / 3525M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 1
  0:00:20.118    19M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.120    19M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 138 times
  0:00:20.120    19M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.243    26M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 1117 times
  0:00:20.243    26M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.244    26M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 2 times
  0:00:20.245    26M / 3525M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 2
  0:00:20.245    26M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.245    26M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 2 times
  0:00:20.245    26M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.245    26M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 0 times
  0:00:20.245    26M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.248    26M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 224 times
  0:00:20.248    26M / 3525M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 3
  0:00:20.248    26M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.249    26M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 7 times
  0:00:20.249    26M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.265    26M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 146 times
  0:00:20.265    26M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.266    26M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 112 times
  0:00:20.266    26M / 3525M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 4
  0:00:20.266    26M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.267    26M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 9 times
  0:00:20.267    26M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.275    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 77 times
  0:00:20.275    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.275    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 70 times
  0:00:20.275    27M / 3525M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 5
  0:00:20.275    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.276    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 5 times
  0:00:20.276    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.281    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 53 times
  0:00:20.281    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.282    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 37 times
  0:00:20.282    27M / 3525M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 6
  0:00:20.282    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.282    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 1 times
  0:00:20.282    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.285    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 24 times
  0:00:20.285    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.285    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 27 times
  0:00:20.285    27M / 3525M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 7
  0:00:20.285    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.285    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 4 times
  0:00:20.285    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.286    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 10 times
  0:00:20.287    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.287    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 22 times
  0:00:20.287    27M / 3525M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 8
  0:00:20.287    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.287    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 4 times
  0:00:20.288    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.288    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 7 times
  0:00:20.288    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.289    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 18 times
  0:00:20.289    27M / 3525M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 9
  0:00:20.289    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.289    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 1 times
  0:00:20.289    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.290    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 10 times
  0:00:20.290    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.290    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 8 times
  0:00:20.290    27M / 3525M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 10
  0:00:20.291    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.291    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 0 times
  0:00:20.291    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.291    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 0 times
  0:00:20.291    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.291    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 9 times
  0:00:20.291    27M / 3525M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 11
  0:00:20.291    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.292    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 1 times
  0:00:20.292    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.293    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 0 times
  0:00:20.293    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.293    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 0 times
  0:00:20.293    27M / 3525M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 12
  0:00:20.293    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.294    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 0 times
  0:00:20.294    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.294    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 0 times
  0:00:20.294    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.294    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 0 times
  0:00:20.294    27M / 3525M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Gap Closer (id: late_gapcloser)
  0:00:20.298    40M / 3525M INFO    General                 (edge_index.hpp            : 132)   Size of edge index entries: 12/8
  0:00:20.298    40M / 3525M INFO    General                 (gap_closer.cpp            : 102)   Total edges in tip neighborhood: 265 out of 1768, length: 4113019
  0:00:20.299    40M / 3525M INFO    General                 (edge_index.hpp            : 196)   Using small index (max_id = 186878)
  0:00:20.316    40M / 3525M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 410)   Building perfect hash indices
  0:00:20.477    43M / 3525M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 446)   Index built. Total 3131990 kmers, 2262848 bytes occupied (5.77996 bits per kmer).
  0:00:20.494    67M / 3525M INFO    General                 (edge_index_builders.hpp   : 266)   Collecting edge information from graph, this takes a while.
  0:00:20.559    67M / 3525M INFO    General                 (sequence_mapper_notifier.h:  64)   Starting sequence mapping
  0:00:22.603    67M / 3525M INFO    General                 (sequence_mapper_notifier.h: 103)   Total 1089342 reads processed
  0:00:22.607    41M / 3525M INFO    General                 (gap_closer.cpp            : 491)   Initializing gap closer
  0:00:22.611    41M / 3525M INFO   GapCloser                (gap_closer.cpp            : 406)   Collecting gap candidates
  0:00:22.611    41M / 3525M INFO   GapCloser                (gap_closer.cpp            : 410)   Total 84 tips collected, total 91 connection candidates
  0:00:22.612    41M / 3525M INFO   GapCloser                (gap_closer.cpp            : 431)   Closing short gaps complete: filled 4 gaps after checking 86 candidates
  0:00:22.613    41M / 3525M INFO    General                 (gap_closer.cpp            : 495)   Gap closer done
  0:00:22.614    27M / 3525M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Simplification Cleanup (id: simplification_cleanup)
  0:00:22.614    27M / 3525M INFO    General                 (simplification.cpp        : 189)   PROCEDURE == Post simplification
  0:00:22.614    27M / 3525M INFO    General                 (graph_simplification.hpp  : 455)   Disconnection of relatively low covered edges disabled
  0:00:22.614    27M / 3525M INFO    General                 (graph_simplification.hpp  : 494)   Complex tip clipping disabled
  0:00:22.614    27M / 3525M INFO    General                 (graph_simplification.hpp  : 646)   Creating parallel br instance
  0:00:22.615    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:22.615    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 0 times
  0:00:22.615    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:22.616    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 0 times
  0:00:22.616    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:22.617    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 0 times
  0:00:22.617    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:22.618    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 0 times
  0:00:22.618    27M / 3525M INFO    General                 (simplification.cpp        : 348)   Disrupting self-conjugate edges
  0:00:22.618    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Removing isolated edges
  0:00:22.619    27M / 3525M INFO   Simplification           (parallel_processing.hpp   : 171)   Removing isolated edges triggered 0 times
  0:00:22.619    27M / 3525M INFO    General                 (simplification.cpp        : 508)   After simplification:
  0:00:22.619    27M / 3525M INFO    General                 (simplification.cpp        : 509)     Average coverage = 39.6118
  0:00:22.620    27M / 3525M INFO    General                 (simplification.cpp        : 510)     Total length = 5015454
  0:00:22.620    27M / 3525M INFO    General                 (simplification.cpp        : 511)     Median edge length: 43382
  0:00:22.620    27M / 3525M INFO    General                 (simplification.cpp        : 512)     Edges: 1760
  0:00:22.620    27M / 3525M INFO    General                 (simplification.cpp        : 513)     Vertices: 1478
  0:00:22.620    27M / 3525M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Contig Output (id: contig_output)
  0:00:22.624    27M / 3525M INFO    General                 (read_converter.cpp        : 135)   Outputting contigs to "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K55/simplified_contigs"
  0:00:22.633    31M / 3525M INFO    General                 (binary_converter.cpp      : 143)   880 reads written
  0:00:22.649    27M / 3525M INFO    General                 (pipeline.cpp              : 292)   SPAdes finished
  0:00:22.663     1M / 3525M INFO    General                 (main.cpp                  : 131)   Assembling time: 0 hours 0 minutes 22 seconds

===== K55 finished. 


===== K77 started. 


== Running: /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/spades/4.0.0/bin/spades-core /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/configs/config.info

  0:00:00.000     1M / 16M   INFO    General                 (main.cpp                  :  94)   Loaded config from "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/configs/config.info"
  0:00:00.022     1M / 16M   INFO    General                 (memory_limit.cpp          :  55)   Memory limit set to 187 Gb
  0:00:00.032     1M / 16M   INFO    General                 (main.cpp                  : 102)   Starting SPAdes, built from N/A, git revision N/A
  0:00:00.041     1M / 16M   INFO    General                 (main.cpp                  : 103)   Maximum k-mer length: 128
  0:00:00.051     1M / 16M   INFO    General                 (main.cpp                  : 104)   Assembling dataset ("/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/dataset.info") with K=77
  0:00:00.060     1M / 16M   INFO    General                 (main.cpp                  : 105)   Maximum # of threads to use (adjusted due to OMP capabilities): 16
  0:00:00.070     1M / 16M   INFO    General                 (pipeline.cpp              : 212)   SPAdes started
  0:00:00.077     1M / 16M   INFO    General                 (pipeline.cpp              : 225)   Starting from stage: read_conversion
  0:00:00.088     1M / 16M   INFO    General                 (pipeline.cpp              : 234)   Two-step repeat resolution disabled
  0:00:00.098     1M / 16M   INFO   GraphCore                (graph_core.hpp            : 689)   Graph created, vertex min_id: 3, edge min_id: 3
  0:00:00.107     1M / 16M   INFO   GraphCore                (graph_core.hpp            : 690)   Vertex size: 48, edge size: 40
  0:00:00.116     1M / 16M   INFO    General                 (edge_index.hpp            : 132)   Size of edge index entries: 12/8
  0:00:00.125     1M / 16M   INFO    General                 (pipeline.cpp              : 245)   Will need read mapping, kmer mapper will be attached
  0:00:00.135     1M / 16M   INFO   StageManager             (stage.cpp                 : 189)   STAGE == Binary Read Conversion (id: read_conversion)
  0:00:00.145     1M / 16M   INFO    General                 (read_converter.cpp        :  57)   Binary reads detected
  0:00:00.156     1M / 16M   INFO   StageManager             (stage.cpp                 : 189)   STAGE == de Bruijn graph construction (id: construction)
  0:00:00.170     1M / 16M   INFO    General                 (construction.cpp          : 115)   Contigs from previous K will be used: /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K55/simplified_contigs
  0:00:00.232     1M / 16M   INFO    General                 (construction.cpp          : 150)   Max read length 151
  0:00:00.250     1M / 16M   INFO    General                 (construction.cpp          : 156)   Average read length 146.498
  0:00:00.260     1M / 16M   INFO    General                 (stage.cpp                 : 121)   PROCEDURE == k+1-mer counting (id: construction:kpomer_counting)
  0:00:00.269     1M / 16M   INFO    General                 (kmer_index_builder.hpp    : 258)   Splitting kmer instances into 160 files using 16 threads. This might take a while.
  0:00:00.293     1M / 16M   INFO    General                 (file_limit.hpp            :  43)   Open file limit set to 51200
  0:00:00.314     1M / 16M   INFO    General                 (kmer_splitter.hpp         :  94)   Memory available for splitting buffers: 3.89582 Gb
  0:00:00.324     1M / 16M   INFO    General                 (kmer_splitter.hpp         : 102)   Using cell size of 139810
  0:00:03.681  9602M / 9602M INFO    General                 (kmer_splitters.hpp        : 128)   Processed 4362166 reads
  0:00:03.707     2M / 4078M INFO    General                 (kmer_splitters.hpp        : 134)   Used 4362166 reads
  0:00:03.890     2M / 4078M INFO    General                 (kmer_index_builder.hpp    : 264)   Starting k-mer counting.
  0:00:05.156     2M / 4078M INFO    General                 (kmer_index_builder.hpp    : 275)   K-mer counting done. There are 5778873 kmers in total.
  0:00:05.166     1M / 4078M INFO    General                 (stage.cpp                 : 121)   PROCEDURE == Extension index construction (id: construction:extension_index_construction)
  0:00:05.272     2M / 4078M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 453)   Building kmer index
  0:00:05.279     2M / 4078M INFO    General                 (kmer_index_builder.hpp    : 258)   Splitting kmer instances into 160 files using 16 threads. This might take a while.
  0:00:05.292     2M / 4078M INFO    General                 (file_limit.hpp            :  43)   Open file limit set to 51200
  0:00:05.306     2M / 4078M INFO    General                 (kmer_splitter.hpp         :  94)   Memory available for splitting buffers: 3.89581 Gb
  0:00:05.314     2M / 4078M INFO    General                 (kmer_splitter.hpp         : 102)   Using cell size of 139810
  0:00:07.391  9602M / 9602M INFO    General                 (kmer_splitters.hpp        : 197)   Processed 5778873 kmers
  0:00:07.405  9602M / 9602M INFO    General                 (kmer_splitters.hpp        : 202)   Used 5778873 kmers.
  0:00:07.416     2M / 4078M INFO    General                 (kmer_index_builder.hpp    : 264)   Starting k-mer counting.
  0:00:08.672     2M / 4078M INFO    General                 (kmer_index_builder.hpp    : 275)   K-mer counting done. There are 5778495 kmers in total.
  0:00:08.683     2M / 4078M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 410)   Building perfect hash indices
  0:00:08.828     7M / 4078M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 446)   Index built. Total 5778495 kmers, 4294808 bytes occupied (5.94592 bits per kmer).
  0:00:08.828     7M / 4078M INFO    General                 (kmer_index_builder.hpp    : 168)   Merging final buckets.
  0:00:10.114    12M / 4078M INFO   DeBruijnExtensionIndexBu (kmer_extension_index_build: 101)   Building k-mer extensions from k+1-mers
  0:00:10.327    12M / 4078M INFO   DeBruijnExtensionIndexBu (kmer_extension_index_build: 106)   Building k-mer extensions from k+1-mers finished.
  0:00:10.337    12M / 4078M INFO    General                 (stage.cpp                 : 121)   PROCEDURE == Condensing graph (id: construction:graph_condensing)
  0:00:10.407    12M / 4078M INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 381)   Extracting unbranching paths
  0:00:10.570    18M / 4078M INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 400)   Extracting unbranching paths finished. 85131 sequences extracted
  0:00:10.711    17M / 4078M INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 336)   Collecting perfect loops
  0:00:10.810    18M / 4078M INFO   UnbranchingPathExtractor (debruijn_graph_constructor: 369)   Collecting perfect loops finished. 0 loops collected
  0:00:10.824    17M / 4078M INFO   DeBruijnGraphConstructor (debruijn_graph_constructor: 586)   Sorting edges...
  0:00:10.827    17M / 4078M INFO   DeBruijnGraphConstructor (debruijn_graph_constructor: 588)   Edges sorted
  0:00:10.827    17M / 4078M INFO    General                 (debruijn_graph_constructor: 516)   Total 170262 edges to create
  0:00:10.827    24M / 4078M INFO    General                 (debruijn_graph_constructor: 519)   Collecting link records
  0:00:10.991    27M / 4078M INFO    General                 (debruijn_graph_constructor: 521)   Ordering link records
  0:00:10.993    27M / 4078M INFO    General                 (debruijn_graph_constructor: 524)   Sorting done
  0:00:10.995    28M / 4078M INFO    General                 (debruijn_graph_constructor: 537)   Sorting LinkRecords...
  0:00:10.996    28M / 4078M INFO    General                 (debruijn_graph_constructor: 540)   LinkRecords sorted
  0:00:10.996    28M / 4078M INFO    General                 (debruijn_graph_constructor: 542)   Total 84753 vertices to create
  0:00:10.996    36M / 4078M INFO    General                 (debruijn_graph_constructor: 545)   Connecting the graph
  0:00:11.157    31M / 4078M INFO    General                 (stage.cpp                 : 121)   PROCEDURE == Filling coverage indices (PHM) (id: construction:coverage_filling_phm)
  0:00:11.157    31M / 4078M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 410)   Building perfect hash indices
  0:00:11.280    36M / 4078M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 446)   Index built. Total 5778873 kmers, 4295048 bytes occupied (5.94586 bits per kmer).
  0:00:11.294    59M / 4078M INFO    General                 (coverage_hash_map_builder.:  49)   Collecting k-mer coverage information from reads, this takes a while.
  0:00:13.667    59M / 4078M INFO    General                 (construction.cpp          : 427)   Filling coverage and flanking coverage from PHM
  0:00:13.799    59M / 4078M INFO    General                 (coverage_filling.hpp      :  83)   Processed 170261 edges
  0:00:14.281    20M / 4078M INFO   StageManager             (stage.cpp                 : 189)   STAGE == EC Threshold Finding (id: ec_threshold_finder)
  0:00:14.282    20M / 4078M INFO    General                 (kmer_coverage_model.cpp   : 182)   Kmer coverage valley at: 2
  0:00:14.282    20M / 4078M INFO    General                 (kmer_coverage_model.cpp   : 202)   K-mer histogram maximum: 24
  0:00:14.282    20M / 4078M INFO    General                 (kmer_coverage_model.cpp   : 238)   Estimated median coverage: 24. Coverage mad: 11.8608
  0:00:14.282    20M / 4078M INFO    General                 (kmer_coverage_model.cpp   : 260)   Fitting coverage model
  0:00:14.312    20M / 4078M INFO    General                 (kmer_coverage_model.cpp   : 296)   ... iteration 2
  0:00:14.395    20M / 4078M INFO    General                 (kmer_coverage_model.cpp   : 296)   ... iteration 4
  0:00:14.766    20M / 4078M INFO    General                 (kmer_coverage_model.cpp   : 296)   ... iteration 8
  0:00:15.589    20M / 4078M INFO    General                 (kmer_coverage_model.cpp   : 296)   ... iteration 16
  0:00:16.720    20M / 4078M INFO    General                 (kmer_coverage_model.cpp   : 310)   Fitted mean coverage: 25.27. Fitted coverage std. dev: 10.0533
  0:00:16.729    20M / 4078M INFO    General                 (kmer_coverage_model.cpp   : 335)   Probability of erroneous kmer at valley: 0.792056
  0:00:16.737    20M / 4078M INFO    General                 (kmer_coverage_model.cpp   : 359)   Preliminary threshold calculated as: 12
  0:00:16.746    20M / 4078M INFO    General                 (kmer_coverage_model.cpp   : 363)   Threshold adjusted to: 12
  0:00:16.754    20M / 4078M INFO    General                 (kmer_coverage_model.cpp   : 376)   Estimated genome size (ignoring repeats): 4512326
  0:00:16.758    20M / 4078M INFO    General                 (genomic_info_filler.cpp   :  56)   Mean coverage was calculated as 25.27
  0:00:16.766    20M / 4078M INFO    General                 (genomic_info_filler.cpp   :  71)   EC coverage threshold value was calculated as 12
  0:00:16.775    20M / 4078M INFO    General                 (genomic_info_filler.cpp   :  72)   Trusted kmer low bound: 0
  0:00:16.783    20M / 4078M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Gap Closer (id: early_gapcloser)
  0:00:16.797    33M / 4078M INFO    General                 (edge_index.hpp            : 132)   Size of edge index entries: 12/8
  0:00:16.816    34M / 4078M INFO    General                 (gap_closer.cpp            : 102)   Total edges in tip neighborhood: 43707 out of 170261, length: 809360
  0:00:16.830    35M / 4078M INFO    General                 (edge_index.hpp            : 196)   Using small index (max_id = 171116)
  0:00:16.854    35M / 4078M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 410)   Building perfect hash indices
  0:00:17.110    35M / 4078M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 446)   Index built. Total 742600 kmers, 537096 bytes occupied (5.78611 bits per kmer).
  0:00:17.121    41M / 4078M INFO    General                 (edge_index_builders.hpp   : 266)   Collecting edge information from graph, this takes a while.
  0:00:17.149    41M / 4078M INFO    General                 (sequence_mapper_notifier.h:  64)   Starting sequence mapping
  0:00:19.503    41M / 4078M INFO    General                 (sequence_mapper_notifier.h: 103)   Total 1089342 reads processed
  0:00:19.520    35M / 4078M INFO    General                 (gap_closer.cpp            : 491)   Initializing gap closer
  0:00:19.524    35M / 4078M INFO   GapCloser                (gap_closer.cpp            : 406)   Collecting gap candidates
  0:00:19.543    35M / 4078M INFO   GapCloser                (gap_closer.cpp            : 410)   Total 162 tips collected, total 176 connection candidates
  0:00:19.552    35M / 4078M INFO   GapCloser                (gap_closer.cpp            : 431)   Closing short gaps complete: filled 0 gaps after checking 176 candidates
  0:00:19.568    35M / 4078M INFO    General                 (gap_closer.cpp            : 495)   Gap closer done
  0:00:19.577    20M / 4078M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Raw Simplification (id: raw_simplification)
  0:00:19.584    20M / 4078M INFO    General                 (simplification.cpp        : 129)   PROCEDURE == Initial cleaning
  0:00:19.593    20M / 4078M INFO    General                 (graph_simplification.hpp  : 674)   Flanking coverage based disconnection disabled
  0:00:19.597    20M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Self conjugate edge remover
  0:00:19.615    20M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Self conjugate edge remover triggered 0 times
  0:00:19.624    20M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Initial tip clipper
  0:00:19.939    18M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Initial tip clipper triggered 38912 times
  0:00:19.951    18M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Initial ec remover
  0:00:19.986    18M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Initial ec remover triggered 871 times
  0:00:19.990    18M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Initial isolated edge remover
  0:00:20.009    18M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Initial isolated edge remover triggered 1761 times
  0:00:20.016    17M / 4078M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Simplification (id: simplification)
  0:00:20.024    17M / 4078M INFO    General                 (simplification.cpp        : 397)   Graph simplification started
  0:00:20.033    17M / 4078M INFO    General                 (graph_simplification.hpp  : 646)   Creating parallel br instance
  0:00:20.041    17M / 4078M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 1
  0:00:20.050    17M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.068    17M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 189 times
  0:00:20.076    17M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.212    27M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 849 times
  0:00:20.224    27M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.244    27M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 0 times
  0:00:20.253    27M / 4078M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 2
  0:00:20.261    27M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.270    27M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 7 times
  0:00:20.278    27M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.287    27M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 1 times
  0:00:20.296    27M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.305    27M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 30 times
  0:00:20.312    27M / 4078M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 3
  0:00:20.321    27M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.330    27M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 5 times
  0:00:20.338    27M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.351    28M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 23 times
  0:00:20.351    28M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.353    28M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 59 times
  0:00:20.353    28M / 4078M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 4
  0:00:20.353    28M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.353    28M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 6 times
  0:00:20.353    28M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.361    28M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 45 times
  0:00:20.361    28M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.361    28M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 42 times
  0:00:20.361    28M / 4078M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 5
  0:00:20.361    28M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.361    28M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 5 times
  0:00:20.361    28M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.366    28M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 27 times
  0:00:20.366    28M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.367    28M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 29 times
  0:00:20.367    28M / 4078M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 6
  0:00:20.367    28M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.367    28M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 2 times
  0:00:20.367    28M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.369    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 12 times
  0:00:20.369    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.370    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 22 times
  0:00:20.370    29M / 4078M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 7
  0:00:20.370    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.370    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 2 times
  0:00:20.370    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.371    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 7 times
  0:00:20.371    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.371    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 10 times
  0:00:20.371    29M / 4078M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 8
  0:00:20.371    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.371    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 0 times
  0:00:20.371    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.371    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 2 times
  0:00:20.371    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.372    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 6 times
  0:00:20.372    29M / 4078M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 9
  0:00:20.372    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.372    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 0 times
  0:00:20.372    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.372    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 4 times
  0:00:20.372    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.373    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 14 times
  0:00:20.373    29M / 4078M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 10
  0:00:20.373    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.373    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 0 times
  0:00:20.373    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.373    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 0 times
  0:00:20.373    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.374    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 9 times
  0:00:20.374    29M / 4078M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 11
  0:00:20.374    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.374    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 1 times
  0:00:20.374    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.375    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 1 times
  0:00:20.375    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.375    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 0 times
  0:00:20.375    29M / 4078M INFO    General                 (simplification.cpp        : 402)   PROCEDURE == Simplification cycle, iteration 12
  0:00:20.375    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:20.375    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 0 times
  0:00:20.375    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:20.375    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 0 times
  0:00:20.375    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Low coverage edge remover
  0:00:20.375    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Low coverage edge remover triggered 0 times
  0:00:20.376    29M / 4078M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Gap Closer (id: late_gapcloser)
  0:00:20.382    42M / 4078M INFO    General                 (edge_index.hpp            : 132)   Size of edge index entries: 12/8
  0:00:20.382    42M / 4078M INFO    General                 (gap_closer.cpp            : 102)   Total edges in tip neighborhood: 243 out of 1352, length: 4639905
  0:00:20.382    42M / 4078M INFO    General                 (edge_index.hpp            : 196)   Using small index (max_id = 171116)
  0:00:20.407    42M / 4078M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 410)   Building perfect hash indices
  0:00:20.624    44M / 4078M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 446)   Index built. Total 3512423 kmers, 2537632 bytes occupied (5.77979 bits per kmer).
  0:00:20.641    71M / 4078M INFO    General                 (edge_index_builders.hpp   : 266)   Collecting edge information from graph, this takes a while.
  0:00:20.724    71M / 4078M INFO    General                 (sequence_mapper_notifier.h:  64)   Starting sequence mapping
  0:00:22.225    72M / 4078M INFO    General                 (sequence_mapper_notifier.h: 103)   Total 1089342 reads processed
  0:00:22.232    42M / 4078M INFO    General                 (gap_closer.cpp            : 491)   Initializing gap closer
  0:00:22.236    42M / 4078M INFO   GapCloser                (gap_closer.cpp            : 406)   Collecting gap candidates
  0:00:22.237    42M / 4078M INFO   GapCloser                (gap_closer.cpp            : 410)   Total 79 tips collected, total 81 connection candidates
  0:00:22.238    42M / 4078M INFO   GapCloser                (gap_closer.cpp            : 431)   Closing short gaps complete: filled 0 gaps after checking 81 candidates
  0:00:22.238    42M / 4078M INFO    General                 (gap_closer.cpp            : 495)   Gap closer done
  0:00:22.240    29M / 4078M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Simplification Cleanup (id: simplification_cleanup)
  0:00:22.240    29M / 4078M INFO    General                 (simplification.cpp        : 189)   PROCEDURE == Post simplification
  0:00:22.240    29M / 4078M INFO    General                 (graph_simplification.hpp  : 455)   Disconnection of relatively low covered edges disabled
  0:00:22.240    29M / 4078M INFO    General                 (graph_simplification.hpp  : 494)   Complex tip clipping disabled
  0:00:22.240    29M / 4078M INFO    General                 (graph_simplification.hpp  : 646)   Creating parallel br instance
  0:00:22.240    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:22.240    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 0 times
  0:00:22.240    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:22.241    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 0 times
  0:00:22.241    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Tip clipper
  0:00:22.242    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Tip clipper triggered 0 times
  0:00:22.242    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Bulge remover
  0:00:22.242    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Bulge remover triggered 0 times
  0:00:22.243    29M / 4078M INFO    General                 (simplification.cpp        : 348)   Disrupting self-conjugate edges
  0:00:22.243    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 168)   Running Removing isolated edges
  0:00:22.243    29M / 4078M INFO   Simplification           (parallel_processing.hpp   : 171)   Removing isolated edges triggered 0 times
  0:00:22.244    29M / 4078M INFO    General                 (simplification.cpp        : 508)   After simplification:
  0:00:22.244    29M / 4078M INFO    General                 (simplification.cpp        : 509)     Average coverage = 30.0773
  0:00:22.244    29M / 4078M INFO    General                 (simplification.cpp        : 510)     Total length = 5021745
  0:00:22.245    29M / 4078M INFO    General                 (simplification.cpp        : 511)     Median edge length: 45206
  0:00:22.245    29M / 4078M INFO    General                 (simplification.cpp        : 512)     Edges: 1352
  0:00:22.245    29M / 4078M INFO    General                 (simplification.cpp        : 513)     Vertices: 1200
  0:00:22.245    29M / 4078M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Mismatch Correction (id: mismatch_correction)
  0:00:22.245    29M / 4078M INFO    General                 (graph_pack_helpers.cpp    :  44)   Index refill
  0:00:22.245    29M / 4078M INFO    General                 (edge_index.hpp            : 175)   Using small index (max_id = 171116)
  0:00:22.276    29M / 4078M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 410)   Building perfect hash indices
  0:00:22.577    32M / 4078M INFO   K-mer Index Building     (kmer_index_builder.hpp    : 446)   Index built. Total 5021745 kmers, 3627760 bytes occupied (5.77928 bits per kmer).
  0:00:22.600    72M / 4078M INFO    General                 (edge_index_builders.hpp   : 253)   Collecting edge information from graph, this takes a while.
  0:00:22.891    72M / 4078M INFO    General                 (graph_pack_helpers.cpp    :  54)   Normalizing k-mer map. Total 161036 kmers to process
  0:00:22.920    72M / 4078M INFO    General                 (graph_pack_helpers.cpp    :  56)   Normalizing done
  0:00:22.920    72M / 4078M INFO    General                 (mismatch_correction.cpp   : 392)   Collect potential mismatches
  0:00:22.940    72M / 4078M INFO    General                 (mismatch_correction.cpp   : 193)   Total 298 edges (out of 1352) with 4999 potential mismatch positions (16.7752 positions per edge)
  0:00:22.940    72M / 4078M INFO    General                 (mismatch_correction.cpp   : 394)   Potential mismatches collected
  0:00:22.953    73M / 4078M INFO    General                 (sequence_mapper_notifier.h:  64)   Starting sequence mapping
  0:00:23.558    73M / 4078M INFO    General                 (sequence_mapper_notifier.h:  85)   Processed 200000 reads
  0:00:23.560    73M / 4078M INFO    General                 (sequence_mapper_notifier.h:  85)   Processed 400000 reads
  0:00:23.562    74M / 4078M INFO    General                 (sequence_mapper_notifier.h:  85)   Processed 600000 reads
  0:00:23.563    74M / 4078M INFO    General                 (sequence_mapper_notifier.h:  85)   Processed 800000 reads
  0:00:23.636    74M / 4078M INFO    General                 (sequence_mapper_notifier.h:  85)   Processed 1000000 reads
  0:00:23.647    75M / 4078M INFO    General                 (sequence_mapper_notifier.h:  85)   Processed 1200000 reads
  0:00:23.660    75M / 4078M INFO    General                 (sequence_mapper_notifier.h:  85)   Processed 2200000 reads
  0:00:23.925    79M / 4078M INFO    General                 (sequence_mapper_notifier.h: 103)   Total 4360406 reads processed
  0:00:23.938    73M / 4078M INFO    General                 (mismatch_correction.cpp   : 387)   All edges processed
  0:00:23.939    72M / 4078M INFO    General                 (mismatch_correction.cpp   : 442)   Corrected 1 nucleotides
  0:00:23.939    72M / 4078M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Contig Output (id: contig_output)
  0:00:23.939    72M / 4078M INFO    General                 (contig_output.hpp         :  20)   Outputting contigs to "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/before_rr.fasta"
  0:00:24.006    72M / 4078M INFO    General                 (contig_output_stage.cpp   : 155)   Writing GFA graph to "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/assembly_graph_after_simplification.gfa"
  0:00:24.053    72M / 4078M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Paired Information Counting (id: late_pair_info_count)
  0:00:24.055    73M / 4078M INFO    General                 (graph_pack_helpers.cpp    :  54)   Normalizing k-mer map. Total 161036 kmers to process
  0:00:24.082    73M / 4078M INFO    General                 (graph_pack_helpers.cpp    :  56)   Normalizing done
  0:00:24.082    73M / 4078M INFO    General                 (pair_info_count.cpp       : 157)   Min edge length for estimation: 45206
  0:00:24.083    73M / 4078M INFO    General                 (pair_info_count.cpp       : 168)   Estimating insert size for library #0
  0:00:24.083    73M / 4078M INFO    General                 (pair_info_count.cpp       :  41)   Selecting usual mapper
  0:00:24.083    73M / 4078M INFO    General                 (paired_info_utils.cpp     :  87)   Estimating insert size (takes a while)
  0:00:24.237   345M / 4078M INFO    General                 (sequence_mapper_notifier.h:  64)   Starting sequence mapping
  0:00:25.073   345M / 4078M INFO    General                 (sequence_mapper_notifier.h: 103)   Total 1089342 reads processed
  0:00:25.146   345M / 4078M INFO    General                 (paired_info_utils.cpp     : 104)   Edge pairs: 6660
  0:00:25.156   345M / 4078M INFO    General                 (paired_info_utils.cpp     : 106)   445187 paired reads (40.8675% of all) aligned to long edges
  0:00:25.165    73M / 4078M INFO    General                 (pair_info_count.cpp       : 188)     Insert size = 563.483, deviation = 96.5037, left quantile = 464, right quantile = 685, read length = 151
  0:00:25.177    73M / 4078M INFO    General                 (pair_info_count.cpp       : 203)   Filtering data for library #0
  0:00:25.186    73M / 4078M INFO    General                 (pair_info_count.cpp       :  41)   Selecting usual mapper
  0:00:25.208    73M / 4078M INFO    General                 (sequence_mapper_notifier.h:  64)   Starting sequence mapping
  0:00:25.719    73M / 4078M INFO    General                 (sequence_mapper_notifier.h: 103)   Total 1089342 reads processed
  0:00:25.728    73M / 4078M INFO    General                 (pair_info_count.cpp       : 207)   Mapping library #0
  0:00:25.736    73M / 4078M INFO    General                 (pair_info_count.cpp       : 209)   Mapping paired reads (takes a while)
  0:00:25.745    73M / 4078M INFO    General                 (pair_info_count.cpp       :  41)   Selecting usual mapper
  0:00:25.753    73M / 4078M INFO    General                 (paired_info_utils.cpp     : 138)   Left insert size quantile 464, right insert size quantile 685, filtering threshold 2, rounding threshold 0
  0:00:25.775    86M / 4078M INFO    General                 (sequence_mapper_notifier.h:  64)   Starting sequence mapping
  0:00:26.557    88M / 4078M INFO    General                 (sequence_mapper_notifier.h: 103)   Total 1089342 reads processed
  0:00:26.572    75M / 4078M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Distance Estimation (id: distance_estimation)
  0:00:26.572    75M / 4078M INFO    General                 (distance_estimation.cpp   :  42)   Processing library #0
  0:00:26.572    75M / 4078M INFO    General                 (distance_estimation_utils.: 133)   Weight Filter Done
  0:00:26.573    75M / 4078M INFO   DistanceEstimator        (distance_estimation.hpp   : 116)   Using SIMPLE distance estimator
  0:00:26.576    75M / 4078M INFO    General                 (distance_estimation_utils.:  28)   Filtering info
  0:00:26.577    75M / 4078M INFO    General                 (pair_info_filters.hpp     : 243)   Start filtering; library index size: 6678
  0:00:26.578    75M / 4078M INFO    General                 (pair_info_filters.hpp     : 264)   Done filtering; library index size: 5126
  0:00:26.578    75M / 4078M INFO    General                 (distance_estimation_utils.: 139)   Refining clustered pair information
  0:00:26.578    75M / 4078M INFO    General                 (distance_estimation_utils.: 141)   The refining of clustered pair information has been finished
  0:00:26.578    75M / 4078M INFO    General                 (distance_estimation_utils.: 143)   Improving paired information
  0:00:26.586    75M / 4078M INFO   PairInfoImprover         (pair_info_improver.hpp    :  65)   Paired info stats: missing = 699; contradictional = 94
  0:00:26.592    75M / 4078M INFO   PairInfoImprover         (pair_info_improver.hpp    :  65)   Paired info stats: missing = 27; contradictional = 32
  0:00:26.592    75M / 4078M INFO    General                 (distance_estimation_utils.:  86)   Filling scaffolding index
  0:00:26.593    75M / 4078M INFO   DistanceEstimator        (distance_estimation.hpp   : 116)   Using SMOOTHING distance estimator
  0:00:26.979    75M / 4078M INFO    General                 (distance_estimation_utils.:  28)   Filtering info
  0:00:26.979    75M / 4078M INFO    General                 (pair_info_filters.hpp     : 243)   Start filtering; library index size: 1994
  0:00:26.980    75M / 4078M INFO    General                 (pair_info_filters.hpp     : 264)   Done filtering; library index size: 1994
  0:00:26.980    75M / 4078M INFO    General                 (distance_estimation.cpp   :  51)   Clearing raw paired index
  0:00:26.981    73M / 4078M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Repeat Resolving (id: repeat_resolving)
  0:00:26.982    73M / 4078M INFO    General                 (repeat_resolving.cpp      :  88)   Using Path-Extend repeat resolving
  0:00:26.982    73M / 4078M INFO    General                 (launcher.cpp              : 603)   ExSPAnder repeat resolving tool started
  0:00:26.986    73M / 4078M INFO    General                 (coverage_uniformity_analyz:  25)   genomic coverage is 25.7252 calculated of length 4948592
  0:00:26.989    76M / 4078M INFO    General                 (launcher.cpp              : 420)   Creating main extenders, unique edge length = 2000
  0:00:26.989    76M / 4078M INFO    General                 (extenders_logic.cpp       : 341)   Estimated coverage of library #0 is 30.0773
  0:00:26.989    76M / 4078M INFO    General                 (extenders_logic.cpp       : 352)   Creating extender; library index size: 6452
  0:00:26.990    76M / 4078M INFO    General                 (extenders_logic.cpp       : 341)   Estimated coverage of library #0 is 30.0773
  0:00:26.990    76M / 4078M INFO    General                 (extenders_logic.cpp       : 352)   Creating extender; library index size: 6452
  0:00:26.991    76M / 4078M INFO    General                 (extenders_logic.cpp       : 550)   Using 1 paired-end library
  0:00:26.991    76M / 4078M INFO    General                 (extenders_logic.cpp       : 551)   Using 1 paired-end scaffolding library
  0:00:26.991    76M / 4078M INFO    General                 (extenders_logic.cpp       : 552)   Using 0 single read libraries
  0:00:26.991    76M / 4078M INFO    General                 (launcher.cpp              : 449)   Total number of extenders is 3
  0:00:26.991    76M / 4078M INFO    General                 (path_extenders.cpp        :  36)   Processed 0 paths from 659 (0%)
  0:00:27.003    76M / 4078M INFO    General                 (path_extenders.cpp        :  36)   Processed 66 paths from 659 (10%)
  0:00:27.010    77M / 4078M INFO    General                 (path_extenders.cpp        :  34)   Processed 128 paths from 659 (19%)
  0:00:27.010    77M / 4078M INFO    General                 (path_extenders.cpp        :  36)   Processed 132 paths from 659 (20%)
  0:00:27.018    78M / 4078M INFO    General                 (path_extenders.cpp        :  36)   Processed 198 paths from 659 (30%)
  0:00:27.022    78M / 4078M INFO    General                 (path_extenders.cpp        :  34)   Processed 256 paths from 659 (38%)
  0:00:27.023    78M / 4078M INFO    General                 (path_extenders.cpp        :  36)   Processed 264 paths from 659 (40%)
  0:00:27.025    78M / 4078M INFO    General                 (path_extenders.cpp        :  36)   Processed 330 paths from 659 (50%)
  0:00:27.026    78M / 4078M INFO    General                 (path_extenders.cpp        :  36)   Processed 396 paths from 659 (60%)
  0:00:27.027    78M / 4078M INFO    General                 (path_extenders.cpp        :  36)   Processed 462 paths from 659 (70%)
  0:00:27.028    78M / 4078M INFO    General                 (path_extenders.cpp        :  34)   Processed 512 paths from 659 (77%)
  0:00:27.029    78M / 4078M INFO    General                 (path_extenders.cpp        :  36)   Processed 528 paths from 659 (80%)
  0:00:27.029    78M / 4078M INFO    General                 (path_extenders.cpp        :  36)   Processed 594 paths from 659 (90%)
  0:00:27.030    76M / 4078M INFO    General                 (launcher.cpp              : 252)   Finalizing paths
  0:00:27.030    76M / 4078M INFO    General                 (launcher.cpp              : 254)   Deduplicating paths
  0:00:27.033    76M / 4078M INFO    General                 (launcher.cpp              : 258)   Paths deduplicated
  0:00:27.033    76M / 4078M INFO   PEResolver               (pe_resolver.cpp           :  60)   Removing overlaps
  0:00:27.033    76M / 4078M INFO   PEResolver               (pe_resolver.cpp           :  63)   Sorting paths
  0:00:27.033    76M / 4078M INFO   PEResolver               (pe_resolver.cpp           :  70)   Marking overlaps
  0:00:27.033    76M / 4078M INFO   OverlapRemover           (overlap_remover.hpp       : 117)   Marking start/end overlaps
  0:00:27.038    76M / 4078M INFO   OverlapRemover           (overlap_remover.hpp       : 120)   Marking remaining overlaps
  0:00:27.043    76M / 4078M INFO   PEResolver               (pe_resolver.cpp           :  73)   Splitting paths
  0:00:27.044    77M / 4078M INFO   PEResolver               (pe_resolver.cpp           :  78)   Deduplicating paths
  0:00:27.046    76M / 4078M INFO   PEResolver               (pe_resolver.cpp           :  80)   Overlaps removed
  0:00:27.046    76M / 4078M INFO    General                 (launcher.cpp              : 275)   Paths finalized
  0:00:27.046    76M / 4078M INFO    General                 (launcher.cpp              : 456)   Closing gaps in paths
  0:00:27.047    77M / 4078M INFO    General                 (launcher.cpp              : 486)   Gap closing completed
  0:00:27.048    74M / 4078M INFO    General                 (launcher.cpp              : 304)   Traversing tandem repeats
  0:00:27.053    74M / 4078M INFO    General                 (launcher.cpp              : 314)   Traversed 0 loops
  0:00:27.054    74M / 4078M INFO    General                 (launcher.cpp              : 252)   Finalizing paths
  0:00:27.054    74M / 4078M INFO    General                 (launcher.cpp              : 254)   Deduplicating paths
  0:00:27.054    74M / 4078M INFO    General                 (launcher.cpp              : 258)   Paths deduplicated
  0:00:27.054    74M / 4078M INFO   PEResolver               (pe_resolver.cpp           :  60)   Removing overlaps
  0:00:27.054    74M / 4078M INFO   PEResolver               (pe_resolver.cpp           :  63)   Sorting paths
  0:00:27.054    74M / 4078M INFO   PEResolver               (pe_resolver.cpp           :  70)   Marking overlaps
  0:00:27.054    74M / 4078M INFO   OverlapRemover           (overlap_remover.hpp       : 117)   Marking start/end overlaps
  0:00:27.054    74M / 4078M INFO   OverlapRemover           (overlap_remover.hpp       : 120)   Marking remaining overlaps
  0:00:27.055    74M / 4078M INFO   PEResolver               (pe_resolver.cpp           :  73)   Splitting paths
  0:00:27.055    74M / 4078M INFO   PEResolver               (pe_resolver.cpp           :  78)   Deduplicating paths
  0:00:27.055    74M / 4078M INFO   PEResolver               (pe_resolver.cpp           :  80)   Overlaps removed
  0:00:27.055    74M / 4078M INFO    General                 (launcher.cpp              : 275)   Paths finalized
  0:00:27.056    74M / 4078M INFO    General                 (launcher.cpp              : 666)   ExSPAnder repeat resolving tool finished
  0:00:27.056    74M / 4078M INFO   StageManager             (stage.cpp                 : 189)   STAGE == Contig Output (id: contig_output)
  0:00:27.056    74M / 4078M INFO    General                 (contig_output.hpp         :  20)   Outputting contigs to "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/before_rr.fasta"
  0:00:27.086    74M / 4078M INFO    General                 (contig_output_stage.cpp   : 155)   Writing GFA graph to "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/assembly_graph_with_scaffolds.gfa"
  0:00:27.103    74M / 4078M INFO    General                 (contig_output_stage.cpp   : 169)   Outputting FastG graph to "/scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/assembly_graph.fastg"
  0:00:27.174    74M / 4078M INFO    General                 (contig_output_stage.cpp   : 200)   Breaking scaffolds
  0:00:27.186    81M / 4078M INFO    General                 (contig_output_stage.cpp   : 101)   Outputting contigs to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/final_contigs.fasta
  0:00:27.217    81M / 4078M INFO    General                 (contig_output_stage.cpp   : 107)   Outputting FastG paths to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/final_contigs.paths
  0:00:27.246    79M / 4078M INFO    General                 (contig_output_stage.cpp   : 101)   Outputting contigs to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/scaffolds.fasta
  0:00:27.276    79M / 4078M INFO    General                 (contig_output_stage.cpp   : 107)   Outputting FastG paths to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/scaffolds.paths
  0:00:27.283    79M / 4078M INFO    General                 (contig_output_stage.cpp   : 114)   Populating GFA with scaffold paths
  0:00:27.295    74M / 4078M INFO    General                 (pipeline.cpp              : 292)   SPAdes finished
  0:00:27.318     1M / 4078M INFO    General                 (main.cpp                  : 131)   Assembling time: 0 hours 0 minutes 27 seconds

===== K77 finished. 


===== Copy files started. 


== Running: /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/python/3.11.5/bin/python3 /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/spades/4.0.0/share/spades/spades_pipeline/scripts/copy_files.py /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/before_rr.fasta /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/before_rr.fasta /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/assembly_graph_after_simplification.gfa /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/assembly_graph_after_simplification.gfa /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/final_contigs.fasta /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/contigs.fasta /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/first_pe_contigs.fasta /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/first_pe_contigs.fasta /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/strain_graph.gfa /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/strain_graph.gfa /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/scaffolds.fasta /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/scaffolds.fasta /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/scaffolds.paths /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/scaffolds.paths /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/assembly_graph_with_scaffolds.gfa /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/assembly_graph_with_scaffolds.gfa /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/assembly_graph.fastg /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/assembly_graph.fastg /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/K77/final_contigs.paths /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/contigs.paths


===== Copy files finished. 


===== Assembling finished. 


===== Breaking scaffolds started. 


== Running: /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/python/3.11.5/bin/python3 /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/spades/4.0.0/share/spades/spades_pipeline/scripts/breaking_scaffolds_script.py --result_scaffolds_filename /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/scaffolds.fasta --misc_dir /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/misc --threshold_for_breaking_scaffolds 3


===== Breaking scaffolds finished. 


===== Terminate started. 


===== Terminate finished. 

 * Corrected reads are in /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/corrected/
 * Assembled contigs are in /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/contigs.fasta
 * Assembled scaffolds are in /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/scaffolds.fasta
 * Paths in the assembly graph corresponding to the contigs are in /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/contigs.paths
 * Paths in the assembly graph corresponding to the scaffolds are in /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/scaffolds.paths
 * Assembly graph is in /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/assembly_graph.fastg
 * Assembly graph in GFA format is in /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/assembly_graph_with_scaffolds.gfa

======= SPAdes pipeline finished.

SPAdes log can be found here: /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/Genome_assembly/spades.log

Thank you for using SPAdes! If you use it in your research, please cite:

  Prjibelski, A., Antipov, D., Meleshko, D., Lapidus, A. and Korobeynikov, A., 2020. Using SPAdes de novo assembler. Current protocols in bioinformatics, 70(1), p.e102.
  doi.org/10.1002/cpbi.102



# QUAST - use Quast to retrieve and compare the genome assemblies between Spades, Flye and Unicycler assemblies

# command:
# Load the following modules for Quast
module load StdEnv/2020
module load gcc/9.3.0
module load python
module load quast
module load abricate

# output:
quast -o ../analysis/quast_output/ ../analysis/spades_genome_assembly/contigs.fasta \
./Flye_136x2/30-contigger/contigs.fasta \
./Unicycler_136x2/Unicycler_assembly.fasta >> ../analysis/analysis_log.txt

/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/quast/5.2.0/bin/quast -o ../analysis/quast_output/ ../analysis/spades_genome_assembly/contigs.fasta ./Flye_136x2/30-contigger/contigs.fasta ./Unicycler_136x2/Unicycler_assembly.fasta

Version: 5.2.0

System information:
  OS: Linux-3.10.0-1160.80.1.el7.x86_64-x86_64-with-glibc2.30 (linux_64)
  Python version: 3.10.2
  CPUs number: 40

Started: 2025-03-23 21:22:06

Logging to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/quast.log
NOTICE: Maximum number of threads is set to 10 (use --threads option to set it manually)

CWD: /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/data
Main parameters: 
  MODE: default, threads: 10, min contig length: 500, min alignment length: 65, min alignment IDY: 95.0, \
  ambiguity: one, min local misassembly length: 200, min extensive misassembly length: 1000

Contigs:
  Pre-processing...
  1  ../analysis/spades_genome_assembly/contigs.fasta ==> spades_genome_assembly_contigs
  2  ./Flye_136x2/30-contigger/contigs.fasta ==> 30-contigger_contigs
  3  ./Unicycler_136x2/Unicycler_assembly.fasta ==> Unicycler_assembly

2025-03-23 21:22:56
Running Basic statistics processor...
  Contig files: 
    1  spades_genome_assembly_contigs
    2  30-contigger_contigs
    3  Unicycler_assembly
  Calculating N50 and L50...
    1  spades_genome_assembly_contigs, N50 = 56712, L50 = 29, auN = 64861.9, Total length = 5028197, GC % = 50.76, # N's per 100 kbp =  0.00
    2  30-contigger_contigs, N50 = 1028955, L50 = 3, auN = 843510.8, Total length = 5112588, GC % = 50.80, # N's per 100 kbp =  0.00
    3  Unicycler_assembly, N50 = 1297187, L50 = 2, auN = 1297287.2, Total length = 5142940, GC % = 50.75, # N's per 100 kbp =  0.00
  Drawing Nx plot...
    saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/basic_stats/Nx_plot.pdf
  Drawing cumulative plot...
    saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/basic_stats/cumulative_plot.pdf
  Drawing GC content plot...
    saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/basic_stats/GC_content_plot.pdf
  Drawing spades_genome_assembly_contigs GC content plot...
    saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/basic_stats/spades_genome_assembly_contigs_GC_content_plot.pdf
  Drawing 30-contigger_contigs GC content plot...
    saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/basic_stats/30-contigger_contigs_GC_content_plot.pdf
  Drawing Unicycler_assembly GC content plot...
    saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/basic_stats/Unicycler_assembly_GC_content_plot.pdf
  Drawing Coverage histogram (bin size: 1x)...
    saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/basic_stats/coverage_histogram.pdf
  Drawing spades_genome_assembly_contigs coverage histogram (bin size: 1x)...
    saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/basic_stats/spades_genome_assembly_contigs_coverage_histogram.pdf
Done.

NOTICE: Genes are not predicted by default. Use --gene-finding or --glimmer option to enable it.

2025-03-23 21:22:59
Creating large visual summaries...
This may take a while: press Ctrl-C to skip this step..
  1 of 2: Creating PDF with all tables and plots...
  2 of 2: Creating Icarus viewers...
Done

2025-03-23 21:23:00
RESULTS:
  Text versions of total report are saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/report.txt, report.tsv, and report.tex
  Text versions of transposed total report are saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/transposed_report.txt, transposed_report.tsv, and transposed_report.tex
  HTML version (interactive tables and plots) is saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/report.html
  PDF version (tables and plots) is saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/report.pdf
  Icarus (contig browser) is saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/icarus.html
  Log is saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/quast.log

Finished: 2025-03-23 21:23:00
Elapsed time: 0:00:53.899811
NOTICEs: 2; WARNINGs: 0; non-fatal ERRORs: 0

Thank you for using QUAST!


# ABRICATE - Use abricate to compare genome assemblies in terms of plasmids, virulence genes and resistance genes

# command: 
# Create naming conventions
Spades_assembly="/home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/analysis/spades_genome_assembly/scaffolds.fasta"
Flye_assembly="/home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/data/Flye_136x2/Flye_assembly.fasta"
Unicycler_assembly="/home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/data/Unicycler_136x2/Unicycler_assembly.fasta"

# Make output directory to hold analysis results
mkdir -p ../analysis/abricate_analysis_outputs
output_dir=../analysis/abricate_analysis_outputs

# Make subdirectories for each assembly
mkdir -p $output_dir/spades_abricate_output
mkdir -p $output_dir/flye_abricate_output
mkdir -p $output_dir/unicycler_abricate_output


# Run Abricate for SPADES
abricate --db plasmidfinder $Spades_assembly > $output_dir/spades_abricate_output/Spades_plasmidfinder_output.txt 2>> ../analysis/analysis_log.txt
abricate --db resfinder $Spades_assembly > $output_dir/spades_abricate_output/Spades_resfinder_output.txt 2>> ../analysis/analysis_log.txt
abricate --db vfdb $Spades_assembly > $output_dir/spades_abricate_output/Spades_vfdb_output.txt 2>> ../analysis/analysis_log.txt

# Run Abricate for FLYE
abricate --db plasmidfinder $Flye_assembly > $output_dir/flye_abricate_output/Flye_plasmidfinder_output.txt 2>> ../analysis/analysis_log.txt
abricate --db resfinder $Flye_assembly > $output_dir/flye_abricate_output/Flye_resfinder_output.txt 2>> ../analysis/analysis_log.txt
abricate --db vfdb $Flye_assembly > $output_dir/flye_abricate_output/Flye_vfdb_output.txt 2>> ../analysis/analysis_log.txt

# Run Abricate for UNICYCLER
abricate --db plasmidfinder $Unicycler_assembly > $output_dir/unicycler_abricate_output/Unicycler_plasmidfinder_output.txt 2>> ../analysis/analysis_log.txt
abricate --db resfinder $Unicycler_assembly > $output_dir/unicycler_abricate_output/Unicycler_resfinder_output.txt 2>> ../analysis/analysis_log.txt
abricate --db vfdb $Unicycler_assembly > $output_dir/unicycler_abricate_output/Unicycler_vfdb_output.txt 2>> ../analysis/analysis_log.txt

# output:
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/quast/5.2.0/bin/quast -o ../analysis/quast_output/ ../analysis/spades_genome_assembly/contigs.fasta ./Flye_136x2/30-contigger/contigs.fasta ./Unicycler_136x2/Unicycler_assembly.fasta

Version: 5.2.0

System information:
  OS: Linux-3.10.0-1160.80.1.el7.x86_64-x86_64-with-glibc2.30 (linux_64)
  Python version: 3.10.2
  CPUs number: 44

Started: 2025-03-23 21:52:10

Logging to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/quast.log
NOTICE: Output directory already exists and looks like a QUAST output dir. Existing results can be reused (e.g. previously generated alignments)!
NOTICE: Maximum number of threads is set to 11 (use --threads option to set it manually)

CWD: /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/data
Main parameters: 
  MODE: default, threads: 11, min contig length: 500, min alignment length: 65, min alignment IDY: 95.0, \
  ambiguity: one, min local misassembly length: 200, min extensive misassembly length: 1000

Contigs:
  Pre-processing...
  1  ../analysis/spades_genome_assembly/contigs.fasta ==> spades_genome_assembly_contigs
  2  ./Flye_136x2/30-contigger/contigs.fasta ==> 30-contigger_contigs
  3  ./Unicycler_136x2/Unicycler_assembly.fasta ==> Unicycler_assembly

2025-03-23 21:52:14
Running Basic statistics processor...
  Contig files: 
    1  spades_genome_assembly_contigs
    2  30-contigger_contigs
    3  Unicycler_assembly
  Calculating N50 and L50...
    1  spades_genome_assembly_contigs, N50 = 56712, L50 = 29, auN = 64861.9, Total length = 5028197, GC % = 50.76, # N's per 100 kbp =  0.00
    2  30-contigger_contigs, N50 = 1028955, L50 = 3, auN = 843510.8, Total length = 5112588, GC % = 50.80, # N's per 100 kbp =  0.00
    3  Unicycler_assembly, N50 = 1297187, L50 = 2, auN = 1297287.2, Total length = 5142940, GC % = 50.75, # N's per 100 kbp =  0.00
  Drawing Nx plot...
    saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/basic_stats/Nx_plot.pdf
  Drawing cumulative plot...
    saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/basic_stats/cumulative_plot.pdf
  Drawing GC content plot...
    saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/basic_stats/GC_content_plot.pdf
  Drawing spades_genome_assembly_contigs GC content plot...
    saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/basic_stats/spades_genome_assembly_contigs_GC_content_plot.pdf
  Drawing 30-contigger_contigs GC content plot...
    saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/basic_stats/30-contigger_contigs_GC_content_plot.pdf
  Drawing Unicycler_assembly GC content plot...
    saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/basic_stats/Unicycler_assembly_GC_content_plot.pdf
  Drawing Coverage histogram (bin size: 1x)...
    saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/basic_stats/coverage_histogram.pdf
  Drawing spades_genome_assembly_contigs coverage histogram (bin size: 1x)...
    saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/basic_stats/spades_genome_assembly_contigs_coverage_histogram.pdf
Done.

NOTICE: Genes are not predicted by default. Use --gene-finding or --glimmer option to enable it.

2025-03-23 21:52:21
Creating large visual summaries...
This may take a while: press Ctrl-C to skip this step..
  1 of 2: Creating PDF with all tables and plots...
  2 of 2: Creating Icarus viewers...
Done

2025-03-23 21:52:22
RESULTS:
  Text versions of total report are saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/report.txt, report.tsv, and report.tex
  Text versions of transposed total report are saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/transposed_report.txt, transposed_report.tsv, and transposed_report.tex
  HTML version (interactive tables and plots) is saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/report.html
  PDF version (tables and plots) is saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/report.pdf
  Icarus (contig browser) is saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/icarus.html
  Log is saved to /scratch/mirza/Genomics_Assignment_3/Assignment3_W25_files/analysis/quast_output/quast.log

Finished: 2025-03-23 21:52:22
Elapsed time: 0:00:11.803156
NOTICEs: 3; WARNINGs: 0; non-fatal ERRORs: 0

Thank you for using QUAST!
Using nucl database plasmidfinder:  460 sequences -  2021-Oct-12
Processing: /home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/analysis/spades_genome_assembly/scaffolds.fasta
Found 12 genes in /home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/analysis/spades_genome_assembly/scaffolds.fasta
Tip: the abricate manual is at https://github.com/tseemann/abricate/blob/master/README.md
Done.
Using nucl database resfinder:  3077 sequences -  2021-Oct-12
Processing: /home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/analysis/spades_genome_assembly/scaffolds.fasta
Found 5 genes in /home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/analysis/spades_genome_assembly/scaffolds.fasta
Tip: abricate can also find virulence factors; use --list to see all supported databases.
Done.
Using nucl database vfdb:  2597 sequences -  2021-Oct-12
Processing: /home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/analysis/spades_genome_assembly/scaffolds.fasta
Found 65 genes in /home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/analysis/spades_genome_assembly/scaffolds.fasta
Tip: you can use the --summary option to combine reports in a presence/absence matrix.
Done.
Using nucl database plasmidfinder:  460 sequences -  2021-Oct-12
Processing: /home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/data/Flye_136x2/Flye_assembly.fasta
Found 13 genes in /home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/data/Flye_136x2/Flye_assembly.fasta
Tip: remember that abricate is unable to find AMR-mediated SNPs.
Done.
Using nucl database resfinder:  3077 sequences -  2021-Oct-12
Processing: /home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/data/Flye_136x2/Flye_assembly.fasta
Found 5 genes in /home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/data/Flye_136x2/Flye_assembly.fasta
Tip: have a suggestion for abricate? Tell me at https://github.com/tseemann/abricate/issues
Done.
Using nucl database vfdb:  2597 sequences -  2021-Oct-12
Processing: /home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/data/Flye_136x2/Flye_assembly.fasta
Found 66 genes in /home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/data/Flye_136x2/Flye_assembly.fasta
Tip: it's important to realise abricate uses DNA matching, not AA.
Done.
Using nucl database plasmidfinder:  460 sequences -  2021-Oct-12
Processing: /home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/data/Unicycler_136x2/Unicycler_assembly.fasta
Found 11 genes in /home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/data/Unicycler_136x2/Unicycler_assembly.fasta
Tip: the abricate manual is at https://github.com/tseemann/abricate/blob/master/README.md
Done.
Using nucl database resfinder:  3077 sequences -  2021-Oct-12
Processing: /home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/data/Unicycler_136x2/Unicycler_assembly.fasta
Found 5 genes in /home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/data/Unicycler_136x2/Unicycler_assembly.fasta
Tip: remember that abricate is unable to find AMR-mediated SNPs.
Done.
Using nucl database vfdb:  2597 sequences -  2021-Oct-12
Processing: /home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/data/Unicycler_136x2/Unicycler_assembly.fasta
Found 66 genes in /home/mirza/scratch/Genomics_Assignment_3/Assignment3_W25_files/data/Unicycler_136x2/Unicycler_assembly.fasta
Tip: have a suggestion for abricate? Tell me at https://github.com/tseemann/abricate/issues
Done.
