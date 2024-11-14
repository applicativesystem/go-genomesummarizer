# go-genomesummarizer
- golang genome summarizer
- allows you to do all calculative stats on the genome and metagenome annotation. done
- estimates all descriptive and variance and all linear algebra options included. partial code completed
- writes a linear regression and a decision tree. partial code completed.
- one stop solution for the given gff and the fasta file, will summarize the genome, and run the machine learning models on the genome and metagenome annotations.

```
➜  go-genomesummarizer git:(main) ✗ go run main.go -h
Run the genome analyzer function and summarize the genome

Usage:
  genomesummarize [flags]

Flags:
  -G, --genomesummarize string   summarize the genome (default "genomeannotate")
  -h, --help                     help for genomesummarize
exit status
```

Gaurav Sablok 
