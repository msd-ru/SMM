#1. Identification of UMI families. Command line parameters (positional, whithout keys)
main.py umi <S file> <Chrom (use '-' for None)>

A report is generated for each UMI family:
UMI, Start, End, strang count +, strang count -.

#2. Variant calling. Command line parameters (positional, whithout keys)
main.py vc <S file> <Chrom (use '-' for None)> <G file> <Reference> <SNP database> <G VCF records file>

Result formats are presented in the project description.

The additional parameters are listed in the file Params.ini. Their meaning is clear enough from the names and values.
Some explanations are given below. 

The family includes all reads with a common UMI, the areas [nStart, nEnd) of which intersect at least at maxPosDist positions.
A related pair read1, read2 is united if they contents at least one common position.

Qualified families are determined from the size of each strand +,- (each >= minStrandCount)
or the total size of the family ( >= minFamilyCount, if this parameter is positive).

File names and chromosomes can also be specified by parameters in Params.ini (see examples there). But these parameters can
"overlapped" by command line parameters. You can set up to 6 of them, see the main.py file and comments there.






