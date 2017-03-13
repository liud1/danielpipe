package Daniel::Blast;
use strict;
use warnings;
use threads;
use Cwd  qw(abs_path);
use File::Basename qw(dirname);
use Exporter qw(import);
use Bio::Search::Result::BlastResult;
use Bio::SearchIO; 
our @EXPORT_OK = qw(create_blast_shell parseblast);
 
sub create_blast_shell {
  my $blast_root = shift // '/export/arrayPRO1/commonSoftwares/blast+'; # the location of blast bin with default /export/arrayPRO1/commonSoftwares/blast+
  my $shell_fn    = shift; # the file name of the blast shell 
  my $type = shift;  # type of blast [blastp|blastn|blastx]
  my $query = shift; # Query file name
  my $db    = shift; # database_name  
  my $outfmt = shift // '0'; # alignment view options:0 = pairwise,6 = tabular,8 = Text ASN.1
  my $evalue = shift // '1e-15' ; # Expect value (E) for saving hits with default 1e-15  
  my $num_alignments = shift // '1'; # Show alignments for this number of database sequences with default 1   
  my $num_descriptions = shift // '1'; # Show one-line descriptions for this number of database sequences with default 1    
  my $query_gencode = shift // '1' ; # Genetic code to translate query [ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt] 1: Standard, 11: Bacterial, Archaeal and Plant Plastid with default 1  
  my $num_threads = shift // '10'; # the thread number with default 10
 
 #!/bin/sh
#$ -N blast
#$ -e /export/EC1680U/Daniel_DeNovo/HDARES_YCTsaiLab/B01/output/3.blastp/B01.gene.aa_9.sh.err
#$ -o /export/EC1680U/Daniel_DeNovo/HDARES_YCTsaiLab/B01/output/3.blastp/B01.gene.aa_9.sh.out
#/export/arrayPRO1/commonSoftwares/blast+/blastp -query /export/EC1680U/Daniel_DeNovo/HDARES_YCTsaiLab/B01/output/3.blastp/B01.gene.aa_9.fa -db /export/arrayPRO1/blastdb/refSeq_bacteria_20150710/RefSeq_bacteriaAllaa.db -outfmt 0 -out /export/EC1680U/Daniel_DeNovo/HDARES_YCTsaiLab/B01/output/3.blastp/B01.gene.aa_9.blast -task blastp -num_threads 10 -evalue 1e-15 -num_alignments 1 -num_descriptions 1 -seg no
  
  my $out   = $query; # Output file name
  $out=~s/.fa/.blast/;
  my $cmd;  
  if($type eq "blastp"){
    $cmd="$blast_root/blastp -query $query -db $db -outfmt $outfmt -out $out -task blastp -num_threads $num_threads -evalue $evalue -num_alignments $num_alignments -num_descriptions $num_descriptions ";
  }
  if($type eq "blastn"){
    $cmd="$blast_root/blastn -query $query -db $db -outfmt $outfmt -out $out -task blastn -num_threads $num_threads -evalue $evalue -num_alignments $num_alignments -num_descriptions $num_descriptions ";
  }
  if($type eq "blastx"){
    $cmd="$blast_root/blastx -query $query -db $db -query_gencode $query_gencode -outfmt $outfmt -out $out -num_threads $num_threads -evalue $evalue -num_alignments $num_alignments -num_descriptions $num_descriptions ";
  }
  
  unless (-e $shell_fn) {
    open (OUT, ">", $shell_fn) or die "Can't create $shell_fn file";
    print OUT "#!/bin/sh\n";
    print OUT "#\$ -N blast\n";
    print OUT "#\$ -e $shell_fn.err\n";
    print OUT "#\$ -o $shell_fn.out\n";
    print OUT "$cmd\n";
    close OUT;
  }
}

sub parseblast {
  my $blast_dir = shift;# the directory that have blast result *.blast
  my $output    = shift;# the output file of the parse blast results
  my @fn;               # blast files
  my $file;             # blast file
  my $report;           # the report of blast file that may have more than one query
  my $result;           # the result is the entire analysis for a single query sequence
  my $algorithm = "";        # algorithm
  my $query_name = "";       # query name
  my $query_length = "";     # query length
  my $hits = "";             # the number of hit
  my $hit = "";              # hits are sequences in the searched database which could be aligned to the query sequence and met the minimal search parameters
  my $hit_name = "";         # hit name 
  my $hsp = "";              # high-scoring segment Pairs 
  my $query_start = "";      # query start of hsp
  my $query_end = "";        # query end pf hsp
  my $percent_identity = ""; # % identical 
  my $num_identical = "";    # number of identical residues 
  my $hsp_length = "";       # Length of the HSP (including gaps) alias for length('total') 
  my $frac_conserved = "";   # fraction conserved (conservative and identical replacements aka "fraction similar") (only valid for Protein alignments will be same as frac_identical) 
  my $num_conserved = "";    # number of conserved (conservative replacements, aka "similar") residues 
  my $gaps = "";             # number of gaps 
  my $bits = "";             # score in bits 
  my $expect = "";           # alias for evalue() 
  my $hit_strand = "";       # strand of the hit 
  my $query_strand = "";     # strand of the query    
  my $hit_start = "";        # hit start of hsp
  my $hit_end = "";          # hit end pf hsp    
  my $hit_length = "";    	# length of the hit sequence 
  my $description = "";      # hit description 
  my $query_string = "";     # query string from alignment  
  my @info_items;     # the output info item
  
  @fn = glob("$blast_dir/*.blast");
  unless (-e $output) {
    open (OUT, '>',$output) or die "Can't create $output file\n";
	print OUT "Type\tQ_Name\tQ_Length\tQ_Start\tQ_End\tId%\tIdentity\tPos%\tPositives\tGaps\tScore\tExpectValue\tHit Strand\tQuery Strand\tS_Start\tS_End\tS_Length\tS_Namen\ttranslatedAA\n";
    foreach $file (@fn) {    	  
      $report = Bio::SearchIO->new( -file=>$file, -format => 'blast');  
      $algorithm = ""; $query_name = ""; $query_length = ""; $query_start = ""; $query_end = ""; $percent_identity = ""; $num_identical = ""; $hsp_length = ""; $frac_conserved = "";
	  $num_conserved = ""; $hsp_length = "";  $gaps = ""; $hsp_length = ""; $bits = ""; $expect = ""; $hit_strand = ""; $query_strand = "";$hit_start = ""; $hit_end = ""; $hit_length = ""; 
	  $hit_name = ""; $description = ""; $query_string = ""; 	  
      while ($result = $report->next_result) {
        $query_name = $result->query_name;
        $algorithm = $result->algorithm;
        $query_length = $result->query_length ;
        $hits = $result->hits ;            
        if ($hits == 0) {          
        } else {
          $hit               = $result->next_hit;
          $description       = $hit->description;
          $hit_length        = $hit->length;
          $hit_name          = $hit->name;
          $hsp               = $hit->next_hsp;      	
          $query_start       = $hsp->start('query');
          $query_end         = $hsp->end('query');
          $percent_identity  = sprintf("%.f", $hsp->percent_identity);
          $num_identical     = $hsp->num_identical  ;
          $hsp_length        = $hsp->hsp_length ; 
          $frac_conserved    = sprintf("%.f", 100 * $hsp->frac_conserved);
          $num_conserved     = $hsp->num_conserved;
          $gaps              = $hsp->gaps ;
          $bits              = $hsp->bits;
          $expect            = $hsp->expect;
          $hit_strand        = $hsp->strand('hit') ;
          $query_strand      = $hsp->strand('query') ;
          $hit_start         = $hsp->start('hit');
          $hit_end           = $hsp->end('hit');   
          $query_string      = $hsp->query_string();
          $query_string =~s/-//g; 	
          $query_string =~s/ //g; 	                          
        }    
      }
	  
	  @info_items=($algorithm, $query_name, $query_length, $query_start, $query_end, "$percent_identity\%", "$num_identical\/$hsp_length", "$frac_conserved\%", "$num_conserved\/$hsp_length", "$gaps\/$hsp_length", $bits, $expect, $hit_strand, $query_strand, $hit_start, $hit_end, $hit_length, "$hit_name$description", $query_string);
	  print OUT join("\t", @info_items)."\n";
    }	
  }
  close OUT;
}
