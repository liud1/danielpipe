package Daniel::Report;
use strict;
use warnings;
use threads;
 
use Exporter qw(import);
use lib '/export/EC1680U/perl/lib'; # Global
use Daniel::Utility qw(checkjobid create_fullshell process_cmd);
 
our @EXPORT_OK = qw(fastqstats samstats procFastqStats procSamStats createFastqStatsHtml createSamStatsHtml procTopHatAlignStats createTopHatAlignStatsHtml createEnrichHtml procRSEMAlignStats createRSEMAlignStatsHtml procHisatStats createHisatAlignStatsHtml procGffcompareStats createGffcompareStatsHtml targetcoverage createBarcodeHtml createQCHtml createOtuHtml createTaxonomyHtml);
$ENV{PATH} = "$ENV{PATH}:/opt/sge/bin/lx-amd64/"; 

sub fastqstats {
  my $fastqstats_root = shift ;
  my $sample = shift;
  my $r1 = shift;
  my $r2 = shift;
  my $stats_out = shift;
  my $shell;
  my $err;
  my $log;      
  my $cmd;
  my $qname = "fastqstats";  	
  $shell = "$stats_out/$sample-fastqstats.sh";
  $err = "$stats_out/$sample-fastqstats.err";
  $log = "$stats_out/$sample-fastqstats.log";
  $cmd=$cmd."export PATH=/opt/sge/bin/lx-amd64:\$PATH \n";
  $cmd=$cmd."export LD_LIBRARY_PATH=/opt/zlib-1.2.8/lib:\$LD_LIBRARY_PATH \n";  
  $cmd=$cmd."export LD_LIBRARY_PATH=/opt/glibc-2.14/lib:\$LD_LIBRARY_PATH \n";  
  $cmd = $cmd."$fastqstats_root/fastq-stats $r1 2>&1 | tee $stats_out/$sample.R1.fq-stats \n";  
  $cmd = $cmd."$fastqstats_root/fastq-stats $r2 2>&1 | tee $stats_out/$sample.R2.fq-stats \n";  
  unless (-e $err) {create_fullshell ($shell, $qname, $err, $log, $cmd); `qsub $shell >> ./jobs`;}    
  checkjobid("");
}
sub procFastqStats {
  my $samples = shift;
  my $stats_out = shift;
  my $data;  # data structure for fastq-stats output 
  my $ln;
  my $sample;
  my $type;
  my @types = ("R1", "R2"); # It means R1 & R2 tag
   
  foreach $sample (@$samples) {    
    foreach $type (@types) {	      
      open (IN, '<', "$stats_out/$sample.$type.fq-stats") or die "Can't open $stats_out/$sample.$type.fq-stats file\n";
      while (<IN>) {
        $ln = $_;
	    chomp($ln);
	    if ($ln=~/reads\s(\d.*)/) {
	      $data->{reads}->{$sample}->{$type} = $1;
	    }
	    if ($ln=~/len\s(\d.*)/) {
	      $data->{len}->{$sample}->{$type} = $1;
	    }	
	    if ($ln=~/\%dup\s(\d.*)/) {
	      $data->{pct_dup}->{$sample}->{$type} = $1;
	    }	
	    if ($ln=~/qual mean\s(\d.*)/) {
	      $data->{qv_mean}->{$sample}->{$type} = $1;
	    }			
	    if ($ln=~/total bases\s(\d.*)/) {
	      $data->{total_base}->{$sample}->{$type} = $1;
	    }				
      }
      close IN;
	}
  }
  return $data;
}
sub createFastqStatsHtml {
  my $samples = shift;
  my $data = shift;
  my $fastqstats_template = shift;
  my $output = shift;
  my $cmd;
  my $content;
  my @items = ("reads","len","pct_dup","qv_mean","total_base");
  my $sample;
  my $r1_item;
  my $r2_item;
  my $title = "Fastq Statistics";
  $content = $content."<table>";
  foreach $sample (@$samples) {    
    $content = $content."<tr><td colspan=3>$sample</td></tr>";
	$content = $content."<tr><td></td><td>$sample\-R1</td><td>$sample\-R2</td></tr>";
	foreach my $item (@items) {
      $content = $content."<tr><td>$item</td>";
      $r1_item =  $data->{$item}->{$sample}->{"R1"};
	  $r2_item =  $data->{$item}->{$sample}->{"R2"};		  
      $content = $content."<td>$r1_item</td><td>$r2_item</td></tr>";		  
	}		
  } 
  $content = $content."</table>";  
  $cmd = "sed 's|Title|$title|' $fastqstats_template | sed 's|Content|$content|' > $output/fastqstats.html";
  process_cmd ($cmd); 
}
sub samstats {
  my $fastqstats_root = shift ;
  my $sample = shift;
  my $bam = shift;
  my $stats_out = shift;
  my $shell;
  my $err;
  my $log;      
  my $cmd;
  my $qname = "samstats";  	
  $shell = "$stats_out/$sample-samstats.sh";
  $err = "$stats_out/$sample-samstats.err";
  $log = "$stats_out/$sample-samstats.log";
  $cmd=$cmd."export PATH=/opt/sge/bin/lx-amd64:\$PATH \n";
  $cmd=$cmd."export LD_LIBRARY_PATH=/opt/zlib-1.2.8/lib:\$LD_LIBRARY_PATH \n";  
  $cmd=$cmd."export LD_LIBRARY_PATH=/opt/glibc-2.14/lib:\$LD_LIBRARY_PATH \n";  
  $cmd = $cmd."$fastqstats_root/sam-stats $bam  2>&1 | tee $stats_out/$sample-sam-stats.log \n";    
  unless (-e $err) {create_fullshell ($shell, $qname, $err, $log, $cmd); `qsub $shell >> ./jobs`;}    
  checkjobid("");
}
sub procSamStats {
  my $samples = shift;
  my $stats_out = shift;
  my $data;  # data structure for fastq-stats output 
  my $ln;
  my $sample;
  foreach $sample (@$samples) {  
    open (IN, '<', "$stats_out/$sample-sam-stats.log") or die "Can't open $stats_out/$sample-sam-stats.log\n";
    while (<IN>) {
      $ln = $_;chomp($ln);
	  
	  if ($ln=~/reads\s(\d.*)/) {
	    $data->{$sample}->{reads} = $1;
	  }
	  if ($ln=~/mapped reads\s(\d.*)/) {
	    $data->{$sample}->{mapped_reads} = $1;
	  }	
	  if ($ln=~/pct align\s(\d.*)/) {
	    $data->{$sample}->{pct_align} = $1;
	  }	
	  if ($ln=~/mapped bases\s(\d.*)/) {
	    $data->{$sample}->{mapped_bases} = $1;
	  }	
	  if ($ln=~/mapq mean\s(\d.*)/) {
	    $data->{$sample}->{mapq_mean} = $1;
	  }					  	  
    }
    close IN;
  }
  return $data;
}
sub createSamStatsHtml {
  my $samples = shift;
  my $data = shift;
  my $samstats_template = shift;
  my $output = shift;
  my $cmd;
  my @items = ("mapped_reads","pct_align","mapped_bases","mapq_mean"); 
  my $sample;
  my $item_value;
  my $title = "Alignment Statistics";
  my $content = "<table><tr><td>Sample</td><td>Metrics</td></tr>";  
  foreach $sample (@$samples) {    
    $content = $content."<tr><td colspan=2>$sample</td></tr>";	
	foreach my $item (@items) {
      $item_value  = $data->{$sample}->{$item} ;
      $content=$content."<tr><td>$item</td><td>$item_value</td></tr>";
    }
  } 
  $content = $content."</table>";  
  $cmd = "sed 's|Title|$title|' $samstats_template | sed 's|Content|$content|' > $output/samstats.html";
  process_cmd ($cmd); 
}
sub procTopHatAlignStats {
  my $samples = shift;
  my $tophat_out = shift;
  my $data;  # data structure for tophat align_summary.txt 
  my $ln;
  my $sample;
  my $type; # 1 (R1) or 2 (R2)
  foreach $sample (@$samples) { 
    open (IN, '<', "$tophat_out/$sample/align_summary.txt") or die "Can't open $tophat_out/$sample/align_summary.txt file\n";
	while (<IN>) {
	  $ln = $_; 
	  chomp($ln);
	  if ($ln=~/^(.*) reads/) {if ($1 eq "Left") {$type = "R1"} else {$type = "R2"};}
	  if ($ln=~/Input\s+:\s(.*)/) {$data->{$sample}->{$type}->{input} = $1;} # Input Reads
	  if ($ln=~/Mapped\s+:\s(.*)/) {$data->{$sample}->{$type}->{mapped} = $1;} # Mapped Reads
	  if ($ln=~/(.*) overall read mapping rate./) {$data->{$sample}->{mappingrate} = $1;} # Mapping Rate
	  if ($ln=~/(.*) concordant pair alignment rate./) {$data->{$sample}->{concordant} = $1;} # Concordant Rate
	}
	close IN;	
  }
  return $data;
}
sub createTopHatAlignStatsHtml {
  # create align stats of tophat align_summary.txt
  my $samples = shift;
  my $data = shift;
  my $samstats_template = shift;
  my $output = shift;
  my $cmd;
  my $concordant;
  my $content;  
  my $item;
  my @items = ("input","mapped");
  my $item_value;
  my $mappingrate;
  my $sample;
  my $title = "Alignment Statistics";
  my $type;
  my @types = ("R1", "R2");
  $content = $content."<table>";
  foreach $sample (@$samples) { 
    $mappingrate = $data->{$sample}->{mappingrate} ;
	$concordant = $data->{$sample}->{concordant} ;
    $content = $content."<tr><td colspan=3>$sample</td></tr>";
	$content = $content."<tr><td>Type</td><td>input</td><td>mapped</td></tr>";  	
    foreach $type (@types) {
	  $content = $content."<tr><td>$type</td>";  	
	  foreach $item (@items) {
	    $item_value =  $data->{$sample}->{$type}->{$item};	    
		$content = $content."<td>$item_value</td>"; 
	  }
	  $content = $content."</tr>";  	  
	}
	$content = $content."<tr><td>overall read mapping rate</td><td colspan=2>$mappingrate</td></tr>";
	$content = $content."<tr><td> concordant pair alignment rate</td><td colspan=2>$concordant</td></tr>";
  }
  $content = $content."</table>";  
  $cmd = "sed 's|Title|$title|' $samstats_template | sed 's|Content|$content|' > $output/samstats.html";
  process_cmd ($cmd); 
}
sub procRSEMAlignStats {
  my $samples = shift;
  my $calculateexpression = shift;
  my $data;  # data structure for tophat align_summary.txt 
  my $ln;
  my $sample;
  foreach $sample (@$samples) { 
    open (IN, '<', "$calculateexpression/$sample.CalculateExpression.err") or die "Can't open $calculateexpression/$sample.CalculateExpression.err file\n";
	while (<IN>) {
	  $ln = $_; 
	  chomp($ln);
	  if ($ln=~/^(.*) reads/) {$data->{$sample}->{input} = $1;} # Total Pair reads
	  if ($ln=~/(.*) \(.*\) aligned concordantly 0 times/) {$data->{$sample}->{unmapped} = $1;} # (unmapped) aligned concordantly 0 times
	  if ($ln=~/(.*) \(.*\) aligned concordantly exactly 1 time/) {$data->{$sample}->{mappedonce} = $1;} # (mappedonce) aligned concordantly exactly 1 time
	  if ($ln=~/(.*) \(.*\) aligned concordantly >1 times/) {$data->{$sample}->{mappedmore} = $1;} # (mappedmore) aligned concordantly >1 times
	  if ($ln=~/(.*) overall alignment rate/) {$data->{$sample}->{maprate} = $1;} # (maprate) overall alignment rate
	}
	close IN;	
  }
  return $data;
}
sub createRSEMAlignStatsHtml {
# create align stats of tophat align_summary.txt
  my $samples = shift;
  my $data = shift;
  my $samstats_template = shift;
  my $output = shift;
  my $cmd;
  my $content;  
  my $input;
  my $mappedonce;
  my $mappedmore;
  my $maprate;
  my $sample;
  my $title = "Alignment Statistics";
  my $unmapped;
  $content = $content."<table>";
  foreach $sample (@$samples) { 
    $input      = $data->{$sample}->{input} ;
	$unmapped   = $data->{$sample}->{unmapped};
	$mappedonce = $data->{$sample}->{mappedonce};
	$mappedmore = $data->{$sample}->{mappedmore};
	$maprate    = $data->{$sample}->{maprate};
    $content = $content."<tr><td colspan=2>$sample</td></tr>";
	$content = $content."<tr><td>Total Reads Pair</td><td>$input</td></tr>";  	
    $content = $content."<tr><td>aligned concordantly 0 times</td><td>$unmapped</td></tr>";  	
    $content = $content."<tr><td>aligned concordantly exactly 1 time</td><td>$mappedonce</td></tr>";  	
    $content = $content."<tr><td>aligned concordantly >1 times</td><td>$mappedmore</td></tr>";  	
    $content = $content."<tr><td>overall alignment rate</td><td>$maprate</td></tr>";  		
  }
  $content = $content."</table>";  
  $cmd = "sed 's|Title|$title|' $samstats_template | sed 's|Content|$content|' > $output/alignstats.html";
  process_cmd ($cmd); 
}
sub procHisatStats {
  my $samples = shift;
  my $hisat = shift;
  my $data;  # data structure for hisat align_summary.txt 
  my $ln;
  my $sample;
  my $sentinel = "N";
  foreach $sample (@$samples) { 
    $sentinel = "N";
    open (IN, '<', "$hisat/$sample.hisat.err") or die "Can't open $hisat/$sample.hisat.err file\n";
	while (<IN>) {
	  $ln = $_; 
	  chomp($ln);
	  if  ($ln=~/overall alignment rate/) {
	    $data->{$sample}->{alncontent} = $data->{$sample}->{alncontent}."$ln"."AAA";
		$sentinel = "Y";
	  }
	  if ($sentinel eq "N" ) {
	    if (defined $data->{$sample}->{alncontent}) {
		  $data->{$sample}->{alncontent} = $data->{$sample}->{alncontent}."$ln"."AAA";
		} else {
		  $data->{$sample}->{alncontent} = "$ln"."AAA";
		}
	  }
	}
	close IN;	
  }
  return $data;
}
sub createHisatAlignStatsHtml {
# create align stats of hisat align_summary.txt
  my $samples = shift;
  my $data = shift;
  my $samstats_template = shift;
  my $output = shift;
  my $alncontent;
  my $cmd;
  my $content;  
  my $sample;
  my $title = "Alignment Statistics";
  $content = $content."<table>";
  foreach $sample (@$samples) { 
    $alncontent = $data->{$sample}->{alncontent} ;
    $content = $content."<tr><td>$sample</td></tr>";
	$content = $content."<tr><td><pre>$alncontent</pre></td></tr>";  	 		
  }
  $content = $content."</table>";  
  $cmd = "sed 's|Title|$title|' $samstats_template | sed 's|Content|$content|' | sed 's|AAA|\\n|g'> $output/alignstats.html";
  process_cmd ($cmd); 
}
sub procGffcompareStats {
  my $gffcompare = shift;
  my $data;  # data structure for hisat align_summary.txt 
  my $ln;
  my $sentinel = "N";
  open (IN, '<', "$gffcompare/gffcompare.stats") or die "Can't open $gffcompare/gffcompare.stats file\n";
  while (<IN>) {
	$ln = $_; 
	chomp($ln);
	$ln=~s/,/ /;
	if  ($ln=~/Summary for dataset/) {
	  $sentinel = "Y";
	}
    if ($sentinel eq "Y" ) {
	  if (defined $data->{gffcomparecontent}) {
	    $data->{gffcomparecontent} = $data->{gffcomparecontent}."$ln"."AAA";
	  } else {
		$data->{gffcomparecontent} = "$ln"."AAA";
	  }
	}
  }
  close IN;	  
  return $data;  
}
sub createGffcompareStatsHtml {
  my $data = shift;
  my $gffcompare_template = shift;
  my $output = shift;
  my $gffcomparecontent;
  my $cmd;
  my $content;  
  my $sample;
  my $title = "GffCompare Statistics";
  $gffcomparecontent = $data->{gffcomparecontent};
  $content = $content."<table>";
  $content = $content."<tr><td><pre>$gffcomparecontent</pre></td></tr>";  	 		
  
  $content = $content."</table>";  
  #$cmd = "sed 's|Title|$title|' $gffcompare_template | sed 's|Content|$content|' | sed 's|AAA|\\n|g'> $output/gffcompare.html";
  $cmd = "sed 's,Title,$title,' $gffcompare_template | sed 's,Content,$content,' | sed 's,AAA,\\n,g'> $output/gffcompare.html";  
  process_cmd ($cmd); 
}
sub createEnrichHtml {
  # create GO & KEGG enrichment
  my $enrich_template = shift;
  my $output = shift;
  my $path = shift;
  my $cmd;
  $cmd = "sed 's|Path|$path|' $enrich_template > $output/enrich.html";  
  process_cmd ($cmd); 
}
sub targetcoverage {
  my $sample = shift;
  my $bam   = shift; 
  my $rscript = shift;
  my $bedtools_root = shift;
  my $bed = shift;
  my $stats_out = shift;
  my $output = shift;
  my $type = shift;
  my $cmd;
  
  if ($type eq "wes") {
    unless (-e "$output/$sample.exome.coverage.hist.txt") {
      $cmd = "$bedtools_root/bedtools coverage -hist -abam $bam -b $bed  > $stats_out/$sample.exome.coverage.hist.txt";
      process_cmd($cmd);
	  $cmd = "grep all $stats_out/$sample.exome.coverage.hist.txt > $stats_out/$sample.all.exome.coverage.hist.txt";
	  process_cmd($cmd);
    }  
    $cmd = "Rscript $rscript/coverage-wes.r $stats_out/$sample.all.exome.coverage.hist.txt $output/$sample-covreage.jpg";
	process_cmd($cmd); 
  } elsif ($type eq "amplicon") {
    unless (-e "$output/$sample.amplicon.coverage.hist.txt") {
      $cmd = "$bedtools_root/bedtools coverage -hist -abam $bam -b $bed  > $stats_out/$sample.amplicon.coverage.hist.txt";
      process_cmd($cmd);
	  $cmd = "grep all $stats_out/$sample.amplicon.coverage.hist.txt > $stats_out/$sample.all.amplicon.coverage.hist.txt";
	  process_cmd($cmd);
    }  
    $cmd = "Rscript $rscript/coverage-amplicon.r $stats_out/$sample.all.amplicon.coverage.hist.txt $output/$sample-covreage.jpg";
	process_cmd($cmd);     
  }
  
    
}
sub createBarcodeHtml {
  # This is for 16S barcode statistics
  my $barcodelist = shift;
  my $barcode_template = shift;  
  my $output = shift;
  my $content;
  my $count = 1;
  my $cmd;
  my $ln;
  my $idx1;
  my $idx2;
  my @tab;
  my $title = "Barcode Stats";
  my $sample;
  
  $content = $content."<table>";
  $content = $content."<tr><td>#</td><td>sampleName</td><td>index1</td><td>index2</td></tr>";
  open (IN, '<', $barcodelist) or die "Can't create $barcodelist file\n";
  while (<IN>) {
    $ln = $_; chomp($ln);
	if ($ln=~/#dualbarcodes/) {next;}
	@tab = split(/\t/, $ln);
	$content = $content."<tr>";
	$idx1 = $tab[1];
	$idx2 = $tab[2];
	$sample = $tab[3];
	$content = $content."<tr><td>$count</td><td>$sample</td><td>$idx1</td><td>$idx2</td></tr>";
	$count = $count + 1;
  }
  close IN;
  $content = $content."</table>";
  $cmd = "sed 's|Title|$title|' $barcode_template | sed 's|Content|$content|'> $output/barcode.html";  
  process_cmd ($cmd); 
}
sub createQCHtml {
  # This is for 16S QC
  my $readstat = shift;
  my $readstat_template = shift;  
  my $output = shift;
  my $content;
  my $cmd;
  my $ln;
  my @tab;
  my $title = "QC Stats";
  my $sample;
  
  $content = $content."<table>";
  $content = $content."<tr><td>SampleID</td><td>CleanPE(#)</td><td>Joined(#)</td><td>Qualified(#)</td><td>NoChime(#)</td><td>AvgQ</td><td>Q20\%</td><td>Q30%</td><td>GC\%</td><td>Effective\%</td></tr>";
  open (IN, '<', $readstat) or die "Can't create $readstat file\n";
  while (<IN>) {
    $ln = $_; chomp($ln);
	if ($ln=~/SampleID/) {next;}
	@tab = split(/\t/, $ln);	
	$tab[1]=~s/\,//;
	if ($tab[1] < 50000) {
	  warn "$tab[0] CleanPE below 50,000!";     # does show up
	  #exit();
	}
	$content = $content."<tr><td>$tab[0]</td><td>$tab[1]</td><td>$tab[2]</td><td>$tab[3]</td><td>$tab[4]</td><td>$tab[7]</td><td>$tab[8]</td><td>$tab[9]</td><td>$tab[10]</td><td>$tab[11]</td></tr>";
  }
  close IN;
  $content = $content."</table>";
  $cmd = "sed 's|Title|$title|' $readstat_template | sed 's|Content|$content|'> $output/readstat.html";  
  process_cmd ($cmd); 
}
sub createOtuHtml {
  # This is for 16S OTU
  my $otustat = shift;
  my $otustat_template = shift;  
  my $output = shift;
  my $content;
  my $cmd;
  my $ln;
  my @tab;
  my $title = "OTU Stats";
  my $sample;
  
  $content = $content."<table>";
  $content = $content."<tr><td>SampleID</td><td>EffectiveRead(#)</td><td>EffectiveBases(bp)</td><td>TaxonReads</td><td>NrReads(#)</td><td>OTUs(#)</td></tr>";
  open (IN, '<', $otustat) or die "Can't create $otustat file\n";
  while (<IN>) {
    $ln = $_; chomp($ln);
	if ($ln=~/SampleID/) {next;}
	@tab = split(/\t/, $ln);	
	$content = $content."<tr><td>$tab[0]</td><td>$tab[1]</td><td>$tab[2]</td><td>$tab[3]</td><td>$tab[4]</td><td>$tab[5]</td></tr>";
  }
  close IN;
  $content = $content."</table>";
  $cmd = "sed 's|Title|$title|' $otustat_template | sed 's|Content|$content|'> $output/otustat.html";  
  process_cmd ($cmd); 
}
sub createTaxonomyHtml {
  # This is for 16S Taxonomy
  my $taxonomystat = shift;
  my $taxonomystat_template = shift;  
  my $output = shift;
  my $content;
  my $cmd;
  my $ln;
  my @tab;
  my $title = "Taxonomy Stats";
  my $sample;
  
  $content = $content."<table>";
  $content = $content."<tr><td>SampleID</td><td></td><td>Kindom</td><td>Phylum</td><td>Class</td><td>Order</td><td>Family</td><td>Genus</td><td>Species</td></tr>";
  open (IN, '<', $taxonomystat) or die "Can't create $taxonomystat file\n";
  while (<IN>) {
    $ln = $_; chomp($ln);
	if ($ln=~/SampleID/) {next;}
	@tab = split(/\t/, $ln);	
	$content = $content."<tr><td>$tab[0]</td><td>$tab[1]</td><td>$tab[2]</td><td>$tab[3]</td><td>$tab[4]</td><td>$tab[5]</td><td>$tab[6]</td><td>$tab[7]</td><td>$tab[8]</td></tr>";
  }
  close IN;
  $content = $content."</table>";
  $cmd = "sed 's|Title|$title|' $taxonomystat_template | sed 's|Content|$content|'> $output/taxonomy.html";  
  process_cmd ($cmd); 
}