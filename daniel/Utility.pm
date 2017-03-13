package Daniel::Utility;
use strict;
use warnings;
use threads;
 
use Exporter qw(import);
our @EXPORT_OK = qw(process_cmds_parallel process_cmds_serial process_cmd splitfasta create_fullshell checkjob extract_seq_by_id checkjobid);
 
$ENV{PATH} = "$ENV{PATH}:/opt/sge/bin/lx-amd64/"; 

sub getCurrentDateTime  {
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
  my $nice_timestamp = sprintf ( "%04d-%02d-%02d %02d:%02d:%02d",
  $year+1900,$mon+1,$mday,$hour,$min,$sec);
  return "[$nice_timestamp]";
}
sub printlog {
  my $msg  = shift;
  open (OUT, ">> my.log") or die "Can't create my.log file \n";
  print OUT "$msg\n";
  close OUT;
}
sub process_cmds_parallel {
  my @cmds = @_;
  my @threads;
  my $msg;
  foreach my $cmd (@cmds) { 
    my $thread = threads->create('process_cmd', $cmd);
    push (@threads, $thread);
  }
                
  my $ret = 0;   
  foreach my $thread (@threads) {
    $thread->join();
    if (my $error = $thread->error()) {
      $msg = "ERROR, thread exited with error $error\n";
      printlog($msg);
      $ret++;
    }
  }
  if ($ret) {
    $msg = "ERROR, $ret threads errored out";
    printlog($msg);
  }
  return;
}
sub process_cmds_serial {
  my @cmds = @_;
  foreach my $cmd (@cmds) {
    process_cmd($cmd);
  }
  return;
}
sub process_cmd {
  my ($cmd) = @_;
  my $start_time = getCurrentDateTime();
  my $ret = system("bash", "-c", $cmd);
  my $end_time = getCurrentDateTime();
  my $msg = "CMD FINISHED: '$cmd' $start_time - $end_time \n";
  if ($ret) {
    $msg = "CMD ERROR: '$cmd' died with ret $ret";
    printlog($msg);
    return;
  }    
  printlog($msg);
  return;
}
sub splitfasta {
  my $fa = shift;     # big fasta file
  my $output = shift; # the output location of tidy fasta files
  my $prefix = shift; # the prefix of header of tidy fasta files
  my $count = 1;# the count variable
  my $ln;     # the content of each line
  my $header =""; # the header of the tidy fasta
  my $seq = "";# the sequence of the tidy fasta
  my $sentinel = "N";
  my $outfn;     # file name of tidy fasta
  unless (-e $output) {system ("mkdir -p $output");}
  open (IN, '<', $fa) or die "Can't read $fa file\n";
  while (<IN>) {
    $ln = $_; chomp($ln);
	if ($ln=~/>(.*)/) {	  
	  if ($sentinel eq "Y") {	   
	    open (OUT, '>', $outfn) or die "Can't write $outfn file\n";
		print OUT ">$header\n$seq\n";
	    close OUT;				
	  }
	  $header = $1;
	  $outfn = "$output/$prefix\_$count.fa";	
	  $sentinel = "Y";
	  $count = $count + 1;	  
	  $seq = "";
	} else {
	  $seq = $seq.$ln."\n";
	}
  }
  close IN;
  unless (-e $outfn) {
    open (OUT, '>', $outfn) or die "Can't write $outfn file\n";
    print OUT ">$header\n$seq\n";
    close OUT;
  }
}
sub create_fullshell {
# create_fullshell (ShellScript, QueueName, ErrFile, LogFile, Cmd);
  my $shell = shift;
  my $qname = shift;
  my $err = shift;
  my $log = shift;
  my $cmd = shift;  
  open (OUT, '>',$shell) or die "Can't create $shell file";
  print OUT "#!/bin/sh\n";
  print OUT "#\$ -N $qname\n";
  print OUT "#\$ -e $err\n";
  print OUT "#\$ -o $log \n";
  print OUT "$cmd\n";
  close OUT;
}
sub checkjob {
  # This needs to more special
  # I will check the jobid for the job completely
  # use `qsub $shell` instead to write down jobid in a temp file and new subroute will check jobid exists
  my $qname = shift; # Queue Name
  my $cmd;
  my $value = 1;
  sleep 60;
  while(1) {    
    $cmd = `qstat -f | grep $qname`;
	if ($cmd eq "") {$value = 0;}
    if ($value eq 0) {last;} else {sleep 60;}
  }   
}

sub checkjobid {
  # This is more special
  # I will check the jobid for the job completely
  # use `qsub $shell` instead to write down jobid in a temp file and new subroute will check jobid exists
  my $qname = shift; # Queue Name
  my $cmd;
  my $value = 1;
  my $ln;
  my @jobs;
  my $string;
  sleep 60;
  if (-e "./jobs") {
    open (IN, '<', "./jobs") or die "Can't read jobs file\n";;
	while (<IN>) {
	  $ln = $_; chomp($ln);
	  if ($ln=~/Your job (\d+) \((.*)\) has been submitted/) {
	    push (@jobs, $1);
	  }
	}
	close IN;
	$string =  join("|", @jobs);
    while(1) {    	  
      $cmd = `qstat -f | grep -E \'$string\'`;
	  if ($cmd eq "") {$value = 0;}
      if ($value eq 0) {last;} else {sleep 60;}
    } 		
  } else { 
    while(1) {    	  
      $cmd = `qstat -f | grep $qname`;
	  if ($cmd eq "") {$value = 0;}
      if ($value eq 0) {last;} else {sleep 60;}
    }   
  }
}

sub extract_seq_by_id {
  my $list = shift; # list of sequence id, the word behind >
  my $fa = shift; # the multiple fasta file
  my $ln;
  my $data;
  open (IN, '<', $list) or die "Can't read $list file\n";
  while (<IN>) {
    $ln= $_;chomp($ln);
    $data->{$ln} = $ln;
  }
  close IN;

  my $sentinel = "N";
  open (IN, '<', $fa) or die "Can't read $fa file\n";
  while (<IN>) {
    $ln = $_;chomp($ln);
    if ($ln=~/>(.*)/) {
      $sentinel = "N";
      if (defined $data->{$1}) {$sentinel = "Y";print "$ln\n";}
    } else {
      if ($sentinel eq "Y") {print "$ln\n";}
    }
  }
  close IN;
}

1;
