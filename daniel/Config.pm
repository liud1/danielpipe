package Daniel::Config;
use strict;
use warnings;
use threads;

use Exporter qw(import);

our @EXPORT_OK = qw(init_cfg);

sub init_cfg {
  my $config = shift;
  my $type;
  my $key;
  my $value;
  my %cfg;
  open (IN, '< ',$config) or die "Can't read $config file";
  while (<IN>) { 
    my $ln = $_;
    chomp($ln);
    if ($ln=~/\[(.*)\]/) {
      $type = $1;
      #print "$type\n";
    } elsif ($ln=~/(.*)=(.*)/) {
      $key = $1;
      $value = $2;
      $cfg{$type}{$key} = $value;
    }
  }
  close IN;
  return %cfg;
}
