use strict;
use warnings;
use XML::Simple;

my $xml = new XML::Simple;

my $data = $xml->XMLin($blast);

for(@{$data->{BlastOutput_iterations}->{Iteration}->{Iteration_hits}->{Hit}}){
	print OUT ">$_->{Hit_id}\n";
	print OUT "$_->{Hit_hsps}->{Hsp}->{Hsp_hseq}\n";
}
