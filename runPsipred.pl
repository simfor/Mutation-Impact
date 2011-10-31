use File::Temp qw/ tempfile tempdir /;

my $basename = $ARGV[0];
$basename =~ s/^.*\/(.*)\..*/$1/;

($fh, $tmpfile) = tempfile();

#Runs psiblast
`psiblast -db $db_path -query $ARGV[0] -inclusion_ethresh 0.001 -out_pssm $tmpfile.chk -num_iterations 3 -num_alignments 0 > $tmpfile.blast`;

#Psipred pipeline
`chkparse $tmpfile.chk > $tmpfile.mtx`;
`psipred $tmpfile.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 > $tmpfile.ss`;
`psipass2 $datadir/weights_p2.dat 1 1.0 1.0 $tmpfile.ss2 $tmpfile.ss > $basename.horiz`;
