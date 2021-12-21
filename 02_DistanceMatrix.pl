#!/usr/bin/perl
use warnings;
use strict;

my @headers;
my %sequences;

my $fileName = shift or die "No file provided\n";
my @data = readFasta($fileName);
@headers = @{ shift @data };
%sequences = %{ shift @data };

#
# initialize and fill distance matrix
#
my $n = scalar(@headers);
my @dist_mtx = ( map( [(0)x$n], (0..$n-1) ) );

foreach my $i (0..$n-1) {
    foreach my $j ($i+1..$n-1) {
        my $seq1 = $sequences{$headers[$i]};
        my $seq2 = $sequences{$headers[$j]};
        my $dist = hammingDist($seq1, $seq2);
        $dist_mtx[$i][$j] = $dist;
        $dist_mtx[$j][$i] = $dist;
    }
}

print "# ".join(" | ", @headers)."\n";
print out2DArray(\@dist_mtx);

exit;


# read in sequences in fasta format
sub readFasta {
    my ( $fileName ) = @_;
    my %sequences;
    my @headers;
    my $header;

    open(F, $fileName) or die "File not found\n";
    while( <F> ) {
        chomp( $_ );
        next if( $_ =~ /^\s*$/ ); # skip empty lines
        if ($_ =~ /^>/ ) {
            $header = $_;
            $header =~ s/>//;
            push( @headers, $header);
            if ( exists $sequences{$header} ) {
                die "Duplicated entries in file: $header\n";
            }
        } else {
            next if( !$header );
            $sequences{$header} .= $_; # concatenate sequence fragments
        }
    }
    close(F);
    return(\@headers, \%sequences);
}

# calculates differences between a pair of sequences
sub hammingDist {
    my ($seq1, $seq2) = @_;
    my $diff;

    if (length($seq1) != length($seq2)) {
        die "Input sequences have different length\n";
    }

    foreach my $i (0..length($seq1)-1) {
        if ((substr($seq1, $i, 1)) ne (substr($seq2, $i, 1))) {
            $diff++;
        }
    }
    return($diff);
}

# outputs 2D array as string
sub out2DArray {
    my @mtx = @{ shift(@_) };
    my $n = scalar(@mtx);
    my $out = '';

    foreach my $i (0..$n-1) {
        foreach my $j (0..$n-1) {
            $out .= "$mtx[$i][$j]\t";
        }
        $out .= "\n";
    }
    return($out);
}
