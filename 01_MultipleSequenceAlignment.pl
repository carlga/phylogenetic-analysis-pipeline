#!/usr/bin/perl
use warnings;
use strict;

my @headers;
my %sequences;

foreach my $fileName ( @ARGV ) {
    my @refs = readFasta($fileName);
    @headers = ( @headers, @{$refs[0]} );
    %sequences = ( %sequences, %{$refs[1]} );
}

#
# initialize tracking arrays
#
my $n = scalar(@headers);
my @alignment = ('') x $n;
my @toalign = (0..$n-1);

#
# find best scoring pair and align
#
my @score_mtx = ( map( [(0)x$n], (0..$n-1) ) );
foreach my $i (0..$n-1) {
    foreach my $j ($i+1..$n-1) {
        my $seq1 = $sequences{$headers[$i]};
        my $seq2 = $sequences{$headers[$j]};
        $score_mtx[$i][$j] = (scoreAlign($seq1, $seq2))[0];
        $score_mtx[$j][$i] = $score_mtx[$i][$j];
    }
}
my ($i, $j) = maxVal(\@score_mtx);
my $seq1 = $sequences{$headers[$i]};
my $seq2 = $sequences{$headers[$j]};
my $trace = (scoreAlign($seq1, $seq2))[1];

# update tracking
($alignment[$i], $alignment[$j]) = seqAlign($seq1, $seq2, $trace);
splice(@toalign, $i, 1);
my $idx = 0;
$idx++ until $toalign[$idx] == $j;
splice(@toalign, $idx, 1);

#
# align remaining sequences by scores
#
while (scalar(@toalign) > 0) {
    my @score_mtx = ( map( [(0)x$n], (0..$n-1) ) );
    my $i = $toalign[0];

    foreach my $j (0..$#alignment) {
        if ($alignment[$j] ne '') {
            $score_mtx[$i][$j] = (scoreAlign($sequences{$headers[$i]}, $alignment[$j]))[0];
        }
    }
    my $j = (maxVal(\@score_mtx))[1];
    my $trace = (scoreAlign($sequences{$headers[$i]}, $alignment[$j]))[1];

    # update tracking
    ($alignment[$i], $alignment[$j]) = seqAlign($sequences{$headers[$i]}, $alignment[$j], $trace);
    shift(@toalign);
}

my $score = 0;
foreach my $i (0..$#alignment) {
    foreach my $j ($i+1..$#alignment) {
        $score += (scoreAlign($alignment[$i], $alignment[$j]))[0];
    }
}

foreach (0..$#headers) {
    print ">$headers[$_]\n$alignment[$_]\n";
}

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

# scores alignment of a pair of sequences
sub scoreAlign {
    my ( $seq1, $seq2 ) = @_;
    my $n1 = length($seq1);
    my $n2 = length($seq2);

    # initialize matrices
    my @dist = ( map( [(0)x($n2+1)], (0..$n1) ) );
    my @trace = ( map( [(0)x($n2+1)], (0..$n1) ) );
    foreach (1..$n1) { $dist[$_][0] = -$_; }
    foreach (1..$n2) { $dist[0][$_] = -$_; }

    # populate with scores and movements
    foreach my $i (1..$n1) {
        foreach my $j (1..$n2) {
            my $m = $dist[$i-1][$j-1] - (substr($seq1,$i-1,1) ne substr($seq2,$j-1,1));
            my $ins = $dist[$i-1][$j] - 1;
            my $del = $dist[$i][$j-1] - 1;
            $dist[$i][$j] = (sort {$a <=> $b} ($m, $ins, $del))[-1];
            if ($m == $dist[$i][$j]) {$trace[$i][$j] = 0;}
            elsif ($ins == $dist[$i][$j]) {$trace[$i][$j] = 1;}
            else {$trace[$i][$j] = 2;}
        }
    }
    my $score = $dist[-1][-1];
    return($score, \@trace);
}

# finds highest value in matrix
sub maxVal {
	my @mtx = @{shift(@_)};
	my $max;
    my $i; my $j;

    foreach my $n (0..$#mtx) {
        foreach my $m (0..$#mtx) {
            next if ( $mtx[$n][$m] == 0);
            if (!defined $max || $max < $mtx[$n][$m]) {
                $max = $mtx[$n][$m];
                $i = $n;
                $j = $m;
            }
        }
    }
    return($i,$j);
}

# aligns a pair of sequences
sub seqAlign {
    my $seq1 = shift;
    my $seq2 = shift;
    my @trace = @{shift(@_)};

    my ($seq1_al, $seq2_al) = ($seq1, $seq2);
    my ($i, $j) = (length($seq1), length($seq2));

    # trace movements from end of each axis
    while ($i>0 && $j>0) {
        if ($trace[$i][$j] == 0) {
            $i -= 1;
            $j -= 1;
        } elsif ($trace[$i][$j] == 1) {
            $i -= 1;
            $seq2_al = substr($seq2_al,0,$j).'-'.substr($seq2_al,$j);
        } else {
            $j -= 1;
            $seq1_al = substr($seq1_al,0,$i).'-'.substr($seq1_al,$i);
        }
    }
    
    # append missing indels
    $seq1_al = ('-' x $j).$seq1_al;
    $seq2_al = ('-' x $i).$seq2_al;
    return($seq1_al, $seq2_al);
}
