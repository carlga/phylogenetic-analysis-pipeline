#!/usr/bin/perl
use warnings;
use strict;

my @headers;
my @dist_mtx;

#
# read in ids and distance matrix
#
my $fileName = shift or die "No file provided\n";
open(F, $fileName) or die "File not found\n";
while( <F> ) {
    chomp( $_ );
    next if( $_ =~ /^\s*$/ ); # skip empty lines
    if ( $_ =~ /^#/ ) {
        my $header = $_;
        $header =~ s/# //;
        @headers = split( /\s*\|\s*/, $header );
    } else {
        next if( !@headers );
        push (@dist_mtx, [split(/\s+/, $_)]);
    }
}
close(F);

my $n = scalar(@headers);
if ($n < 3) {
    die "Header error: At least 3 elements required\n";
}
if ($n != scalar(@dist_mtx)) {
    die "Dimentions of headers and distance matrix must match\n";
}

#
# initialize variables
#
my @tree = ( map( [(-1)x(2*$n-1)], (0..2*$n-2) ) );  # total nodes = 2n-1
my @toblock;
my $next_node = $n;
my $total_groups = $n-2;    # total groups = n-2

#
# create clusters
#
while ($total_groups > 0) {
    my @Sx = (0) x scalar(@dist_mtx);
    my @Mij = ( map( [(0) x scalar(@dist_mtx)], (0..$#dist_mtx) ) );

    # Set average distances (Sx)
    foreach my $i (0..$#dist_mtx) {
        my $cum_dist = 0;
        my $block = 0;
		foreach (@toblock) {if ($_ == $i) {$block = 1;}}
		next if ($block == 1);
        foreach my $j (0..$#dist_mtx) {
            next if ( $i == $j);
            $block = 0;
		    foreach (@toblock) {if ($_ == $j) {$block = 1;}}
		    next if ($block == 1);
            $cum_dist += $dist_mtx[$i][$j];
        }
        $Sx[$i]= $cum_dist/$total_groups;
    }
	
    # Set Mij matrix
    foreach my $i (0..$#dist_mtx) {
        my $block = 0;
		foreach (@toblock) {if ($_ == $i) {$block = 1;}}
		next if ($block == 1);
        foreach my $j ($i+1..$#dist_mtx) {
            $block = 0;
		    foreach (@toblock) {if ($_ == $j) {$block = 1;}}
		    next if ($block == 1);
            $Mij[$i][$j] = $dist_mtx[$i][$j]-$Sx[$i]-$Sx[$j];
			$Mij[$j][$i] = $Mij[$i][$j];
        }
    }

    # find minimal value
	my ($node_a, $node_b)= minVal(\@Mij, \@toblock);

    # update distances in tree
	$tree[$next_node][$node_a] = ($dist_mtx[$node_a][$node_b]+$Sx[$node_a]-$Sx[$node_b])/2;;
	$tree[$next_node][$node_b] = ($dist_mtx[$node_a][$node_b]+$Sx[$node_b]-$Sx[$node_a])/2;
	$tree[$node_a][$next_node] = $tree[$next_node][$node_a];
	$tree[$node_b][$next_node] = $tree[$next_node][$node_b];
	
    # keep track of connected nodes
	push (@toblock, $node_a);
	push (@toblock, $node_b);
	
    # register distance to new node
    foreach my $j (0..2*$n-2-$total_groups-1) {
        my $block = 0;
		foreach (@toblock) {if ($_ == $j) {$block = 1;}}
		next if ($block == 1);
        $dist_mtx[$next_node][$j] = ($dist_mtx[$node_a][$j]+$dist_mtx[$node_b][$j]-$dist_mtx[$node_a][$node_b])/2;
		$dist_mtx[$j][$next_node] = $dist_mtx[$next_node][$j];
		if ($total_groups == 1){
			$tree[$next_node][$j] = $dist_mtx[$next_node][$j];
			$tree[$j][$next_node ]= $dist_mtx[$next_node][$j];
        }
    }
	$next_node++;
    $total_groups--;
}

#
# output distances between nodes
#
print "# ".join(" | ", @headers)."\n";
foreach my $i (0..$#tree) {
    foreach my $j (0..$#tree) {
        if ($tree[$i][$j] != -1) {
            printf "%d->%d:%0.3f\n", $i, $j, $tree[$i][$j];
        }
    }
}

exit;


# finds lowest value in distance matrix
sub minVal {
    my @mtx = @{ shift(@_) };
    my @toblock = @{ shift(@_) };
	my $min;
    my $i; my $j;
	
    foreach my $n (0..$#mtx) {
        my $block = 0;
		foreach (@toblock) {if ($_ == $n) {$block = 1;}}
		next if ($block == 1);
        foreach my $m ($n+1..$#mtx) {
            next if ( $mtx[$n][$m] == 0);
            $block = 0;
		    foreach (@toblock) {if ($_ == $m) {$block = 1;}}
		    next if ($block == 1);
            if (!defined $min || $min > $mtx[$n][$m]) {
                $min = $mtx[$n][$m];
                $i = $n;
                $j = $m;
            }
        }
    }
	return ($i, $j);
}
