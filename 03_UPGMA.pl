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
my @node_dist = ( (0)x$n );
my @node_groups = (0..$n-1);
my $next_node = $n;
my $total_groups = $n-1;    # total groups = n-1

#
# create clusters
#
while ($total_groups > 0) {
	my ($node_a, $node_b)= minVal(\@dist_mtx, \@toblock);

    # update connections in tree
	$tree[$next_node][$node_a] = 1;
	$tree[$next_node][$node_b] = 1;
	$tree[$node_a][$next_node] = 1;
	$tree[$node_b][$next_node] = 1;
	
    # keep track of connected nodes
	push (@toblock, $node_a);
	push (@toblock, $node_b);
	
    # register distance to new node
	$node_dist[$next_node] = $dist_mtx[$node_a][$node_b]/2;
	
    # update grouping for new parent node
	$node_groups[$next_node] = ("$node_groups[$node_a] $node_groups[$node_b]");

    # add parent node to distance matrix
	@dist_mtx = groupNodes(\@dist_mtx, $node_groups[$next_node], $next_node);

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
            $tree[$i][$j] = abs($node_dist[$i]-$node_dist[$j]);
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

# creates parent node to group child nodes
sub groupNodes {
    my @mtx = @{ shift(@_) };
    my @nodes = split(' ', shift);
    my $parent_idx = shift;

    # initialize parent node
	$mtx[$parent_idx][$parent_idx] = 0;
	
    foreach my $i (0..$#mtx-1) {
        my $cum_dist = 0;
        foreach my $j (0..$#nodes) {
            if ($mtx[$nodes[$j]][$i] == 0) {
                $mtx[$parent_idx][$i] = 0;
                $mtx[$i][$parent_idx] = 0;
                last;
            } else {
                $cum_dist += $mtx[$nodes[$j]][$i];
                if ( $j == $#nodes ) {
                    $mtx[$parent_idx][$i] = $cum_dist/scalar(@nodes);
                    $mtx[$i][$parent_idx] = $cum_dist/scalar(@nodes);
                }
            }
        }
    }
    return(@mtx);
}
