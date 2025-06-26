#################################################################
# Copyright(c) 2004 Whitehead Institute for Biomedical Research.
#              All Right Reserve
#
# Created:     1/23/2004
# author:      Bingbing Yuan
#################################################################
 
# ============================================================================
# purpose: calculate the thermodynamic energy of the 3'end of dangling strand
# need 2 files: dangling_3_energy.txt
#               nearest_neighbor_energy.txt
# ============================================================================


package Thermodynamics;

use strict;
use Class::Struct;

struct Thermodynamics => {
    seq                   => '$',    # dangling sequence: 7nt from pos 15 to 21 in 5'->3' of a strand
    base                  => '$',    # the dangling base in in the opposite strand
    nearest_neighbor_file => '$',
    dangling_file         => '$',
    energy                => '$',     # float
};


sub cal_energy {

    my $self = shift;

    # convert DNA to RNA; lowercase to uppercase
    my $seq = $self->seq;
    $seq =~ tr/a-z/A-Z/;

    # calculate the nearest_neighbor score
    my $ns_seq = substr($self->seq, 0, 5);
    my $nb_hash_ref = file_to_hash($self->nearest_neighbor_file, 0);
    my $sum_nb = cal_nb($ns_seq, $nb_hash_ref);
   
    # calculate the dangling score
    my $dangling_seq = substr($self->seq, 4, 2);
    my $dangling_hash_ref = file_to_hash($self->dangling_file, 1);
    my $dangling_score = cal_dangling($dangling_seq, $self->base, $dangling_hash_ref);

    $self->energy($sum_nb + $dangling_score);
}


# calculate the total nearest_neighbor value
sub cal_nb {
    my ($seq, $nb_hash_ref) = @_;
    my $score = 0;
    
    my @nts = split(//, $seq);

    for my $i(0..$#nts-1) {
 	my $nb = $nts[$i] . $nts[$i+1];
#	print "pair=$nb, score=", $nb_hash_ref->{$nb}, "<br>"; 
 	$score += $nb_hash_ref->{$nb};
    }
    
    return $score;
}

# calculate the dangling pair value
sub cal_dangling {

    my ($seq, $base, $dangling_hash_ref) = @_;
    my $score = 0;

#    print "base=$base, dangling=$seq, ";
    foreach my $i(0..$#{$dangling_hash_ref->{$seq}} ) {

	if ($base eq $dangling_hash_ref->{$seq}->[$i]) {
	    $score = $dangling_hash_ref->{$seq}->[1];
#	    print "score=$score\n<br>"; 
	    return $score;
	}
    }
#    print "WRONG\n<br>";
    return $score;
}


# ===================================== #
# extract hash from a file
# the key is counted start=0
# only 2 fields: result is simple hash
# more than 3 fields: hash of arrays
# ===================================== #
sub file_to_hash {
    my ($file, $key_pos) = @_;
    my %hash = ();
    
    open (FL, $file) || die "Can not open $file: $!";
    while (<FL>) {
	chomp;
	next if ($_ =~ /^\s*$/);

        my @array = split(/\t/, $_);
	
	for my $i(0..$#array) { 
	    if ($i != $key_pos) {
		
		# only two element: simple hash
		if ($#array == 1) {
		    $hash{$array[$key_pos]} = $array[$i]; 
		}

		# hash of arrays with same order as original file
		else {
		    push @{ $hash{$array[$key_pos]} }, $array[$i];
		}
	    }
	}
    }
    close (FL);	
    
    return \%hash;
    
}

1;
