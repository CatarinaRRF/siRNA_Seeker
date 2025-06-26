#!/usr/bin/perl
use strict;
use warnings;

use lib '.';
use Thermodynamics;

# Argumentos: sequência, base, arquivo NN, arquivo Dangling
my ($seq, $base, $nn_file, $dangling_file) = @ARGV;

# Verificar entradas
warn "SEQ: $seq\n";
warn "BASE: $base\n";
warn "NN_FILE: $nn_file\n";
warn "DANGLING_FILE: $dangling_file\n";

my $thermo = Thermodynamics->new();
$thermo->seq($seq);
$thermo->base($base);
$thermo->nearest_neighbor_file($nn_file);
$thermo->dangling_file($dangling_file);

$thermo->cal_energy();

my $energia = $thermo->energy;
warn "Energia calculada (perl): $energia\n";

print $energia, "\n";
