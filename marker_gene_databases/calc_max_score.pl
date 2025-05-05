#! /usr/bin/perl
use strict;
use warnings;


# A quick and dirty script modified from a previously existing script to get the maximum alignment score for a sequence (a self hit)
#
# input: protein fasta file
# input: matrix (choose between BLOSUM45 and BLOSUM62)
# output: tab delimited file with max score and 
#
# Daan Speth, 2021

# requirements

if (scalar @ARGV != 3){
	die "requires 3 arguments:\n\n1) protein fasta file\n2) matrix name (BLOSUM45 or BLOSUM62)\n3) user-defined name for outfile\n\n";
} 

# files

my $fasta = $ARGV[0];
my $matrix = $ARGV[1];
my $out = $ARGV[2];


# read fasta into hash

my %sequences;
my @header_full;
my $header;
open FASTA, $fasta or die "can not open fasta";
while (my $line = <FASTA>){
	chomp $line;
	$line =~ s/\r//;
	my $first_char = substr($line,0,1);
	if ($first_char eq ">"){
		@header_full = split(' ', $line);
		$header = substr($header_full[0],1);
		}
	else {
		$sequences{$header} .= $line;
	}
}
close FASTA;


# print hash id and max_score
# create file and set header
my %max_score;
open OUT, (">> $out") or die "can not create file";
print OUT "Protein_id\tmax_score\n";

# loop over hash and get AA counts
foreach $header(keys(%sequences)) {
	my $A_count = () = $sequences{$header} =~ /\QA/g;
	my $R_count = () = $sequences{$header} =~ /\QR/g;
	my $N_count = () = $sequences{$header} =~ /\QN/g;
	my $D_count = () = $sequences{$header} =~ /\QD/g;
	my $C_count = () = $sequences{$header} =~ /\QC/g;
	my $Q_count = () = $sequences{$header} =~ /\QQ/g;
	my $E_count = () = $sequences{$header} =~ /\QE/g;
	my $G_count = () = $sequences{$header} =~ /\QG/g;
	my $H_count = () = $sequences{$header} =~ /\QH/g;
	my $I_count = () = $sequences{$header} =~ /\QI/g;
	my $L_count = () = $sequences{$header} =~ /\QL/g;
	my $K_count = () = $sequences{$header} =~ /\QK/g;
	my $M_count = () = $sequences{$header} =~ /\QM/g;		
	my $F_count = () = $sequences{$header} =~ /\QF/g;
	my $P_count = () = $sequences{$header} =~ /\QP/g;
	my $S_count = () = $sequences{$header} =~ /\QS/g;
	my $T_count = () = $sequences{$header} =~ /\QT/g;		
	my $W_count = () = $sequences{$header} =~ /\QW/g;
	my $Y_count = () = $sequences{$header} =~ /\QY/g;
	my $V_count = () = $sequences{$header} =~ /\QV/g;
	my $B_count = () = $sequences{$header} =~ /\QB/g;
	my $J_count = () = $sequences{$header} =~ /\QJ/g;
	my $Z_count = () = $sequences{$header} =~ /\QZ/g;
	my $X_count = () = $sequences{$header} =~ /\QX/g;

# within loop, turn AA counts into matrix values and sum
# matrices from here:
# https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM45
# https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM62
	if ($matrix eq "BLOSUM45"){
		$max_score{$header} = $A_count * 5 + $R_count * 7 + $N_count * 6 + $D_count * 7 + $C_count * 12 + $Q_count * 6 + $E_count * 6 + $G_count * 7 + $H_count * 10 + $I_count * 5 + $L_count * 5 + $K_count * 5 + $M_count * 6 + $F_count * 8 + $P_count * 9 + $S_count * 4 + $T_count * 5 + $W_count * 15 + $Y_count * 8 + $V_count * 5 + $B_count * 5 + $J_count * 4 + $Z_count * 5
	}

	if ($matrix eq "BLOSUM62"){
		$max_score{$header} = $A_count * 4 + $R_count * 5 + $N_count * 6 + $D_count * 6 + $C_count * 9 + $Q_count * 5 + $E_count * 5 + $G_count * 6 + $H_count * 8 + $I_count * 4 + $L_count * 4 + $K_count * 5 + $M_count * 5 + $F_count * 6 + $P_count * 7 + $S_count * 4 + $T_count * 5 + $W_count * 11 + $Y_count * 7 + $V_count * 4 + $B_count * 4 + $J_count * 3 + $Z_count * 4
	}
	print OUT "$header\t$max_score{$header}\n";
}
close OUT;
