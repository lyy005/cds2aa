#!/usr/bin/perl
use strict;

die "perl $0 [cds fasta file] [Nuclear codon table (NUC) | Invertebrate mitochondrial codon table (MT)] [output file]\n Translate based on phase=0, no gap will be removed\n" unless (@ARGV == 3);
my %CODE_NUC = ( # Nuclear condon table !!! Not invertebrate condon table
                                'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',                               # Alanine
                                'TGC' => 'C', 'TGT' => 'C',                                                           # Cysteine
                                'GAC' => 'D', 'GAT' => 'D',                                                           # Aspartic Acid
                                'GAA' => 'E', 'GAG' => 'E',                                                           # Glutamic Acid
                                'TTC' => 'F', 'TTT' => 'F',                                                           # Phenylalanine
                                'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',                               # Glycine
                                'CAC' => 'H', 'CAT' => 'H',                                                           # Histidine
                                'ATC' => 'I', 'ATT' => 'I', 'ATA' => 'I',                                           # Isoleucine
                                'AAA' => 'K', 'AAG' => 'K',                                                           # Lysine
                                'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', 'TTA' => 'L', 'TTG' => 'L',   # Leucine
                                'ATG' => 'M',                                                                         # Methionine
                                'AAC' => 'N', 'AAT' => 'N',                                                           # Asparagine
                                'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',                               # Proline
                                'CAA' => 'Q', 'CAG' => 'Q',                                                           # Glutamine
                                'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 'AGG' => 'R', 'AGA' => 'R', 'AGG' => 'R',   # Arginine
                                'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 'AGC' => 'S', 'AGT' => 'S',   # Serine
                                'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',                               # Threonine
                                'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',                               # Valine
                                'TGG' => 'W',                                                                         # Tryptophan
                                'TAC' => 'Y', 'TAT' => 'Y',                                                           # Tyrosine
                                'TAA' => 'U', 'TAG' => 'U', 'TGA' => 'U',                                             # Stop
				'---' => '-'
        );

my %CODE_MT = (
                                'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',                               # Alanine
                                'TGC' => 'C', 'TGT' => 'C',                                                           # Cysteine
                                'GAC' => 'D', 'GAT' => 'D',                                                           # Aspartic Acid
                                'GAA' => 'E', 'GAG' => 'E',                                                           # Glutamic Acid
                                'TTC' => 'F', 'TTT' => 'F',                                                           # Phenylalanine
                                'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',                               # Glycine
                                'CAC' => 'H', 'CAT' => 'H',                                                           # Histidine
                                'ATC' => 'I', 'ATT' => 'I',                                                           # Isoleucine
                                'AAA' => 'K', 'AAG' => 'K',                                                           # Lysine
                                'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', 'TTA' => 'L', 'TTG' => 'L',   # Leucine
                                'ATG' => 'M', 'ATA' => 'M',                                                           # Methionine
                                'AAC' => 'N', 'AAT' => 'N',                                                           # Asparagine
                                'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',                               # Proline
                                'CAA' => 'Q', 'CAG' => 'Q',                                                           # Glutamine
                                'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R',                               # Arginine
                                'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 'AGC' => 'S', 'AGT' => 'S', 'AGA' => 'S', 'AGG' => 'S', # Serine
                                'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',                               # Threonine
                                'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',                               # Valine
                                'TGG' => 'W', 'TGA' => 'W',                                                           # Tryptophan
                                'TAC' => 'Y', 'TAT' => 'Y',                                                           # Tyrosine
                                'TAA' => 'U', 'TAG' => 'U',                                                           # Stop
				'---' => '-'
        );

open (FA, $ARGV[0]) or die "$ARGV[0] $!\n";
open (OT, ">$ARGV[2]") or die "$ARGV[2] $!\n";

$/=">"; 
<FA>;
while (<FA>) {
	chomp;
        my @line = split/\n+/;
	my $name = shift @line;
	my $seq = join "", @line;	
	$seq =~ tr/atcg/ATCG/;
	
	my $len = length ($seq);
	my $aa;
	for(my $c=0; $c<$len; $c+=3){
		my $codon = substr($seq, $c, 3);
		
		if($ARGV[1] eq "NUC"){
			if($CODE_NUC{$codon}){
				$aa .= $CODE_NUC{$codon};
			}else{
				$aa .= "X";
				print "Unknown codon is found: $codon\n";
			}
		}elsif($ARGV[1] eq "MT"){
			if($CODE_MT{$codon}){
				$aa .= $CODE_MT{$codon};
			}else{
				$aa .= "X";
				print "Unknown codon is found: $codon\n";
			}
		}else{
			print "Input code error: $ARGV[1]\n";
		}
	}
	print OT ">$name\n$aa\n";
}
close IN;
print "Done!\n";
