#!/usr/bin/perl
use strict;
use Cwd;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;

# Script for making SNP substitutions into fasta files. Requires fasta format for individual chrmosomes to be named chr?.fa. Will use last line of fasta for sequence
# Use bedtools getfasta -fi [GENOME_FA] -bed [CHR_BED] -fo [OUT.FA] with beds for each chromosome (or all chromosomes in order for next command)
# Use "head -n 4 | tail -n 2 > chr2.fa" to split directly from fasta


&GetOptions
(
    "bed=s"=>\ my $input,   # --bed             chr,start,stop,rsID,ref,var  (Tab/space separated)           
    "name=s"=>\ my $name,   # --name            Output name
);


open (INPUT, $input);

my $old_chr= "NIL";
my $counter = 0;
my $fasta;
my $last_start;

open (OUTPUT, ">$name\_Modified_chromosomes.fa");
open (REPORT, ">$name\_Report.bed");
print REPORT "chr\tstart\tstop\trsID\tref\tvar\tnew\n";
open (ERROR, ">$name\_Error_list.bed");
print ERROR "chr\tstart\tstop\trsID\tref\tvar\tbase_found\n";

while (my $line = <INPUT>)
            {
            chomp $line;
            my ($chr, $start, $stop,$rsID, $ref, $var) = split(' ',$line);
            my $ref_lc = lc($ref);
            my $var_lc = lc($var);         
            # Convert to lowercase to match
            if ($chr ne $old_chr)
                {
                        if ($counter == 0)      #Put in header of first chromosome.
                            {print OUTPUT ">$chr\n";}
                        if ($counter != 0)      #Print the end of the old chr, and header for new chromosome.
                            {my $rest = substr $fasta, $last_start;print OUTPUT "$rest\n>$chr\n";}
                        $counter++;
                        close FASTA;
                        open (FASTA, "$chr.fa");        #Chromosome fasta file (2 lines).                           
                        while (my $line2 = <FASTA>)
                            {
                                    chomp $line2;
                                    $fasta = $line2;
                            }
                        my $before_length = $start-1;
                        my $base_coord = $start-1;
                        my $before = substr $fasta, 0, $before_length;# Before: 0, start-2; (adjust to 0-based and 1 before)
                        print OUTPUT "$before";
                        my $base = substr $fasta, $base_coord, 1;                 # Base: start-1,1: (adjust to 0-based take 1 digit)
                        my $base_lc = lc($base);
                        if($base_lc eq $ref_lc) {print OUTPUT $var_lc;print REPORT "$line\tVAR\n";}
                        if($base_lc eq $var_lc) {print OUTPUT $ref_lc;print REPORT "$line\tREF\n";}
                        if($base_lc ne $ref_lc && $base_lc ne $var_lc){print OUTPUT $base_lc;print ERROR "$line\t$base_lc\n";print REPORT "$line\tUNMATCHED\n";}
                        $old_chr=$chr;
                        $last_start = $start;
                        next;                    
                }
            if ($chr eq $old_chr)
                {
                    my $space_length = $start - $last_start -1;
                    my $space_seq = substr $fasta, $last_start, $space_length;
                    print OUTPUT "$space_seq";
                    my $base_coord = $start-1;
                    my $base = substr $fasta, $base_coord, 1;
                    my $base_lc = lc($base);  
                    if($base_lc eq $ref_lc) {print OUTPUT $var_lc;print REPORT "$line\tVAR\n";}
                    if ($base_lc eq $var_lc) {print OUTPUT $ref_lc;print REPORT "$line\tREF\n";}
                    if($base_lc ne $ref_lc && $base_lc ne $var_lc){print OUTPUT $base_lc;print ERROR "$line\t$base_lc\n";print REPORT "$line\tUNMATCHED\n";}
                    $old_chr=$chr;
                    $last_start = $start;                      
                }
            }
my $rest = substr $fasta, $last_start;
print OUTPUT "$rest\n";
close OUTPUT;
close INPUT;

exit;

