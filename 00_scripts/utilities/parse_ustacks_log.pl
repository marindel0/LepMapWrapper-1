#!/usr/bin/perl -w

use strict;
use File::Basename;


#########################################################################################################
#
# This is a little script to parse logged output from ustacks to generate a tabbed report of
# the stack building process for each sample
# 
# Written by Zophonías O. Jónsson - July 3rd 2020
#########################################################################################################

my @arguments = @ARGV;
my $outfilename = 'parsed_ustackslog.tsv';
my %datahash;

my $filename = shift(@arguments);

sub parsethefile(){

    my $text ="";
    my $sampleid ="";
    my $samplename="";
    
    
    print "Reading file: $filename\n";
    open (INFILE, $filename) or die "can't find file $filename: $!";
    while (<INFILE>){
        chomp;
        s/^\s*//;
        
        if ($_ =~ m/Input file:/){
            ($text, $samplename) = split(": ");
            ($samplename=basename($samplename)) =~ s/\.1.fq.gz//;
            next;
        }
        if ($_ =~ m/Sample ID:/){
            #print "$_\n";
            ($text, $sampleid) = split(": ");
            #print "Sample ID: $value\n";
            $datahash{$sampleid}{'sampleid'} = $sampleid;
            $datahash{$sampleid}{'samplename'} = $samplename;
            next;
        }
        if ($_ =~ m/^Loaded /){
            my @fields = split(" ");
            my $reads = $fields[1];
            $datahash{$sampleid}{'loadedreads'} = $reads;
            next;
        }
        if ($_ =~ m/^Final coverage:/){
            my @fields = split(" ");
            my $meancoverage = $fields[2];
            my ( $number ) = ( $meancoverage =~ /(\d+(?:\.\d+)?)/ );
            $datahash{$sampleid}{'meancoverage'} = $number;
            my $stdev = $fields[3];
            my ( $number ) = ( $stdev =~ /(\d+(?:\.\d+)?)/ );
            $datahash{$sampleid}{'stdev'} = $number;
            my $max = $fields[4];
            my ( $number ) = ( $max =~ /(\d+(?:\.\d+)?)/ );
            $datahash{$sampleid}{'max'} = $number;
            my $n_reads = $fields[5];
            my ( $number ) = ( $n_reads =~ /(\d+(?:\.\d+)?)/ );
            $datahash{$sampleid}{'n_reads'} = $number;
        }
        next;
  
  
    }

}

sub writetheoutputfile (){

	print "writing output data to $outfilename\n";
	unless (open (OUTFILE, ">$outfilename")){
		die "Can't write to $outfilename: $!";
		}
    
    print OUTFILE "Sample ID\tSamplename\tLoaded reads\tMean coverage\tStandard deviation\tMax\tFinal nr of reads\n";
	foreach my $key (sort {$a <=> $b} keys %datahash) {
        print OUTFILE "$key\t$datahash{$key}{'samplename'}\t";
        print OUTFILE "$datahash{$key}{'loadedreads'}\t";
        print OUTFILE "$datahash{$key}{'meancoverage'}\t";
        print OUTFILE "$datahash{$key}{'stdev'}\t";
        print OUTFILE "$datahash{$key}{'max'}\t";
        print OUTFILE "$datahash{$key}{'n_reads'}\n";
    }
	

	close (OUTFILE);
	print "Successfully wrote analysis to $outfilename\n";
}


sub write2stdout(){
    foreach my $key (sort {$a <=> $b} keys %datahash) {
        print "$key => $datahash{$key}{'samplename'}\n";
        print "Loaded reads $datahash{$key}{'loadedreads'}\n";
        print "Mean coverage $datahash{$key}{'meancoverage'}\n";
        print "Standard dev $datahash{$key}{'stdev'}\n";
        print "Max $datahash{$key}{'max'}\n";
        print "Final n reads $datahash{$key}{'n_reads'}\n\n";
    }

}

parsethefile();
write2stdout();
writetheoutputfile();
print "done";

