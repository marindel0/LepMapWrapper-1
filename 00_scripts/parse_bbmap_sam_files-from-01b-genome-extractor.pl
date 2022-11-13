#!/usr/bin/perl -w

use strict;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use Cwd qw(getcwd);


#########################################################################################################
#
# BBmap produces rather confusing SAM files.  This is an attempt to tidy those up and extract relevant
# information
# 
# Written by Zophonías O. Jónsson - 7. Nov. 2022
#########################################################################################################

### TODO fasta deflines from 01b_genome_extractor, 02_catalog_extractor and 02b_denovo_catalog extrctor
### are not identical.  Should ask or recognize them and deal with all versions.

#########################################################################################################
# Set parameters - initialize variables
#########################################################################################################

my @arguments = @ARGV; #save arguments as a precaution

### The get options block

my ($help, $samfile, $outfile, $outpath );

#Let's be friendly
sub Usage {
  print "Unknown option: @_\n" if ( @_ );
  print "\nThis is how you use this script: $0

An imput file must be provided and optionally an outfilename and outfilepath

-s/--samfile <file>\t:\tThe sam file from bbmap. \n\n";

print "Optional parameters:

-o/--outfile <filename>\t:\tFilname for outputfile (if not provided name of samfile will be used as template).

-p/--path-out <path>\t:\tPath to write output files to (defaults to current directory).

-h/--help\t\t:\tPrint this message.

";
  exit;
}

GetOptions('h|help' => \$help, 's|samfile=s' => \$samfile, 'o|outfile:s' => \$outfile, 'p|path-out:s' => \$outpath,) or die ("Arguments in error!\n");
Usage() if (@arguments == 0 || $help || not ($samfile));

(my $strippedbasename = basename($samfile)) =~ s/\.[^.]+$//; 
#(my $outfilebasename = $strippedbasename) =~ s/^contig_order_//;
unless (defined($outfile))  { $outfile = $strippedbasename . ".tsv" };

if (!defined($outpath)){
    $outpath=getcwd();  
}
elsif (-d $outpath) {
    print "Outputfiles will be written to $outpath\n";
    $outpath =~ s/\/$//;    # remove trailing slash if there is one and add back later (may be uneccesary)
}
elsif (-e $outpath) {
    die "$outpath is a file, not a directory";
}
else { 
    die "Couldn't find directory $outpath: $!";
} 


print "\nInput SAM: $samfile \n";
print "Output filename: $outfile \n";
print "Output path: $outpath\n";

my %datahash;


sub parsethefile(){

    my $text ="";
    my $sampleid ="";
    my $samplename="";
    
    
    print "Reading file: $samfile\n";
    open (INFILE, '<', $samfile) or die "can't find file $samfile: $!";
    while (<INFILE>){
        chomp;
        s/^\s*//;    #strip leading whitespace
        s/\s+$//;    #strip trailing whitespace
        
        if ($_ =~ m/^@/){      #lose all the comment lines
            next;
        }
        
### defline from 01b_genome_extractor: >1:LG=LG1:female_pos=0:marker=224:region=ptg000072l:283487-283737:marker_pos=283612
### defline from 02b_denovo_catalog_extractor: >1:LG=LG1:female_pos=0:marker=1301:contig=158754:pos=181

## this version is for the 01b_case !!!

        my @linepieces=split("\t", $_);
        my ($defline, $match, $map_chrom, $map_pos, $map_qual, $cigar, $star, $null, $nill, $seq) = @linepieces;
        my @deflinepieces=split(":", $defline);
        my $map_order=$deflinepieces[0];
        $datahash{$map_order}{'LG'}=(split("=", $deflinepieces[1]))[1];
        $datahash{$map_order}{'Female_pos'}=(split("=", $deflinepieces[2]))[1];
        $datahash{$map_order}{'marker'}=(split("=", $deflinepieces[3]))[1];
        $datahash{$map_order}{'contig'}=(split("=", $deflinepieces[4]))[1];
        $datahash{$map_order}{'contig_region'}=($deflinepieces[5]);
        $datahash{$map_order}{'marker_pos'}=(split("=", $deflinepieces[6]))[1];
        $datahash{$map_order}{'match'}=$match;
        my @map_chrompieces=split(" ", $map_chrom);
        $datahash{$map_order}{'Map_chrom'}=$map_chrompieces[0];
        $datahash{$map_order}{'Map_chrom_pos'}=$map_pos;
        $datahash{$map_order}{'Map_qual'}=$map_qual;
        $datahash{$map_order}{'Cigar'}=$cigar;
        $datahash{$map_order}{'seq'}=$seq;        
        next;
    }
    close(INFILE);
}

sub writetheoutputfile (){

	print "writing output data to $outfile\n";
	unless (open (OUTFILE, ">", "$outpath/$outfile")){
		die "Can't write to $outfile: $!";
		}
    
    print OUTFILE "Map_order\tLG\tFemale_pos\tMarker\tContig\tContig_region\tMarker_pos\tMatch_qual\tMap_CHR\tMap_CHR_pos\tMap_qual\tCigar\tSequence\n";
	foreach my $key (sort {$a <=> $b} keys %datahash) {
        print OUTFILE "$key\t$datahash{$key}{'LG'}\t";
        print OUTFILE "$datahash{$key}{'Female_pos'}\t";
        print OUTFILE "$datahash{$key}{'marker'}\t";
        print OUTFILE "$datahash{$key}{'contig'}\t";
        print OUTFILE "$datahash{$key}{'contig_region'}\t";
        print OUTFILE "$datahash{$key}{'marker_pos'}\t";
        print OUTFILE "$datahash{$key}{'match'}\t";
        print OUTFILE "$datahash{$key}{'Map_chrom'}\t";
        print OUTFILE "$datahash{$key}{'Map_chrom_pos'}\t";
        print OUTFILE "$datahash{$key}{'Map_qual'}\t";
        print OUTFILE "$datahash{$key}{'Cigar'}\t";
        print OUTFILE "$datahash{$key}{'seq'}\n";
    }
	close(OUTFILE);
	print "Successfully wrote analysis to $outfile\n";
}


sub write2stdout(){
	foreach my $key (sort {$a <=> $b} keys %datahash) {
        print "$key\t$datahash{$key}{'LG'}\t";
        print "$datahash{$key}{'Female_pos'}\t";
        print "$datahash{$key}{'marker'}\t";
        print "$datahash{$key}{'stack'}\t";
        print "$datahash{$key}{'stack_pos'}\t";
        print "$datahash{$key}{'match'}\t";
        print "$datahash{$key}{'Map_chrom'}\t";
        print "$datahash{$key}{'Map_chrom_pos'}\t";
        print "$datahash{$key}{'Map_qual'}\t";
        print "$datahash{$key}{'Cigar'}\t";
        print "$datahash{$key}{'seq'}\n";
    }

}

parsethefile();
#write2stdout();
writetheoutputfile();
print "done";

