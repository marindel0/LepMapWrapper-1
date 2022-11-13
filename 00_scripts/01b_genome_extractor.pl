#!/usr/bin/perl 

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use Cwd qw(getcwd);


#########################################################################################################
#
# This little script parses the markerlist output files from LepMap3 given the genome used to 
# call the SNPs and retrieves 250 bp of surrounding sequence
# Meant to replace the R and .py scripts in 01_Extract_loci
#
# Output formats: various types of lists ... For example
# 
# No header
# Species: <string>
# Family: <number>
# Female pos: <float>
# Phenotype: <number?> Not sure about this column - have to look it up in documentation or code.
# Sequence: <string>
# 
# 
# Written by Zophonías O. Jónsson - Oct 25th 2022
#########################################################################################################

#########################################################################################################
# Set parameters - initialize variables
#########################################################################################################

my @arguments = @ARGV; #save arguments as a precaution

### The get options block

my ($help, $markerlist_tsv, $genome_fasta, $outpath, $specieslabel);

#Let's be friendly
sub Usage {
  print "Unknown option: @_\n" if ( @_ );
  print "\nThis is how you use this script: $0

Two input files must be provided with the correct option switches

-m/--markerlist <file>\t:\tThe markerlist map file (from step 8).

-g/--genome <file>\t:\tThe catalog file (gzipped)\n\n";

print "Optional parameters:

-o/--outpath <path>\t:\tPath to write output files to (defaults to current directory)

-l/--label\t:\tA short label included at the beginning of each line in mapcomp files.

-h/--help\t\t:\tPrint this message.

";
  exit;
}

GetOptions('h|help' => \$help, 'm|markerlist=s' => \$markerlist_tsv, 'g|genome=s' => \$genome_fasta, 'o|outpath:s' => \$outpath, 'l|label:s' => \$specieslabel) or die ("Arguments in error!\n");
Usage() if (@arguments == 0 || $help || not ($markerlist_tsv) || not ($genome_fasta));


# TODO provide option to enter select whick types of output files to write ... or maybe not.

#(print $markerlist_tsv . " " . $genome_fasta) or Usage();
### Set global variables
my %markerlist_hash;  
my %catalog_hash;     
my %catalog_cache; 

#my $markerlist_tsv = shift(@arguments) or die "no markerfile $1";
#my $genome_fasta = shift(@arguments)or die "no catalog file $1"; 

(my $strippedbasename = basename($markerlist_tsv)) =~ s/\.[^.]+$//; 
(my $outfilebasename = $strippedbasename) =~ s/^contig_order_//;
#my $workingdir=getcwd();
my $counter=0;
my @fileheader;   #This should be removed and made local TODO


print "\nMarkerlist order file: $markerlist_tsv \n";
print "Genome: $genome_fasta \n";

if (defined($outpath)){
    chomp($outpath);
    if (-d $outpath) {
        print "Outputfiles will be written to $outpath\n";
        $outpath =~ s/\/$//;    # remove trailing slash if there is one and add back later (may be uneccesary)
    }
    elsif (-e $outpath) {
        die "$outpath is a file, not a directory";
    }
    else { 
        die "Couldn't find directory $outpath: $!";
    } 
}
else{
    $outpath=getcwd();
}

### TODO Check for number of arguments passed.

print "Output filename(s): $outfilebasename + extension \n";

sub get_specieslabel {
    print "\n\nMapcomp files have a \"SpeciesName\" identifier field. What do you want to use ?\n";
    print "Do not use commas and avoid spaces or other special characters.\n";
    print "Unless you enter something <" . basename($markerlist_tsv) . "> will be used.\n";
    print "Enter a name: ";
    $specieslabel = <STDIN>;
    chomp($specieslabel);
    $specieslabel ||= $strippedbasename;
}

unless (defined($specieslabel))  { get_specieslabel() };


#########################################################################################################
# Parse the markerlist from LepMap3 and stuff it into a hash
# 
# Passing filname as an argument may seem uneccessarily complicated but I want this to be reusable
#########################################################################################################


sub parse_markerlist_file($){

   my $filename = shift;    #get the filename passed to sub 

   # The format of the markerlist is: (Header + 1st line) -> 7 columns, tab separated 
   # No empty lines 
   #
   # CHR	marker_id	male_pos	female_pos	contig	pos	contig_pos    
   # LG1	1738	84.32	10.668	NC_036874.1	7219174	NC_036874.1_7219174
   #
   # Location: LepMapWrap dir 08_analyze_maps/01_maps
   
   my $text ="";
   my $sampleid ="";
   my $samplename="";
   
   print "Reading markerlist file: $filename\n";
   open (MARKERLISTFILE, $filename) or die "can't find file $filename: $!";
   
   # TODO Check for correct header (now we will just print it out )
   my $firstline = <MARKERLISTFILE>;
   chomp($firstline);
   @fileheader = split /\t+/, $firstline;
   my ($CHR, $marker_id, $male_pos, $female_pos, $contig, $pos, $contig_pos) = @fileheader;
 
   foreach (@fileheader){
        print "Column " . ++$counter . ": $_\n";
   }
   
   my $linecounter=0;
   
   while (<MARKERLISTFILE>){
       chomp;
       s/^\s*|\s+$//g;  #remove leading and trailing whitespace - an old habbit
 
       my @column = split /\t+/, $_;
       
       #We use contig:pos as key becaust that's what binds the markerlist and loci
       my $CHROM = $column[4];
       $CHROM  =~ s/(_)(\d)$/.$2/ ;
       my $POS = $column[5];
       my $marker= join ('_', $CHROM, $POS);
       $markerlist_hash{$marker}{'CHR'} = $column[0];
       $markerlist_hash{$marker}{'marker_id'} = $column[1];  
       $markerlist_hash{$marker}{'male_pos'} = $column[2];
       $markerlist_hash{$marker}{'female_pos'} = $column[3];
       $markerlist_hash{$marker}{'contig'} = $CHROM;
       $markerlist_hash{$marker}{'pos'} = $column[5];
       $markerlist_hash{$marker}{'contig_pos'} = $column[6];  #redundant but let's keep it
       $markerlist_hash{$marker}{'line_number'} = ++$linecounter;
       next;
                       
   }                   
   close MARKERLISTFILE ;
}

parse_markerlist_file($markerlist_tsv);

sub add_sequences {

    print "Fetching 250 bp of sequence surrounding marker on chromosome.\n";
    foreach my $marker (keys %markerlist_hash) {
        my $contig = $markerlist_hash{$marker}{'contig'};
        my $from = ($markerlist_hash{$marker}{'pos'} - 125);
        my $to = ($markerlist_hash{$marker}{'pos'} + 125);
        if ($from <=0) {   #if marker is close to beginning of contig
            $to += (abs($from) +1);
            $from = 1;
        }
        my $command = "samtools faidx $genome_fasta $contig:$from-$to";
        #print "Commandline: $command\n"; 
        print "*";
        
        sub handle_exit_gracefully(\$) {
            my $com = shift;
            die "Something went wrong in when executing $com" ;
            # Not so graceful after all but can be improved if needed
        }
        
        open(INPIPE, "$command |") or handle_exit_gracefully($command) ; #($exit=$? and handle_exit_gracefully($exit) ) );
        while( my $line = <INPIPE> ){
            chomp $line;
            if ( $line =~ /^>(.*)$/ ){
               $markerlist_hash{$marker}{'fasta_defline'} = $1;
            }
            elsif ( $line !~ /^\s*$/ ){
               $markerlist_hash{$marker}{'sequence'} .= $line;
            }
        }
        close INPIPE;
    }
    print "Finished fetching sequences.\n\n";
}
add_sequences();

sub sort_markers_by_line(\%) {  #make sure that we have the same order as in the markerlist input file
    my $hashref = shift;
    my %sortinghash = %$hashref;
    my @sorted_keys = (sort { "$sortinghash{$a}{'line_number'}" <=> "$sortinghash{$b}{'line_number'}"} keys %sortinghash);
    return @sorted_keys;
}   

sub write_output_files(\@;\@;\@;$;$) {      #could pass hash by reference but let's keeep it simple for now         
    my ($ref1, $ref2, $ref3, $sep, $filename) = @_;
    my @fileheader = @$ref1;
    my @columnames = @$ref2;     # just swap out the old columnames   
    my @markerorder = @$ref3;

    my @first_five=@markerorder[0..5];  
    # First print something to screen to show that all is well
    print "\nBeginning of custom outfile: $filename \n";
    if (@fileheader) { print join("$sep", @fileheader) . "\n" };
    foreach my $marker (@first_five) {
        my @line= ();    #clear preceding line
        foreach my $column (@columnames) {           
            push(@line, $markerlist_hash{$marker}{$column});
        }
        print join("$sep", @line) . "\n";
     }


    unless (open (OUT, ">$filename")){
        die "Can't write to $filename $!";
    }
    (print OUT join("$sep", @fileheader) . "\n") unless (not @fileheader); 
    foreach my $marker (@markerorder) {
        my @line= ();    #clear preceding line
        foreach my $column (@columnames) {           
            push(@line, $markerlist_hash{$marker}{$column});
        }
        print OUT join("$sep", @line) . "\n";
    }
    close OUT;
}

sub write_custom_output_file {   #quick & sloppy but works first some custom format output fil with most of the info gathered
    my @fileheader = ("CHR","marker_id","male_pos","female_pos","contig","pos","fasta_defline","sequence");
    my @columnames = ('CHR','marker_id','male_pos','female_pos','contig','pos','fasta_defline','sequence'); 
    my @markers_by_line = sort_markers_by_line(%markerlist_hash);
    write_output_files(@fileheader, @columnames, @markers_by_line, "\t", "$outpath/flank_seq_$outfilebasename.tsv");
}
write_custom_output_file();

sub write_mapcomp_file {   
    #my @columnames = ('CHR','female_pos','marker_id','sequence');  #not needed
    my @markerorder = sort_markers_by_line(%markerlist_hash);
    
    my $filename ="$outpath/mapcomp_gen_$outfilebasename.csv";
    unless (open (OUT, ">$filename")){
        die "Can't write to $filename $!";
    }
    foreach my $marker (@markerorder) {
        my @line= ();    #clear preceding line
        (my $linkage_group = $markerlist_hash{$marker}{'CHR'}) =~ s/LG// ;  #strip LG from linkage group as in mapcomp files
        push(@line, $specieslabel, $linkage_group);
        push(@line, $markerlist_hash{$marker}{'female_pos'}, "0", $markerlist_hash{$marker}{'marker_id'}, $markerlist_hash{$marker}{'sequence'});
        print OUT join(",", @line) . "\n";
        }
    close OUT;
}
write_mapcomp_file();

sub write_fasta_file {
    my $filename = "$outpath/flank_seq_$outfilebasename.fasta";
    my @markerorder = sort_markers_by_line(%markerlist_hash);
    unless (open (OUT, ">$filename")){
        die "Can't write to $filename $!";
    }
    
    foreach my $marker (@markerorder) {
        my $defline = ">" . $markerlist_hash{$marker}{'line_number'} . ":" ;
        $defline .= "LG=". $markerlist_hash{$marker}{'CHR'} . ":" ;
        $defline .= "female_pos=" . $markerlist_hash{$marker}{'female_pos'} . ":" ;
        $defline .= "marker=" . $markerlist_hash{$marker}{'marker_id'} . ":" ;
        $defline .= "region=" . $markerlist_hash{$marker}{'fasta_defline'} . ":" ;
        $defline .= "marker_pos=" . $markerlist_hash{$marker}{'pos'} . "\n";
        print OUT $defline;
        my $seq = $markerlist_hash{$marker}{'sequence'};
        $seq =~ s/(.{0,80})/$1\n/g;
        chomp($seq);
        print OUT $seq;
    }
}

write_fasta_file();

#write_mapcomp_format_markerlist();
print "\nDone whith this one !\n\n";




