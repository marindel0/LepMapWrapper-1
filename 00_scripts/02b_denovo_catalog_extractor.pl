#! /usr/bin/perl 

use warnings;
use strict;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use Cwd qw(getcwd);


#########################################################################################################
#
# This little script parses the gzipped catalog files from stacks and markerlist output files from
# LepMap3 and appends the DNA sequence of the stacks locus to the markerlist file
# Output format: list as for prepare for chromonomer
# 
# No header
# Species: <string>
# Family: <number>
# Female pos: <float>
# Phenotype: <number?> Not sure about this column - have to look it up in documentation or code.
# Sequence: <string>
# 
# 
# Written by Zophonías O. Jónsson - Oct 18th 2022
#########################################################################################################

#########################################################################################################
# Set parameters - initialize variables
#########################################################################################################

my @arguments = @ARGV; #save arguments as a precaution

### The get options block

my ($help, $markerlist_tsv, $catalog_fasta_gz, $outpath, $discard_below);

#Let's be friendly
sub Usage {
  print "Unknown option: @_\n" if ( @_ );
  print "\nThis is how you use this script: $0

Two input files must be provided with the correct option switches

-m/--markerlist <file>\t:\tThe markerlist file.

-c/--catalog <file>\t:\tThe catalog file (gzipped)\n\n";

print "Optional parameters:

-o/--outpath <path>\t:\tPath to write output files to (defaults to current directory)

-d/--discard <N>\t:\tDiscard stacks with fewer than <N> samples in stacks database (for speedup) defaults to 10

-h/--help\t\t:\tPrint this message.

";
  exit;
}

GetOptions('h|help' => \$help, 'd|discard:i' => \$discard_below, 'm|markerlist=s' => \$markerlist_tsv, 'c|catalog=s' => \$catalog_fasta_gz, 'o|outpath:s' => \$outpath) or die ("Arguments in error!\n");
Usage() if (@arguments == 0 || $help || not ($markerlist_tsv) || not ($catalog_fasta_gz));

$discard_below = 10 if not ($discard_below); # set default to 10 if not provided

# TODO provide option to enter outfilename and select various types of output files

# TODO pass in abbreviation of name for mapcomp file

#(print $markerlist_tsv . " " . $catalog_fasta_gz) or Usage();
### Set global variables
my %markerlist_hash;  
my %catalog_hash;     
my %catalog_cache; 

#my $markerlist_tsv = shift(@arguments) or die "no markerfile $1";
#my $catalog_fasta_gz = shift(@arguments)or die "no catalog file $1"; 

(my $strippedbasename = basename($markerlist_tsv)) =~ s/\.[^.]+$//; 
(my $outfilebasename = $strippedbasename) =~ s/^contig_order_//;
#my $workingdir=getcwd();
my $counter=0;
my @fileheader;   #TODO make local for individual subs 


print "\nMarkerlist: $markerlist_tsv \n";
print "catalog: $catalog_fasta_gz \n";

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




print "Output filename(s): $outfilebasename + extension \n";

print "disregarding stacks with less than $discard_below sample hits in catalog\n";

#   ## begin by figuring out what type of catalog we are dealing with
#   my $firstline =  `zcat $catalog_fasta_gz | head -n 1` or die "can't find the file";
#   print "First line of catalog file: $firstline  => ";
#   
#   my $catalogtype = "denovo" if ( $firstline =~ m/^>\d+\sNS=\d+\s\S+.$/);
#   $catalogtype = "reference_aligned" if ( $firstline =~ m/^>\d+\spos=\S+\sNS=\d+$/);
#   $catalogtype = "undefined" unless defined($catalogtype);
#   #should turn this into a case statement but this works.
#   
#   print "Catalog_type detected as: $catalogtype\n\n";
#   
#   die "This won't work -> play with the source code or talk to the doctor and ask for a new pattern !\n" if ($catalogtype = "undefined");



#########################################################################################################
# Parse the markerlist from LepMap3 and stuff it into a hash
# 
# Passing filname as an argument may seem uneccessarily complicated but I want this to be reusable
#########################################################################################################



sub parse_denovo_markerlist_file($){

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
       
       #Use contig (stack) as nested key (column 2 => 1 in perl) because that's shared with the catalog
       my $marker = $column[4];
       $markerlist_hash{$marker}{'CHR'} = $column[0];
       $markerlist_hash{$marker}{'marker_id'} = $column[1];  
       $markerlist_hash{$marker}{'male_pos'} = $column[2];
       $markerlist_hash{$marker}{'female_pos'} = $column[3];
       $markerlist_hash{$marker}{'contig'} = $column[4];    #redundant but I'll leve it for now
       $markerlist_hash{$marker}{'pos'} = $column[5];
       $markerlist_hash{$marker}{'contig_pos'} = $column[6];  
       $markerlist_hash{$marker}{'line_number'} = ++$linecounter;
       next;
                       
   }                   
   close MARKERLISTFILE ;
}

parse_denovo_markerlist_file($markerlist_tsv);

sub read_fasta_file {

    my $marker_ID='';
    my $stack_counter=0;
    my $discard_counter=0;
    
    
    if ($catalog_fasta_gz =~ /.gz$/) {
        open(INPIPE, "zcat $catalog_fasta_gz |") || die "can’t open pipe to $catalog_fasta_gz"; #or gunzip -c
    }
    else {
        open(INPIPE, $catalog_fasta_gz) || die "can’t open $catalog_fasta_gz";
    }
    
    print "\n* Please be patient. Reading the catalog takes a while *\n\n";
        
    
        my $discard=0;
        NEXTLINE: while( my $line = <INPIPE> ){
            chomp $line;

            if ( $line =~ /^(>.*)$/ ){
                my $defline = $1;
        #print "DEF: $defline\n";
                $stack_counter++;
                $discard=0;
                $defline =~ m/^>(\d+)\sNS=(\d+)/ ;
                my ($stacks_no, $nr_of_samples) = ($1, $2);  #just stuff the rest into $chr_pos for now - should work
        #print "No = $stacks_no, Number of samples = $nr_of_samples";
                if ($nr_of_samples < $discard_below) {
        #die if (! defined($nr_of_samples));
                    $discard_counter++;
                    $discard=1;  #discard until next match
                    next NEXTLINE ;
                }
                $marker_ID = $stacks_no;
                $defline =~ s/^>//;
                $catalog_hash{$marker_ID}{'defline'} = $defline;
                #$catalog_hash{$marker_ID}{'catalog_ID'} = $stacks_no;
                $catalog_hash{$marker_ID}{'nr_of_samples'} = $nr_of_samples;  #number of samples hitting catalog (useful for calculating missingless)
                #$catalog_hash{$marker_ID}{'contiginfo'} = $contiginfo;
            }
            elsif ( $line !~ /^\s*$/ and not $discard ){
            $catalog_hash{$marker_ID}{'sequence'} .= $line;
            }
#die if ($stack_counter > 5);
        }      
    
    
    close(INPIPE);
    print "Done reading fasta file \n";
    print "Discarded $discard_counter and kept " . ($stack_counter - $discard_counter) . " out of $stack_counter stacks in catalog file.\n";
}
read_fasta_file();  #No need to be fancy and pass reference to catalog_hash (maybe later if this script grows)

my ($matches, $misses) = (0,0);


sub sequences_2_markerlisthash {  #this is the brute force approach - optimize later by breaking up into LG's

    print "Start cross referencing catalogs and markerlist\n";
    
    foreach my $stack_key ( keys %catalog_hash ) {
       # my $stack_contig = $catalog_hash{$stack_key}{'contig'};
       # my $stack_pos = $catalog_hash{$stack_key}{'pos'};

        foreach my $markerlist_key (keys %markerlist_hash) {     #move comparison here - save some time
            #my $markerlist_contig = $markerlist_hash{$markerlist_key}{'contig'};            
            #my $marker_hits=0;

            if ($markerlist_key eq $stack_key) { 
                $matches++;
                $markerlist_hash{$markerlist_key}{'catalog_defline'} = $catalog_hash{$stack_key}{'defline'};
                $markerlist_hash{$markerlist_key}{'nr_of_samples'} = $catalog_hash{$stack_key}{'nr_of_samples'};
                $markerlist_hash{$markerlist_key}{'sequence'} = $catalog_hash{$stack_key}{'sequence'} ;
                #$markerlist_hash{$markerlist_key}{'catalog_ID'} = $catalog_hash{$stack_key}{'catalog_ID'};
                #$markerlist_hash{$markerlist_key}{'strandpos'} = $catalog_hash{$stack_key}{'strandpos'};
            }
        }
    }
    
    print "\nPadding gaps in catalog (in case cutoff was to aggressife) \n";
    
    foreach my $markerlist_key (keys %markerlist_hash) {      #iterate through markerlist hash 

        if (not defined ($markerlist_hash{$markerlist_key}{'sequence'}) ) {                                          # if no entry in catalog cache then there is no info
            $misses++;
            $markerlist_hash{$markerlist_key}{'catalog_defline'} = "NA";
            $markerlist_hash{$markerlist_key}{'nr_of_samples'} = "NA";
            $markerlist_hash{$markerlist_key}{'sequence'} = "NA";
            #$markerlist_hash{$markerlist_key}{'catalog_ID'} = "NA";
            #$markerlist_hash{$markerlist_key}{'strandpos'} = "NA";
            #$markerlist_hash{$markerlist_key}{'distance'} = "NA";
            #$markerlist_hash{$markerlist_key}{'alternate_stacks'} = "NA";
        }
    }  
}

sequences_2_markerlisthash();

print "\n";
print "Found: $matches matching catalog stacks. $misses markers were missing from the catalog (truncated below $discard_below samples).\n";

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
    print "\nBeginning of outfile\n";
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
    my @fileheader = ("CHR","marker_id","male_pos","female_pos","contig","pos","catalog_defline","nr_of_samples","sequence");
    my @columnames = ('CHR','marker_id','male_pos','female_pos','contig','pos','catalog_defline','nr_of_samples','sequence'); 
    my @markers_by_line = sort_markers_by_line(%markerlist_hash);
    write_output_files(@fileheader, @columnames, @markers_by_line, "\t", "$outpath/union_denovo_$outfilebasename.tsv");
}
write_custom_output_file();

sub write_mapcomp_file {   
    #my @columnames = ('CHR','female_pos','marker_id','sequence');  #not needed
    my @markerorder = sort_markers_by_line(%markerlist_hash);
    my $filename ="$outpath/mapcomp_cat_$outfilebasename.csv";
    unless (open (OUT, ">$filename")){
        die "Can't write to $filename $!";
    }
    foreach my $marker (@markerorder) {
        my @line= ();    #clear preceding line
        (my $linkage_group = $markerlist_hash{$marker}{'CHR'}) =~ s/LG// ;  #strip LG from linkage group as in mapcomp files
        push(@line, basename($markerlist_tsv), $linkage_group);   #TODO fix this like in the genome extractor.
        push(@line, $markerlist_hash{$marker}{'female_pos'}, "0", $markerlist_hash{$marker}{'marker_id'}, $markerlist_hash{$marker}{'sequence'});
        print OUT join(",", @line) . "\n";
        }
    close OUT;
}
write_mapcomp_file();

sub write_fasta_file {
    my $filename = "$outpath/stacks_from_$outfilebasename.fasta";
    my @markerorder = sort_markers_by_line(%markerlist_hash);
    unless (open (OUT, ">$filename")){
        die "Can't write to $filename $!";
    }
    foreach my $marker (@markerorder) {
        my $defline = ">". $markerlist_hash{$marker}{'CHR'} . ":" ;
        $defline .= "female_pos=" . $markerlist_hash{$marker}{'female_pos'} . ":" ;
        $defline .= "marker=" . $markerlist_hash{$marker}{'marker_id'} . ":" ;
        $defline .= "contig=" . $markerlist_hash{$marker}{'contig'} . ":" ;
        $defline .= $markerlist_hash{$marker}{'pos'} . "\n";
        print OUT $defline;
        my $seq = $markerlist_hash{$marker}{'sequence'};
        $seq =~ s/(.{0,80})/$1\n/g;
        chomp($seq);
        print OUT $seq;
    }
}

write_fasta_file();

#write_mapcomp_format_markerlist();
print "\nDone!\n\n";
exit 1;


