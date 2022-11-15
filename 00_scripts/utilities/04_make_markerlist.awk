#!/usr/bin/awk -f

BEGIN {
    POSFILE=ARGV[2];
    line=0;
    while (getline < POSFILE) {
       ++line;
       pos[line]=$0;   
    }
    close(POSFILE);
    line=0
    print "marker_id contig pos contig_pos"
}
(FNR==NR)&&(!/#java|CHR/) {
       line++;
       a = line +7;
       print line, $0, pos[a], $0"_"pos[a];
}
        
