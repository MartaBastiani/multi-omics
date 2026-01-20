#!/home/marta/miniconda3/bin/perl

use strict;
use warnings;
use Data::Dumper;
use FindBin;
use Getopt::Std;

use lib "$FindBin::Bin/./";
use DbUtils;

getopts("f:d:h", \my %opts);

my $file = $opts{f};
my $db = $opts{d};

if ($opts{h}) {
    print_usage();
    exit;
}

unless($opts{f} && $opts{d}) {
    print_usage();
    exit;
}

my %expData = (); 

open(my $fhf, "<", $file)
    or die "Could not open file '$file': $!";

while(my $gene = <$fhf>) {
    chomp $gene;
    my @gene = split(/\t/, $gene);
    $gene[0] =~ s/\.\d+//;

    if (!$expData{$gene[0]}) {
        %{$expData{$gene[0]}} = (id=>$gene[0], dev=>$gene[1], pval=>$gene[2]);
    }
}
close $fhf;

open(my $fhdb, "<", $db)
    or die "Could not open file '$db': $!";

# locate geneID basing on the db and store it in a trueID key
if ($db =~ /hsa.gene_info/) {
    %expData = checkGeneInfo(
        fileHandle => $fhdb,
        expData => \%expData
    );

} elsif ($db =~ /hsa_mirtarbase.tsv/) {
    %expData = checkmiRNInfo(
        fileHandle => $fhdb,
        expData => \%expData
    );

} elsif ($db =~ /hsa_uniprot.tsv/) {
    %expData = checkProtInfo(
        fileHandle => $fhdb,
        expData => \%expData
    );

} else {
    die "Please choose a proper database file\n";
}
close $fhdb;


my $file_new = $file . "_newId";

open(my $fhfw, ">", $file_new)
    or die "Could not open output file '$file_new': $!";

foreach my $key (sort keys %expData) {
    # if geneId exists print it
    if (${$expData{$key}}{trueId}) {
        print $fhfw "${$expData{$key}}{trueId}\t${$expData{$key}}{dev}\t${$expData{$key}}{pval}\n";
    } else {
    # if it doesn't, use the gene symbol
        print $fhfw "${$expData{$key}}{id}\t${$expData{$key}}{dev}\t${$expData{$key}}{pval}\n";
    }
}

close $fhfw;




sub print_usage{
    print << "HOW TO";
    
      USAGE: ./main.pl -f <file_path> -d <database_path>



<file_path>         The path to the tsv file containing gene 
                    identifications  (gene symbols or Ensembl gene 
                    IDs), log2 fold changes and adjusted p-values.
            
<database_path>     The path to a database (either hsa.gene_info,
                    hsa_mirtarbase.tsv or hsa_uniprot.tsv).   

HOW TO
}

exit;