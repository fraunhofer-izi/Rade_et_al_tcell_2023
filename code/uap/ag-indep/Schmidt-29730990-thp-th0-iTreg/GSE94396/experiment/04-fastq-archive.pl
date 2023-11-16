#/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use feature ':5.18';
use File::Basename;
use File::Find;

my $queryDIR = "/data/RNA-Seq/2018-MAVO-MCF/Schmidt-29730990-thp-th0-iTreg/GSE94396/Unaligned/Project_Schmidt-29730990-GSE94396-thp-th0-iTreg/";
my $metaDataFILE = "SraRunTable";

my @raw_files;
my @metaData;

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Push fastq file paths to array
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
find sub {
    return unless -f;
    return unless /\.*fastq.gz$/;
    # print "$File::Find::name\n";
    push @raw_files, &trim($File::Find::name);
}, $queryDIR;
# print Dumper(\@raw_files);

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## make GSM.* folder and move SRR.* files match with GSM.* from metadata there.
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
open(METAFILE, "$metaDataFILE") or die "FEHLER: $!";
while(<METAFILE>){
    chomp;
    next if $. < 2; # Skip first line
    push @metaData, &trim($_);
}
close METAFILE;
# print Dumper(\@metaData);

foreach my $x (@metaData) {
    chomp $x;

    my @col = split(/\t/, $x);
    my $srrIDMeta = &trim($col[5]);
    $srrIDMeta =~ s/"//g;
    my $gsmIDMeta = &trim($col[7]);
    $gsmIDMeta =~ s/"//g;

    # say $srrIDMeta;
    # say $gsmIDMeta;

    my $newFolderName = "Sample_".$gsmIDMeta;

    foreach my $fileNameToParse (@raw_files) {

        $fileNameToParse = &trim(basename($fileNameToParse));
        $fileNameToParse =~ /(.*)(_[12]).fastq.gz/;
        # say $fileNameToParse;
        my $sampleID = &trim($1);
        my $read = &trim($2);

        if ($sampleID eq $srrIDMeta) {
            if (!-d "$queryDIR/$newFolderName") {
            # say "mkdir $queryDIR/$newFolderName";
            # say "touch $queryDIR/$newFolderName/SampleSheet.csv";
            `mkdir $queryDIR/$newFolderName`;
            `touch $queryDIR/$newFolderName/SampleSheet.csv`;
            }
        if ($sampleID eq $srrIDMeta) {
            # say "mv $queryDIR/$fileNameToParse $queryDIR/$newFolderName";
            `mv $queryDIR/$fileNameToParse $queryDIR/$newFolderName`;
        }
        }
    }
}

sub trim(){
  my $string=shift;
  $string=~s/^\s+//;
  $string=~s/\s+$//;
  return $string;
}
