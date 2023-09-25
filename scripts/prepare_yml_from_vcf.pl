#!/usr/bin/env perl 

use warnings; 
use strict;

use YAML::Tiny;
use Excel::Writer::XLSX;



# Parse the input file and write YAML output 
# MIE: Mendelian Inheritance Error in vcf format
# Also save the result into Excel format with http links to local files.
# use igv_case_bed_2_yaml.pl as a reference

# perl src/prepare_mie.pl Files/Meredith_Case1/t0008c1.mie.pt.filt.mdnm.relaxed.txt Files/Meredith_Case1/trios_final_IDs.txt /DCEG/Scimentis/DNM/data/BATCH2_b38 my.yaml  test.xlsx 

@ARGV == 6 or die "$0 <VCF input> <Ped file> <BAM_DIR> <YAML output file> <Excel Table> <snapshot output folder>\n"; 

my $vcffn = shift;
my $pedfn = shift;
my $bam_dir = shift; # SC742321.dedup.bam
my $yaml_outfn = shift;
my $excel_outfn = shift;
my $snapshot_dir = shift;

### Extract kid, dad, mom ids 
open PED, "<", $pedfn or die "$pedfn: $!";
my @ids=();
my $fam;

while(<PED>){
# t0007c1 SC499427        SC499423        SC499428        1       1
    my @e = split /\t/;
    if($e[2] ne "0"){
        @ids=@e[1..3];
        $fam = $e[0];
        print join("\t", @ids)."\n";
        last;
    }

}
close PED;

### Extract the DNMs from the VCF file
my $cat_cmd = $vcffn =~ /gz$/ ?'zcat':'cat';
open (VCF, "$cat_cmd $vcffn |") or die "$vcffn: $!";

my @hdr=();
while(<VCF>){
    chomp;
    if(/^#CHROM/){
        @hdr = split /\t/;
        last;
    }
}



### Get sample position from the header of vcf
# #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SC109353        SC109336        SC109341
my @idxs;
foreach my $sample (@ids) {
    
    my ($index) = grep { $hdr[$_] eq $sample } (0 .. @hdr-1);
    push @idxs, defined $index ? $index : -1;
}
# print "idx: ".join("\t", @idxs)."\n";

my @new_hdr= (@hdr[0..8,@idxs], "IGV_snapshot", "IGV_script");

# print "header: ".join("\t", @new_hdr)."\n";

### SC742321.dedup.bam 
my @bam_files = map{$bam_dir.'/'.$_.".dedup.bam"} @ids; #needs to be refined in the future



### Create an Excel file
my $workbook  = Excel::Writer::XLSX->new($excel_outfn);
my $worksheet = $workbook->add_worksheet();

### Create a new object with a single hashref document
my $yaml = YAML::Tiny->new();

# The whole YAML output is an array with only one element. 

&write_line_excel($worksheet, 0, \@new_hdr);
my $row_id=1;

my $entry = {};
my $snapshots = [];

### add bam files and names to yaml
$entry->{name}=$fam;
$entry->{bam_files}= \@bam_files;


### continue to process VCF files 
while(<VCF>){
    # print $_;

    chomp;
    # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SC109353        SC109336        SC109341
    
    my @v = split /\t/;

    my $snapshot_name = sprintf("%s_MIE_%05d", $fam, $row_id);
    my $sp_hash = {};
    $sp_hash->{name} = $snapshot_name;
    $sp_hash->{chr} = $v[0];
    $sp_hash->{start} = int($v[1]); # 1-based cooordinates
    $sp_hash->{stop} = int($v[1])+length($v[3])-1;
    push @$snapshots, $sp_hash;

    my @row_data = @v[0..8,@idxs];
    # print join("\t", $cid, $mani_hash{$cid}, $mani_hash{$fid}, $mani_hash{$mid})."\n";

    # Besides, each row together with the UID (and the additional http links to the local scriptsa and pngs ) will also be saved as Excel file,
    &write_line_excel($worksheet, $row_id, \@row_data);

    # Add links to igv script
    $worksheet->write_url($row_id, $#row_data + 1, $snapshot_dir."/".$fam."/".$sp_hash->{name}.".png",  undef, $sp_hash->{name}.".png" );
    $worksheet->write_url($row_id, $#row_data + 2, $snapshot_dir."/".$fam."/".$sp_hash->{name}.".bat",  undef, $sp_hash->{name}.".bat" );
    $row_id++;
}

close VCF;
$workbook->close();

# Create the only entry for the whole yaml file
$entry->{snapshots} = $snapshots;

my $entries = [];
push @$entries, $entry;
push @$yaml, $entries; 
$yaml->write($yaml_outfn);

exit 0;

###########
sub write_line_excel{
    my $ws = shift;
    my $row = shift;
    my $data_ref = shift;

    for my $col (0  .. @{$data_ref} - 1) {
        my $cell_data = $data_ref->[$col];
        $ws->write( $row, $col, $cell_data);
    }
}