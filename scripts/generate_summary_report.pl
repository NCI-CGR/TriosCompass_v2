#!/usr/bin/env perl 

### ref: ../IGV_Curation//src/prepare_mie.pl
use strict;
use Excel::Writer::XLSX;

@ARGV==2 or die "$0 <output_xlsx><Pedigress list file>\n";

my $output_xlsx=shift;
my $workbook  = Excel::Writer::XLSX->new($output_xlsx);
my $worksheet = $workbook->add_worksheet();

my $trio_hash={};
while(<>){
    chomp;
    my $fn = $_;
    open FIN, $fn or die $!;
    while(<FIN>){
        # only take the line with full columns
        # t0766c2 SC736756        SC736703        SC730933        2       1
        chomp;
        my @e = split(/\t/);
        #print $e[0]."\n";
        next if $e[3] eq "0";
        my $gender=($e[4]==1)?'M':'F';
        $trio_hash->{$e[0]} = [@e[1..3], $gender];
        
    }
    close FIN;
}


my @hdr=("FamilyID", "ID", "SampleID","FatherSampleID", "MotherSampleID", "Gender", "DNM_Count_DG", "DNM_Count_Strelka", "JIGV_DG", "JIGV_Strelka");

&write_line_excel($worksheet, 0, \@hdr);
### print out hash
my $row_id=1;
for my $k (keys(%{$trio_hash})){
    # We got indel/snp type from the truth file before. Here, we may just have the variant counts instead.
    my $DG_DNM_fn=sprintf("output/GATK_DV/D_and_G.%s.dnm.vcf.gz", $k);
    my $Strelka_DNM_fn=sprintf("output/slivar/strelka_%s.dnm.vcf.gz", $k); 
    my $DG_JIGV_link=sprintf("./output/call_JIGV/D_and_G_%s.JIGV.html", $k);
    my $Strelka_JIGV_link=sprintf("./output/call_JIGV/strelka_%s.JIGV.html", $k);

    ### ref: 38M.process_new_CGR_data.md
    my $DG_DNM_cnt = `zgrep -v -c "^#" $DG_DNM_fn`;
    my $Strelka_DNM_cnt = `zgrep -v -c "^#" $Strelka_DNM_fn`;
    chomp($DG_DNM_cnt);
    chomp($Strelka_DNM_cnt);
    print join("\t", $k, @{$trio_hash->{$k}}, $DG_DNM_cnt, $Strelka_DNM_cnt)."\n";
    my $fam=substr($k, 0,5);
    my @row_data = ($fam, $k, @{$trio_hash->{$k}}, $DG_DNM_cnt, $Strelka_DNM_cnt);
    &write_line_excel($worksheet, $row_id, \@row_data);

    # Add links to igv script
    $worksheet->write_url($row_id, $#row_data + 1, $DG_JIGV_link,  undef, "JIGV HTML" );
    $worksheet->write_url($row_id, $#row_data + 2, $Strelka_JIGV_link,  undef, "JIGV HTML" );
    $row_id++;
}

############################################################################
sub write_line_excel{
    my $ws = shift;
    my $row = shift;
    my $data_ref = shift;

    for my $col (0  .. @{$data_ref} - 1) {
        my $cell_data = $data_ref->[$col];
        $worksheet->write( $row, $col, $cell_data);
    }
}