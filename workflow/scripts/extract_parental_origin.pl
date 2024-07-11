#!/usr/bin/env perl 

use warnings; 
use strict;
use Data::Dumper;

# Assume bcftools are available and both vcf input files have been compressed and indexed properly.
### output
# $dnm, GT of DNM (in phase child), $FM_cnt, $MF_cnt, $status, $origin
# chr2_119386326_G_A  0|1     0       3       MF      paternal
@ARGV == 5 or die "$0  <DNM VID list> <Ped> <Child Phase VCF> <Trios Phase VCF>  <OUTFN>\n";

# my $WIN=shift;
my $dnm_vid_fn = shift;
my $ped_fn = shift;
my $child_phase_fn = shift;
my $trios_phase_fn = shift;
my $outfn=shift;


open FOUT, ">$outfn" or die $!;

### get ped sample information ###
open PED, $ped_fn or die $!;
my @sample_ids = (); #child, father, mother
while(<PED>){
    chomp;
    my @e = split /\t/;
    if($e[2] ne "0"){
        @sample_ids = @e[1..3];
        last;
    }

}
close PED;
print join("\t", @sample_ids)."\n";


### Process DNM one by one
open DNM, $dnm_vid_fn or die $!;

my @DNMs = ();
while(<DNM>){
    chomp;
    my $dnm=$_;
    # $dnm =~ s/:/_/g; # should be comment out once id is fixed
    
    print ("Processing DNM $dnm ...\n");

    open my $child_fh, "zcat $child_phase_fn | " or die $!;
    my ($variant_hash, $ps_hash) = &process_vcf($child_fh, @sample_ids);
    close $child_fh;

    open my $trios_fh, " zcat $trios_phase_fn | " or die $!;
    my ($variant_hash2, $ps_hash2) = &process_vcf($trios_fh, @sample_ids);
    close $trios_fh;
    
    # print Dumper $variant_hash2; die;

    my $rv = &lookup_dnm($ps_hash, $dnm);
    my $ps = $rv->{$dnm};
    print "\nProcessing the block $ps with the DNM $dnm ...\n";
    
    my ($FM_cnt, $MF_cnt, $status, $origin)= &identify_dnm_phase($dnm, $ps_hash->{$ps}, $variant_hash, $variant_hash2, $ps_hash2);
    my ($gt2, $status2, $origin2, $num_inform_sites) = &is_FM($dnm, $ps_hash->{$ps},$variant_hash);
    # print "Prediction: ".join("\t", $dnm, $gt, $status, $origin)."\n"; 
    # size of the block (excdue DNM ; ie, -1)
    my $origin_change = sprintf("%s=>%s",$origin2, $origin );
    my $block_size = ($ps eq "")?0:scalar @{$ps_hash->{$ps}} - 1;
    print FOUT join("\t", $dnm, $origin, $variant_hash->{$dnm}->{"child"}->{"GT"}, $block_size , $num_inform_sites, $FM_cnt, $MF_cnt, $status,$origin_change)."\n";

}

close FOUT;
exit;




### print results
# foreach my $k (keys %{$ps_hash}){
#     print join("\t", $k.":", @{$ps_hash->{$k}})."\n";
# }

### for each DNM, we need to know the genotype of DNMs
# and the genotype, GQ of the variatns other than DNMs in the block


########
### get the phase of variants from trios
sub identify_dnm_phase{
    my $dnm = shift;
    my $child_block = shift;
    my $variant_hash = shift;
    my $variant_hash2 = shift;
    my $ps_hash2 = shift;

    # get informative variatns in the block (remove dnm)
    # my @vars=();

    my $FM_cnt = 0;
    my $MF_cnt = 0;

    foreach my $v (@{$child_block}){
        my @gt = map{$variant_hash->{$v}->{$_}->{"GT"} } ("child", "father", "mother");
        next if $gt[1] eq "0/0" && $gt[2] eq "0/0";
        next if $gt[1] eq "1/1" && $gt[2] eq "1/1";
        print "$v...\n";
        
        # get the phase of the variant in trios phased vcf
        if(defined $variant_hash2->{$v} && $variant_hash2->{$v}->{"child"}->{"GT"} =~ /\|/){
            if ($variant_hash2->{$v}->{"child"}->{"GT"} eq $variant_hash->{$v}->{"child"}->{"GT"}){
                $FM_cnt++;
            }else{
                $MF_cnt++;
            }

            print "Trios phased: ".join("\t", $v, $variant_hash2->{$v}->{"child"}->{"GT"}, $variant_hash2->{$v}->{"child"}->{"PS"}, $variant_hash->{$v}->{"child"}->{"GT"})."\n";
        }else{
            # print "no $v\n"; die;
        }
        
        
            
    }


    my $status = "ND";

    if($FM_cnt>0 && $MF_cnt==0){
        $status = "FM";
    }elsif($FM_cnt==0 && $MF_cnt>0){
        $status = "MF";
    }
    
    my $origin = &get_origin($variant_hash->{$dnm}->{"child"}->{"GT"}, $status);
    #print "Prediction: ".join("\t", $dnm, $variant_hash->{$dnm}->{"child"}->{"GT"}, $FM_cnt, $MF_cnt, $status, $origin)."\n";
    return($FM_cnt, $MF_cnt, $status, $origin);
}

### To determine the genotype of child is F|M or not
sub is_FM{
    my $dnm = shift;
    my $variants_in_block = shift;
    my $variant_hash = shift;

    my $FM_incompatible_cnt=0;
    my $MF_incompatible_cnt=0;

    foreach my $v (@{$variants_in_block}){
        next if $v eq $dnm;
        my @gt = map{$variant_hash->{$v}->{$_}->{"GT"} } ("child", "father", "mother");
        next if $gt[1] eq "0/0" && $gt[2] eq "0/0";
        next if $gt[1] eq "1/1" && $gt[2] eq "1/1";

        if($gt[0] eq "0|1" || $gt[0] eq "1|0"){
            my @haps = split /\|/, $gt[0];
            print join("\t", $v, @gt)."\n";

            $FM_incompatible_cnt += 1 if &is_incompatible($haps[0],$gt[1]) > 0 ||  &is_incompatible($haps[1],$gt[2]) >0 ;
            $MF_incompatible_cnt += 1 if &is_incompatible($haps[0],$gt[2]) >0 ||  &is_incompatible($haps[1],$gt[1]) >0 ;
        }
        
    }
    print join("\t", $FM_incompatible_cnt, $MF_incompatible_cnt)."\n";

    my $status = "ND";
    my $origin = "ND";
    my $num_inform_sites=0;

    if($FM_incompatible_cnt ==0 && $MF_incompatible_cnt>0 ){
        $status = "FM";
        $num_inform_sites = $MF_incompatible_cnt;
    }
    if($FM_incompatible_cnt >0 && $MF_incompatible_cnt==0 ){
        $status = "MF";
        $num_inform_sites = $FM_incompatible_cnt;
    }

    my $dnm_gt = $variant_hash->{$dnm}->{"child"}->{"GT"};

    if($status ne "ND"){
        if($dnm_gt eq "0|1"){
            $origin = ($status eq "FM")?"maternal":"paternal";
        }else{
            $origin = ($status eq "FM")?"paternal":"maternal";
        }
    }
    return $dnm_gt, $status, $origin, $num_inform_sites;
}

sub get_origin{
    my $gt = shift;
    my $status = shift;

    my $origin = "ND";

    return $origin if !($gt eq "0|1" || $gt eq "1|0");

    if($status ne "ND"){
        if($gt eq "0|1"){
            $origin = ($status eq "FM")?"maternal":"paternal";
        }else{
            $origin = ($status eq "FM")?"paternal":"maternal";
        }
    }
    return $origin;

}

sub is_incompatible{
    my $child_hap = shift;
    my $parent_gt = shift;
    
    return 0 if($parent_gt eq "./.");
    
    return 0 if $parent_gt =~ /$child_hap/;
    return 1;
}


### Search DNM in the ps blocks

sub lookup_dnm{
    my $ps_hash = shift;
    my @dnms = @_;

    my %rv = map{ $_ => "" } @dnms;

    
    foreach my $k (keys %$ps_hash){
        for my $e (@{$ps_hash->{$k}}){
            if(defined $rv{$e}){
                # $e is dnm
                $rv{$e}=$k;
            }
        }  
    }
    return \%rv;
}

### Process the VCF file and return two hashes
sub process_vcf{
    my $fh = shift;
    my @ids = @_;

    my @hdr=();
    while(<$fh>){
        chomp;
        if(/^#CHROM/){
            @hdr = split /\t/;
            last;
        }
    }
    
    ### Get sample position from the header of vcf
    # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SC109353        SC109336        SC109341
    # my @idxs;
    # foreach my $sample (@ids) {
        
    #     my ($index) = grep { $hdr[$_] eq $sample } (0 .. @hdr-1);
    #     push @idxs, defined $index ? $index : -1;
    # }
    
    my ($child_ind, $father_ind, $mother_ind)  = &get_index(\@hdr, \@ids);
    #print "idx: ".join("\t", $child_ind, $father_ind, $mother_ind)."\n";
    
    my @targeted_tags=("GT", "GQ", "PS");
    my $ps_hash = {};
    my $variant_hash = {};

    while(<$fh>){
        # print $_;

        chomp;
        # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SC109353        SC109336        SC109341
        # chr1    109944636       chr1_109944636_A_G      A       G       32      .       AF=0.166667;AQ=32       GT:DP:AD:GQ:PL:RNC:PS   0/0:102:102,0:50:0,300,2999:..:.        0|1:151:80,71:32:32,0,41:..:109942755   0/0:120:120,0:50:0,300,2999:..:.
        
        my @v = split /\t/;
        my @format_tags = split /:/, $v[8];
        my @tags_index = &get_index(\@format_tags, \@targeted_tags);
        # print "Format: ".join("\t", @tags_index)."\n";

        # save variant information any way
        $variant_hash->{$v[2]} = {};
        
        $variant_hash->{$v[2]}->{child} = &get_format_values($v[8], $v[$child_ind], ("GT", "GQ", "PS"));
        $variant_hash->{$v[2]}->{father} = &get_format_values($v[8], $v[$father_ind], ("GT", "GQ", "PS"));
        $variant_hash->{$v[2]}->{mother} = &get_format_values($v[8], $v[$mother_ind], ("GT", "GQ", "PS"));

        if ($variant_hash->{$v[2]}->{child}->{PS} ne "."){
            if(!defined $ps_hash->{$variant_hash->{$v[2]}->{child}->{PS}}){
                $ps_hash->{$variant_hash->{$v[2]}->{child}->{PS}} = [];
            }
            push @{$ps_hash->{$variant_hash->{$v[2]}->{child}->{PS}}},$v[2];
        } 

    }
    close $fh;
    return $variant_hash, $ps_hash;
}

### Extract format information and return as a hash
# variant -> 
#    child -> GQ/GT

sub get_format_values{
    my $format_str = shift;
    my $value_str = shift;
    my @selected_tags = @_;

    my $rv = {};

    my @format_tags = split /:/, $format_str;
    my @values = split /:/, $value_str;

    my @tags_index = &get_index(\@format_tags, \@selected_tags);

    foreach my $i (0..$#tags_index){
        $rv->{$selected_tags[$i]} = ($tags_index[$i]>=0)?$values[$tags_index[$i]]:".";
    }    
    return $rv;
}

### Given a reference tags array and return the indexes of the values in the reference array
# return -1 if no exits
sub get_index{
    my $tags_ref = shift;
    my $values_ref = shift;

    my @rv = ();
    foreach my $v (@$values_ref) {
        
        my ($index) = grep { $tags_ref->[$_] eq $v } (0 .. @{$tags_ref}-1);
        push @rv, defined $index ? $index : -1;
    }
    return @rv;
}