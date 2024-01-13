use warnings;
use strict;

my $abs_dir="/data1/juns/phase2/";
my $choice="";
my $clust=$ARGV[1];
my $re=$ARGV[0];

my $out_dir=$abs_dir."ML/scaffold/$choice"."_repeat"."$clust/$re/";
my $fea_dir="./data_raw/";
#age
my %age=();
my $sum=0;
my $age_avg=0;
my $an=0;
open(FB, "$fea_dir"."HUMAN_PPODv4_Jaccard_wagner1.0_age-time.protein_list") or die "$!";
while (<FB>)
{
    chomp;
    next if /^#/;
    my @item=split/\t/;
    $age{$item[1]}=$item[-1];
    $sum+=$item[-1];
    $an++;
}
close FB;
$age_avg=$sum/$an;

#dnds
my %dnds=();
$sum=0;
my $dnds_avg=0;
$an=0;
open(FB, "$fea_dir"."uniprot_dN_dS_ratio.txt") or die "$!";
while (<FB>)
{
    chomp;
    my @item=split/\t/;
    $dnds{$item[0]}=$item[5];
    $sum+=$item[5];
    $an++;
}
close FB;
$dnds_avg=$sum/$an;

my %e=();
open(IN, "$fea_dir"."Essential.txt") or die "$!";
while (<IN>)
{
    chomp;
    $e{$_}++
}
close IN;

my %h=();
open(IN, "$fea_dir"."housekeeping.txt") or die "$!";
while (<IN>)
{
    chomp;
    $h{$_}++
}
close IN;

#expression
my %pax=();
$sum=0;
my $pax_avg=0;
$an=0;
open(IN,"/data1/juns/phase/PAX/pax.txt") or die "$!";
while (<IN>)
{
    chomp;
    next if /^uniprot/;
    my @item=split/\t/;
    $pax{$item[0]}=$item[1];
    $an++;
    $sum+=$item[1];
}
$pax_avg=$sum/$an;
close IN;


my $f_dir=$out_dir;
my $o_dir=$f_dir."ExpEnsProdnds/";
mkdir $o_dir;
open(OUT, ">$o_dir"."CVinput.txt") or die "$!";
open(FA, "$f_dir"."CV.txt") or die "$!";
while (<FA>)
{
    chomp;
    my ($uni,$tag)=split/\s+/,$_;
    my @fea=();
    if (exists $e{$uni})
    {
        push(@fea,1);
    }else{
         push(@fea,0);
    }
    
     if (exists $h{$uni})
    {
        push(@fea,1);
    }else{
         push(@fea,0);
    }

    if (exists $age{$uni})
    {
        push(@fea,$age{$uni})
    }else{
        push(@fea,$age_avg)
    }
    
    if (exists $dnds{$uni})
    {
        push(@fea,$dnds{$uni})
    }else{
        push(@fea,$dnds_avg)
    }
    if (exists $pax{$uni})
    {
        push(@fea,$pax{$uni});
    }else{
        push(@fea,$pax_avg)
    }
    
    print OUT "$tag\t",join "\t",@fea,"\n";
}
close OUT;
close FA;

open(OUT, ">$o_dir"."Intestinput.txt") or die "$!";
open(FA, "$f_dir"."Intest.txt") or die "$!";
while (<FA>)
{
   chomp;
    my ($uni,$tag)=split/\s+/,$_;
    my @fea=();
    
    if (exists $e{$uni})
    {
        push(@fea,1);
    }else{
         push(@fea,0);
    }
    
     if (exists $h{$uni})
    {
        push(@fea,1);
    }else{
         push(@fea,0);
    }

    if (exists $age{$uni})
    {
        push(@fea,$age{$uni})
    }else{
        push(@fea,$age_avg)
    }
    
    if (exists $dnds{$uni})
    {
        push(@fea,$dnds{$uni})
    }else{
        push(@fea,$dnds_avg)
    }
    if (exists $pax{$uni})
    {
        push(@fea,$pax{$uni});
    }else{
        push(@fea,$pax_avg)
    }
    print OUT "$tag\t",join "\t",@fea,"\n";
}
close OUT;
close FA;

