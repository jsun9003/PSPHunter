use warnings;
use strict;

my $abs_dir="/data3/juns/phase2/";
my $fea_dir="/data3/juns/phase/";
my $choice="scaffold";
my $clust="30";
my $re=$ARGV[0];

my $out_dir=$abs_dir."ML/$choice/repeat"."$clust/$re/";

my $data_dir="/data3/juns/PathHost/data_raw/";


my $f_dir=$out_dir;
my $o_dir=$f_dir."FunSite/";
mkdir $o_dir;

my %hash=();
open(FA, "$data_dir"."Pathhost_ASA_tag.txt") or die "$!";
while (<FA>)
{
    chomp;
    my ($uni,$pert1,$pert2,$pert3,$pert4)=split/\s+/;
    my @ref=($pert1,$pert2,$pert3,$pert4);
    $hash{$uni}=\@ref;
}

my %c=();
my $disorder_avg=0;
my $an=0;
my $sum=0;
open(FA, "$fea_dir"."IDR/SPINE-D-local.txt") or die "$!";
while (<FA>)
{
    chomp;
    next if /^uniprot/;
    my ($uni,$dis,$l,$count,$pert)=split/\s+/;
    next if exists $c{$uni};
    $c{$uni}=$pert;
    $sum+=$pert;
    $an++;
}
$disorder_avg=$sum/$an;
my %rna=();
my $rna_avg=0;
$an=0;
$sum=0;
open(FA, "$fea_dir"."NAR_binding/RNA_binding.txt") or die "$!";
while (<FA>)
{
    chomp;
    next if /^uniprot/;
    my ($uni,$dis,$l,$count,$pert)=split/\s+/;
    next if exists $rna{$uni};
    $rna{$uni}=$pert;
    $sum+=$pert;
    $an++;
}
$rna_avg=$sum/$an;

my %dna=();
my $dna_avg=0;
$an=0;
$sum=0;
open(FA, "$fea_dir"."NAR_binding/DNA_binding.txt") or die "$!";
while (<FA>)
{
    chomp;
    next if /^uniprot/;
    my ($uni,$dis,$l,$count,$pert)=split/\s+/;
    next if exists $dna{$uni};
    $dna{$uni}=$pert;
    $sum+=$pert;
    $an++;
}
$dna_avg=$sum/$an;

#Spscore
my %ps=();
my $ps_avg=0;
$an=0;
$sum=0;
open(IN,"$fea_dir"."PScore/pscore.txt") or die "$!";
while (<IN>)
{
    chomp;
    next if /^uniprot/;
    my @item=split/\t/;
    $ps{$item[0]}=$item[1];
    $sum+=$item[1];
    $an++;
}
close IN;
$ps_avg=$sum/$an;

#predicted mutation
my %mu=();
my $dis_avg=0;
my $neu_avg=0;
$an=0;
my $sum1=0;
my $sum2=0;
open(IN,"$fea_dir"."Disease/rhapsody_length.txt") or die "$!";
while (<IN>)
{
    chomp;
    next if /^uniprot/;
    my @item=split/\t/;
    $mu{$item[0]}=$item[1];
    $sum1+=$item[1];
    $an++;
}
close IN;
$dis_avg=$sum1/$an;

my %ptm=();
my $phos_avg=0;
my $meth_avg=0;
my $nitro_avg=0;
my $palm_avg=0;
$an=0;
$sum1=0;
$sum2=0;
my $sum3=0;
my $sum4=0;
open(FA, "$fea_dir"."PTM/predict/PTM_length.txt") or die "$!";
while (<FA>)
{
    chomp;
    next if /^Uniprot/;
    my ($uni,$type,$ace,$meth,$phos,$sumo,$nitro,$palm)=split/\s+/;
    next if exists $ptm{$uni};
    my @vec=($phos,$meth,$nitro,$palm);
    $ptm{$uni}=\@vec;
    $sum1+=$phos;
    $sum2+=$meth;
    $sum3+=$nitro;
    $sum4+=$palm;
    $an++;
}
$phos_avg=$sum1/$an;
$meth_avg=$sum2/$an;
$nitro_avg=$sum3/$an;
$palm_avg=$sum4/$an;

my %mem=();
open(IN, "$data_dir"."mem.txt") or die "$!";
while (<IN>)
{
    chomp;
    my @item=split/\t/;
    my ($name,$mnum)=split/=/,$item[4];
    $mem{$item[0]}=$mnum;
}
close IN;
$an=0;
$sum=0;
my $le_avg=0;
my %le=();
open(IN,"$fea_dir"."proteinLength/proteinLength.txt") or die "$!";
while (<IN>)
{
    chomp;
    next if /^uniprot/;
    my @item=split/\t/;
    $le{$item[0]}=$item[1];
    $sum+=$item[1];
    $an++;
}
close IN;
$le_avg=$sum/$an;

open(OUT, ">$o_dir"."CVinput.txt") or die "$!";
open(FA, "$f_dir"."CV.txt") or die "$!";
while (<FA>)
{
    chomp;
    my ($uni,$tag)=split/\s+/,$_;
    my @fea=();
    if (exists $hash{$uni})
    {
       push (@fea,@{$hash{$uni}})
    }else{
       push (@fea,qw(0 0 0 0))
    }
    
    if (exists $c{$uni})
    {
        push(@fea,$c{$uni});
    }else{
        push(@fea,$disorder_avg);
    }
    
     if (exists $rna{$uni})
    {
        push(@fea,$rna{$uni});
    }else{
        push(@fea,$rna_avg);
    }
    
    if (exists $dna{$uni})
    {
        push(@fea,$dna{$uni});
    }else{
        push(@fea,$dna_avg);
    }
    if (exists $ps{$uni})
    {
        push(@fea,$ps{$uni});
    }else{
        push(@fea,$ps_avg);
    }
    
    my $vec1="";
    if (exists $mu{$uni})
    {
        $vec1=$mu{$uni};
    }else{
        $vec1=$dis_avg;
    }
    my @vec1=();
    if (exists $ptm{$uni})
    {
        @vec1=@{$ptm{$uni}};
    }else{
        @vec1=($phos_avg,$meth_avg,$nitro_avg,$palm_avg);
    }
    push(@fea,$vec1,@vec1);
    
    if (exists $mem{$uni})
    {
        push(@fea,$mem{$uni});
    }else{
        push(@fea,0);
    }
    if (exists $le{$uni})
    {
        push(@fea,$le{$uni});
    }else{
        push(@fea,$le_avg);
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
    if (exists $hash{$uni})
    {
       push (@fea,@{$hash{$uni}})
    }else{
       push (@fea,qw(0 0 0 0))
    }
    
    if (exists $c{$uni})
    {
        push(@fea,$c{$uni});
    }else{
        push(@fea,$disorder_avg);
    }
    
     if (exists $rna{$uni})
    {
        push(@fea,$rna{$uni});
    }else{
        push(@fea,$rna_avg);
    }
    
    if (exists $dna{$uni})
    {
        push(@fea,$dna{$uni});
    }else{
        push(@fea,$dna_avg);
    }
    if (exists $ps{$uni})
    {
        push(@fea,$ps{$uni});
    }else{
        push(@fea,$ps_avg);
    }
    
    my $vec1="";
    if (exists $mu{$uni})
    {
        $vec1=$mu{$uni};
    }else{
        $vec1=$dis_avg;
    }
    my @vec1=();
    if (exists $ptm{$uni})
    {
        @vec1=@{$ptm{$uni}};
    }else{
        @vec1=($phos_avg,$meth_avg,$nitro_avg,$palm_avg);
    }
    push(@fea,$vec1,@vec1);
    
    if (exists $mem{$uni})
    {
        push(@fea,$mem{$uni});
    }else{
        push(@fea,0);
    }
    if (exists $le{$uni})
    {
        push(@fea,$le{$uni});
    }else{
        push(@fea,$le_avg);
    }
    print OUT "$tag\t",join "\t",@fea,"\n";
}
close OUT;
close FA;

