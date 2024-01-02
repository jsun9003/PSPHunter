use warnings;
use strict;

my $abs_dir="/data1/juns/phase2/";
my $data_dir="/data1/juns/PathHost/data_raw/";
my $choice="$ARGV[0]";#"scaffold"
#age
my %age=();
my $sum=0;
my $age_avg=0;
my $an=0;
open(FB, "$data_dir"."HUMAN_PPODv4_Jaccard_wagner1.0_age-time.protein_list") or die "$!";
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
open(FB, "$data_dir"."uniprot_dN_dS_ratio.txt") or die "$!";
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
open(IN, "$data_dir"."Essential.txt") or die "$!";
while (<IN>)
{
    chomp;
    $e{$_}++
}
close IN;

my %h=();
open(IN, "$data_dir"."housekeeping.txt") or die "$!";
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

#mutation information
my %mu=();
my $dis_avg=0;
my $neu_avg=0;
$an=0;
my $sum11=0;
my $sum22=0;
open(IN,"/data1/juns/phase/Disease/HumanVarBase_phase_length.txt") or die "$!";
while (<IN>)
{
    chomp;
    next if /^uniprot/;
    my @item=split/\t/;
    my @ref=($item[1],$item[2]);
    $mu{$item[0]}=$item[1];
    $sum11+=$item[1];
    $sum22+=$item[2];
    $an++;
}
close IN;
$dis_avg=$sum11/$an;
$neu_avg=$sum22/$an;
#PTM information
my %ptm=();
$an=0;
my ($avg11,$avg22,$avg33,$avg44,$avg55,$avg66)=(0,0,0,0.0.0);
my ($sum111,$sum222,$sum33,$sum44,$sum55,$sum66)=(0,0,0,0,0,0);
open(IN,"/data1/juns/phase/PTM/PTM_length.txt") or die "$!";
while (<IN>)
{
    chomp;
    next if /^Uniprot/;
    my @item=split/\t/;
    my @ref=($item[2],$item[3],$item[4],$item[5],$item[6]);
    $ptm{$item[0]}=\@ref;
    $sum111+=$item[2];
    $sum222+=$item[3];
    $sum33+=$item[4];
    $sum44+=$item[5];
    $sum55+=$item[6];
    $sum66+=$item[7];
    $an++;
}
close IN;
($avg11,$avg22,$avg33,$avg44,$avg55,$avg66)=($sum111/$an,$sum222/$an,$sum33/$an,$sum44/$an,$sum55/$an,$sum66/$an);

#topo
my %topo=();
open(IN, "/data1/juns/phase2/data_raw/ppi/feature_net.txt") or die "$!";
my ($avg1,$avg2,$avg3,$avg4)=(0,0,0,0);
my ($sum1,$sum2,$sum3,$sum4)=(0,0,0,0);
my $c=0;
while (<IN>)
{
    chomp;
    next if /^uniprot/;
    my @item=split/\t/;
    my @ref=@item[1..$#item];
    $topo{$item[0]}=\@ref;
    $sum1+=$item[1];
    $sum2+=$item[2];
    $sum3+=$item[3];
    $sum4+=$item[4];
    $c++;
}
close IN;
($avg1,$avg2,$avg3,$avg4)=($sum1/$c,$sum2/$c,$sum3/$c,$sum4/$c);

open(IN, "/data1/juns/phase2/data_raw/ppi_string/ppi_uniprot/Mixmerge.txt") or die "$!";
while (<IN>)
{
    chomp;
    next if /^uniprot/;
    my @item=split/\t/;
    my @ref=@item[1..$#item];
    $topo{$item[0]}=\@ref;
}
close IN;


my $f_dir=$abs_dir."ML/PHApply/scaffold/$choice/test/";
my $o_dir=$f_dir."mergeFun/";
mkdir $o_dir;

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
    
    if (exists $mu{$uni})
    {
        my $vec=$mu{$uni};
        push(@fea,$vec);
    }else{
        my $vec=$dis_avg;
        push(@fea,$vec);
    }

    if (exists $ptm{$uni})
    {
        my @vec=@{$ptm{$uni}};
        push(@fea,@vec);
    }else{
        my @vec=($avg11,$avg22,$avg33,$avg44,$avg55);
        push(@fea,@vec);
    }
    
    if (exists $topo{$uni})
    {
        my @vec=@{$topo{$uni}};
        push(@fea,@vec);
    }else{
        my @vec=($avg1,$avg2,$avg3,$avg4);
        push(@fea,@vec);
    }
    print OUT "$tag\t",join "\t",@fea,"\n";
}
close OUT;
close FA;

