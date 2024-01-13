use warnings;
use strict;

my $abs_dir="/data1/juns/phase2/";
my $choice="";
my $clust=$ARGV[1];
my $re=$ARGV[0];

my $out_dir=$abs_dir."ML/scaffold/$choice"."_repeat"."$clust/$re/";
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


my $f_dir=$out_dir;
my $o_dir=$f_dir."Topo/";
mkdir $o_dir;
open(OUT, ">$o_dir"."CVinput.txt") or die "$!";
open(FA, "$f_dir"."CV.txt") or die "$!";
while (<FA>)
{
    chomp;
    my ($uni,$tag)=split/\s+/,$_;
    my @fea=();
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

open(OUT, ">$o_dir"."Intestinput.txt") or die "$!";
open(FA, "$f_dir"."Intest.txt") or die "$!";
while (<FA>)
{
   chomp;
    my ($uni,$tag)=split/\s+/,$_;
    my @fea=();
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

