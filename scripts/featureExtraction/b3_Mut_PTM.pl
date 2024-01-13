use warnings;
use strict;

my $abs_dir="/data1/juns/phase2/";
my $choice="";
my $clust=$ARGV[1];
my $re=$ARGV[0];

my $out_dir=$abs_dir."ML/scaffold/$choice"."_repeat"."$clust/$re/";
my $fea_dir="/data1/juns/PathHost/data_raw/";


#mutation information
my %mu=();
my $dis_avg=0;
my $neu_avg=0;
my $an=0;
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

my $f_dir=$out_dir;
my $o_dir=$f_dir."MutPTM/";
mkdir $o_dir;
open(OUT, ">$o_dir"."CVinput.txt") or die "$!";
open(FA, "$f_dir"."CV.txt") or die "$!";
while (<FA>)
{
    chomp;
    my ($uni,$tag)=split/\s+/,$_;
    my @fea=();
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
    print OUT "$tag\t",join "\t",@fea,"\n";
}
close OUT;
close FA;

