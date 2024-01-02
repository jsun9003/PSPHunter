use warnings;
use strict;

my $abs_dir="/data1/juns/phase2/";
my $choice="$ARGV[0]";
my $o_dir="/data1/juns/phase2/ML/PHApply/scaffold/$choice/mergeFeature/";
mkdir $o_dir;

open(OUT,">$o_dir"."InProbAVG.txt") or die "$!";

my $out_dir=$abs_dir."ML/PHApply/scaffold/$choice/test/";
open(IN,"$out_dir"."Intest.txt") or die "$!";
my $n=0;
my %tag=();
while (<IN>)
{
    chomp;
    $n++;
    my ($uni,$tag)=split;
    $tag{$n}=$uni;
}

close IN;


my %scaffold=();
for(my $i=1;$i<=100;$i++)
{
    open(IN,"$out_dir"."mergeFeature/InProbphase$i.txt") or die "$!";
    my $k=0;
    while (<IN>)
    {
        chomp;
        $k++;
        if (exists $scaffold{$k})
        {
            $scaffold{$k}+=$_;
        }else{
            $scaffold{$k}=$_;
        }
    }
    close IN;
}

my %map=();
open(IN, "/data1/juns/data/uniprot/uniprot_mod.txt") or die "$!";
while (<IN>)
{
    chomp;
    next if /uniprot/;
    my @item=split;
    $map{$item[0]}=$item[1];
}
close IN;

#known
my %known=();
open(IN, "/data1/juns/phase2/data_raw/integrate/anno_list.txt") or die "$!";
print OUT "uniprot\ttype\tgene\tscaffold\n";
$n=0;
while (<IN>)
{
    chomp;
    my ($uni,$tag)=split;
    print OUT $_,"\t";
    $n++;
    
    my $index=$tag{$n};
    if ($uni eq $index)
    {
        
    }else{
        print $_,"\n";
        die;
    }
    
    if (exists $map{$uni})
    {
        print OUT $map{$uni},"\t";
    }else{
        print OUT "NA\t";
    }
    
    if (exists $scaffold{$n})
    {
        my $prob=sprintf "%0.3f",$scaffold{$n}/100;
        print OUT $prob,"\t";
    }else{
        print OUT "NA\t";
    }
    
    print OUT "\n";
    
}
close IN;
close OUT;
