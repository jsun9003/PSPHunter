use warnings;
use strict;

my $abs_dir="/data1/juns/phase2/";
my $choice="PhaSePredMix";
my $clust=$ARGV[1];
my $re=$ARGV[0];

my $out_dir=$abs_dir."ML/scaffold/".$choice."_repeat"."$clust/$re/";

my $f_dir=$out_dir;
my $o_dir=$f_dir."mergeSeq/";
mkdir $o_dir;
my %hash=();
#777
my $in_dir=$f_dir."ResComposition/";
open(FA, "$in_dir"."CVinput.txt") or die "$!";
my $n=0;
while (<FA>)
{
    chomp;
    s/\s+$//g;
    $n++;
    my @item=split/\s+/;
    my @ref=@item[1..$#item];
    $hash{$n}=\@ref;
}
close FA;

$in_dir=$f_dir."SeqConservation/";
open(FA, "$in_dir"."CVinput.txt") or die "$!";
$n=0;
while (<FA>)
{
    chomp;
    s/\s+$//g;
    $n++;
    my @item=split/\s+/;
    if (exists $hash{$n})
    {
        my @ref=@{$hash{$n}};
        push(@ref,@item[1..$#item]);
        $hash{$n}=\@ref;
    }else{
        print $n,"\tSeqConservation\n";
        die;
    }
}
close FA;

$in_dir=$f_dir."FunSite/";
open(FA, "$in_dir"."CVinput.txt") or die "$!";
$n=0;
while (<FA>)
{
    chomp;
    s/\s+$//g;
    $n++;
    my @item=split/\s+/;
    if (exists $hash{$n})
    {
        my @ref=@{$hash{$n}};
        push(@ref,@item[1..$#item]);
        $hash{$n}=\@ref;
    }else{
        print $n,"\tReEntropyDis\n";
        die;
    }
}
close FA;

$in_dir=$f_dir."word2vec70_60/";
open(OUT, ">$o_dir"."CVinput.txt") or die "$!";
open(FA, "$in_dir"."CVinput.txt") or die "$!";
$n=0;
while (<FA>)
{
    chomp;
    s/\s+$//g;
    $n++;
    my @item=split/\s+/;
    if (exists $hash{$n})
    {
        my @ref=@{$hash{$n}};
        push(@ref,@item[1..$#item]);
        print OUT $item[0],"\t",join "\t",@ref,"\n";
    }else{
        print "\tDisPTM\n";
        die;
    }
}
close FA;
close OUT;

%hash=();
$in_dir=$f_dir."ResComposition/";
open(FA, "$in_dir"."Intestinput.txt") or die "$!";
$n=0;
while (<FA>)
{
    chomp;
    s/\s+$//g;
    $n++;
    my @item=split/\s+/;
    my @ref=@item[1..$#item];
    $hash{$n}=\@ref;
}
close FA;

$in_dir=$f_dir."SeqConservation/";
open(FA, "$in_dir"."Intestinput.txt") or die "$!";
$n=0;
while (<FA>)
{
    chomp;
    s/\s+$//g;
    $n++;
    my @item=split/\s+/;
    if (exists $hash{$n})
    {
        my @ref=@{$hash{$n}};
        push(@ref,@item[1..$#item]);
        $hash{$n}=\@ref;
    }else{
        print $n,"\tSeqConservation\n";
        die;
    }
}
close FA;

$in_dir=$f_dir."FunSite/";
open(FA, "$in_dir"."Intestinput.txt") or die "$!";
$n=0;
while (<FA>)
{
    chomp;
    s/\s+$//g;
    $n++;
    my @item=split/\s+/;
    if (exists $hash{$n})
    {
        my @ref=@{$hash{$n}};
        push(@ref,@item[1..$#item]);
        $hash{$n}=\@ref;
    }else{
        print $n,"\tSsDisNARPs\n";
        die;
    }
}
close FA;

$in_dir=$f_dir."word2vec70_60/";
open(OUT, ">$o_dir"."Intestinput.txt") or die "$!";
open(FA, "$in_dir"."Intestinput.txt") or die "$!";
$n=0;
while (<FA>)
{
    chomp;
    s/\s+$//g;
    $n++;
    my @item=split/\s+/;
    if (exists $hash{$n})
    {
        my @ref=@{$hash{$n}};
        push(@ref,@item[1..$#item]);
        print OUT $item[0],"\t",join "\t",@ref,"\n";
    }else{
        print "\tDisPTM\n";
        die;
    }
}
close FA;
close OUT;
