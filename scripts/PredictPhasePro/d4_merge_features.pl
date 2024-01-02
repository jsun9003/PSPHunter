use warnings;
use strict;

my $abs_dir="/data1/juns/phase2/";
my $choice="$ARGV[0]";#"scaffold";
my $f_dir=$abs_dir."ML/PHApply/scaffold/$choice/test/";
my $o_dir=$f_dir."mergeFeature/";
mkdir $o_dir;
my %hash=();

my $in_dir=$f_dir."mergeSeq/";
open(FA, "$in_dir"."Intestinput.txt") or die "$!";
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

#FunA
$in_dir=$f_dir."mergeFun/";
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
        print "\tFunA\n";
        die;
    }
}
close FA;
close OUT;

