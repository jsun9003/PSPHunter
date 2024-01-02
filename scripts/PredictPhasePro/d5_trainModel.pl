use warnings;
use strict;

my $abs_dir="/data1/juns/phase2/";
my $re=$ARGV[0];
my $choice="$ARGV[1]";#"scaffold";
my $model="mergeFeature";#mergeFeatureNew
my $dir=$abs_dir."ML/scaffold/$choice"."_repeat100/$re/";

my $o_dir=$abs_dir."ML/PHApply/scaffold/$choice/train/$re/";
system "mkdir -p $o_dir";

my $train_dir=$o_dir."$model/";
mkdir $train_dir;

my %hash=();
my $n=0;
open(IN,"$dir"."CV.txt") or die "$!";
while (<IN>)
{
    chomp;
    $n++;
    $hash{$n}=$_;
}
close IN;

open(IN,"$dir"."Intest.txt") or die "$!";
while (<IN>)
{
    chomp;
    $n++;
    $hash{$n}=$_;
}
close IN;

open(OUT,">$train_dir"."list.txt") or die "$!";
foreach my $key (sort {$a<=>$b} keys %hash)
{
    print OUT $hash{$key},"\n";
}
close OUT;

open(OUT,">$train_dir"."CVinput.txt") or die "$!";
open(IN,"$dir"."$model/CVinput.txt") or die "$!";
while (<IN>)
{
    chomp;
    print OUT $_,"\n";
}
close IN;

open(IN,"$dir"."$model/Intestinput.txt") or die "$!";
while (<IN>)
{
    chomp;
    print OUT $_,"\n";
}
close IN;
close OUT;

my $python=$abs_dir."program/PHApply/Intest-apply-train.py";
my $f_dir=$abs_dir."ML/PHApply/scaffold/$choice/train/$re/$model/";

system "python $python $f_dir";
