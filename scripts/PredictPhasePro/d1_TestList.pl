use warnings;
use strict;

my $abs_dir="/data1/juns/phase2/";
my $choice="$ARGV[0]";#"scaffold";
my $out_dir=$abs_dir."ML/PHApply/scaffold/$choice/test/";
mkdir $out_dir;

my @human=();
open(IN, "$abs_dir"."data_raw/integrate/anno_list.txt") or die "$!";
while (<IN>)
{
    chomp;
    my ($uni,$tag)=split;
    push(@human,$uni) unless $_~~@human;
}
close IN;

print scalar @human,"\n";
open(OUT,">$out_dir"."Intest.txt") or die "$!";
foreach my $human (@human)
{
    print OUT "$human\t0\n";
}
close OUT;
