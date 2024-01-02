use warnings;
use strict;

my $abs_dir="/data1/juns/phase2/";
my $re=$ARGV[0];
my $choice="$ARGV[1]";
my $python=$abs_dir."program/PHApply/Intest-apply-test.py";
#my $f1_dir=$abs_dir."ML/PHApply/scaffold/$choice/train/$re/mergeSeq/";
#my $f2_dir=$abs_dir."ML/PHApply/scaffold/$choice/test/mergeSeq/";
#system "python $python $f1_dir $f2_dir $re";

#$f1_dir=$abs_dir."ML/PHApply/$choice/train/$re/mergeFun/";
#$f2_dir=$abs_dir."ML/PHApply/$choice/test/mergeFun/";
#system "python $python $f1_dir $f2_dir $re";

my $f1_dir=$abs_dir."ML/PHApply/scaffold/$choice/train/$re/mergeFeature/";
my $f2_dir=$abs_dir."ML/PHApply/scaffold/$choice/test/mergeFeature/";
system "python $python $f1_dir $f2_dir $re";
