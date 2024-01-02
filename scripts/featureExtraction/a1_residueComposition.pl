use warnings;
use strict;

my $abs_dir="/data1/juns/phase2/";
my $choice="PhaSePredMix";
my $clust=$ARGV[1];
my $re=$ARGV[0];

my $out_dir=$abs_dir."ML/scaffold/"."$choice"."_repeat"."$clust/$re/";
my $data_dir="/data1/juns/PathHost/data_raw/";
my $fasta_dir="/data1/juns/PathHost/data_raw/fasta_pathhost/";
my $hhm_dir="/data1/juns/PathHost/data_raw/hhm_pathhost/";
my $pssm_dir="/data1/juns/PathHost/data_raw/pssm_pathhost/";
my %class=(
            "C"=>"1","W"=>"1","H"=>"1","M"=>"1","Y"=>"1","F"=>"1","I"=>"1","V"=>"1","A"=>"1","L"=>"1",
            "G"=>"2","P"=>"2",
            "S"=>"3","Q"=>"3","T"=>"3","N"=>"3",
            "E"=>"4","R"=>"4","K"=>"4","D"=>"4"
           );
my $f_dir=$out_dir;
my $o_dir=$f_dir."ResComposition/";
mkdir $o_dir;

open(OUT, ">$o_dir"."CVinput.txt") or die "$!";
open(FA, "$f_dir"."CV.txt") or die "$!";
while (<FA>)
{
    chomp;
    my ($uni,$tag)=split/\s+/,$_;
    my @vec=();
    my $bacfasta="$fasta_dir"."$uni.fasta";
    if (!-s $bacfasta) {
        print $bacfasta,"\t","no bac protein\n";
        for(my $j=1;$j<=20;$j++) #39
        {
            push(@vec,0);
        }
    }else{
        open(F,$bacfasta) or die "$!";
        chomp(my @bacteria_seq=<F>);
        close F;
        my $vector1=&sequence_feature1($bacteria_seq[1]);
        $vector1=~s/^\s+//;
        my @vec1=split/\t/,$vector1;
        my $vector2=&sequence_feature2($bacteria_seq[1]);
        $vector2=~s/^\s+//;
        my @vec2=split/\t/,$vector2;

        push(@vec,@vec1,@vec2);
    }
    my @fea=();
    push(@fea,@vec);
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
    my @vec=();
    my $bacfasta="$fasta_dir"."$uni.fasta";
    if (!-s $bacfasta) {
        print $bacfasta;
        for(my $j=1;$j<=20;$j++)
        {
            push(@vec,0);
        }
    }else{
        open(F,$bacfasta) or die "$!";
        chomp(my @bacteria_seq=<F>);
        close F;
        my $vector1=&sequence_feature1($bacteria_seq[1]);
        $vector1=~s/^\s+//;
        my @vec1=split/\t/,$vector1;
        my $vector2=&sequence_feature2($bacteria_seq[1]);
        $vector2=~s/^\s+//;
        my @vec2=split/\t/,$vector2;
        push(@vec,@vec1,@vec2);
    }

    my @fea=();
    push(@fea,@vec);
    print OUT "$tag\t",join "\t",@fea,"\n";
}
close OUT;
close FA;

#feature
sub sequence_feature2
{
    my %freq=();
    for(my $i=1;$i<=4;$i++)
    {
        for(my $j=1;$j<=4;$j++)
        {
            my $ref=$i.$j;
            $freq{$ref}=0;
        }
    }
    my $sequnce=$_[0];
    my @elements=split//,$sequnce;
    my $new_seq="";
    for(my $i=0;$i<@elements;$i++)
    {
        if (exists $class{$elements[$i]})
        {
            $new_seq=$new_seq.$class{$elements[$i]};
        }
    }
    
    
    for(my $i=0;$i<(length $new_seq)-1;$i++)
    {
        my $triple=substr($new_seq,$i,2);
        $freq{$triple}++;
    }
    
    #my @vector=();
    my $max=0;
    my $min=0;
    for(my $i=1;$i<=4;$i++)
    {
        for(my $j=1;$j<=4;$j++)
        {
            my $ref=$i.$j;
            my $frequncy=$freq{$ref};
            if ($frequncy>$max)
            {
                $max=$frequncy;
            }
            if ($frequncy<$min)
            {
                $min=$frequncy;
            }
        }
    }
    #print $new_seq,"\t$min\t$max\n";
    my $result="";
    for(my $i=1;$i<=4;$i++)
    {
        for(my $j=1;$j<=4;$j++)
        {
            my $ref=$i.$j;
            my $frequncy=sprintf "%0.3f",($freq{$ref}-$min)/$max;
            $result=$result."\t".$frequncy;
        }
    }
    return $result;
}

sub sequence_feature1
{
    my %freq=();
    for(my $i=1;$i<=4;$i++)
    {
        $freq{$i}=0;
    }
    my $sequnce=$_[0];
    my @elements=split//,$sequnce;
    my $new_seq="";
    for(my $i=0;$i<@elements;$i++)
    {
        if (exists $class{$elements[$i]})
        {
            $new_seq=$new_seq.$class{$elements[$i]};
        }
    }
    for(my $i=0;$i<(length $new_seq);$i++)
    {
        my $triple=substr($new_seq,$i,1);
        $freq{$triple}++;
    }
    my $len=length $sequnce;
    my $result="";
    for(my $i=1;$i<=4;$i++)
    {
        my $ref=$i;
        my $frequncy=sprintf "%0.3f",$freq{$ref}/$len;
        $result=$result."\t".$frequncy;
    }
    return $result;
}
