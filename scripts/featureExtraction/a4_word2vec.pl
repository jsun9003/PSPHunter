use warnings;
use strict;

my $abs_dir="/data1/juns/phase2/";
my $choice="PhaSePredMix";
my $clust=$ARGV[1];
my $re=$ARGV[0];

my $out_dir=$abs_dir."ML/scaffold/$choice"."_repeat"."$clust/$re/";
my $fasta_dir="/data1/juns/PathHost/data_raw/fasta_pathhost/";

#topo
my %wv=();
my $win="70";
my $size="60";
open(IN, "/data1/juns/word2v/uniprot_sprot"."$win"."_size$size.txt") or die "$!";
my $c=0;
while (<IN>)
{
    chomp;
    next if /^Uniprot/;
    my @item=split/\t/;
    my @ref=@item[1..$#item];
    $wv{$item[0]}=\@ref;
}
close IN;

my $f_dir=$out_dir;
my $o_dir=$f_dir."word2vec$win"."_$size/";
print $o_dir,"\n";

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
        for(my $j=1;$j<=$size;$j++) #39
        {
            push(@vec,0);
        }
    }else{
        open(F,$bacfasta) or die "$!";
        chomp(my @bacteria_seq=<F>);
        close F;
        my @aa=split//,$bacteria_seq[1];
        for(my $i=0;$i<(scalar @aa)-2;$i++)
        {
            my $tri=$aa[$i].$aa[$i+1].$aa[$i+2];
            if (exists $wv{$tri})
            {
                my @tvec=@{$wv{$tri}};
                for(my $j=0;$j<@tvec;$j++)
                {
                    $vec[$j]+=$tvec[$j];
                }
            }
            
        }
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
        print $bacfasta,"\t","no bac protein\n";
        for(my $j=1;$j<=$size;$j++) #39
        {
            push(@vec,0);
        }
    }else{
        open(F,$bacfasta) or die "$!";
        chomp(my @bacteria_seq=<F>);
        close F;
        my @aa=split//,$bacteria_seq[1];
        for(my $i=0;$i<(scalar @aa)-2;$i++)
        {
            my $tri=$aa[$i].$aa[$i+1].$aa[$i+2];
            if (exists $wv{$tri})
            {
                my @tvec=@{$wv{$tri}};
                for(my $j=0;$j<@tvec;$j++)
                {
                    $vec[$j]+=$tvec[$j];
                }
            }
            
        }
    }
    my @fea=();
    push(@fea,@vec);
    print OUT "$tag\t",join "\t",@fea,"\n";
}
close OUT;
close FA;

