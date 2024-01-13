#!/usr/bin/perl
my ($id,$email)=@ARGV;
#my $identi=0.3;

my $home_dir="/home/PSPHunter/";
my $dir="/var/www/psphunter.stemcellding.org/job/$id/";
#my $dir="$home_dir"."catalytic/job/";
my $progam_dir=$home_dir."scripts/";
#my $data_raw="$home_dir"."catalytic/data_raw/";
my $sendmail_file="/var/www/psphunter.stemcellding.org/sendmail.pl";
my $sendmail_filetome="/var/www/psphunter.stemcellding.org/sendmailtome.pl";
my $emailsend_file_check ="$dir"."emailcheck.txt";

my $jsun_email="o.sj\@qq.com";
my $own_flag=0;
if (-e $emailsend_file_check)
{
    print "you will get an email latter;";
    $own_flag++;
}else{
    system "/home/PSPHunter/software/miniconda3/bin/perl $sendmail_filetome $jsun_email $id 1";
}

my %wv=();
open(IN, "$home_dir"."wordvec/uniprot_sprot70_size60.txt") or die "$!";
while (<IN>)
{
    chomp;
    next if /^Uniprot/;
    my @item=split/\t/;
    my @ref=@item[1..$#item];
    $wv{$item[0]}=\@ref;
}
close IN;

####Split fasta into different truncations
my $seq_fa=$dir."temp.fasta";
my $fasta_dir=$dir."fasta1/";
mkdir $fasta_dir;

open(IN,"$seq_fa") or die "$!";
chomp(my @fastas=<IN>);
close IN;

my %fa=();
my $pro="";
my $flag=0;
foreach my $fasta (@fastas)
{
    $fasta=~s/\r//g;
    if ($fasta=~/^>(.*)/)
    {
        $pro=$1;
        $pro=~s/\s+|\*|\.|\[|\]|\"|\||\(|\)|\=//g;
        $pro=substr($pro,0,15);
        open(OUT,">>$fasta_dir"."$pro.fasta") or die "$!";
        print OUT ">",$pro,"\n";
        $fa{$pro}++;
        $flag++;
        close OUT;
    }else{
        open(OUT,">>$fasta_dir"."$pro.fasta") or die "$!";
        $fasta=uc $fasta;
        print OUT "$fasta";
        close OUT;
    }
    
}
print scalar keys %fa,"\n";
my $proNo=scalar keys %fa;

my $o_dir=$dir."ML1/";
mkdir $o_dir;

my $size=60;

if ($flag==0)
{
    $proNo=1;
    $fa{"query"}++;
    open(OUT,">$fasta_dir"."query.fasta") or die "$!";
    print OUT ">query\n";
    foreach my $fasta (@fastas)
    {
        $fasta=uc $fasta;
        print OUT "$fasta";
    }
    close OUT;
}


if ($proNo == 1)
{
    $fa{"testSun"}++;
    my $seqtest="/var/www/psphunter.stemcellding.org/job/testSun.fasta";
    system("cp $seqtest $fasta_dir");
    open(OUT, ">$o_dir"."Intestinput.txt") or die "$!";
    open(OU, ">$o_dir"."Intest.txt") or die "$!";
    foreach my $uni (sort keys %fa)
    {
        my @fea=();
        my @vec=();
        my $bacfasta="$fasta_dir"."$uni.fasta";
         if (!-s $bacfasta) {
            print $bacfasta,"\t","no bac protein\n";
            for(my $j=1;$j<=$size;$j++) 
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
        push(@fea,@vec);
        print OUT "0\t",join "\t",@fea,"\n";
        print OU "$uni\n";
    }
    close OUT;
    close OU;
}else{
    open(OUT, ">$o_dir"."Intestinput.txt") or die "$!";
    open(OU, ">$o_dir"."Intest.txt") or die "$!";
    foreach my $uni (sort keys %fa)
    {
        my @fea=();
        my @vec=();
        my $bacfasta="$fasta_dir"."$uni.fasta";
         if (!-s $bacfasta) {
            print $bacfasta,"\t","no bac protein\n";
            for(my $j=1;$j<=$size;$j++) 
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
        push(@fea,@vec);
        print OUT "0\t",join "\t",@fea,"\n";
        print OU "$uni\n";
    }
    close OUT;
    close OU;
}

my $python=$progam_dir."Intest-apply-test.py";
for (my $re=1;$re<=100;$re++)
{
    my $f1_dir=$home_dir."/train/$re/word2vec70_60/";
    my $f2_dir=$o_dir;
    system "/home/PSPHunter/software/miniconda3/bin/python $python $f1_dir $f2_dir $re";
}

my $n=0;
my %name=();
open(IN,"$o_dir"."Intest.txt") or die "$!";
while (<IN>)
{
    chomp;
    $n++;
    $name{$n}=$_;
}
close IN;

my %prob=();
for (my $re=1;$re<=100;$re++)
{
    open(IN,"$o_dir"."InProbphase$re.txt") or die "$!";
    my $m=0;
    while (<IN>)
    {
        chomp;
        $m++;
        $prob{$m}+=$_;
    }
    close IN;
}

if ($proNo==1)
{
    open(OUT,">$dir"."Avg.txt") or die "$!";
    print OUT "ProteinName\tProb\n";
    foreach my $name (sort {$a<=>$b} keys %name)
    {
        next if $name{$name} eq "testSun";
        print OUT $name{$name},"\t",$prob{$name}/100,"\n";
    }
    close OUT; 
}else{
    open(OUT,">$dir"."Avg.txt") or die "$!";
    print OUT "ProteinName\tProb\n";
    foreach my $name (sort {$a<=>$b} keys %name)
    {
        print OUT $name{$name},"\t",$prob{$name}/100,"\n";
    }
    close OUT; 
}



my $sendmail=$progam_dir."sendmail_myemail.py";

open(DO,">$dir"."done_predictProProb.txt")||die"$!";
foreach my $key(keys %done){
  print DO $key,",",$done{$key},"\n";
}
close DO;

print $email,"\tEmail To\n";
if (-e $emailsend_file_check)
{
    print "you will get an email latter;"
}else{
    system "/home/PSPHunter/software/miniconda3/bin/perl $sendmail_file $email $id 1";
}

system "rm -rf $fasta_dir";
