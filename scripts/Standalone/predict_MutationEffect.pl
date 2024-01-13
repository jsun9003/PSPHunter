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
    system "/home/PSPHunter/software/miniconda3/bin/perl $sendmail_filetome $jsun_email $id 3";
}

my @aa=qw(C W H M Y F I V A L G P S Q T N E R K D);

####Split fasta into different truncations
my $seq_fa=$dir."temp.fasta";

my $fasta_dir=$dir."fasta3/";
mkdir $fasta_dir;

my $seq="";
open(IN,$seq_fa) or die "$!";
while (<IN>)
{
    chomp;
    print $_,"\tRaw seq Each Line\n";
    if($_ =~/>/)
    {
       print $_,"\tHeader sequence Infor\n";
       next;
    }else{
       print $seq,"\tSeq\n";
       $_=~s/\n|\r|\s//g;
       $seq=$seq."".$_;
    }
}

close IN;
#my $seq=$lines[1];
$seq=~s/-|\*//g;

my $len=length $seq;
$seq=uc $seq;
print $seq,"\tfinal seq\n";

my @list=();
my @seq=split//,$seq;
my %index=();
for(my $i=0;$i<@seq;$i++)
{
    if ($seq[$i]~~@aa)
    {
         $index{$i+1}=$seq[$i];
        for(my $j=0;$j<@aa;$j++)
        {
            my $raw=$seq;
            substr($raw,$i,1,$aa[$j]);
            my $temp="query"."_$i"."_$aa[$j]";
            push(@list,$temp);
            open(OUT,">$fasta_dir"."$temp.fasta") or die "$!";
            print OUT ">query"."_$i"."_$aa[$j]","\n$raw\\n";
            close OUT;
        }
    }
}

my $ML_dir=$dir."ML3/";
mkdir $ML_dir;

open(OUT,">$ML_dir"."Intest.txt") or die "$!";
foreach my $human (@list)
{
    print OUT "$human\t0\n";
}
close OUT;

###word2vec hash
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

my $o_dir=$ML_dir."word2vec/";
mkdir $o_dir;

my $size=60;
open(OUT, ">$o_dir"."Intestinput.txt") or die "$!";
open(FA, "$ML_dir"."Intest.txt") or die "$!";
while (<FA>)
{
    chomp;
    my ($uni,$tag)=split/\s+/,$_;
    my @fea=();
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
    push(@fea,@vec);
    print OUT "$tag\t",join "\t",@fea,"\n";
}
close OUT;
close FA;

###Model Test
my $python=$progam_dir."Intest-apply-test.py";
for(my $i=1;$i<=100;$i++)
{
    my $f1_dir=$home_dir."/train/$i/word2vec70_60/";
    my $f2_dir=$o_dir;
    system "/home/PSPHunter/software/miniconda3/bin/python $python $f1_dir $f2_dir $i";
}

#Avg Prob

open(IN,"$ML_dir"."Intest.txt") or die "$!";
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

open(OUT,">$dir"."InProbAVG.txt") or die "$!";
my %hash=();
for(my $i=1;$i<=100;$i++)
{
    open(IN,"$o_dir"."InProbphase$i.txt") or die "$!";
    my $k=0;
    while (<IN>)
    {
        chomp;
        $k++;
        if (exists $hash{$k})
        {
            $hash{$k}+=$_;
        }else{
            $hash{$k}=$_;
        }
    }
    close IN;
}

my %prob=();
foreach my $key (sort {$hash{$b}<=>$hash{$a}} keys %hash)
{
    my $prob=sprintf "%0.3f",$hash{$key}/100;
    my ($query,$index,$aa)=split/_/,$tag{$key};
    $prob{$index+1}{$aa}=$prob;
    #print OUT $key,"\t$prob\n";
}

print OUT "SeqIndex\tAA\t";
foreach my $aa (sort @aa)
{
    print OUT $aa,"\t";
}
print OUT "\n";

foreach my $key1(sort {$a<=>$b} keys %prob)
{
    print OUT $key1,"\t$index{$key1}\t";
    foreach my $key2 (sort keys %{$prob{$key1}} )
    {
        print OUT $prob{$key1}{$key2},"\t";
    }
    print OUT "\n";
}
close OUT;

my $sendmail=$progam_dir."sendmail_myemail.py";

open(DO,">$dir"."done_predictMutEffect.txt")||die"$!";
foreach my $key(keys %done){
  print DO $key,",",$done{$key},"\n";
}
close DO;

print $email,"\tEmail To\n";
if (-e $emailsend_file_check)
{
    print "you will get an email latter;"
}else{
    system "/home/PSPHunter/software/miniconda3/bin/perl $sendmail_file $email $id 3";
}

system "rm -rf $fasta_dir";
system "rm -rf $ML_dir";

#system "zip -j $dir/PSPHunter_results.zip $dir/*.out $dir/*.pdb";
