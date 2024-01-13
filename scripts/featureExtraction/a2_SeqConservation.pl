use warnings;
use strict;

my $abs_dir="/data1/juns/phase2/";
my $choice="";
my $clust=$ARGV[1];
my $re=$ARGV[0];

my $out_dir=$abs_dir."ML/scaffold/$choice"."_repeat"."$clust/$re/";
my $data_dir="./data_raw/";
my $hhm_dir="./data_raw/hhm/";
my $pssm_dir="./data_raw/pssm/";

my $f_dir=$out_dir;
my $o_dir=$f_dir."SeqConservation/";
mkdir $o_dir;


my %cs=();
open(FB, "$data_dir"."cscore_quantile.txt") or die "$!";
while (<FB>)
{
    chomp;
    my @item=split/\t/;
    my @ref=@item[1..$#item];
    $cs{$item[0]}=\@ref;
}
close FB;

open(FB, "$abs_dir"."Mut/NonPS/cscore_quantile.txt") or die "$!";
while (<FB>)
{
    chomp;
    my @item=split/\t/;
    my @ref=@item[1..$#item];
    $cs{$item[0]}=\@ref;
}
close FB;

open(FB, "$abs_dir"."Mut/PhaSePred/cscore_quantile.txt") or die "$!";
while (<FB>)
{
	chomp;
	my @item=split/\t/;
	my @ref=@item[1..$#item];
	$cs{$item[0]}=\@ref;
}
close FB;

open(OUT, ">$o_dir"."CVinput.txt") or die "$!";
open(FA, "$f_dir"."CV.txt") or die "$!";
while (<FA>)
{
    chomp;
    my ($uni,$tag)=split/\s+/,$_;
    my @fea=();
    my @vec1=();
    my $flag=0;
    my $bacpssm="$hhm_dir"."$uni.hhm";
    if (!-s $bacpssm) {
        print $bacpssm,"\tno bac pssm\n";
        for(my $j=1;$j<=4;$j++)
        {
            push(@vec1,0);
        }
    }else{
        @vec1=@{&hhm_feature($bacpssm)};
    }
    my @vec2=();
    $bacpssm="$pssm_dir"."$uni.pssm";
    if (!-s $bacpssm) {
        print $bacpssm,"\tno bac pssm\n";
        for(my $j=1;$j<=4;$j++)
        {
            push(@vec2,0);
        }
    }else{
        @vec2=@{&pssm_feature($bacpssm)};
    }
    my @vec3=();
    if (exists $cs{$uni})
    {
        @vec3=@{$cs{$uni}};
    }else{
        for(my $j=1;$j<=5;$j++)
        {
            push(@vec3,0);
        }
    }
    push(@fea,@vec1,@vec2,@vec3);
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
    my @fea=();
    my @vec1=();
    my $bacpssm="$hhm_dir"."$uni.hhm";
    if (!-s $bacpssm) {
        print $bacpssm,"\tno bac pssm\n";
        for(my $j=1;$j<=4;$j++)
        {
            push(@vec1,0);
        }
    }else{
        @vec1=@{&hhm_feature($bacpssm)};
    }
    my @vec2=();
    $bacpssm="$pssm_dir"."$uni.pssm";
    if (!-s $bacpssm) {
        print $bacpssm,"\tno bac pssm\n";
        for(my $j=1;$j<=4;$j++)
        {
            push(@vec2,0);
        }
    }else{
        @vec2=@{&pssm_feature($bacpssm)};
    }
    
    my @vec3=();
    if (exists $cs{$uni})
    {
        @vec3=@{$cs{$uni}};
    }else{
        for(my $j=1;$j<=5;$j++)
        {
            push(@vec3,0);
        }
    }
    push(@fea,@vec1,@vec2,@vec3);
    print OUT "$tag\t",join "\t",@fea,"\n";
}
close OUT;
close FA;

sub hhm_feature
{
    my $file=$_[0];
    #print $file,"\n";
    my @twenty=();
    my %mtx=();
    for(my $i=0;$i<20;$i++)
    {
        $twenty[$i]=0;
    }
    my $l=0;
    open(INS,"$file") or die "$!";
    chomp(my @hhm=<INS>);
	close INS;
    for(my $i=0;$i<@hhm;$i++)
    {
        last if $hhm[$i]=~/^\/\//;
        next if $i<=10;
        next if $hhm[$i]=~/^>/;
        my @item=split/\s+/,$hhm[$i];
        if ($item[0] && $item[1])
        {
            next if $item[0] eq "NULL";
            if ($item[0]=~/[A-Z]/ && $item[1]=~/[0-9]/)
            {
                my @hhm_line=@item[2..21];
                for(my $k=0;$k<=$#hhm_line;$k++)
                {
                    if($hhm_line[$k] ne "*")
                    {
                        $hhm_line[$k]=2**(-($hhm_line[$k]/1000));
                    }else
                    {
                        $hhm_line[$k]=0;
                    }
                }
                my @n_val=@{&zscore(@hhm_line)};
                $l++;
                for(my $j=0;$j<@n_val;$j++)
                {
                    $twenty[$j]+=$n_val[$j];
                }
                $i=$i+2;
            }
        }
    }
    close INS;
    
    for(my $j=0;$j<@twenty;$j++)
    {
        $twenty[$j]=sprintf "%0.3f",$twenty[$j]/$l;
    }
    my $h=$twenty[1]+$twenty[18]+$twenty[6]+$twenty[10]+$twenty[19]+$twenty[4]+$twenty[7]+$twenty[17]+$twenty[0]+$twenty[9];
    my $c=$twenty[3]+$twenty[14]+$twenty[8]+$twenty[2];
    my $p=$twenty[15]+$twenty[13]+$twenty[3]+$twenty[16]+$twenty[11];
    my $g=$twenty[5]+$twenty[12];
    my @res=($h,$c,$p,$g);
    return \@res;
}

sub pssm_feature
{
    my $file=$_[0];
    #print $file,"\n";
    my @twenty=();
    my %mtx=();
    for(my $i=0;$i<20;$i++)
    {
        $twenty[$i]=0;
    }
    my $l=0;
    open(INS,"$file") or die "$!";
    while (<INS>)
    {
        #chomp;
        my $raw=$_;
        $raw=~s/^\s+//g;
        my @raw_inf=split/ +/,$raw;
        my $num=scalar(@raw_inf);
        #print $num,"\n";
        if($num>42)
        {
            my @val=@raw_inf[2..21];
            my @n_val=@{&zscore(@val)};
            $l++;
            for(my $j=0;$j<@n_val;$j++)
            {
                $twenty[$j]+=$n_val[$j];
            }
        }
    }
    close INS;
    
    for(my $j=0;$j<@twenty;$j++)
    {
        $twenty[$j]=sprintf "%0.3f",$twenty[$j]/$l;
    }
    my $h=$twenty[4]+$twenty[17]+$twenty[8]+$twenty[12]+$twenty[18]+$twenty[13]+$twenty[9]+$twenty[19]+$twenty[0]+$twenty[10];
    my $c=$twenty[6]+$twenty[1]+$twenty[11]+$twenty[3];
    my $p=$twenty[15]+$twenty[5]+$twenty[16]+$twenty[2];
    my $g=$twenty[7]+$twenty[14];
    my @res=($h,$c,$p,$g);
    return \@res;
}

sub zscore
{
    my @tempsub=@_;
    my $total=0;
    foreach my $item (@tempsub)
    {
        $total+=$item;
    }
    my $avgsub=$total/scalar @tempsub;
    
    my $sdsub=0;
    foreach my $item (@tempsub){
        $sdsub+=($item-$avgsub)**2;
    }
    my $stdsub=sqrt($sdsub/($#tempsub));
   # print $avgsub,"\t$stdsub\n";
    my @norm=();
    foreach my $item (@tempsub)
    {
        my $norm=0;
        if ($stdsub ==0)
        {
            $norm=0;
        }else{
            $norm=($item-$avgsub)/$stdsub;
            #$norm=sprintf"%.4f",1/(1+exp(-$norm));
        }
        push(@norm,$norm);
    }
    return \@norm;
}



