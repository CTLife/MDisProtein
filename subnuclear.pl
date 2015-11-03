
#最后一次修改于2012年12月11日。

#严格筛选，分为9类。
#忽略核外位置。               
#建立亚核蛋白数据集，选出单定位和多定位。
#要求准备好从UniProt下载的flat text与fasta格式的序列，2个文件只有后缀不同。

#步骤如下：
##############################################################################################################################
# 01  去掉注释信息中含有多个“SUBCELLULAR LOCATION:”的序列。
#     (此步会把亚细胞定位信息中含有ISOFORM的部分序列除去，只含一个‘SUBCELLULAR LOCATION: Isoform’不会被去掉。) 
#     (此步去掉的序列需要单独处理。)                                                                                                               
# 02  去除序列的多余注释信息。(序列数不会改变)    (留下了5种信息)                                                  
# 03  去掉NOTE信息。        (序列数不会改变)       
# 04  去掉亚细胞定位信息模糊的序列。（nucleus([PB])）                                                                                                                                                                                                                                                          
# 05  去掉片段。                                                                                               
# 06  去掉序列中含有UOBZJX的序列。                                                                                  
# 07  去掉短序列。（小于50AA的去掉）                                                                                                                                                                                                         
# 08  进一步去掉多余注释信息。(留下了3种信息：ID，CC，//)   (序列数不会改变)                                                                                                                                                   
# 09  选出各位置的序列(单定位和多定位)，分别输出位于各亚细胞核位置的序列。（注释了定位在核外的序列一律不要）                                                                                                                                                                        
# 10  提取出含全部注释信息的序列。       
# 11  提取fasta格式的序列。                                                                                                                                                                      
###############################################################################################################################
#所有结果都放在了一个新建的文件夹下，各步骤的结果放在各个子文件夹下。




#!/usr/bin/env perl
use strict;
use warnings;

my $script_name = $0;                            #此perl脚本的名字。
$script_name =~ s/\.pl//   or die  "$!";         #去掉扩展名。
mkdir  "$script_name"   or  die "$!" ;           #在当前目录(此程序文件所在目录)下创建一个文件夹，用于存放此程序的所有输出结果。
  


#以下是4个函数，供其它函数调用，用于统计序列条数或某种模式的出现次数。
##########################################################################################################################################
sub NUMBER_TEXT    #统计序列条数(flat text format)
##########################################################################################################################################
{
my $name = $_[0];          #第1个参数,待统计文件的名称(含路径)。
my $path1 = $_[1];         #第2个参数,结果文件存放的路径。
open(FILE1, "<", "$name")   or die "$!";                            #以读的方式打开待处理文件。
open(FILE2, ">>", "$path1/NUMBER_TEXT.txt")   or die "$!";          #以追加的方式打开存放结果的文件。
my $n1=0;
my $n2=0;

while (my $line=<FILE1>) {
  if ($line=~m/^ID/) {$n1++;}
  if ($line=~m"^//") {$n2++;}
}

print  FILE2  "
文件$name中的序列条数:以ID为标志有$n1条序列。
文件$name中的序列条数:以//为标志有$n2条序列。
\n\n";
close FILE1;
close FILE2;
}#################################################### END NUMBER_TEXT





##############################################################################################################
sub NUMBER_FASTA      #统计序列条数(fasta format)
##############################################################################################################
{ 
my $name = $_[0];          #第1个参数,待统计文件的名称(含路径)。
my $path1 = $_[1];         #第2个参数,结果文件存放的路径。
open(FILE1, "<", "$name")   or die "$!";                             #以读的方式打开待处理文件。
open(FILE2, ">>", "$path1/NUMBER_FASTA.txt")   or die "$!";          #以追加的方式打开存放结果的文件。
my $n1=0;

while (my $line=<FILE1>) {
  if ($line=~m/^>/) {$n1++;}
}

print  FILE2  "文件$name中的序列条数:以>为标志有$n1条序列。 \n\n";
close FILE1;
close FILE2;
}#################################################### END NUMBER_FASTA 





##############################################################################################################
sub PATTERN_SEQ      #统计某种模式的出现次数，一条序列统计一次。
##############################################################################################################
{
my $name=$_[0];     #第一个参数，待处理文件的名字(含路径)。
my $path1 = $_[1];  #第2个参数,结果文件存放的路径。
my $pattern=$_[2];  #第3个参数为待统计的模式。(传递参数时，用单引号把内容括起来)
open(FILE1, "<", "$name") or die "$!";                        #以读的方式打开待处理文件。
open(FILE2, ">>", "$path1/PATTERN_SEQ.txt") or die "$!";      #以追加的方式打开存放结果的文件。
my $n1=0;

while (my $line=<FILE1>) {
  my $one='';
  $line=~m/^ID/  or die "$!";        #把一条序列的所有信息放在一个标量变量里。
  do {
    $one=$one.$line;
    $line=<FILE1>;  #读取下一行。
  }until $line=~m"^//";
  $one=$one.$line;        
  if ($one=~m/$pattern/) {$n1++;}
}

print  FILE2  "文件$name中模式$pattern的出现次数（一条序列统计一次）：$n1  \n\n";
close FILE1;
close FILE2;
}##################################### END PATTERN_SEQ





##############################################################################################################
sub PATTERN_LINE  #统计某种模式的出现次数，一行统计一次。
##############################################################################################################
{
my $name=$_[0];     #第一个参数，待处理文件的名字(含路径)。
my $path1 = $_[1];  #第2个参数,结果文件存放的路径。
my $pattern=$_[2];  #第3个参数为待统计的模式。(传递参数时，用单引号把内容括起来)
open(FILE1, "<", "$name") or die "$!";                     #以读的方式打开待处理文件。
open(FILE2, ">>", "$path1/PATTERN_LINE.txt") or die "$!";  #以追加的方式打开存放结果的文件。
my $n1=0;

while (my $line=<FILE1>) {
  if ($line=~m/$pattern/) {$n1++;}
}

print  FILE2  "文件$name中模式$pattern的出现次数（一行统计一次。）：$n1  \n\n\n";
close FILE1;
close FILE2;
}########################################### END PATTERN_LINE 









#正式开始处理文件的11个步骤：
################################################################################################################################################
sub ISOFORM_1    #去掉注释信息中含有多个(>=2)“SUBCELLULAR LOCATION: Isoform”的序列。 亚细胞定位信息中含有ISOFORM的序列会被部分除去。
################################################################################################################################################
{
my $name1=$_[0];                        #待处理文件的名字(含路径)。
my $name2="1.txt";                      #用于存放结果的文件的名称。
my $name3="no_1.txt";                   #用于存放其它结果的文件的名称。
my $file_n = "1_ISOFORM";
mkdir  "$script_name/$file_n"  or die;  #建立一个子文件夹，存放这个子程序的输出结果。  

open(FILE1, "<", "$name1") or die "$!";                        #以读的方式打开待处理的文件。
open(FILE2, ">", "$script_name/$file_n/$name2") or die "$!";   #以写的方式打开存放结果的文件。
open(FILE3, ">", "$script_name/$file_n/$name3") or die "$!";   #以写的方式打开存放其它结果的文件。

while (my $line1=<FILE1>) {
  my $one='';   #用来存储一条序列的所有信息。
  $line1=~m/^ID/  or die "$!";
  do {
    $one=$one.$line1;
    $line1=<FILE1>;
  }until ($line1=~m"^//");
  $one=$one.$line1;
  if ($one=~m/(SUBCELLULAR LOCATION:)([\d\D]+)(SUBCELLULAR LOCATION:)/) {
    print FILE3 $one; 
  }else{ 
    print FILE2 $one; 
  } 
}

close FILE1;
close FILE2;
close FILE3;
&NUMBER_TEXT("$name1",                          "$script_name/$file_n");
&NUMBER_TEXT("$script_name/$file_n/$name2",     "$script_name/$file_n");
&NUMBER_TEXT("$script_name/$file_n/$name3",     "$script_name/$file_n");

&PATTERN_SEQ("$name1",                          "$script_name/$file_n",     '(SUBCELLULAR LOCATION:)([\d\D]+)(SUBCELLULAR LOCATION:)');
&PATTERN_SEQ("$script_name/$file_n/$name2",     "$script_name/$file_n",     '(SUBCELLULAR LOCATION:)([\d\D]+)(SUBCELLULAR LOCATION:)');
&PATTERN_SEQ("$script_name/$file_n/$name3",     "$script_name/$file_n",     '(SUBCELLULAR LOCATION:)([\d\D]+)(SUBCELLULAR LOCATION:)');

&PATTERN_LINE("$name1",                          "$script_name/$file_n",     '(SUBCELLULAR LOCATION:)');
&PATTERN_LINE("$script_name/$file_n/$name2",     "$script_name/$file_n",     '(SUBCELLULAR LOCATION:)');
&PATTERN_LINE("$script_name/$file_n/$name3",     "$script_name/$file_n",     '(SUBCELLULAR LOCATION:)');
}######################################### END ISOFORM_1






##########################################################################################################################################
sub FIVE_2    #只输出5类信息。剔除多余信息，只输出ID,DE,CC(SUBCELLULAR LOCATION),SQ,// 这5种信息。
##########################################################################################################################################
{ 
my $name1=$_[0];                        #待处理文件的名字(含路径)。
my $name2="2.txt";                      #用于存放结果的文件的名称。
my $name3="no_2.txt";                   #用于存放其它结果的文件的名称。
my $file_n = "2_FIVE";
mkdir  "$script_name/$file_n"  or die;  #建立一个子文件夹，存放这个子程序的输出结果。  

open(FILE1, "<", "$name1") or die "$!";                        #以读的方式打开待处理的文件。
open(FILE2, ">", "$script_name/$file_n/$name2") or die "$!";   #以写的方式打开存放结果的文件。
open(FILE3, ">", "$script_name/$file_n/$name3") or die "$!";   #以写的方式打开存放其它结果的文件。

while (my $line=<FILE1>) {
  my $decide=0;  #用于判断的开关。
  if ($line=~m/^ID/) {print FILE2 $line; $decide=1; }
  if ($line=~m/^DE/) {print FILE2 $line; $decide=1; }  
  if ($line=~m/^CC\s{3}-!-\sSUBCELLULAR\sLOCATION:/) {     #输出亚细胞定位信息。
    do {
        print FILE2 $line;
        $decide=1; 
        $line=<FILE1>; #读取下一行。
    } until (($line =~ m/^CC\s{3}-!-/) or ($line =~ m/^CC\s{3}--------/));
  }
  if ($line=~m/^SQ/) {   #输出序列。
    do {
      print FILE2 $line;
      $decide=1; 
      $line=<FILE1>; #读取下一行。
     } until ($line=~m"^//");
  }  
  if ($line=~m"^//") {print FILE2 $line; $decide=1; }    #输出一条序列的结束标志。
  if ($decide==0) {print FILE3 $line;}                   #把上述信息以外的信息输出到另外一个文件中。
} 

close FILE1;
close FILE2;
close FILE3;  #必须先关闭文件，再统计序列条数，否则内容可能还没有写入文件，还在缓存中。
&NUMBER_TEXT("$name1",                          "$script_name/$file_n");
&NUMBER_TEXT("$script_name/$file_n/$name2",     "$script_name/$file_n");
&NUMBER_TEXT("$script_name/$file_n/$name3",     "$script_name/$file_n");

&PATTERN_LINE("$name1",                          "$script_name/$file_n",     '(SUBCELLULAR LOCATION:)');
&PATTERN_LINE("$script_name/$file_n/$name2",     "$script_name/$file_n",     '(SUBCELLULAR LOCATION:)');
&PATTERN_LINE("$script_name/$file_n/$name3",     "$script_name/$file_n",     '(SUBCELLULAR LOCATION:)');
}####################################################### END FIVE_2



##########################################################################################################################################
sub NOTE_3   #去掉信息：note=
##########################################################################################################################################
{ 
my $name1=$_[0];                        #待处理文件的名字(含路径)。
my $name2="3.txt";                      #用于存放结果的文件的名称。
my $name3="no_3.txt";                   #用于存放其它结果的文件的名称。
my $file_n = "3_NOTE";
mkdir  "$script_name/$file_n"  or die;  #建立一个子文件夹，存放这个子程序的输出结果。  

open(FILE1, "<", "$name1") or die "$!";                        #以读的方式打开待处理的文件。
open(FILE2, ">", "$script_name/$file_n/$name2") or die "$!";   #以写的方式打开存放结果的文件。
open(FILE3, ">", "$script_name/$file_n/$name3") or die "$!";   #以写的方式打开存放其它结果的文件。

while (my $line=<FILE1>) {  
  if ($line=~m/^ID/) {print FILE2 $line;  }
  if ($line=~m/^DE/) {print FILE2 $line;  }
  my $subcell='';
  if ($line=~m/^CC\s{3}-!-\sSUBCELLULAR\sLOCATION:/) {
    do {
      $subcell=$subcell.$line;  #把亚细胞定位信息放在一个变量里。    
      $line=<FILE1>;
    } until $line=~m/^SQ/ ;
  }
  if ($subcell=~m/(^CC\s{3}-!-\sSUBCELLULAR\sLOCATION:[\d\D]+)(Note=[\d\D]+)(\n)/) {
    print  FILE2 $1,$3;   #去除note信息。
    print  FILE3 "$2 \n\n\n\n";
  }else{
    print FILE2 $subcell;
  }
  if ($line =~ m/^SQ/) {   #输出序列。
    do {
      print FILE2 $line;
      $line=<FILE1>;
    } until ($line=~m"^//");
  }
  if ($line=~m"^//") {print FILE2 $line; }    #输出一条序列的结束标志。
} 

close FILE1;
close FILE2;
close FILE3;  #必须先关闭文件，再统计序列条数，否则内容可能还没有写入文件，还在缓存中。
&NUMBER_TEXT("$name1",                          "$script_name/$file_n");
&NUMBER_TEXT("$script_name/$file_n/$name2",     "$script_name/$file_n");
&NUMBER_TEXT("$script_name/$file_n/$name3",     "$script_name/$file_n");

&PATTERN_LINE("$name1",                          "$script_name/$file_n",     'Note=');
&PATTERN_LINE("$script_name/$file_n/$name2",     "$script_name/$file_n",     'Note=');
&PATTERN_LINE("$script_name/$file_n/$name3",     "$script_name/$file_n",     'Note=');
}############################################ END NOTE_3



################################################################################################################################################
sub FUZZY_4    #此步不进行筛选。 #子程序，亚细胞定位信息中含有(By similarity)、(Potential)、(Probable)的序列一律不要。
################################################################################################################################################
{
my $name1=$_[0];                        #待处理文件的名字(含路径)。
my $name2="4.txt";                      #用于存放结果的文件的名称。
my $name3="no_4.txt";                   #用于存放其它结果的文件的名称。
my $file_n = "4_FUZZY";
mkdir  "$script_name/$file_n"  or die;  #建立一个子文件夹，存放这个子程序的输出结果。  

open(FILE1, "<", "$name1") or die "$!";                        #以读的方式打开待处理的文件。
open(FILE2, ">", "$script_name/$file_n/$name2") or die "$!";   #以写的方式打开存放结果的文件。
open(FILE3, ">", "$script_name/$file_n/$name3") or die "$!";   #以写的方式打开存放其它结果的文件。

while (my $line1=<FILE1>) {
  my $one='';   #用来存储一条序列的所有信息。
  $line1=~m/^ID/  or die;
  do {
    $one=$one.$line1;
    $line1=<FILE1>;
  }until ($line1=~m"^//");
  $one=$one.$line1;
  $one=~m/\nCC\s{3}-!-\sSUBCELLULAR\sLOCATION:([\d\D]+)\nSQ\s/  or die;
  my $subcell1=$1;     #把亚细胞定位信息放在一个变量里。   
  if ($subcell1=~m/Nucleus[\sC]+\([PB]/) {
    print FILE3 $one; 
  }else{ 
    print FILE2 $one; 
  }
}

close FILE1;
close FILE2;
close FILE3;  #必须先关闭文件，再统计序列条数，否则内容可能还没有写入文件，还在缓存中。
&NUMBER_TEXT("$name1",                          "$script_name/$file_n");
&NUMBER_TEXT("$script_name/$file_n/$name2",     "$script_name/$file_n");
&NUMBER_TEXT("$script_name/$file_n/$name3",     "$script_name/$file_n");


&PATTERN_SEQ("$name1",                          "$script_name/$file_n",      'SUBCELLULAR LOCATION:');
&PATTERN_SEQ("$script_name/$file_n/$name2",     "$script_name/$file_n",      'SUBCELLULAR LOCATION:');
&PATTERN_SEQ("$script_name/$file_n/$name3",     "$script_name/$file_n",      'SUBCELLULAR LOCATION:');

&PATTERN_LINE("$name1",                          "$script_name/$file_n",     'SUBCELLULAR LOCATION:');
&PATTERN_LINE("$script_name/$file_n/$name2",     "$script_name/$file_n",     'SUBCELLULAR LOCATION:');
&PATTERN_LINE("$script_name/$file_n/$name3",     "$script_name/$file_n",     'SUBCELLULAR LOCATION:');
}################################# END FUZZY_4 







##########################################################################################################################################
sub FRAGMENT_5  #去掉片段
##########################################################################################################################################
{ 
my $name1=$_[0];                        #待处理文件的名字(含路径)。
my $name2="5.txt";                      #用于存放结果的文件的名称。
my $name3="no_5.txt";                   #用于存放其它结果的文件的名称。
my $file_n = "5_FRAGMENT";
mkdir  "$script_name/$file_n"  or die;  #建立一个子文件夹，存放这个子程序的输出结果。  

open(FILE1, "<", "$name1") or die "$!";                        #以读的方式打开待处理的文件。
open(FILE2, ">", "$script_name/$file_n/$name2") or die "$!";   #以写的方式打开存放结果的文件。
open(FILE3, ">", "$script_name/$file_n/$name3") or die "$!";   #以写的方式打开存放其它结果的文件。

while (my $line=<FILE1>) {
  my $one='';
  $line=~m/^ID/  or die;
  do {
    $one=$one.$line;
    $line=<FILE1>;
  }until ($line=~m"^//");
  $one=$one.$line;
  if ($one=~m/\nDE\s{3}Flags:[^\n]+Fragment[s;]/) {  #此时DE已不位于字符窜的开头。
    print FILE3 $one;
  }else{
    print FILE2 $one;
  }
}   

close FILE1;
close FILE2;
close FILE3; #必须先关闭文件，再统计序列条数，否则内容可能还没有写入文件，还在缓存中。
&NUMBER_TEXT("$name1",  "$script_name/$file_n");
&NUMBER_TEXT("$script_name/$file_n/$name2",     "$script_name/$file_n");
&NUMBER_TEXT("$script_name/$file_n/$name3",     "$script_name/$file_n");

&PATTERN_SEQ("$name1",                          "$script_name/$file_n",      '\nDE\s{3}Flags:[^\n]+Fragment[s;]');
&PATTERN_SEQ("$script_name/$file_n/$name2",     "$script_name/$file_n",      '\nDE\s{3}Flags:[^\n]+Fragment[s;]');
&PATTERN_SEQ("$script_name/$file_n/$name3",     "$script_name/$file_n",      '\nDE\s{3}Flags:[^\n]+Fragment[s;]');

&PATTERN_LINE("$name1",                          "$script_name/$file_n",     'DE\s{3}Flags:[^\n]+Fragment[s;]');
&PATTERN_LINE("$script_name/$file_n/$name2",     "$script_name/$file_n",     'DE\s{3}Flags:[^\n]+Fragment[s;]');
&PATTERN_LINE("$script_name/$file_n/$name3",     "$script_name/$file_n",     'DE\s{3}Flags:[^\n]+Fragment[s;]');
}############################################ END FRAGMENT_5





##########################################################################################################################################
sub UOBZJX_6  #去掉序列中含有UOBZJX的序列。
##########################################################################################################################################
{ 
my $name1=$_[0];                        #待处理文件的名字(含路径)。
my $name2="6.txt";                      #用于存放结果的文件的名称。
my $name3="no_6.txt";                   #用于存放其它结果的文件的名称。
my $file_n = "6_UOBZJX";
mkdir  "$script_name/$file_n"  or die;  #建立一个子文件夹，存放这个子程序的输出结果。  

open(FILE1, "<", "$name1") or die "$!";                        #以读的方式打开待处理的文件。
open(FILE2, ">", "$script_name/$file_n/$name2") or die "$!";   #以写的方式打开存放结果的文件。
open(FILE3, ">", "$script_name/$file_n/$name3") or die "$!";   #以写的方式打开存放其它结果的文件。

while (my $line=<FILE1>) {
  my $one='';
  $line=~m/^ID/  or die;
  do {
    $one=$one.$line;
    $line=<FILE1>;
  }until ($line=~m"^//");
  $one=$one.$line;  
  if ($one=~m/SQ\s+SEQUENCE.+AA;.+MW;.+;[\d\D]+[UOBZJX]+/) {
    print FILE3 $one;
  }else{
    print FILE2 $one;
  }
} 

close FILE1;
close FILE2;
close FILE3; #必须先关闭文件，再统计序列条数，否则内容可能还没有写入文件，还在缓存中。

&NUMBER_TEXT("$name1",  "$script_name/$file_n");
&NUMBER_TEXT("$script_name/$file_n/$name2",     "$script_name/$file_n");
&NUMBER_TEXT("$script_name/$file_n/$name3",     "$script_name/$file_n");

&PATTERN_SEQ("$name1",                          "$script_name/$file_n",      'SQ\s+SEQUENCE.+AA;.+MW;.+;[\d\D]+[UOBZJX]+');
&PATTERN_SEQ("$script_name/$file_n/$name2",     "$script_name/$file_n",      'SQ\s+SEQUENCE.+AA;.+MW;.+;[\d\D]+[UOBZJX]+');
&PATTERN_SEQ("$script_name/$file_n/$name3",     "$script_name/$file_n",      'SQ\s+SEQUENCE.+AA;.+MW;.+;[\d\D]+[UOBZJX]+');

&PATTERN_LINE("$name1",                          "$script_name/$file_n",     '\s{5}[A-Z]+[UOBZJX]');
&PATTERN_LINE("$script_name/$file_n/$name2",     "$script_name/$file_n",     '\s{5}[A-Z]+[UOBZJX]');
&PATTERN_LINE("$script_name/$file_n/$name3",     "$script_name/$file_n",     '\s{5}[A-Z]+[UOBZJX]');
}############################################ END UOBZJX_6





##########################################################################################################################################
sub SHORT_7   #去掉短序列。
##########################################################################################################################################
{ 
my $name1=$_[0];                        #待处理文件的名字(含路径)。
my $name2="7.txt";                      #用于存放结果的文件的名称。
my $name3="no_7.txt";                   #用于存放其它结果的文件的名称。
my $file_n = "7_SHORT";
mkdir  "$script_name/$file_n"  or die;  #建立一个子文件夹，存放这个子程序的输出结果。  

open(FILE1, "<", "$name1") or die "$!";                        #以读的方式打开待处理的文件。
open(FILE2, ">", "$script_name/$file_n/$name2") or die "$!";   #以写的方式打开存放结果的文件。
open(FILE3, ">", "$script_name/$file_n/$name3") or die "$!";   #以写的方式打开存放其它结果的文件。

my $length=50;  #设置序列的最小长度。

while (my $line=<FILE1>) {
  my $one='';
  my $sequence='';
  $line=~m/^ID/  or die;
  do {
    $one=$one.$line;
    $line=<FILE1>;
  }until ($line=~m"^//");
  $one=$one.$line;  
  $one=~m"SQ\s+SEQUENCE.+AA;.+MW;.+;([\d\D]+)//"  or die;
  $sequence=$1;
  my $i1=0;
  for ($i1=1,$i1<=10,$i1++) {  #去掉序列中的空格和换行符
    my @tem=split(/\s/,$sequence);
    $sequence=join('',@tem);
  }  
  if ((length($sequence))>=$length) {print FILE2 $one;}
  if (((length($sequence))<$length) and ((length($sequence))>=1)) {print FILE3 $one;}
}   


close FILE1;
close FILE2;
close FILE3; #必须先关闭文件，再统计序列条数，否则内容可能还没有写入文件，还在缓存中。
&NUMBER_TEXT("$name1",  "$script_name/$file_n");
&NUMBER_TEXT("$script_name/$file_n/$name2",     "$script_name/$file_n");
&NUMBER_TEXT("$script_name/$file_n/$name3",     "$script_name/$file_n");

&PATTERN_SEQ("$name1",  "$script_name/$file_n",      'SUBCELLULAR LOCATION:');
&PATTERN_SEQ("$script_name/$file_n/$name2",     "$script_name/$file_n",      'SUBCELLULAR LOCATION:');
&PATTERN_SEQ("$script_name/$file_n/$name3",     "$script_name/$file_n",      'SUBCELLULAR LOCATION:');

&PATTERN_LINE("$name1",  "$script_name/$file_n",     'SUBCELLULAR LOCATION:');
&PATTERN_LINE("$script_name/$file_n/$name2",     "$script_name/$file_n",     'SUBCELLULAR LOCATION:');
&PATTERN_LINE("$script_name/$file_n/$name3",     "$script_name/$file_n",     'SUBCELLULAR LOCATION:');
}############################################  SHORT_7





##########################################################################################################################################
sub THREE_8    #只输出3类信息：ID,CC,//。
##########################################################################################################################################
{ 
my $name1=$_[0];                        #待处理文件的名字(含路径)。
my $name2="8.txt";                      #用于存放结果的文件的名称。
my $name3="no_8.txt";                   #用于存放其它结果的文件的名称。
my $file_n = "8_THREE";
mkdir  "$script_name/$file_n"  or die;  #建立一个子文件夹，存放这个子程序的输出结果。  

open(FILE1, "<", "$name1") or die "$!";                        #以读的方式打开待处理的文件。
open(FILE2, ">", "$script_name/$file_n/$name2") or die "$!";   #以写的方式打开存放结果的文件。
open(FILE3, ">", "$script_name/$file_n/$name3") or die "$!";   #以写的方式打开存放其它结果的文件。

while (my $line=<FILE1>) {
  my $decide=0;  #用于判断的开关。
  if ($line=~m/^ID/) {print FILE2 $line; $decide=1; }
  if ($line=~m/^CC\s{3}-!-\sSUBCELLULAR\sLOCATION:/) {     #输出亚细胞定位信息。
    do {
        print FILE2 $line;
        $decide=1; 
        $line=<FILE1>; #读取下一行。
    } until $line =~ m/^SQ/;
  }
  if ($line=~m"^//") {print FILE2 $line; $decide=1; }    #输出一条序列的结束标志。
  if ($decide==0) {print FILE3 $line;}  #把上述信息以外的信息输出到另外一个文件中。
  if ($line=~m"^SQ") {print FILE3 $line; }
} 

close FILE1;
close FILE2;
close FILE3; #必须先关闭文件，再统计序列条数，否则内容可能还没有写入文件，还在缓存中。
&NUMBER_TEXT("$name1",                          "$script_name/$file_n");
&NUMBER_TEXT("$script_name/$file_n/$name2",     "$script_name/$file_n");
&NUMBER_TEXT("$script_name/$file_n/$name3",     "$script_name/$file_n");

&PATTERN_LINE("$name1",                          "$script_name/$file_n",     '^SQ');
&PATTERN_LINE("$script_name/$file_n/$name2",     "$script_name/$file_n",     '^SQ');
&PATTERN_LINE("$script_name/$file_n/$name3",     "$script_name/$file_n",     '^SQ');
}############################################ END THREE_8






################################################################################################################################################
sub LOCATION_9   #找出所有单定位和多定位(6类)。
################################################################################################################################################
{
my @locations=(   #每一个元素为一个亚核位置模式。
'',               #0 
'(chromosome)',    #1
'((nucleus[^.]*membrane)|(nucleus[^.]*envelope))',   #2
'(nucleus[^.]*matrix)', #3
'(nuclear[^.]*pore)',  #4
'(speckle)',  #5
'(nucleolus)',  #6
'(Nucleus,[\sc]*nucleoplasm)',  #7
'(PML)',  #8    
'((cajal)|(Nucleus,[\sc]*gem))', #9                                
);

my $file_n = "9_LOCATION";
mkdir  "$script_name/$file_n"  or die;  #建立一个子文件夹，存放这个子程序的输出结果。 
 
my $name1=$_[0];                        #待处理文件的名字(含路径)。
open(FILE1, "<", "$name1") or die "$!";    #以读的方式打开待处理的文件。

while (my $line1=<FILE1>) {
  my $one='';   #用来存储一条序列的所有信息。考虑一个文件中的每一条序列属于哪一类。
  $line1=~m/^ID/  or die;
  do {
    $one=$one.$line1;
    $line1=<FILE1>;
  }until ($line1=~m"^//");
  $one=$one.$line1;
  $one=~m{\nCC\s{3}-!-\sSUBCELLULAR\sLOCATION(:[\d\D]+)//}  or die;
  my $subcell = $1;     #把亚细胞定位信息放在一个变量里。 
  $subcell =~ m/^:[\d\D]+\n$/   or die;

  for (my $i1=1;$i1<=$#locations;$i1++) {
    if ($subcell =~ m/$locations[$i1](?![^.,;]*\([pb])/i) { 
      my $p = 0;
      for (my $j=1;$j<=$#locations;$j++) {
        if(($j != $i1)   and  ($subcell =~ m/$locations[$j]/i)) {$p = 1; last;}
      }
      if($p==0) {
        open(FILE2, ">>", "$script_name/$file_n/$i1.txt") or die;   #以追加的方式打开文件。文件名与位置编号相对应。
        print FILE2  $one;
        goto label;  #一条序列一旦输出，就跳到最后去，考察下一条序列。
      }
    }
  }

  for (my $i1=1;$i1<=$#locations;$i1++) {       #考虑了所有组合，不必担心注释信息中各亚核位置的顺序问题。
    for (my $i2=1;$i2<=$#locations;$i2++) {
      if ($subcell =~ m/$locations[$i1](?![^.,;]*\([pb])([\d\D]*)$locations[$i2](?![^.,;]*\([pb])/i) {   
        my $p = 0;
        for (my $j=1;$j<=$#locations;$j++) {
          if(($j != $i1)   and  ($j != $i2)   and  ($subcell =~ m/$locations[$j]/i)) {$p = 1; last;}
        }
        if($p==0) {
          open(FILE2, ">>", "$script_name/$file_n/$i1=$i2.txt") or die;  
          print FILE2  $one;
          goto label;  #一条序列一旦输出，就跳到最后去，考察下一条序列。
        }
      }
    }
  }


  for (my $i1=1;$i1<=$#locations;$i1++) {   #三定位情况。
    for (my $i2=1;$i2<=$#locations;$i2++) {
      for (my $i3=1;$i3<=$#locations;$i3++) {
        if ($subcell =~ m/$locations[$i1](?![^.,;]*\([pb])([\d\D]*)$locations[$i2](?![^.,;]*\([pb])([\d\D]*)$locations[$i3](?![^.,;]*\([pb])/i) {   
          my $p = 0;
          for (my $j=1;$j<=$#locations;$j++) { #不能匹配其它位置。
            if(($j != $i1)   and  ($j != $i2)   and  ($j != $i3)   and  ($subcell =~ m/$locations[$j]/i)) {$p = 1; last;}
          }
          if($p==0) {
            open(FILE2, ">>", "$script_name/$file_n/$i1=$i2=$i3.txt") or die;
            print FILE2  $one;
            goto label;  #一条序列一旦输出，就跳到最后去，考察下一条序列。
          }  
        }
      }
    }
  }
 

  for (my $i1=1;$i1<=$#locations;$i1++) {   #4定位情况。
    for (my $i2=1;$i2<=$#locations;$i2++) {
      for (my $i3=1;$i3<=$#locations;$i3++) {
        for (my $i4=1;$i4<=$#locations;$i4++) {

          if ($subcell =~ m/$locations[$i1](?![^.,;]*\([pb])([\d\D]*)$locations[$i2](?![^.,;]*\([pb])([\d\D]*)$locations[$i3](?![^.,;]*\([pb])([\d\D]*)$locations[$i4](?![^.,;]*\([pb])/i) {   
            my $p = 0;
            for (my $j=1;$j<=$#locations;$j++) { #不能匹配其它位置。
              if(($j != $i1)   and  ($j != $i2)   and  ($j != $i3)   and  ($j != $i4)  and  ($subcell =~ m/$locations[$j]/i)) {$p = 1; last;}
            }
            if($p==0) {
              open(FILE2, ">>", "$script_name/$file_n/$i1=$i2=$i3=$i4.txt") or die;
              print FILE2  $one;
              goto label;  #一条序列一旦输出，就跳到最后去，考察下一条序列。
            }
          } 
 
        }
      }
    }
  }


  label:{$one='';}  #大括号里面的内容可以为空。
} 
 
my $dir = "$script_name/$file_n";
opendir DIRHANDLE,$dir or die;         #打开目录句柄，读取目录里的文件名。
while (my $file=readdir DIRHANDLE) {   #读取的仅仅是文件名称，不含路径。
  next unless $file=~m/[0-9 -]+\.txt$/;
  &NUMBER_TEXT("$dir/$file", "$script_name/$file_n");          #传递的参数是含有路径的文件名。
}

}###################################################  LOCATION_9




##############################################################################################################################
sub   TEXT_TEXT_10       #以一个TXT文件中序列的ID号为标准，找出另外一个TXT文件中的序列。
###############################################################################################################################
{

my $name1=$_[0];                #待处理的文件（不含路径），首先需要其ID号。
my $name2='subnuclear.txt';     #最终要输出此文件中的序列。
my $name3='a'.$name1;           #输出文件的名字。

open(FILE1, "<", "$script_name/9_LOCATION/$name1") or die "Can't open $name1 : $!";    #以读的方式打开文件。
open(FILE2, "<", "$name2") or die "Can't open $name2 : $!";                   #以读的方式打开文件。
open(FILE3, ">", "$script_name/10_TXT/$name3") or die "Can't open $name3 : $!";   #以写的方式打开文件。

my $IDs='';
while (my $line1=<FILE1>) {
  if ($line1=~m/^ID\s+( [A-Z_0-9]+ )\s/) {
    $IDs=$IDs.$1.'  ';
  }
} 

while (my $line2=<FILE2>) {
  if ($line2=~m/^ID\s{3}([A-Z_0-9]+)\s/) {
    if ($IDs=~m/ $1 /) {
      do {
        print FILE3 $line2;
        $line2=<FILE2>;
      } until $line2=~m"^//";
      print FILE3 $line2;
    }
  }
}

close FILE1;
close FILE2;
close FILE3;

&NUMBER_TEXT("$script_name/10_TXT/$name3",     "$script_name/10_TXT");

}############################################  END TEXT_TEXT_10 





##############################################################################################################################
sub TEXT_FASTA_11    #以一个TXT文件中序列的ID号为标准，找出另外一个FASTA文件中的序列。
###############################################################################################################################
{
my $name1=$_[0];                    #待处理的文件（不含路径），首先需要其ID号。
my $name2='subnuclear.fasta';       #最终要输出此文件中的序列。
my $name3=$name1;   
$name3=~s/\.txt/.fasta/ or die;          #输出文件的名字。
open(FILE1, "<", "$script_name/9_LOCATION/$name1") or die "Can't open $name1 : $!";   #以读的方式打开文件。
open(FILE2, "<", "$name2") or die "Can't open $name2 : $!";   #以读的方式打开文件。
open(FILE3, ">", "$script_name/11_FASTA/$name3") or die "Can't open $name3 : $!";   #以写的方式打开文件。


my $IDs='';
while (my $line1=<FILE1>) {
  if ($line1=~m/^ID\s+( [A-Z_0-9]+ )\s/) {
    $IDs=$IDs.$1.'  ';
  }
} 

my $i1=0;
my @lines=<FILE2>;
while ($lines[$i1]) {
  if ($lines[$i1]=~m/^>[a-z]{2}\|.+\|([A-Z_0-9]+)\s/) {
    my $i2=$i1;
    if ($IDs=~m/ $1 /) {
      do {
        print FILE3 $lines[$i2];
        $i2++;
      } until (($lines[$i2]=~m"^>") or ($lines[$i2] eq ''));
    }
  }
  $i1++;
}

close FILE1;
close FILE2;
close FILE3;
&NUMBER_FASTA("$script_name/11_FASTA/$name3",     "$script_name/11_FASTA");
}######################################################  END TEXT_FASTA_11 







&ISOFORM_1("subnuclear.txt");
&FIVE_2(     "$script_name/1_ISOFORM/1.txt");
&NOTE_3(     "$script_name/2_FIVE/2.txt");
&FUZZY_4(    "$script_name/3_NOTE/3.txt");
&FRAGMENT_5( "$script_name/4_FUZZY/4.txt");
&UOBZJX_6(   "$script_name/5_FRAGMENT/5.txt");
&SHORT_7(    "$script_name/6_UOBZJX/6.txt");
&THREE_8(    "$script_name/7_SHORT/7.txt");
&LOCATION_9( "$script_name/8_THREE/8.txt");


{
mkdir  "$script_name/10_TXT"    or  die ;     #建立一个子文件夹，存放这个子程序的输出结果。 
mkdir  "$script_name/11_FASTA"    or  die ;     #建立一个子文件夹，存放这个子程序的输出结果。 
my $dir = "$script_name/9_LOCATION";
opendir DIRHANDLE,$dir or die;            #打开目录句柄，读取目录里的文件名。
while (my $file=readdir DIRHANDLE) {      #读取的仅仅是文件名称，不含路径。
  next unless $file=~m/[0-9 -]+\.txt$/;
  &TEXT_TEXT_10("$file");                 #传递的参数是不含路径的文件名。
  &TEXT_FASTA_11("$file");
}
}



