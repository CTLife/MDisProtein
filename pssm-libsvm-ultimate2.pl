#最后一次修改于2012年03月26日。

#                                         说明
####################################################################################################################################################
#处理的是单定位蛋白或(和)多定位蛋白的PSSM矩阵。
#处理PSI-BLAST输出的PSSM矩阵,把PSSM矩阵过滤(其特殊情况为不过滤)以后再转化为向量，再转化为LibSVM需要的格式。
#一个文件存放一条序列的PSSM,文件名含序列的ID号。
#转化为向量时，可能不会用到PSSM矩阵的元素，而只需要其序列。
#多定位蛋白的特征向量转化为LibSVM格式时，类别标签取其中一个便可，这样在处理预测结果时，需要自己编程重新统计哪些预测正确。



#要求：
#待处理文件的名字的格式为：ID号.pssm
#不同类别序列的PSSM放在不同文件夹下，文件夹的名字分别为类别标签1、2、3、4、5、......，再把这些文件夹放在一个名为pssm的文件夹(与此程序位于同一目录)下。
#对于多定位蛋白，类别号用短线隔开作为文件名，如1-2，1-2-3,4-9-10-13.



#包括10步:
#1  分割: 输出原始PSSM矩阵。（1个文件对应一条序列的PSSM）（1个子程序）(之后的文件名为序列的ID号)
#2  检查: 检查输出格式是否正确，并去掉含有UOBZJX的行。考虑相邻行，把他们相加。
#3  统计: 统计PSSM矩阵各行的最大值的分布情况。（所有序列的，每一类序列的，每一条序列的。下同）
#        统计PSSM矩阵所有矩阵元的分布。
#        统计20种氨基酸的分布。(此步的结果，后面各步骤不会使用到。此步共9个子程序。) 
#4  缩放: 把PSSM矩阵的元素归一化到[0，1]或不归一，5种情况。(第一次缩放)(处理第2步的结果)
#5  过滤: 过滤掉PSSM矩阵的特定行或特定元素，14种方法。                     (此步涉及到4个命令行参数)
#6  转化: 提取特征向量，目前是5种方法。（400D向量是一个矩阵的格式）（特征信息提取）(此步涉及到6个命令行参数)
#7  一行: 进一步处理向量(把向量变成一行，一条序列对应一行，对应一个文件)。
#8  合并: 把各类的向量分别放在不同的文件里，同一类的放在同一个文件。还要把所有向量放在一个文件里。（4个子程序）
#9  缩放: 把PSSM矩阵的元素归一化到[0，1]或不归一，4种情况。(第二次缩放)
#10 格式: 把向量转化为LibSVM需要的格式。 ####此步由分类器所需要的格式决定。





#12个命令行参数：
#
#1  -scale1  n1   
#    (第一次缩放)（对应第4步）
#    可以为0、1、2、3、4，表示选择哪一种归一化方法： 
#       0:表示不归一化。
#       1:表示1/(1+e(-x))。
#       2:表示(x-min)/(max-min)（按行）。
#       3:表示(x-min)/(max-min)（按列，20种氨基酸分别考虑）。  
#       4:表示(x-min)/(max-min)（按列）。  
#
#2  -filter  n2 
#   （对应第5步）
#    表示过滤方法：
#        1:把一列的某些元素过滤掉，看其是否大于阈值，等价于单独考虑每一个元素的值是否大于阈值。(过滤以后的元素按0输出，转化为向量时不计入。) 此时最好先归一化（按列，20种氨基酸分别考虑）。
#        2:把一列的某些元素过滤掉，看其是否小于阈值，等价于单独考虑每一个元素的值是否小于阈值。(过滤以后的元素按0输出，转化为向量时不计入。) 此时最好先归一化（按列，20种氨基酸分别考虑）。
#        3:把某些行过滤掉，考虑每一行的最大值(max)，在长度为m的窗口内若有n个大于阈值，则留下。(过滤以后的行不再输出）（n可以等于m,下同）
#        4:把某些行过滤掉，考虑每一行的最大值(max)，在长度为m的窗口内若有n个小于阈值，则留下。(过滤以后的行不再输出）                                                                     
#        5:把某些行过滤掉，考虑连续多少行的最大值的平均值，若大于阈值，则留下。             (过滤以后的行不再输出）
#        6:把某些行过滤掉，考虑连续多少行的最大值的平均值，若小于阈值，则留下。             (过滤以后的行不再输出）  
#        7:把某些行过滤掉，考虑20种氨基酸的背景频率，不同类氨基酸应该采用不同阈值，在长度为m的窗口内若有n个大于阈值，则留下。(过滤以后的行不再输出）
#        8:把某些行过滤掉，考虑20种氨基酸的背景频率，不同类氨基酸应该采用不同阈值，在长度为m的窗口内若有n个小于阈值，则留下。(过滤以后的行不再输出） 
#        9:不考虑PSSM矩阵的元素，根据序列过滤，在长度为m的窗口内若有n个氨基酸满足条件，则留下。              (过滤以后的行不再输出） 
#       10:不考虑PSSM矩阵的元素，根据序列过滤，在长度为m的窗口内若有n个氨基酸不满足条件，则留下。             (过滤以后的行不再输出） 
#       11:不考虑PSSM矩阵的元素，根据序列过滤，在长度为m的窗口内若其某种AAindex的平均值大于阈值，则留下。     (过滤以后的行不再输出） 
#       12:不考虑PSSM矩阵的元素，根据序列过滤，在长度为m的窗口内若其某种AAindex的平均值小于阈值，则留下。     (过滤以后的行不再输出）     
#       13:不考虑PSSM矩阵的元素，根据序列过滤，在长度为m的窗口内若有n个氨基酸的某种AAindex大于阈值，则留下。   (过滤以后的行不再输出） 
#       14:不考虑PSSM矩阵的元素，根据序列过滤，在长度为m的窗口内若有n个氨基酸的某种AAindex小于阈值，则留下。   (过滤以后的行不再输出）                 
#
#3  -window1  n3 
#   （对应第5步）
#    为自然数，表示窗口大小。第1、2种过滤方法与这个参数无关。（分别考虑窗口中的所有AA，所有AA的平均值，部分AA。）
#
#4  -window2  n4
#   （对应第5步）
#    考虑在长度为n2的窗口内的n3个氨基酸.(n3<=n2)
#
#5  -threshold  n5
#   （对应第5步）
#    取值范围为[-13,13],是设定的阈值。 第7至第10种过滤方法与这个参数无关。 
#                  
#6  -vector      n6
#   （对应第6步）
#    表示把PSSM矩阵转或序列转化为向量的方法，即特征信息提取方法：
#        1：各行相加，成为20D的向量。
#        2：同种氨基酸的行相加，成为400D的向量。 
#        3：把序列分成N段，每一段转化成一个20维的向量（同1），20*N D。(N=1时与第一种情况相同。)
#        4：序列长度。（1D）
#        5：SAAC。（把20种氨基酸约化成几类或只考虑20种氨基酸中的部分氨基酸，考虑K肽组分）
#
#7  -split       n7
#   （对应第6步）
#    把序列等分成n6段。
#
#8  -whether  n8
#   （对应第6步）
#    为0或1：
#           0：把20种氨基酸约化成几类。      （此时-kinds n9有效）
#           1：只考虑20种氨基酸中的部分氨基酸。（此时-AAs  n10有效）
#
#9  -kinds n9
#   （对应第6步）
#    把20种氨基酸约化成n7种。（约化成几种）
#
#10 -AAs  n10
#   （对应第6步）
#    只考虑20种之中的n8种。  （只考虑几种）
#
#11  -kmer n11
#   （对应第6步）
#    考虑K肽组分.
#
#12  -scale2  n12  (第二次缩放)
#   （对应第9步）
#    可以为0、1、2、3，表示选择哪一种归一化方法： 
#       0:表示不归一化。
#       1:表示1/(1+e(-x))。
#       2:表示(x-min)/(max-min)（按行）。 
#       3:表示(x-min)/(max-min)（按列）。 
###############################################################################################################################################





#!/usr/bin/env perl 
use strict;
use warnings;
use v5.14.2;  #为了使用given-when语句。



($ARGV[0]   eq '-scale1'    )   and  ($ARGV[1]  =~ m/[0-4]/    )     or die "The  parameter is wrong. : $!";
($ARGV[2]   eq '-filter'    )   and  ($ARGV[3]  =~ m/[0-9]+/   )     or die "The  parameter is wrong. : $!";
($ARGV[4]   eq '-window1'   )   and  ($ARGV[5]  =~ m/[0-9]+/   )     or die "The  parameter is wrong. : $!";
($ARGV[6]   eq '-window2'   )   and  ($ARGV[7]  =~ m/[0-9]+/   )     or die "The  parameter is wrong. : $!";
($ARGV[8]   eq '-threshold' )   and  ($ARGV[9]  =~ m/[-0-9.]+/ )     or die "The  parameter is wrong. : $!";
($ARGV[10]  eq '-vector'    )   and  ($ARGV[11] =~ m/[0-9]+/   )     or die "The  parameter is wrong. : $!";
($ARGV[12]  eq '-split'     )   and  ($ARGV[13] =~ m/[0-9]+/   )     or die "The  parameter is wrong. : $!";
($ARGV[14]  eq '-whether'   )   and  ($ARGV[15] =~ m/[01]/     )     or die "The  parameter is wrong. : $!";
($ARGV[16]  eq '-kinds'     )   and  ($ARGV[17] =~ m/[0-9]+/   )     or die "The  parameter is wrong. : $!";
($ARGV[18]  eq '-AAs'       )   and  ($ARGV[19] =~ m/[0-9]+/   )     or die "The  parameter is wrong. : $!";
($ARGV[20]  eq '-kmer'      )   and  ($ARGV[21] =~ m/[0-9]+/   )     or die "The  parameter is wrong. : $!";
($ARGV[22]  eq '-scale2'    )   and  ($ARGV[23] =~ m/[0-3]/    )     or die "The  parameter is wrong. : $!";
($ARGV[24]  eq '-width'    )   and  ($ARGV[25] =~ m/[0-9]+/    )     or die "The  parameter is wrong. : $!";

#故执行此程序的命令的格式是一定的，例如：
#perl   pssm-libsvm-ultimate2.pl     -scale1  0     -filter 1     -window1  5      -window2  5      -threshold  4      -vector 1     -split 1   -whether 0      -kinds 20     -AAs 2    -kmer 2      -scale2 4    -width 0                                                                                                                           


my $whole = "$ARGV[1]"."-$ARGV[3]"."-$ARGV[5]"."-$ARGV[7]"."-$ARGV[9]"."-$ARGV[11]"."-$ARGV[13]"."-$ARGV[15]"."-$ARGV[17]"."-$ARGV[19]"."-$ARGV[21]"."-$ARGV[23]"."-$ARGV[25]";   #是全局变量，每个子程序都会用到。
mkdir "$whole"   or die;   #创建一个文件夹，存放这个程序的所有输出结果。





###################################################
sub MAX  #求任意一组数的最大值的子程序
###################################################
{ 
my $max = shift @_;  #暂时把数组的第一个元素当作最大值。其自变量为任一n维向量(n可以为任一自然数)，函数值为这组数的最大值。
foreach (@_) {       #遍历数组中的其它元素。
  if ($_ > $max) {
    $max = $_;
  }
}
$max;   #返回值
}###################################### END MAX





###################################################
sub MIN   #求任意一组数的最小值的子程序
###################################################
{ 
my $min = shift @_;  #暂时把数组的第一个元素当作最小值。其自变量为任一n维向量(n可以为任一自然数)，函数值为这组数的最小值。
foreach (@_) {       #遍历数组中的其它个元素。
  if ($_ < $min) {
    $min = $_;
  }
}
$min;   #返回值
}###################################### END MIN






###################################################
sub SEQUENCE   #接受一个PSSM文件，返回其序列。
###################################################
{ 
my $name1 = $_[0];                   #接受参数，待处理文件的名字(含路径)。
open(FILE1, "<", "$name1") or die;   #以读的方式打开待处理文件。
my @lines = <FILE1>;
my $seq='';

for(my $i=1; $i<=$#lines; $i++){  #从第2行开始考虑。
  $lines[$i] =~ m/^[ 0-9]{6}([A-Z])\s\s/ or die; 
  my $aa = $1;
  $aa =~ m/^[A-Z]$/ or die;
  $aa !~ m/[UOBZJX]/ or die;
  $seq=$seq.$aa;  
}

length($seq)==$#lines or die;
$seq;   #返回值
}###################################### END SEQUENCE





###################################################
sub LOGXY 
###################################################
{      #以x为底,y的对数.
my $x = $_[0];
my $y = $_[1];
my $m = log($y)/log($x);
return($m);
}############################### END LOGXY








my @mer = ();             #最后，其元素个数应该为:$kinds**$k_m 或 $AAs**$k_m
my $kinds = $ARGV[17];    #约化成多少类。为正整数，从2到20
my $AAss = $ARGV[19];     #只考虑20种氨基酸中的部分氨基酸。为正整数，从2到20
my $k_m = $ARGV[21];      #考虑的是几肽组分。为正整数，从1到10
########################################################
{ #括起来，限制变量作用域。
my @amino_acids=();

if($ARGV[15]==1) {      #只考虑20种氨基酸中的部分氨基酸。
  given($AAss) {
    when(20) {
      my $x1 ="[A]";
      my $x2 ="[C]";
      my $x3 ="[D]";
      my $x4 ="[E]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[I]";
      my $x9 ="[K]";
      my $x10="[L]";
      my $x11="[M]";
      my $x12="[N]";
      my $x13="[P]";
      my $x14="[Q]";
      my $x15="[R]";
      my $x16="[S]";
      my $x17="[T]";
      my $x18="[V]";
      my $x19="[W]";
      my $x20="[Y]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11","$x12","$x13","$x14","$x15","$x16","$x17","$x18","$x19","$x20");                                      
    }
    when(19) {
      my $x1 ="[A]";
      my $x2 ="[C]";
      my $x3 ="[D]";
      my $x4 ="[E]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[I]";
      my $x9 ="[K]";
      my $x10="[L]";
      my $x11="[M]";
      my $x12="[N]";
      my $x13="[P]";
      my $x14="[Q]";
      my $x15="[R]";
      my $x16="[S]";
      my $x17="[T]";
      my $x18="[W]";
      my $x19="[Y]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11","$x12","$x13","$x14","$x15","$x16","$x17","$x18","$x19");                                      
    }
    when(18) {
      my $x1 ="[A]";
      my $x2 ="[C]";
      my $x3 ="[D]";
      my $x4 ="[E]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[I]";
      my $x9 ="[K]";
      my $x10="[M]";
      my $x11="[N]";
      my $x12="[P]";
      my $x13="[Q]";
      my $x14="[R]";
      my $x15="[S]";
      my $x16="[T]";
      my $x17="[W]";
      my $x18="[Y]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11","$x12","$x13","$x14","$x15","$x16","$x17","$x18");                                      
    }
    when(17) {
      my $x1 ="[A]";
      my $x2 ="[C]";
      my $x3 ="[DE]";
      my $x4 ="[Y]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[I]";
      my $x9 ="[K]";
      my $x10="[M]";
      my $x11="[N]";
      my $x12="[P]";
      my $x13="[Q]";
      my $x14="[R]";
      my $x15="[S]";
      my $x16="[T]";
      my $x17="[W]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11","$x12","$x13","$x14","$x15","$x16","$x17");                                      
    }
    when(16) {
      my $x1 ="[A]";
      my $x2 ="[C]";
      my $x3 ="[DE]";
      my $x4 ="[Y]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[I]";
      my $x9 ="[K]";
      my $x10="[M]";
      my $x11="[N]";
      my $x12="[P]";
      my $x13="[Q]";
      my $x14="[R]";
      my $x15="[ST]";
      my $x16="[W]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11","$x12","$x13","$x14","$x15","$x16");                                      
    }
    when(15) {
      my $x1 ="[A]";
      my $x2 ="[C]";
      my $x3 ="[DE]";
      my $x4 ="[Y]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[I]";
      my $x9 ="[K]";
      my $x10="[M]";
      my $x11="[N]";
      my $x12="[P]";
      my $x13="[W]";
      my $x14="[R]";
      my $x15="[S]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11","$x12","$x13","$x14","$x15");                                      
    }
    when(14) {
      my $x1 ="[A]";
      my $x2 ="[C]";
      my $x3 ="[D]";
      my $x4 ="[S]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[I]";
      my $x9 ="[K]";
      my $x10="[M]";
      my $x11="[N]";
      my $x12="[P]";
      my $x13="[W]";
      my $x14="[R]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11","$x12","$x13","$x14");                                      
    }
    when(13) {
      my $x1 ="[R]";
      my $x2 ="[C]";
      my $x3 ="[D]";
      my $x4 ="[S]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[A]";
      my $x9 ="[K]";
      my $x10="[M]";
      my $x11="[N]";
      my $x12="[P]";
      my $x13="[W]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11","$x12","$x13");                                      
    }
    when(12) {
      my $x1 ="[W]";
      my $x2 ="[C]";
      my $x3 ="[D]";
      my $x4 ="[S]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[A]";
      my $x9 ="[K]";
      my $x10="[M]";
      my $x11="[N]";
      my $x12="[P]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11","$x12");                                      
    }
    when(11) {
      my $x1 ="[P]";
      my $x2 ="[C]";
      my $x3 ="[D]";
      my $x4 ="[S]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[A]";
      my $x9 ="[K]";
      my $x10="[M]";
      my $x11="[N]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11");                                      
    }
    when(10) {
      my $x1 ="[P]";
      my $x2 ="[C]";
      my $x3 ="[D]";
      my $x4 ="[S]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[N]";
      my $x8 ="[A]";
      my $x9 ="[H]";
      my $x10="[M]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10");                                      
    }
    when(9) {
      my $x1 ="[P]";
      my $x2 ="[C]";
      my $x3 ="[D]";
      my $x4 ="[S]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[N]";
      my $x8 ="[A]";
      my $x9 ="[H]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9");                                      
    }
    when(8) {
      my $x1 ="[P]";
      my $x2 ="[C]";
      my $x3 ="[D]";
      my $x4 ="[S]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[A]";
      my $x8 ="[H]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8");                                      
    }
    when(7) {
      my $x1="[D]";
      my $x2="[A]";
      my $x3="[S]";
      my $x4="[P]";
      my $x5="[G]";
      my $x6="[C]";
      my $x7="[K]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7");                                      
    }
    when(6) {
      my $x1="[D]";
      my $x2="[A]";
      my $x3="[S]";
      my $x4="[P]";
      my $x5="[G]";
      my $x6="[C]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6");                                      
    }
    when(5) {
      my $x1="[D]";
      my $x2="[A]";
      my $x3="[S]";
      my $x4="[P]";
      my $x5="[C]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5");                                      
    }
    when(4) {
      my $x1="[D]";
      my $x2="[A]";
      my $x3="[S]";
      my $x4="[C]";
      @amino_acids=("$x1","$x2","$x3","$x4");                                      
    }
    when(3) {
      my $x1="[K]";
      my $x2="[R]";
      my $x3="[H]";
      @amino_acids=("$x1","$x2","$x3");                                      
    }
    when(2) {
      my $x1="[K]";
      my $x2="[R]";
      @amino_acids=("$x1","$x2");                                      
    }
  }
}  

if($ARGV[15]==0) {    #下面考虑约化成多少类:
  given($kinds) {
    when(20) {
      my $x1 ="[A]";
      my $x2 ="[C]";
      my $x3 ="[D]";
      my $x4 ="[E]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[I]";
      my $x9 ="[K]";
      my $x10="[L]";
      my $x11="[M]";
      my $x12="[N]";
      my $x13="[P]";
      my $x14="[Q]";
      my $x15="[R]";
      my $x16="[S]";
      my $x17="[T]";
      my $x18="[V]";
      my $x19="[W]";
      my $x20="[Y]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11","$x12","$x13","$x14","$x15","$x16","$x17","$x18","$x19","$x20");                                      
    }
    when(19) {
      my $x1 ="[A]";
      my $x2 ="[C]";
      my $x3 ="[D]";
      my $x4 ="[E]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[I]";
      my $x9 ="[K]";
      my $x10="[LV]";
      my $x11="[M]";
      my $x12="[N]";
      my $x13="[P]";
      my $x14="[Q]";
      my $x15="[R]";
      my $x16="[S]";
      my $x17="[T]";
      my $x18="[W]";
      my $x19="[Y]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11","$x12","$x13","$x14","$x15","$x16","$x17","$x18","$x19");                                      
    }
    when(18) {
      my $x1 ="[A]";
      my $x2 ="[C]";
      my $x3 ="[D]";
      my $x4 ="[E]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[ILV]";
      my $x9 ="[K]";
      my $x10="[M]";
      my $x11="[N]";
      my $x12="[P]";
      my $x13="[Q]";
      my $x14="[R]";
      my $x15="[S]";
      my $x16="[T]";
      my $x17="[W]";
      my $x18="[Y]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11","$x12","$x13","$x14","$x15","$x16","$x17","$x18");                                      
    }
    when(17) {
      my $x1 ="[A]";
      my $x2 ="[C]";
      my $x3 ="[DE]";
      my $x4 ="[Y]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[ILV]";
      my $x9 ="[K]";
      my $x10="[M]";
      my $x11="[N]";
      my $x12="[P]";
      my $x13="[Q]";
      my $x14="[R]";
      my $x15="[S]";
      my $x16="[T]";
      my $x17="[W]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11","$x12","$x13","$x14","$x15","$x16","$x17");                                      
    }
    when(16) {
      my $x1 ="[A]";
      my $x2 ="[C]";
      my $x3 ="[DE]";
      my $x4 ="[Y]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[ILV]";
      my $x9 ="[K]";
      my $x10="[M]";
      my $x11="[N]";
      my $x12="[P]";
      my $x13="[Q]";
      my $x14="[R]";
      my $x15="[ST]";
      my $x16="[W]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11","$x12","$x13","$x14","$x15","$x16");                                      
    }
    when(15) {
      my $x1 ="[A]";
      my $x2 ="[C]";
      my $x3 ="[DE]";
      my $x4 ="[Y]";
      my $x5 ="[F]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[ILV]";
      my $x9 ="[K]";
      my $x10="[M]";
      my $x11="[NQ]";
      my $x12="[P]";
      my $x13="[W]";
      my $x14="[R]";
      my $x15="[ST]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11","$x12","$x13","$x14","$x15");                                      
    }
    when(14) {
      my $x1 ="[A]";
      my $x2 ="[C]";
      my $x3 ="[DE]";
      my $x4 ="[ST]";
      my $x5 ="[FY]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[ILV]";
      my $x9 ="[K]";
      my $x10="[M]";
      my $x11="[NQ]";
      my $x12="[P]";
      my $x13="[W]";
      my $x14="[R]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11","$x12","$x13","$x14");                                      
    }
    when(13) {
      my $x1 ="[R]";
      my $x2 ="[C]";
      my $x3 ="[DE]";
      my $x4 ="[ST]";
      my $x5 ="[FY]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[AILV]";
      my $x9 ="[K]";
      my $x10="[M]";
      my $x11="[NQ]";
      my $x12="[P]";
      my $x13="[W]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11","$x12","$x13");                                      
    }
    when(12) {
      my $x1 ="[W]";
      my $x2 ="[C]";
      my $x3 ="[DE]";
      my $x4 ="[ST]";
      my $x5 ="[FY]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[AILV]";
      my $x9 ="[KR]";
      my $x10="[M]";
      my $x11="[NQ]";
      my $x12="[P]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11","$x12");                                      
    }
    when(11) {
      my $x1 ="[P]";
      my $x2 ="[C]";
      my $x3 ="[DE]";
      my $x4 ="[ST]";
      my $x5 ="[FYW]";
      my $x6 ="[G]";
      my $x7 ="[H]";
      my $x8 ="[AILV]";
      my $x9 ="[KR]";
      my $x10="[M]";
      my $x11="[NQ]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10","$x11");                                      
    }
    when(10) {
      my $x1 ="[P]";
      my $x2 ="[C]";
      my $x3 ="[DE]";
      my $x4 ="[ST]";
      my $x5 ="[FYW]";
      my $x6 ="[G]";
      my $x7 ="[NQ]";
      my $x8 ="[AILV]";
      my $x9 ="[HKR]";
      my $x10="[M]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9","$x10");                                      
    }
    when(9) {
      my $x1 ="[P]";
      my $x2 ="[C]";
      my $x3 ="[DE]";
      my $x4 ="[ST]";
      my $x5 ="[FYW]";
      my $x6 ="[G]";
      my $x7 ="[NQ]";
      my $x8 ="[AILV]";
      my $x9 ="[HKRM]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8","$x9");                                      
    }
    when(8) {
      my $x1 ="[P]";
      my $x2 ="[C]";
      my $x3 ="[DENQ]";
      my $x4 ="[ST]";
      my $x5 ="[FYW]";
      my $x6 ="[G]";
      my $x7 ="[AILV]";
      my $x8 ="[HKRM]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7","$x8");                                      
    }
    when(7) {
      my $x1="[DEHNQ]";
      my $x2="[AFILMV]";
      my $x3="[STWY]";
      my $x4="[P]";
      my $x5="[G]";
      my $x6="[C]";
      my $x7="[KR]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6","$x7");                                      
    }
    when(6) {
      my $x1="[DEHNQKR]";
      my $x2="[AFILMV]";
      my $x3="[STWY]";
      my $x4="[P]";
      my $x5="[G]";
      my $x6="[C]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5","$x6");                                      
    }
    when(5) {
      my $x1="[DEHNQKRG]";
      my $x2="[AFILMV]";
      my $x3="[STWY]";
      my $x4="[P]";
      my $x5="[C]";
      @amino_acids=("$x1","$x2","$x3","$x4","$x5");                                      
    }
    when(4) {
      my $x1="[DEHNQKRGP]";
      my $x2="[AFILMV]";
      my $x3="[STWY]";
      my $x4="[C]";
      @amino_acids=("$x1","$x2","$x3","$x4");                                      
    }
    when(3) {
      my $x1="[NQAFILMVSTWYPGC]";
      my $x2="[KRH]";
      my $x3="[DE]";
      @amino_acids=("$x1","$x2","$x3");                                      
    }
    when(2) {
      my $x1="[GAVLIMFWP]";
      my $x2="[STCYNQDEKRH]";
      @amino_acids=("$x1","$x2");                                      
    }
  }
}

  #下面考虑K肽(k-mer):
  given($k_m){
    when(10){
      for my $aa1(@amino_acids) {
        for my $aa2(@amino_acids) {
          for my $aa3(@amino_acids) {
            for my $aa4(@amino_acids) {
  	      for my $aa5(@amino_acids) {
  	        for my $aa6(@amino_acids) {
  	          for my $aa7(@amino_acids) {
  	            for my $aa8(@amino_acids) {
  	              for my $aa9(@amino_acids) {
  	                for my $aa10(@amino_acids) {
                          $mer[++$#mer]=$aa1.$aa2.$aa3.$aa4.$aa5.$aa6.$aa7.$aa8.$aa9.$aa10;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    when(9){
      for my $aa1(@amino_acids) {
        for my $aa2(@amino_acids) {
          for my $aa3(@amino_acids) {
            for my $aa4(@amino_acids) {
  	      for my $aa5(@amino_acids) {
  	        for my $aa6(@amino_acids) {
  	          for my $aa7(@amino_acids) {
  	            for my $aa8(@amino_acids) {
  	              for my $aa9(@amino_acids) {
                        $mer[++$#mer]=$aa1.$aa2.$aa3.$aa4.$aa5.$aa6.$aa7.$aa8.$aa9;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    when(8){
      for my $aa1(@amino_acids) {
        for my $aa2(@amino_acids) {
          for my $aa3(@amino_acids) {
            for my $aa4(@amino_acids) {
  	      for my $aa5(@amino_acids) {
  	        for my $aa6(@amino_acids) {
  	          for my $aa7(@amino_acids) {
  	            for my $aa8(@amino_acids) {
                      $mer[++$#mer]=$aa1.$aa2.$aa3.$aa4.$aa5.$aa6.$aa7.$aa8;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    when(7){
      for my $aa1(@amino_acids) {
        for my $aa2(@amino_acids) {
          for my $aa3(@amino_acids) {
            for my $aa4(@amino_acids) {
  	      for my $aa5(@amino_acids) {
  	        for my $aa6(@amino_acids) {
  	          for my $aa7(@amino_acids) {
                    $mer[++$#mer]=$aa1.$aa2.$aa3.$aa4.$aa5.$aa6.$aa7;
                  }
                }
              }
            }
          }
        }
      }
    }
    when(6){
      for my $aa1(@amino_acids) {
        for my $aa2(@amino_acids) {
          for my $aa3(@amino_acids) {
            for my $aa4(@amino_acids) {
  	      for my $aa5(@amino_acids) {
  	        for my $aa6(@amino_acids) {
                  $mer[++$#mer]=$aa1.$aa2.$aa3.$aa4.$aa5.$aa6;
                }
              }
            }
          }
        }
      }
    }
    when(5){
      for my $aa1(@amino_acids) {
        for my $aa2(@amino_acids) {
          for my $aa3(@amino_acids) {
            for my $aa4(@amino_acids) {
  	      for my $aa5(@amino_acids) {
                $mer[++$#mer]=$aa1.$aa2.$aa3.$aa4.$aa5;
              }
            }
          }
        }
      }
    }
    when(4){
      for my $aa1(@amino_acids) {
        for my $aa2(@amino_acids) {
          for my $aa3(@amino_acids) {
            for my $aa4(@amino_acids) {
              $mer[++$#mer]=$aa1.$aa2.$aa3.$aa4;
            }
          }
        }
      }
    }
    when(3){
      for my $aa1(@amino_acids) {
        for my $aa2(@amino_acids) {
          for my $aa3(@amino_acids) {
            $mer[++$#mer]=$aa1.$aa2.$aa3;
          }
        }
      }
    }
    when(2){
      for my $aa1(@amino_acids) {
        for my $aa2(@amino_acids) {
          $mer[++$#mer]=$aa1.$aa2;
        }
      }
    }
    when(1){
      for my $aa1(@amino_acids) {
        $mer[++$#mer]=$aa1;
      }
    }
  }

my $num_mer = @mer;
($num_mer == $kinds**$k_m) or ($num_mer == $AAss**$k_m) or die;  #元素个数必须正确。约化成几种或只考虑几种。
open(FILE2, ">", "$whole/info_of_kmer") or die;    #以写的方式打开存放结果的文件。

for (@mer) {
  print FILE2 ;
  print FILE2 "\n";
}
print FILE2 "\n  kmer数组的元素个数：$num_mer \n"; 

} ########################################################




 


##############################################################################################
sub SEPARATE_1    #把初始文件一分为二,只输出PSSM矩阵。
##############################################################################################
{ 
my $name1 = $_[0];                                                 #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^pssm) \/ ([0-9-]+) \/ ([^\s]+)\.pssm$ /x    or die; #可随意加入空白，空白会被忽略。
my $folder = $2;                                                   #把类别标签取出来。
my $name2 = $3;                                                    #把ID号取出来。

if (!(-e "$whole/1_SEPARATE"))         {mkdir "$whole/1_SEPARATE" or die;}           #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/1_SEPARATE/$folder")) {mkdir "$whole/1_SEPARATE/$folder" or die;}   #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。
                                           
open(FILE1, "<", "$name1") or die;                              #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/1_SEPARATE/$folder/$name2") or die;    #以写的方式打开存放结果的文件。

while (my $line=<FILE1>) {   #一行一行地依次读取。
  $line=<FILE1>;
  $line=<FILE1>;    #前2行不要。下一行的格式是一定的。
  $line=~m/(^ {9})(  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V)(   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V)/ or die;    #依次为：9个空格，每个字母前2个空格（20个字母），每个字母前3个空格（20个字母）。                                                                                                                            
  my $g = ' 'x19;   #19个空格。

  $line = $g.' '.$g.'A'.$g.'R'.$g.'N'.$g.'D'.$g.'C'.$g.'Q'.$g.'E'.$g.'G'.$g.'H'.$g.'I'.$g.'L'.$g.'K'.$g.'M'.$g.'F'.$g.'P'.$g.'S'.$g.'T'.$g.'W'.$g.'Y'.$g.'V';                  #依次为：先是20个空格，然后每个字母前19个空格（20个字母）。
  $line=~m/ (^\s{20}) (\s{19}A  \s{19}R  \s{19}N  \s{19}D  \s{19}C  \s{19}Q  \s{19}E  \s{19}G  \s{19}H  \s{19}I  \s{19}L  \s{19}K  \s{19}M  \s{19}F  \s{19}P  \s{19}S  \s{19}T  \s{19}W  \s{19}Y  \s{19}V)/x  or  die; 
  $line =~ m/ ^ ([A-Z\s]{20}){21} $ /x   or die;   #连续这3行是同一个意思。

  print FILE2 $line,"\n";  #输出第一行。
  
  $line=<FILE1>;           #取下一行。
  do{
    $line=~m/(^[ 0-9A-Z*]{9})([ 0-9-]{60})([-. 0-9]{81})/ or die "\n $name1 \n";    #接下来的每一行都应当匹配这个模式。
    print FILE2 $1,' ' x 11;                                        #一行输出的开始，先输出20个字符。
    my $mid1 = $2;                                                  #取出60个字符：20个字母，每个字母前2个空格。  
    for (my $index=0;$index<=59;$index=$index+3) {
      my $mid3 = substr($mid1,$index,3);  #取出子串。
      printf FILE2 "%20s",$mid3;	
    }
    print FILE2 "\n"; #一行输出的结束
    $line=<FILE1>;  
  }until ($line eq "\n");
  
  do {
    $line=<FILE1>;
  }until $line=~m/^PSI Gapped/;

}                              

close FILE1;
close FILE2;
}######################################## END SEPARATE_1





##############################################################################################
sub CHECK_2     #检查输出格式是否正确，并去掉含有UOBZJX的行。
##############################################################################################
{ 
my $name1 = $_[0];                                                         #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/1_SEPARATE) \/ ([0-9-]+) \/ ([^\s]+)$ /x    or die; #可随意加入空白，空白会被忽略。
my $folder = $2;                                                           #把类别标签取出来。
my $name2 = $3;                                                            #把ID号取出来。

if (!(-e "$whole/2_CHECK"))         {mkdir "$whole/2_CHECK" or die;}           #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/2_CHECK/$folder")) {mkdir "$whole/2_CHECK/$folder" or die;}   #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。
                                           
open(FILE1, "<", "$name1") or die;                           #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/2_CHECK/$folder/$name2") or die;    #以写的方式打开存放结果的文件。
open(FILE3, ">>", "$whole/2_CHECK/num_uobzjx") or die;       #以追加的方式打开存放结果的文件

my @lines = <FILE1>;
my @array2d = ();
my @temp2d = ();

for(my $i=0; $i<=$#lines; $i++) {
  $lines[$i] =~ m/^[^\n]{420}\n$/ or die;  #每一行都是421个字符，20+400+1，包括最后的换行符。
  if ($lines[$i] =~ m/^\s{20}/) {
    $i==0  or die;
    $lines[$i] =~ m/ ^\s{20}   ([\sA-Z]{20}){20}   \n$ /x   or die;
    $lines[$i] =~ m/ (^\s{20}) (\s{19}A  \s{19}R  \s{19}N  \s{19}D  \s{19}C  \s{19}Q  \s{19}E  \s{19}G  \s{19}H  \s{19}I  \s{19}L  \s{19}K  \s{19}M  \s{19}F  \s{19}P  \s{19}S  \s{19}T  \s{19}W  \s{19}Y  \s{19}V) \n$ /x    or  die; 
    print FILE2  $lines[$i] ;
  }else{
    $i>0  or die;
    $lines[$i] =~ m/ ^[0-9\sA-Z*]{20}   ([-\s0-9]{20}){20}   \n$ /x   or die;
    my $temp = $lines[$i];
    $temp =~ s/^\s+([0-9])/$1/  or  die;
    $temp =~ s/([0-9])\n$/$1/  or  die;
    $array2d[$i] = [split(/\s+/,$temp)];
    $#{$array2d[$i]} == 21 or die "$name1,$i";
    $temp2d[$i] = [split(/\s+/,$temp)];
    $#{$temp2d[$i]} == 21 or die "$name1,$i";
  }
}                              

my $up_down = $ARGV[25];
for(my $i=1; $i<=$#lines; $i++) {
  my $start = $i - $up_down;
  if($start < 1) {$start=1;}
  my $end = $i + $up_down;
  if($end > $#lines) {$end=$#lines;} 

  for(my $j=2; $j<=21; $j++) {
    for(my $i1=$start; $i1<=$end; $i1++) {
      $temp2d[$i][$j] =~ m/^[-0-9]+$/  or die "$name1, $i \n";
      $array2d[$i1][$j] =~ m/^[-0-9]+$/  or die "\n $name1, $i,$j, \n $array2d[$i1][$j] \n";
      if($i1 != $i) {$temp2d[$i][$j] = $temp2d[$i][$j] + $array2d[$i1][$j]; }
    }
  }

  if($lines[$i] =~ m/[UOBZJX*]/) {
    print FILE3  "$name2: $lines[$i]";
  }else{
    my $pre20 = substr($lines[$i],0,20);
    print FILE2 "$pre20";
    for(my $j1=2; $j1<=21; $j1++) {
      my $temp1 =$temp2d[$i][$j1];  
      printf FILE2 "%20s",$temp1;
    }
    print FILE2 "\n";
  }

}




close FILE1;
close FILE2;
}############################################# END CHECK_2









##############################################################################################
sub STATISTIC_3_1     #统计每一条序列的PSSM矩阵的每一行的最大值的分布。
##############################################################################################
{ 
my $name1 = $_[0];                                                    #接受参数，待处理文件的名字(含路径)，一个文件一条序列。
$name1 =~ m/ (^$whole\/2_CHECK) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  #可随意加入空白，空白会被忽略。
my $folder = $2;                                                      #把类别标签取出来。
my $name2 = $3;                                                       #把不带路径的文件名取出来。

if (!(-e "$whole/3_STATISTIC_1"))          {mkdir "$whole/3_STATISTIC_1" or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/3_STATISTIC_1/$folder"))  {mkdir "$whole/3_STATISTIC_1/$folder" or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。   
                                        
open(FILE1, "<", "$name1") or die;                                #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/3_STATISTIC_1/$folder/$name2") or die;   #以写的方式打开存放结果的文件。

my @max = ();                              #存储1条序列的PSSM矩阵的每一行的最大值的分布。
for(my $i=0;$i<=100;$i++) {$max[$i]=0;}    #初始化数组。

my $seq_len = 0;                           #序列长度。
my $panduan = 0;                           #起判断作用。 
  
while (my $line=<FILE1>) {
  $line =~ m/^[^\n]{420}\n$/ or die;                                     #每一行都匹配这种模式。
  if ($line =~ m/^\s{20}/) {                                             #只有第1行匹配这种模式。故下面的语句只应被执行一次，否则停止执行。
    $panduan==0 or die;
    $panduan=1;
    $line =~ m/ ^\s{20}   ([\sA-Z]{20}){20}   \n$ /x   or die;           #只有第1行匹配这种模式。
  }else{
    $seq_len++;
    $line =~ m/ ^[0-9\sA-Z]{20}   ([-\s0-9]{20}){20}   \n$ /x   or die;  #除了第1行以外，剩下的所有行都应匹配这种模式。
    $line=~m/(^[ 0-9A-Z]{20})([ 0-9-]{400})\n$/ or die;                  #除了第1行以外，剩下的所有行都应匹配这种模式。
    my $mid1=$2;
    $mid1 =~ s/ ^\s* ([-0-9][^\n]+[0-9]) \s*$ /$1/x or die;              #去掉首尾的空格。
    my @array1=split(/\s+/,$mid1);                                       #把字符串变成数组。
    $#array1==19 or die;                                                 #数组的元素个数应该为20，即最大下标为19.
    my $temp = &MAX(@array1);                                            #求出数组的最大值。
    given ($temp) {                                                      #统计各行的最大值，各出现了多少次。
      when ($temp <= -20)  {$max[0]++;}
      when ($temp == -19)  {$max[1]++;}
      when ($temp == -18)  {$max[2]++;}
      when ($temp == -17)  {$max[3]++;}
      when ($temp == -16)  {$max[4]++;}
      when ($temp == -15)  {$max[5]++;}
      when ($temp == -14)  {$max[6]++;}
      when ($temp == -13)  {$max[7]++;}
      when ($temp == -12)  {$max[8]++;}
      when ($temp == -11)  {$max[9]++;}
      when ($temp == -10)  {$max[10]++;}
      when ($temp ==  -9)  {$max[11]++;}
      when ($temp ==  -8)  {$max[12]++;}    
      when ($temp ==  -7)  {$max[13]++;}
      when ($temp ==  -6)  {$max[14]++;}
      when ($temp ==  -5)  {$max[15]++;}
      when ($temp ==  -4)  {$max[16]++;}
      when ($temp ==  -3)  {$max[17]++;}
      when ($temp ==  -2)  {$max[18]++;}
      when ($temp ==  -1)  {$max[19]++;}
      when ($temp ==   0)  {$max[20]++;}
      when ($temp ==   1)  {$max[21]++;}
      when ($temp ==   2)  {$max[22]++;}
      when ($temp ==   3)  {$max[23]++;}
      when ($temp ==   4)  {$max[24]++;}
      when ($temp ==   5)  {$max[25]++;}
      when ($temp ==   6)  {$max[26]++;}
      when ($temp ==   7)  {$max[27]++;}
      when ($temp ==   8)  {$max[28]++;}
      when ($temp ==   9)  {$max[29]++;}
      when ($temp ==  10)  {$max[30]++;}
      when ($temp ==  11)  {$max[31]++;}
      when ($temp ==  12)  {$max[32]++;}    
      when ($temp ==  13)  {$max[33]++;}
      when ($temp ==  14)  {$max[34]++;}
      when ($temp ==  15)  {$max[35]++;}
      when ($temp ==  16)  {$max[36]++;}
      when ($temp ==  17)  {$max[37]++;}
      when ($temp ==  18)  {$max[38]++;}
      when ($temp ==  19)  {$max[39]++;}
      when ($temp >=  20)  {$max[40]++;}
    }                             
  }
}

my $sum = 0;
for (my $i=0; $i<=40; $i++) {
  $sum = $sum + $max[$i];
}
($sum == $seq_len)  or die;    #数组@max的所有元素之和应为该序列的长度。


my @ratio = ();                   #频数数组@max所对应的频率数组。          
my $percent = 0;                  #频率之和，应该为1。
for (my $i=0; $i<=40; $i++) {     #由频数求频率。  
  $ratio[$i] = $max[$i]/$sum;
  $percent = $percent + $ratio[$i];
}


my @a = ();                              #累积频数。
for(my $i=0;$i<=100;$i++) {$a[$i]=0;}    #初始化数组。

for(my $i=0;$i<=$#max;$i++) {            #求累积频数分布。
  for(my $j=0;$j<=$i;$j++){
    $a[$i]=$a[$i]+$max[$j];
  }
}


my @b = ();                              #累积频率。
for(my $i=0;$i<=100;$i++) {$b[$i]=0;}    #初始化数组。

for(my $i=0;$i<=$#ratio;$i++) {          #求累积频率。
  for(my $j=0;$j<=$i;$j++){
    $b[$i]=$b[$i]+$ratio[$j];
  }
}

#以下是输出：
for(my $i1 = -20; $i1<=20;$i1++) {
  my $i2 = $i1 + 20;
  printf  FILE2  "max == %3s :     %-10s     %-10s     %-10.7f     %-10.7f  \n", $i1, $max[$i2], $a[$i2], $ratio[$i2], $b[$i2];
}

print FILE2 "
频数的和为$sum ,  频率的和为$percent
序列长度(以氨基酸为单位): $seq_len  AA.
\n\n\n\n\n\n\n";    #输出中的8个变量依次表示：最大值，频数，累积频数,频率，累积频率. 频数的和，频率的和，序列长度(多少个氨基酸)。

close FILE1;
close FILE2;
} ############################################# END STATISTIC_3_1






##############################################################################################
sub STATISTIC_3_2     #统计每一类序列的PSSM矩阵的每一行的最大值的分布情况。
##############################################################################################
{ 
my $dir1 = "$whole/2_CHECK";
opendir(my $dh1, $dir1)  ||  die "can't opendir $dir1: $!"; 

while (my $folder = readdir $dh1) {  #读取此目录下的所有目录。一个文件夹一个文件夹地处理，即一类一类地处理。
  my $num_seq = 0;                   #每一类所含的序列条数。
  my $aa = 0;                        #每一类的氨基酸总数。
  my @max = ();
  for(my $i=0;$i<=100;$i++) {$max[$i]=0;}    #初始化数组。

  next unless $folder =~ m/^[0-9-]+$/;                               #目录名只含类别标签。
  open(FILE1, ">", "$whole/3_STATISTIC_1/$folder-max.dis1") or die;  #以写的方式打开存放结果的文件。
  my $dir2 = "$dir1/$folder";               
  opendir(my $dh2, $dir2)  ||  die "can't opendir $dir2: $!";        #打开目录句柄，读取目录里的文件名。

  while (my $file=readdir $dh2) {               #读取的仅仅是文件名称，不含路径。一个目录里是一类序列。
    next unless $file !~ m/^[.]/; 
    next unless $file !~ m/[~]$/;
    $num_seq++; #每一类所含的序列条数。
    open(FILE, "<", "$dir2/$file") or die ;     #以读的方式打开待处理文件。
    while (my $line=<FILE>) {
      $line =~ m/^[^\n]{420}\n$/ or die;
      if ($line =~ m/^\s{20}/) {
        $line =~ m/ ^\s{20}   ([\sA-Z]{20}){20}   \n$ /x   or die;
      }else{
        $aa++;  #每一类的氨基酸总数。
        $line =~ m/ ^[0-9\sA-Z]{20}   ([-\s0-9]{20}){20}   \n$ /x   or die;
        $line=~m/(^[ 0-9A-Z]{20})([ 0-9-]{400})\n$/ or die; 
        my $mid1=$2;
        $mid1 =~ s/ ^\s* ([-0-9][^\n]+[0-9]) \s*$ /$1/x or die;   #去掉首尾的空格。
        my @array1=split(/\s+/,$mid1);
        $#array1==19 or die;   #数组的元素个数应该为20，即最大下标为19.
        my $temp = &MAX(@array1);
        given ($temp) {
          when ($temp <= -20)  {$max[0]++;}
          when ($temp == -19)  {$max[1]++;}
          when ($temp == -18)  {$max[2]++;}
          when ($temp == -17)  {$max[3]++;}
          when ($temp == -16)  {$max[4]++;}
          when ($temp == -15)  {$max[5]++;}
          when ($temp == -14)  {$max[6]++;}
          when ($temp == -13)  {$max[7]++;}
          when ($temp == -12)  {$max[8]++;}
          when ($temp == -11)  {$max[9]++;}
          when ($temp == -10)  {$max[10]++;}
          when ($temp ==  -9)  {$max[11]++;}
          when ($temp ==  -8)  {$max[12]++;}    
          when ($temp ==  -7)  {$max[13]++;}
          when ($temp ==  -6)  {$max[14]++;}
          when ($temp ==  -5)  {$max[15]++;}
          when ($temp ==  -4)  {$max[16]++;}
          when ($temp ==  -3)  {$max[17]++;}
          when ($temp ==  -2)  {$max[18]++;}
          when ($temp ==  -1)  {$max[19]++;}
          when ($temp ==   0)  {$max[20]++;}
          when ($temp ==   1)  {$max[21]++;}
          when ($temp ==   2)  {$max[22]++;}
          when ($temp ==   3)  {$max[23]++;}
          when ($temp ==   4)  {$max[24]++;}
          when ($temp ==   5)  {$max[25]++;}
          when ($temp ==   6)  {$max[26]++;}
          when ($temp ==   7)  {$max[27]++;}
          when ($temp ==   8)  {$max[28]++;}
          when ($temp ==   9)  {$max[29]++;}
          when ($temp ==  10)  {$max[30]++;}
          when ($temp ==  11)  {$max[31]++;}
          when ($temp ==  12)  {$max[32]++;}    
          when ($temp ==  13)  {$max[33]++;}
          when ($temp ==  14)  {$max[34]++;}
          when ($temp ==  15)  {$max[35]++;}
          when ($temp ==  16)  {$max[36]++;}
          when ($temp ==  17)  {$max[37]++;}
          when ($temp ==  18)  {$max[38]++;}
          when ($temp ==  19)  {$max[39]++;}
          when ($temp >=  20)  {$max[40]++;}
        }
      }
    }                             
  }#一个文件夹处理完毕。即一类序列处理完毕。

  my $sum = 0; #频数的和
  for (my $i=0; $i<=40; $i++) {
    $sum = $sum + $max[$i];
  }
  ($sum == $aa)  or die;  #频数的和应该等于这一类序列的氨基酸总数。


  my @ratio = ();  #频数数组@max所对应的频率数组。
  my $percent = 0; #频率之和，应该为1。
  for (my $i=0; $i<=40; $i++) {
    $ratio[$i] = $max[$i]/$sum; #由频数求频率。
    $percent = $percent + $ratio[$i];
  }

  my @a = ();                              #累积频数。
  for(my $i=0;$i<=100;$i++) {$a[$i]=0;}    #初始化数组。

  for(my $i=0;$i<=$#max;$i++) {            #求累积频数分布。
    for(my $j=0;$j<=$i;$j++){
      $a[$i]=$a[$i]+$max[$j];
    }
  }

  my @b = ();
  for(my $i=0;$i<=100;$i++) {$b[$i]=0;}    #初始化数组。

  for(my $i=0;$i<=$#ratio;$i++) {   #求累积频率。
    for(my $j=0;$j<=$i;$j++){
      $b[$i]=$b[$i]+$ratio[$j];
    }
  }

  #以下是输出：
  for(my $i1 = -20; $i1<=20;$i1++) {
    my $i2 = $i1 + 20;
    printf  FILE1  "max == %3s :     %-10s     %-10s     %-10.7f     %-10.7f  \n", $i1, $max[$i2], $a[$i2], $ratio[$i2], $b[$i2];
  }

  print FILE1 "
  频数的和为$sum , 频率的和为$percent
  这一类的氨基酸数: $aa
  这一类的序列数:$num_seq
  \n\n\n\n\n\n\n";      #输出中的8个变量依次表示：最大值，频数，累积频数，频率，累积频率。频数的和，频率的和，每一类的氨基酸数，每一类的序列数。
}

close FILE;
close FILE1;
}############################################# END STATISTIC_3_2







##############################################################################################
sub STATISTIC_3_3     #统计所有序列的PSSM矩阵的每一行的最大值的分布情况。
##############################################################################################
{ 
my @max = ();
for(my $i=0;$i<=100;$i++) {$max[$i]=0;}    #初始化数组。
open(FILE1, ">", "$whole/3_STATISTIC_1/all-max.dis1") or die;       #以写的方式打开存放结果的文件。

my $dir1 = "$whole/2_CHECK";
opendir(my $dh1, $dir1)  ||  die "can't opendir $dir1: $!"; 
my $num_seq = 0;     #序列条数。
my $aa = 0;          #氨基酸总数。

while (my $folder = readdir $dh1) {                             #读取此目录下的所有目录。
  next unless $folder =~ m/^[0-9-]+$/;                           #目录名只含类别标签。
  my $dir2 = "$dir1/$folder";               
  opendir(my $dh2, $dir2)  ||  die "can't opendir $dir2: $!";   #打开目录句柄，读取目录里的文件名。

  while (my $file=readdir $dh2) {                               #读取的仅仅是文件名称，不含路径。
    next unless $file !~ m/^[.]/; 
    next unless $file !~ m/[~]$/;
    $num_seq++; #序列条数。
    open(FILE, "<", "$dir2/$file") or die ;                     #以读的方式打开待处理文件。
    while (my $line=<FILE>) {
      $line =~ m/^[^\n]{420}\n$/ or die;
      if ($line =~ m/^\s{20}/) {
        $line =~ m/ ^\s{20}   ([\sA-Z]{20}){20}   \n$ /x   or die;
      }else{
        $aa++; #氨基酸总数。
        $line =~ m/ ^[0-9\sA-Z]{20}   ([-\s0-9]{20}){20}   \n$ /x   or die;
        $line=~m/(^[ 0-9A-Z]{20})([ 0-9-]{400})\n$/ or die; 
        my $mid1=$2;
        $mid1 =~ s/ ^\s* ([-0-9][^\n]+[0-9]) \s*$ /$1/x or die;   #去掉首尾的空格。
        my @array1=split(/\s+/,$mid1);
        $#array1==19 or die;                                      #数组的元素个数应该为20，即最大下标为19.
        my $temp = &MAX(@array1);
        given ($temp) {
          when ($temp <= -20)  {$max[0]++;}
          when ($temp == -19)  {$max[1]++;}
          when ($temp == -18)  {$max[2]++;}
          when ($temp == -17)  {$max[3]++;}
          when ($temp == -16)  {$max[4]++;}
          when ($temp == -15)  {$max[5]++;}
          when ($temp == -14)  {$max[6]++;}
          when ($temp == -13)  {$max[7]++;}
          when ($temp == -12)  {$max[8]++;}
          when ($temp == -11)  {$max[9]++;}
          when ($temp == -10)  {$max[10]++;}
          when ($temp ==  -9)  {$max[11]++;}
          when ($temp ==  -8)  {$max[12]++;}    
          when ($temp ==  -7)  {$max[13]++;}
          when ($temp ==  -6)  {$max[14]++;}
          when ($temp ==  -5)  {$max[15]++;}
          when ($temp ==  -4)  {$max[16]++;}
          when ($temp ==  -3)  {$max[17]++;}
          when ($temp ==  -2)  {$max[18]++;}
          when ($temp ==  -1)  {$max[19]++;}
          when ($temp ==   0)  {$max[20]++;}
          when ($temp ==   1)  {$max[21]++;}
          when ($temp ==   2)  {$max[22]++;}
          when ($temp ==   3)  {$max[23]++;}
          when ($temp ==   4)  {$max[24]++;}
          when ($temp ==   5)  {$max[25]++;}
          when ($temp ==   6)  {$max[26]++;}
          when ($temp ==   7)  {$max[27]++;}
          when ($temp ==   8)  {$max[28]++;}
          when ($temp ==   9)  {$max[29]++;}
          when ($temp ==  10)  {$max[30]++;}
          when ($temp ==  11)  {$max[31]++;}
          when ($temp ==  12)  {$max[32]++;}    
          when ($temp ==  13)  {$max[33]++;}
          when ($temp ==  14)  {$max[34]++;}
          when ($temp ==  15)  {$max[35]++;}
          when ($temp ==  16)  {$max[36]++;}
          when ($temp ==  17)  {$max[37]++;}
          when ($temp ==  18)  {$max[38]++;}
          when ($temp ==  19)  {$max[39]++;}
          when ($temp >=  20)  {$max[40]++;}
        }
      }
    }                             
  }
}

my $sum = 0;  #频数的和
for (my $i=0; $i<=40; $i++) {
  $sum = $sum + $max[$i];
}
($sum == $aa)  or die;  #频数的和应该等于氨基酸总数。


my @ratio = ();   #频数数组@max所对应的频率数组。
my $percent = 0;  #频率的和。
for (my $i=0; $i<=40; $i++) {
  $ratio[$i] = $max[$i]/$sum;
  $percent = $percent + $ratio[$i];
}

my @a = ();                              #累积频数。
for(my $i=0;$i<=100;$i++) {$a[$i]=0;}    #初始化数组。

for(my $i=0;$i<=$#max;$i++) {            #求累积频数分布。
  for(my $j=0;$j<=$i;$j++){
    $a[$i]=$a[$i]+$max[$j];
  }
}

my @b = ();
for(my $i=0;$i<=100;$i++) {$b[$i]=0;}    #初始化数组。

for(my $i=0;$i<=$#ratio;$i++) {   #求累积频率分布。
  for(my $j=0;$j<=$i;$j++){
    $b[$i]=$b[$i]+$ratio[$j];
  }
}

#以下是输出：
for(my $i1 = -20; $i1<=20;$i1++) {
  my $i2 = $i1 + 20;
  printf  FILE1  "max == %3s :     %-10s     %-10s     %-10.7f     %-10.7f  \n", $i1, $max[$i2], $a[$i2], $ratio[$i2], $b[$i2];
}

print FILE1 "
频数的和:$sum , 频率的和：$percent
所有序列的氨基酸数: $aa
序列数:$num_seq
\n\n\n\n\n\n\n";                 #输出中的8个变量依次表示：最大值，频数，累积频数，频率，累积频率。频数的和，频率的和，所有序列的氨基酸数，序列数。

close FILE;
close FILE1;
} ############################################# END STATISTIC_3_3






##############################################################################################
sub STATISTIC_3_4    #统计每一条序列的PSSM矩阵的每一个元素的分布。
##############################################################################################
{ 
my $name1 = $_[0];                                                    #接受参数，待处理文件的名字(含路径)，一个文件一条序列。
$name1 =~ m/ (^$whole\/2_CHECK) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  #可随意加入空白，空白会被忽略。
my $folder = $2;                                                      #把类别标签取出来。
my $name2 = $3;                                                       #把不带路径的文件名取出来。

if (!(-e "$whole/3_STATISTIC_3"))          {mkdir "$whole/3_STATISTIC_2" or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/3_STATISTIC_3/$folder"))  {mkdir "$whole/3_STATISTIC_2/$folder" or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。   
                                        
open(FILE1, "<", "$name1") or die;                                #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/3_STATISTIC_2/$folder/$name2") or die;   #以写的方式打开存放结果的文件。

my @max = ();                              #存储1条序列的PSSM矩阵的每一行的最大值的分布。
for(my $i=0;$i<=100;$i++) {$max[$i]=0;}    #初始化数组。

my $seq_len = 0;                           #最终表示序列长度。
my $element = 0;                           #PSSM矩阵的元素个数。
my $panduan = 0;                           #起判断作用。 
  
while (my $line=<FILE1>) {
  $line =~ m/^[^\n]{420}\n$/ or die;                                     #每一行都匹配这种模式。
  if ($line =~ m/^\s{20}/) {                                             #只有第1行匹配这种模式。故下面的语句只应被执行一次，否则停止执行。
    $panduan==0 or die;
    $panduan=1;
    $line =~ m/ ^\s{20}   ([\sA-Z]{20}){20}   \n$ /x   or die;           #只有第1行匹配这种模式。
  }else{
    $seq_len++; #最终表示序列长度。
    $line =~ m/ ^[0-9\sA-Z]{20}   ([-\s0-9]{20}){20}   \n$ /x   or die;  #除了第1行以外，剩下的所有行都应匹配这种模式。
    $line=~m/(^[ 0-9A-Z]{20})([ 0-9-]{400})\n$/ or die;                  #除了第1行以外，剩下的所有行都应匹配这种模式。
    my $mid1=$2;
    $mid1 =~ s/ ^\s* ([-0-9][^\n]+[0-9]) \s*$ /$1/x or die;              #去掉首尾的空格。
    my @array1=split(/\s+/,$mid1);                                       #把字符串变成数组。
    $#array1==19 or die;                                                 #数组的元素个数应该为20，即最大下标为19.
    for (my $j=0;$j<=$#array1;$j++)  {
      $element++; #PSSM矩阵的元素个数。
      my $temp = $array1[$j];                                             #求出数组的最大值。
      given ($temp) {                                                     #统计各行的最大值，各出现了多少次。
        when ($temp <= -20)  {$max[0]++;}
        when ($temp == -19)  {$max[1]++;}
        when ($temp == -18)  {$max[2]++;}
        when ($temp == -17)  {$max[3]++;}
        when ($temp == -16)  {$max[4]++;}
        when ($temp == -15)  {$max[5]++;}
        when ($temp == -14)  {$max[6]++;}
        when ($temp == -13)  {$max[7]++;}
        when ($temp == -12)  {$max[8]++;}
        when ($temp == -11)  {$max[9]++;}
        when ($temp == -10)  {$max[10]++;}
        when ($temp ==  -9)  {$max[11]++;}
        when ($temp ==  -8)  {$max[12]++;}    
        when ($temp ==  -7)  {$max[13]++;}
        when ($temp ==  -6)  {$max[14]++;}
        when ($temp ==  -5)  {$max[15]++;}
        when ($temp ==  -4)  {$max[16]++;}
        when ($temp ==  -3)  {$max[17]++;}
        when ($temp ==  -2)  {$max[18]++;}
        when ($temp ==  -1)  {$max[19]++;}
        when ($temp ==   0)  {$max[20]++;}
        when ($temp ==   1)  {$max[21]++;}
        when ($temp ==   2)  {$max[22]++;}
        when ($temp ==   3)  {$max[23]++;}
        when ($temp ==   4)  {$max[24]++;}
        when ($temp ==   5)  {$max[25]++;}
        when ($temp ==   6)  {$max[26]++;}
        when ($temp ==   7)  {$max[27]++;}
        when ($temp ==   8)  {$max[28]++;}
        when ($temp ==   9)  {$max[29]++;}
        when ($temp ==  10)  {$max[30]++;}
        when ($temp ==  11)  {$max[31]++;}
        when ($temp ==  12)  {$max[32]++;}    
        when ($temp ==  13)  {$max[33]++;}
        when ($temp ==  14)  {$max[34]++;}
        when ($temp ==  15)  {$max[35]++;}
        when ($temp ==  16)  {$max[36]++;}
        when ($temp ==  17)  {$max[37]++;}
        when ($temp ==  18)  {$max[38]++;}
        when ($temp ==  19)  {$max[39]++;}
        when ($temp >=  20)  {$max[40]++;}
      }
    }                             
  }
}

my $sum = 0; #频数的和
for (my $i=0; $i<=40; $i++) {
  $sum = $sum + $max[$i];
}
($sum == $element)  or die;    #数组@max的所有元素之和应为该序列的PSSM矩阵的元素个数。
($seq_len * 20 == $element)  or die;

my @ratio = ();                   #各最大值的所占的比例。          
my $percent = 0;                  #数组@ratio各元素之和，应该为1。
for (my $i=0; $i<=40; $i++) {     #由频数求频率。  
  $ratio[$i] = $max[$i]/$sum;
  $percent = $percent + $ratio[$i];
}


my @a = ();                              #累积频数。
for(my $i=0;$i<=100;$i++) {$a[$i]=0;}    #初始化数组。

for(my $i=0;$i<=$#max;$i++) {            #求累积频数分布。
  for(my $j=0;$j<=$i;$j++){
    $a[$i]=$a[$i]+$max[$j];
  }
}


my @b = ();                              #累积频率。
for(my $i=0;$i<=100;$i++) {$b[$i]=0;}    #初始化数组。

for(my $i=0;$i<=$#ratio;$i++) {          #求累积分布。
  for(my $j=0;$j<=$i;$j++){
    $b[$i]=$b[$i]+$ratio[$j];
  }
}

#以下是输出：
for(my $i1 = -20; $i1<=20;$i1++) {
  my $i2 = $i1 + 20;
  printf  FILE2  "element == %3s :     %-10s     %-10s     %-10.7f     %-10.7f  \n", $i1, $max[$i2], $a[$i2], $ratio[$i2], $b[$i2];
}

print FILE2 "
频数的和: $sum ,     频率的和:$percent
序列长度: $seq_len
PSSM矩阵的元素个数: $element
\n\n\n\n\n\n\n";    #输出中的7个变量依次表示：最大值，频数，累积频数，频率，累积频率。频数的和，频率的和，序列长度(多少个氨基酸)。

close FILE1;
close FILE2;
} ############################################# END STATISTIC_3_4






##############################################################################################
sub STATISTIC_3_5     #统计每一类序列的PSSM矩阵的每一个元素的分布情况。
##############################################################################################
{ 
my $dir1 = "$whole/2_CHECK";
opendir(my $dh1, $dir1)  ||  die "can't opendir $dir1: $!"; 

while (my $folder = readdir $dh1) {  #读取此目录下的所有目录。
  my $num_seq = 0;                   #每一类所含的序列条数。
  my $aa = 0;                        #每一类的氨基酸总数。
  my $element;                       #每一类PSSM矩阵的元素个数。                       
  my @max = ();
  for(my $i=0;$i<=100;$i++) {$max[$i]=0;}    #初始化数组。

  next unless $folder =~ m/^[0-9-]+$/;                                   #目录名只含类别标签。
  open(FILE1, ">", "$whole/3_STATISTIC_2/$folder-element.dis2") or die;  #以写的方式打开存放结果的文件。
  my $dir2 = "$dir1/$folder";               
  opendir(my $dh2, $dir2)  ||  die "can't opendir $dir2: $!";     #打开目录句柄，读取目录里的文件名。

  while (my $file=readdir $dh2) {               #读取的仅仅是文件名称，不含路径。
    next unless $file !~ m/^[.]/; 
    next unless $file !~ m/[~]$/;
    $num_seq++; #每一类所含的序列条数。
    open(FILE, "<", "$dir2/$file") or die ;     #以读的方式打开待处理文件。
    while (my $line=<FILE>) {
      $line =~ m/^[^\n]{420}\n$/ or die;
      if ($line =~ m/^\s{20}/) {
        $line =~ m/ ^\s{20}   ([\sA-Z]{20}){20}   \n$ /x   or die;
      }else{
        $aa++; #每一类的氨基酸总数。
        $line =~ m/ ^[0-9\sA-Z]{20}   ([-\s0-9]{20}){20}   \n$ /x   or die;
        $line=~m/(^[ 0-9A-Z]{20})([ 0-9-]{400})\n$/ or die; 
        my $mid1=$2;
        $mid1 =~ s/ ^\s* ([-0-9][^\n]+[0-9]) \s*$ /$1/x or die;   #去掉首尾的空格。
        my @array1=split(/\s+/,$mid1);
        $#array1==19 or die;   #数组的元素个数应该为20，即最大下标为19.
        for (my $j=0;$j<=$#array1;$j++)  {
          $element++; #每一类PSSM矩阵的元素个数。
          my $temp = $array1[$j]; 
          given ($temp) {
            when ($temp <= -20)  {$max[0]++;}
            when ($temp == -19)  {$max[1]++;}
            when ($temp == -18)  {$max[2]++;}
            when ($temp == -17)  {$max[3]++;}
            when ($temp == -16)  {$max[4]++;}
            when ($temp == -15)  {$max[5]++;}
            when ($temp == -14)  {$max[6]++;}
            when ($temp == -13)  {$max[7]++;}
            when ($temp == -12)  {$max[8]++;}
            when ($temp == -11)  {$max[9]++;}
            when ($temp == -10)  {$max[10]++;}
            when ($temp ==  -9)  {$max[11]++;}
            when ($temp ==  -8)  {$max[12]++;}    
            when ($temp ==  -7)  {$max[13]++;}
            when ($temp ==  -6)  {$max[14]++;}
            when ($temp ==  -5)  {$max[15]++;}
            when ($temp ==  -4)  {$max[16]++;}
            when ($temp ==  -3)  {$max[17]++;}
            when ($temp ==  -2)  {$max[18]++;}
            when ($temp ==  -1)  {$max[19]++;}
            when ($temp ==   0)  {$max[20]++;}
            when ($temp ==   1)  {$max[21]++;}
            when ($temp ==   2)  {$max[22]++;}
            when ($temp ==   3)  {$max[23]++;}
            when ($temp ==   4)  {$max[24]++;}
            when ($temp ==   5)  {$max[25]++;}
            when ($temp ==   6)  {$max[26]++;}
            when ($temp ==   7)  {$max[27]++;}
            when ($temp ==   8)  {$max[28]++;}
            when ($temp ==   9)  {$max[29]++;}
            when ($temp ==  10)  {$max[30]++;}
            when ($temp ==  11)  {$max[31]++;}
            when ($temp ==  12)  {$max[32]++;}    
            when ($temp ==  13)  {$max[33]++;}
            when ($temp ==  14)  {$max[34]++;}
            when ($temp ==  15)  {$max[35]++;}
            when ($temp ==  16)  {$max[36]++;}
            when ($temp ==  17)  {$max[37]++;}
            when ($temp ==  18)  {$max[38]++;}
            when ($temp ==  19)  {$max[39]++;}
            when ($temp >=  20)  {$max[40]++;}
          }
        }
      }
    }                             
  }#一个文件夹处理完毕。即一类序列处理完毕。

  my $sum = 0;
  for (my $i=0; $i<=40; $i++) {
    $sum = $sum + $max[$i];
  }
  ($sum == $element)  or die;
  ($element == $aa * 20) or die;

  my @ratio = ();
  my $percent = 0;
  for (my $i=0; $i<=40; $i++) {
    $ratio[$i] = $max[$i]/$sum;
    $percent = $percent + $ratio[$i];
  }

  my @a = ();                              #累积频数。
  for(my $i=0;$i<=100;$i++) {$a[$i]=0;}    #初始化数组。

  for(my $i=0;$i<=$#max;$i++) {            #求累积频数分布。
    for(my $j=0;$j<=$i;$j++){
      $a[$i]=$a[$i]+$max[$j];
    }
  }

  my @b = ();
  for(my $i=0;$i<=100;$i++) {$b[$i]=0;}    #初始化数组。

  for(my $i=0;$i<=$#ratio;$i++) {   #求累积分布。
    for(my $j=0;$j<=$i;$j++){
      $b[$i]=$b[$i]+$ratio[$j];
    }
  }

  #以下是输出：
  for(my $i1 = -20; $i1<=20;$i1++) {
    my $i2 = $i1 + 20;
    printf  FILE1  "element == %3s :     %-10s     %-10s     %-10.7f     %-10.7f  \n", $i1, $max[$i2], $a[$i2], $ratio[$i2], $b[$i2];
  }

  print FILE1 "
  频数的和:$sum ,  频率的和:$percent
  每一类的氨基酸数: $aa
  每一类的序列数:$num_seq
  每一类PSSM矩阵的元素个数:$element
  \n\n\n\n\n\n\n";      #输出中的8个变量依次表示：最大值，频数，累积频数，频率，累积频率。频数的和，频率的和，每一类的氨基酸数，每一类的序列数。
}

close FILE;
close FILE1;
}############################################# END STATISTIC_3_5








##############################################################################################
sub  STATISTIC_3_6    #统计所有序列的PSSM矩阵的每一个元素的分布情况。
##############################################################################################
{ 
my @max = ();
for(my $i=0;$i<=100;$i++) {$max[$i]=0;}    #初始化数组。
open(FILE1, ">", "$whole/3_STATISTIC_2/all-element.dis2") or die;   #以写的方式打开存放结果的文件。

chomp(my $dir1 = "$whole/2_CHECK");
opendir DIRHANDLE1,$dir1 or die; 
my $all_element = 0;

while (my $folder = readdir DIRHANDLE1) { #读取此目录下的所有目录。
  next unless $folder =~ m/^[0-9-]+$/;    #目录名只含类别标签。
  my $dir2 = "$dir1/$folder";               
  opendir DIRHANDLE2,$dir2 or die $!;     #打开目录句柄，读取目录里的文件名。
  while (my $file=readdir DIRHANDLE2) {   #读取的仅仅是文件名称，不含路径。
    next unless $file !~ m/^[.]/; 
    next unless $file !~ m/[~]$/;
    open(FILE, "<", "$dir2/$file") or die ;     #以读的方式打开待处理文件。
    while (my $line=<FILE>) {
      $line =~ m/^[^\n]{420}\n$/ or die;
      if ($line =~ m/^\s{20}/) {
        $line =~ m/ ^\s{20}   ([\sA-Z]{20}){20}   \n$ /x   or die;
      }else{
        $line =~ m/ ^[0-9\sA-Z]{20}   ([-\s0-9]{20}){20}   \n$ /x   or die;
        $line=~m/(^[ 0-9A-Z]{20})([ 0-9-]{400})\n$/ or die; 
        my $mid1=$2;
        $mid1 =~ s/ ^\s* ([-0-9][^\n]+[0-9]) \s*$ /$1/x or die;   #去掉首尾的空格。
        my @array1=split(/\s+/,$mid1);
        $#array1==19 or die;   #数组的元素个数应该为20，即最大下标为19.
        for (my $j=0;$j<=$#array1;$j++)  {
          $all_element++;
          my $temp = $array1[$j]; 
          given ($temp) {
            when ($temp <= -20)  {$max[0]++;}
            when ($temp == -19)  {$max[1]++;}
            when ($temp == -18)  {$max[2]++;}
            when ($temp == -17)  {$max[3]++;}
            when ($temp == -16)  {$max[4]++;}
            when ($temp == -15)  {$max[5]++;}
            when ($temp == -14)  {$max[6]++;}
            when ($temp == -13)  {$max[7]++;}
            when ($temp == -12)  {$max[8]++;}
            when ($temp == -11)  {$max[9]++;}
            when ($temp == -10)  {$max[10]++;}
            when ($temp ==  -9)  {$max[11]++;}
            when ($temp ==  -8)  {$max[12]++;}    
            when ($temp ==  -7)  {$max[13]++;}
            when ($temp ==  -6)  {$max[14]++;}
            when ($temp ==  -5)  {$max[15]++;}
            when ($temp ==  -4)  {$max[16]++;}
            when ($temp ==  -3)  {$max[17]++;}
            when ($temp ==  -2)  {$max[18]++;}
            when ($temp ==  -1)  {$max[19]++;}
            when ($temp ==   0)  {$max[20]++;}
            when ($temp ==   1)  {$max[21]++;}
            when ($temp ==   2)  {$max[22]++;}
            when ($temp ==   3)  {$max[23]++;}
            when ($temp ==   4)  {$max[24]++;}
            when ($temp ==   5)  {$max[25]++;}
            when ($temp ==   6)  {$max[26]++;}
            when ($temp ==   7)  {$max[27]++;}
            when ($temp ==   8)  {$max[28]++;}
            when ($temp ==   9)  {$max[29]++;}
            when ($temp ==  10)  {$max[30]++;}
            when ($temp ==  11)  {$max[31]++;}
            when ($temp ==  12)  {$max[32]++;}    
            when ($temp ==  13)  {$max[33]++;}
            when ($temp ==  14)  {$max[34]++;}
            when ($temp ==  15)  {$max[35]++;}
            when ($temp ==  16)  {$max[36]++;}
            when ($temp ==  17)  {$max[37]++;}
            when ($temp ==  18)  {$max[38]++;}
            when ($temp ==  19)  {$max[39]++;}
            when ($temp >=  20)  {$max[40]++;}
          }
        }
      }
    }                             
  }
}

my $sum = 0;
for (my $i=0; $i<=40; $i++) {
  $sum = $sum + $max[$i];
}
($sum == $all_element)  or die;


my @ratio = ();
my $percent = 0;
for (my $i=0; $i<=40; $i++) {
  $ratio[$i] = $max[$i]/$sum;
  $percent = $percent + $ratio[$i];
}

my @a = ();                              #累积频数。
for(my $i=0;$i<=100;$i++) {$a[$i]=0;}    #初始化数组。

for(my $i=0;$i<=$#max;$i++) {            #求累积频数分布。
  for(my $j=0;$j<=$i;$j++){
    $a[$i]=$a[$i]+$max[$j];
  }
}

my @b = ();
for(my $i=0;$i<=100;$i++) {$b[$i]=0;}    #初始化数组。

for(my $i=0;$i<=$#ratio;$i++) {   #求累积分布。
  for(my $j=0;$j<=$i;$j++){
    $b[$i]=$b[$i]+$ratio[$j];
  }
}

#以下是输出：
for(my $i1 = -20; $i1<=20;$i1++) {
  my $i2 = $i1 + 20;
  printf  FILE1  "element == %3s :     %-10s     %-10s     %-10.7f     %-10.7f  \n", $i1, $max[$i2], $a[$i2], $ratio[$i2], $b[$i2];
}

print FILE1 "
频数的和:   $sum ,   频率的和:$percent
所有PSSM矩阵的元素个数:  $all_element
\n\n\n\n\n\n\n";                  #输出中的8个变量依次表示：最大值，频数，累积频数，频率，累积频率。频数的和，频率的和，所有PSSM矩阵的元素个数。

close FILE;
close FILE1;
}############################################# END  STATISTIC_3_6







##############################################################################################
sub STATISTIC_3_7     #统计每一条序列的20种氨基酸的分布。
##############################################################################################
{ 
my $name1 = $_[0];                                                    #接受参数，待处理文件的名字(含路径)，一个文件一条序列。
$name1 =~ m/ (^$whole\/2_CHECK) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  #可随意加入空白，空白会被忽略。
my $folder = $2;                                                      #把类别标签取出来。
my $name2 = $3;                                                       #把不带路径的文件名取出来。

if (!(-e "$whole/3_STATISTIC_3"))          {mkdir "$whole/3_STATISTIC_3" or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/3_STATISTIC_3/$folder"))  {mkdir "$whole/3_STATISTIC_3/$folder" or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。   
                                        
open(FILE1, "<", "$name1") or die;                                #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/3_STATISTIC_3/$folder/$name2") or die;   #以写的方式打开存放结果的文件。

my @max = ();                              #频数数组。
for(my $i=0;$i<=100;$i++) {$max[$i]=0;}    #初始化数组。
my $seq_len = 0;                           #最终表示序列长度。
my $panduan = 0;                           #起判断作用。 
  
while (my $line=<FILE1>) {
  $line =~ m/^[^\n]{420}\n$/ or die;                                     #每一行都匹配这种模式。
  if ($line =~ m/^\s{20}/) {                                             #只有第1行匹配这种模式。故下面的语句只应被执行一次，否则停止执行。
    $panduan==0 or die;
    $panduan=1;
    $line =~ m/ ^\s{20}   ([\sA-Z]{20}){20}   \n$ /x   or die;           #只有第1行匹配这种模式。
  }else{
    $seq_len++;  #频数数组。
    $line =~ m/ ^[0-9\sA-Z]{20}   ([-\s0-9]{20}){20}   \n$ /x   or die;  #除了第1行以外，剩下的所有行都应匹配这种模式。
    $line=~m/^[ 0-9]{6}([A-Z])\s{13}([ 0-9-]{400})\n$/ or die;           #除了第1行以外，剩下的所有行都应匹配这种模式。
    my $mid1=$1;
    $mid1 =~ m/^[A-Z]$/ or die;              #去掉首尾的空格。
    my $temp = $mid1;
    given ($temp) {                          #统计各行的最大值，各出现了多少次。
            when ('A')  {$max[0]++;}
            when ('B')  {$max[1]++;}
            when ('C')  {$max[2]++;}
            when ('D')  {$max[3]++;}
            when ('E')  {$max[4]++;}
            when ('F')  {$max[5]++;}
            when ('G')  {$max[6]++;}
            when ('H')  {$max[7]++;}
            when ('I')  {$max[8]++;}
            when ('J')  {$max[9]++;}
            when ('K')  {$max[10]++;}
            when ('L')  {$max[11]++;}
            when ('M')  {$max[12]++;}    
            when ('N')  {$max[13]++;}
            when ('O')  {$max[14]++;}
            when ('P')  {$max[15]++;}
            when ('Q')  {$max[16]++;}
            when ('R')  {$max[17]++;}
            when ('S')  {$max[18]++;}
            when ('T')  {$max[19]++;}
            when ('U')  {$max[20]++;}
            when ('V')  {$max[21]++;}
            when ('W')  {$max[22]++;}
            when ('X')  {$max[23]++;}
            when ('Y')  {$max[24]++;}
            when ('Z')  {$max[25]++;}
    }                             
  }
}

my $sum = 0;
for (my $i=0; $i<=25; $i++) {
  $sum = $sum + $max[$i];
}
($sum == $seq_len)  or die;    #数组@max的所有元素之和应为该序列的长度。


my @ratio = ();                   #各最大值的所占的比例。          
my $percent = 0;                  #数组@ratio各元素之和，应该为1。
for (my $i=0; $i<=25; $i++) {     #由频数求频率。  
  $ratio[$i] = $max[$i]/$sum;
  $percent = $percent + $ratio[$i];
}


my @a = ();                              #累积频数。
for(my $i=0;$i<=100;$i++) {$a[$i]=0;}    #初始化数组。

for(my $i=0;$i<=$#max;$i++) {            #求累积频数分布。
  for(my $j=0;$j<=$i;$j++){
    $a[$i]=$a[$i]+$max[$j];
  }
}


my @b = ();                              #累积频率。
for(my $i=0;$i<=100;$i++) {$b[$i]=0;}    #初始化数组。

for(my $i=0;$i<=$#ratio;$i++) {          #求累积分布。
  for(my $j=0;$j<=$i;$j++){
    $b[$i]=$b[$i]+$ratio[$j];
  }
}

#以下是输出：
my @aa_temp = qw(A B C D E F G H I J K L M N O P Q R S T U V W X Y Z);
for(my $i1 = 0; $i1<=25;$i1++) {
  my $i2 = $i1;
  printf  FILE2  "AA == %3s :     %-10s     %-10s     %-10.7f     %-10.7f  \n", $aa_temp[$i2], $max[$i2], $a[$i2], $ratio[$i2], $b[$i2];
}

print FILE2 "
频数的和:$sum  ,     频率的和：$percent
序列长度: $seq_len
\n\n\n\n\n\n\n";    #输出中的7个变量依次表示：最大值，频数，累积频数，频率，累积频率。频数的和，频率的和，序列长度(多少个氨基酸)。

close FILE1;
close FILE2;
} ############################################# END STATISTIC_3_7






##############################################################################################
sub STATISTIC_3_8     #统计每一类序列的20种氨基酸的分布情况。
##############################################################################################
{ 
my $dir1 = "$whole/2_CHECK";
opendir(my $dh1, $dir1)  ||  die "can't opendir $dir1: $!"; 

while (my $folder = readdir $dh1) {  #读取此目录下的所有目录。
  my $num_seq = 0;                   #每一类所含的序列条数。
  my $aa = 0;                        #每一类的氨基酸总数。
  my @max = ();
  for(my $i=0;$i<=100;$i++) {$max[$i]=0;}    #初始化数组。
  next unless $folder =~ m/^[0-9-]+$/;                               #目录名只含类别标签。
  open(FILE1, ">", "$whole/3_STATISTIC_3/$folder-aas.dis3") or die;  #以写的方式打开存放结果的文件。
  my $dir2 = "$dir1/$folder";               
  opendir(my $dh2, $dir2)  ||  die "can't opendir $dir2: $!";     #打开目录句柄，读取目录里的文件名。

  while (my $file=readdir $dh2) {               #读取的仅仅是文件名称，不含路径。
    next unless $file !~ m/^[.]/; 
    next unless $file !~ m/[~]$/;
    $num_seq++;  #每一类所含的序列条数。
    open(FILE, "<", "$dir2/$file") or die ;     #以读的方式打开待处理文件。
    while (my $line=<FILE>) {
      $line =~ m/^[^\n]{420}\n$/ or die;
      if ($line =~ m/^\s{20}/) {
        $line =~ m/ ^\s{20}   ([\sA-Z]{20}){20}   \n$ /x   or die;
      }else{
        $aa++; #每一类的氨基酸总数。
        $line =~ m/ ^[0-9\sA-Z]{20}   ([-\s0-9]{20}){20}   \n$ /x   or die;
        $line=~m/^[ 0-9]{6}([A-Z])\s{13}([ 0-9-]{400})\n$/ or die;                  #除了第1行以外，剩下的所有行都应匹配这种模式。
        my $mid1=$1;
        $mid1 =~ m/^[A-Z]$/ or die;              #去掉首尾的空格。
        my $temp = $mid1;
        given ($temp) {                                                      #统计各行的最大值，各出现了多少次。
            when ('A')  {$max[0]++;}
            when ('B')  {$max[1]++;}
            when ('C')  {$max[2]++;}
            when ('D')  {$max[3]++;}
            when ('E')  {$max[4]++;}
            when ('F')  {$max[5]++;}
            when ('G')  {$max[6]++;}
            when ('H')  {$max[7]++;}
            when ('I')  {$max[8]++;}
            when ('J')  {$max[9]++;}
            when ('K')  {$max[10]++;}
            when ('L')  {$max[11]++;}
            when ('M')  {$max[12]++;}    
            when ('N')  {$max[13]++;}
            when ('O')  {$max[14]++;}
            when ('P')  {$max[15]++;}
            when ('Q')  {$max[16]++;}
            when ('R')  {$max[17]++;}
            when ('S')  {$max[18]++;}
            when ('T')  {$max[19]++;}
            when ('U')  {$max[20]++;}
            when ('V')  {$max[21]++;}
            when ('W')  {$max[22]++;}
            when ('X')  {$max[23]++;}
            when ('Y')  {$max[24]++;}
            when ('Z')  {$max[25]++;}
          }
        }
      }
    }                             
  my $sum = 0;
  for (my $i=0; $i<=25; $i++) {
    $sum = $sum + $max[$i];
  }
  ($sum == $aa)  or die;


  my @ratio = ();
  my $percent = 0;
  for (my $i=0; $i<=25; $i++) {
    $ratio[$i] = $max[$i]/$sum;
    $percent = $percent + $ratio[$i];
  }

  my @a = ();                              #累积频数。
  for(my $i=0;$i<=100;$i++) {$a[$i]=0;}    #初始化数组。

  for(my $i=0;$i<=$#max;$i++) {            #求累积频数分布。
    for(my $j=0;$j<=$i;$j++){
      $a[$i]=$a[$i]+$max[$j];
    }
  }

  my @b = ();
  for(my $i=0;$i<=100;$i++) {$b[$i]=0;}    #初始化数组。

  for(my $i=0;$i<=$#ratio;$i++) {   #求累积分布。
    for(my $j=0;$j<=$i;$j++){
      $b[$i]=$b[$i]+$ratio[$j];
    }
  }

  #以下是输出：
  my @aa_temp = qw(A B C D E F G H I J K L M N O P Q R S T U V W X Y Z);
  for(my $i1 = 0; $i1<=25;$i1++) {
    my $i2 = $i1;
    printf  FILE1  "AA == %3s :     %-10s     %-10s     %-10.7f     %-10.7f  \n", $aa_temp[$i2], $max[$i2], $a[$i2], $ratio[$i2], $b[$i2];
  }

  print FILE1 "
  频数的和: $sum     , 频率的和：$percent
  每一类的氨基酸数: $aa
  每一类的序列数:$num_seq
  \n\n\n\n\n\n\n";      #输出中的8个变量依次表示：最大值，频数，累积频数，频率，累积频率。频数的和，频率的和，每一类的氨基酸数，每一类的序列数。
}

close FILE;
close FILE1;
}############################################# END STATISTIC_3_8






##############################################################################################
sub  STATISTIC_3_9    #统计所有序列的20种氨基酸的分布情况。
##############################################################################################
{ 
my @max = ();
for(my $i=0;$i<=100;$i++) {$max[$i]=0;}    #初始化数组。
open(FILE1, ">", "$whole/3_STATISTIC_3/all-aas.dis3") or die;   #以写的方式打开存放结果的文件。

chomp(my $dir1 = "$whole/2_CHECK");
opendir DIRHANDLE1,$dir1 or die; 
my $aa = 0;  #氨基酸总数

while (my $folder = readdir DIRHANDLE1) { #读取此目录下的所有目录。
  next unless $folder =~ m/^[0-9-]+$/;    #目录名只含类别标签。
  my $dir2 = "$dir1/$folder";               
  opendir DIRHANDLE2,$dir2 or die $!;     #打开目录句柄，读取目录里的文件名。
  while (my $file=readdir DIRHANDLE2) {   #读取的仅仅是文件名称，不含路径。
    next unless $file !~ m/^[.]/; 
    next unless $file !~ m/[~]$/;
    open(FILE, "<", "$dir2/$file") or die ;     #以读的方式打开待处理文件。
    while (my $line=<FILE>) {
      $line =~ m/^[^\n]{420}\n$/ or die;
      if ($line =~ m/^\s{20}/) {
        $line =~ m/ ^\s{20}   ([\sA-Z]{20}){20}   \n$ /x   or die;
      }else{
        $aa++; #氨基酸总数
        $line =~ m/ ^[0-9\sA-Z]{20}   ([-\s0-9]{20}){20}   \n$ /x   or die;
        $line=~m/^[ 0-9]{6}([A-Z])\s{13}([ 0-9-]{400})\n$/ or die;                  #除了第1行以外，剩下的所有行都应匹配这种模式。
        my $mid1=$1;
        $mid1 =~ m/^[A-Z]$/ or die;              #去掉首尾的空格。
       my $temp = $mid1;
       given ($temp) {                                                      #统计各行的最大值，各出现了多少次。
            when ('A')  {$max[0]++;}
            when ('B')  {$max[1]++;}
            when ('C')  {$max[2]++;}
            when ('D')  {$max[3]++;}
            when ('E')  {$max[4]++;}
            when ('F')  {$max[5]++;}
            when ('G')  {$max[6]++;}
            when ('H')  {$max[7]++;}
            when ('I')  {$max[8]++;}
            when ('J')  {$max[9]++;}
            when ('K')  {$max[10]++;}
            when ('L')  {$max[11]++;}
            when ('M')  {$max[12]++;}    
            when ('N')  {$max[13]++;}
            when ('O')  {$max[14]++;}
            when ('P')  {$max[15]++;}
            when ('Q')  {$max[16]++;}
            when ('R')  {$max[17]++;}
            when ('S')  {$max[18]++;}
            when ('T')  {$max[19]++;}
            when ('U')  {$max[20]++;}
            when ('V')  {$max[21]++;}
            when ('W')  {$max[22]++;}
            when ('X')  {$max[23]++;}
            when ('Y')  {$max[24]++;}
            when ('Z')  {$max[25]++;}
        }
      }
    }                             
  }
}

my $sum = 0;
for (my $i=0; $i<=25; $i++) {
  $sum = $sum + $max[$i];
}
($sum == $aa)  or die;


my @ratio = ();
my $percent = 0;
for (my $i=0; $i<=25; $i++) {
  $ratio[$i] = $max[$i]/$sum;
  $percent = $percent + $ratio[$i];
}

my @a = ();                              #累积频数。
for(my $i=0;$i<=100;$i++) {$a[$i]=0;}    #初始化数组。

for(my $i=0;$i<=$#max;$i++) {            #求累积频数分布。
  for(my $j=0;$j<=$i;$j++){
    $a[$i]=$a[$i]+$max[$j];
  }
}

my @b = ();
for(my $i=0;$i<=100;$i++) {$b[$i]=0;}    #初始化数组。

for(my $i=0;$i<=$#ratio;$i++) {   #求累积分布。
  for(my $j=0;$j<=$i;$j++){
    $b[$i]=$b[$i]+$ratio[$j];
  }
}

#以下是输出：
my @aa_temp = qw(A B C D E F G H I J K L M N O P Q R S T U V W X Y Z);
for(my $i1 = 0; $i1<=25;$i1++) {
  my $i2 = $i1;
  printf  FILE1  "AA == %3s :     %-10s     %-10s     %-10.7f     %-10.7f  \n", $aa_temp[$i2], $max[$i2], $a[$i2], $ratio[$i2], $b[$i2];
}

print FILE1 "
频数的和: $sum , 频率的和：$percent
氨基酸总数:  $aa
\n\n\n\n\n\n\n";                  #输出中的8个变量依次表示：最大值，频数，累积频数，频率，累积频率。频数的和，频率的和，氨基酸总数。

close FILE;
close FILE1;
}############################################# END  STATISTIC_3_9








##############################################################################################
sub NORMALIZE_4_0     #原样输出，不对数值进行缩放。
##############################################################################################
{ 
my $name1 = $_[0];                                                       #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/2_CHECK) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;     #可随意加入空白，空白会被忽略。 
my $folder = $2;                                                         #把类别标签取出来。
my $name2 = $3;                                                          #把不带路径的文件名取出来。

if (!(-e "$whole/4_NORMALIZE"))         {mkdir "$whole/4_NORMALIZE" or die;}          #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/4_NORMALIZE/$folder")) {mkdir "$whole/4_NORMALIZE/$folder" or die;}  #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。   
                                        
open(FILE1, "<", "$name1") or die;                              #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/4_NORMALIZE/$folder/$name2") or die;   #以写的方式打开存放结果的文件。

while (my $line=<FILE1>) {
  $line =~ m/^[^\n]{420}\n$/ or die;
  if ($line =~ m/^\s{20}/) {
    $line =~ m/ ^\s{20}   ([\sA-Z]{20}){20}   \n$ /x   or die;
  }else{
    $line =~ m/ ^[0-9\sA-Z*]{20}   ([-\s.0-9]{20}){20}   \n$ /x   or die;
  }
  print FILE2 $line;
}                              

close FILE1;
close FILE2;
} ############################################# END NORMALIZE_4_0 






##############################################################################################
sub NORMALIZE_4_1     #缩放，1/(1 + e^(-x))
##############################################################################################
{ 
my $name1 = $_[0];                                                       #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/2_CHECK) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;     #可随意加入空白，空白会被忽略。 
my $folder = $2;                                                         #把类别标签取出来。
my $name2 = $3;                                                          #把不带路径的文件名取出来。

if (!(-e "$whole/4_NORMALIZE"))         {mkdir "$whole/4_NORMALIZE" or die;}          #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/4_NORMALIZE/$folder")) {mkdir "$whole/4_NORMALIZE/$folder" or die;}  #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。   
                                        
open(FILE1, "<", "$name1") or die;                              #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/4_NORMALIZE/$folder/$name2") or die;   #以写的方式打开存放结果的文件。
  
while (my $line=<FILE1>) {
  if ($line=~m/^ {11}/) {
    $line =~ m/^ ([A-Z\s]{20}){21} $/x or die;
    print FILE2 $line;  #输出第一行
  }else{
    $line=~m/ (^[ 0-9A-Z]{20}) ([ 0-9-]{400}) /x or die; 
    print FILE2 $1;  #先输出前20个字符。
    my $mid1=$2;
    $mid1 =~ s/ ^\s* ([-0-9][^\n]+[0-9]) \s*$ /$1/x or die;   #去掉首尾的空格。
    my @array1=split(/\s+/,$mid1);
    $#array1==19 or die;   #数组的元素个数应该为20，即最大下标为19.
    for(my $i1=0;$i1<=$#array1;$i1++){
      $array1[$i1]=1/(1+exp(-$array1[$i1]));  #归一化。
      printf FILE2 "%20.10f",$array1[$i1];
    }
    print FILE2 "\n";
  } 
}                              

close FILE1;
close FILE2;
} ############################################# END NORMALIZE_4_1 







##############################################################################################
sub NORMALIZE_4_2    #缩放，(x - min)/(max - min),按行。
##############################################################################################
{ 
my $name1 = $_[0];                                                       #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/2_CHECK) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;     #可随意加入空白，空白会被忽略。 
my $folder = $2;                                                         #把类别标签取出来。
my $name2 = $3;                                                          #把不带路径的文件名取出来。

if (!(-e "$whole/4_NORMALIZE"))         {mkdir "$whole/4_NORMALIZE" or die;}          #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/4_NORMALIZE/$folder")) {mkdir "$whole/4_NORMALIZE/$folder" or die;}  #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。   
                                        
open(FILE1, "<", "$name1") or die;                              #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/4_NORMALIZE/$folder/$name2") or die;   #以写的方式打开存放结果的文件。

while (my $line=<FILE1>) {
  if ($line=~m/^ {11}/) {
    $line =~ m/^ ([A-Z\s]{20}){21} $/x or die;
    print FILE2 $line;  #输出第一行
  }else{
    $line=~m/ (^[ 0-9A-Z]{20}) ([ 0-9-]{400}) /x or die; 
    print FILE2 $1;  #先输出前20个字符。
    my $mid1=$2;
    $mid1 =~ s/ ^\s* ([-0-9][^\n]+[0-9]) \s*$ /$1/x or die;   #去掉首尾的空格。
    my @array1=split(/ +/,$mid1);
    $#array1==19 or die; #数组的元素个数应该为20，即最大下标为19.
    my $max1=&MAX(@array1);
    my $min1=&MIN(@array1);
    my $i1=0;
    ($max1 >= $min1) or die;
    for($i1=0;$i1<=$#array1;$i1++){
      $array1[$i1] = ($array1[$i1] - $min1)/($max1 - $min1 + 10**(-20));  #归一化,把分母加上一个很小的数，防止除数为0.   
      printf FILE2 "%20.10f",$array1[$i1];	
    }
    print FILE2 "\n";
  } 
}                              

close FILE1;
close FILE2;
} ######################################## END NORMALIZE_4_2







##############################################################################################
sub NORMALIZE_4_3    #缩放，(x - min)/(max - min),按列。20中氨基酸分别考虑。
##############################################################################################
{ 
my $name1 = $_[0];                                                       #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/2_CHECK) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;     #可随意加入空白，空白会被忽略。 
my $folder = $2;                                                         #把类别标签取出来。
my $name2 = $3;                                                          #把不带路径的文件名取出来。

if (!(-e "$whole/4_NORMALIZE"))         {mkdir "$whole/4_NORMALIZE" or die;}          #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/4_NORMALIZE/$folder")) {mkdir "$whole/4_NORMALIZE/$folder" or die;}  #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。   
                                        
open(FILE1, "<", "$name1") or die;                              #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/4_NORMALIZE/$folder/$name2") or die;   #以写的方式打开存放结果的文件。
                      
my @lines = <FILE1>;
my @array2d = '';                                                  #声明一个数组,此时数组的第一个元素是有值的，为空值。若使用“my @array2d = ();”，则$array2d[0]的值也为undef,即所有元素都是没有值的。
for (my $i1=1;$i1<=$#lines;$i1++) {                                #数组@lines的第一个元素不做处理。
  $lines[$i1] =~ s/ ^\s* ([-0-9][^\n]+[0-9]) \s*$ /$1/x or die;    #去掉首尾的空格。
  $array2d[$i1] = [split(/ +/,$lines[$i1])];                       #用方括号生成匿名引用。
  $#{$array2d[$i1]} == 21 or die;                                  #看成一维数组时，其每个元素均为22维的数组。
}	                                                           #现在@array2d是一个L*22维的2维数组。L是序列长度。

my @table = qw(A R N D C Q E G H I L K M F P S T W Y V);	#20种氨基酸。
my @temp1 = ();                                                 #所有元素都是没有值的。
my @temp2 = ();                                                 #所有元素都是没有值的。

for my $letter (@table) {                              # $letter的一个值对应PSSM中的某些行。
  @temp1 = ();                                         # @temp1的所有元素决定了一种氨基酸在蛋白质序列中的那些位置上出现了。
  for (my $i=1;$i<=$#lines;$i++) {                     # @temp1是从第1个元素开始赋值的，即从$temp1[0]开始。
    if ($array2d[$i][1] eq $letter) {                  #把不同位置的同种氨基酸提出来。
      $temp1[++$#temp1] = $i;                          #把行号（位置号）存起来,蛋白质序列的这些位置上出现了字母$letter。  
    }                                                  #@temp1可能只有一个元素，这样按列归一化时，每个元素均为0.即某个字母在氨基酸序列中只出现了一次，很有可能是半胱氨酸(C)。
  }                                                    #此时，$temp1[0]的值也还有可能为undef,这样运行时会出现警告信息。即某个字母在氨基酸序列中一次也没有出现。  
  for (my $j=2;$j<=21;$j++) {
    @temp2 = ();
    foreach my $k (@temp1){                            #@temp1没有初始化时，此语句不会被执行。
      $temp2[++$#temp2] = $array2d[$k][$j];            #同一列的同种氨基酸。
    }
    my $min = &MIN(@temp2);                            #@temp1没有初始化时，@temp2也没有初始化，故$min也没有初始化，但运行时不会显示警告信息，因为还没有使用它。
    my $max = &MAX(@temp2);                            #@temp1没有初始化时，@temp2也没有初始化，故$max也没有初始化，但运行时不会显示警告信息，因为还没有使用它。  
    foreach my $n (@temp1){                            #@temp1没有初始化时，此语句不会被执行。故不会因为$min和$max没有初始化而显示警告信息。
      $array2d[$n][$j] = ($array2d[$n][$j] - $min) / ($max -$min+10**(-20));     #归一化。
    }
  }
}

#归一化完毕，下面是输出:
$lines[0] =~ m/^ ([A-Z\s]{20}){21} $/x or die;
print FILE2 $lines[0];  #输出第一行

for (my $i=1;$i<=$#lines;$i++) {  #从第2行开始，下标为1.
  printf FILE2 "%5s",$array2d[$i][0];
  printf FILE2 "%2s             ",$array2d[$i][1];	#占用的是15个位置。
  for (my $j=2;$j<=21;$j++) {
    printf FILE2 "%20.10f",$array2d[$i][$j];
  }
  print FILE2 "\n";
}

close FILE1;
close FILE2;
}##################################################### END NORMALIZE_4_3







##############################################################################################
sub NORMALIZE_4_4    #缩放，(x - min)/(max - min),按列。#最好是按这种方法归一化。
##############################################################################################
{ 
my $name1 = $_[0];                                                          #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/2_CHECK) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;        #可随意加入空白，空白会被忽略。 
my $folder = $2;                                                            #把类别标签取出来。
my $name2 = $3;                                                             #把不带路径的文件名取出来。

if (!(-e "$whole/4_NORMALIZE"))         {mkdir "$whole/4_NORMALIZE" or die;}          #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/4_NORMALIZE/$folder")) {mkdir "$whole/4_NORMALIZE/$folder" or die;}  #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。   
                                        
open(FILE1, "<", "$name1") or die;                              #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/4_NORMALIZE/$folder/$name2") or die;   #以写的方式打开存放结果的文件。
                      
my @lines = <FILE1>;
my @array2d = '';                                                  #声明一个数组,此时数组的第一个元素是有值的，为空值。若使用“my @array2d = ();”，则$array2d[0]的值也为undef,即所有元素都是没有值的。
for (my $i1=1;$i1<=$#lines;$i1++) {                                #数组@lines的第一个元素不做处理。
  $lines[$i1] =~ s/ ^\s* ([-0-9][^\n]+[0-9]) \s*$ /$1/x or die;    #去掉首尾的空格。
  $array2d[$i1] = [split(/ +/,$lines[$i1])];                       #用方括号生成匿名引用。
  $#{$array2d[$i1]} == 21 or die;                                  #看成一维数组时，其每个元素均为22维的数组。
}	                                                           #现在@array2d是一个L*22维的2维数组。L是序列长度。

my @table = qw(A R N D C Q E G H I L K M F P S T W Y V);	#20种氨基酸。

for (my $j=2;$j<=21;$j++) {
  my @temp2 = ();
  for(my $k=1; $k<=$#lines; $k++){                            #@temp1没有初始化时，此语句不会被执行。
    $temp2[++$#temp2] = $array2d[$k][$j];            #同一列的氨基酸。
  }
  my $min = &MIN(@temp2);                            #@temp1没有初始化时，@temp2也没有初始化，故$min也没有初始化，但运行时不会显示警告信息，因为还没有使用它。
  my $max = &MAX(@temp2);                            #@temp1没有初始化时，@temp2也没有初始化，故$max也没有初始化，但运行时不会显示警告信息，因为还没有使用它。  
  for(my $n=1; $n<=$#lines; $n++){                   #@temp1没有初始化时，此语句不会被执行。故不会因为$min和$max没有初始化而显示警告信息。
    $array2d[$n][$j] = ($array2d[$n][$j] - $min) / ($max -$min+10**(-20));     #归一化。
  }
}

#归一化完毕，下面是输出:
$lines[0] =~ m/^ ([A-Z\s]{20}){21} $/x or die;
print FILE2 $lines[0];  #输出第一行

for (my $i=1;$i<=$#lines;$i++) {  #从第2行开始，下标为1.
  printf FILE2 "%5s",$array2d[$i][0];
  printf FILE2 "%2s             ",$array2d[$i][1];	#占用的是15个位置。
  for (my $j=2;$j<=21;$j++) {
    printf FILE2 "%20.10f",$array2d[$i][$j];
  }
  print FILE2 "\n";
}

close FILE1;
close FILE2;
}##################################################### END NORMALIZE_4_4








open(AVE_FILE, ">", "$whole/Ave-Ratio.txt") or die ;
my $ave_ratio = 0;   #所有序列平均被保留了多少行。
my $all_seq = 0;     #总共多少条序列。
my $all_line = 0;    #总共多少行。
my $delete_line = 0; #总共过滤掉了多少行。
my $remain_line = 0; #总共保留了多少行。
##############################################################################################
sub FILTER_5_1  #按列过滤（即按单个元素过滤），大于阈值的留下。
##############################################################################################
{ 
my $name1 = $_[0];               #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/4_NORMALIZE) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;                 #把类别标签取出来。
my $name2 = $3;                  #把不带路径的文件名取出来。

if (!(-e "$whole/5_FILTER"))          {mkdir "$whole/5_FILTER"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/5_FILTER/$folder"))  {mkdir "$whole/5_FILTER/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";                             #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/5_FILTER/$folder/$name2") or die "Can't open  > $name2 : $!";    #以写的方式打开存放结果的文件。
open(FILE3, ">>", "$whole/5_FILTER/ratio") or die "Can't open  >> ratio : $!";            #存储每一条序列保留了多少行。  
  
my $threshold = $ARGV[7];  ########################################设定的阈值。

my $x1 = 0;   #总共多少个元素。
my $x2 = 0;   #保留了多少个元素。
my $x3 = 0;   #过滤掉了多少个元素。

my @lines = <FILE1>;

for (my $i=0;$i<=$#lines;$i++) {
  if ($lines[$i]=~m/^ {11}/) {
    $i==0 or die;
    $lines[$i] =~ m/^ ([A-Z\s]{20}){21} $/x or die;
    print FILE2 $lines[$i];  #输出第一行
  }else{
    $lines[$i]=~m/ (^[ 0-9A-Z*]{20}) ([ \.0-9-]{400}) /x or die; 
    my $temp = $1;
    my $mid1=$2;
    $mid1 =~ s/ ^\s* ([-0-9][^\n]+[0-9]) \s*$/$1/x or die;   #去掉首尾的空格。
    my @array1=split(/ +/,$mid1);
    $#array1==19 or die;  #应为20个元素。
    print FILE2 $temp;    #输出一行的前20个字符。
    for (my $j=0;$j<=$#array1;$j++) {
      $x1++; #总共多少个元素。
      $all_line++;
      if ($array1[$j] >= $threshold) { 
        $x2++;  #保留了多少个元素。
        $remain_line++;
        printf FILE2 "%20s", $array1[$j] ;   #元素的值大于阈值则输出，否则输出0。
      }else{
        $array1[$j] =~ m/[-0-9]+/ or die ;
        $array1[$j] =~ s/[0-9-.]+/0.0/ or die ;
        $array1[$j] =~ m/[0.]{3}/ or die ;
        $x3++; #过滤掉了多少个元素。
        $delete_line++;
        printf FILE2 "%20.1f", $array1[$j] ; 
      }
    }
    print FILE2 "\n";
  }
}                             
 
my $ratio = $x2/$x1 ;
$ratio = $ratio * 100;
($x1 == $x2+$x3)   or die;              
print FILE3  "$name1:\n总共$x1个元素。\n保留了$x2个元素。\n过滤掉了$x3个元素。\n保留比例为：$ratio%  \n\n ";

$ave_ratio = $ave_ratio + $ratio;  #各条序列保留比例的和。
$all_seq++; #总共多少条序列。

close FILE1;
close FILE2;
} ################################################### END FILTER_5_1






##############################################################################################
sub FILTER_5_2    #按列过滤（即按单个元素过滤），小于阈值的留下。
##############################################################################################
{ 
my $name1 = $_[0];               #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/4_NORMALIZE) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;                 #把类别标签取出来。
my $name2 = $3;                  #把不带路径的文件名取出来。

if (!(-e "$whole/5_FILTER"))          {mkdir "$whole/5_FILTER"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/5_FILTER/$folder"))  {mkdir "$whole/5_FILTER/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";                             #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/5_FILTER/$folder/$name2") or die "Can't open  > $name2 : $!";    #以写的方式打开存放结果的文件。
open(FILE3, ">>", "$whole/5_FILTER/ratio") or die "Can't open  >> ratio : $!";            #存储每一条序列保留了多少行。
 
my $threshold = $ARGV[7];  ########################################设定的阈值。

my $x1 = 0;   #总共多少个元素。
my $x2 = 0;   #保留了多少个元素。
my $x3 = 0;   #过滤掉了多少个元素。

my @lines = <FILE1>;

for (my $i=0;$i<=$#lines;$i++) {
  if ($lines[$i]=~m/^ {11}/) {
    $i==0 or die;
    $lines[$i] =~ m/^ ([A-Z\s]{20}){21} $/x or die;
    print FILE2 $lines[$i]++;  #输出第一行
  }else{
    $lines[$i]=~m/(^[ 0-9A-Z]{20})([ \.0-9-]{400})/ or die; 
    my $temp = $1;
    my $mid1=$2;
    $mid1 =~ s/^\s*([-0-9][^\n]+[0-9])\s*$/$1/ or die;   #去掉首尾的空格。
    my @array1=split(/ +/,$mid1);
    $#array1==19 or die;  #应为20个元素。
    print FILE2 $temp;  #输出一行的前20个字符。
    for (my $j=0;$j<=$#array1;$j++) {
      $x1++; #总共多少个元素。
      $all_line++;
      if ($array1[$j] <= $threshold) { 
        $x2++; #保留了多少个元素。
        $remain_line++;
        printf FILE2 "%20s", $array1[$j] ;   #元素的值小于阈值则输出，否则输出0。
      }else{
        $array1[$j] =~ m/[0-9-]+/ or die ;
        $array1[$j] =~ s/[0-9-.]+/0.0/ or die ;
        $array1[$j] =~ m/[0.]{3}/ or die ;
        $x3++; #过滤掉了多少个元素。
        $delete_line++;
        printf FILE2 "%20.1f", $array1[$j] ; 
      }
    }
    print FILE2 "\n";
  }
}                             
 
my $ratio = $x2/$x1 ;
$ratio = $ratio * 100;  
($x1 == $x2+$x3)   or die;            
print FILE3  "$name1:\n总共$x1个元素。\n保留了$x2个元素。\n过滤掉了$x3个元素。\n保留比例为：$ratio%  \n\n ";

$ave_ratio = $ave_ratio + $ratio;  #各条序列保留比例的和。
$all_seq++;  #总共多少条序列。

close FILE1;
close FILE2;
}################################################### END FILTER_5_2






##############################################################################################
sub FILTER_5_3   #把某些行过滤掉，考虑每一行的最大值(max)，在长度为m的窗口内若有n个大于阈值，则留下。(过滤以后的行不再输出）（n可以等于m,下同）
##############################################################################################
{ 
my $name1 = $_[0];               #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/4_NORMALIZE) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;                 #把类别标签取出来。
my $name2 = $3;                  #把不带路径的文件名取出来。

if (!(-e "$whole/5_FILTER"))          {mkdir "$whole/5_FILTER"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/5_FILTER/$folder"))  {mkdir "$whole/5_FILTER/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";                             #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/5_FILTER/$folder/$name2") or die "Can't open  > $name2 : $!";    #以写的方式打开存放结果的文件。
open(FILE3, ">>", "$whole/5_FILTER/ratio") or die "Can't open  >> ratio : $!";            #存储每一条序列保留了多少行。

my $window1 = $ARGV[3];
my $window2 = $ARGV[5];
my $threshold = $ARGV[7];  ########################################设定的阈值。


my $x1 = 0;   #总共多少行。
my $x2 = 0;   #保留了多少行。
my $x3 = 0;   #过滤掉了多少行。

my @lines = <FILE1>;
my @array_max = ();

for (my $i=0;$i<=$#lines;$i++) {
  if ($lines[$i]=~m/^ {11}/) {
    $i==0 or die;
    $lines[$i] =~ m/^ ([A-Z\s]{20}){21} $/x or die;
    print FILE2 $lines[$i];  #输出第一行
  }else{
    $lines[$i]=~m/ (^[ 0-9A-Z]{20}) ([ 0-9-]{400}) /x or die; 
    my $mid1=$2;
    $mid1 =~ s/ ^\s* ([-0-9][^\n]+[0-9]) \s*$ /$1/x or die;   #去掉首尾的空格。
    my @array1=split(/ +/,$mid1);
    $#array1==19 or die;  #应为20个元素。
    my $max = &MAX(@array1);
    $array_max[$i] = $max;     #各行的最大值都在这里面。$array_max[0]的值应为undef.
  }
}                             
$#lines == $#array_max  or  die;  # $array_max[0]的值应为undef。

my @sites = ();  #保留下来的位点其元素值变为1.
for(my $i=0; $i<=$#array_max; $i++) {$sites[$i] = 0; }

for (my $i1=1; $i1<=$#array_max - $window1 + 1; $i1++) { #每一次循环考虑一个窗口
  my $p = 0;
  for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) {  #每一次循环考虑一个位点。
    if($array_max[$i2]>=$threshold) { $p++; }
  }
  if($p >= $window2) { #若满足条件则保留。
    for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) { 
      $sites[$i2] = 1;
    } 
  }
}

#输出:
for(my $i=1;$i<=$#array_max;$i++) {
  $x1++; #总共多少行。
  $all_line++;
  if($sites[$i] == 1){
    print FILE2 $lines[$i];
    $x2++; #保留了多少行。
    $remain_line++; 
  }else{
    $x3++;#过滤掉了多少行。
    $delete_line++;
  }
}
       
my $ratio = $x2/$x1 ;
$ratio = $ratio * 100;   
($x1 == $x2+$x3)   or die;           
print FILE3  "$name1:\n总共$x1行。\n保留了$x2行。\n过滤掉了$x3行。\n保留比例为：$ratio%  \n\n ";

$ave_ratio = $ave_ratio + $ratio;
$all_seq++;

close FILE1;
close FILE2;
}################################################### END FILTER_5_3





##############################################################################################
sub FILTER_5_4   #把某些行过滤掉，考虑每一行的最大值(max)，在长度为m的窗口内若有n个小于阈值，则留下。
##############################################################################################
{ 
my $name1 = $_[0];               #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/4_NORMALIZE) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;                 #把类别标签取出来。
my $name2 = $3;                  #把不带路径的文件名取出来。

if (!(-e "$whole/5_FILTER"))          {mkdir "$whole/5_FILTER"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/5_FILTER/$folder"))  {mkdir "$whole/5_FILTER/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";                             #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/5_FILTER/$folder/$name2") or die "Can't open  > $name2 : $!";    #以写的方式打开存放结果的文件。
open(FILE3, ">>", "$whole/5_FILTER/ratio") or die "Can't open  >> ratio : $!";            #存储每一条序列保留了多少行。

my $window1 = $ARGV[3];
my $window2 = $ARGV[5];
my $threshold = $ARGV[7];  ########################################设定的阈值。


my $x1 = 0;   #总共多少行。
my $x2 = 0;   #保留了多少行。
my $x3 = 0;   #过滤掉了多少行。

my @lines = <FILE1>;
my @array_max = ();

for (my $i=0;$i<=$#lines;$i++) {
  if ($lines[$i]=~m/^ {11}/) {
    $i==0 or die;
    $lines[$i] =~ m/^ ([A-Z\s]{20}){21} $/x or die;
    print FILE2 $lines[$i];  #输出第一行
  }else{
    $lines[$i]=~m/ (^[ 0-9A-Z]{20}) ([ 0-9-]{400}) /x or die; 
    my $mid1=$2;
    $mid1 =~ s/ ^\s* ([-0-9][^\n]+[0-9]) \s*$ /$1/x or die;   #去掉首尾的空格。
    my @array1=split(/ +/,$mid1);
    $#array1==19 or die;  #应为20个元素。
    my $max = &MAX(@array1);
    $array_max[$i] = $max;     #各行的最大值都在这里面。$array_max[0]的值应为undef.
  }
}                             
$#lines == $#array_max  or  die;  # $array_max[0]的值应为undef。

my @sites = ();  #保留下来的位点其元素值变为1.
for(my $i=0; $i<=$#array_max; $i++) {$sites[$i] = 0; }

for (my $i1=1; $i1<=$#array_max - $window1 + 1; $i1++) { #每一次循环考虑一个窗口
  my $p = 0;
  for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) {  #每一次循环考虑一个位点。
    if($array_max[$i2]<$threshold) { $p++; }
  }
  if($p >= $window2) { #若满足条件则保留。
    for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) { 
      $sites[$i2] = 1;
    } 
  }
}

#输出:
for(my $i=1;$i<=$#array_max;$i++) {
  $x1++; #总共多少行。
  $all_line++;
  if($sites[$i] == 1){
    print FILE2 $lines[$i];
    $x2++; #保留了多少行。
    $remain_line++; 
  }else{
    $x3++;#过滤掉了多少行。
    $delete_line++;
  }
}
 
       
my $ratio = $x2/$x1 ;
$ratio = $ratio * 100;   
($x1 == $x2+$x3)   or die;           
print FILE3  "$name1:\n总共$x1行。\n保留了$x2行。\n过滤掉了$x3行。\n保留比例为：$ratio%  \n\n ";

$ave_ratio = $ave_ratio + $ratio;
$all_seq++;

close FILE1;
close FILE2;
}################################################### END FILTER_5_4







##############################################################################################
sub FILTER_5_5   #把某些行过滤掉，考虑连续多少行的最大值的平均值，若大于阈值，则留下。
##############################################################################################
{ 
my $name1 = $_[0];               #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/4_NORMALIZE) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;                 #把类别标签取出来。
my $name2 = $3;                  #把不带路径的文件名取出来。

if (!(-e "$whole/5_FILTER"))          {mkdir "$whole/5_FILTER"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/5_FILTER/$folder"))  {mkdir "$whole/5_FILTER/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";                             #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/5_FILTER/$folder/$name2") or die "Can't open  > $name2 : $!";    #以写的方式打开存放结果的文件。
open(FILE3, ">>", "$whole/5_FILTER/ratio") or die "Can't open  >> ratio : $!";            #存储每一条序列保留了多少行。
open(FILE8, ">>", "$whole/5_FILTER/max-ave") or die ;                                     #存储每一条序列的PSSM的某些行的最大值的平均值。

my $window1 = $ARGV[3];
my $threshold = $ARGV[7];  ########################################设定的阈值。

my $x1 = 0;   #总共多少行。
my $x2 = 0;   #保留了多少行。
my $x3 = 0;   #过滤掉了多少行。

my @lines = <FILE1>;
print FILE2 $lines[0];
my @array_max = ();

for (my $i=1;$i<=$#lines;$i++) {   #第一行$lines[0]不用处理。
    $lines[$i]=~m/ (^[ 0-9A-Z]{20}) ([ .0-9-]{400}) /x or die; 
    my $mid1=$2;
    $mid1 =~ s/ ^\s* ([-0-9][^\n]+[0-9]) \s*$ /$1/x or die;   #去掉首尾的空格。
    my @array1=split(/\s+/,$mid1);
    $#array1==19 or die;  #应为20个元素。
    my $max = &MAX(@array1);
    $array_max[$i] = $max;     #各行的最大值都在这里面。$array_max[0]的值应为undef.
}                         
$#lines == $#array_max  or  die;  # $array_max[0]的值应为undef.

my @sites = ();  #保留下来的位点其元素值变为1.
for(my $i=0; $i<=$#array_max; $i++) {$sites[$i] = 0; }

for (my $i1=1; $i1<=$#array_max - $window1 + 1; $i1++) { #每一次循环考虑一个窗口
  my $p = 0;
  for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) {  #每一次循环考虑一个位点。
    $p = $p + $array_max[$i2]; 
  }
  $p = $p/$window1 ;
  if($p >= $threshold) { #若满足条件则保留。
    for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) { 
      $sites[$i2] = 1;
    } 
  }
}

#输出:
for(my $i=1;$i<=$#array_max;$i++) {
  $x1++; #总共多少行。
  $all_line++;
  if($sites[$i] == 1){
    print FILE2 $lines[$i];
    $x2++; #保留了多少行。
    $remain_line++; 
  }else{
    $x3++;#过滤掉了多少行。
    $delete_line++;
  }
}
 

my $ratio = $x2/$x1 ;
$ratio = $ratio * 100; 
($x1 == $x2+$x3)   or die;             
print FILE3  "$name1:\n总共$x1行。\n保留了$x2行。\n过滤掉了$x3行。\n保留比例为：$ratio%  \n\n ";

$ave_ratio = $ave_ratio + $ratio;
$all_seq++;

close FILE1;
close FILE2;
close FILE3;
} ################################################### END FILTER_5_5




##############################################################################################
sub FILTER_5_6   #把某些行过滤掉，考虑连续多少行的最大值的平均值，若大于阈值，则留下。
##############################################################################################
{ 
my $name1 = $_[0];               #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/4_NORMALIZE) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;                 #把类别标签取出来。
my $name2 = $3;                  #把不带路径的文件名取出来。

if (!(-e "$whole/5_FILTER"))          {mkdir "$whole/5_FILTER"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/5_FILTER/$folder"))  {mkdir "$whole/5_FILTER/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";                             #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/5_FILTER/$folder/$name2") or die "Can't open  > $name2 : $!";    #以写的方式打开存放结果的文件。
open(FILE3, ">>", "$whole/5_FILTER/ratio") or die "Can't open  >> ratio : $!";            #存储每一条序列保留了多少行。
open(FILE8, ">>", "$whole/5_FILTER/max-ave") or die ;                                     #存储每一条序列的PSSM的某些行的最大值的平均值。

my $window1 = $ARGV[3];
my $threshold = $ARGV[7];  ########################################设定的阈值。

my $x1 = 0;   #总共多少行。
my $x2 = 0;   #保留了多少行。
my $x3 = 0;   #过滤掉了多少行。

my @lines = <FILE1>;
print FILE2 $lines[0];
my @array_max = ();

for (my $i=1;$i<=$#lines;$i++) {   #第一行$lines[0]不用处理。
    $lines[$i]=~m/ (^[ 0-9A-Z]{20}) ([ .0-9-]{400}) /x or die; 
    my $mid1=$2;
    $mid1 =~ s/ ^\s* ([-0-9][^\n]+[0-9]) \s*$ /$1/x or die;   #去掉首尾的空格。
    my @array1=split(/\s+/,$mid1);
    $#array1==19 or die;  #应为20个元素。
    my $max = &MAX(@array1);
    $array_max[$i] = $max;     #各行的最大值都在这里面。$array_max[0]的值应为undef.
}                         
$#lines == $#array_max  or  die;  # $array_max[0]的值应为undef.

my @sites = ();  #保留下来的位点其元素值变为1.
for(my $i=0; $i<=$#array_max; $i++) {$sites[$i] = 0; }

for (my $i1=1; $i1<=$#array_max - $window1 + 1; $i1++) { #每一次循环考虑一个窗口
  my $p = 0;
  for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) {  #每一次循环考虑一个位点。
    $p = $p + $array_max[$i2]; 
  }
  $p = $p/$window1 ;
  printf  FILE8  "%-5s  \n", $p;         #输出平均值。
  if($p < $threshold) { #若满足条件则保留。
    for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) { 
      $sites[$i2] = 1;
    } 
  }
}

#输出:
for(my $i=1;$i<=$#array_max;$i++) {
  $x1++; #总共多少行。
  $all_line++;
  if($sites[$i] == 1){
    print FILE2 $lines[$i];
    $x2++; #保留了多少行。
    $remain_line++; 
  }else{
    $x3++;#过滤掉了多少行。
    $delete_line++;
  }
}
 
my $ratio = $x2/$x1 ;
$ratio = $ratio * 100; 
($x1 == $x2+$x3)   or die;             
print FILE3  "$name1:\n总共$x1行。\n保留了$x2行。\n过滤掉了$x3行。\n保留比例为：$ratio%  \n\n ";

$ave_ratio = $ave_ratio + $ratio;
$all_seq++;

close FILE1;
close FILE2;
close FILE3;
}################################################### END FILTER_5_6








##############################################################################################
sub FILTER_5_7   #把某些行过滤掉，考虑20种氨基酸的背景频率，不同类氨基酸应该采用不同阈值，在长度为m的窗口内若有n个大于阈值，则留下。
##############################################################################################
{ 
my $name1 = $_[0];               #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/4_NORMALIZE) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;                 #把类别标签取出来。
my $name2 = $3;                  #把不带路径的文件名取出来。

if (!(-e "$whole/5_FILTER"))          {mkdir "$whole/5_FILTER"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/5_FILTER/$folder"))  {mkdir "$whole/5_FILTER/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";                             #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/5_FILTER/$folder/$name2") or die "Can't open  > $name2 : $!";    #以写的方式打开存放结果的文件。
open(FILE3, ">>", "$whole/5_FILTER/ratio") or die "Can't open  >> ratio : $!";            #存储每一条序列保留了多少行。

my $other = 0.00;
my %hash = (A => 0.0826+$other,
            C => 0.0136+$other,
            D => 0.0545+$other,
            E => 0.0675+$other,
            F => 0.0386+$other,
            G => 0.0708+$other,
            H => 0.0227+$other,
            I => 0.0597+$other,
            K => 0.0584+$other,
            L => 0.0966+$other,
            M => 0.0242+$other,
            N => 0.0406+$other,
            P => 0.0470+$other,
            Q => 0.0393+$other,
            R => 0.0553+$other,
            S => 0.0655+$other,
            T => 0.0534+$other,
            V => 0.0687+$other,
            W => 0.0108+$other,
            Y => 0.0292+$other,
           );
for my $key (keys(%hash)) {
  $hash{$key} = $hash{$key} / 0.01;
  $hash{$key} = 10 * &LOGXY(10,$hash{$key});           #需要把背景频率值转化一下，根据PSSM的构造规则转化。

}


my $window1 = $ARGV[3];
my $window2 = $ARGV[5];
my $threshold = 0;  #阈值随位点的改变而改变。


my $x1 = 0;   #总共多少行。
my $x2 = 0;   #保留了多少行。
my $x3 = 0;   #过滤掉了多少行。

my @lines = <FILE1>;
my @array_max = ();
my @all_sites = ();

for (my $i=0;$i<=$#lines;$i++) {
  if ($lines[$i]=~m/^ {11}/) {
    $i==0 or die;
    $lines[$i] =~ m/^ ([A-Z\s]{20}){21} $/x or die;
    print FILE2 $lines[$i];  #输出第一行
  }else{
    $lines[$i]=~m/ ^[\s0-9]{6}([A-Z])\s{13} ([ 0-9-]{400}) /x or die; 
    $all_sites[$i] = $1;
    my $mid1=$2;
    $all_sites[$i] =~ m/^[A-Z]$/ or die; 
    $mid1 =~ s/ ^\s* ([-0-9][^\n]+[0-9]) \s*$ /$1/x or die;   #去掉首尾的空格。
    my @array1=split(/ +/,$mid1);
    $#array1==19 or die;  #应为20个元素。
    my $max = &MAX(@array1);
    $array_max[$i] = $max;     #各行的最大值都在这里面。$array_max[0]的值应为undef.
  }
}                             
$#lines == $#array_max  or  die;  # $array_max[0]的值应为undef。

my @sites = ();  #保留下来的位点其元素值变为1.
for(my $i=0; $i<=$#array_max; $i++) {$sites[$i] = 0; }

for (my $i1=1; $i1<=$#array_max - $window1 + 1; $i1++) { #每一次循环考虑一个窗口
  given($all_sites[$i1]) {  #阈值随位点的改变而改变。
    when('A') {$threshold = $hash{A}; }
    when('C') {$threshold = $hash{C}; }
    when('D') {$threshold = $hash{D}; }
    when('E') {$threshold = $hash{E}; }
    when('F') {$threshold = $hash{F}; }
    when('G') {$threshold = $hash{G}; }
    when('H') {$threshold = $hash{H}; }
    when('I') {$threshold = $hash{I}; }
    when('K') {$threshold = $hash{K}; }
    when('L') {$threshold = $hash{L}; }
    when('M') {$threshold = $hash{M}; }
    when('N') {$threshold = $hash{N}; }
    when('P') {$threshold = $hash{P}; }
    when('Q') {$threshold = $hash{Q}; }
    when('R') {$threshold = $hash{R}; }
    when('S') {$threshold = $hash{S}; }
    when('T') {$threshold = $hash{T}; }
    when('V') {$threshold = $hash{V}; }
    when('W') {$threshold = $hash{W}; }
    when('Y') {$threshold = $hash{Y}; }
    when($all_sites[$i1] =~ m/[UOBZJX]/) {$threshold = 0; }
    default   {die;        }
  }
  my $p = 0;
  for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) {  #每一次循环考虑一个位点。
    if($array_max[$i2]>=$threshold) { $p++; }
  }
  if($p >= $window2) { #若满足条件则保留。
    for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) { 
      $sites[$i2] = 1;
    } 
  }
}

#输出:
for(my $i=1;$i<=$#array_max;$i++) {
  $x1++; #总共多少行。
  $all_line++;
  if($sites[$i] == 1){
    print FILE2 $lines[$i];
    $x2++; #保留了多少行。
    $remain_line++; 
  }else{
    $x3++;#过滤掉了多少行。
    $delete_line++;
  }
}
       
my $ratio = $x2/$x1 ;
$ratio = $ratio * 100;   
($x1 == $x2+$x3)   or die;           
print FILE3  "$name1:\n总共$x1行。\n保留了$x2行。\n过滤掉了$x3行。\n保留比例为：$ratio%  \n\n ";

$ave_ratio = $ave_ratio + $ratio;
$all_seq++;

close FILE1;
close FILE2;
}################################################### END FILTER_5_7








##############################################################################################
sub FILTER_5_8   #把某些行过滤掉，考虑20种氨基酸的背景频率，不同类氨基酸应该采用不同阈值，在长度为m的窗口内若有n个小于阈值，则留下。
##############################################################################################
{ 
my $name1 = $_[0];               #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/4_NORMALIZE) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;                 #把类别标签取出来。
my $name2 = $3;                  #把不带路径的文件名取出来。

if (!(-e "$whole/5_FILTER"))          {mkdir "$whole/5_FILTER"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/5_FILTER/$folder"))  {mkdir "$whole/5_FILTER/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";                             #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/5_FILTER/$folder/$name2") or die "Can't open  > $name2 : $!";    #以写的方式打开存放结果的文件。
open(FILE3, ">>", "$whole/5_FILTER/ratio") or die "Can't open  >> ratio : $!";            #存储每一条序列保留了多少行。

my $other = 0.00;
my %hash = (A => 0.0826+$other,
            C => 0.0136+$other,
            D => 0.0545+$other,
            E => 0.0675+$other,
            F => 0.0386+$other,
            G => 0.0708+$other,
            H => 0.0227+$other,
            I => 0.0597+$other,
            K => 0.0584+$other,
            L => 0.0966+$other,
            M => 0.0242+$other,
            N => 0.0406+$other,
            P => 0.0470+$other,
            Q => 0.0393+$other,
            R => 0.0553+$other,
            S => 0.0655+$other,
            T => 0.0534+$other,
            V => 0.0687+$other,
            W => 0.0108+$other,
            Y => 0.0292+$other,
           );
for my $key (keys(%hash)) {
  $hash{$key} = $hash{$key} / 0.01;
  $hash{$key} = 10 * &LOGXY(10,$hash{$key});           #需要把背景频率值转化一下，根据PSSM的构造规则转化。

}


my $window1 = $ARGV[3];
my $window2 = $ARGV[5];
my $threshold = 0;  #阈值随位点的改变而改变。


my $x1 = 0;   #总共多少行。
my $x2 = 0;   #保留了多少行。
my $x3 = 0;   #过滤掉了多少行。

my @lines = <FILE1>;
my @array_max = ();
my @all_sites = ();

for (my $i=0;$i<=$#lines;$i++) {
  if ($lines[$i]=~m/^ {11}/) {
    $i==0 or die;
    $lines[$i] =~ m/^ ([A-Z\s]{20}){21} $/x or die;
    print FILE2 $lines[$i];  #输出第一行
  }else{
    $lines[$i]=~m/ ^[\s0-9]{6}([A-Z])\s{13} ([ 0-9-]{400}) /x or die; 
    $all_sites[$i] = $1;
    my $mid1=$2;
    $all_sites[$i] =~ m/^[A-Z]$/ or die; 
    $mid1 =~ s/ ^\s* ([-0-9][^\n]+[0-9]) \s*$ /$1/x or die;   #去掉首尾的空格。
    my @array1=split(/ +/,$mid1);
    $#array1==19 or die;  #应为20个元素。
    my $max = &MAX(@array1);
    $array_max[$i] = $max;     #各行的最大值都在这里面。$array_max[0]的值应为undef.
  }
}                             
$#lines == $#array_max  or  die;  # $array_max[0]的值应为undef。

my @sites = ();  #保留下来的位点其元素值变为1.
for(my $i=0; $i<=$#array_max; $i++) {$sites[$i] = 0; }

for (my $i1=1; $i1<=$#array_max - $window1 + 1; $i1++) { #每一次循环考虑一个窗口
  given($all_sites[$i1]) {  #阈值随位点的改变而改变。
    when('A') {$threshold = $hash{A}; }
    when('C') {$threshold = $hash{C}; }
    when('D') {$threshold = $hash{D}; }
    when('E') {$threshold = $hash{E}; }
    when('F') {$threshold = $hash{F}; }
    when('G') {$threshold = $hash{G}; }
    when('H') {$threshold = $hash{H}; }
    when('I') {$threshold = $hash{I}; }
    when('K') {$threshold = $hash{K}; }
    when('L') {$threshold = $hash{L}; }
    when('M') {$threshold = $hash{M}; }
    when('N') {$threshold = $hash{N}; }
    when('P') {$threshold = $hash{P}; }
    when('Q') {$threshold = $hash{Q}; }
    when('R') {$threshold = $hash{R}; }
    when('S') {$threshold = $hash{S}; }
    when('T') {$threshold = $hash{T}; }
    when('V') {$threshold = $hash{V}; }
    when('W') {$threshold = $hash{W}; }
    when('Y') {$threshold = $hash{Y}; }
    when($all_sites[$i1] =~ m/[UOBZJX]/) {$threshold = 0; }
    default   {die;        }
  }
  my $p = 0;
  for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) {  #每一次循环考虑一个位点。
    if($array_max[$i2]<$threshold) { $p++; }
  }
  if($p >= $window2) { #若满足条件则保留。
    for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) { 
      $sites[$i2] = 1;
    } 
  }
}

#输出:
for(my $i=1;$i<=$#array_max;$i++) {
  $x1++; #总共多少行。
  $all_line++;
  if($sites[$i] == 1){
    print FILE2 $lines[$i];
    $x2++; #保留了多少行。
    $remain_line++; 
  }else{
    $x3++;#过滤掉了多少行。
    $delete_line++;
  }
}
       
my $ratio = $x2/$x1 ;
$ratio = $ratio * 100;   
($x1 == $x2+$x3)   or die;           
print FILE3  "$name1:\n总共$x1行。\n保留了$x2行。\n过滤掉了$x3行。\n保留比例为：$ratio%  \n\n ";

$ave_ratio = $ave_ratio + $ratio;
$all_seq++;

close FILE1;
close FILE2;
}################################################### END FILTER_5_8







##############################################################################################
sub FILTER_5_9   #不考虑PSSM矩阵的元素，根据序列过滤，在长度为m的窗口内若有n个氨基酸满足条件，则留下。
##############################################################################################
{ 
my $name1 = $_[0];               #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/4_NORMALIZE) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;                 #把类别标签取出来。
my $name2 = $3;                  #把不带路径的文件名取出来。

if (!(-e "$whole/5_FILTER"))          {mkdir "$whole/5_FILTER"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/5_FILTER/$folder"))  {mkdir "$whole/5_FILTER/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";                             #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/5_FILTER/$folder/$name2") or die "Can't open  > $name2 : $!";    #以写的方式打开存放结果的文件。
open(FILE3, ">>", "$whole/5_FILTER/ratio") or die "Can't open  >> ratio : $!";            #存储每一条序列保留了多少行。

my $window1 = $ARGV[3];
my $window2 = $ARGV[5];

my $x1 = 0;   #总共多少行。
my $x2 = 0;   #保留了多少行。
my $x3 = 0;   #过滤掉了多少行。

my @lines = <FILE1>;
my @all_sites = ();

for (my $i=0;$i<=$#lines;$i++) {
  if ($lines[$i]=~m/^ {11}/) {
    $i==0 or die;
    $lines[$i] =~ m/^ ([A-Z\s]{20}){21} $/x or die;
    print FILE2 $lines[$i];  #输出第一行
  }else{
    $lines[$i]=~m/ ^[\s0-9]{6}([A-Z])\s{13} ([ 0-9-]{400}) /x or die; 
    $all_sites[$i] = $1;
    $all_sites[$i] =~ m/^[A-Z]$/ or die; 
  }
}                             
$#lines == $#all_sites  or  die;  # $array_max[0]的值应为undef。

my @sites = ();  #保留下来的位点其元素值变为1.
for(my $i=0; $i<=$#all_sites; $i++) {$sites[$i] = 0; }

for (my $i1=1; $i1<=$#all_sites - $window1 + 1; $i1++) { #每一次循环考虑一个窗口
  my $p = 0;
  for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) {  #每一次循环考虑一个位点。
    if($all_sites[$i2] =~ m/[KR]/) { $p++; }
  }
  if($p >= $window2) { #若满足条件则保留。
    for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) { 
      $sites[$i2] = 1;
    } 
  }
}

#输出:
for(my $i=1; $i<=$#all_sites; $i++) {
  $x1++; #总共多少行。
  $all_line++;
  if($sites[$i] == 1){
    print FILE2 $lines[$i];
    $x2++; #保留了多少行。
    $remain_line++; 
  }else{
    $x3++;#过滤掉了多少行。
    $delete_line++;
  }
}
       
my $ratio = $x2/$x1 ;
$ratio = $ratio * 100;   
($x1 == $x2+$x3)   or die;           
print FILE3  "$name1:\n总共$x1行。\n保留了$x2行。\n过滤掉了$x3行。\n保留比例为：$ratio%  \n\n ";

$ave_ratio = $ave_ratio + $ratio;
$all_seq++;

close FILE1;
close FILE2;
}################################################### END FILTER_5_9








##############################################################################################
sub FILTER_5_10   #不考虑PSSM矩阵的元素，根据序列过滤，在长度为m的窗口内若有n个氨基酸满足条件，则留下。
##############################################################################################
{ 
my $name1 = $_[0];               #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/4_NORMALIZE) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;                 #把类别标签取出来。
my $name2 = $3;                  #把不带路径的文件名取出来。

if (!(-e "$whole/5_FILTER"))          {mkdir "$whole/5_FILTER"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/5_FILTER/$folder"))  {mkdir "$whole/5_FILTER/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";                             #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/5_FILTER/$folder/$name2") or die "Can't open  > $name2 : $!";    #以写的方式打开存放结果的文件。
open(FILE3, ">>", "$whole/5_FILTER/ratio") or die "Can't open  >> ratio : $!";            #存储每一条序列保留了多少行。

my $window1 = $ARGV[3];
my $window2 = $ARGV[5];

my $x1 = 0;   #总共多少行。
my $x2 = 0;   #保留了多少行。
my $x3 = 0;   #过滤掉了多少行。

my @lines = <FILE1>;
my @all_sites = ();

for (my $i=0;$i<=$#lines;$i++) {
  if ($lines[$i]=~m/^ {11}/) {
    $i==0 or die;
    $lines[$i] =~ m/^ ([A-Z\s]{20}){21} $/x or die;
    print FILE2 $lines[$i];  #输出第一行
  }else{
    $lines[$i]=~m/ ^[\s0-9]{6}([A-Z])\s{13} ([ 0-9-]{400}) /x or die; 
    $all_sites[$i] = $1;
    $all_sites[$i] =~ m/^[A-Z]$/ or die; 
  }
}                             
$#lines == $#all_sites  or  die;  # $array_max[0]的值应为undef。

my @sites = ();  #保留下来的位点其元素值变为1.
for(my $i=0; $i<=$#all_sites; $i++) {$sites[$i] = 0; }

for (my $i1=1; $i1<=$#all_sites - $window1 + 1; $i1++) { #每一次循环考虑一个窗口
  my $p = 0;
  for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) {  #每一次循环考虑一个位点。
    if($all_sites[$i2] =~ m/[KR]/) { $p++; }
  }
  if($p < $window2) { #若满足条件则保留。
    for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) { 
      $sites[$i2] = 1;
    } 
  }
}

#输出:
for(my $i=1; $i<=$#all_sites; $i++) {
  $x1++; #总共多少行。
  $all_line++;
  if($sites[$i] == 1){
    print FILE2 $lines[$i];
    $x2++; #保留了多少行。
    $remain_line++; 
  }else{
    $x3++;#过滤掉了多少行。
    $delete_line++;
  }
}
       
my $ratio = $x2/$x1 ;
$ratio = $ratio * 100;   
($x1 == $x2+$x3)   or die;           
print FILE3  "$name1:\n总共$x1行。\n保留了$x2行。\n过滤掉了$x3行。\n保留比例为：$ratio%  \n\n ";

$ave_ratio = $ave_ratio + $ratio;
$all_seq++;

close FILE1;
close FILE2;
}################################################### END FILTER_5_10









##############################################################################################
sub FILTER_5_11   #不考虑PSSM矩阵的元素，根据序列过滤，在长度为m的窗口内若其某种AAindex的平均值大于阈值，则留下。
##############################################################################################
{ 
my $name1 = $_[0];               #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/4_NORMALIZE) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;                 #把类别标签取出来。
my $name2 = $3;                  #把不带路径的文件名取出来。

if (!(-e "$whole/5_FILTER"))          {mkdir "$whole/5_FILTER"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/5_FILTER/$folder"))  {mkdir "$whole/5_FILTER/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";                             #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/5_FILTER/$folder/$name2") or die "Can't open  > $name2 : $!";    #以写的方式打开存放结果的文件。
open(FILE3, ">>", "$whole/5_FILTER/ratio") or die "Can't open  >> ratio : $!";            #存储每一条序列保留了多少行。

my $other = 0.00;
my %hash = (A =>  1.10+$other,  #某种AAindex，Modified Kyte-Doolittle hydrophobicity scale (Juretic et al., 1998)
            C =>  2.50+$other,
            D => -3.60+$other,
            E => -3.20+$other,
            F =>  2.80+$other,
            G => -0.64+$other,
            H => -3.20+$other,
            I =>  4.50+$other,
            K => -4.11+$other,
            L =>  3.80+$other,
            M =>  1.90+$other,
            N => -3.50+$other,
            P => -1.90+$other,
            Q => -3.68+$other,
            R => -5.10+$other,
            S => -0.50+$other,
            T => -0.70+$other,
            V =>  4.20+$other,
            W => -0.46+$other,
            Y => -1.30+$other,
           );



my $window1 = $ARGV[3];
my $threshold = $ARGV[7];  #阈值随位点的改变而改变。


my $x1 = 0;   #总共多少行。
my $x2 = 0;   #保留了多少行。
my $x3 = 0;   #过滤掉了多少行。

my @lines = <FILE1>;
print FILE2 $lines[0];
my @all_sites = ();

for (my $i=0;$i<=$#lines;$i++) {
  if ($lines[$i]=~m/^ {11}/) {
    $i==0 or die;
  }else{
    $i>=1 or die;
    $lines[$i]=~m/ ^[\s0-9]{6}([A-Z])\s{13} ([ 0-9-]{400}) /x or die; 
    $all_sites[$i] = $1;
    $all_sites[$i] =~ m/^[A-Z]$/ or die; 
  }
}                             
$#lines == $#all_sites  or  die;  # $array_max[0]的值应为undef。

my @sites = ();  #保留下来的位点其元素值变为1.
for(my $i=0; $i<=$#all_sites; $i++) {$sites[$i] = 0; }

for (my $i1=1; $i1<=$#all_sites - $window1 + 1; $i1++) { #每一次循环考虑一个窗口
  my $p = 0;
  for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) {  #每一次循环考虑一个位点。
    if($all_sites[$i2] !~ m/^[0-9.-]+$/) {
    given($all_sites[$i2]) {  
      when('A') {$all_sites[$i2] = $hash{A}; }
      when('C') {$all_sites[$i2] = $hash{C}; }
      when('D') {$all_sites[$i2] = $hash{D}; }
      when('E') {$all_sites[$i2] = $hash{E}; }
      when('F') {$all_sites[$i2] = $hash{F}; }
      when('G') {$all_sites[$i2] = $hash{G}; }
      when('H') {$all_sites[$i2] = $hash{H}; }
      when('I') {$all_sites[$i2] = $hash{I}; }
      when('K') {$all_sites[$i2] = $hash{K}; }
      when('L') {$all_sites[$i2] = $hash{L}; }
      when('M') {$all_sites[$i2] = $hash{M}; }
      when('N') {$all_sites[$i2] = $hash{N}; }
      when('P') {$all_sites[$i2] = $hash{P}; }
      when('Q') {$all_sites[$i2] = $hash{Q}; }
      when('R') {$all_sites[$i2] = $hash{R}; }
      when('S') {$all_sites[$i2] = $hash{S}; }
      when('T') {$all_sites[$i2] = $hash{T}; }
      when('V') {$all_sites[$i2] = $hash{V}; }
      when('W') {$all_sites[$i2] = $hash{W}; }
      when('Y') {$all_sites[$i2] = $hash{Y}; }
      when($all_sites[$i2] =~ m/[UOBZJX]/) {$all_sites[$i2] = 0; }
      default   {print  "$all_sites[$i2]\n"; die;        }
    }
    }
    $all_sites[$i2] =~ m/^[0-9.-]+$/  or die ;
    $p = $p + $all_sites[$i2]; 
  }
  $p = $p / $window1 ; 
  if($p >= $threshold) { #若满足条件则保留。
    for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) { 
      $sites[$i2] = 1;
    } 
  }
}

#输出:
for(my $i=1;$i<=$#all_sites;$i++) {
  $x1++; #总共多少行。
  $all_line++;
  if($sites[$i] == 1){
    print FILE2 $lines[$i];
    $x2++; #保留了多少行。
    $remain_line++; 
  }else{
    $x3++;#过滤掉了多少行。
    $delete_line++;
  }
}
       
my $ratio = $x2/$x1 ;
$ratio = $ratio * 100;   
($x1 == $x2+$x3)   or die;           
print FILE3  "$name1:\n总共$x1行。\n保留了$x2行。\n过滤掉了$x3行。\n保留比例为：$ratio%  \n\n ";

$ave_ratio = $ave_ratio + $ratio;
$all_seq++;

close FILE1;
close FILE2;
}################################################### END FILTER_5_11







##############################################################################################
sub FILTER_5_12   #不考虑PSSM矩阵的元素，根据序列过滤，在长度为m的窗口内若其某种AAindex的平均值小于阈值，则留下。
##############################################################################################
{ 
my $name1 = $_[0];               #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/4_NORMALIZE) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;                 #把类别标签取出来。
my $name2 = $3;                  #把不带路径的文件名取出来。

if (!(-e "$whole/5_FILTER"))          {mkdir "$whole/5_FILTER"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/5_FILTER/$folder"))  {mkdir "$whole/5_FILTER/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";                             #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/5_FILTER/$folder/$name2") or die "Can't open  > $name2 : $!";    #以写的方式打开存放结果的文件。
open(FILE3, ">>", "$whole/5_FILTER/ratio") or die "Can't open  >> ratio : $!";            #存储每一条序列保留了多少行。

my $other = 0.00;
my %hash = (A =>  1.10+$other,  #某种AAindex，Modified Kyte-Doolittle hydrophobicity scale (Juretic et al., 1998)
            C =>  2.50+$other,
            D => -3.60+$other,
            E => -3.20+$other,
            F =>  2.80+$other,
            G => -0.64+$other,
            H => -3.20+$other,
            I =>  4.50+$other,
            K => -4.11+$other,
            L =>  3.80+$other,
            M =>  1.90+$other,
            N => -3.50+$other,
            P => -1.90+$other,
            Q => -3.68+$other,
            R => -5.10+$other,
            S => -0.50+$other,
            T => -0.70+$other,
            V =>  4.20+$other,
            W => -0.46+$other,
            Y => -1.30+$other,
           );



my $window1 = $ARGV[3];
my $threshold = $ARGV[7];  #阈值随位点的改变而改变。


my $x1 = 0;   #总共多少行。
my $x2 = 0;   #保留了多少行。
my $x3 = 0;   #过滤掉了多少行。

my @lines = <FILE1>;
print FILE2 $lines[0];
my @all_sites = ();

for (my $i=0;$i<=$#lines;$i++) {
  if ($lines[$i]=~m/^ {11}/) {
    $i==0 or die;
  }else{
    $i>=1 or die;
    $lines[$i]=~m/ ^[\s0-9]{6}([A-Z])\s{13} ([ 0-9-]{400}) /x or die; 
    $all_sites[$i] = $1;
    $all_sites[$i] =~ m/^[A-Z]$/ or die; 
  }
}                             
$#lines == $#all_sites  or  die;  # $array_max[0]的值应为undef。

my @sites = ();  #保留下来的位点其元素值变为1.
for(my $i=0; $i<=$#all_sites; $i++) {$sites[$i] = 0; }

for (my $i1=1; $i1<=$#all_sites - $window1 + 1; $i1++) { #每一次循环考虑一个窗口
  my $p = 0;
  for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) {  #每一次循环考虑一个位点。
    if($all_sites[$i2] !~ m/^[0-9.-]+$/) {
    given($all_sites[$i2]) {  
      when('A') {$all_sites[$i2] = $hash{A}; }
      when('C') {$all_sites[$i2] = $hash{C}; }
      when('D') {$all_sites[$i2] = $hash{D}; }
      when('E') {$all_sites[$i2] = $hash{E}; }
      when('F') {$all_sites[$i2] = $hash{F}; }
      when('G') {$all_sites[$i2] = $hash{G}; }
      when('H') {$all_sites[$i2] = $hash{H}; }
      when('I') {$all_sites[$i2] = $hash{I}; }
      when('K') {$all_sites[$i2] = $hash{K}; }
      when('L') {$all_sites[$i2] = $hash{L}; }
      when('M') {$all_sites[$i2] = $hash{M}; }
      when('N') {$all_sites[$i2] = $hash{N}; }
      when('P') {$all_sites[$i2] = $hash{P}; }
      when('Q') {$all_sites[$i2] = $hash{Q}; }
      when('R') {$all_sites[$i2] = $hash{R}; }
      when('S') {$all_sites[$i2] = $hash{S}; }
      when('T') {$all_sites[$i2] = $hash{T}; }
      when('V') {$all_sites[$i2] = $hash{V}; }
      when('W') {$all_sites[$i2] = $hash{W}; }
      when('Y') {$all_sites[$i2] = $hash{Y}; }
      when($all_sites[$i2] =~ m/[UOBZJX]/) {$all_sites[$i2] = 0; }
      default   {print  "$all_sites[$i2]\n"; die;        }
    }
    }
    $all_sites[$i2] =~ m/^[0-9.-]+$/  or die ;
    $p = $p + $all_sites[$i2]; 
  }
  $p = $p / $window1 ; 
  if($p < $threshold) { #若满足条件则保留。
    for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) { 
      $sites[$i2] = 1;
    } 
  }
}

#输出:
for(my $i=1;$i<=$#all_sites;$i++) {
  $x1++; #总共多少行。
  $all_line++;
  if($sites[$i] == 1){
    print FILE2 $lines[$i];
    $x2++; #保留了多少行。
    $remain_line++; 
  }else{
    $x3++;#过滤掉了多少行。
    $delete_line++;
  }
}
       
my $ratio = $x2/$x1 ;
$ratio = $ratio * 100;   
($x1 == $x2+$x3)   or die;           
print FILE3  "$name1:\n总共$x1行。\n保留了$x2行。\n过滤掉了$x3行。\n保留比例为：$ratio%  \n\n ";

$ave_ratio = $ave_ratio + $ratio;
$all_seq++;

close FILE1;
close FILE2;
}################################################### END FILTER_5_12







##############################################################################################
sub FILTER_5_13   #不考虑PSSM矩阵的元素，根据序列过滤，在长度为m的窗口内若有n个氨基酸的某种AAindex大于阈值，则留下。
##############################################################################################
{ 
my $name1 = $_[0];               #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/4_NORMALIZE) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;                 #把类别标签取出来。
my $name2 = $3;                  #把不带路径的文件名取出来。

if (!(-e "$whole/5_FILTER"))          {mkdir "$whole/5_FILTER"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/5_FILTER/$folder"))  {mkdir "$whole/5_FILTER/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";                             #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/5_FILTER/$folder/$name2") or die "Can't open  > $name2 : $!";    #以写的方式打开存放结果的文件。
open(FILE3, ">>", "$whole/5_FILTER/ratio") or die "Can't open  >> ratio : $!";            #存储每一条序列保留了多少行。

my $other = 0.00;
my %hash = (A =>  1.10+$other,  #某种AAindex，Modified Kyte-Doolittle hydrophobicity scale (Juretic et al., 1998)
            C =>  2.50+$other,
            D => -3.60+$other,
            E => -3.20+$other,
            F =>  2.80+$other,
            G => -0.64+$other,
            H => -3.20+$other,
            I =>  4.50+$other,
            K => -4.11+$other,
            L =>  3.80+$other,
            M =>  1.90+$other,
            N => -3.50+$other,
            P => -1.90+$other,
            Q => -3.68+$other,
            R => -5.10+$other,
            S => -0.50+$other,
            T => -0.70+$other,
            V =>  4.20+$other,
            W => -0.46+$other,
            Y => -1.30+$other,
           );



my $other1 = 0.00;
my %hash1 =(A => 0.0+$other1,  #某种AAindex的阈值。
            C => 0.0+$other1,
            D => 0.0+$other1,
            E => 0.0+$other1,
            F => 0.0+$other1,
            G => 0.0+$other1,
            H => 0.0+$other1,
            I => 0.0+$other1,
            K => 0.0+$other1,
            L => 0.0+$other1,
            M => 0.0+$other1,
            N => 0.0+$other1,
            P => 0.0+$other1,
            Q => 0.0+$other1,
            R => 0.0+$other1,
            S => 0.0+$other1,
            T => 0.0+$other1,
            V => 0.0+$other1,
            W => 0.0+$other1,
            Y => 0.0+$other1,
           );

my $window1 = $ARGV[3];
my $window2 = $ARGV[5];
my $threshold = $ARGV[7];  #阈值随位点的改变而改变。


my $x1 = 0;   #总共多少行。
my $x2 = 0;   #保留了多少行。
my $x3 = 0;   #过滤掉了多少行。

my @lines = <FILE1>;
print FILE2 $lines[0];
my @all_sites = ();

for (my $i=0;$i<=$#lines;$i++) {
  if ($lines[$i]=~m/^ {11}/) {
    $i==0 or die;
  }else{
    $i>=1 or die;
    $lines[$i]=~m/ ^[\s0-9]{6}([A-Z])\s{13} ([ 0-9-]{400}) /x or die; 
    $all_sites[$i] = $1;
    $all_sites[$i] =~ m/^[A-Z]$/ or die; 
  }
}                             
$#lines == $#all_sites  or  die;  # $array_max[0]的值应为undef。

my @sites = ();  #保留下来的位点其元素值变为1.
for(my $i=0; $i<=$#all_sites; $i++) {$sites[$i] = 0; }

for (my $i1=1; $i1<=$#all_sites - $window1 + 1; $i1++) { #每一次循环考虑一个窗口,步长为1.
  my $p = 0;
  given($all_sites[$i1]) {  #阈值随位点而改变。
    when('A') {$threshold = $hash1{A}; }
    when('C') {$threshold = $hash1{C}; }
    when('D') {$threshold = $hash1{D}; }
    when('E') {$threshold = $hash1{E}; }
    when('F') {$threshold = $hash1{F}; }
    when('G') {$threshold = $hash1{G}; }
    when('H') {$threshold = $hash1{H}; }
    when('I') {$threshold = $hash1{I}; }
    when('K') {$threshold = $hash1{K}; }
    when('L') {$threshold = $hash1{L}; }
    when('M') {$threshold = $hash1{M}; }
    when('N') {$threshold = $hash1{N}; }
    when('P') {$threshold = $hash1{P}; }
    when('Q') {$threshold = $hash1{Q}; }
    when('R') {$threshold = $hash1{R}; }
    when('S') {$threshold = $hash1{S}; }
    when('T') {$threshold = $hash1{T}; }
    when('V') {$threshold = $hash1{V}; }
    when('W') {$threshold = $hash1{W}; }
    when('Y') {$threshold = $hash1{Y}; }
    when($all_sites[$i1] =~ m/[UOBZJX]/) {$threshold = 0; }
    default   {print  "$all_sites[$i1]\n"; die;        }
  } 
  for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) {  #把氨基酸序列数值化。
    my $mid = $all_sites[$i2];
    if($mid !~ m/^[0-9.-]+$/) {
    given($mid) {  
      when('A') {$mid = $hash{A}; }
      when('C') {$mid = $hash{C}; }
      when('D') {$mid = $hash{D}; }
      when('E') {$mid = $hash{E}; }
      when('F') {$mid = $hash{F}; }
      when('G') {$mid = $hash{G}; }
      when('H') {$mid = $hash{H}; }
      when('I') {$mid = $hash{I}; }
      when('K') {$mid = $hash{K}; }
      when('L') {$mid = $hash{L}; }
      when('M') {$mid = $hash{M}; }
      when('N') {$mid = $hash{N}; }
      when('P') {$mid = $hash{P}; }
      when('Q') {$mid = $hash{Q}; }
      when('R') {$mid = $hash{R}; }
      when('S') {$mid = $hash{S}; }
      when('T') {$mid = $hash{T}; }
      when('V') {$mid = $hash{V}; }
      when('W') {$mid = $hash{W}; }
      when('Y') {$mid = $hash{Y}; }
      when($mid =~ m/[UOBZJX]/) {$mid = 0; }
      default   {print  "$mid \n"; die;        }
    }
    }
    $mid =~ m/^[0-9.-]+$/  or die ;
    if($mid >= $threshold) { $p++; } 
  }
  if($p >= $window2) { #若满足条件则保留。
    for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) { 
      $sites[$i2] = 1;
    } 
  }
}

#输出:
for(my $i=1;$i<=$#all_sites;$i++) {
  $x1++; #总共多少行。
  $all_line++;
  if($sites[$i] == 1){
    print FILE2 $lines[$i];
    $x2++; #保留了多少行。
    $remain_line++; 
  }else{
    $x3++;#过滤掉了多少行。
    $delete_line++;
  }
}
       
my $ratio = $x2/$x1 ;
$ratio = $ratio * 100;   
($x1 == $x2+$x3)   or die;           
print FILE3  "$name1:\n总共$x1行。\n保留了$x2行。\n过滤掉了$x3行。\n保留比例为：$ratio%  \n\n ";

$ave_ratio = $ave_ratio + $ratio;
$all_seq++;

close FILE1;
close FILE2;
}################################################### END FILTER_5_13







##############################################################################################
sub FILTER_5_14   #不考虑PSSM矩阵的元素，根据序列过滤，在长度为m的窗口内若有n个氨基酸的某种AAindex小于阈值，则留下。
##############################################################################################
{ 
my $name1 = $_[0];               #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/4_NORMALIZE) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;                 #把类别标签取出来。
my $name2 = $3;                  #把不带路径的文件名取出来。

if (!(-e "$whole/5_FILTER"))          {mkdir "$whole/5_FILTER"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (!(-e "$whole/5_FILTER/$folder"))  {mkdir "$whole/5_FILTER/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";                             #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/5_FILTER/$folder/$name2") or die "Can't open  > $name2 : $!";    #以写的方式打开存放结果的文件。
open(FILE3, ">>", "$whole/5_FILTER/ratio") or die "Can't open  >> ratio : $!";            #存储每一条序列保留了多少行。

my $other = 0.00;
my %hash = (A =>  1.10+$other,  #某种AAindex，Modified Kyte-Doolittle hydrophobicity scale (Juretic et al., 1998)
            C =>  2.50+$other,
            D => -3.60+$other,
            E => -3.20+$other,
            F =>  2.80+$other,
            G => -0.64+$other,
            H => -3.20+$other,
            I =>  4.50+$other,
            K => -4.11+$other,
            L =>  3.80+$other,
            M =>  1.90+$other,
            N => -3.50+$other,
            P => -1.90+$other,
            Q => -3.68+$other,
            R => -5.10+$other,
            S => -0.50+$other,
            T => -0.70+$other,
            V =>  4.20+$other,
            W => -0.46+$other,
            Y => -1.30+$other,
           );



my $other1 = 0.00;
my %hash1 =(A => 0.0+$other1,  #某种AAindex的阈值。
            C => 0.0+$other1,
            D => 0.0+$other1,
            E => 0.0+$other1,
            F => 0.0+$other1,
            G => 0.0+$other1,
            H => 0.0+$other1,
            I => 0.0+$other1,
            K => 0.0+$other1,
            L => 0.0+$other1,
            M => 0.0+$other1,
            N => 0.0+$other1,
            P => 0.0+$other1,
            Q => 0.0+$other1,
            R => 0.0+$other1,
            S => 0.0+$other1,
            T => 0.0+$other1,
            V => 0.0+$other1,
            W => 0.0+$other1,
            Y => 0.0+$other1,
           );

my $window1 = $ARGV[3];
my $window2 = $ARGV[5];
my $threshold = $ARGV[7];  #阈值随位点的改变而改变。


my $x1 = 0;   #总共多少行。
my $x2 = 0;   #保留了多少行。
my $x3 = 0;   #过滤掉了多少行。

my @lines = <FILE1>;
print FILE2 $lines[0];
my @all_sites = ();

for (my $i=0;$i<=$#lines;$i++) {
  if ($lines[$i]=~m/^ {11}/) {
    $i==0 or die;
  }else{
    $i>=1 or die;
    $lines[$i]=~m/ ^[\s0-9]{6}([A-Z])\s{13} ([ 0-9-]{400}) /x or die; 
    $all_sites[$i] = $1;
    $all_sites[$i] =~ m/^[A-Z]$/ or die; 
  }
}                             
$#lines == $#all_sites  or  die;  # $array_max[0]的值应为undef。

my @sites = ();  #保留下来的位点其元素值变为1.
for(my $i=0; $i<=$#all_sites; $i++) {$sites[$i] = 0; }

for (my $i1=1; $i1<=$#all_sites - $window1 + 1; $i1++) { #每一次循环考虑一个窗口,步长为1.
  my $p = 0;
  given($all_sites[$i1]) {  #阈值随位点而改变。
    when('A') {$threshold = $hash1{A}; }
    when('C') {$threshold = $hash1{C}; }
    when('D') {$threshold = $hash1{D}; }
    when('E') {$threshold = $hash1{E}; }
    when('F') {$threshold = $hash1{F}; }
    when('G') {$threshold = $hash1{G}; }
    when('H') {$threshold = $hash1{H}; }
    when('I') {$threshold = $hash1{I}; }
    when('K') {$threshold = $hash1{K}; }
    when('L') {$threshold = $hash1{L}; }
    when('M') {$threshold = $hash1{M}; }
    when('N') {$threshold = $hash1{N}; }
    when('P') {$threshold = $hash1{P}; }
    when('Q') {$threshold = $hash1{Q}; }
    when('R') {$threshold = $hash1{R}; }
    when('S') {$threshold = $hash1{S}; }
    when('T') {$threshold = $hash1{T}; }
    when('V') {$threshold = $hash1{V}; }
    when('W') {$threshold = $hash1{W}; }
    when('Y') {$threshold = $hash1{Y}; }
    when($all_sites[$i1] =~ m/[UOBZJX]/) {$threshold = 0; }
    default   {print  "$all_sites[$i1]\n"; die;        }
  } 
  for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) {  #把氨基酸序列数值化。
    my $mid = $all_sites[$i2];
    if($mid !~ m/^[0-9.-]+$/) {
    given($mid) {  
      when('A') {$mid = $hash{A}; }
      when('C') {$mid = $hash{C}; }
      when('D') {$mid = $hash{D}; }
      when('E') {$mid = $hash{E}; }
      when('F') {$mid = $hash{F}; }
      when('G') {$mid = $hash{G}; }
      when('H') {$mid = $hash{H}; }
      when('I') {$mid = $hash{I}; }
      when('K') {$mid = $hash{K}; }
      when('L') {$mid = $hash{L}; }
      when('M') {$mid = $hash{M}; }
      when('N') {$mid = $hash{N}; }
      when('P') {$mid = $hash{P}; }
      when('Q') {$mid = $hash{Q}; }
      when('R') {$mid = $hash{R}; }
      when('S') {$mid = $hash{S}; }
      when('T') {$mid = $hash{T}; }
      when('V') {$mid = $hash{V}; }
      when('W') {$mid = $hash{W}; }
      when('Y') {$mid = $hash{Y}; }
      when($mid =~ m/[UOBZJX]/) {$mid = 0; }
      default   {print  "$mid \n"; die;        }
    }
    }
    $mid =~ m/^[0-9.-]+$/  or die ;
    if($mid < $threshold) { $p++; } 
  }
  if($p >= $window2) { #若满足条件则保留。
    for (my $i2=$i1; $i2<=$i1 + $window1 - 1; $i2++) { 
      $sites[$i2] = 1;
    } 
  }
}

#输出:
for(my $i=1;$i<=$#all_sites;$i++) {
  $x1++; #总共多少行。
  $all_line++;
  if($sites[$i] == 1){
    print FILE2 $lines[$i];
    $x2++; #保留了多少行。
    $remain_line++; 
  }else{
    $x3++;#过滤掉了多少行。
    $delete_line++;
  }
}
       
my $ratio = $x2/$x1 ;
$ratio = $ratio * 100;   
($x1 == $x2+$x3)   or die;           
print FILE3  "$name1:\n总共$x1行。\n保留了$x2行。\n过滤掉了$x3行。\n保留比例为：$ratio%  \n\n ";

$ave_ratio = $ave_ratio + $ratio;
$all_seq++;

close FILE1;
close FILE2;
}################################################### END FILTER_5_14












##############################################################################################
sub VECTOR_6_1    #各行相加，成为20D的向量。
##############################################################################################
{ 
my $name1 = $_[0];   #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/5_FILTER) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;    #把类别标签取出来。
my $name2 = $3;     #把不带路径的文件名取出来。

if (-e "$whole/6_VECTOR")         { }else{mkdir "$whole/6_VECTOR"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (-e "$whole/6_VECTOR/$folder") { }else{mkdir "$whole/6_VECTOR/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";  #以读的方式打开文件。
open(FILE2, ">", "$whole/6_VECTOR/$folder/$name2") or die "Can't open  > $name2 : $!";   

my @lines = <FILE1>;
my @sub = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);

$lines[0]=~m/(^ {20})( {19}A {19}R {19}N {19}D {19}C {19}Q {19}E {19}G {19}H {19}I {19}L {19}K {19}M {19}F {19}P {19}S {19}T {19}W {19}Y {19}V)/ or die;
#第一行不输出。

my $remain = 0;  #记录保留了多少行。
for (my $i=1;$i<=$#lines;$i++) { 
    $remain++;
    $lines[$i]=~m/(^[ 0-9]{5} [A-Z*] {13})([ .0-9-]{400})/ or die; 
    my $mid1=$2;
    $mid1 =~ s/^\s*([-0-9][^\n]+[0-9A-Z])\s*$/$1/ or die;   #去掉首尾的空格。
    my @array1=split(/ +/,$mid1);
    $#array1==19 or die;
    for (my $j=0;$j<=19;$j++) {
      $sub[$j] = $sub[$j] + $array1[$j];
    }
}

$remain==$#lines   or die;
my @average = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
for (my $j=0;$j<=19;$j++) {
  $average[$j] = $sub[$j]/($remain+10**(-20));
}

#输出结果:
for (my $j=0;$j<=19;$j++) {
  printf FILE2 "%20.10f",  $average[$j];
}                           

close FILE1;
close FILE2;
} #################################### END VECTOR_6_1






##############################################################################################
sub VECTOR_6_2     #同种氨基酸的行相加，成为400D的向量。
##############################################################################################
{ 
my $name1 = $_[0];   #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/5_FILTER) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;    #把类别标签取出来。
my $name2 = $3;     #把不带路径的文件名取出来。

if (-e "$whole/6_VECTOR")         { }else{mkdir "$whole/6_VECTOR"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (-e "$whole/6_VECTOR/$folder") { }else{mkdir "$whole/6_VECTOR/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";  #以读的方式打开文件。
open(FILE2, ">", "$whole/6_VECTOR/$folder/$name2") or die "Can't open  > $name2 : $!";   

my @lines = <FILE1>;

$lines[0]=~m/(^ {20})( {19}A {19}R {19}N {19}D {19}C {19}Q {19}E {19}G {19}H {19}I {19}L {19}K {19}M {19}F {19}P {19}S {19}T {19}W {19}Y {19}V)/ or die;
print FILE2 $lines[0];  #输出第一行
    
my @array2d = ();  #声明一个数组。
for (my $i1=1;$i1<=$#lines;$i1++) {  #数组的第一个元素不做处理。
  $lines[$i1] =~ s/^\s*([-0-9][^\n]+[0-9A-Z])\s*$/$1/ or die;     #把每个元素(不含第一个，下同)开头和结尾的空格去掉。
  $array2d[$i1] = [split(/ +/,$lines[$i1])];  ##用方括号生成匿名引用。
  $#{$array2d[$i1]} == 21 or $#{$array2d[$i1]} == 1 or die; #看成一维数组时，其每个元素均为22维或2维的数组。
}	
#现在@array2d是一个L*22维的2维数组。L是序列长度。

my $protein = '';
for (my $i=1;$i <= $#lines; $i++ ) {
	$protein = $protein."$array2d[$i][1]";   #把蛋白质序列放在这个变量里。
}
(length($protein) == $#lines)  or die;

my @table = qw( A R N D C Q E G H I L K M F P S T W Y V);	
my @temp1 = ();
my @average = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);

for my $letter (@table) {
  if ($protein =~ m/$letter/) {  #蛋白质序列含有这种氨基酸时才执行下面的语句。
    @temp1 = ();
    for (my $i=1;$i<=$#lines;$i++) {
      if ($array2d[$i][1] eq $letter) {          #把不同位置的同种氨基酸提出来。
        $temp1[++$#temp1] = $i;                  #把行号（位置号）存起来。         
      }
    }
    my @sub = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
    foreach my $k (@temp1){
      for (my $j=2;$j<=21;$j++) { 
  	$sub[$j] = $sub[$j] + $array2d[$k][$j];
      }
    }
    @average = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
    $#temp1 >= 0 or die;
    for (my $j=2;$j<=21;$j++) {
      $average[$j] = $sub[$j]/($#temp1 + 1);
    }
  }else{
    for (my $j=2;$j<=21;$j++) { 
      $average[$j] = 0;
    }
  }
  #输出：
  printf FILE2 "%20s",$letter;
  for (my $j=2;$j<=21;$j++) {
    printf FILE2 "%20.10f",  $average[$j];
  }	
  print FILE2 "\n";	 
}

close FILE1;
close FILE2;
}############################################### END VECTOR_6_2






##############################################################################################
sub VECTOR_6_3    #把序列平均分成N段，每一段的各行相加，成为20*N D的向量。
##############################################################################################
{ 
my $name1 = $_[0];   #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/5_FILTER) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;    #把类别标签取出来。
my $name2 = $3;     #把不带路径的文件名取出来。

if (-e "$whole/6_VECTOR")         { }else{mkdir "$whole/6_VECTOR"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (-e "$whole/6_VECTOR/$folder") { }else{mkdir "$whole/6_VECTOR/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";  #以读的方式打开文件。
open(FILE2, ">", "$whole/6_VECTOR/$folder/$name2") or die "Can't open  > $name2 : $!";   

my @lines = <FILE1>;


$lines[0]=~m/(^ {20})( {19}A {19}R {19}N {19}D {19}C {19}Q {19}E {19}G {19}H {19}I {19}L {19}K {19}M {19}F {19}P {19}S {19}T {19}W {19}Y {19}V)/ or die;
#第一行不输出。

my $leng = $#lines/$ARGV[11];            #把序列分成N段，每一段的长度。
my $l=0;
if($leng =~ m/^[0-9]+$/){  
  $l=$leng;   
}else{
  $leng =~ m/^([0-9]+)\.[0-9]+$/ or die;  #小数一律按去尾法处理。
  $l = $1;
}
$l =~ m/^[0-9]+$/ or die;
my $remain = 0;

for (my $m=1; $m<=$ARGV[11] - 1; $m++) {   #一段一段地处理，前N-1段，最后一段特殊处理。
  my @sub = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
  my @average = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
  for (my $i=1+($m-1)*$l; $i<=1+$m*$l-1; $i++) { 
    $remain++;
    $lines[$i]=~m/(^[ 0-9]{5} [A-Z] {13})([ .0-9-]{400})/ or die; 
    my $mid1=$2;
    $mid1 =~ s/^\s*([-0-9][^\n]+[0-9A-Z])\s*$/$1/ or die;   #去掉首尾的空格。
    my @array1=split(/ +/,$mid1);
    $#array1==19 or die;
    for (my $j=0;$j<=19;$j++) {
      $sub[$j] = $sub[$j] + $array1[$j];
    }
  }
  for (my $j=0;$j<=19;$j++) {
    $average[$j] = $sub[$j]/($l+10**(-20));
  }
  #输出结果:
  for (my $j=0;$j<=$#average;$j++) {
    printf  FILE2  "%20.10f",  $average[$j];
  }  
}
 
{
  my @sub = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
  my @average = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
  for (my $i=1+($ARGV[11]-1)*$l; $i<=$#lines; $i++) { 
    $remain++;
    $lines[$i]=~m/(^[ 0-9]{5} [A-Z] {13})([ .0-9-]{400})/ or die; 
    my $mid1=$2;
    $mid1 =~ s/^\s*([-0-9][^\n]+[0-9A-Z])\s*$/$1/ or die;   #去掉首尾的空格。
    my @array1=split(/ +/,$mid1);
    $#array1==19 or die;
    for (my $j=0;$j<=19;$j++) {
      $sub[$j] = $sub[$j] + $array1[$j];
    }
  }
  for (my $j=0;$j<=19;$j++) {
    $average[$j] = $sub[$j]/($l+10**(-20));
  }
  #输出结果:
  for (my $j=0;$j<=$#average;$j++) {
    printf  FILE2  "%20.10f",  $average[$j];
  } 
}

$remain==$#lines  or die;
close FILE1;
close FILE2;
}#################################### END VECTOR_6_3






##############################################################################################
sub VECTOR_6_4    #求序列长度，成为1D的向量。
##############################################################################################
{ 
my $name1 = $_[0];   #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/5_FILTER) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;    #把类别标签取出来。
my $name2 = $3;     #把不带路径的文件名取出来。

if (-e "$whole/6_VECTOR")         { }else{mkdir "$whole/6_VECTOR"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (-e "$whole/6_VECTOR/$folder") { }else{mkdir "$whole/6_VECTOR/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";  #以读的方式打开文件。
open(FILE2, ">", "$whole/6_VECTOR/$folder/$name2") or die "Can't open  > $name2 : $!";   

my @lines = <FILE1>;
my @sub = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);

$lines[0]=~m/(^ {20})( {19}A {19}R {19}N {19}D {19}C {19}Q {19}E {19}G {19}H {19}I {19}L {19}K {19}M {19}F {19}P {19}S {19}T {19}W {19}Y {19}V)/ or die;
#第一行不输出。

my $sequence = &SEQUENCE($_[0]);
my $seq_length = length($sequence);

#输出结果:
printf FILE2 "%20.10f",  $seq_length;
                         

close FILE1;
close FILE2;
} #################################### END VECTOR_6_4








##############################################################################################
sub VECTOR_6_5    #分段氨基酸组成.SAAC.
##############################################################################################
{ 
my $name1 = $_[0];   #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/5_FILTER) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;    #把类别标签取出来。
my $name2 = $3;     #把不带路径的文件名取出来。

if (-e "$whole/6_VECTOR")         { }else{mkdir "$whole/6_VECTOR"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (-e "$whole/6_VECTOR/$folder") { }else{mkdir "$whole/6_VECTOR/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";  #以读的方式打开文件。
open(FILE2, ">", "$whole/6_VECTOR/$folder/$name2") or die "Can't open  > $name2 : $!";   

my @lines = <FILE1>;
my $sequence = &SEQUENCE($_[0]); 
$sequence =~ m/^[A-Z]+$/ or die "$sequence"; #首先保证序列不为空。
$sequence =~ s/[UOBZJX]//g ;  #去掉非标准氨基酸。
$sequence !~ m/[UOBZJX]/ or die "$sequence"; 
my $remain = length($sequence); 


#把序列平均分成N段：
my $leng = $#lines/$ARGV[13];            #把序列分成N段，每一段的长度。
my $l=0;
if($leng =~ m/^[0-9]+$/){  
  $l=$leng;   
}else{
  $leng =~ m/^([0-9]+)\.[0-9]+$/ or die;  #小数一律按去尾法处理。
  $l = $1;
}
$l =~ m/^[0-9]+$/ or die;

my @seq = ();
for (my $m=1; $m<=$ARGV[13] - 1; $m++) {   #一段一段地处理，前N-1段，最后一段特殊处理。
  my $i=0+($m-1)*$l; 
  $seq[$m] = substr($sequence,$i,$l);
}

my $last = $#lines - ($ARGV[13]-1)*$l +1;
$seq[$ARGV[13]] = substr($sequence, ($ARGV[13]-1)*$l, $last); #最后一段特殊处理。

my $all_aa=0;
for (my $m=1; $m<=$ARGV[13]; $m++) { 
  $all_aa = $all_aa + length($seq[$m]);
}
($all_aa== $#lines)  or die;  #各片段长度之和应该等于序列长度。


#一段一段地处理:
for (my $m=1; $m<=$ARGV[13]; $m++) { 
  #全序列的k肽数组：
  my @all=();  #其每一个元素表示一个k肽。
  my $frag = length($seq[$m]) ;                        
  for (my $i=0; $i<=$frag-$k_m; $i++) {    #若线性序列的长度为N，则k肽的数目应该为N-(k-1).
    $all[++$#all]=substr($seq[$m],$i,$k_m);
  }
     
  #频数数组:
  my @number=();      
  for (my $i1=0;$i1<=$#mer;$i1++) {      #$i1是循环变量。
    my $p1=0;                              
    for (my $i2=0;$i2<=$#all;$i2++){     #注意区分K-mer有多少种与K-mer有多少个.
      if ($all[$i2]=~m/$mer[$i1]/) {     #统计频数
        $number[$i1]++;                  #数组@number的第i个元素表示@mer的第i个元素出现的频数。
        $p1=1;
      }
    }
    if ($p1==0) { $number[$i1]=0;}
  }  

  #求频数的和:
  my $add1=0;  
  foreach my $num(@number) {    
    $add1=$add1+$num;
  }
  $add1==($frag-$k_m+1 ) or die;

  #频率数组:
  my @frequency=();
  for (my $j=0;$j<=$#number;$j++) {               
    $frequency[$j]=$number[$j]/($add1+10**(-20));
  }

  ($#mer==$kinds**$k_m -1 and  $#number==$kinds**$k_m -1 and $#frequency==$kinds**$k_m -1) or die ;

  #输出结果:
  for (my $j=0;$j<=$kinds**$k_m -1;$j++) {
    printf FILE2 "%20.10f",  $frequency[$j];
  }                           
}

close FILE1;
close FILE2;
} #################################### END VECTOR_6_5 









##############################################################################################
sub CONVERT_7_1    #对向量的进一步处理，适用于VECTOR_5_1，VECTOR_5_3的结果。
##############################################################################################
{ 
my $name1 = $_[0];   #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/6_VECTOR) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;    #把类别标签取出来。
my $name2 = $3;     #把不带路径的文件名取出来。

if (-e "$whole/7_CONVERT")         { }else{mkdir "$whole/7_CONVERT"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (-e "$whole/7_CONVERT/$folder") { }else{mkdir "$whole/7_CONVERT/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";  #以读的方式打开文件。
open(FILE2, ">", "$whole/7_CONVERT/$folder/$name2") or die "Can't open  > $name2 : $!";   

my @lines = <FILE1>;
$#lines == 0 or die;
print FILE2   $lines[0];  

close FILE1;
close FILE2;
} ############################################ END CONVERT_7_1






##############################################################################################
sub CONVERT_7_2     #对向量的进一步处理，适用于VECTOR_5_2的结果。
##############################################################################################
{ 
my $name1 = $_[0];   #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/6_VECTOR) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;    #把类别标签取出来。
my $name2 = $3;     #把不带路径的文件名取出来。

if (-e "$whole/7_CONVERT")         { }else{mkdir "$whole/7_CONVERT"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。
if (-e "$whole/7_CONVERT/$folder") { }else{mkdir "$whole/7_CONVERT/$folder"   or die;}      #若此目录不存在，则创建。用于存放这个子程序的相应类别的输出结果。                                           

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";  #以读的方式打开文件。
open(FILE2, ">", "$whole/7_CONVERT/$folder/$name2") or die "Can't open  > $name2 : $!";   

my @lines = <FILE1>;
$#lines == 20 or die;
$lines[0]=~m/(^ {20})( {19}A {19}R {19}N {19}D {19}C {19}Q {19}E {19}G {19}H {19}I {19}L {19}K {19}M {19}F {19}P {19}S {19}T {19}W {19}Y {19}V)/ or die;

for (my $i=1;$i<=$#lines;$i++) {
  $lines[$i] =~ m/^\s{19}[A-Z]([- .0-9]{400})\n$/ or die;  
  printf FILE2 $1;
}

close FILE1;
close FILE2;
} ################################# END CONVERT_7_2







##############################################################################################
sub MERGE_8_1   #合并: 把各类的向量分别放在不同的文件里，同一类的放在同一个文件。还要把所有向量放在一个文件里。
##############################################################################################
{ 
my $name1 = $_[0];   #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/7_CONVERT) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;    #把类别标签取出来。

if (-e "$whole/8_1_MERGE")         { }else{mkdir "$whole/8_1_MERGE"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";  #以读的方式打开文件。
open(FILE2, ">>", "$whole/8_1_MERGE/$folder") or die "Can't open >> $folder: $!";      
open(FILE3, ">>", "$whole/8_1_MERGE/subnuc1") or die "Can't open  >> subnuc : $!";   

my @lines = <FILE1>;    
$#lines == 0 or die;
$lines[0]=~m/^([0-9 .-]{20})+$/ or die;
$folder =~ s/^([0-9]+)-([0-9]+)/$1/  ;
printf FILE2 "+%-8s $lines[0]\n",$folder;
if($folder !~ m/-/) {printf FILE3 "+%-8s $lines[0]\n",$folder;}  #多定位的不用合并到这个文件。

close FILE1;
close FILE2;
close FILE3;
}################################# END MERGE_8_1 










##############################################################################################
sub MERGE_8_2   #合并，并把ID号加在前面。
##############################################################################################
{ 
my $name1 = $_[0];   #接受参数，待处理文件的名字(含路径)。
$name1 =~ m/ (^$whole\/7_CONVERT) \/ ([0-9-]+) \/ ([^\s]+)$ /x or die;  
my $folder = $2;    #把类别标签取出来。
my $id = $3;

if (-e "$whole/8_2_MERGE")         { }else{mkdir "$whole/8_2_MERGE"   or die;}              #若此目录不存在，则创建。用于存放这个子程序的所有输出结果。

open(FILE1, "<", "$name1") or die "Can't open < $name1 : $!";  #以读的方式打开文件。
open(FILE2, ">>", "$whole/8_2_MERGE/$folder") or die "Can't open  > $folder : $!";      
open(FILE3, ">>", "$whole/8_2_MERGE/subnuc2") or die "Can't open  >> subnuc : $!";   

my @lines = <FILE1>;    
$#lines == 0 or die;
$lines[0]=~m/^([0-9 .-]{20})+$/ or die;
$folder =~ s/^([0-9]+)-([0-9]+)/$1/  ;
printf FILE2 "%-19s +%-8s $lines[0]\n",$id,$folder;
if($folder !~ m/-/) {printf FILE3 "%-19s +%-8s $lines[0]\n",$id,$folder;}  #多定位的不用合并到这个文件。

close FILE1;
close FILE2;
close FILE3;
}################################# END MERGE_8_2 







##############################################################################################
sub NORMALIZE_9_0A     #原样输出，不对数值进行缩放。(第二次缩放)
##############################################################################################
{ 
my $name1 = $_[0];              #接受参数，待处理文件的名字(含路径)。
$name1 =~ m|/([^/]+)$| or die;   
my $name2 = $1;                 #输出文件的名字。

if (-e "$whole/9_1_NORMALIZE")         { }else{mkdir "$whole/9_1_NORMALIZE"   or die;}  
                                       
open(FILE1, "<", "$name1") or die;                        #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/9_1_NORMALIZE/$name2") or die;   #以写的方式打开存放结果的文件。

while (my $line=<FILE1>) {
  $line =~ m/^\+/ or die;
  print FILE2 $line;
}                              

close FILE1;
close FILE2;
} ############################################# END NORMALIZE_9_0A 





##############################################################################################
sub NORMALIZE_9_0B     #原样输出，不对数值进行缩放。(第二次缩放)
##############################################################################################
{ 
my $name1 = $_[0];              #接受参数，待处理文件的名字(含路径)。
$name1 =~ m|/([^/]+)$| or die;   
my $name2 = $1;                 #输出文件的名字。

if (-e "$whole/9_2_NORMALIZE")  { }else{mkdir "$whole/9_2_NORMALIZE"   or die;}  
                                       
open(FILE1, "<", "$name1") or die;                        #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/9_2_NORMALIZE/$name2") or die;   #以写的方式打开存放结果的文件。

while (my $line=<FILE1>) {
  $line =~ m/^[A-Z_0-9]+/ or die;
  print FILE2 $line;
}                              

close FILE1;
close FILE2;
} ############################################# END NORMALIZE_9_0B 










##############################################################################################
sub NORMALIZE_9_1A     #缩放，1/(1 + e^(-x))    (第二次缩放)
##############################################################################################
{ 
my $name1 = $_[0];              #接受参数，待处理文件的名字(含路径)。
$name1 =~ m|/([^/]+)$| or die;   
my $name2 = $1;                 #输出文件的名字。

if (-e "$whole/9_1_NORMALIZE")         { }else{mkdir "$whole/9_1_NORMALIZE"   or die;}                                         
open(FILE1, "<", "$name1") or die;                      #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/9_1_NORMALIZE/$name2") or die;   #以写的方式打开存放结果的文件。

while (my $line=<FILE1>) {
    $line =~ s/ ^\s* ([^\s][^\n]+[^\s]) \s*$ /$1/x or die;   #去掉首尾的空格。
    my @array=split(/\s+/,$line);
    for(my $i=0;$i<=$#array;$i++){
      if($array[$i] =~ m/^[0-9.-]+$/) {
        $i >= 1  or die;
        $array[$i]=1/(1+exp(-$array[$i]));  #归一化。
        printf FILE2 "%20.10f",$array[$i];
      }else{
        $i==0  or die;
        printf FILE2 "%-20s",$array[$i];
      }
    }
    print FILE2 "\n";
}                              

close FILE1;
close FILE2;
} ############################################# END NORMALIZE_9_1A 






##############################################################################################
sub NORMALIZE_9_1B     #缩放，1/(1 + e^(-x))    (第二次缩放)
##############################################################################################
{ 
my $name1 = $_[0];              #接受参数，待处理文件的名字(含路径)。
$name1 =~ m|/([^/]+)$| or die;   
my $name2 = $1;                 #输出文件的名字。

if (-e "$whole/9_2_NORMALIZE")         { }else{mkdir "$whole/9_2_NORMALIZE"   or die;}                                         
open(FILE1, "<", "$name1") or die;                      #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/9_2_NORMALIZE/$name2") or die;   #以写的方式打开存放结果的文件。

while (my $line=<FILE1>) {
    $line =~ s/ ^\s* ([^\s][^\n]+[^\s]) \s*$ /$1/x or die;   #去掉首尾的空格。
    my @array=split(/\s+/,$line);
    for(my $i=0;$i<=$#array;$i++){
      if($array[$i] =~ m/^[0-9.-]+$/) {
        $i >= 2  or die;
        $array[$i]=1/(1+exp(-$array[$i]));  #归一化。
        printf FILE2 "%20.10f",$array[$i];
      }else{
        $i==0  or $i==1 or die;
        printf FILE2 "%-20s",$array[$i];
      }
    }
    print FILE2 "\n";
}                              

close FILE1;
close FILE2;
} ############################################# END NORMALIZE_9_1B 






##############################################################################################
sub NORMALIZE_9_2A    #缩放，(x - min)/(max - min),按行。
##############################################################################################
{ 
my $name1 = $_[0];              #接受参数，待处理文件的名字(含路径)。
$name1 =~ m|/([^/]+)$| or die;   
my $name2 = $1;                 #输出文件的名字。

if (-e "$whole/9_1_NORMALIZE")         { }else{mkdir "$whole/9_1_NORMALIZE"   or die;}                                         
open(FILE1, "<", "$name1") or die;                        #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/9_1_NORMALIZE/$name2") or die;   #以写的方式打开存放结果的文件。

while (my $line=<FILE1>) {
    $line =~ s/ ^\s* ([^\s][^\n]+[^\s]) \s*$ /$1/x or die;   #去掉首尾的空格。
    my @array1=split(/ +/,$line);
    my $t1 = shift(@array1);    #移去第一个元素，并返回它。
    printf FILE2 "%-20s",$t1;
    my $max1=&MAX(@array1);
    my $min1=&MIN(@array1);
    ($max1 >= $min1) or die;
    for(my $i1=1;$i1<=$#array1;$i1++){
      $array1[$i1] = ($array1[$i1] - $min1)/($max1 - $min1 + 10**(-20));  #归一化,把分母加上一个很小的数，防止除数为0.   
      printf FILE2 "%20.10f",$array1[$i1];	
    }
    print FILE2 "\n";
}                              

close FILE1;
close FILE2;
} ######################################## END NORMALIZE_9_2A








##############################################################################################
sub NORMALIZE_9_2B    #缩放，(x - min)/(max - min),按行。
##############################################################################################
{ 
my $name1 = $_[0];              #接受参数，待处理文件的名字(含路径)。
$name1 =~ m|/([^/]+)$| or die;   
my $name2 = $1;                 #输出文件的名字。

if (-e "$whole/9_2_NORMALIZE")         { }else{mkdir "$whole/9_2_NORMALIZE"   or die;}                                         
open(FILE1, "<", "$name1") or die;                        #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/9_2_NORMALIZE/$name2") or die;   #以写的方式打开存放结果的文件。

while (my $line=<FILE1>) {
    $line =~ s/ ^\s* ([^\s][^\n]+[^\s]) \s*$ /$1/x or die;   #去掉首尾的空格。
    my @array1=split(/ +/,$line);
    my $t1 = shift(@array1);    #移去第一个元素，并返回它。
    my $t2 = shift(@array1);    #移去第一个元素，并返回它。
    printf FILE2 "%-20s",$t1;
    printf FILE2 "%-20s",$t2;
    my $max1=&MAX(@array1);
    my $min1=&MIN(@array1);
    ($max1 >= $min1) or die;
    for(my $i1=2;$i1<=$#array1;$i1++){
      $array1[$i1] = ($array1[$i1] - $min1)/($max1 - $min1 + 10**(-20));  #归一化,把分母加上一个很小的数，防止除数为0.   
      printf FILE2 "%20.10f",$array1[$i1];	
    }
    print FILE2 "\n";
}                              

close FILE1;
close FILE2;
} ######################################## END NORMALIZE_9_2B










##############################################################################################
sub NORMALIZE_9_3A    #缩放，(x - min)/(max - min),按列。#最好是按这种方法归一化。
##############################################################################################
{ 
my $name1 = $_[0];              #接受参数，待处理文件的名字(含路径)。
$name1 =~ m|/([^/]+)$| or die;   
my $name2 = $1;                 #输出文件的名字。

if (-e "$whole/9_1_NORMALIZE")         { }else{mkdir "$whole/9_1_NORMALIZE"   or die;}                                         
open(FILE1, "<", "$name1") or die;                        #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/9_1_NORMALIZE/$name2") or die;   #以写的方式打开存放结果的文件。
                    
my @lines = <FILE1>;
my @array2d = '';      
my @column = ();
                                            #声明一个数组,此时数组的第一个元素是有值的，为空值。若使用“my @array2d = ();”，则$array2d[0]的值也为undef,即所有元素都是没有值的。
for (my $i1=0;$i1<=$#lines;$i1++) {                                #数组@lines的第一个元素不做处理。
  $lines[$i1] =~ s/ ^\s* ([^\s][^\n]+[^\s]) \s*$ /$1/x or die;     #去掉首尾的空格。
  $array2d[$i1] = [split(/ +/,$lines[$i1])];                       #用方括号生成匿名引用。
  $#{$array2d[$i1]} >= 1 or die;                                   #看成一维数组时，其每个元素均为数组。
  $column[$i1] = $#{$array2d[$i1]};
}	                                                           #现在@array2d是一个2维数组。

for (my $i1=0; $i1<=$#column - 1; $i1++) {  
  $column[$i1] == $column[$i1+1]   or die;
}


for (my $j=1;$j<=$column[0];$j++) {    #一次处理一列。
  my @temp2 = ();
  for(my $k=0; $k<=$#lines; $k++){                   
    $temp2[++$#temp2] = $array2d[$k][$j];            #同一列的氨基酸。
  }
  my $min = &MIN(@temp2);                            
  my $max = &MAX(@temp2);                            
  for(my $n=0; $n<=$#lines; $n++){                  
    $array2d[$n][$j] = ($array2d[$n][$j] - $min) / ($max -$min+10**(-20));     #归一化。
  }
}

#归一化完毕，下面是输出:

for (my $i=0;$i<=$#lines;$i++) {  
  printf FILE2 "%20s",$array2d[$i][0];	
  for (my $j=1;$j<=$column[0];$j++) {
    printf FILE2 "%20.10f",$array2d[$i][$j];
  }
  print FILE2 "\n";
}

close FILE1;
close FILE2;
}##################################################### END NORMALIZE_9_3A









##############################################################################################
sub NORMALIZE_9_3B    #缩放，(x - min)/(max - min),按列。#最好是按这种方法归一化。
##############################################################################################
{ 
my $name1 = $_[0];              #接受参数，待处理文件的名字(含路径)。
$name1 =~ m|/([^/]+)$| or die;   
my $name2 = $1;                 #输出文件的名字。

if (-e "$whole/9_2_NORMALIZE")         { }else{mkdir "$whole/9_2_NORMALIZE"   or die;}                                         
open(FILE1, "<", "$name1") or die;                        #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/9_2_NORMALIZE/$name2") or die;   #以写的方式打开存放结果的文件。
                    
my @lines = <FILE1>;
my @array2d = ''; 
my @column = ();
                                                 #声明一个数组,此时数组的第一个元素是有值的，为空值。若使用“my @array2d = ();”，则$array2d[0]的值也为undef,即所有元素都是没有值的。
for (my $i1=0;$i1<=$#lines;$i1++) {                                #数组@lines的第一个元素不做处理。
  $lines[$i1] =~ s/ ^\s* ([^\s][^\n]+[^\s]) \s*$ /$1/x or die;    #去掉首尾的空格。
  $array2d[$i1] = [split(/ +/,$lines[$i1])];                       #用方括号生成匿名引用。
  $#{$array2d[$i1]} >= 1 or die;                                   #看成一维数组时，其每个元素均为数组。
  $column[$i1] = $#{$array2d[$i1]};
}	                                                           #现在@array2d是一个2维数组。

for (my $i1=0; $i1<=$#column - 1; $i1++) {  
  $column[$i1] == $column[$i1+1]   or die;
}

for (my $j=2;$j<=$column[0];$j++) {    #一次处理一列。
  my @temp2 = ();
  for(my $k=0; $k<=$#lines; $k++){                   
    $temp2[++$#temp2] = $array2d[$k][$j];            #同一列的氨基酸。
  }
  my $min = &MIN(@temp2);                            
  my $max = &MAX(@temp2);                            
  for(my $n=0; $n<=$#lines; $n++){                  
    $array2d[$n][$j] = ($array2d[$n][$j] - $min) / ($max -$min+10**(-20));     #归一化。
  }
}

#归一化完毕，下面是输出:

for (my $i=0;$i<=$#lines;$i++) {  
  printf FILE2 "%20s",$array2d[$i][0];	
  printf FILE2 "%20s",$array2d[$i][1];
  for (my $j=2; $j<=$column[0]; $j++) {
    printf FILE2 "%20.10f",$array2d[$i][$j];
  }
  print FILE2 "\n";
}

close FILE1;
close FILE2;
}##################################################### END NORMALIZE_9_3B









##############################################################################################
sub FORMAT_10_1    #格式: 把向量转化为LibSVM需要的格式。 ####此步由分类器所需要的格式决定。
##############################################################################################
{ 
my $name1 = $_[0];              #接受参数，待处理文件的名字(含路径)。
$name1 =~ m|/([^/]+)$| or die;   
my $name2 = $1;                 #输出文件的名字。

if (-e "$whole/10_1_FORMAT") { }else{mkdir "$whole/10_1_FORMAT"   or die;}                                         
open(FILE1, "<", "$name1") or die;                      #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/10_1_FORMAT/$name2") or die;   #以写的方式打开存放结果的文件。

while(my $line=<FILE1>) {
  $line =~ s/^\s*([^\s][^\n]+[^\s])\s*$/$1/ or die;   
  my @array1=split(/ +/,$line);
  $#array1 >= 2  or die;
  printf FILE2 "%20s" ,$array1[0];
  for (my $i=1;$i<=$#array1;$i++) {
    my $n=$i;
    printf FILE2 "%4s:%-20.10f",$n,$array1[$i];
  }
  print FILE2 "\n";
}

close FILE1;
close FILE2;
}################################################# END FORMAT_10_1








##############################################################################################
sub FORMAT_10_2    #格式: 把向量转化为LibSVM需要的格式。 ####此步由分类器所需要的格式决定。
##############################################################################################
{ 
my $name1 = $_[0];              #接受参数，待处理文件的名字(含路径)。
$name1 =~ m|/([^/]+)$| or die;   
my $name2 = $1;                 #输出文件的名字。

if (-e "$whole/10_2_FORMAT") { }else{mkdir "$whole/10_2_FORMAT"   or die;}                                         
open(FILE1, "<", "$name1") or die;                        #以读的方式打开待处理文件。
open(FILE2, ">", "$whole/10_2_FORMAT/$name2") or die;   #以写的方式打开存放结果的文件。

while(my $line=<FILE1>) {
  $line =~ s/^\s*([^\s][^\n]+[^\s])\s*$/$1/ or die;   
  my @array1=split(/ +/,$line);
  $#array1 >= 3  or die;
  printf FILE2 "%20s", $array1[0];
  printf FILE2 "%20s", $array1[1];
  for (my $i=2;$i<=$#array1;$i++) {
    my $n=$i-1;
    printf FILE2 "%4s:%-20.10f",$n,$array1[$i];
  }
  print FILE2 "\n";
}

close FILE1;
close FILE2;
}################################################# END FORMAT_10_2













#调用子程序。共10步：

#第1步，分割:
##############################################################################################################
{   #括起来，限制变量的作用域。
my $dir1 = 'pssm';
opendir(my $dh1, $dir1)  ||  die "can't opendir $dir1: $!";

while (my $folder = readdir $dh1) {         #读取此目录下的所有目录。
  next unless $folder =~ m/^[0-9-]+$/;      #目录名只含类别标签。
  my $dir2 = "$dir1/$folder";
  opendir(my $dh2, $dir2)  ||  die "can't opendir $dir2: $!";        #打开目录句柄，读取目录里的文件名。
  while (my $file=readdir $dh2) {                                    #读取的仅仅是文件名称，不含路径。
    next unless $file !~ m/^[.]/; 
    next unless $file !~ m/[~]$/;  
    &SEPARATE_1("$dir2/$file");                                      #传递的参数是含有路径的文件名。
  }
}

}
##################################################







#第2步，检查:
##############################################################################################################
{   
my $dir1 = "$whole/1_SEPARATE";
opendir(my $dh1, $dir1)  ||  die "can't opendir $dir1: $!";

while (my $folder = readdir $dh1) {         #读取此目录下的所有目录。
  next unless $folder =~ m/^[0-9-]+$/;      #目录名只含类别标签。
  my $dir2 = "$dir1/$folder";
  opendir(my $dh2, $dir2)  ||  die "can't opendir $dir2: $!";        #打开目录句柄，读取目录里的文件名。
  while (my $file=readdir $dh2) {                                    #读取的仅仅是文件名称，不含路径。
    next unless $file !~ m/^[.]/; 
    next unless $file !~ m/[~]$/;  
    &CHECK_2("$dir2/$file");                                      #传递的参数是含有路径的文件名。
  }
}

}
##################################################






#第3步，统计:
my $stat = 12;
if ($stat == 1) {    #并不需要每次都去运行此子程序。
##############################################################################################################
{  
my $dir1 = "$whole/2_CHECK";
opendir(my $dh1, $dir1)  ||  die "can't opendir $dir1: $!";

while (my $folder = readdir $dh1) {        #读取此目录下的所有目录。
  next unless $folder =~ m/^[0-9-]+$/;     #目录名只含类别标签。
  my $dir2 = "$dir1/$folder";
  opendir(my $dh2, $dir2)  ||  die "can't opendir $dir2: $!";        #打开目录句柄，读取目录里的文件名。
  while (my $file=readdir $dh2) {                                    #读取的仅仅是文件名称，不含路径。
    next unless $file !~ m/^[.]/; 
    next unless $file !~ m/[~]$/;
    &STATISTIC_3_1("$dir2/$file"); 
    &STATISTIC_3_4("$dir2/$file");
    &STATISTIC_3_7("$dir2/$file");
  }
}

&STATISTIC_3_2;
&STATISTIC_3_3;
&STATISTIC_3_5;
&STATISTIC_3_6;
&STATISTIC_3_8;
&STATISTIC_3_9;
}
##################################################
}






#第4步，第1次归一化:
##############################################################################################################
{
my $dir1 = "$whole/2_CHECK";
opendir(my $dh1, $dir1)  ||  die "can't opendir $dir1: $!";

while (my $folder = readdir $dh1) {                                  #读取此目录下的所有目录。
  next unless $folder =~ m/^[0-9-]+$/;                               #目录名只含类别标签。
  my $dir2 = "$dir1/$folder";
  opendir(my $dh2, $dir2)  ||  die "can't opendir $dir2: $!";        #打开目录句柄，读取目录里的文件名。
  while (my $file=readdir $dh2) {                                    #读取的仅仅是文件名称，不含路径。
    next unless $file !~ m/^[.]/; 
    next unless $file !~ m/[~]$/;
    if($ARGV[1]==0) { &NORMALIZE_4_0("$dir2/$file");} 
    if($ARGV[1]==1) { &NORMALIZE_4_1("$dir2/$file");}               #传递的参数是含有路径的文件名。
    if($ARGV[1]==2) { &NORMALIZE_4_2("$dir2/$file");}
    if($ARGV[1]==3) { &NORMALIZE_4_3("$dir2/$file");}
    if($ARGV[1]==4) { &NORMALIZE_4_4("$dir2/$file");}
  }
}

}
#######################################################







#第5步，过滤:
##############################################################################################################
{
chomp(my $dir1 = "$whole/4_NORMALIZE");
opendir(my $dh1, $dir1)  ||  die "can't opendir $dir1: $!";

while (my $folder = readdir $dh1) {        #读取此目录下的所有目录。
  next unless $folder =~ m/^[0-9-]+$/;     #目录名只含类别标签。
  my $dir2 = "$dir1/$folder";
  opendir(my $dh2, $dir2)  ||  die "can't opendir $dir2: $!";        #打开目录句柄，读取目录里的文件名。
  while (my $file=readdir $dh2) {          #读取的仅仅是文件名称，不含路径。
    next unless $file !~ m/^[.]/; 
    next unless $file !~ m/[~]$/;
    given($ARGV[3]) {
      when( 1) {&FILTER_5_1 ("$dir2/$file"); }    #传递的参数是含有路径的文件名。
      when( 2) {&FILTER_5_2 ("$dir2/$file"); }
      when( 3) {&FILTER_5_3 ("$dir2/$file"); }
      when( 4) {&FILTER_5_4 ("$dir2/$file"); }
      when( 5) {&FILTER_5_5 ("$dir2/$file"); }
      when( 6) {&FILTER_5_6 ("$dir2/$file"); }
      when( 7) {&FILTER_5_7 ("$dir2/$file"); }
      when( 8) {&FILTER_5_8 ("$dir2/$file"); }
      when( 9) {&FILTER_5_9 ("$dir2/$file"); }
      when(10) {&FILTER_5_10("$dir2/$file"); }
      when(11) {&FILTER_5_11("$dir2/$file"); }
      when(12) {&FILTER_5_12("$dir2/$file"); }
      when(13) {&FILTER_5_13("$dir2/$file"); }
      when(14) {&FILTER_5_14("$dir2/$file"); }
    }
  }
}
}
########################################################
$ave_ratio = $ave_ratio/$all_seq;
my $all_ratio = $remain_line/$all_line;

print  AVE_FILE  "
总共$all_seq条序列。
平均保留比例为：$ave_ratio% 

总共$all_line行(或元素)。
总共保留了$remain_line行(或元素)。
总共过滤掉了$delete_line行(或元素)。
总的保留比例为$all_ratio。
\n";






#第6步，转化为向量:
##############################################################################################################
{
chomp(my $dir1 = "$whole/5_FILTER");
opendir(my $dh1, $dir1)  ||  die "can't opendir $dir1: $!";
while (my $folder = readdir $dh1) {        #读取此目录下的所有目录。
  next unless $folder =~ m/^[0-9-]+$/;     #目录名只含类别标签。
  my $dir2 = "$dir1/$folder";
  opendir(my $dh2, $dir2)  ||  die "can't opendir $dir2: $!";        #打开目录句柄，读取目录里的文件名。
  while (my $file=readdir $dh2) {          #读取的仅仅是文件名称，不含路径。。
    next unless $file !~ m/^[.]/; 
    next unless $file !~ m/[~]$/;
    if ($ARGV[11]==1) {&VECTOR_6_1("$dir2/$file");}     #传递的参数是含有路径的文件名。
    if ($ARGV[11]==2) {&VECTOR_6_2("$dir2/$file");} 
    if ($ARGV[11]==3) {&VECTOR_6_3("$dir2/$file");} 
    if ($ARGV[11]==4) {&VECTOR_6_4("$dir2/$file");} 
    if ($ARGV[11]==5) {&VECTOR_6_5("$dir2/$file");}    
  }
}
}
########################################################







#第7步，进一步处理向量:
##############################################################################################################
{
chomp(my $dir1 = "$whole/6_VECTOR");
opendir(my $dh1, $dir1)  ||  die "can't opendir $dir1: $!";
while (my $folder = readdir $dh1) {        #读取此目录下的所有目录。
  next unless $folder =~ m/^[0-9-]+$/;     #目录名只含类别标签。
  my $dir2 = "$dir1/$folder";
  opendir(my $dh2, $dir2)  ||  die "can't opendir $dir2: $!";        #打开目录句柄，读取目录里的文件名。
  while (my $file=readdir $dh2) {          #读取的仅仅是文件名称，不含路径。
    next unless $file !~ m/^[.]/; 
    next unless $file !~ m/[~]$/;
    if ($ARGV[11] =~ m/^[1345]$/ ) {&CONVERT_7_1("$dir2/$file");}     #传递的参数是含有路径的文件名。
    if ($ARGV[11] =~ m/^[2]$/    ) {&CONVERT_7_2("$dir2/$file");}  
  }
}
}
########################################################






#第8步，合并: 把各类的向量分别放在不同的文件里，同一类的放在同一个文件。还要把所有向量放在一个文件里:
##############################################################################################################
{
chomp(my $dir1 = "$whole/7_CONVERT");
opendir(my $dh1, $dir1)  ||  die "can't opendir $dir1: $!";
while (my $folder = readdir $dh1) {        #读取此目录下的所有目录。
  next unless $folder =~ m/^[0-9-]+$/;     #目录名只含类别标签。
  my $dir2 = "$dir1/$folder";
  opendir(my $dh2, $dir2)  ||  die "can't opendir $dir2: $!";        #打开目录句柄，读取目录里的文件名。
  while (my $file=readdir $dh2) {          #读取的仅仅是文件名称，不含路径。
    next unless $file !~ m/^[.]/; 
    next unless $file !~ m/[~]$/;
    &MERGE_8_1("$dir2/$file");
    &MERGE_8_2("$dir2/$file");  
  }
}
}
########################################################








#第9步，缩放: 把PSSM矩阵的元素归一化到[0，1]或不归一，4种情况。(第二次缩放):
##############################################################################################################
{
chomp(my $dir1 = "$whole/8_1_MERGE");
opendir(my $dh1, $dir1)  ||  die "can't opendir $dir1: $!";
while (my $file=readdir $dh1) {          #读取的仅仅是文件名称，不含路径。
  next unless $file !~ m/^[.]/; 
  next unless $file !~ m/[~]$/;
  if($ARGV[23]==0) {&NORMALIZE_9_0A("$dir1/$file");} 
  if($ARGV[23]==1) {&NORMALIZE_9_1A("$dir1/$file");} 
  if($ARGV[23]==2) {&NORMALIZE_9_2A("$dir1/$file");} 
  if($ARGV[23]==3) {&NORMALIZE_9_3A("$dir1/$file");} 
}
}


{
chomp(my $dir1 = "$whole/8_2_MERGE");
opendir(my $dh1, $dir1)  ||  die "can't opendir $dir1: $!";
while (my $file=readdir $dh1) {          #读取的仅仅是文件名称，不含路径。
  next unless $file !~ m/^[.]/; 
  next unless $file !~ m/[~]$/;
  if($ARGV[23]==0) {&NORMALIZE_9_0B("$dir1/$file");} 
  if($ARGV[23]==1) {&NORMALIZE_9_1B("$dir1/$file");} 
  if($ARGV[23]==2) {&NORMALIZE_9_2B("$dir1/$file");} 
  if($ARGV[23]==3) {&NORMALIZE_9_3B("$dir1/$file");} 
}
}


########################################################






#第10步，格式: 把向量转化为LibSVM需要的格式:
##############################################################################################################
{
chomp(my $dir1 = "$whole/9_1_NORMALIZE");
opendir(my $dh1, $dir1)  ||  die "can't opendir $dir1: $!";
while (my $file=readdir $dh1) {          #读取的仅仅是文件名称，不含路径。
  next unless $file !~ m/^[.]/; 
  next unless $file !~ m/[~]$/;
  &FORMAT_10_1("$dir1/$file");
}
}


{
chomp(my $dir1 = "$whole/9_2_NORMALIZE");
opendir(my $dh1, $dir1)  ||  die "can't opendir $dir1: $!";
while (my $file=readdir $dh1) {          #读取的仅仅是文件名称，不含路径。
  next unless $file !~ m/^[.]/; 
  next unless $file !~ m/[~]$/;
  &FORMAT_10_2("$dir1/$file");
}
}
########################################################







