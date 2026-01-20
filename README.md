## 基于三代测序数据的遗传变异识别及新生抗原筛选研究
### 一、从2022年 HCC1395金标准里看二代三代结构变异结果衍生新抗原的关系
#### 1.筛选SV
* 三代检测出的（不考虑在二代中有没有检测到）
  
        # 进入文件目录
        cd /data/renweijie/data/HCC1395/HCC1395_truth/
        # 1. Pacbio和nanopore都为1的结果（置信度更高一点）
        awk -F'\t' 'NR==1 || ($9==1 && $10==1)' 13059_2022_2816_MOESM3_ESM.txt > pacbio_nanopore_both_1.txt
        #查看总共有多少个（1777）
        wc -l /data/renweijie/data/HCC1395/HCC1395_truth/pacbio_nanopore_both_1.txt
        # 2. pacbio nanopore至少有一个为1的结果
        awk -F'\t' 'NR==1 || ($9==1 || $10==1)' 13059_2022_2816_MOESM3_ESM.txt > pacbio_nanopore_at_least_one_1.txt
        #查看总共有多少个（3552）
        wc -l /data/renweijie/data/HCC1395/HCC1395_truth/pacbio_nanopore_at_least_one_1.txt
* 二代检测出的（不考虑在二代中有没有检测到）

        # 3. illumina 10x至少有一个为1的结果
        awk -F'\t' 'NR==1 || ($12==1 || $13==1)' cleaned_format.txt > illumina_10x_at_least_one_1.txt
        #查看总共有多少个（2152）
        wc -l /data/renweijie/data/HCC1395/HCC1395_truth/illumina_10x_at_least_one_1.txt
        # 4. illumina和10x都为1的结果
        awk -F'\t' 'NR==1 || ($12==1 && $13==1)' cleaned_correct.txt > illumina_10x_both_1.txt
        #查看总共有多少个（321）
        wc -l /data/renweijie/data/HCC1395/HCC1395_truth/illumina_10x_both_1.txt
* 二代、三代都能检测出来的

        # 5. pacbio nanopore至少有一个为1，且illumina 10x至少有一个为1
        awk -F'\t' 'NR==1 || (($9==1 || $10==1) && ($12==1 || $13==1))' cleaned_correct.txt > test_both_techs.txt
        #查看总共有多少个（1042）
        wc -l /data/renweijie/data/HCC1395/HCC1395_truth/test_both_techs.txt
* 只有三代能检测出来的
        
        # 6. pacbio nanopore至少有一个为1，且illumina 10x必须都为0
        awk -F'\t' 'NR==1 || (($9==1 || $10==1) && $12==0 && $13==0)' cleaned_correct.txt > test_pacbio_nanopore_only.txt
        #查看总共有多少个（2510）
        wc -l /data/renweijie/data/HCC1395/HCC1395_truth/test_pacbio_nanopore_only.txt
* 只有二代能检测出来的

        # 7. illumina 10x至少有一个为1，且pacbio nanopore必须都为0
        awk -F'\t' 'NR==1 || (($12==1 || $13==1) && $9==0 && $10==0)' cleaned_correct.txt  > illumina_10x_only.txt
        #查看总共有多少个（1110）
        wc -l /data/renweijie/data/HCC1395/HCC1395_truth/illumina_10x_only.txt
#### 2.从转换后的文件中提取neosv需要的列

       # 1. 创建一个名为 "sv_analysis" 的环境，并指定 python 版本
       conda create -n sv_analysis python=3.9 -y
       # 2. 激活这个环境
       conda activate sv_analysis
       # 3. 安装脚本所需的 pandas 库
       conda install pandas -y
       # 4.执行脚本（修改输入输出文件路径）
       python /data/renweijie/data/HCC1395/HCC1395_truth/txt_to_neosv_bedpe.py
#### 3.运行neosv得到二代三代结构变异衍生新抗原的数量关系
* 三代only SV 2510 新抗原 24
* 二代only SV 1110 新抗原 110
* 二代三代共同 SV 2510 新抗原 85
#### 4.绘图
##### （1）二代三代结构变异数量及新抗原韦恩图
* 创建绘图环境并加载库


        # 1. 创建环境（指定 Python 3.9）
        conda create -n plot_env python=3.9 -y
        # 2. 激活环境
        conda activate plot_env
        # 3. 安装绘图必需的库
        pip install matplotlib matplotlib-venn
* 绘图


        cd /data/renweijie/python_plots/2022_HCC1395_4662_rawCallsSV
        python draw_sv_venn.py
        python draw_neoantigen_venn.py
##### （2）不同类型SV衍生新抗原数量图
* 加载库


        pip install pandas seaborn matplotlib
* 绘图
### 二、检查SV结果的准确性
### 1.处理2022年 4662个变异结果
* 筛选进用PacBio数据检测出来的变异


        #将需要的数据提取出来并变成tsv格式
        python /data/renweijie/data/HCC1395/HCC1395_truth/2022_HCC1395_truth_filter/PacBio_only/extract_pacbio_specific.py
        #将tsv格式转换成vcf
        python /data/renweijie/data/HCC1395/HCC1395_truth/2022_HCC1395_truth_filter/PacBio_only/tsv_to_vcf.py
        # 处理severus输出中的BND类型
        python /data/renweijie/Softwares/SV_tools/severus/HCC1395_Somatic_SV_output/preprocess_for_truvari/severus_bnd_converter.py
* SV类型及数量比对


        # 绘制金标准 severus结果 自己结果的柱状图
        cd /data/renweijie/python_plots/Severus_myResults
        python compare_sv_count.py
* 使用truvari前对文件进行处理


        cd /data/renweijie/data/HCC1395/HCC1395_truvari_output/truvari_vcfs
        chmod +x prepare_truvari.sh
        ./prepare_truvari.sh
* 使用Truvari进行比对

        cd /data/renweijie/data/HCC1395/HCC1395_truvari_output
        truvari bench \
        -b /data/renweijie/data/HCC1395/HCC1395_truvari_output/truvari_vcfs/truth_std_final.vcf.gz \
        -c /data/renweijie/data/HCC1395/HCC1395_truvari_output/truvari_vcfs/paper_result_final.vcf.gz \
        --typeignore \
        --dup-to-ins \
        -p 0 \
        -s 30 \
        -S 0 \
        --sizemax 100000000 \
        --passonly \
        -o PacBio_truvari_output

        


        

      

