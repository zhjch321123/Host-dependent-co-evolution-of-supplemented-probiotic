#绘制基因组圈图加SNP位点
方法如下：
首先，需要有组装成一条基因组的文件 “line.fasta”
而后运行基因预测，获得".gff"文件
prodigal -a line.pep -d line.cds -f gff -g 11  -o line.gff -p single -s sample1.stat -i line.fasta > ref.log

1.将参考基因组原始测序序列文件与组装后的line.fasta进行bowtie2比对，获得“.bam“文件，命令类似如下：
bowtie2 -p 24 -x /userdatanode4/data_z/ISP_metagenome/mdis_snp/verify/ref/ref --no-mixed --very-sensitive --n-ceil 0,0.01 -1 /userdatanode4/data_z/ISP_metagenome/mdis_snp/verify/VA2_1.fq -2 /userdatanode4/data_z/ISP_metagenome/mdis_snp/verify/VA2_2.fq | samtools sort -O bam -@ 12 -o - > /userdatanode4/data_z/ISP_metagenome/mdis_snp/verify/VA2.bam
2. 用samtools统计
samtools depth ref.bam > samtools_depth
3.#再使用自编脚本 depth_base_stat.py，读取上一步统计的文件 samtools_depth，以及参考基因组 fasta 文件
#根据设定滑窗大小（2000bp），统计参考序列各滑窗区间中的各碱基占比以及平均 reads 覆盖深度
python3 depth_base_stat.py -g line.fasta -d samtools_depth -s HNU082.depth_base.txt -l 2000

加载R包，指定文件,主要以托R语言的circlize包完成，下面进入正式程序：
setwd("C:/Users/zhjch/Desktop")
library(stringr) #方便处理字符串
library(circlize)      #绘制圈图
library(grid)    #组合画板，圈图 + 图例

##命令传递
sample_name <- 'HNU082' #测序样本名称
ref_name <- 'LP_HNU082'    #参考基因组名称

genome_gff <- 'LP_HNU082.gff'  #参考基因组 gff 注释文件
snp_vcf <- 'snp_line.vcf'    #SNP 检测结果 vcf 文件
indel_vcf <- 'snp_selevt.vcf' #用筛选验证过的SNP代替indel
depth_base_stat <- 'HNU082.depth_base.txt'    #测序深度、碱基含量统计结果文件
seq_split <- 2000    #滑窗大小，与 depth_base_stat 中使用的滑窗大小对应
 
out_dir = 'output'    #生成一个目录，用于存放结果文件
if (!file.exists(out_dir)) dir.create(out_dir)

##参考基因组长度、GC 统计 & 测序覆盖度、深度统计
depth_base <- read.delim(depth_base_stat, stringsAsFactors = FALSE)
genome_size <- sum(depth_base$seq_end - depth_base$seq_start + 1) 
genome_GC <- round(mean(depth_base$GC), 2)
 
depth_exist <- subset(depth_base, depth != 0)
coverage <- round(100 * sum(depth_exist$seq_end - depth_exist$seq_start + 1) / genome_size, 2)
average_depth <- round(mean(depth_base$depth), 0)
 
seq_stat <- NULL
for (seq_id in unique(depth_base$seq_ID)) seq_stat <- rbind(seq_stat, c(seq_id, 1, max(subset(depth_base, seq_ID == seq_id)$seq_end)))
seq_stat <- data.frame(seq_stat, stringsAsFactors = FALSE)
colnames(seq_stat) <- c('seq_ID', 'seq_start', 'seq_end')
rownames(seq_stat) <- seq_stat$seq_ID
seq_stat$seq_start <- as.numeric(seq_stat$seq_start)
seq_stat$seq_end <- as.numeric(seq_stat$seq_end)
 
write.table(seq_stat, str_c(out_dir, '/', sample_name, '.genome_stat.txt'), row.names = FALSE, sep = '\t', quote = FALSE)

#读取 vcf 文件，统计 SNP 类型
snp <- read.delim(snp_vcf, header = FALSE, colClasses = 'character', comment.char = '#')[c(1, 2, 4, 5)]
snp$V2 <- as.numeric(snp$V2)
snp$change <- str_c(snp$V4, snp$V5)
 
change <- which(snp$change == 'AT')
snp[change,'type1'] <- 'A>T|T>A'; snp[change,'type2'] <- 'tv'
change <- which(snp$change == 'AG')
snp[change,'type1'] <- 'A>G|T>C'; snp[change,'type2'] <- 'ti'
change <- which(snp$change == 'AC')
snp[change,'type1'] <- 'A>C|T>G'; snp[change,'type2'] <- 'tv'
 
change <- which(snp$change == 'TA')
snp[change,'type1'] <- 'A>T|T>A'; snp[change,'type2'] <- 'tv'
change <- which(snp$change == 'TG')
snp[change,'type1'] <- 'A>C|T>G'; snp[change,'type2'] <- 'tv'
change <- which(snp$change == 'TC')
snp[change,'type1'] <- 'A>G|T>C'; snp[change,'type2'] <- 'ti'
 
change <- which(snp$change == 'GA')
snp[change,'type1'] <- 'G>A|C>T'; snp[change,'type2'] <- 'ti'
change <- which(snp$change == 'GT')
snp[change,'type1'] <- 'G>T|C>A'; snp[change,'type2'] <- 'tv'
change <- which(snp$change == 'GC')
snp[change,'type1'] <- 'G>C|C>G'; snp[change,'type2'] <- 'tv'
 
change <- which(snp$change == 'CA')
snp[change,'type1'] <- 'G>T|C>A'; snp[change,'type2'] <- 'tv'
change <- which(snp$change == 'CT')
snp[change,'type1'] <- 'G>A|C>T'; snp[change,'type2'] <- 'ti'
change <- which(snp$change == 'CG')
snp[change,'type1'] <- 'G>C|C>G'; snp[change,'type2'] <- 'tv'
 
snp_ti <- length(which(snp$type2 == 'ti'))
snp_tv <- length(which(snp$type2 == 'tv'))
 
snp_at <- length(which(snp$type1 == 'A>T|T>A'))
snp_ag <- length(which(snp$type1 == 'A>G|T>C'))
snp_ac <- length(which(snp$type1 == 'A>C|T>G'))
snp_ga <- length(which(snp$type1 == 'G>A|C>T'))
snp_gt <- length(which(snp$type1 == 'G>T|C>A'))
snp_gc <- length(which(snp$type1 == 'G>C|C>G'))

#统计 SNP 密度
snp <- snp[c(1, 2, 5, 6, 7)]
colnames(snp)[1:2] <- c('seq_ID', 'seq_site')
 
snp_stat <- NULL
seq_ID <- unique(snp$seq_ID)
 
for (seq_ID_n in seq_ID) {
    snp_subset <- subset(snp, seq_ID == seq_ID_n)
    seq_end <- seq_split
    snp_num <- 0
    
    for (i in 1:nrow(snp_subset)) {
        if (snp_subset[i,'seq_site'] <= seq_end) snp_num <- snp_num + 1
        else {
            snp_stat <- rbind(snp_stat, c(seq_ID_n, seq_end - seq_split + 1, seq_end, snp_num))
            
            seq_end <- seq_end + seq_split
            snp_num <- 0
            while (snp_subset[i,'seq_site'] > seq_end) {
                snp_stat <- rbind(snp_stat, c(seq_ID_n, seq_end - seq_split + 1, seq_end, snp_num))
                seq_end <- seq_end + seq_split
            }
            snp_num <- snp_num + 1
        }
    }
    
    while (seq_end < seq_stat[seq_ID_n,'seq_end']) {
        snp_stat <- rbind(snp_stat, c(seq_ID_n, seq_end - seq_split + 1, seq_end, snp_num))
        seq_end <- seq_end + seq_split
        snp_num <- 0
    }
    snp_stat <- rbind(snp_stat, c(seq_ID_n, seq_end - seq_split + 1, seq_stat[seq_ID_n,'seq_end'], snp_num))
}
 
snp_stat <- data.frame(snp_stat, stringsAsFactors = FALSE)
names(snp_stat) <- c('seq_ID', 'seq_start', 'seq_end', 'snp_num')
snp_stat$seq_start <- as.numeric(snp_stat$seq_start)
snp_stat$seq_end <- as.numeric(snp_stat$seq_end)
snp_stat$snp_num <- as.numeric(snp_stat$snp_num)
 
write.table(snp_stat, str_c(out_dir, '/', sample_name, '.snp_stat.txt'), row.names = FALSE, sep = '\t', quote = FALSE)

#读取 vcf 文件，统计 InDel 长度
indel <- read.delim(indel_vcf, header = FALSE, colClasses = 'character', comment.char = '#')[c(1, 2, 4, 5)]
indel$V2 <- as.numeric(indel$V2)
indel$length <- str_length(indel[ ,4]) - str_length(indel[ ,3])
indel_insert <- length(which(indel$length > 0))
indel_delet <- length(which(indel$length < 0))

#统计 InDel 密度
indel <- indel[c(1, 2, 5)]
colnames(indel)[1:2] <- c('seq_ID', 'seq_site')
 
indel_stat <- NULL
seq_ID <- unique(indel$seq_ID)
for (seq_ID_n in seq_ID) {
    indel_subset <- subset(indel, seq_ID == seq_ID_n)
    seq_end <- seq_split
    indel_num <- 0
    
    for (i in 1:nrow(indel_subset)) {
        if (indel_subset[i,'seq_site'] <= seq_end) indel_num <- indel_num + 1
        else {
            indel_stat <- rbind(indel_stat, c(seq_ID_n, seq_end - seq_split + 1, seq_end, indel_num))
            
            seq_end <- seq_end + seq_split
            indel_num <- 0
            while (indel_subset[i,'seq_site'] > seq_end) {
                indel_stat <- rbind(indel_stat, c(seq_ID_n, seq_end - seq_split + 1, seq_end, indel_num))
                seq_end <- seq_end + seq_split
            }
            indel_num <- indel_num + 1
        }
    }
    
    while (seq_end < seq_stat[seq_ID_n,'seq_end']) {
        indel_stat <- rbind(indel_stat, c(seq_ID_n, seq_end - seq_split + 1, seq_end, indel_num))
        seq_end <- seq_end + seq_split
        indel_num <- 0
    }
    indel_stat <- rbind(indel_stat, c(seq_ID_n, seq_end - seq_split + 1, seq_stat[seq_ID_n,'seq_end'], indel_num))
}
 
indel_stat <- data.frame(indel_stat, stringsAsFactors = FALSE)
names(indel_stat) <- c('seq_ID', 'seq_start', 'seq_end', 'indel_num')
indel_stat$seq_start <- as.numeric(indel_stat$seq_start)
indel_stat$seq_end <- as.numeric(indel_stat$seq_end)
indel_stat$indel_num <- as.numeric(indel_stat$indel_num)

##然后接下来就是使用circlize包绘制circos圈图了

pdf(str_c(out_dir, '/', sample_name, '.circlize.pdf'), width = 14, height = 8)
circle_size = unit(1, "snpc")
circos.par(gap.degree = 2)
circos.genomicInitialize(seq_stat, plotType = 'axis')
 
circos.track(
    ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = '#BCEE68',
    panel.fun = function(x, y) {
        xlim = CELL_META$xlim
        ylim = CELL_META$ylim
        seq_ID = CELL_META$sector.index
    } )

#GC% 含量图
circos.genomicTrack(
    depth_base[c(1:3, 5)], track.height = 0.08, bg.col = '#EEEEEE6E', bg.border = NA,
    panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, col = 'blue', lwd = 0.35, ...)
        circos.lines(c(0, max(region)), c(genome_GC, genome_GC), col = 'blue2', lwd = 0.15, lty = 2)
        circos.yaxis(labels.cex = 0.2, lwd = 0.1, tick.length = convert_x(0.15, 'mm'))
    } )
	
#覆盖度 & 深度
circos.genomicTrack(
    depth_base[1:4], track.height = 0.08, ylim = c(0, (max(depth_base$depth) + 1)), bg.col = '#EEEEEE6E', bg.border = NA,
    panel.fun = function(region,value, ...) {
        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, border = 'red2', lwd = 0.02, col = 'red', ...)
        circos.lines(c(0, max(region)), c(average_depth, average_depth), col = 'red3', lwd = 0.15, lty = 2)
        circos.yaxis(labels.cex = 0.2, lwd = 0.1, tick.length = convert_x(0.15, 'mm'))
    } )

#SNP 密度
value_max <- max(snp_stat$snp_num)
colorsChoice <- colorRampPalette(c('white', '#245B8E'))
color_assign <- colorRamp2(breaks = c(0:value_max), col = colorsChoice(value_max + 1))
 
circos.genomicTrackPlotRegion(
    snp_stat, track.height = 0.12, stack = TRUE, bg.border = NA,
    panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
    } )

#InDel 密度
value_max <- max(indel_stat$indel_num)
colorsChoice <- colorRampPalette(c('white', 'red2'))
color_assign <- colorRamp2(breaks = c(0:value_max), col = colorsChoice(value_max + 1))
 
circos.genomicTrackPlotRegion(
    indel_stat, track.height = 0.18, stack = TRUE, bg.border = NA,
    panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
    } )

dev.off()

#获得完成图
