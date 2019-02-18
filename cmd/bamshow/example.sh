curl -O http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam
curl -O http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam.bai
./bamshow -r chr1 -s 1112420 -e 1112478 wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam
