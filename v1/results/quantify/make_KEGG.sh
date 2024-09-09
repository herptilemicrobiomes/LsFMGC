perl -p -e 's/^/over1717.1710_/' LsFMGC_DESeq_Basid1717_vs_Basid1710_95_over.KEGG.tsv > Compare.over_under_KEGG.tsv
perl -p -e 's/^/under1717.1710_/' LsFMGC_DESeq_Basid1717_vs_Basid1710_95_under.KEGG.tsv > Compare.over_under_KEGG.tsv
perl -p -e 's/^/over1717.1710_/' LsFMGC_DESeq_Basid1717_vs_Basid1710_95_over.KEGG.tsv > Compare.over_under_KEGG.tsv
perl -p -e 's/^/under1717.1710_/' LsFMGC_DESeq_Basid1717_vs_Basid1710_95_under.KEGG.tsv >> Compare.over_under_KEGG.tsv
perl -p -e 's/^/underctl.1710_/' LsFMGC_DESeq_Ctl_vs_Basid1710_95_under.KEGG.tsv >> Compare.over_under_KEGG.tsv
perl -p -e 's/^/overctl.1710_/' LsFMGC_DESeq_Ctl_vs_Basid1710_95_over.KEGG.tsv >> Compare.over_under_KEGG.tsv
perl -p -e 's/^/overctl.1717_/' LsFMGC_DESeq_Ctl_vs_Basid1717_95_over.KEGG.tsv >> Compare.over_under_KEGG.tsv
perl -p -e 's/^/underctl.1717_/' LsFMGC_DESeq_Ctl_vs_Basid1717_95_under.KEGG.tsv >> Compare.over_under_KEGG.tsv

grep -v Gene_ID Compare.over_under_KEGG.tsv > x
mv x Compare.over_under_KEGG.tsv

