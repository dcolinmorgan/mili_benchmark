recap_ttest <- read.delim("~/recap_ttest.txt")
library(dplyr)

wgplusP=filter(wgplus,region=='promoter')
wgplusB=filter(wgplus,region=='body')
wgminusB=filter(wgminus,region=='body')
wgminusP=filter(wgminus,region=='promoter')
arrayplusP=filter(arrayplus,region=='promoter')
arrayplusB=filter(arrayplus,region=='body')
arrayminusB=filter(arrayminus,region=='body')
arrayminusP=filter(arrayminus,region=='promoter')

t.test(arrayminusP$AUROC,arrayplusP$AUROC,alternative = c("greater"))
t.test(arrayminusB$AUROC,arrayplusB$AUROC,alternative = c("greater"))

t.test(wgminusP$AUROC,wgplusP$AUROC,alternative = c("greater"))
t.test(wgminusB$AUROC,wgplusB$AUROC,alternative = c("greater"))



