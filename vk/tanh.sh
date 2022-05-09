echo "" | awk '{n=1000; for(i=0;i<=n;i++){y=-100.0+i*200.0/n; print y,"\t", 1.0, "\t",(1+.5*((exp(y)-exp(-y))/(exp(y)+exp(-y)))),"\t",1.0;}}'
