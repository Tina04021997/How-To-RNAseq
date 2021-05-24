#!/bin/bash
# Author: Tina Yang
# Date: Apr 21, 2021
# Download data from GEO

for i in $(seq 79 1 81)
do
  /home/ldap_tina/sratoolkit.2.11.0-ubuntu64/bin/fasterq-dump SRR99733$i -O /LVM_data/tina/RNAseq/data/
done



for i in $(seq 85 1 87)
do 
  /home/ldap_tina/sratoolkit.2.11.0-ubuntu64/bin/fasterq-dump SRR99733$i -O /LVM_data/tina/RNAseq/data/
done
