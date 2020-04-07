#!/bin/sh

DATA_PATH="../data/"

docker run -v $DATA_PATH:/data rmats:turbo01\
	--b1 /data/raw/rmats_groups/RPL22_oe_control.txt\
	--b2 /data/raw/rmats_groups/RPL22_oe_treatment.txt\
	--gtf /data/raw/Homo_sapiens.GRCh37.75.gtf\
	--od /data/raw/rmats_output/RPL22_oe\
	-t paired\
	--nthread 6\
	--tstat 6

docker run -v $DATA_PATH:/data rmats:turbo01\
	--b1 /data/raw/rmats_groups/RPL22L1_oe_control.txt\
	--b2 /data/raw/rmats_groups/RPL22L1_oe_treatment.txt\
	--gtf /data/raw/Homo_sapiens.GRCh37.75.gtf\
	--od /data/raw/rmats_output/RPL22L1_oe\
	-t paired\
	--nthread 6\
	--tstat 6

docker run -v $DATA_PATH:/data rmats:turbo01\
	--b1 /data/raw/rmats_groups/shluc_control.txt\
	--b2 /data/raw/rmats_groups/sh704_treatment.txt\
	--gtf /data/raw/Homo_sapiens.GRCh37.75.gtf\
	--od /data/raw/rmats_output/sh704\
	-t paired\
	--nthread 6\
	--tstat 6

docker run -v $DATA_PATH:/data rmats:turbo01\
	--b1 /data/raw/rmats_groups/shluc_control.txt\
	--b2 /data/raw/rmats_groups/sh705_treatment.txt\
	--gtf /data/raw/Homo_sapiens.GRCh37.75.gtf\
	--od /data/raw/rmats_output/sh705\
	-t paired\
	--nthread 6\
	--tstat 6
