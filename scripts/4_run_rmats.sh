#!/bin/sh

docker run -v /Volumes/aulyd_xk/RPL22L1_analysis/data/:/data rmats:turbo01\
	--b1 /data/RPL22_oe_control.txt\
	--b2 /data/RPL22_oe_treatment.txt\
	--gtf /data/Homo_sapiens.GRCh37.75.gtf\
	--od /data/rmats_output/RPL22_oe\
	-t paired\
	--nthread 6\
	--tstat 6

docker run -v /Volumes/aulyd_xk/RPL22L1_analysis/data/:/data rmats:turbo01\
	--b1 /data/RPL22L1_oe_control.txt\
	--b2 /data/RPL22L1_oe_treatment.txt\
	--gtf /data/Homo_sapiens.GRCh37.75.gtf\
	--od /data/rmats_output/RPL22L1_oe\
	-t paired\
	--nthread 6\
	--tstat 6

docker run -v /Volumes/aulyd_xk/RPL22L1_analysis/data/:/data rmats:turbo01\
	--b1 /data/shluc_control.txt\
	--b2 /data/sh704_treatment.txt\
	--gtf /data/Homo_sapiens.GRCh37.75.gtf\
	--od /data/rmats_output/sh704\
	-t paired\
	--nthread 6\
	--tstat 6

docker run -v /Volumes/aulyd_xk/RPL22L1_analysis/data/:/data rmats:turbo01\
	--b1 /data/shluc_control.txt\
	--b2 /data/sh705_treatment.txt\
	--gtf /data/Homo_sapiens.GRCh37.75.gtf\
	--od /data/rmats_output/sh705\
	-t paired\
	--nthread 6\
	--tstat 6