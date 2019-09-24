rule CalculatedNdSNormal:
    input:
        oesophaguspatientinfo="data/oesophagus/patient_info.xlsx",
        oesophagusdata="data/oesophagus/esophagus.csv",
        skindata="data/skin/skin_mutation.csv"
    output:
        oesophagusdnds="results/oesophagus/dnds.csv",
        oesophagusdndsgenes="results/oesophagus/dnds_genes.csv",
        skindnds="results/skin/dnds.csv",
        skindndsgenes="results/skin/dnds_genes.csv",
        oesophagusdndsneutral="results/oesophagus/dnds_neutral.csv",
        oesophagusdndsgenesneutral="results/oesophagus/dnds_genes_neutral.csv",
        singlepatientdnds="results/dataforfigures/singlepatient_bins.csv",
        singlepatientdndsgenes="results/dataforfigures/singlepatient_bins_genes.csv"
    params:
        singlepatient=config["patient"],
        step=config["idndslimits"]["step"],
        minarea=config["idndslimits"]["minarea"],
        maxarea=config["idndslimits"]["maxarea"],
        singularityimage=config["stansingularity"]
    shell:
        """
        module unload python
        module load singularity
        singularity exec {params.singularityimage} \
        Rscript R/CalculatedNdS-normal.R \
            --patientinfo {input.oesophaguspatientinfo} \
            --oesophagusdata {input.oesophagusdata} \
            --skindata {input.skindata} \
            --oesophagusdnds {output.oesophagusdnds} \
            --oesophagusdndsgenes {output.oesophagusdndsgenes} \
            --skindnds {output.skindnds} \
            --skindndsgenes {output.skindndsgenes} \
            --oesophagusdndsneutral {output.oesophagusdndsneutral} \
            --oesophagusdndsgenesneutral {output.oesophagusdndsgenesneutral} \
            --singlepatient {params.singlepatient} \
            --singlepatientdnds {output.singlepatientdnds} \
            --singlepatientdndsgenes {output.singlepatientdndsgenes} \
            --step {params.step} \
            --minarea {params.minarea} \
            --maxarea {params.maxarea}
        """

rule CalculatedNdSNormalSNV:
    input:
        oesophaguspatientinfo="data/oesophagus/patient_info.xlsx",
        oesophagusdata="data/oesophagus/esophagus.csv",
        skindata="data/skin/skin_mutation.csv"
    output:
        oesophagusdnds="results/oesophagus/dnds_snv.csv",
        oesophagusdndsgenes="results/oesophagus/dnds_genes_snv.csv",
        skindnds="results/skin/dnds_snv.csv",
        skindndsgenes="results/skin/dnds_genes_snv.csv",
        oesophagusdndsneutral="results/oesophagus/dnds_neutral_snv.csv",
        oesophagusdndsgenesneutral="results/oesophagus/dnds_genes_neutral_snv.csv",
        singlepatientdnds="results/dataforfigures/singlepatient_bins_snv.csv",
        singlepatientdndsgenes="results/dataforfigures/singlepatient_bins_genes_snv.csv"
    params:
        singlepatient=config["patient"],
        step=config["idndslimits"]["step"],
        minarea=config["idndslimits"]["minarea"],
        maxarea=config["idndslimits"]["maxarea"],
        singularityimage=config["stansingularity"]
    shell:
        """
        module unload python
        module load singularity
        singularity exec {params.singularityimage} \
        Rscript R/CalculatedNdS-normal-snv.R \
            --patientinfo {input.oesophaguspatientinfo} \
            --oesophagusdata {input.oesophagusdata} \
            --skindata {input.skindata} \
            --oesophagusdnds {output.oesophagusdnds} \
            --oesophagusdndsgenes {output.oesophagusdndsgenes} \
            --skindnds {output.skindnds} \
            --skindndsgenes {output.skindndsgenes} \
            --oesophagusdndsneutral {output.oesophagusdndsneutral} \
            --oesophagusdndsgenesneutral {output.oesophagusdndsgenesneutral} \
            --singlepatient {params.singlepatient} \
            --singlepatientdnds {output.singlepatientdnds} \
            --singlepatientdndsgenes {output.singlepatientdndsgenes} \
            --step {params.step} \
            --minarea {params.minarea} \
            --maxarea {params.maxarea}
        """

rule MakeFilesSSB:
    output:
        directory("results/oesophagus/SSBfiles/{oes_sample}/")
    input:
        oesophaguspatientinfo="data/oesophagus/patient_info.xlsx",
        oesophagusdata="data/oesophagus/esophagus.csv",
    params:
        singlepatient=config["patient"],
        step=config["idndslimits"]["step"],
        minarea=config["idndslimits"]["minarea"],
        maxarea=config["idndslimits"]["maxarea"]
    shell:
        """
        module load R
        module load gcc
        module load julia
        Rscript R/SSB-dNdS-files.R \
            --patientinfo {input.oesophaguspatientinfo} \
            --oesophagusdata {input.oesophagusdata} \
            --outputdir {output} \
            --sample {wildcards.oes_sample} \
            --step {params.step} \
            --minarea {params.minarea} \
            --maxarea {params.maxarea}
        """

rule CalculateSitedNdSNormal:
    input:
        oesophaguspatientinfo="data/oesophagus/patient_info.xlsx",
        oesophagusdata="data/oesophagus/esophagus.csv",
    output:
        oesophagusdnds="results/oesophagus/sitednds.csv",
        oesophagusdndsgenes="results/oesophagus/sitednds_genes.csv",
        oesophagushotspots="results/oesophagus/sitednds_genes_hotspots.csv",
    params:
        step=config["idndslimits"]["step"],
        minarea=config["idndslimits"]["minarea"],
        maxarea=config["idndslimits"]["maxarea"],
        singularityimage=config["stansingularity"]
    shell:
        """
        module unload python
        module load singularity
        singularity exec {params.singularityimage} \
        Rscript R/CalculatesitedNdS-normal.R \
            --patientinfo {input.oesophaguspatientinfo} \
            --oesophagusdata {input.oesophagusdata} \
            --oesophagusdnds {output.oesophagusdnds} \
            --oesophagusdndsgenes {output.oesophagusdndsgenes} \
            --hotspots {output.oesophagushotspots} \
            --step {params.step} \
            --minarea {params.minarea} \
            --maxarea {params.maxarea}
        """
