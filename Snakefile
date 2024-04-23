
rule download_spine_d:
	message:
		"downloading spine D"
	output:
		".dependencies/downloads/spinedlocal_v2.0.tar.gz"
	shell:
		"""
		wget "https://apisz.sparks-lab.org:8443/downloads/Resource/Protein/2_Protein_local_structure_prediction/old_versions/spinedlocal_v2.0.tar.gz" -nc
		"""

rule unzip_spine_d:
    input: ".dependencies/downloads/spinedlocal_v2.0.tar.gz"
    output:
        tar=temp(".dependencies/downloads/spinedlocal_v2.0.tar"),
        directory=directory(".dependencies/programmes/SPINE-D-local")
    shell:
        """
        gunzip -f {input}
        tar -xf {output.tar}
        mv "SPINE-D-local" {output.directory}
        """
        
