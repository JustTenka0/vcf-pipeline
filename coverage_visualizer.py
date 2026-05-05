import webbrowser
import urllib.parse

def open_ucsc_trio(case_id, github_user, repo_name, branch="main", genome="hg38"):
    base_raw_url = f"https://raw.githubusercontent.com/{github_user}/{repo_name}/{branch}"
    folder = f"/coverage_data"
    
    files = {
        "child": f"{base_raw_url}/{folder}/{case_id}_child.bg",
        "father": f"{base_raw_url}/{folder}/{case_id}_father.bg",
        "mother": f"{base_raw_url}/{folder}/{case_id}_mother.bg",
        "vcf": f"{base_raw_url}/vcf_data/{case_id}.vep_filtered.vcf"
    }

    # Definiamo le tracce per UCSC
    # Nota: il file VCF su UCSC richiede spesso l'indice (.tbi) nella stessa cartella
    tracks = [
        f"track type=bedGraph name='{case_id}_Child' visibility=full color=0,0,255 bigDataUrl={files['child']}",
        f"track type=bedGraph name='{case_id}_Father' visibility=full color=255,0,0 bigDataUrl={files['father']}",
        f"track type=bedGraph name='{case_id}_Mother' visibility=full color=0,255,0 bigDataUrl={files['mother']}",
        f"track type=vcf name='{case_id}_Variants' visibility=pack bigDataUrl={files['vcf']}"
    ]

    # Codifica l'URL
    combined_tracks = "\n".join(tracks)
    encoded_tracks = urllib.parse.quote(combined_tracks)
    
    # Crea l'URL finale di UCSC
    ucsc_url = f"https://genome.ucsc.edu/cgi-bin/hgTracks?db={genome}&hgct_customText={encoded_tracks}"
    
    print(f"Opening UCSC for Case: {case_id}")
    webbrowser.open(ucsc_url)

# --- CONFIGURAZIONE ---
USER = "JustTenka0"
REPO = "vcf-pipeline"
CASES = ["trio_1", "trio_2", "trio_3", "trio_4", "trio_5"] # Aggiungi qui i tuoi ID

for case in CASES:
    open_ucsc_trio(case, USER, REPO)