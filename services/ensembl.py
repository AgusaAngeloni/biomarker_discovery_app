import requests

def get_sequence(
    chrom,
    start,
    end
):

    url = (
        "https://rest.ensembl.org/"
        f"sequence/region/human/"
        f"{chrom}:{start}..{end}:1"
    )

    r = requests.get(
        url,
        headers={
            "Content-Type":"text/plain"
        }
    )

    r.raise_for_status()

    return r.text.upper()
